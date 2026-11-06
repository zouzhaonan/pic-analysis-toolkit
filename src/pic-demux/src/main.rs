//! pic-demux: PIC の 1 つの FASTQ の組 (R1/R2) を 1 回だけ読み、barcode でサンプルへ振り分ける。
//!
//! 旧 scripts/awk/demux_prefix.awk (`paste <(R1) <(R2) | awk`) と同じ中身を書く:
//!   - 読みの名前は「R1 と R2 の行を TAB でつないだ行を空白 (space/TAB) で分けた 3 番目の欄
//!     + "_" + barcode + "_" + UMI」、空白、4 番目の欄。配列・+・品質の行はつないだ行の 2 番目の欄。
//!     (普通の Illumina の FASTQ では、R2 の名前 + "_BC_UMI" + " " + R2 の注記、R2 の配列・+・品質)
//!   - barcode / UMI は R1 の配列 (1 番目の欄) の 7〜12 文字目 / 1〜6 文字目。--barcode-first で逆。
//!   - 同じ barcode がシートに 2 回あれば後のほうのサンプル。
//!
//! --threads n: 1 本で使うコアの目安 (圧縮のスレッドは n-2、1 以上。ほかに解凍 2 本と突き合わせ 1 本)。
//! --level: libdeflate の圧縮レベル (既定 6。大きさは pigz/gzip の既定 -6 とほぼ同じ)。
//! 出力: <out-dir>/<sample><suffix>.tmp (読みが 1 つも無いサンプルは作らない)。
//!   --block-kb (既定 1024 KiB) ごとに独立した gzip のメンバーにして、スレッドで並べて圧縮する
//!   (libdeflate)。メンバーをつないだ gzip は RFC 1952 どおりで、gzip/zcat/pigz/seqkit/hisat2 が
//!   そのまま読める (bcl-convert の FASTQ も同じ形)。
//! --nread-file: sample<TAB>読みの数 (シートの全サンプル、最初に出てきた順、0 も書く)。最後に書く。
//!
//! 入力の不正 (gzip が壊れている・途中で切れている、R1 と R2 の行の数が違う、
//! 4 行単位でない、名前の行が @ でない、3 行目が + でない、R2 の配列と品質の長さが違う) は
//! 終了コード 2。

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom, Write};
use std::os::raw::c_void;
use std::path::{Path, PathBuf};
use std::process::exit;
use std::sync::mpsc::{sync_channel, Receiver, SyncSender};
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};

use libdeflate_sys as ld;

type Res<T> = Result<T, String>;

const VERSION: &str = env!("CARGO_PKG_VERSION");
const USAGE: &str = "usage: pic-demux --r1 <R1.fastq.gz> --r2 <R2.fastq.gz> --barcodes <b1,b2,...> \
--samples <s1,s2,...> --out-dir <dir> --suffix <suffix> --nread-file <file> \
[--threads <n>] [--level <1-12>] [--barcode-first] [--block-kb <kb>]";

struct Args {
    r1: PathBuf,
    r2: PathBuf,
    barcodes: Vec<String>,
    samples: Vec<String>,
    out_dir: PathBuf,
    suffix: String,
    nread_file: PathBuf,
    threads: usize,
    level: i32,
    barcode_first: bool,
    block: usize,
}

fn parse_args() -> Res<Args> {
    let mut it = std::env::args().skip(1);
    let (mut r1, mut r2, mut barcodes, mut samples) = (None, None, None, None);
    let (mut out_dir, mut suffix, mut nread_file) = (None, None, None);
    let (mut threads, mut level, mut barcode_first, mut block_kb) = (4usize, 6i32, false, 1024usize);
    while let Some(a) = it.next() {
        let mut val = || it.next().ok_or_else(|| format!("missing value for {a}"));
        match a.as_str() {
            "--version" => {
                println!("pic-demux {VERSION}");
                exit(0);
            }
            "-h" | "--help" => {
                println!("{USAGE}");
                exit(0);
            }
            "--r1" => r1 = Some(PathBuf::from(val()?)),
            "--r2" => r2 = Some(PathBuf::from(val()?)),
            "--barcodes" => barcodes = Some(val()?),
            "--samples" => samples = Some(val()?),
            "--out-dir" => out_dir = Some(PathBuf::from(val()?)),
            "--suffix" => suffix = Some(val()?),
            "--nread-file" => nread_file = Some(PathBuf::from(val()?)),
            "--threads" => threads = val()?.parse().map_err(|_| "--threads: not a number")?,
            "--level" => level = val()?.parse().map_err(|_| "--level: not a number")?,
            "--block-kb" => block_kb = val()?.parse().map_err(|_| "--block-kb: not a number")?,
            "--barcode-first" => barcode_first = true,
            _ => return Err(format!("unknown argument: {a}\n{USAGE}")),
        }
    }
    let need = |o: Option<String>, n: &str| o.ok_or_else(|| format!("{n} is required\n{USAGE}"));
    let barcodes: Vec<String> = need(barcodes, "--barcodes")?.split(',').map(String::from).collect();
    let samples: Vec<String> = need(samples, "--samples")?.split(',').map(String::from).collect();
    if barcodes.len() != samples.len() {
        return Err("--barcodes and --samples must have the same number of items".into());
    }
    if !(1..=12).contains(&level) {
        return Err("--level must be 1-12".into());
    }
    Ok(Args {
        r1: r1.ok_or("--r1 is required")?,
        r2: r2.ok_or("--r2 is required")?,
        barcodes,
        samples,
        out_dir: out_dir.ok_or("--out-dir is required")?,
        suffix: need(suffix, "--suffix")?,
        nread_file: nread_file.ok_or("--nread-file is required")?,
        threads: threads.max(1),
        level,
        barcode_first,
        block: block_kb.max(16) * 1024,
    })
}

// ---------------------------------------------------------------------------
// 解凍 (1 ファイル 1 スレッド)。
// bcl-convert の FASTQ は 1 MiB ずつのメンバーの gzip なので、メンバーごとに libdeflate で
// 一気に解く。メンバーが大きい (pigz などで作った 1 メンバーの gzip) ときは、そこから先を
// zlib-rs のストリームで解く。どちらも CRC と長さを確かめる。
// ---------------------------------------------------------------------------

const MEMBER_OUT_CAP: usize = 4 << 20;
const READ_AHEAD: usize = 4 << 20;
const READ_AHEAD_MAX: usize = 64 << 20;
const STREAM_CHUNK: usize = 1 << 20;

struct Decompressor(*mut ld::libdeflate_decompressor);
impl Drop for Decompressor {
    fn drop(&mut self) {
        unsafe { ld::libdeflate_free_decompressor(self.0) }
    }
}

fn spawn_reader(path: PathBuf) -> (Receiver<Res<Vec<u8>>>, JoinHandle<()>) {
    let (tx, rx) = sync_channel::<Res<Vec<u8>>>(32);
    let h = thread::spawn(move || {
        if let Err(e) = read_gzip(&path, &tx) {
            let _ = tx.send(Err(format!("{}: {}", path.display(), e)));
        }
    });
    (rx, h)
}

fn fill(f: &mut File, buf: &mut Vec<u8>, want: usize, eof: &mut bool) -> Res<()> {
    while buf.len() < want && !*eof {
        let old = buf.len();
        buf.resize(want, 0);
        match f.read(&mut buf[old..]) {
            Ok(0) => {
                buf.truncate(old);
                *eof = true;
            }
            Ok(n) => buf.truncate(old + n),
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => buf.truncate(old),
            Err(e) => return Err(format!("read error: {e}")),
        }
    }
    Ok(())
}

fn read_gzip(path: &Path, tx: &SyncSender<Res<Vec<u8>>>) -> Res<()> {
    let mut f = File::open(path).map_err(|e| format!("cannot open: {e}"))?;
    let d = Decompressor(unsafe { ld::libdeflate_alloc_decompressor() });
    if d.0.is_null() {
        return Err("libdeflate_alloc_decompressor failed".into());
    }
    let mut buf: Vec<u8> = Vec::with_capacity(READ_AHEAD * 2);
    let mut start = 0usize; // buf のうち解き終えたところ
    let mut file_off = 0u64; // buf[0] のファイル上の位置
    let mut eof = false;
    let mut want = READ_AHEAD;
    loop {
        if buf.len() - start < want && !eof {
            buf.drain(..start);
            file_off += start as u64;
            start = 0;
            fill(&mut f, &mut buf, want, &mut eof)?;
        }
        if start == buf.len() {
            return Ok(()); // 終わり
        }
        let mut out: Vec<u8> = Vec::with_capacity(MEMBER_OUT_CAP);
        let (mut used_in, mut used_out) = (0usize, 0usize);
        let r = unsafe {
            ld::libdeflate_gzip_decompress_ex(
                d.0,
                buf.as_ptr().add(start) as *const c_void,
                buf.len() - start,
                out.as_mut_ptr() as *mut c_void,
                MEMBER_OUT_CAP,
                &mut used_in,
                &mut used_out,
            )
        };
        if r == ld::libdeflate_result_LIBDEFLATE_SUCCESS {
            unsafe { out.set_len(used_out) };
            start += used_in;
            want = READ_AHEAD;
            if used_out > 0 && tx.send(Ok(out)).is_err() {
                return Ok(()); // 受け手が止まった
            }
            continue;
        }
        if r != ld::libdeflate_result_LIBDEFLATE_INSUFFICIENT_SPACE && !eof && want < READ_AHEAD_MAX {
            want *= 4; // メンバーの終わりがまだ読めていないかもしれない
            continue;
        }
        // 大きいメンバー、または壊れている: ここから先はストリームで解く (壊れていればエラーになる)。
        drop(buf);
        return stream_rest(f, file_off + start as u64, tx);
    }
}

fn stream_rest(mut f: File, off: u64, tx: &SyncSender<Res<Vec<u8>>>) -> Res<()> {
    f.seek(SeekFrom::Start(off)).map_err(|e| format!("seek error: {e}"))?;
    let mut dec = flate2::read::MultiGzDecoder::new(BufReader::with_capacity(1 << 20, f));
    loop {
        let mut out = vec![0u8; STREAM_CHUNK];
        let mut n = 0;
        while n < STREAM_CHUNK {
            match dec.read(&mut out[n..]) {
                Ok(0) => break,
                Ok(k) => n += k,
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
                Err(e) => return Err(format!("gzip error (corrupt or truncated?): {e}")),
            }
        }
        if n == 0 {
            return Ok(());
        }
        out.truncate(n);
        if tx.send(Ok(out)).is_err() {
            return Ok(());
        }
    }
}

// ---------------------------------------------------------------------------
// 行を 4 行 (1 リード) ずつ取り出す。
// ---------------------------------------------------------------------------

struct Lines {
    rx: Receiver<Res<Vec<u8>>>,
    buf: Vec<u8>,
    pos: usize,
    done: bool,
}

type Spans = [(usize, usize); 4];

impl Lines {
    fn new(rx: Receiver<Res<Vec<u8>>>) -> Self {
        Lines { rx, buf: Vec::with_capacity(8 << 20), pos: 0, done: false }
    }

    /// 次の 4 行の位置を spans に入れ、取れた行の数 (0〜4) を返す。
    /// 返した位置は、次に呼ぶまで有効。最後の行に改行が無くても 1 行と数える。
    fn next_record(&mut self, spans: &mut Spans) -> Res<usize> {
        let mut n = 0;
        let mut scan = self.pos;
        loop {
            while n < 4 {
                match memchr::memchr(b'\n', &self.buf[scan..]) {
                    Some(i) => {
                        spans[n] = (scan, scan + i);
                        scan += i + 1;
                        n += 1;
                    }
                    None => break,
                }
            }
            if n == 4 {
                self.pos = scan;
                return Ok(4);
            }
            if !self.done {
                // 足りない: 前を詰めて次の塊を足す
                let shift = self.pos;
                if shift > 0 {
                    self.buf.drain(..shift);
                    for s in spans.iter_mut().take(n) {
                        s.0 -= shift;
                        s.1 -= shift;
                    }
                    scan -= shift;
                    self.pos = 0;
                }
                match self.rx.recv() {
                    Ok(Ok(chunk)) => {
                        self.buf.extend_from_slice(&chunk);
                        continue;
                    }
                    Ok(Err(e)) => return Err(e),
                    Err(_) => self.done = true,
                }
            }
            // 入力の終わり: 改行の無い最後の行も 1 行と数える
            if scan < self.buf.len() {
                spans[n] = (scan, self.buf.len());
                scan = self.buf.len();
                n += 1;
                continue;
            }
            self.pos = scan;
            return Ok(n);
        }
    }
}

// ---------------------------------------------------------------------------
// awk の既定の欄の分け方 (space と TAB の並びで分け、前後の空白は無視) を、
// R1 の行と R2 の行を TAB でつないだ行に当てたときの n 番目 (0 始まり) の欄。
// ---------------------------------------------------------------------------

#[inline]
fn is_ws(c: u8) -> bool {
    c == b' ' || c == b'\t'
}

#[inline]
fn nth_field<'a>(a: &'a [u8], b: &'a [u8], mut n: usize) -> &'a [u8] {
    for line in [a, b] {
        let len = line.len();
        let mut i = 0;
        loop {
            while i < len && is_ws(line[i]) {
                i += 1;
            }
            if i >= len {
                break;
            }
            let st = i;
            while i < len && !is_ws(line[i]) {
                i += 1;
            }
            if n == 0 {
                return &line[st..i];
            }
            n -= 1;
        }
    }
    &[]
}

// ---------------------------------------------------------------------------
// 圧縮 (スレッドを並べる) と書き出し (1 スレッド、サンプルごとに順を守る)。
// ---------------------------------------------------------------------------

struct Job {
    sample: usize,
    seq: u64,
    data: Vec<u8>,
}

struct Compressor(*mut ld::libdeflate_compressor);
impl Drop for Compressor {
    fn drop(&mut self) {
        unsafe { ld::libdeflate_free_compressor(self.0) }
    }
}

fn spawn_compressor(level: i32, rx: Arc<Mutex<Receiver<Job>>>, tx: SyncSender<Job>) -> JoinHandle<Res<()>> {
    thread::spawn(move || {
        let c = Compressor(unsafe { ld::libdeflate_alloc_compressor(level) });
        if c.0.is_null() {
            return Err("libdeflate_alloc_compressor failed".into());
        }
        loop {
            let job = match rx.lock().unwrap().recv() {
                Ok(j) => j,
                Err(_) => return Ok(()),
            };
            let bound = unsafe { ld::libdeflate_gzip_compress_bound(c.0, job.data.len()) };
            let mut out: Vec<u8> = Vec::with_capacity(bound);
            let n = unsafe {
                ld::libdeflate_gzip_compress(
                    c.0,
                    job.data.as_ptr() as *const c_void,
                    job.data.len(),
                    out.as_mut_ptr() as *mut c_void,
                    bound,
                )
            };
            if n == 0 {
                return Err("libdeflate_gzip_compress failed".into());
            }
            unsafe { out.set_len(n) };
            if tx.send(Job { sample: job.sample, seq: job.seq, data: out }).is_err() {
                return Ok(()); // 書き手が止まった (エラーは書き手が返す)
            }
        }
    })
}

fn spawn_writer(paths: Vec<PathBuf>, rx: Receiver<Job>) -> JoinHandle<Res<()>> {
    thread::spawn(move || {
        let ns = paths.len();
        let mut files: Vec<Option<File>> = (0..ns).map(|_| None).collect();
        let mut next = vec![0u64; ns];
        let mut pending: Vec<BTreeMap<u64, Vec<u8>>> = (0..ns).map(|_| BTreeMap::new()).collect();
        for job in rx {
            let s = job.sample;
            pending[s].insert(job.seq, job.data);
            while let Some(data) = pending[s].remove(&next[s]) {
                if files[s].is_none() {
                    let f = File::create(&paths[s])
                        .map_err(|e| format!("cannot create {}: {e}", paths[s].display()))?;
                    files[s] = Some(f);
                }
                files[s]
                    .as_mut()
                    .unwrap()
                    .write_all(&data)
                    .map_err(|e| format!("write error {}: {e}", paths[s].display()))?;
                next[s] += 1;
            }
        }
        if pending.iter().any(|p| !p.is_empty()) {
            return Err("internal error: compressed blocks left unwritten".into());
        }
        Ok(())
    })
}

// ---------------------------------------------------------------------------

fn run(args: Args) -> Res<()> {
    // サンプル (最初に出てきた順) と barcode -> サンプル
    let mut sample_names: Vec<String> = Vec::new();
    let mut sample_index: HashMap<String, usize> = HashMap::new();
    let mut bc_to_sample: HashMap<Vec<u8>, usize> = HashMap::new();
    for (bc, s) in args.barcodes.iter().zip(args.samples.iter()) {
        let idx = *sample_index.entry(s.clone()).or_insert_with(|| {
            sample_names.push(s.clone());
            sample_names.len() - 1
        });
        bc_to_sample.insert(bc.as_bytes().to_vec(), idx);
    }
    let ns = sample_names.len();
    let paths: Vec<PathBuf> = sample_names
        .iter()
        .map(|s| args.out_dir.join(format!("{s}{}.tmp", args.suffix)))
        .collect();

    let (rx1, h1) = spawn_reader(args.r1.clone());
    let (rx2, h2) = spawn_reader(args.r2.clone());
    let mut l1 = Lines::new(rx1);
    let mut l2 = Lines::new(rx2);

    // --threads は 1 本で使うコアの目安。解凍 (R1・R2) と突き合わせで約 2 コアを使うので、
    // 圧縮のスレッドはその残り (1 以上)。
    let nworkers = args.threads.saturating_sub(2).max(1);
    let (job_tx, job_rx) = sync_channel::<Job>(nworkers * 2);
    let job_rx = Arc::new(Mutex::new(job_rx));
    let (done_tx, done_rx) = sync_channel::<Job>(nworkers * 4);
    let workers: Vec<_> = (0..nworkers)
        .map(|_| spawn_compressor(args.level, job_rx.clone(), done_tx.clone()))
        .collect();
    drop(done_tx);
    let writer = spawn_writer(paths, done_rx);

    let block = args.block;
    let mut bufs: Vec<Vec<u8>> = (0..ns).map(|_| Vec::new()).collect();
    let mut seqs = vec![0u64; ns];
    let mut nreads = vec![0u64; ns];
    let mut s1: Spans = [(0, 0); 4];
    let mut s2: Spans = [(0, 0); 4];
    let mut nrec: u64 = 0;
    let mut send_failed = false;

    let main_result: Res<()> = (|| loop {
        let n1 = l1.next_record(&mut s1)?;
        let n2 = l2.next_record(&mut s2)?;
        if n1 == 0 && n2 == 0 {
            return Ok(());
        }
        let line_no = nrec * 4;
        if n1 != n2 {
            let count = |n: usize| if n < 4 { format!("{}", line_no + n as u64) } else { format!("at least {}", line_no + 4) };
            return Err(format!(
                "R1 and R2 have different numbers of lines (R1: {} lines, R2: {} lines)",
                count(n1),
                count(n2)
            ));
        }
        if n1 < 4 {
            return Err(format!("truncated FASTQ input ({} lines, not a multiple of 4)", line_no + n1 as u64));
        }
        nrec += 1;
        let a = |i: usize| &l1.buf[s1[i].0..s1[i].1];
        let b = |i: usize| &l2.buf[s2[i].0..s2[i].1];
        let (a0, a1, a2, a3) = (a(0), a(1), a(2), a(3));
        let (b0, b1, b2, b3) = (b(0), b(1), b(2), b(3));
        if a0.first() != Some(&b'@') || b0.first() != Some(&b'@') || a2.first() != Some(&b'+') || b2.first() != Some(&b'+') {
            return Err(format!(
                "not a FASTQ record at read {nrec} (line {}): the 1st line must start with @ and the 3rd with +",
                line_no + 1
            ));
        }
        if b1.len() != b3.len() {
            return Err(format!("R2 sequence and quality lengths differ at read {nrec} (line {})", line_no + 1));
        }
        let f1 = nth_field(a1, b1, 0);
        let k6 = f1.len().min(6);
        let k12 = f1.len().min(12);
        let (bc, umi) = if args.barcode_first { (&f1[..k6], &f1[k6..k12]) } else { (&f1[k6..k12], &f1[..k6]) };
        let Some(&s) = bc_to_sample.get(bc) else { continue };
        nreads[s] += 1;
        let out = &mut bufs[s];
        if out.capacity() == 0 {
            out.reserve(block + 4096);
        }
        out.extend_from_slice(nth_field(a0, b0, 2));
        out.push(b'_');
        out.extend_from_slice(bc);
        out.push(b'_');
        out.extend_from_slice(umi);
        out.push(b' ');
        out.extend_from_slice(nth_field(a0, b0, 3));
        out.push(b'\n');
        out.extend_from_slice(nth_field(a1, b1, 1));
        out.push(b'\n');
        out.extend_from_slice(nth_field(a2, b2, 1));
        out.push(b'\n');
        out.extend_from_slice(nth_field(a3, b3, 1));
        out.push(b'\n');
        if out.len() >= block {
            let data = std::mem::replace(out, Vec::with_capacity(block + 4096));
            if job_tx.send(Job { sample: s, seq: seqs[s], data }).is_err() {
                send_failed = true;
                return Ok(());
            }
            seqs[s] += 1;
        }
    })();

    if main_result.is_ok() && !send_failed {
        for s in 0..ns {
            if !bufs[s].is_empty() {
                let data = std::mem::take(&mut bufs[s]);
                if job_tx.send(Job { sample: s, seq: seqs[s], data }).is_err() {
                    send_failed = true;
                    break;
                }
                seqs[s] += 1;
            }
        }
    }
    drop(job_tx);
    // 読み手を止める (エラーのときは受け手が無くなって止まる)
    drop(l1);
    drop(l2);
    let _ = h1.join();
    let _ = h2.join();
    let mut errs: Vec<String> = Vec::new();
    match writer.join() {
        Ok(Ok(())) => {}
        Ok(Err(e)) => errs.push(e),
        Err(_) => errs.push("writer thread panicked".into()),
    }
    for w in workers {
        match w.join() {
            Ok(Ok(())) => {}
            Ok(Err(e)) => errs.push(e),
            Err(_) => errs.push("compressor thread panicked".into()),
        }
    }
    main_result?;
    if let Some(e) = errs.into_iter().next() {
        return Err(e);
    }
    if send_failed {
        return Err("output pipeline stopped".into());
    }

    let mut text = String::new();
    for (i, name) in sample_names.iter().enumerate() {
        text.push_str(&format!("{name}\t{}\n", nreads[i]));
    }
    std::fs::write(&args.nread_file, text).map_err(|e| format!("cannot write {}: {e}", args.nread_file.display()))?;
    Ok(())
}

fn main() {
    let args = match parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("[ERROR] pic-demux: {e}");
            exit(2);
        }
    };
    if let Err(e) = run(args) {
        eprintln!("[ERROR] pic-demux: {e}");
        exit(2);
    }
}
