//! pic-aggregate: pic aggregate のメタジーン集計 (サンプルごとの全遺伝子平均プロファイル)。
//!
//! 旧来の道 (scripts/r/cmd_aggregate.R) は、ビンの BED に対してサンプルごとに UCSC
//! bigWigAverageOverBed (bwaob, v482) を回し、その出力 (%g の文字列) を R の data.table::fread
//! で読み、GForce の mean で集計していた。ここではそれを 1 本で行い、同じ値 (ビット単位) を返す:
//!
//!  1. ビンごとの和 (bwaob の averageFetchingEachChrom / addBufIntervalInfo と同じ):
//!     double sum1 = 0; for (i=start; i<end; ++i) sum1 += valBuf[i];
//!     valBuf は染色体ごとに 0 で埋めた double の配列に bigWig の float 値を書いたもの。
//!     値の無いベースは +0.0 を足すだけ (和は変わらない) なので、値のある区間だけを左から
//!     1 ベースずつ double で足す。染色体全長の配列は持たない (区間の配列を二分探索)。
//!     bigWig に無い染色体のビンは 0。
//!  2. mean0 = sum/size を C の "%g" で書き (bwaob の fprintf と同じ snprintf)、それを fread と同じ
//!     やり方 (仮数を整数に貯め、long double で 10^k を掛け/割り、double に落とす) で読み戻す。
//!  3. 行の順 = bwaob の出力順: BED の順を、染色体名で cmpStringsWithEmbeddedNumbers (kent) を
//!     比較に使った安定ソートで並べ替えたもの。
//!  4. 遺伝子ごとの mean0 の和 > 0 の遺伝子の行だけで、globalpos ごとに、その行順で double で
//!     足し、件数で割る (data.table 1.17 の gmean と同じ)。行が無い globalpos は NA。
//!
//! 使い方:
//!   pic-aggregate --bins <bins.tsv> --bigwigs <list.tsv> --out <out.tsv> [--threads N]
//!                 [--total N] [--dump-mean0 <dir>]
//!   bins.tsv : chr, start0, end, gene_index, globalpos (タブ区切り、見出し無し、BED と同じ行の順)
//!   list.tsv : sample, bigwig のパス (タブ区切り、見出し無し)
//!   out.tsv  : 見出し "sample globalpos value mant exp2"。value は %.17g、value = mant * 2^exp2
//!              (R はこちらで読むと丸め無しで同じ double になる)。値の無い globalpos は NA。
//!   --threads: 同時に処理するサンプルの数 (1 サンプルは 1 スレッド)。
//!   --dump-mean0: 確かめ用。<dir>/<sample>.tsv に bwaob の出力の 1 列目と 5 列目と同じもの
//!              ("<gene>|<globalpos>\t<mean0 の %g>") を bwaob と同じ行順で書く。
//! 終了コード: 0 = 済んだ (読めないサンプルは標準エラーに [WARN] を書いて飛ばす。out.tsv に行が無い)、
//!            1 = 入力・出力の誤り、2 = 引数の誤り。

use std::collections::HashMap;
use std::ffi::{c_char, c_int};
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

const VERSION: &str = env!("CARGO_PKG_VERSION");

type Res<T> = Result<T, String>;

// ---------------------------------------------------------------- C の printf

extern "C" {
    fn snprintf(buf: *mut c_char, n: usize, fmt: *const c_char, ...) -> c_int;
}

/// C の snprintf(buf, n, fmt, x) (fmt は NUL 終端の "%g" など)。
fn c_fmt<'a>(fmt: &[u8], x: f64, buf: &'a mut [u8; 64]) -> &'a [u8] {
    let n = unsafe { snprintf(buf.as_mut_ptr() as *mut c_char, buf.len(), fmt.as_ptr() as *const c_char, x) };
    let n = (n.max(0) as usize).min(buf.len() - 1);
    &buf[..n]
}

// ---------------------------------------------------------------- fread の double の読み

/// data.table (1.17) fread の parse_double_regular_core と同じ読み。
/// acc (18 桁まで) と 10 の指数 e を作り、r = (long double)acc; r = e<0 ? r/10^-e : r*10^e; (double)r。
struct FreadParser {
    #[cfg(target_arch = "x86_64")]
    p10: Vec<[u8; 10]>, // 10^0..10^300 を long double (x87 80 bit) にしたもの ("1.0E<k>L" と同じ丸め)
    #[cfg(not(target_arch = "x86_64"))]
    p10: Vec<f64>,
}

impl FreadParser {
    fn new() -> Self {
        #[cfg(target_arch = "x86_64")]
        {
            let p = FreadParser { p10: (0..=300).map(pow10_f80).collect() };
            p.self_test();
            p
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            // macOS arm64 では long double = double。そのほか (Linux aarch64 の 128 bit) では近似になる。
            FreadParser { p10: (0..=300).map(|k| format!("1e{}", k).parse::<f64>().unwrap()).collect() }
        }
    }

    #[cfg(target_arch = "x86_64")]
    fn self_test(&self) {
        // 割り算の向きと丸めの確かめ (x87 の命令の向きを取り違えていないか)。
        assert_eq!(x87::scale(1, &self.p10[1], None, true), 0.1);
        assert_eq!(x87::scale(123457, &self.p10[10], None, true), 1.23457e-05);
        assert_eq!(x87::scale(5, &self.p10[2], None, false), 500.0);
    }

    /// 文字列 s (fread が double として読む形) を読む。読めなければ NaN (fread の NA)。
    fn parse(&self, s: &[u8]) -> f64 {
        let mut i = 0usize;
        let n = s.len();
        let neg = n > 0 && s[0] == b'-';
        if n > 0 && (s[0] == b'-' || s[0] == b'+') {
            i += 1;
        }
        let mut acc: u64 = 0;
        let mut e: i32 = 0;
        let mut sflimit = 18i32;
        while i < n && s[i] == b'0' {
            i += 1;
        }
        while i < n && s[i].is_ascii_digit() && sflimit > 0 {
            acc = 10 * acc + (s[i] - b'0') as u64;
            sflimit -= 1;
            i += 1;
        }
        if sflimit == 0 && i < n && s[i].is_ascii_digit() {
            while i < n && s[i].is_ascii_digit() {
                i += 1;
                e += 1;
            }
        }
        if i < n && s[i] == b'.' {
            i += 1;
            if i < n && s[i] == b'0' && acc == 0 {
                let mut k = 0;
                while i + k < n && s[i + k] == b'0' {
                    k += 1;
                }
                i += k;
                e = -(k as i32);
            }
            let mut k = 0usize;
            while i + k < n && s[i + k].is_ascii_digit() && sflimit > 0 {
                acc = 10 * acc + (s[i + k] - b'0') as u64;
                k += 1;
                sflimit -= 1;
            }
            i += k;
            e -= k as i32;
            while i < n && s[i].is_ascii_digit() {
                i += 1; // 18 桁を超えた小数の桁は捨てる
            }
        }
        if i < n && (s[i] == b'e' || s[i] == b'E') {
            i += 1;
            let eneg = i < n && s[i] == b'-';
            if i < n && (s[i] == b'-' || s[i] == b'+') {
                i += 1;
            }
            let mut ee: i32 = 0;
            while i < n && s[i].is_ascii_digit() {
                ee = ee.saturating_mul(10).saturating_add((s[i] - b'0') as i32);
                i += 1;
            }
            e += if eneg { -ee } else { ee };
        }
        if i != n || !(-350..=350).contains(&e) {
            return f64::NAN;
        }
        let r = self.scale(acc, e);
        if neg { -r } else { r }
    }

    fn scale(&self, acc: u64, e: i32) -> f64 {
        if acc == 0 {
            return 0.0;
        }
        let (extra, e2) = if e < -300 || e > 300 {
            let extra = if e < 0 { e + 300 } else { e - 300 };
            (Some(extra), e - extra)
        } else {
            (None, e)
        };
        #[cfg(target_arch = "x86_64")]
        {
            let pe = &self.p10[e2.unsigned_abs() as usize];
            let px = extra.map(|x| &self.p10[x.unsigned_abs() as usize]);
            // extra と e2 は同じ符号 (どちらも割るか、どちらも掛ける)。
            x87::scale(acc, pe, px, e2 < 0)
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            let mut r = acc as f64;
            if let Some(x) = extra {
                let p = self.p10[x.unsigned_abs() as usize];
                r = if x < 0 { r / p } else { r * p };
            }
            let p = self.p10[e2.unsigned_abs() as usize];
            if e2 < 0 { r / p } else { r * p }
        }
    }
}

/// 10^k を x87 の 80 bit (64 bit の仮数、最も近い偶数への丸め) にする。
#[cfg(target_arch = "x86_64")]
fn pow10_f80(k: usize) -> [u8; 10] {
    // 10^k を 32 bit の桁の多倍長整数で作る (下位から)。
    let mut v: Vec<u32> = vec![1];
    for _ in 0..k {
        let mut carry: u64 = 0;
        for d in v.iter_mut() {
            let t = (*d as u64) * 10 + carry;
            *d = t as u32;
            carry = t >> 32;
        }
        if carry > 0 {
            v.push(carry as u32);
        }
    }
    let bit = |i: usize| -> u64 { ((v[i / 32] >> (i % 32)) & 1) as u64 };
    let nbits = 32 * (v.len() - 1) + (32 - v.last().unwrap().leading_zeros() as usize);
    let mut mant: u64 = 0;
    let take = nbits.min(64);
    for j in 0..take {
        mant = (mant << 1) | bit(nbits - 1 - j);
    }
    let mut exp = nbits as i64 - 1; // 値 = 1.xxx * 2^exp
    if nbits <= 64 {
        mant <<= 64 - nbits;
    } else {
        let rest = nbits - 64; // 捨てる下位の桁数
        let guard = bit(rest - 1);
        let sticky = (0..rest - 1).any(|i| bit(i) == 1);
        if guard == 1 && (sticky || (mant & 1) == 1) {
            let (m2, of) = mant.overflowing_add(1);
            if of {
                mant = 1u64 << 63;
                exp += 1;
            } else {
                mant = m2;
            }
        }
    }
    let mut out = [0u8; 10];
    out[..8].copy_from_slice(&mant.to_le_bytes());
    let se = (exp + 16383) as u16;
    out[8..].copy_from_slice(&se.to_le_bytes());
    out
}

#[cfg(target_arch = "x86_64")]
mod x87 {
    use std::arch::asm;
    /// (long double)acc を p1 (と p0) で割る (div) か掛け、double に丸めて返す。p0 があれば先に p0。
    /// x87 の精度は Linux の既定 (拡張 64 bit) のまま、C の long double の演算と同じ。
    pub fn scale(acc: u64, p1: &[u8; 10], p0: Option<&[u8; 10]>, div: bool) -> f64 {
        let a: i64 = acc as i64; // acc < 10^18 < 2^63
        let mut out: f64 = 0.0;
        let one = super::pow10_f80(0);
        let p0 = p0.unwrap_or(&one);
        unsafe {
            if div {
                asm!(
                    "fild qword ptr [{a}]",
                    "fld tbyte ptr [{p0}]",
                    "fdivp st(1), st",
                    "fld tbyte ptr [{p1}]",
                    "fdivp st(1), st",
                    "fstp qword ptr [{o}]",
                    a = in(reg) &a, p0 = in(reg) p0.as_ptr(), p1 = in(reg) p1.as_ptr(),
                    o = in(reg) &mut out as *mut f64,
                    out("st(0)") _, out("st(1)") _, out("st(2)") _, out("st(3)") _,
                    out("st(4)") _, out("st(5)") _, out("st(6)") _, out("st(7)") _,
                    options(nostack),
                );
            } else {
                asm!(
                    "fild qword ptr [{a}]",
                    "fld tbyte ptr [{p0}]",
                    "fmulp st(1), st",
                    "fld tbyte ptr [{p1}]",
                    "fmulp st(1), st",
                    "fstp qword ptr [{o}]",
                    a = in(reg) &a, p0 = in(reg) p0.as_ptr(), p1 = in(reg) p1.as_ptr(),
                    o = in(reg) &mut out as *mut f64,
                    out("st(0)") _, out("st(1)") _, out("st(2)") _, out("st(3)") _,
                    out("st(4)") _, out("st(5)") _, out("st(6)") _, out("st(7)") _,
                    options(nostack),
                );
            }
        }
        out
    }
}

// ---------------------------------------------------------------- 染色体名の比較 (kent)

/// glibc の atoi (= (int)strtol(s, NULL, 10))。
fn c_atoi(s: &[u8]) -> i32 {
    let mut i = 0;
    while i < s.len() && matches!(s[i], b' ' | b'\t' | b'\n' | b'\x0b' | b'\x0c' | b'\r') {
        i += 1;
    }
    let mut neg = false;
    if i < s.len() && (s[i] == b'+' || s[i] == b'-') {
        neg = s[i] == b'-';
        i += 1;
    }
    let mut v: i64 = 0;
    let mut over = false;
    while i < s.len() && s[i].is_ascii_digit() {
        let d = (s[i] - b'0') as i64;
        if !over {
            match v.checked_mul(10).and_then(|x| x.checked_add(d)) {
                Some(x) => v = x,
                None => over = true,
            }
        }
        i += 1;
    }
    let l: i64 = if over {
        if neg { i64::MIN } else { i64::MAX }
    } else if neg {
        -v
    } else {
        v
    };
    l as i32
}

/// kent の cmpStringsWithEmbeddedNumbers (bwaob の bedCmpChrom が使う)。
fn cmp_embedded(mut a: &[u8], mut b: &[u8]) -> std::cmp::Ordering {
    use std::cmp::Ordering::*;
    let lead = |s: &[u8], dig: bool| s.iter().take_while(|c| c.is_ascii_digit() == dig).count();
    loop {
        let an = lead(a, true);
        let bn = lead(b, true);
        let diff = c_atoi(a).wrapping_sub(c_atoi(b));
        if diff != 0 {
            return diff.cmp(&0);
        }
        a = &a[an..];
        b = &b[bn..];
        let ann = lead(a, false);
        let bnn = lead(b, false);
        if ann != bnn {
            return a.cmp(b); // strcmp
        } else if ann == 0 {
            return Equal;
        } else {
            let c = a[..ann].cmp(&b[..ann]);
            if c != Equal {
                return c;
            }
            a = &a[ann..];
            b = &b[bnn..];
        }
    }
}

// ---------------------------------------------------------------- ビンの表

struct Bins {
    chrom_names: Vec<String>,
    chrom: Vec<u32>,
    start: Vec<u32>,
    end: Vec<u32>,
    gene: Vec<u32>,
    gpos: Vec<u32>,
    name: Vec<String>,  // 確かめ用 ("<gene>|<globalpos>"、--dump-mean0 のときだけ作る)
    order: Vec<u32>,    // bwaob の出力の行順
    by_chrom: Vec<Vec<u32>>,
    ngenes: usize,
    total: usize,
}

fn parse_u32(f: &[u8], what: &str, line: usize) -> Res<u32> {
    if f.is_empty() {
        return Err(format!("bins line {}: empty {}", line, what));
    }
    let mut v: u64 = 0;
    for &c in f {
        if !c.is_ascii_digit() {
            return Err(format!("bins line {}: bad {} '{}'", line, what, String::from_utf8_lossy(f)));
        }
        v = v * 10 + (c - b'0') as u64;
        if v > u32::MAX as u64 {
            return Err(format!("bins line {}: {} too large", line, what));
        }
    }
    Ok(v as u32)
}

fn read_bins(path: &Path, want_names: bool) -> Res<Bins> {
    let mut data = Vec::new();
    File::open(path)
        .and_then(|mut f| f.read_to_end(&mut data))
        .map_err(|e| format!("cannot read bins {}: {}", path.display(), e))?;
    let mut b = Bins {
        chrom_names: vec![], chrom: vec![], start: vec![], end: vec![], gene: vec![], gpos: vec![],
        name: vec![], order: vec![], by_chrom: vec![], ngenes: 0, total: 0,
    };
    let mut ids: HashMap<Vec<u8>, u32> = HashMap::new();
    for (ln, line) in data.split(|&c| c == b'\n').enumerate() {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.is_empty() {
            continue;
        }
        let mut it = line.split(|&c| c == b'\t');
        let mut next = |w: &str| it.next().ok_or_else(|| format!("bins line {}: missing {}", ln + 1, w));
        let c = next("chr")?;
        let s = parse_u32(next("start0")?, "start0", ln + 1)?;
        let e = parse_u32(next("end")?, "end", ln + 1)?;
        let g = parse_u32(next("gene_index")?, "gene_index", ln + 1)?;
        let p = parse_u32(next("globalpos")?, "globalpos", ln + 1)?;
        if e <= s {
            return Err(format!("bins line {}: end <= start0", ln + 1));
        }
        if p == 0 {
            return Err(format!("bins line {}: globalpos must be >= 1", ln + 1));
        }
        let id = match ids.get(c) {
            Some(&id) => id,
            None => {
                let id = b.chrom_names.len() as u32;
                b.chrom_names.push(String::from_utf8_lossy(c).into_owned());
                ids.insert(c.to_vec(), id);
                id
            }
        };
        b.chrom.push(id);
        b.start.push(s);
        b.end.push(e);
        b.gene.push(g);
        b.gpos.push(p);
        if want_names {
            b.name.push(format!("{}|{}", g, p));
        }
        b.ngenes = b.ngenes.max(g as usize + 1);
        b.total = b.total.max(p as usize);
    }
    drop(data);
    // bwaob の行順: 染色体名で安定ソート (同じと比べられる名前は同じ順位)。
    let nc = b.chrom_names.len();
    let mut cidx: Vec<usize> = (0..nc).collect();
    cidx.sort_by(|&x, &y| cmp_embedded(b.chrom_names[x].as_bytes(), b.chrom_names[y].as_bytes()));
    let mut rank = vec![0usize; nc];
    let mut r = 0;
    for (k, &ci) in cidx.iter().enumerate() {
        if k > 0 && cmp_embedded(b.chrom_names[cidx[k - 1]].as_bytes(), b.chrom_names[ci].as_bytes()) != std::cmp::Ordering::Equal {
            r += 1;
        }
        rank[ci] = r;
    }
    let mut buckets: Vec<Vec<u32>> = vec![vec![]; r + 1];
    b.by_chrom = vec![vec![]; nc];
    for (i, &c) in b.chrom.iter().enumerate() {
        buckets[rank[c as usize]].push(i as u32);
        b.by_chrom[c as usize].push(i as u32);
    }
    b.order = buckets.concat();
    Ok(b)
}

// ---------------------------------------------------------------- bigWig

struct BigWig {
    data: Vec<u8>,
    chroms: HashMap<String, (u32, u32)>, // name -> (id, size)
    leaves: Vec<(u32, u32, u64, u64)>,   // (startChromIx, endChromIx, offset, size)
    compressed: bool,
}

fn rd_u16(d: &[u8], o: usize) -> Res<u16> {
    d.get(o..o + 2).map(|s| u16::from_le_bytes(s.try_into().unwrap())).ok_or_else(|| "truncated bigWig".to_string())
}
fn rd_u32(d: &[u8], o: usize) -> Res<u32> {
    d.get(o..o + 4).map(|s| u32::from_le_bytes(s.try_into().unwrap())).ok_or_else(|| "truncated bigWig".to_string())
}
fn rd_u64(d: &[u8], o: usize) -> Res<u64> {
    d.get(o..o + 8).map(|s| u64::from_le_bytes(s.try_into().unwrap())).ok_or_else(|| "truncated bigWig".to_string())
}

impl BigWig {
    fn open(path: &Path) -> Res<BigWig> {
        let mut data = Vec::new();
        File::open(path)
            .and_then(|mut f| f.read_to_end(&mut data))
            .map_err(|e| format!("cannot read {}: {}", path.display(), e))?;
        if rd_u32(&data, 0)? != 0x888F_FC26 {
            return Err(format!("not a (little-endian) bigWig: {}", path.display()));
        }
        let chrom_tree = rd_u64(&data, 8)? as usize;
        let index = rd_u64(&data, 24)? as usize;
        let compressed = rd_u32(&data, 52)? > 0;
        // 染色体の B+ 木
        if rd_u32(&data, chrom_tree)? != 0x78CA_8C91 {
            return Err("bad chromosome tree magic".into());
        }
        let key_size = rd_u32(&data, chrom_tree + 8)? as usize;
        let mut chroms = HashMap::new();
        let mut stack = vec![chrom_tree + 32];
        let mut guard = 0usize;
        while let Some(off) = stack.pop() {
            guard += 1;
            if guard > 10_000_000 {
                return Err("chromosome tree loop".into());
            }
            let is_leaf = *data.get(off).ok_or("truncated bigWig")? != 0;
            let count = rd_u16(&data, off + 2)? as usize;
            let mut p = off + 4;
            for _ in 0..count {
                let key = data.get(p..p + key_size).ok_or("truncated bigWig")?;
                let klen = key.iter().position(|&c| c == 0).unwrap_or(key_size);
                if is_leaf {
                    let id = rd_u32(&data, p + key_size)?;
                    let size = rd_u32(&data, p + key_size + 4)?;
                    chroms.insert(String::from_utf8_lossy(&key[..klen]).into_owned(), (id, size));
                    p += key_size + 8;
                } else {
                    stack.push(rd_u64(&data, p + key_size)? as usize);
                    p += key_size + 8;
                }
            }
        }
        // R 木 (unzoomed)
        if rd_u32(&data, index)? != 0x2468_ACE0 {
            return Err("bad R-tree magic".into());
        }
        let mut leaves = vec![];
        let mut stack = vec![index + 48];
        guard = 0;
        while let Some(off) = stack.pop() {
            guard += 1;
            if guard > 100_000_000 {
                return Err("R-tree loop".into());
            }
            let is_leaf = *data.get(off).ok_or("truncated bigWig")? != 0;
            let count = rd_u16(&data, off + 2)? as usize;
            let mut p = off + 4;
            for _ in 0..count {
                let sc = rd_u32(&data, p)?;
                let ec = rd_u32(&data, p + 8)?;
                if is_leaf {
                    leaves.push((sc, ec, rd_u64(&data, p + 16)?, rd_u64(&data, p + 24)?));
                    p += 32;
                } else {
                    stack.push(rd_u64(&data, p + 16)? as usize);
                    p += 24;
                }
            }
        }
        leaves.sort_by_key(|l| l.2); // ファイルの順 (bwaob が値を書く順)
        Ok(BigWig { data, chroms, leaves, compressed })
    }

    /// 染色体 id の区間 (start, end, val) をファイルの順で。[0, size) に切る。
    fn intervals(&self, id: u32, size: u32) -> Res<Vec<(u32, u32, f32)>> {
        let mut out = vec![];
        let mut buf = Vec::new();
        for &(sc, ec, off, sz) in &self.leaves {
            if id < sc || id > ec {
                continue;
            }
            let raw = self.data.get(off as usize..(off + sz) as usize).ok_or("truncated bigWig block")?;
            let blk: &[u8] = if self.compressed {
                buf.clear();
                flate2::read::ZlibDecoder::new(raw).read_to_end(&mut buf).map_err(|e| format!("bad bigWig block: {}", e))?;
                &buf
            } else {
                raw
            };
            let cid = rd_u32(blk, 0)?;
            if cid != id {
                continue;
            }
            let sec_start = rd_u32(blk, 4)?;
            let step = rd_u32(blk, 12)?;
            let span = rd_u32(blk, 16)?;
            let typ = *blk.get(20).ok_or("truncated bigWig block")?;
            let n = rd_u16(blk, 22)? as usize;
            let mut push = |s: u32, e: u32, v: f32| {
                let s = s.min(size);
                let e = e.min(size);
                if s < e {
                    out.push((s, e, v));
                }
            };
            let fv = |o: usize| -> Res<f32> { Ok(f32::from_bits(rd_u32(blk, o)?)) };
            match typ {
                1 => {
                    for i in 0..n {
                        let o = 24 + 12 * i;
                        push(rd_u32(blk, o)?, rd_u32(blk, o + 4)?, fv(o + 8)?);
                    }
                }
                2 => {
                    for i in 0..n {
                        let o = 24 + 8 * i;
                        let s = rd_u32(blk, o)?;
                        push(s, s.saturating_add(span), fv(o + 4)?);
                    }
                }
                3 => {
                    for i in 0..n {
                        let s = sec_start.saturating_add((i as u32).saturating_mul(step));
                        push(s, s.saturating_add(span), fv(24 + 4 * i)?);
                    }
                }
                t => return Err(format!("unknown bigWig section type {}", t)),
            }
        }
        Ok(out)
    }
}

/// 1 染色体の値。普通は並んだ重ならない区間。重なり・逆順があれば bwaob と同じく後の値で上書きした
/// 密な配列にする (bedGraphToBigWig の出力では起きない)。
enum ChromVals {
    Sparse(Vec<(u32, u32, f32)>),
    Dense(Vec<f32>),
}

impl ChromVals {
    fn new(iv: Vec<(u32, u32, f32)>, size: u32) -> ChromVals {
        let sorted = iv.windows(2).all(|w| w[1].0 >= w[0].1);
        if sorted {
            ChromVals::Sparse(iv)
        } else {
            let mut d = vec![0f32; size as usize];
            for &(s, e, v) in &iv {
                for x in &mut d[s as usize..e as usize] {
                    *x = v;
                }
            }
            ChromVals::Dense(d)
        }
    }

    /// bwaob の addBufIntervalInfo の和: [s, e) を左から 1 ベースずつ double で足す
    /// (値の無いベースの +0.0 は和を変えないので飛ばす)。
    fn sum(&self, s: u32, e: u32) -> f64 {
        let mut sum = 0.0f64;
        match self {
            ChromVals::Sparse(iv) => {
                let mut k = iv.partition_point(|x| x.1 <= s);
                while k < iv.len() && iv[k].0 < e {
                    let (a, b, v) = iv[k];
                    let lo = a.max(s);
                    let hi = b.min(e);
                    let v = v as f64;
                    for _ in lo..hi {
                        sum += v;
                    }
                    k += 1;
                }
            }
            ChromVals::Dense(d) => {
                let hi = (e as usize).min(d.len());
                for &x in d.get(s as usize..hi).unwrap_or(&[]) {
                    sum += x as f64;
                }
            }
        }
        sum
    }
}

// ---------------------------------------------------------------- 1 サンプル

fn sample_profile(bins: &Bins, fp: &FreadParser, bw_path: &Path, dump: Option<&Path>) -> Res<Vec<Option<f64>>> {
    let bw = BigWig::open(bw_path)?;
    let n = bins.chrom.len();
    let mut mean0 = vec![0f64; n];
    let mut gbuf = [0u8; 64];
    let mut text: Vec<Vec<u8>> = if dump.is_some() { vec![b"0".to_vec(); n] } else { vec![] };
    for (ci, idx) in bins.by_chrom.iter().enumerate() {
        let Some(&(id, size)) = bw.chroms.get(&bins.chrom_names[ci]) else { continue };
        if size == 0 {
            continue;
        }
        let vals = ChromVals::new(bw.intervals(id, size)?, size);
        for &b in idx {
            let b = b as usize;
            let (s, e) = (bins.start[b], bins.end[b]);
            let size_b = (e - s) as i32 as f64; // bwaob: int size; sum/size
            let sum = vals.sum(s, e.min(size).max(s));
            let g = c_fmt(b"%g\0", sum / size_b, &mut gbuf);
            mean0[b] = fp.parse(g);
            if dump.is_some() {
                text[b] = g.to_vec();
            }
        }
    }
    if let Some(p) = dump {
        let mut w = BufWriter::new(File::create(p).map_err(|e| format!("{}: {}", p.display(), e))?);
        for &b in &bins.order {
            let b = b as usize;
            w.write_all(bins.name[b].as_bytes()).and_then(|_| w.write_all(b"\t")).and_then(|_| w.write_all(&text[b])).and_then(|_| w.write_all(b"\n"))
                .map_err(|e| e.to_string())?;
        }
        w.flush().map_err(|e| e.to_string())?;
    }
    // 遺伝子ごとの mean0 の和 (R: gtot := sum(mean0), by = gi)。値は 0 以上なので > 0 は「どれかが > 0」。
    let mut gtot = vec![0f64; bins.ngenes];
    for &b in &bins.order {
        gtot[bins.gene[b as usize] as usize] += mean0[b as usize];
    }
    // R: tab[gtot > 0, .(value = mean(mean0)), by = binpos] (GForce の gmean: 行順に double で足して件数で割る)
    let mut acc = vec![0f64; bins.total + 1];
    let mut cnt = vec![0u64; bins.total + 1];
    for &b in &bins.order {
        let b = b as usize;
        if gtot[bins.gene[b] as usize] > 0.0 {
            let p = bins.gpos[b] as usize;
            acc[p] += mean0[b];
            cnt[p] += 1;
        }
    }
    Ok((1..=bins.total).map(|p| if cnt[p] > 0 { Some(acc[p] / cnt[p] as f64) } else { None }).collect())
}

/// x = mant * 2^exp2 (mant は奇数か 0、|mant| < 2^53)。
fn decompose(x: f64) -> (i64, i32) {
    if x == 0.0 {
        return (0, 0);
    }
    let bits = x.to_bits();
    let ex = ((bits >> 52) & 0x7ff) as i32;
    let fr = (bits & ((1u64 << 52) - 1)) as i64;
    let (mut m, mut e) = if ex == 0 { (fr, -1074) } else { (fr | (1i64 << 52), ex - 1075) };
    let tz = m.trailing_zeros() as i32;
    m >>= tz;
    e += tz;
    if x < 0.0 { (-m, e) } else { (m, e) }
}

// ---------------------------------------------------------------- main

fn usage() -> ! {
    eprintln!(
        "usage: pic-aggregate --bins <bins.tsv> --bigwigs <list.tsv> --out <out.tsv> [--threads N] [--total N] [--dump-mean0 <dir>]\n\
         \x20 bins.tsv: chr start0 end gene_index globalpos (tab, no header; same rows/order as the BED)\n\
         \x20 list.tsv: sample bigwig_path (tab, no header)\n\
         \x20 out.tsv : sample globalpos value mant exp2 (value = mant * 2^exp2; NA if no row)"
    );
    std::process::exit(2)
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (mut bins_p, mut list_p, mut out_p, mut dump_p) = (None, None, None, None);
    let mut threads = 1usize;
    let mut total_arg: Option<usize> = None;
    let mut i = 0;
    while i < args.len() {
        let k = args[i].as_str();
        let val = |i: usize| args.get(i + 1).cloned().unwrap_or_else(|| usage());
        match k {
            "--version" | "-V" => {
                println!("pic-aggregate {}", VERSION);
                return;
            }
            "--help" | "-h" => usage(),
            "--bins" => { bins_p = Some(PathBuf::from(val(i))); i += 2; }
            "--bigwigs" => { list_p = Some(PathBuf::from(val(i))); i += 2; }
            "--out" => { out_p = Some(PathBuf::from(val(i))); i += 2; }
            "--dump-mean0" => { dump_p = Some(PathBuf::from(val(i))); i += 2; }
            "--threads" => { threads = val(i).parse().unwrap_or_else(|_| usage()); i += 2; }
            "--total" => { total_arg = Some(val(i).parse().unwrap_or_else(|_| usage())); i += 2; }
            _ => { eprintln!("pic-aggregate: unknown option {}", k); usage() }
        }
    }
    let (Some(bins_p), Some(list_p), Some(out_p)) = (bins_p, list_p, out_p) else { usage() };
    if let Err(e) = run(&bins_p, &list_p, &out_p, dump_p.as_deref(), threads.max(1), total_arg) {
        eprintln!("[ERROR] pic-aggregate: {}", e);
        std::process::exit(1);
    }
}

fn run(bins_p: &Path, list_p: &Path, out_p: &Path, dump: Option<&Path>, threads: usize, total_arg: Option<usize>) -> Res<()> {
    let t0 = Instant::now();
    let list = std::fs::read_to_string(list_p).map_err(|e| format!("cannot read {}: {}", list_p.display(), e))?;
    let samples: Vec<(String, PathBuf)> = list
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let mut it = l.splitn(2, '\t');
            let s = it.next().unwrap().to_string();
            let p = it.next().ok_or_else(|| format!("bigwig list: no path for '{}'", s))?;
            Ok((s, PathBuf::from(p)))
        })
        .collect::<Res<_>>()?;
    if let Some(d) = dump {
        std::fs::create_dir_all(d).map_err(|e| format!("{}: {}", d.display(), e))?;
    }
    let mut bins = read_bins(bins_p, dump.is_some())?;
    if let Some(t) = total_arg {
        if t < bins.total {
            return Err(format!("--total {} < max globalpos {}", t, bins.total));
        }
        bins.total = t;
    }
    let fp = FreadParser::new();
    eprintln!("[INFO] pic-aggregate: {} bins, {} genes, {} samples, threads {} (read {:.1}s)",
              bins.chrom.len(), bins.ngenes, samples.len(), threads, t0.elapsed().as_secs_f64());

    let next = AtomicUsize::new(0);
    let results: Mutex<Vec<Option<Vec<Option<f64>>>>> = Mutex::new(vec![None; samples.len()]);
    std::thread::scope(|sc| {
        for _ in 0..threads.min(samples.len().max(1)) {
            sc.spawn(|| loop {
                let k = next.fetch_add(1, Ordering::SeqCst);
                if k >= samples.len() {
                    break;
                }
                let (name, path) = &samples[k];
                let t = Instant::now();
                let dpath = dump.map(|d| d.join(format!("{}.tsv", name)));
                match sample_profile(&bins, &fp, path, dpath.as_deref()) {
                    Ok(v) => {
                        eprintln!("[INFO] pic-aggregate: {} ({:.1}s)", name, t.elapsed().as_secs_f64());
                        results.lock().unwrap()[k] = Some(v);
                    }
                    Err(e) => eprintln!("[WARN] pic-aggregate: {}: {}", name, e),
                }
            });
        }
    });
    let results = results.into_inner().unwrap();

    let tmp = out_p.with_extension("tmp");
    let write = || -> std::io::Result<()> {
        let mut w = BufWriter::new(File::create(&tmp)?);
        writeln!(w, "sample\tglobalpos\tvalue\tmant\texp2")?;
        let mut buf = [0u8; 64];
        for (k, r) in results.iter().enumerate() {
            let Some(v) = r else { continue };
            for (p, x) in v.iter().enumerate() {
                match x {
                    Some(x) => {
                        let (m, e) = decompose(*x);
                        w.write_all(samples[k].0.as_bytes())?;
                        write!(w, "\t{}\t", p + 1)?;
                        w.write_all(c_fmt(b"%.17g\0", *x, &mut buf))?;
                        writeln!(w, "\t{}\t{}", m, e)?;
                    }
                    None => writeln!(w, "{}\t{}\tNA\tNA\tNA", samples[k].0, p + 1)?,
                }
            }
        }
        w.flush()
    };
    write().map_err(|e| format!("cannot write {}: {}", tmp.display(), e))?;
    std::fs::rename(&tmp, out_p).map_err(|e| format!("cannot rename to {}: {}", out_p.display(), e))?;
    let ok = results.iter().filter(|r| r.is_some()).count();
    eprintln!("[INFO] pic-aggregate: {}/{} samples, {:.1}s", ok, samples.len(), t0.elapsed().as_secs_f64());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn embedded() {
        use std::cmp::Ordering::*;
        assert_eq!(cmp_embedded(b"GL456210.1", b"1"), Less);
        assert_eq!(cmp_embedded(b"X", b"1"), Less);
        assert_eq!(cmp_embedded(b"MT", b"X"), Less);
        assert_eq!(cmp_embedded(b"9", b"10"), Less);
        assert_eq!(cmp_embedded(b"chr4", b"chr14"), Less);
        assert_eq!(cmp_embedded(b"chr1", b"chr1"), Equal);
    }
    #[test]
    fn fread_parse() {
        let fp = FreadParser::new();
        assert_eq!(fp.parse(b"0"), 0.0);
        assert_eq!(fp.parse(b"0.5"), 0.5);
        assert_eq!(fp.parse(b"1.23457e-05"), 1.23457e-05);
        assert_eq!(fp.parse(b"12345.6"), 12345.6);
    }
    #[test]
    fn decomp() {
        for x in [0.1f64, 3.0, 1e-300, 123456.789] {
            let (m, e) = decompose(x);
            assert_eq!((m as f64) * 2f64.powi(e), x);
        }
    }
}
