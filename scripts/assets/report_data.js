// 埋め込みデータの読み出し (report.js より前に読み込む)。
// R (report_json.R の pic_z_*) が zlib (RFC1950) で縮めて base64 にした塊を、実行されない
// <script type="application/x-pic-z" id="..."> に入れている。ここで標準の DecompressionStream("deflate")
// で展開する。外部のライブラリは使わない。
//   pic-z-plotly         plotly.min.js           → 開いたときに展開して実行する
//   pic-z-PIC_PLOTS など  図・発現・PCA の JSON    → 開いたときに展開して window.PIC_PLOTS などに入れる
//   pic-z-frag-<k>       大きな表の HTML          → 開いたときに展開して <template class="pic-zfrag"> の場所に戻す
//   pic-f-<id>           ダウンロード用のファイル → ダウンロードのときに展開する (data-n がファイル名)
// 公開する関数:
//   picData(name)    → Promise (JSON の値。name は "PIC_PLOTS" など。1 回だけ展開して覚える)
//   picFileBlob(id)  → Promise ({name, blob}。blob は展開済みの元のバイト列)
//   picReady(fn)     → DOMContentLoaded と上の展開がすべて終わってから fn を呼ぶ
(function () {
  var HAS_DS = typeof DecompressionStream !== "undefined";
  var DATA_NAMES = ["PIC_PLOTS", "PIC_EXPR", "PIC_PCA"];

  function block(id) { return document.getElementById(id); }

  // base64 → バイト列 (Uint8Array.fromBase64 があれば使う)
  function b64ToBytes(b64) {
    if (typeof Uint8Array.fromBase64 === "function") {
      try { return Uint8Array.fromBase64(b64); } catch (e) {}
    }
    var bin = atob(b64), n = bin.length, a = new Uint8Array(n);
    for (var i = 0; i < n; i++) a[i] = bin.charCodeAt(i);
    return a;
  }
  // zlib のバイト列を展開した Response (呼ぶ側で .text() / .blob())
  function inflate(bytes) {
    if (!HAS_DS) return Promise.reject(new Error("DecompressionStream is not supported"));
    return Promise.resolve(new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream("deflate"))));
  }
  // 塊を展開した文字列 (UTF-8)。drop なら展開のあと塊を DOM から外す (メモリを空ける)
  function blockText(id, drop) {
    var el = block(id);
    if (!el) return Promise.resolve(null);
    var bytes = b64ToBytes(el.textContent);
    if (drop) el.remove();
    return inflate(bytes).then(function (r) { return r.text(); });
  }

  var dataCache = {};
  window.picData = function (name) {
    if (!dataCache[name]) {
      dataCache[name] = blockText("pic-z-" + name, true).then(function (t) { return t == null ? null : JSON.parse(t); });
    }
    return dataCache[name];
  };

  window.picFileBlob = function (id) {
    var el = block("pic-f-" + id);
    if (!el) return Promise.reject(new Error("no embedded file: " + id));
    var name = el.dataset.n || "download";
    return inflate(b64ToBytes(el.textContent))
      .then(function (r) { return r.blob(); })
      .then(function (b) { return { name: name, blob: b }; });
  };
  // 展開できない環境用: 縮めたままのバイト列 (zlib)
  window.picFileRaw = function (id) {
    var el = block("pic-f-" + id);
    if (!el) return null;
    return { name: el.dataset.n || "download", bytes: b64ToBytes(el.textContent) };
  };

  // --- 開いたときの展開 ------------------------------------------------------
  var waiting = [], booted = false, domReady = document.readyState !== "loading";
  function flush() {
    if (!booted || !domReady) return;
    if (window.performance && performance.mark) { try { performance.mark("pic-ready"); } catch (e) {} }
    var q = waiting; waiting = [];
    q.forEach(function (fn) { try { fn(); } catch (e) { setTimeout(function () { throw e; }); } });
  }
  window.picReady = function (fn) { waiting.push(fn); flush(); };
  if (!domReady) document.addEventListener("DOMContentLoaded", function () { domReady = true; flush(); });

  function showError(msg) {
    var d = document.createElement("div");
    d.style.cssText = "background:#fde8e8;color:#8a1c1c;padding:10px 16px;font:14px sans-serif;border-bottom:1px solid #e0b4b4";
    d.textContent = msg;
    document.body.insertBefore(d, document.body.firstChild);
  }

  function boot() {
    if (!HAS_DS) {
      showError("このブラウザは DecompressionStream に対応していないため、図と表を表示できません " +
        "(Chrome / Edge 80 以降、Safari 16.4 以降、Firefox 113 以降で開いてください)。");
      return Promise.resolve();
    }
    var jobs = [];
    // plotly: 展開して script 要素で実行する (window.Plotly ができる)
    jobs.push(blockText("pic-z-plotly", true).then(function (code) {
      if (code == null) return;
      var s = document.createElement("script");
      s.text = code;
      document.head.appendChild(s);
    }));
    DATA_NAMES.forEach(function (name) {
      jobs.push(window.picData(name).then(function (v) { if (v != null) window[name] = v; }));
    });
    // 大きな表: 展開して元の場所に戻す
    document.querySelectorAll("template.pic-zfrag").forEach(function (tpl) {
      jobs.push(blockText(tpl.dataset.z, true).then(function (html) {
        if (html == null) return;
        tpl.insertAdjacentHTML("afterend", html);
        tpl.remove();
      }));
    });
    return Promise.all(jobs).then(function () {
      if (window.performance && performance.mark) { try { performance.mark("pic-data-ready"); } catch (e) {} }
    }, function (e) {
      showError("埋め込みデータの展開に失敗しました: " + e);
    });
  }
  boot().then(function () { booted = true; flush(); });
})();
