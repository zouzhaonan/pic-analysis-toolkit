// plotly のレンダリングを遅延実行する。
// 多数のプロットを最初に描くと重いので、details が開いた時 / details 外なら即時に描く。
(function () {
  // 文言 (差し替えるときはここだけ直す)
  // Expression の最初のグループ (pvalue 最小の 5 遺伝子) のタイトル
  var PIC_TOP_GENES_TITLE = "Top 5 DEGs";
  // サンプルを外しているときに DEG・Enrichment・Expression のタブの冒頭に出す 1 文
  var PIC_EXCL_NOTE_TEXT = "Reference values computed with all samples. Re-analysis is required for accurate statistics.";
  // 描いた直後に呼ぶフック (サンプル除外を描いたばかりの図に当てる。DOMContentLoaded で設定)
  var picAfterRender = null;
  // サンプル除外の状態 (全タブ共通・1 つだけ)。excl: sample -> 1。残りは木 (.pic-tree) から作る。
  var EX = { excl: {}, samples: [], groupOf: {}, groups: [], members: {}, ready: false };
  function exList() { return EX.samples.filter(function (s) { return EX.excl[s]; }); }
  function exAny() { return EX.ready && exList().length > 0; }
  // URL のハッシュは "#<tab>" か "#<tab>&excl=<s1>,<s2>"。タブとサンプル除外を分けて返す。
  function picParseHash() {
    var h = (location.hash || "").replace(/^#/, ""), parts = h.split("&"), x = null;
    parts.slice(1).forEach(function (q) { if (q.indexOf("excl=") === 0) x = q.slice(5); });
    return { tab: parts[0] || "", excl: x };
  }
  function exHashSuffix() {
    // 状態を読む前 (最初のタブの表示) は、URL にある除外をそのまま残す
    if (!EX.ready) { var h0 = picParseHash().excl; return h0 ? "&excl=" + h0 : ""; }
    var l = exList();
    return l.length ? "&excl=" + l.map(encodeURIComponent).join(",") : "";
  }
  function renderPlot(el) {
    if (!el || el.dataset.rendered) return;
    // 非表示 (別タブ / 閉じた details / 非表示ブロック) 内は描画しない (表示された時に描く)
    if (el.closest) {
      if (el.closest("section.pic-tab:not(.active)")) return;
      if (el.closest(".pic-subpanel:not(.active)")) return;
      if (el.closest("details:not([open])")) return;
      if (el.closest(".pic-select-item[hidden]")) return;
    }
    var spec = window.PIC_PLOTS && window.PIC_PLOTS[el.id];
    if (!spec) return;
    var layout = spec.layout || {};
    layout.autosize = true;
    var cfg = spec.config || { responsive: true, displaylogo: false };
    try {
      Plotly.newPlot(el, spec.data, layout, cfg);
      el.dataset.rendered = "1";
      if (picAfterRender) { try { picAfterRender(el); } catch (e2) {} }
      // 凡例クリックでサンプルを on/off した後、表示中のトレースだけで軸を再スケール
      el.on("plotly_legendclick", function () {
        setTimeout(function () {
          try { Plotly.relayout(el, { "yaxis.autorange": true }); } catch (e) {}
        }, 60);
        return true;
      });
      el.on("plotly_legenddoubleclick", function () {
        setTimeout(function () {
          try { Plotly.relayout(el, { "yaxis.autorange": true }); } catch (e) {}
        }, 60);
        return true;
      });
    } catch (e) {
      el.innerHTML = '<div style="color:#b00;padding:8px">plot 描画エラー: ' + e + '</div>';
    }
  }

  function renderWithin(root) {
    root.querySelectorAll(".pic-plot").forEach(renderPlot);
    root.querySelectorAll(".pic-expr").forEach(function (b) { if (b._picRefit) b._picRefit(); });
  }
  function resizeWithin(root) {
    root.querySelectorAll(".pic-plot").forEach(function (p) {
      if (p.dataset.rendered) { try { Plotly.Plots.resize(p); } catch (e) {} }
    });
  }
  // xenograft: グローバルゲノム切替。currentGeno が設定されている解析タブは
  // 論理キー key -> セクション <geno>__<key> に解決する (通常レポートは currentGeno=null)。
  var currentGeno = null;
  var currentTab = null;
  function genoSection(key) {
    var direct = document.getElementById(key);
    if (direct && direct.classList.contains("pic-tab")) return key;
    if (currentGeno) {
      var g = document.getElementById(currentGeno + "__" + key);
      if (g && g.classList.contains("pic-tab")) return currentGeno + "__" + key;
    }
    return key;
  }
  // タブ切替: 対象タブを表示し、そのタブのプロットを遅延描画 + リサイズ
  function activateTab(key) {
    if (!key) return;
    currentTab = key;
    var target = genoSection(key);
    var found = false;
    document.querySelectorAll("section.pic-tab").forEach(function (s) {
      var on = (s.id === target);
      s.classList.toggle("active", on);
      if (on) found = true;
    });
    if (!found) return;
    document.querySelectorAll(".pic-tabbtn").forEach(function (b) {
      b.classList.toggle("active", b.dataset.target === key);
    });
    var tab = document.getElementById(target);
    if (tab) { renderWithin(tab); resizeWithin(tab); }
    try { if (history.replaceState) history.replaceState(null, "", "#" + key + exHashSuffix()); } catch (e) {}
    window.scrollTo(0, 0);
  }

  // xenograft: ゲノム切替時に、現タブの UI 状態 (サブタブ / チェックボックス / ラジオ) を
  // 旧ゲノムのセクションから新ゲノムのセクションへ引き継ぐ (ゲノム間比較のため)。
  function stripGeno(s) { return s ? s.replace(/^(graft|host)__/, "") : s; }
  function syncGenoState(fromGeno, toGeno, tabKey) {
    var fromSec = document.getElementById(fromGeno + "__" + tabKey);
    var toSec = document.getElementById(toGeno + "__" + tabKey);
    if (!fromSec || !toSec) return;  // 解析タブ以外は対象外
    // 1. サブタブ (Read distribution/Sequencing depth, Heatmap/MA・Volcano, GSEA/ORA)
    var actSub = fromSec.querySelector(".pic-subtabbtn.active");
    if (actSub) {
      var suf = stripGeno(actSub.dataset.sub);
      var dstSub = Array.prototype.filter.call(toSec.querySelectorAll(".pic-subtabbtn"),
        function (x) { return stripGeno(x.dataset.sub) === suf; })[0];
      if (dstSub && !dstSub.classList.contains("active")) dstSub.click();
    }
    // 2. ラジオ (DEG contrast: data-target / GSEA matrix: data-axis+data-key)
    fromSec.querySelectorAll("input[type=radio]:checked").forEach(function (src) {
      var dst = null;
      if (src.dataset.target) {
        var s = stripGeno(src.dataset.target);
        dst = Array.prototype.filter.call(toSec.querySelectorAll("input[type=radio][data-target]"),
          function (r) { return stripGeno(r.dataset.target) === s; })[0];
      } else if (src.dataset.key && src.dataset.axis) {
        dst = toSec.querySelector("input[type=radio][data-axis='" + src.dataset.axis + "'][data-key='" + src.dataset.key + "']");
      }
      if (dst && !dst.checked) dst.click();
    });
    // 3. チェックボックス (group / sample トグル)
    fromSec.querySelectorAll("input[type=checkbox]").forEach(function (src) {
      var sel = src.dataset.group ? "[data-group='" + src.dataset.group + "']"
              : (src.dataset.sample ? "[data-sample='" + src.dataset.sample + "']" : null);
      if (!sel) return;
      toSec.querySelectorAll("input[type=checkbox]" + sel).forEach(function (dst) {
        if (dst.checked !== src.checked) { dst.checked = src.checked; dst.dispatchEvent(new Event("change", { bubbles: true })); }
      });
    });
  }
  // Downloads タブ内の graft/host ブロックを現ゲノムに合わせて表示切替
  function applyGenoDl() {
    document.querySelectorAll(".pic-geno-dl").forEach(function (d) { d.hidden = (d.dataset.geno !== currentGeno); });
  }

  // HTML 要素 (heatmap / QC / correlation テーブル) を PNG 化してダウンロードする。
  // 外部ライブラリ不要: SVG <foreignObject> に要素とレポート CSS を埋めて canvas 描画。
  function capturePng(el, filename, btn) {
    // レポート自身の CSS のみ (plotly が注入する style は < を含み SVG(XML) を壊すため除外)
    var css = "";
    document.querySelectorAll("style").forEach(function (s) {
      var t = s.textContent || "";
      if (t.indexOf(".pic-tabs") !== -1 || t.indexOf(".pic-hm") !== -1) css += t + "\n";
    });
    // スクロール/sticky を解除して全体を描く
    css += ".pic-hm-scroll{max-height:none!important;overflow:visible!important}" +
      ".pic-hm thead th,.pic-hm tbody th,.pic-hm thead th.corner," +
      ".pic-cor thead th,.pic-cor tbody th,.pic-cor .pic-cor-corner{position:static!important}" +
      ".pic-png-btn{display:none!important}";
    // クローンをオフスクリーンに一時配置してスクロール領域を展開し、フルサイズを実測する
    var clone = el.cloneNode(true);
    clone.querySelectorAll(".pic-hm-scroll").forEach(function (n) { n.style.maxHeight = "none"; n.style.overflow = "visible"; });
    clone.style.position = "absolute"; clone.style.left = "-100000px"; clone.style.top = "0"; clone.style.background = "#fff";
    document.body.appendChild(clone);
    var w = Math.ceil(clone.scrollWidth), h = Math.ceil(clone.scrollHeight);
    // 計測用の配置スタイルは描画前に外す
    clone.style.position = ""; clone.style.left = ""; clone.style.top = "";
    var xhtml;
    try { xhtml = new XMLSerializer().serializeToString(clone); }
    catch (e) { document.body.removeChild(clone); if (btn) btn.textContent = "PNG error"; return; }
    document.body.removeChild(clone);
    var svg = '<svg xmlns="http://www.w3.org/2000/svg" width="' + w + '" height="' + h + '">' +
      '<foreignObject x="0" y="0" width="' + w + '" height="' + h + '">' +
      '<div xmlns="http://www.w3.org/1999/xhtml" style="background:#fff;width:' + w + 'px">' +
      '<style><![CDATA[' + css + ']]></style>' + xhtml + '</div></foreignObject></svg>';
    var img = new Image();
    img.onload = function () {
      var scale = 2, c = document.createElement("canvas");
      c.width = Math.ceil(w * scale); c.height = Math.ceil(h * scale);
      var ctx = c.getContext("2d"); ctx.scale(scale, scale);
      ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, w, h);
      ctx.drawImage(img, 0, 0);
      c.toBlob(function (blob) {
        if (!blob) { if (btn) btn.textContent = "PNG error"; return; }
        var a = document.createElement("a");
        a.href = URL.createObjectURL(blob); a.download = filename;
        document.body.appendChild(a); a.click(); a.remove();
        setTimeout(function () { URL.revokeObjectURL(a.href); }, 2000);
        if (btn) { var t = btn.textContent; btn.textContent = "✓ saved"; setTimeout(function () { btn.textContent = t; }, 1200); }
      }, "image/png");
    };
    img.onerror = function () { if (btn) btn.textContent = "PNG error"; };
    img.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(svg);
  }

  // contrast × method の 2 軸チェックボックス。両軸が選択されたセルのみ表示。
  function initMatrix(root) {
    function refresh() {
      var rk = {}, ck = {};
      root.querySelectorAll("input[data-axis=r]").forEach(function (cb) { if (cb.checked) rk[cb.dataset.key] = 1; });
      root.querySelectorAll("input[data-axis=c]").forEach(function (cb) { if (cb.checked) ck[cb.dataset.key] = 1; });
      root.querySelectorAll(".pic-matrix-cells > .pic-select-item").forEach(function (cell) {
        var vis = rk[cell.dataset.r] && ck[cell.dataset.c];
        cell.hidden = !vis;
        if (vis) renderWithin(cell);
      });
    }
    root.querySelectorAll("input[data-axis]").forEach(function (cb) {
      cb.addEventListener("change", refresh);
    });
  }

  // ---- 遺伝子発現 (normalized counts) の box + beeswarm ----
  function hexToRgba(hex, a) {
    var m = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex || "");
    if (!m) return hex;
    return "rgba(" + parseInt(m[1], 16) + "," + parseInt(m[2], 16) + "," + parseInt(m[3], 16) + "," + a + ")";
  }
  function uniqueInOrder(arr) {
    var seen = {}, out = [];
    arr.forEach(function (x) { if (!seen[x]) { seen[x] = 1; out.push(x); } });
    return out;
  }
  function escHtml(s) {
    return String(s == null ? "" : s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
  }
  var _tip = null;
  function ensureTip() {
    if (_tip) return _tip;
    _tip = document.createElement("div");
    _tip.className = "pic-expr-tip";
    _tip.style.display = "none";
    document.body.appendChild(_tip);
    return _tip;
  }
  // 貼り付けた遺伝子リストを要素に分ける。区切りは改行・コンマ・スラッシュ・空白・タブ・セミコロン・
  // 縦棒 (全角の 、，；／ と全角空白も) とその混在。前後の空白・引用符、空の要素、重複 (大小文字を問わず) を除く。
  function picSplitGeneList(text) {
    var seen = {}, out = [];
    String(text == null ? "" : text).split(/[\s,;\/|　、，；／]+/).forEach(function (t) {
      t = t.replace(/^["'`“”‘’]+|["'`“”‘’]+$/g, "").trim();
      if (!t) return;
      var k = t.toLowerCase();
      if (seen[k]) return;
      seen[k] = 1; out.push(t);
    });
    return out;
  }
  function initExpr(box) {
    var id = box.dataset.expr;
    var data = window.PIC_EXPR && window.PIC_EXPR[id];
    if (!data || box.dataset.wired) return;
    box.dataset.wired = "1";
    var msg = box.querySelector(".pic-expr-msg");
    var plots = box.querySelector(".pic-expr-plots");
    var groups = uniqueInOrder(data.groups);

    // lookup (exact) と候補用の配列
    // lut: 小文字 -> index (従来どおり 1 つ。タグ入力の Enter と既定の遺伝子)
    // lutAll: 小文字 -> index の配列 (貼り付け用。1 つの名前に複数の ID があれば全部)。
    //   遺伝子 ID は版 (ENSG...\.12) を除いた形でも引けるようにする。
    var lut = {}, lutAll = {}, nameN = {}, nm = [];
    var cand = [];   // {disp, low, keys, i}
    function addKey(k, i) {
      if (!k) return;
      k = k.toLowerCase();
      var a = lutAll[k] || (lutAll[k] = []);
      if (a.indexOf(i) === -1) a.push(i);
    }
    function noVer(x) { return String(x).replace(/\.\d+$/, ""); }
    for (var i = 0; i < data.ens.length; i++) {
      var ens = data.ens[i] ? String(data.ens[i]) : "";
      var ext = data.ext[i] ? String(data.ext[i]) : "";
      // 名前の無い遺伝子は ID で埋めてある (ext == ens)。そのときは名前なしとして扱う
      if (ext && ens && noVer(ext).toLowerCase() === noVer(ens).toLowerCase()) ext = "";
      nm.push(ext);
      if (ens) lut[ens.toLowerCase()] = i;
      if (ext) lut[ext.toLowerCase()] = i;
      addKey(ens, i); addKey(ext, i);
      if (/\.\d+$/.test(ens)) { addKey(noVer(ens), i); if (lut[noVer(ens).toLowerCase()] === undefined) lut[noVer(ens).toLowerCase()] = i; }
      if (ext) nameN[ext.toLowerCase()] = (nameN[ext.toLowerCase()] || 0) + 1;
      var keys = [ens.toLowerCase(), noVer(ens).toLowerCase()];
      if (ext) keys.unshift(ext.toLowerCase());
      cand.push({ disp: ext ? (ext + " (" + ens + ")") : ens, low: (ext + " " + ens).toLowerCase(), keys: keys, i: i });
    }
    // 候補・タグの hover に出す「名前 (ID)」
    function candLabel(gi) { return nm[gi] ? nm[gi] + " (" + data.ens[gi] + ")" : String(data.ens[gi]); }
    // 貼り付けたリストを index に引く: {idxs, missing, multi:[{tok, idxs}]}
    // 名前でも ID でも (大文字小文字を問わず、ID の版 .N は付いていてもいなくても) 当たる。
    function resolveList(toks) {
      var idxs = [], missing = [], multi = [];
      toks.forEach(function (t) {
        var k = t.toLowerCase();
        var hits = lutAll[k] || (/\.\d+$/.test(k) ? lutAll[k.replace(/\.\d+$/, "")] : null);
        if (!hits || !hits.length) { missing.push(t); return; }
        if (hits.length > 1) multi.push({ tok: t, idxs: hits });
        hits.forEach(function (h) { if (idxs.indexOf(h) === -1) idxs.push(h); });
      });
      return { idxs: idxs, missing: missing, multi: multi };
    }

    // 遺伝子ごとに独立した y 軸 (下端 0、各遺伝子で最適スケール) を持つ小パネルを並べる。
    // 列数はコンテナ幅から決め、幅が変わって列数が変わるときだけ描き直す。
    // 図は遺伝子のグループごとに 1 枚。図ごとに状態 st = {el, idxs, ncol, shown, title} を持つ。
    var PANEL_ROW_H = 260, PANEL_VGAP = 56, PANEL_HGAP = 56;
    var setStates = [];   // 遺伝子のグループの図の状態
    function panelCols(width, n) {
      var minW = Math.max(170, 70 + groups.length * 34);
      return Math.max(1, Math.min(n, Math.floor((width + PANEL_HGAP) / (minW + PANEL_HGAP))));
    }
    function axId(i) { return i === 0 ? "" : String(i + 1); }
    function plotInto(st, idxs) {
      var target = st.el;
      st.idxs = idxs.slice();
      target.innerHTML = "";
      st.ncol = 0;
      if (!idxs.length) return;
      var width = target.clientWidth;
      if (!width) return;   // 非表示タブ: 表示されたら ResizeObserver が描く
      var on = activeSamples();
      var shownGroups = groups.filter(function (g) {
        return data.samples.some(function (sm, s) { return data.groups[s] === g && on[sm]; });
      });
      var n = idxs.length, ncol = panelCols(width, n), nrow = Math.ceil(n / ncol);
      st.ncol = ncol;
      st.shown = shownGroups;
      // タイトル付きの図 (遺伝子のグループ) は、タイトルの分だけ上を空ける
      var mar = { l: 78, r: 12, t: st.title != null ? 64 : 34, b: 50 };
      var H = mar.t + mar.b + nrow * PANEL_ROW_H + (nrow - 1) * PANEL_VGAP;
      var plotW = Math.max(1, width - mar.l - mar.r), plotH = H - mar.t - mar.b;
      var colW = (plotW - (ncol - 1) * PANEL_HGAP) / ncol;
      var div = document.createElement("div");
      div.className = "pic-plot";
      div.style.height = H + "px";
      target.appendChild(div);
      var layout = {
        showlegend: true, hovermode: false,
        legend: { orientation: "h", yanchor: "bottom", y: 1 + 6 / plotH, x: 0 },
        margin: mar,
        // 縦軸名は左端に 1 つだけ
        annotations: [{ text: "Normalized count", xref: "paper", yref: "paper", x: 0, y: 0.5,
                        xanchor: "center", yanchor: "middle", xshift: -62, textangle: -90,
                        showarrow: false, font: { size: 14 } }]
      };
      if (st.title != null) {
        layout.title = { text: "<b>" + escHtml(st.title) + "</b>", x: 0, xref: "container", xanchor: "left",
                         y: 1, yref: "container", yanchor: "top", pad: { t: 8, l: 8 }, font: { size: 16 } };
      }
      var traces = [];
      idxs.forEach(function (gi, k) {
        var r = Math.floor(k / ncol), c = k % ncol, id = axId(k);
        var x0 = c * (colW + PANEL_HGAP) / plotW, x1 = (c * (colW + PANEL_HGAP) + colW) / plotW;
        var yTop = 1 - r * (PANEL_ROW_H + PANEL_VGAP) / plotH, yBot = yTop - PANEL_ROW_H / plotH;
        var gl = shortLabel(gi);
        // x は群 (色で区別できるので目盛りは出さない)。遺伝子名は x 軸タイトルに置く
        layout["xaxis" + id] = { domain: [x0, Math.min(1, x1)], anchor: "y" + id, type: "category",
                                 categoryorder: "array", categoryarray: shownGroups, showticklabels: false,
                                 title: { text: escHtml(gl), standoff: 6, font: { size: 13 } },
                                 fixedrange: true, showline: true, linecolor: "#d7dee6", mirror: true };
        layout["yaxis" + id] = { domain: [Math.max(0, yBot), yTop], anchor: "x" + id,
                                 rangemode: "tozero", zeroline: false,
                                 showline: true, linecolor: "#d7dee6", mirror: true };
        shownGroups.forEach(function (g) {
          var ys = [], counts = data.counts[gi];
          for (var s = 0; s < data.samples.length; s++) if (data.groups[s] === g && on[data.samples[s]]) ys.push(counts[s]);
          var col = data.palette[g] || "#1f77b4";
          traces.push({
            // 凡例テキストも group 色付き。group トグルは legendgroup (素の g) で全パネル一括
            type: "box", name: '<span style="color:' + col + '">' + escHtml(g) + '</span>',
            legendgroup: g, showlegend: k === 0, xaxis: "x" + id, yaxis: "y" + id,
            x: ys.map(function () { return g; }), y: ys, boxpoints: "all", jitter: 0.5, pointpos: 0,
            marker: { color: col, size: 5 }, line: { color: col }, fillcolor: hexToRgba(col, 0.25)
          });
        });
      });
      var cfg = { responsive: true, displaylogo: false };
      if (st.title != null) cfg.toImageButtonOptions = { filename: String(st.title).replace(/[^A-Za-z0-9_.-]+/g, "_") || "genes" };
      Plotly.newPlot(div, traces, layout, cfg).then(function () {
        attachTickTips(div, idxs, st);
        div.on("plotly_afterplot", function () { attachTickTips(div, idxs, st); });
      });
      div.dataset.rendered = "1";
      // 凡例で group を on/off した後、表示中 group だけで各 y を再スケール
      div.on("plotly_legendclick", function () {
        setTimeout(function () {
          var upd = {};
          idxs.forEach(function (gi, k) { upd["yaxis" + axId(k) + ".autorange"] = true; });
          try { Plotly.relayout(div, upd); } catch (e) {}
        }, 60);
        return true;
      });
    }
    // ウィンドウ/分割幅の変化に追従: 列数が変わるときは描き直し、同じなら plotly の resize で足りる
    // 左パネルの group トグルで表示するサンプルを絞る (トグルが無い group は表示)
    function activeSamples() {
      var off = {}, on = {};
      if (EX.ready) {   // 全タブ共通のサンプル除外
        data.samples.forEach(function (sm) { on[sm] = !EX.excl[sm]; });
        return on;
      }
      box.querySelectorAll(".pic-gtoggle").forEach(function (c) { if (!c.checked) off[c.dataset.group] = 1; });
      data.samples.forEach(function (sm, s) { on[sm] = !off[data.groups[s]]; });
      return on;
    }
    function allStates() { return setStates; }
    box._picRefilter = function () { allStates().forEach(function (st) { if (st.idxs.length) plotInto(st, st.idxs); }); };
    box._picRefit = function () {
      allStates().forEach(function (st) {
        if (!st.idxs.length || !st.el.clientWidth) return;
        if (panelCols(st.el.clientWidth, st.idxs.length) !== st.ncol) plotInto(st, st.idxs);
      });
    };

    // 遺伝子名に hover -> contrast ごとの p / padj を -log10 で表示
    function fmtLog(v) {
      if (v == null || isNaN(v)) return "NA";
      if (v <= 0) return "&gt;300";
      var l = -Math.log10(v);
      return (l < 0.05 ? 0 : l).toFixed(1);
    }
    // x 軸の遺伝子ラベルに、確実に hover できる透明オーバーレイを重ねる。
    // (plotly の目盛りテキストは実操作で mouse イベントが届きにくいため)
    function attachTickTips(div, idxs, st) {
      if (!data.contrasts || !data.qval) return;
      var tip = ensureTip();
      div.style.position = "relative";
      div.querySelectorAll(".pic-expr-tickzone").forEach(function (z) { z.remove(); });
      // 各パネルの x 軸タイトル (= 遺伝子名) を拾う
      var ticks = idxs.map(function (gi, k) {
        return div.querySelector("text." + (k === 0 ? "x" : "x" + (k + 1)) + "title");
      });
      var drect = div.getBoundingClientRect();
      var fdr = (typeof data.fdr === "number") ? data.fdr : 0.1;
      // 表示中の群どうしの contrast だけを対象にする
      var shownSet = {};
      ((st && st.shown) || groups).forEach(function (g) { shownSet[g] = 1; });
      var cidx = [];
      data.contrasts.forEach(function (c, ci) {
        var ab = String(c).split(" / ");
        if (ab.length === 2 && shownSet[ab[0]] && shownSet[ab[1]]) cidx.push(ci);
      });
      ticks.forEach(function (t, k) {
        var gi = idxs[k];
        if (gi === undefined || !t) return;
        var q = data.qval[gi] || [];
        var sig = cidx.some(function (ci) { var v = q[ci]; return v != null && !isNaN(v) && v < fdr; });
        if (sig) { t.style.fill = "#d7301f"; t.style.fontWeight = "700"; }
        var r = t.getBoundingClientRect();
        var zone = document.createElement("div");
        zone.className = "pic-expr-tickzone";
        zone.style.cssText = "position:absolute;z-index:5;cursor:pointer;left:" +
          (r.left - drect.left - 8) + "px;top:" + (r.top - drect.top - 5) + "px;width:" +
          (r.width + 16) + "px;height:" + (r.height + 10) + "px;";
        (function (gi) {
          function show(e) {
            var qq = data.qval[gi] || [], pp = (data.pval && data.pval[gi]) || [];
            var rows = cidx.map(function (ci) {
              var q = qq[ci], sigq = q != null && !isNaN(q) && q < fdr;
              return "<tr><td>" + escHtml(data.contrasts[ci]) + "</td><td>" + fmtLog(pp[ci]) + "</td><td" +
                (sigq ? " class='pic-expr-tip-sig'" : "") + ">" + fmtLog(q) + "</td></tr>";
            }).join("");
            var head = "<div class='pic-expr-tip-h'>" + escHtml(candLabel(gi)) + "</div>";
            tip.innerHTML = head + (rows
              ? "<table><tr><th>contrast</th><th>&minus;log<sub>10</sub> p</th><th>&minus;log<sub>10</sub> padj</th></tr>" + rows + "</table>"
              : "<div class='pic-expr-tip-none'>no comparison among the shown groups</div>");
            tip.style.display = "block";
            // カーソルの横 (右、入らなければ左) に、縦はカーソル中心で viewport 内にクランプ
            tip.style.left = "0px"; tip.style.top = "0px";
            var tw = tip.offsetWidth, th = tip.offsetHeight;
            var vw = window.innerWidth, vh = window.innerHeight;
            var left = e.clientX + 18;
            if (left + tw > vw - 6) left = e.clientX - 18 - tw;
            if (left < 6) left = 6;
            var top = e.clientY - th / 2;
            if (top + th > vh - 6) top = vh - th - 6;
            if (top < 6) top = 6;
            tip.style.left = left + "px"; tip.style.top = top + "px";
          }
          zone.addEventListener("mouseenter", show);
          zone.addEventListener("mousemove", show);
          zone.addEventListener("mouseleave", function () { tip.style.display = "none"; });
        })(gi);
        div.appendChild(zone);
      });
    }

    // パネルの見出し。同じ名前に複数の ID があるときは ID も付けて区別する
    function shortLabel(gi) {
      var ext = nm[gi];
      if (ext && nameN[String(ext).toLowerCase()] > 1) return ext + " (" + data.ens[gi] + ")";
      return ext || data.ens[gi];
    }

    var SEP_RE = /[\s,;\/|　、，；／]/;

    // ---- 遺伝子のグループ (グループごとに 1 枚の図) ----
    // 左: グループごとにタイトルと遺伝子の入力欄 (タグ + 補完。区切りの混ざったリストの貼り付けも可)。
    // 右: グループごとの図。最初は既定の遺伝子 (pvalue 最小の 5 つ) を 1 つのグループとして出す。
    // 中身は localStorage にレポートごと (project@out_dir と expr の id) に保存する。
    var gplots = plots;
    var gside = document.createElement("div"); gside.className = "pic-expr-gsets";
    gside.innerHTML = '<h4>Gene groups</h4><div class="pic-expr-gset-list"></div>';
    var gaddBtn = document.createElement("button"); gaddBtn.type = "button"; gaddBtn.className = "pic-excl-btn pic-expr-gset-add";
    gaddBtn.textContent = "Add group"; gside.appendChild(gaddBtn);
    msg.parentNode.insertBefore(gside, msg);
    var glist = gside.querySelector(".pic-expr-gset-list");
    var meta = (typeof PIC_META !== "undefined" && PIC_META) ? PIC_META : null;
    var gbase = (meta && meta.project ? meta.project + "@" + (meta.out_dir || "") : document.title) + ":" + id;
    var gkey = "picExprSets2:" + gbase;      // {title, genes:[ID], miss:[入力のまま]}
    var gkeyOld = "picExprSets:" + gbase;    // 旧版 (Genes の欄 + {title, text} のグループ)
    var sets = [];
    function gsave() {
      try {
        localStorage.setItem(gkey, JSON.stringify(sets.map(function (g) {
          return { title: g.title, genes: g.idxs.map(function (gi) { return data.ens[gi]; }), miss: g.miss };
        })));
      } catch (e) {}
    }
    function gtitle(g) { return g.title && g.title.trim() ? g.title.trim() : "Group " + (sets.indexOf(g) + 1); }
    function gmsg() { msg.textContent = ""; }   // 操作の説明は出さない (2026-11-07 ユーザーの指示)
    // カードの下の 1 行: 遺伝子数 / 見つからない入力 / 複数の ID に当たった入力
    function gsummary(g) {
      g.sum.innerHTML = "";
      var n = g.idxs.length;
      if (n) g.sum.appendChild(document.createTextNode(n + " gene" + (n === 1 ? "" : "s")));
      if (g.miss.length) {
        var m = document.createElement("div"); m.className = "pic-expr-miss";
        m.appendChild(document.createTextNode("Not found: " + g.miss.join(", ")));
        var x = document.createElement("button"); x.type = "button"; x.className = "pic-expr-miss-x"; x.textContent = "×";
        x.title = "Clear this message";
        x.addEventListener("click", function () { g.miss = []; gsave(); gsummary(g); });
        m.appendChild(x);
        g.sum.appendChild(m);
      }
      if (g.multi && g.multi.length) {
        var d = document.createElement("div"); d.className = "pic-expr-multi";
        d.appendChild(document.createTextNode("Several IDs (all added): " + g.multi.map(function (x) {
          return x.tok + " → " + x.idxs.map(function (gi) { return data.ens[gi]; }).join(", ");
        }).join("; ")));
        g.sum.appendChild(d);
      }
    }
    function gplot(g) {
      g.st.title = gtitle(g);
      g.ti.placeholder = "Group " + (sets.indexOf(g) + 1);
      g.blk.hidden = !g.idxs.length;
      plotInto(g.st, g.idxs);
      gmsg();
    }
    function gchips(g) {
      g.tagbox.querySelectorAll(".pic-expr-chip").forEach(function (c) { c.remove(); });
      g.idxs.forEach(function (gi) {
        var chip = document.createElement("span");
        chip.className = "pic-expr-chip";
        chip.title = candLabel(gi);
        chip.appendChild(document.createTextNode(shortLabel(gi)));
        var b = document.createElement("button");
        b.type = "button"; b.textContent = "×"; b.title = "remove";
        b.addEventListener("click", function () { gremove(g, gi); g.input.focus(); });
        chip.appendChild(b);
        g.tagbox.insertBefore(chip, g.input);
      });
    }
    function grefresh(g) { gchips(g); gsummary(g); gplot(g); gsave(); }
    function gaddIdx(g, gi) {
      if (gi === undefined || gi === null || g.idxs.indexOf(gi) !== -1) return;
      g.idxs.push(gi); g.input.value = ""; hideAc(g); grefresh(g); g.input.focus();
    }
    function gremove(g, gi) { g.idxs = g.idxs.filter(function (x) { return x !== gi; }); grefresh(g); }
    // 区切りの混ざったリスト (貼り付け / Enter) をまとめて足す。見つからない入力は覚えておく
    function gaddList(g, text) {
      var r = resolveList(picSplitGeneList(text));
      r.idxs.forEach(function (gi) { if (g.idxs.indexOf(gi) === -1) g.idxs.push(gi); });
      r.missing.forEach(function (t) { if (g.miss.indexOf(t) === -1) g.miss.push(t); });
      g.multi = r.multi;
      g.input.value = ""; hideAc(g);
      grefresh(g);
    }
    // ---- 補完 (名前でも ID でも。選択済みは除く) ----
    function hideAc(g) { g.ac.style.display = "none"; g.ac.innerHTML = ""; }
    function showAc(g) {
      var t = (g.input.value || "").trim().toLowerCase();
      if (t.length < 1) { hideAc(g); return; }
      var pref = [], sub = [];
      for (var k = 0; k < cand.length && pref.length < 12; k++) {
        var c = cand[k];
        if (g.idxs.indexOf(c.i) !== -1) continue;
        if (c.keys.some(function (x) { return x.indexOf(t) === 0; })) pref.push(c);
      }
      if (pref.length < 12) {
        for (var k2 = 0; k2 < cand.length && (pref.length + sub.length) < 12; k2++) {
          var c2 = cand[k2];
          if (g.idxs.indexOf(c2.i) !== -1 || pref.indexOf(c2) !== -1) continue;
          if (c2.low.indexOf(t) !== -1) sub.push(c2);
        }
      }
      var list = pref.concat(sub);
      if (!list.length) { hideAc(g); return; }
      g.ac.innerHTML = "";
      list.forEach(function (c) {
        var it = document.createElement("div");
        it.className = "pic-expr-ac-item"; it.textContent = c.disp; it._gi = c.i;
        it.addEventListener("mousedown", function (e) { e.preventDefault(); gaddIdx(g, c.i); });
        g.ac.appendChild(it);
      });
      g.ac.style.display = "block";
    }
    function gadd(title, idxs, miss, focus) {
      var g = { title: title || "", idxs: (idxs || []).slice(), miss: (miss || []).slice(), multi: [] };
      var card = document.createElement("div"); card.className = "pic-expr-gset";
      var head = document.createElement("div"); head.className = "pic-expr-gset-head";
      var ti = document.createElement("input"); ti.type = "text"; ti.className = "pic-expr-gset-title"; ti.value = g.title;
      ti.setAttribute("aria-label", "Group title");
      var del = document.createElement("button"); del.type = "button"; del.className = "pic-expr-gset-del"; del.textContent = "×";
      del.title = "Remove this group";
      head.appendChild(ti); head.appendChild(del);
      var ta = document.createElement("div"); ta.className = "pic-expr-typeahead";
      var tagbox = document.createElement("div"); tagbox.className = "pic-expr-tagbox";
      var input = document.createElement("input"); input.type = "text"; input.className = "pic-expr-input";
      input.autocomplete = "off"; input.placeholder = "Add or paste genes…"; input.setAttribute("aria-label", "Add genes");
      tagbox.appendChild(input);
      var acel = document.createElement("div"); acel.className = "pic-expr-ac";
      ta.appendChild(tagbox); ta.appendChild(acel);
      tagbox.addEventListener("mousedown", function (e) { if (e.target === tagbox) { e.preventDefault(); input.focus(); } });
      var sum = document.createElement("div"); sum.className = "pic-expr-gset-sum";
      card.appendChild(head); card.appendChild(ta); card.appendChild(sum);
      glist.appendChild(card);
      var blk = document.createElement("div"); blk.className = "pic-expr-gblock"; blk.hidden = true;
      var gp = document.createElement("div"); gp.className = "pic-expr-gplot";
      blk.appendChild(gp); gplots.appendChild(blk);
      g.card = card; g.ti = ti; g.tagbox = tagbox; g.input = input; g.ac = acel; g.sum = sum; g.blk = blk;
      g.st = { el: gp, idxs: [], ncol: 0, title: "" };
      sets.push(g); setStates.push(g.st);
      var tmr = null;
      ti.addEventListener("input", function () {
        g.title = ti.value; gsave();
        clearTimeout(tmr); tmr = setTimeout(function () { gplot(g); }, 350);
      });
      del.addEventListener("click", function () {
        sets.splice(sets.indexOf(g), 1); setStates.splice(setStates.indexOf(g.st), 1);
        card.remove(); blk.remove();
        gsave();
        sets.forEach(gplot);   // 「Group N」の番号を詰める
        gmsg();
      });
      input.addEventListener("input", function () { showAc(g); });
      input.addEventListener("paste", function (e) {
        var t = (e.clipboardData || window.clipboardData);
        t = t ? t.getData("text") : "";
        if (!t || !SEP_RE.test(t.trim())) return;   // 1 つだけなら候補から選ぶ
        e.preventDefault();
        gaddList(g, (input.value ? input.value + " " : "") + t);
      });
      input.addEventListener("keydown", function (e) {
        if (e.key === "Enter") {
          e.preventDefault();
          var t = (input.value || "").trim();
          if (!t) return;
          if (SEP_RE.test(t)) { gaddList(g, t); return; }
          var first = g.ac.querySelector(".pic-expr-ac-item");
          var r = resolveList([t]);
          if (r.idxs.length) gaddList(g, t);
          else if (first && first._gi !== undefined) gaddIdx(g, first._gi);
          else gaddList(g, t);   // 見つからない入力として出す
        } else if (e.key === "Backspace" && !input.value && g.idxs.length) {
          gremove(g, g.idxs[g.idxs.length - 1]);
        } else if (e.key === "Escape") { hideAc(g); }
      });
      input.addEventListener("blur", function () { setTimeout(function () { hideAc(g); }, 150); });
      gchips(g); gsummary(g); gplot(g);
      if (focus) input.focus();
      return g;
    }
    gaddBtn.addEventListener("click", function () { gadd("", [], [], true); gsave(); });

    // 既定: pvalue 最小の 5 遺伝子を 1 つのグループに
    function defaultIdxs() {
      return uniqueInOrder((data.default || []).map(function (e) { return lut[String(e).toLowerCase()]; })
        .filter(function (x) { return x !== undefined; }));
    }
    var saved = null, old = null;
    try { var sv = localStorage.getItem(gkey); if (sv) saved = JSON.parse(sv); } catch (e) {}
    if (Array.isArray(saved)) {
      saved.forEach(function (x) {
        if (!x) return;
        var r = resolveList((x.genes || []).map(String));
        gadd(String(x.title || ""), r.idxs, (x.miss || []).map(String).concat(r.missing));
      });
    } else {
      gadd(PIC_TOP_GENES_TITLE, defaultIdxs(), []);
      // 旧版で作ったグループ (貼り付けた文字) を引き継ぐ
      try { var so = localStorage.getItem(gkeyOld); if (so) old = JSON.parse(so); } catch (e) {}
      if (Array.isArray(old)) old.forEach(function (x) {
        if (!x || !String(x.text || "").trim()) return;
        var r = resolveList(picSplitGeneList(String(x.text || "")));
        gadd(String(x.title || ""), r.idxs, r.missing);
      });
      gsave();
    }
    gmsg();
    if (window.ResizeObserver) {
      var _roTimer = null;
      new ResizeObserver(function () {
        clearTimeout(_roTimer);
        _roTimer = setTimeout(box._picRefit, 120);
      }).observe(gplots);
    }
  }

  // 埋め込みデータ (図・発現・PCA・plotly・大きな表) の展開を待ってから始める (report_data.js)
  picReady(function () {
    // タブ: ボタンで切替。初期は URL ハッシュ or 先頭タブ。
    // xenograft: グローバルゲノム切替 (graft/host)。initial タブ解決の前に設定する。
    var genobtns = document.querySelectorAll(".pic-genobtn");
    if (genobtns.length) {
      currentGeno = (document.querySelector(".pic-genobtn.active") || genobtns[0]).dataset.geno;
      applyGenoDl();  // Downloads の graft/host ブロック初期表示
      genobtns.forEach(function (b) {
        b.addEventListener("click", function () {
          var from = currentGeno;
          currentGeno = b.dataset.geno;
          genobtns.forEach(function (x) { x.classList.toggle("active", x === b); });
          applyGenoDl();
          if (currentTab && from && from !== currentGeno) syncGenoState(from, currentGeno, currentTab);
          if (currentTab) {
            activateTab(currentTab);
            // 新セクションの plot を描画後に group トグル状態を再適用
            var newSec = document.getElementById(currentGeno + "__" + currentTab);
            if (newSec) applyGroupToggles(newSec);
          }
        });
      });
    }
    var tabbtns = document.querySelectorAll(".pic-tabbtn");
    if (tabbtns.length) {
      tabbtns.forEach(function (b) {
        b.addEventListener("click", function () { activateTab(b.dataset.target); });
      });
      var targets = Array.prototype.map.call(tabbtns, function (b) { return b.dataset.target; });
      var initial = picParseHash().tab;
      if (!initial || targets.indexOf(initial) === -1) {
        var activeBtn = document.querySelector(".pic-tabbtn.active") || tabbtns[0];
        initial = activeBtn ? activeBtn.dataset.target : "";
      }
      activateTab(initial);
    }
    // Overview の解説カード / Downloads の「open tab →」でそのタブへ
    document.querySelectorAll(".pic-ov-step[data-target],.pic-ov-jump[data-target]").forEach(function (el) {
      function go() { activateTab(el.dataset.target); }
      el.addEventListener("click", function (e) { e.preventDefault(); e.stopPropagation(); go(); });
      el.addEventListener("keydown", function (e) {
        if (e.key === "Enter" || e.key === " ") { e.preventDefault(); e.stopPropagation(); go(); }
      });
    });
    // アクティブなタブのプロットを描画 (renderPlot が非表示を自己判定)
    renderWithin(document.body);
    // details は開いた時に描画
    document.querySelectorAll("details").forEach(function (d) {
      // 既に open のものは描画
      if (d.open) renderWithin(d);
      d.addEventListener("toggle", function () {
        if (d.open) renderWithin(d);
      });
    });
    // ラジオ/チェックボックス (バー or group 行列) でブロックを表示/非表示
    function showTarget(inp) {
      var t = document.getElementById(inp.dataset.target); if (!t) return;
      if (inp.checked) { t.hidden = false; renderWithin(t); resizeWithin(t); } else { t.hidden = true; }
    }
    document.querySelectorAll("input[data-target]").forEach(function (cb) {
      cb.addEventListener("change", function () {
        if (cb.type === "radio" && cb.name) {
          document.querySelectorAll('input[data-target][name="' + cb.name + '"]').forEach(showTarget);
        } else { showTarget(cb); }
      });
    });
    // セクション内サブタブの切替 (ボタンはサブパネルごとに複製されるため data-sub で照合)
    document.querySelectorAll(".pic-subtabbtn").forEach(function (b) {
      b.addEventListener("click", function () {
        var sec = b.closest("section.pic-tab"); if (!sec) return;
        var key = b.dataset.sub;
        sec.querySelectorAll(".pic-subtabbtn").forEach(function (x) { x.classList.toggle("active", x.dataset.sub === key); });
        sec.querySelectorAll(".pic-subpanel").forEach(function (p) {
          var on = (p.id === key);
          p.classList.toggle("active", on);
          if (on) { renderWithin(p); resizeWithin(p); }
        });
        sec.querySelectorAll(".pic-method-sel").forEach(function (w) { w.hidden = (w.dataset.sub !== key); });
      });
    });
    // タブの見出しの右の method のリスト (Enrichment)。サブパネルごとに 1 つ作り、表示中のものだけ出す。
    // 元のラジオ (.pic-method-src、CSS で隠す) はそのまま残し、リストから操作する
    // (表示の切り替え・ゲノム間の引き継ぎはラジオの change の仕組みのまま)。
    document.querySelectorAll("section.pic-tab").forEach(function (sec) {
      var host = sec.querySelector(":scope > .pic-tab-head > .pic-head-methods");
      if (!host) return;
      sec.querySelectorAll(".pic-subpanel").forEach(function (sp) {
        var src = sp.querySelector(".pic-method-src"); if (!src) return;
        var radios = Array.prototype.slice.call(src.querySelectorAll("input[type=radio]"));
        if (!radios.length) return;
        var wrap = document.createElement("label"); wrap.className = "pic-method-sel"; wrap.dataset.sub = sp.id;
        var t = document.createElement("span"); t.textContent = src.dataset.title || "Method";
        var sel = document.createElement("select");
        radios.forEach(function (r, i) {
          var o = document.createElement("option"); o.value = String(i);
          o.textContent = (r.closest("label") || r.parentNode).textContent.trim();
          sel.appendChild(o);
        });
        function sync() { radios.forEach(function (r, i) { if (r.checked) sel.value = String(i); }); }
        sync();
        sel.addEventListener("change", function () { var r = radios[+sel.value]; if (r && !r.checked) r.click(); });
        radios.forEach(function (r) { r.addEventListener("change", sync); });
        wrap.appendChild(t); wrap.appendChild(sel); host.appendChild(wrap);
        wrap.hidden = !sp.classList.contains("active");
      });
    });
    // GSEA の contrast × method マトリクスを初期化
    document.querySelectorAll(".pic-matrix").forEach(initMatrix);
    // contrast 行列: セルに合わせたら行・列とヘッダを強調
    document.querySelectorAll("table.pic-cmatrix").forEach(function (tbl) {
      var bodyRows = tbl.querySelectorAll("tbody tr");
      var head = tbl.querySelector("thead tr");
      function clear() {
        tbl.querySelectorAll(".pic-cmx-rowhl,.pic-cmx-colhl,.pic-cmx-hl,.pic-cmx-hit").forEach(function (e) {
          e.classList.remove("pic-cmx-rowhl", "pic-cmx-colhl", "pic-cmx-hl", "pic-cmx-hit");
        });
      }
      tbl.querySelectorAll("td.pic-cmx-cell").forEach(function (cell) {
        cell.addEventListener("mouseenter", function () {
          clear();
          var tr = cell.parentElement;
          var idx = Array.prototype.indexOf.call(tr.children, cell);
          Array.prototype.forEach.call(tr.children, function (c) { c.classList.add("pic-cmx-rowhl"); });
          bodyRows.forEach(function (r) { if (r.children[idx]) r.children[idx].classList.add("pic-cmx-colhl"); });
          if (tr.children[0]) tr.children[0].classList.add("pic-cmx-hl");
          if (head && head.children[idx]) head.children[idx].classList.add("pic-cmx-hl");
          cell.classList.add("pic-cmx-hit");
        });
      });
      tbl.addEventListener("mouseleave", clear);
    });
    // 遺伝子発現 UI を初期化
    document.querySelectorAll(".pic-expr").forEach(initExpr);
    // 左右分割のドラッグ・リサイズ。--ctrl-w を :root に設定 → 全タブ・全ウィジェットに反映
    (function initSplitters() {
      var SEL = ".pic-2pane, .pic-matrix, .pic-degsel, .pic-selgrid, .pic-expr, .pic-hmsel";
      var MIN = 170, MAX = 620;
      try {
        var saved = localStorage.getItem("picCtrlW");
        if (saved && /^\d+px$/.test(saved)) document.documentElement.style.setProperty("--ctrl-w", saved);
      } catch (e) {}
      document.querySelectorAll(SEL).forEach(function (grid) {
        if (grid.classList.contains("pic-1pane") || grid.dataset.split) return;
        grid.dataset.split = "1";
        var handle = document.createElement("div");
        handle.className = "pic-splitter";
        handle.title = "Drag to resize (applies to all tabs)";
        grid.appendChild(handle);
        handle.addEventListener("mousedown", function (e) {
          e.preventDefault();
          var rect = grid.getBoundingClientRect();
          handle.classList.add("pic-drag");
          document.body.classList.add("pic-resizing");
          function move(ev) {
            var w = Math.round(Math.max(MIN, Math.min(MAX, ev.clientX - rect.left)));
            document.documentElement.style.setProperty("--ctrl-w", w + "px");
          }
          function up() {
            document.removeEventListener("mousemove", move);
            document.removeEventListener("mouseup", up);
            handle.classList.remove("pic-drag");
            document.body.classList.remove("pic-resizing");
            var cur = getComputedStyle(document.documentElement).getPropertyValue("--ctrl-w").trim();
            try { if (cur) localStorage.setItem("picCtrlW", cur); } catch (e2) {}
            // リサイズ後、表示中のプロットを再描画してレイアウトを合わせる
            var act = document.querySelector("section.pic-tab.active");
            if (act) resizeWithin(act);
          }
          document.addEventListener("mousemove", move);
          document.addEventListener("mouseup", up);
        });
      });
    })();
    // ---- 部分集合 PCA -------------------------------------------------------
    // PCA タブの group トグルは「表示の on/off」ではなく「PCA の計算対象」を選ぶ。
    // 選ばれた群のサンプルだけで分散上位 ntop を選び直し、中心化/標準化して主成分を
    // 引き直す (R の prcomp(t(rld_top), scale.=TRUE, center=TRUE) と同じ手順)。
    // 遺伝子数 p >> サンプル数 n なので、n×n のグラム行列 G = XXt を Jacobi 法で
    // 対角化すれば PC 座標 = U*sqrt(lambda) が厳密に得られる (SVD と等価)。
    // 対称行列の固有値分解 (Numerical Recipes の巡回 Jacobi)。固有値の降順で返す。
    function picJacobiEigen(A, n) {
      var a = [], v = [], i, j, ip, iq;
      for (i = 0; i < n; i++) {
        a.push(A[i].slice()); v.push([]);
        for (j = 0; j < n; j++) v[i].push(i === j ? 1 : 0);
      }
      for (var sweep = 0; sweep < 100; sweep++) {
        var sm = 0;
        for (ip = 0; ip < n - 1; ip++) for (iq = ip + 1; iq < n; iq++) sm += Math.abs(a[ip][iq]);
        if (sm < 1e-13) break;
        for (ip = 0; ip < n - 1; ip++) {
          for (iq = ip + 1; iq < n; iq++) {
            if (Math.abs(a[ip][iq]) < 1e-300) { a[ip][iq] = 0; continue; }
            var g = 100 * Math.abs(a[ip][iq]), h = a[iq][iq] - a[ip][ip], t;
            if (Math.abs(h) + g === Math.abs(h)) { t = a[ip][iq] / h; }
            else {
              var theta = 0.5 * h / a[ip][iq];
              t = 1 / (Math.abs(theta) + Math.sqrt(1 + theta * theta));
              if (theta < 0) t = -t;
            }
            var c = 1 / Math.sqrt(1 + t * t), s = t * c, tau = s / (1 + c);
            h = t * a[ip][iq];
            a[ip][ip] -= h; a[iq][iq] += h; a[ip][iq] = 0;
            var rot = function (m, i1, j1, i2, j2) {
              var g1 = m[i1][j1], h1 = m[i2][j2];
              m[i1][j1] = g1 - s * (h1 + g1 * tau);
              m[i2][j2] = h1 + s * (g1 - h1 * tau);
            };
            for (j = 0; j < ip; j++) rot(a, j, ip, j, iq);
            for (j = ip + 1; j < iq; j++) rot(a, ip, j, j, iq);
            for (j = iq + 1; j < n; j++) rot(a, ip, j, iq, j);
            for (j = 0; j < n; j++) rot(v, j, ip, j, iq);
          }
        }
      }
      var idx = [];
      for (i = 0; i < n; i++) idx.push(i);
      idx.sort(function (p, q) { return a[q][q] - a[p][p]; });
      return {
        values: idx.map(function (p) { return a[p][p]; }),
        vectors: idx.map(function (p) { return v.map(function (row) { return row[p]; }); })
      };
    }
    // keepIdx のサンプルだけで PCA を計算する。戻り: {coords:[PCごとの座標], varpct:[...]}
    function picSubsetPCA(pd, keepIdx, ntop) {
      var G = pd.mat.length, n = keepIdx.length, i, k, kk, q;
      if (n < 3) return null;
      var stat = [];
      for (i = 0; i < G; i++) {
        var row = pd.mat[i], m = 0;
        for (k = 0; k < n; k++) m += row[keepIdx[k]];
        m /= n;
        var ss = 0;
        for (k = 0; k < n; k++) { var dd = row[keepIdx[k]] - m; ss += dd * dd; }
        stat.push({ i: i, mean: m, sd: Math.sqrt(ss / (n - 1)) });
      }
      // 分散 0 の遺伝子は標準化できない (prcomp も同様) ので除外してから上位 ntop
      var usable = stat.filter(function (s) { return s.sd > 1e-9; });
      usable.sort(function (p2, q2) { return q2.sd - p2.sd; });
      var sel = usable.slice(0, Math.min(ntop, usable.length)), p = sel.length;
      if (p < 2) return null;
      var X = [];
      for (k = 0; k < n; k++) X.push(new Array(p));
      for (q = 0; q < p; q++) {
        var st = sel[q], r = pd.mat[st.i];
        for (k = 0; k < n; k++) X[k][q] = (r[keepIdx[k]] - st.mean) / st.sd;
      }
      var Gm = [];
      for (k = 0; k < n; k++) Gm.push(new Array(n));
      for (k = 0; k < n; k++) for (kk = k; kk < n; kk++) {
        var acc = 0;
        for (q = 0; q < p; q++) acc += X[k][q] * X[kk][q];
        Gm[k][kk] = acc; Gm[kk][k] = acc;
      }
      var e = picJacobiEigen(Gm, n), tot = 0;
      for (k = 0; k < n; k++) tot += Math.max(0, e.values[k]);
      // 寄与率も座標も全 PC 分を返す (scree は PC1..PC15、CSV 出力は全 PC を使う)。
      var varpct = [];
      for (k = 0; k < n; k++) varpct.push(tot > 0 ? 100 * Math.max(0, e.values[k]) / tot : 0);
      var npc = n, coords = [];
      for (var pc = 0; pc < npc; pc++) {
        var lam = Math.max(0, e.values[pc]), sq = Math.sqrt(lam), vec = e.vectors[pc], mi = 0;
        // 符号は任意なので、絶対値最大の成分が正になるよう固定する
        // (選択を変えるたびに軸が反転して見えるのを防ぐ)
        for (k = 1; k < n; k++) if (Math.abs(vec[k]) > Math.abs(vec[mi])) mi = k;
        var sg = vec[mi] < 0 ? -1 : 1, col = [];
        for (k = 0; k < n; k++) col.push(sg * vec[k] * sq);
        coords.push(col);
      }
      return { coords: coords, varpct: varpct, n: n, p: p };
    }
    function picPCAWarn(tab, msg) {
      var host = tab.querySelector(".pic-ctrl");
      if (!host) return;
      var el = host.querySelector(".pic-ctrl-warn");
      if (!msg) { if (el && el.parentNode) el.parentNode.removeChild(el); return; }
      if (!el) { el = document.createElement("p"); el.className = "pic-ctrl-warn"; host.appendChild(el); }
      el.textContent = msg;
    }
    // 初期状態 (R の prcomp が出した PCA) を退避しておく。全群を選び直したときは
    // 再計算値ではなくこれを復元する: ブラウザ側の計算は R と数値的に一致するが、
    // 埋め込み行列は小数を丸めてあるため完全に同じ値にはならない。全群のときは
    // 元の図がそのまま出るほうが誤解がない。
    function picStashPCA0(div) {
      if (!div || !div.data || div._picPCA0) return;
      var ax = function (a) { return a && a.title ? (a.title.text != null ? a.title.text : a.title) : null; };
      div._picPCA0 = {
        x: div.data.map(function (t) { return (t.x || []).slice(); }),
        y: div.data.map(function (t) { return (t.y || []).slice(); }),
        text: div.data.map(function (t) { return (t.text || []).slice(); }),          // replicate 番号
        hovertext: div.data.map(function (t) { return (t.hovertext || []).slice(); }), // サンプル名
        xt: ax(div.layout && div.layout.xaxis), yt: ax(div.layout && div.layout.yaxis)
      };
    }
    function picRestorePCA0(div) {
      if (!div || !div._picPCA0 || !div.data) return;
      var s = div._picPCA0, idx = div.data.map(function (_, i) { return i; });
      try {
        Plotly.restyle(div, { x: s.x, y: s.y, text: s.text, hovertext: s.hovertext }, idx);
        if (s.xt != null) Plotly.relayout(div, { "xaxis.title.text": s.xt, "yaxis.title.text": s.yt });
      } catch (e) {}
    }
    // 現在チェックされているサンプルの添字。計算対象はサンプル側のチェックだけで決まる
    // (group はその一括トグル)。sample トグルを持たない旧レイアウトでは group で判定。
    function picPCAKeepIdx(tab, pd) {
      var keepS = {}, hasS = false, keepG = {}, out = [], i;
      // 全タブ共通の除外の状態があればそれに従う (木を持つレポート)
      if (EX.ready) {
        for (i = 0; i < pd.samples.length; i++) if (!EX.excl[pd.samples[i]]) out.push(i);
        return out;
      }
      tab.querySelectorAll(".pic-stoggle").forEach(function (c) { hasS = true; if (c.checked) keepS[c.dataset.sample] = 1; });
      tab.querySelectorAll(".pic-gtoggle").forEach(function (c) { if (c.checked) keepG[c.dataset.group] = 1; });
      for (i = 0; i < pd.samples.length; i++) {
        if (hasS ? !keepS[pd.samples[i]] : !keepG[pd.groups[i]]) continue;
        out.push(i);
      }
      return out;
    }
    // 選択中のサンプルだけで引き直した PCA 座標を CSV で保存する
    // (R が出す PCA_RegLog_*.csv と同じ sample, PC1..PCn, group の並び)。
    function picPCADownloadSelected(tab) {
      if (!tab) return;
      var pd = (typeof PIC_PCA !== "undefined" && PIC_PCA) ? PIC_PCA[tab.id] : null;
      if (!pd || !pd.mat) return;
      var keep = picPCAKeepIdx(tab, pd);
      if (keep.length < 3) { alert("Select at least 3 samples to export a PCA."); return; }
      var res = picSubsetPCA(pd, keep, pd.ntop || 500);
      if (!res) { alert("Not enough usable genes for the current selection."); return; }
      var npc = res.coords.length, head = ["sample"], pc, k;
      var csvq = function (x) { x = String(x); return /[",\r\n]/.test(x) ? '"' + x.replace(/"/g, '""') + '"' : x; };
      for (pc = 0; pc < npc; pc++) head.push("PC" + (pc + 1));
      head.push("group");
      var rows = [head.join(",")];
      for (k = 0; k < keep.length; k++) {
        var row = [csvq(pd.samples[keep[k]])];
        for (pc = 0; pc < npc; pc++) row.push(res.coords[pc][k]);
        row.push(csvq(pd.groups[keep[k]]));
        rows.push(row.join(","));
      }
      saveBlob(new Blob([rows.join("\n") + "\n"], { type: "text/csv" }),
               "PCA_RegLog_" + (pd.project || "selected") + "_" + keep.length + "samples.csv");
    }
    // PCA タブなら再計算して true を返す。PCA データが無ければ false (従来処理へ)。
    function picRecomputePCA(tab) {
      var pd = (typeof PIC_PCA !== "undefined" && PIC_PCA) ? PIC_PCA[tab.id] : null;
      if (!pd || !pd.mat || !pd.mat.length) return false;
      var pairs = [[0, 1], [1, 2]];
      var divs = pairs.map(function (pr) {
        return document.getElementById(tab.id + "_PC" + (pr[0] + 1) + "_PC" + (pr[1] + 1));
      });
      var scree = document.getElementById(tab.id + "_scree");
      divs.concat([scree]).forEach(picStashPCA0);

      var keepIdx = picPCAKeepIdx(tab, pd);

      // 全群を選んでいるときは R が計算した元の PCA に戻す
      if (keepIdx.length === pd.samples.length) {
        divs.concat([scree]).forEach(picRestorePCA0);
        picPCAWarn(tab, "");
        return true;
      }
      var res = keepIdx.length >= 3 ? picSubsetPCA(pd, keepIdx, pd.ntop || 500) : null;
      if (!res) {
        picPCAWarn(tab, "Select at least 3 samples to compute a PCA.");
        return true;
      }
      picPCAWarn(tab, "");
      // sample -> 新しい PC 座標
      var coord = {};
      for (var k = 0; k < res.n; k++) {
        var s = pd.samples[keepIdx[k]], row = [];
        for (var pc = 0; pc < res.coords.length; pc++) row.push(res.coords[pc][k]);
        coord[s] = row;
      }
      var lbl = function (pc) {
        return "PC" + (pc + 1) + (res.varpct[pc] != null ? " (" + res.varpct[pc].toFixed(1) + "%)" : "");
      };
      // 散布図: trace 構造 (色 / batch 分割 / legendgroup) は保ち、座標だけ差し替える
      pairs.forEach(function (pr, pi) {
        var div = divs[pi];
        if (!div || !div.data) return;
        var xs = [], ys = [], ts = [], hs = [], idxs = [];
        div.data.forEach(function (t, ti) {
          var s0 = div._picPCA0;
          var names = (s0 ? s0.hovertext[ti] : t.hovertext) || [];  // サンプル名 (座標の引き当てに使う)
          var reps = (s0 ? s0.text[ti] : t.text) || [];             // replicate 番号 (点のラベル)
          var nx = [], ny = [], nt = [], nh = [];
          for (var m = 0; m < names.length; m++) {
            var nm = names[m];
            if (coord[nm] && coord[nm].length > pr[1]) {
              nx.push(coord[nm][pr[0]]); ny.push(coord[nm][pr[1]]);
              nt.push(reps[m]); nh.push(nm);
            }
          }
          xs.push(nx); ys.push(ny); ts.push(nt); hs.push(nh); idxs.push(ti);
        });
        try {
          Plotly.restyle(div, { x: xs, y: ys, text: ts, hovertext: hs }, idxs);
          Plotly.relayout(div, { "xaxis.title.text": lbl(pr[0]), "yaxis.title.text": lbl(pr[1]) });
        } catch (e) {}
      });
      // scree (寄与率): R 側と同じく最大 15 本まで描く
      if (scree && scree.data) {
        var px = [], py = [], nb = Math.min(res.varpct.length, 15);
        for (var pc2 = 0; pc2 < nb; pc2++) { px.push("PC" + (pc2 + 1)); py.push(res.varpct[pc2]); }
        try { Plotly.restyle(scree, { x: [px], y: [py] }, [0]); } catch (e) {}
      }
      return true;
    }
    // ---- サンプル除外 (全タブ共通) --------------------------------------------
    // どのタブの group -> sample の木 (.pic-tree) を触っても、状態 EX.excl は 1 つ。
    // 変えたら全タブの木を揃え、localStorage と URL (#<tab>&excl=...) に保存し、全タブに当てる。
    // タブごとの扱い (計算し直すのは PCA だけ。ほかは表示から外すだけ):
    //   QC (Read distribution / Sequencing depth / Aggregation) … サンプルごとの値なので行・線を外す
    //       (データバーは残りのサンプルの最大で引き直す)
    //   Correlation … 行・列を外す。値はそのまま、色の幅だけ表示中のサンプル (対角を除く) の最小〜最大に引き直す
    //   PCA … 残りのサンプルで主成分を引き直す (picRecomputePCA)
    //   DEG / Expression / Enrichment … DESeq2 は全サンプルで計算したもの。外したサンプルの列・点を隠し、
    //       丸ごと外した group を含む比較を隠し、参考値である旨を 1 文 (PIC_EXCL_NOTE_TEXT) で注記する。
    // 除外が無いときは何もしない (図・表は元のまま)。
    var PIC_META_ = (typeof PIC_META !== "undefined") ? PIC_META : null;
    function exInitModel() {
      document.querySelectorAll(".pic-tree .pic-stoggle").forEach(function (c) {
        var s = c.dataset.sample, g = c.dataset.sgroup;
        if (s == null || EX.groupOf[s] !== undefined) return;
        EX.samples.push(s); EX.groupOf[s] = g;
        if (!EX.members[g]) { EX.members[g] = []; EX.groups.push(g); }
        EX.members[g].push(s);
      });
      EX.ready = EX.samples.length > 0;
    }
    function exStoreKey() {
      var m = PIC_META_;
      return "picExcl:" + (m && m.project ? m.project + "@" + (m.out_dir || "") : document.title);
    }
    function exLoad() {
      var h = picParseHash(), list = null;
      if (h.excl !== null) {
        list = h.excl ? h.excl.split(",").map(function (x) { try { return decodeURIComponent(x); } catch (e) { return x; } }) : [];
      } else {
        try { var v = localStorage.getItem(exStoreKey()); if (v) list = JSON.parse(v); } catch (e) {}
      }
      EX.excl = {};
      (list || []).forEach(function (s) { if (EX.groupOf[s] !== undefined) EX.excl[s] = 1; });
    }
    function exSave() {
      var l = exList();
      try { if (l.length) localStorage.setItem(exStoreKey(), JSON.stringify(l)); else localStorage.removeItem(exStoreKey()); } catch (e) {}
      try {
        if (history.replaceState) history.replaceState(null, "", "#" + (picParseHash().tab || currentTab || "") + exHashSuffix());
      } catch (e) {}
    }
    // group が丸ごと除外されているか
    function exGroupOut(g) { var m = EX.members[g]; return !!m && m.every(function (s) { return EX.excl[s]; }); }
    function exGroupKey(g) { return String(g == null ? "" : g).toLowerCase(); }
    function exOutGroups() { return EX.groups.filter(exGroupOut); }

    // 木のチェックを状態に合わせる (group は配下から checked / indeterminate)
    function exSyncTrees() {
      document.querySelectorAll(".pic-tree .pic-stoggle").forEach(function (c) { c.checked = !EX.excl[c.dataset.sample]; });
      document.querySelectorAll(".pic-tree .pic-gtoggle").forEach(function (g) {
        var m = EX.members[g.dataset.group] || [];
        var on = m.filter(function (s) { return !EX.excl[s]; }).length;
        g.checked = on > 0;
        g.indeterminate = on > 0 && on < m.length;
      });
      document.querySelectorAll(".pic-tree").forEach(function (t) {
        var foot = t.nextElementSibling && t.nextElementSibling.classList.contains("pic-excl-foot") ? t.nextElementSibling : null;
        var l = exList();
        if (!l.length) { if (foot) foot.parentNode.removeChild(foot); return; }
        if (!foot) {
          foot = document.createElement("div"); foot.className = "pic-excl-foot";
          t.parentNode.insertBefore(foot, t.nextSibling);
        }
        foot.innerHTML = "";
        var txt = document.createElement("span");
        txt.textContent = l.length + " of " + EX.samples.length + " samples excluded";
        foot.appendChild(txt);
        var rb = document.createElement("button"); rb.type = "button"; rb.className = "pic-excl-btn"; rb.textContent = "Reset";
        rb.title = "Include all samples again";
        rb.addEventListener("click", function () { EX.excl = {}; exCommit(); });
        foot.appendChild(rb);
      });
    }
    function exCommit() { exSyncTrees(); exSave(); exApplyAll(); }

    // ---- タブごとの注記 ----
    function exTabKind(tab) {
      var id = (tab.id || "").replace(/^(graft|host)__/, "");
      if (id === "qc" || id === "aggregate") return "qc";
      if (id === "cor" || id === "pca" || id === "deg" || id === "expr" || id === "enrich") return id;
      return "";
    }
    // DESeq2 / enrichment の結果 (DEG・Expression・Enrichment) は全サンプルで計算したもの。
    // サンプルを外しているときは、これらのタブの冒頭に 1 文 (PIC_EXCL_NOTE_TEXT) だけ出す。
    function exNote(tab) {
      var old = tab.querySelector(":scope > .pic-excl-note");
      if (old) old.parentNode.removeChild(old);
      var kind = exTabKind(tab);
      if (!(kind === "deg" || kind === "expr" || kind === "enrich") || !exAny()) return;
      var box = document.createElement("div");
      box.className = "pic-excl-note pic-excl-warn";
      var p = document.createElement("p"); p.textContent = PIC_EXCL_NOTE_TEXT;
      box.appendChild(p);
      var head = tab.querySelector(":scope > .pic-tab-head");
      if (head) head.parentNode.insertBefore(box, head.nextSibling); else tab.insertBefore(box, tab.firstChild);
    }

    // ---- 行 (QC の棒・深さの表・相関の行) ----
    // 区切り線 (group の先頭) は見えている行の中で引き直す。除外が無ければ元のクラスに戻す。
    function exResepRows(rows, cls) {
      var prev = null, first = true;
      rows.forEach(function (r) {
        if (r._picCls0 === undefined) r._picCls0 = r.className;
        if (!exAny()) { r.className = r._picCls0; return; }
        if (r.style.display === "none") return;
        var g = r.dataset.group || r.dataset.grow || "";
        var base = r._picCls0.split(/\s+/).filter(function (c) { return c && c !== cls; });
        if (!first && g !== prev) base.push(cls);
        r.className = base.join(" ");
        prev = g; first = false;
      });
    }
    function exApplyRows(tab) {
      var any = exAny();
      tab.querySelectorAll(".pic-bars").forEach(function (box) {
        var rows = Array.prototype.slice.call(box.querySelectorAll(".pic-bar-row[data-sample]"));
        rows.forEach(function (r) { r.style.display = (any && EX.excl[r.dataset.sample]) ? "none" : ""; });
        exResepRows(rows, "pic-grp-start");
      });
      tab.querySelectorAll("table.pic-databar-table").forEach(function (tbl) {
        var rows = Array.prototype.slice.call(tbl.querySelectorAll("tbody > tr[data-sample]"));
        rows.forEach(function (r) { r.style.display = (any && EX.excl[r.dataset.sample]) ? "none" : ""; });
        exResepRows(rows, "pic-grp-start");
        // データバー: 列ごとに、見えている行の最大で長さを引き直す (R と同じ 100 * v / max)
        var ncol = rows.length ? rows[0].children.length : 0;
        for (var c = 0; c < ncol; c++) {
          var cells = rows.map(function (r) { return r.children[c]; });
          if (!cells.length || !cells[0] || !cells[0].classList.contains("pic-td-bar")) continue;
          var mx = 0, vals = cells.map(function (td, i) {
            var v = parseFloat(((td.querySelector(".pic-databar-val") || {}).textContent || "").replace(/,/g, ""));
            if (!isFinite(v)) v = 0;
            if (rows[i].style.display !== "none" && v > mx) mx = v;
            return v;
          });
          cells.forEach(function (td, i) {
            var bar = td.querySelector(".pic-databar"); if (!bar) return;
            if (bar._picSt0 === undefined) bar._picSt0 = bar.getAttribute("style") || "";
            if (!any) { bar.setAttribute("style", bar._picSt0); return; }
            var pct = (mx > 0 ? 100 * vals[i] / mx : 0).toFixed(3);
            bar.setAttribute("style", bar._picSt0.replace(/[\d.]+%,transparent [\d.]+%/, pct + "%,transparent " + pct + "%"));
          });
        }
      });
    }

    // ---- 相関 (R の cor_color / contrast_text と同じ配色) ----
    function exRamp3(t) {
      var c0 = [0x21, 0x66, 0xac], c1 = [0xf7, 0xf7, 0xf7], c2 = [0xb2, 0x18, 0x2b], a, b, u;
      if (t <= 0.5) { a = c0; b = c1; u = t / 0.5; } else { a = c1; b = c2; u = (t - 0.5) / 0.5; }
      return [0, 1, 2].map(function (k) { return a[k] + (b[k] - a[k]) * u; });
    }
    function exHex(rgb) {
      // R の rgb(..., maxColorValue = 255) は整数に切り捨てる (四捨五入ではない)
      return "#" + rgb.map(function (v) { var x = Math.floor(v); x = Math.max(0, Math.min(255, x)); return (x < 16 ? "0" : "") + x.toString(16).toUpperCase(); }).join("");
    }
    // 相関の色: lo..hi を青-白-赤に写す (R の cor_color は hi = 1 の場合)
    function exCorColor(v, lo, hi) {
      var t = (v - lo) / Math.max(1e-9, hi - lo);
      if (!isFinite(t)) t = 0.5;
      t = Math.max(0, Math.min(1, t));
      return exHex(exRamp3(t));
    }
    function exTextColor(hex) {
      var r = parseInt(hex.substr(1, 2), 16), g = parseInt(hex.substr(3, 2), 16), b = parseInt(hex.substr(5, 2), 16);
      return (0.299 * r + 0.587 * g + 0.114 * b) < 140 ? "#ffffff" : "#1f2933";
    }
    // 列 (見出しの th[data-sample]) を隠し、区切り線のクラスを見えている列で引き直す表の共通部分
    function exColInfo(tbl) {
      var hdr = tbl.querySelector("thead tr");
      var ths = hdr ? Array.prototype.slice.call(hdr.children) : [];
      var cols = [];
      ths.forEach(function (th, k) { if (th.dataset && th.dataset.sample != null) cols.push({ k: k, s: th.dataset.sample, g: th.dataset.gcol || "" }); });
      return cols;
    }
    // 色の幅の凡例 (カラーバーと数字)。表の下 (PNG に写る .pic-capbox の中) に置く。
    // 全部表示のときは R と同じ幅 (全セルの最小 .. 1)、除外したときは表示中のサンプルの
    // 対角を除く値の最小 .. 最大。
    function exCorLegend(tbl, lo, hi, refit) {
      var wrap = tbl.closest(".pic-cor-wrap") || tbl;
      var leg = wrap.parentNode.querySelector(":scope > .pic-cor-legend");
      if (!leg) {
        leg = document.createElement("div"); leg.className = "pic-cor-legend";
        wrap.parentNode.insertBefore(leg, wrap.nextSibling);
      }
      var d = (hi - lo) < 0.02 ? 4 : 3;
      var f = function (x) { return isFinite(x) ? x.toFixed(d) : "NA"; };
      leg.innerHTML = '<span class="pic-cor-leg-t">Pearson r</span>' +
        '<span class="pic-cor-leg-scale"><span class="pic-cor-leg-bar"></span>' +
        '<span class="pic-cor-leg-ticks"><span>' + f(lo) + '</span><span>' + f((lo + hi) / 2) + '</span><span>' + f(hi) + '</span></span></span>' +
        '<span class="pic-cor-leg-n">' + (refit
          ? "color scale refit to the shown samples (off-diagonal min&ndash;max)"
          : "color scale: all samples (minimum&ndash;1)") + '</span>';
      leg.classList.toggle("pic-cor-leg-refit", !!refit);
    }
    function exApplyCor(tab) {
      var any = exAny();
      tab.querySelectorAll("table.pic-cor").forEach(function (tbl) {
        var cols = exColInfo(tbl); if (!cols.length) return;
        var hdr = tbl.querySelector("thead tr");
        var rows = Array.prototype.slice.call(tbl.querySelectorAll("tbody > tr[data-sample]"));
        var offCol = cols.map(function (c) { return any && !!EX.excl[c.s]; });
        var vis = cols.filter(function (c, j) { return !offCol[j]; });
        // 列の区切り (見えている列の中で group が変わる所)
        var cstart = {}, prevG = null;
        vis.forEach(function (c, j) { if (j > 0 && c.g !== prevG) cstart[c.k] = 1; prevG = c.g; });
        var rstart = {}; prevG = null; var firstR = true;
        rows.forEach(function (r) {
          if (any && EX.excl[r.dataset.sample]) return;
          var g = r.dataset.grow || ""; if (!firstR && g !== prevG) rstart[r.dataset.sample] = 1; prevG = g; firstR = false;
        });
        // 色の幅。全部表示: R と同じ (全セルの最小 .. 1)。除外あり: 表示中の対角以外の最小 .. 最大
        var lo = Infinity, hi = -Infinity;
        rows.forEach(function (r) {
          if (any && EX.excl[r.dataset.sample]) return;
          cols.forEach(function (c, j) {
            if (offCol[j]) return;
            if (any && c.s === r.dataset.sample) return;   // 対角
            var v = parseFloat(r.children[c.k].dataset.v);
            if (isFinite(v)) { if (v < lo) lo = v; if (v > hi) hi = v; }
          });
        });
        if (!any) hi = 1;
        if (!isFinite(lo) || !isFinite(hi)) { lo = 0; hi = 1; }
        exCorLegend(tbl, lo, hi, any);
        cols.forEach(function (c, j) {
          var th = hdr.children[c.k];
          if (th._picCls0 === undefined) th._picCls0 = th.className;
          th.style.display = offCol[j] ? "none" : "";
          th.className = !any ? th._picCls0 : ("pic-cor-ch" + (cstart[c.k] ? " gsep-l" : ""));
        });
        rows.forEach(function (r) {
          var roff = any && !!EX.excl[r.dataset.sample];
          r.style.display = roff ? "none" : "";
          var rh = r.children[0];
          if (rh._picCls0 === undefined) rh._picCls0 = rh.className;
          rh.className = !any ? rh._picCls0 : ("pic-cor-rh" + (rstart[r.dataset.sample] ? " gsep-t" : ""));
          cols.forEach(function (c, j) {
            var td = r.children[c.k]; if (!td) return;
            if (td._picSt0 === undefined) { td._picSt0 = td.getAttribute("style") || ""; td._picCls0 = td.className; }
            td.style.display = offCol[j] ? "none" : "";
            if (!any) { td.setAttribute("style", td._picSt0); td.className = td._picCls0; return; }
            if (offCol[j] || roff) return;
            var v = parseFloat(td.dataset.v);
            var bg = isFinite(v) ? exCorColor(v, lo, hi) : "#ffffff";
            td.setAttribute("style", "background:" + bg + ";color:" + exTextColor(bg));
            td.className = "val" + (cstart[c.k] ? " gsep-l" : "") + (rstart[r.dataset.sample] ? " gsep-t" : "");
          });
        });
      });
    }
    // ---- DEG heatmap: 列を隠す (値と並びは全サンプルのまま) ----
    function exApplyHeatmap(tab) {
      var any = exAny();
      tab.querySelectorAll("table.pic-hm").forEach(function (tbl) {
        if (!any && !tbl._picDirty) return;
        var cols = exColInfo(tbl); if (!cols.length) return;
        var hdr = tbl.querySelector("thead tr");
        var offCol = cols.map(function (c) { return any && !!EX.excl[c.s]; });
        var cstart = {}, prevG = null, j0 = 0;
        cols.forEach(function (c, j) { if (offCol[j]) return; if (j0 > 0 && c.g !== prevG) cstart[c.k] = 1; prevG = c.g; j0++; });
        var orig = {};
        cols.forEach(function (c, j) {
          var th = hdr.children[c.k];
          if (th._picCls0 === undefined) th._picCls0 = th.className;
          orig[c.k] = / gsep\b/.test(" " + th._picCls0);
          th.style.display = offCol[j] ? "none" : "";
          th.className = !any ? th._picCls0 : ("s" + (cstart[c.k] ? " gsep" : ""));
        });
        tbl.querySelectorAll("tbody > tr").forEach(function (r) {
          var ch = r.children;
          cols.forEach(function (c, j) {
            var td = ch[c.k]; if (!td) return;
            var dsp = offCol[j] ? "none" : "";
            if (td.style.display !== dsp) td.style.display = dsp;
            var cl = (any ? cstart[c.k] : orig[c.k]) ? "gsep" : "";
            if (td.className !== cl) td.className = cl;
          });
        });
        tbl._picDirty = any;
      });
    }

    // ---- 比較の選択 (group×group 行列 / ラジオの並び): 丸ごと外した group を含む比較を隠す ----
    function exOutKeys() {
      var outK = {};
      if (exAny()) EX.groups.forEach(function (g) { if (exGroupOut(g)) outK[exGroupKey(g)] = 1; });
      return outK;
    }
    // 選んでいた比較が隠れたら、残っている最初の比較に替える
    function exReselect(checkedOff, firstOn) {
      if (checkedOff && firstOn && checkedOff.type === "radio") {
        firstOn.checked = true;
        firstOn.dispatchEvent(new Event("change", { bubbles: true }));
      }
    }
    function exApplyCmatrix(tab) {
      var outK = exOutKeys();
      tab.querySelectorAll("table.pic-cmatrix").forEach(function (tbl) {
        // 外した group の列 (見出しとセル) と行を隠す。上三角の各比較は、どちらかの group が外れると消える
        var hcells = Array.prototype.slice.call(tbl.querySelectorAll("thead tr > th"));
        var colG = hcells.map(function (th) { return exGroupKey(th.textContent.trim()); });
        hcells.forEach(function (th, k) { th.style.display = outK[colG[k]] ? "none" : ""; });
        var checkedOff = null, firstOn = null;
        tbl.querySelectorAll("tbody tr").forEach(function (tr) {
          var rh = tr.querySelector("th.pic-cmx-rh"), rg = rh ? exGroupKey(rh.textContent.trim()) : "";
          tr.style.display = outK[rg] ? "none" : "";
          Array.prototype.forEach.call(tr.children, function (td, k) {
            if (td.tagName === "TD") td.style.display = outK[colG[k]] ? "none" : "";
            if (!td.classList.contains("pic-cmx-cell")) return;
            var off = !!(outK[rg] || outK[colG[k]]);
            var inp = td.querySelector("input"); if (!inp) return;
            inp.disabled = off;
            if (off && inp.checked) checkedOff = inp;
            if (!off && !firstOn) firstOn = inp;
          });
        });
        exReselect(checkedOff, firstOn);
      });
      // 行列でない比較の並び (ラベルが "A / B" のラジオ)。ORA の method などは対象外
      tab.querySelectorAll(".pic-select-bar").forEach(function (bar) {
        var checkedOff = null, firstOn = null;
        bar.querySelectorAll("label.pic-select-chk").forEach(function (lb) {
          var inp = lb.querySelector("input"); if (!inp) return;
          var ab = lb.textContent.trim().split(" / ");
          var off = ab.length === 2 && !!(outK[exGroupKey(ab[0].trim())] || outK[exGroupKey(ab[1].trim())]);
          lb.style.display = off ? "none" : "";
          inp.disabled = off;
          if (off && inp.checked) checkedOff = inp;
          if (!off && !firstOn) firstOn = inp;
        });
        exReselect(checkedOff, firstOn);
      });
    }

    // ---- group ごとの trace (ORA の cluster expression profiles など): 丸ごと外した group を消す ----
    // legendgroup が group 名の trace だけが対象。category の x 軸からもその group を抜く。
    function exApplyGroupTraces(div) {
      if (!div || !div.data || !div.dataset.rendered) return;
      var kind = exTabKind(div.closest("section.pic-tab") || {});
      if (kind !== "deg" && kind !== "enrich") return;
      var any = exAny();
      if (!any && !div._picGt0) return;
      var gk = {}; EX.groups.forEach(function (g) { gk[exGroupKey(g)] = 1; });
      var idx = [];
      div.data.forEach(function (t, i) { if (t.legendgroup != null && gk[exGroupKey(t.legendgroup)] && !(t.meta && t.meta.pic_sample != null)) idx.push(i); });
      if (!idx.length) return;
      if (!div._picGt0) {
        var ax0 = {};
        Object.keys(div.layout || {}).forEach(function (k) {
          var a = div.layout[k];
          if (/^xaxis\d*$/.test(k) && a && a.categoryarray) ax0[k] = a.categoryarray.slice();
        });
        div._picGt0 = { vis: div.data.map(function (t) { return t.visible === undefined ? true : t.visible; }), ax: ax0 };
      }
      var o = div._picGt0, outK = exOutKeys(), vis = [], upd = {};
      idx.forEach(function (i) { vis.push(any && outK[exGroupKey(div.data[i].legendgroup)] ? false : o.vis[i]); });
      Object.keys(o.ax).forEach(function (k) {
        upd[k + ".categoryarray"] = any ? o.ax[k].filter(function (g) { return !outK[exGroupKey(g)]; }) : o.ax[k];
      });
      try {
        Plotly.restyle(div, { visible: vis }, idx);
        if (Object.keys(upd).length) Plotly.relayout(div, upd);
      } catch (e) {}
      if (!any) div._picGt0 = null;
    }

    // ---- plotly のサンプルごとの trace (Aggregation): 外したサンプルの線を消す ----
    // 凡例は group ごとに先頭の trace が持つので、先頭を消したときは残りの先頭に渡す。
    function exApplyTraces(div) {
      if (!div || !div.data || !div.dataset.rendered) return;
      var has = div.data.some(function (t) { return t.meta && t.meta.pic_sample != null; });
      if (!has) return;
      var any = exAny();
      if (!any && !div._picTr0) return;
      if (!div._picTr0) div._picTr0 = div.data.map(function (t) {
        return { name: t.name, showlegend: t.showlegend, visible: t.visible === undefined ? true : t.visible, g: t.legendgroup };
      });
      var o = div._picTr0, legName = {}, vis = [], nm = [], sl = [], idx = [], seen = {};
      o.forEach(function (x) { if (x.showlegend && x.g != null && legName[x.g] === undefined) legName[x.g] = x.name; });
      div.data.forEach(function (t, i) {
        var s = t.meta && t.meta.pic_sample, off = any && s != null && !!EX.excl[s];
        idx.push(i);
        if (!any) { vis.push(o[i].visible); nm.push(o[i].name); sl.push(o[i].showlegend); return; }
        vis.push(off ? false : o[i].visible);
        var lead = !off && o[i].g != null && !seen[o[i].g];
        if (lead) seen[o[i].g] = 1;
        if (o[i].g != null && legName[o[i].g] !== undefined) {
          sl.push(lead);
          nm.push(lead ? legName[o[i].g] : (o[i].showlegend ? escHtml(s) : o[i].name));
        } else { sl.push(o[i].showlegend); nm.push(o[i].name); }
      });
      try {
        Plotly.restyle(div, { visible: vis, name: nm, showlegend: sl }, idx);
        Plotly.relayout(div, { "yaxis.autorange": true });
      } catch (e) {}
      if (!any) div._picTr0 = null;
    }

    function exApplyTab(tab) {
      if (!tab || !EX.ready) return;
      var any = exAny();
      if (!any && !tab._picExDirty) return;
      if (picRecomputePCA(tab)) { /* PCA のタブ */ }
      tab.querySelectorAll(".pic-plot").forEach(function (d) { exApplyTraces(d); exApplyGroupTraces(d); });
      exApplyRows(tab);
      exApplyCor(tab);
      exApplyHeatmap(tab);
      exApplyCmatrix(tab);
      tab.querySelectorAll(".pic-expr").forEach(function (b) { if (b._picRefilter) b._picRefilter(); });
      exNote(tab);
      tab._picExDirty = any;
    }
    function exApplyAll() { document.querySelectorAll("section.pic-tab").forEach(exApplyTab); }

    // 木を持たない旧来の group トグル (Genome split など): 表示の on/off だけ
    function applyLegacyGroupToggles(tab) {
      var tgl = Array.prototype.filter.call(tab.querySelectorAll(".pic-gtoggle"), function (c) { return !c.closest(".pic-tree"); });
      if (!tgl.length) return;
      var exprs = tab.querySelectorAll(".pic-expr");
      if (exprs.length) { exprs.forEach(function (b) { if (b._picRefilter) b._picRefilter(); }); return; }
      var hide = {};
      tgl.forEach(function (c) { if (!c.checked) hide[c.dataset.group] = 1; });
      tab.querySelectorAll(".pic-plot").forEach(function (div) {
        if (!div.dataset.rendered || !div.data) return;
        var groups = {};
        div.data.forEach(function (t, i) { var g = (t.legendgroup != null) ? t.legendgroup : t.name; if (g != null) (groups[g] = groups[g] || []).push(i); });
        Object.keys(groups).forEach(function (g) {
          try { Plotly.restyle(div, { visible: hide[g] ? "legendonly" : true }, groups[g]); } catch (e) {}
        });
      });
      tab.querySelectorAll("[data-group],[data-grow]").forEach(function (r) {
        if (r.classList.contains("pic-gtoggle")) return;
        r.style.display = hide[r.dataset.group || r.dataset.grow] ? "none" : "";
      });
      tab.querySelectorAll("[data-gcol]").forEach(function (c) {
        c.style.display = hide[c.dataset.gcol] ? "none" : "";
      });
    }
    // ゲノム切替などで呼ぶ: 共通の除外 + 旧来のトグル
    function applyGroupToggles(tab) { exApplyTab(tab); applyLegacyGroupToggles(tab); }
    // summary 内のチェックボックスを押しても折りたたみが切り替わらないようにする
    document.querySelectorAll(".pic-tree-g > summary > .pic-gtoggle").forEach(function (cb) {
      cb.addEventListener("click", function (e) { e.stopPropagation(); });
    });
    // 木のチェック (group = 配下の一括 / sample = 1 つ) は共通の状態を変える。木の外の group トグルは旧来どおり。
    document.querySelectorAll(".pic-gtoggle,.pic-stoggle").forEach(function (cb) {
      cb.addEventListener("change", function () {
        if (cb.closest(".pic-tree")) {
          if (cb.classList.contains("pic-gtoggle")) {
            (EX.members[cb.dataset.group] || []).forEach(function (s) { if (cb.checked) delete EX.excl[s]; else EX.excl[s] = 1; });
          } else if (cb.checked) { delete EX.excl[cb.dataset.sample]; } else { EX.excl[cb.dataset.sample] = 1; }
          exCommit();
          return;
        }
        var tab = cb.closest("section.pic-tab"); if (tab) applyLegacyGroupToggles(tab);
      });
    });
    // HTML テーブルの Download PNG ボタン
    document.querySelectorAll(".pic-png-btn").forEach(function (btn) {
      btn.addEventListener("click", function () {
        var el = document.getElementById(btn.dataset.cap);
        if (el) capturePng(el, (btn.dataset.name || "figure") + ".png", btn);
      });
    });
    // 埋め込みファイルのダウンロード。展開は report_data.js (picFileBlob) が DecompressionStream で行う
    function saveBlob(blob, name) {
      var url = URL.createObjectURL(blob);
      var a = document.createElement("a"); a.href = url; a.download = name;
      document.body.appendChild(a); a.click(); a.remove();
      setTimeout(function () { URL.revokeObjectURL(url); }, 1500);
    }
    function downloadEmbedded(id) {
      picFileBlob(id).then(function (f) { saveBlob(f.blob, f.name); }, function () {
        var r = picFileRaw(id);   // 展開できない環境は zlib のまま保存
        if (r) saveBlob(new Blob([r.bytes]), r.name + ".zz");
      });
    }
    document.querySelectorAll("[data-file]").forEach(function (el) {
      el.addEventListener("click", function (e) { e.preventDefault(); downloadEmbedded(el.dataset.file); });
      el.addEventListener("keydown", function (e) {
        if (e.key === "Enter" || e.key === " ") { e.preventDefault(); downloadEmbedded(el.dataset.file); }
      });
    });
    // PCA の「選択中サンプルだけ」の CSV は埋め込みではなくその場で生成する
    document.querySelectorAll(".pic-pca-dlsel").forEach(function (el) {
      var run = function (e) { e.preventDefault(); picPCADownloadSelected(el.closest("section.pic-tab")); };
      el.addEventListener("click", run);
      el.addEventListener("keydown", function (e) { if (e.key === "Enter" || e.key === " ") run(e); });
    });
    // サンプル除外の状態を読み (URL のハッシュ > localStorage)、あれば全タブに当てる。
    // 以後に描く図には picAfterRender で当てる (PCA はタブごとにまとめて引き直す)。
    exInitModel();
    document.querySelectorAll("section.pic-tab").forEach(exApplyCor);   // 相関の凡例 (除外なしは R の色のまま)
    if (EX.ready) {
      exLoad();
      var pcaTimers = {};
      picAfterRender = function (el) {
        var tab = el.closest && el.closest("section.pic-tab");
        if (!tab || !(exAny() || tab._picExDirty)) return;
        if (typeof PIC_PCA !== "undefined" && PIC_PCA && PIC_PCA[tab.id]) {
          clearTimeout(pcaTimers[tab.id]);
          pcaTimers[tab.id] = setTimeout(function () { picRecomputePCA(tab); }, 0);
          return;
        }
        exApplyTraces(el);
        exApplyGroupTraces(el);
      };
      if (exAny()) exCommit();
    }
  });

  // ウィンドウリサイズで再レイアウト
  var rt;
  window.addEventListener("resize", function () {
    clearTimeout(rt);
    rt = setTimeout(function () {
      document.querySelectorAll(".pic-plot").forEach(function (p) {
        if (p.dataset.rendered) {
          try { Plotly.Plots.resize(p); } catch (e) {}
        }
      });
    }, 200);
  });
})();
