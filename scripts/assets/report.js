// plotly のレンダリングを遅延実行する。
// 多数のプロットを最初に描くと重いので、details が開いた時 / details 外なら即時に描く。
(function () {
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
    try { if (history.replaceState) history.replaceState(null, "", "#" + key); } catch (e) {}
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
  function initExpr(box) {
    var id = box.dataset.expr;
    var data = window.PIC_EXPR && window.PIC_EXPR[id];
    if (!data || box.dataset.wired) return;
    box.dataset.wired = "1";
    var input = box.querySelector(".pic-expr-input");
    var tagbox = box.querySelector(".pic-expr-tagbox");
    var msg = box.querySelector(".pic-expr-msg");
    var plots = box.querySelector(".pic-expr-plots");
    var ac = box.querySelector(".pic-expr-ac");
    var groups = uniqueInOrder(data.groups);
    var selected = [];   // 選択中の遺伝子 index (タグの並び順)

    // lookup (exact) と候補用の配列
    var lut = {};
    var cand = [];   // {disp, low, i}
    for (var i = 0; i < data.ens.length; i++) {
      var ens = data.ens[i] ? String(data.ens[i]) : "";
      var ext = data.ext[i] ? String(data.ext[i]) : "";
      if (ens) lut[ens.toLowerCase()] = i;
      if (ext) lut[ext.toLowerCase()] = i;
      var disp = ext ? (ext + " | " + ens) : ens;
      cand.push({ disp: disp, low: (ext + " " + ens).toLowerCase(), tok: (ext || ens), i: i });
    }

    function dispLabel(gi) { return data.ext[gi] ? (data.ext[gi] + " (" + data.ens[gi] + ")") : data.ens[gi]; }

    // 遺伝子ごとに独立した y 軸 (下端 0、各遺伝子で最適スケール) を持つ小パネルを並べる。
    // 列数はコンテナ幅から決め、幅が変わって列数が変わるときだけ描き直す。
    var PANEL_ROW_H = 260, PANEL_VGAP = 56, PANEL_HGAP = 56;
    var plotState = { idxs: [], ncol: 0 };
    function panelCols(width, n) {
      var minW = Math.max(170, 70 + groups.length * 34);
      return Math.max(1, Math.min(n, Math.floor((width + PANEL_HGAP) / (minW + PANEL_HGAP))));
    }
    function axId(i) { return i === 0 ? "" : String(i + 1); }
    function plotGenes(idxs) {
      plotState.idxs = idxs.slice();
      plots.innerHTML = "";
      plotState.ncol = 0;
      if (!idxs.length) return;
      var width = plots.clientWidth;
      if (!width) return;   // 非表示タブ: 表示されたら ResizeObserver が描く
      var n = idxs.length, ncol = panelCols(width, n), nrow = Math.ceil(n / ncol);
      plotState.ncol = ncol;
      var mar = { l: 78, r: 12, t: 34, b: 50 };
      var H = mar.t + mar.b + nrow * PANEL_ROW_H + (nrow - 1) * PANEL_VGAP;
      var plotW = Math.max(1, width - mar.l - mar.r), plotH = H - mar.t - mar.b;
      var colW = (plotW - (ncol - 1) * PANEL_HGAP) / ncol;
      var div = document.createElement("div");
      div.className = "pic-plot";
      div.style.height = H + "px";
      plots.appendChild(div);
      var layout = {
        showlegend: true, hovermode: false,
        legend: { orientation: "h", yanchor: "bottom", y: 1 + 6 / plotH, x: 0 },
        margin: mar,
        // 縦軸名は左端に 1 つだけ
        annotations: [{ text: "Normalized count", xref: "paper", yref: "paper", x: 0, y: 0.5,
                        xanchor: "center", yanchor: "middle", xshift: -62, textangle: -90,
                        showarrow: false, font: { size: 14 } }]
      };
      var traces = [];
      idxs.forEach(function (gi, k) {
        var r = Math.floor(k / ncol), c = k % ncol, id = axId(k);
        var x0 = c * (colW + PANEL_HGAP) / plotW, x1 = (c * (colW + PANEL_HGAP) + colW) / plotW;
        var yTop = 1 - r * (PANEL_ROW_H + PANEL_VGAP) / plotH, yBot = yTop - PANEL_ROW_H / plotH;
        var gl = data.ext[gi] || data.ens[gi];
        // x は群 (色で区別できるので目盛りは出さない)。遺伝子名は x 軸タイトルに置く
        layout["xaxis" + id] = { domain: [x0, Math.min(1, x1)], anchor: "y" + id, type: "category",
                                 categoryorder: "array", categoryarray: groups, showticklabels: false,
                                 title: { text: escHtml(gl), standoff: 6, font: { size: 13 } },
                                 fixedrange: true, showline: true, linecolor: "#d7dee6", mirror: true };
        layout["yaxis" + id] = { domain: [Math.max(0, yBot), yTop], anchor: "x" + id,
                                 rangemode: "tozero", zeroline: false,
                                 showline: true, linecolor: "#d7dee6", mirror: true };
        groups.forEach(function (g) {
          var ys = [], counts = data.counts[gi];
          for (var s = 0; s < data.samples.length; s++) if (data.groups[s] === g) ys.push(counts[s]);
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
      Plotly.newPlot(div, traces, layout, { responsive: true, displaylogo: false }).then(function () {
        attachTickTips(div, idxs);
        div.on("plotly_afterplot", function () { attachTickTips(div, idxs); });
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
    if (window.ResizeObserver) {
      var _roTimer = null;
      new ResizeObserver(function () {
        clearTimeout(_roTimer);
        _roTimer = setTimeout(function () {
          if (!plotState.idxs.length || !plots.clientWidth) return;
          var nc = panelCols(plots.clientWidth, plotState.idxs.length);
          if (nc !== plotState.ncol) plotGenes(plotState.idxs);
        }, 120);
      }).observe(plots);
    }

    // x 軸の遺伝子ラベルに hover -> contrast ごとの padj (q-value) を表示
    function fmtQ(v) {
      if (v == null || isNaN(v)) return "NA";
      if (v !== 0 && (v < 0.001 || v >= 1000)) return v.toExponential(2);
      return String(v);
    }
    // x 軸の遺伝子ラベルに、確実に hover できる透明オーバーレイを重ねる。
    // (plotly の目盛りテキストは実操作で mouse イベントが届きにくいため)
    function attachTickTips(div, idxs) {
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
      ticks.forEach(function (t, k) {
        var gi = idxs[k];
        if (gi === undefined || !t) return;
        var q = data.qval[gi] || [];
        var sig = q.some(function (v) { return v != null && !isNaN(v) && v < fdr; });
        if (sig) { t.style.fill = "#d7301f"; t.style.fontWeight = "700"; }
        var r = t.getBoundingClientRect();
        var zone = document.createElement("div");
        zone.className = "pic-expr-tickzone";
        zone.style.cssText = "position:absolute;z-index:5;cursor:pointer;left:" +
          (r.left - drect.left - 8) + "px;top:" + (r.top - drect.top - 5) + "px;width:" +
          (r.width + 16) + "px;height:" + (r.height + 10) + "px;";
        (function (gi) {
          function show(e) {
            var qq = data.qval[gi] || [];
            var rows = data.contrasts.map(function (c, ci) {
              var v = qq[ci];
              if (v != null && !isNaN(v) && v >= 1) return "";  // padj=1 (有意でない) の行は非表示
              return "<tr><td>" + escHtml(c) + "</td><td>" + fmtQ(qq[ci]) + "</td></tr>";
            }).join("");
            var head = "<div class='pic-expr-tip-h'>" + escHtml(data.ext[gi] || data.ens[gi]) + "</div>";
            var body = rows
              ? "<table><tr><th>contrast</th><th>padj</th></tr>" + rows + "</table>"
              : "<div class='pic-expr-tip-none'>not significant (padj = 1) in any comparison</div>";
            tip.innerHTML = head + body;
            tip.style.display = "block";
            // 描画後の実寸で viewport 内にクランプ (下端の見切れ防止・必要なら上/左へ反転)
            tip.style.left = "0px"; tip.style.top = "0px";
            var tw = tip.offsetWidth, th = tip.offsetHeight;
            var vw = window.innerWidth, vh = window.innerHeight;
            var left = e.clientX + 14;
            if (left + tw > vw - 6) left = e.clientX - 14 - tw;
            if (left < 6) left = 6;
            var top = e.clientY + 14;
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

    function shortLabel(gi) { return data.ext[gi] || data.ens[gi]; }

    // 選択中の遺伝子タグを描画し、プロットを更新する
    function renderChips() {
      tagbox.querySelectorAll(".pic-expr-chip").forEach(function (c) { c.remove(); });
      selected.forEach(function (gi) {
        var chip = document.createElement("span");
        chip.className = "pic-expr-chip";
        chip.appendChild(document.createTextNode(shortLabel(gi)));
        var b = document.createElement("button");
        b.type = "button"; b.textContent = "×"; b.title = "remove";
        b.addEventListener("click", function () { removeGene(gi); });
        chip.appendChild(b);
        tagbox.insertBefore(chip, input);
      });
    }
    function refresh() {
      renderChips();
      if (selected.length) { msg.textContent = ""; plotGenes(selected); }
      else { plots.innerHTML = ""; msg.textContent = "Add a gene above to see its expression across groups."; }
    }
    function addGene(gi) {
      if (gi === undefined || gi === null || selected.indexOf(gi) !== -1) return;
      selected.push(gi); input.value = ""; hideAc(); refresh(); input.focus();
    }
    function removeGene(gi) {
      selected = selected.filter(function (x) { return x !== gi; });
      refresh(); input.focus();
    }

    // ---- autocomplete (選択済みは除外) ----
    function hideAc() { ac.style.display = "none"; ac.innerHTML = ""; }
    function showAc() {
      var t = (input.value || "").trim().toLowerCase();
      if (t.length < 1) { hideAc(); return; }
      var pref = [], sub = [];
      for (var k = 0; k < cand.length && pref.length < 12; k++) {
        var c = cand[k];
        if (selected.indexOf(c.i) !== -1) continue;
        if (c.low.indexOf(t) === 0 || (c.tok && c.tok.toLowerCase().indexOf(t) === 0)) pref.push(c);
      }
      if (pref.length < 12) {
        for (var k2 = 0; k2 < cand.length && (pref.length + sub.length) < 12; k2++) {
          var c2 = cand[k2];
          if (selected.indexOf(c2.i) !== -1 || pref.indexOf(c2) !== -1) continue;
          if (c2.low.indexOf(t) !== -1) sub.push(c2);
        }
      }
      var list = pref.concat(sub);
      if (!list.length) { hideAc(); return; }
      ac.innerHTML = "";
      list.forEach(function (c) {
        var it = document.createElement("div");
        it.className = "pic-expr-ac-item"; it.textContent = c.disp; it._gi = c.i;
        it.addEventListener("mousedown", function (e) { e.preventDefault(); addGene(c.i); });
        ac.appendChild(it);
      });
      ac.style.display = "block";
    }

    input.addEventListener("input", showAc);
    input.addEventListener("keydown", function (e) {
      if (e.key === "Enter") {
        e.preventDefault();
        var t = (input.value || "").trim();
        if (!t) return;
        var gi = lut[t.toLowerCase()];
        if (gi === undefined) {
          var first = ac.querySelector(".pic-expr-ac-item");
          if (first && first._gi !== undefined) addGene(first._gi);
        } else { addGene(gi); }
      } else if (e.key === "Backspace" && !input.value && selected.length) {
        removeGene(selected[selected.length - 1]);
      } else if (e.key === "Escape") { hideAc(); }
    });
    input.addEventListener("blur", function () { setTimeout(hideAc, 150); });

    // 既定: pvalue 最小の 5 遺伝子
    selected = (data.default || []).map(function (e) { return lut[String(e).toLowerCase()]; })
      .filter(function (x) { return x !== undefined; });
    refresh();
  }

  document.addEventListener("DOMContentLoaded", function () {
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
      var initial = (location.hash || "").replace(/^#/, "");
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
      for (pc = 0; pc < npc; pc++) head.push("PC" + (pc + 1));
      head.push("group");
      var rows = [head.join(",")];
      for (k = 0; k < keep.length; k++) {
        var row = [pd.samples[keep[k]]];
        for (pc = 0; pc < npc; pc++) row.push(res.coords[pc][k]);
        row.push(pd.groups[keep[k]]);
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
    // group トグル: plotly の trace を on/off + HTML テーブルの行/列を隠す
    function applyGroupToggles(tab) {
      // PCA タブだけはトグルの意味が違う (計算対象の選択) ので専用処理に委ねる
      if (picRecomputePCA(tab)) return;
      var hide = {};
      tab.querySelectorAll(".pic-gtoggle").forEach(function (c) { if (!c.checked) hide[c.dataset.group] = 1; });
      // plotly: group ごとに trace の visible を設定
      tab.querySelectorAll(".pic-plot").forEach(function (div) {
        if (!div.dataset.rendered || !div.data) return;
        var groups = {};
        div.data.forEach(function (t, i) { var g = (t.legendgroup != null) ? t.legendgroup : t.name; if (g != null) (groups[g] = groups[g] || []).push(i); });
        Object.keys(groups).forEach(function (g) {
          try { Plotly.restyle(div, { visible: hide[g] ? "legendonly" : true }, groups[g]); } catch (e) {}
        });
      });
      // HTML テーブル/バー: 行 (data-group/data-grow) と列 (data-gcol) を隠す。
      // ただし group トグル自身 (.pic-gtoggle) は data-group を持つが隠さない。
      tab.querySelectorAll("[data-group],[data-grow]").forEach(function (r) {
        if (r.classList.contains("pic-gtoggle")) return;
        r.style.display = hide[r.dataset.group || r.dataset.grow] ? "none" : "";
      });
      tab.querySelectorAll("[data-gcol]").forEach(function (c) {
        c.style.display = hide[c.dataset.gcol] ? "none" : "";
      });
    }
    // group -> sample 階層セレクタ: group のチェックは配下サンプルの一括トグル。
    // 逆にサンプルを個別に触ったら、group 側は配下の状態から checked / indeterminate
    // を決め直す (一部だけ選ばれている状態が見て分かるように)。
    function picTreeApplyGroup(cb, tab) {
      tab.querySelectorAll('.pic-stoggle[data-sgroup="' + cb.dataset.group + '"]').forEach(function (s) {
        s.checked = cb.checked;
      });
    }
    function picTreeSyncGroups(tab) {
      tab.querySelectorAll(".pic-gtoggle").forEach(function (g) {
        var kids = tab.querySelectorAll('.pic-stoggle[data-sgroup="' + g.dataset.group + '"]');
        if (!kids.length) return;
        var on = 0;
        kids.forEach(function (k) { if (k.checked) on++; });
        g.checked = on > 0;
        g.indeterminate = on > 0 && on < kids.length;
      });
    }
    // summary 内のチェックボックスを押しても折りたたみが切り替わらないようにする
    document.querySelectorAll(".pic-tree-g > summary > .pic-gtoggle").forEach(function (cb) {
      cb.addEventListener("click", function (e) { e.stopPropagation(); });
    });
    // group トグルと sample トグル (PCA の計算対象選択) の両方を拾う
    document.querySelectorAll(".pic-gtoggle,.pic-stoggle").forEach(function (cb) {
      cb.addEventListener("change", function () {
        var tab = cb.closest("section.pic-tab"); if (!tab) return;
        if (cb.classList.contains("pic-gtoggle")) picTreeApplyGroup(cb, tab);
        picTreeSyncGroups(tab);
        applyGroupToggles(tab);
      });
    });
    // HTML テーブルの Download PNG ボタン
    document.querySelectorAll(".pic-png-btn").forEach(function (btn) {
      btn.addEventListener("click", function () {
        var el = document.getElementById(btn.dataset.cap);
        if (el) capturePng(el, (btn.dataset.name || "figure") + ".png", btn);
      });
    });
    // 埋め込みファイル (gzip+base64) のダウンロード。DecompressionStream で展開
    function b64ToBytes(b64) {
      var bin = atob(b64), n = bin.length, a = new Uint8Array(n);
      for (var i = 0; i < n; i++) a[i] = bin.charCodeAt(i);
      return a;
    }
    function saveBlob(blob, name) {
      var url = URL.createObjectURL(blob);
      var a = document.createElement("a"); a.href = url; a.download = name;
      document.body.appendChild(a); a.click(); a.remove();
      setTimeout(function () { URL.revokeObjectURL(url); }, 1500);
    }
    function downloadEmbedded(id) {
      var f = (window.PIC_FILES || {})[id]; if (!f) return;
      // R の memCompress(type="gzip") は実際には zlib/deflate 形式 (RFC1950) を出力する
      var bytes = b64ToBytes(f.gz), name = f.n || "download";
      if (typeof DecompressionStream !== "undefined") {
        try {
          var ds = new DecompressionStream("deflate");
          new Response(new Blob([bytes]).stream().pipeThrough(ds)).blob()
            .then(function (b) { saveBlob(b, name); })
            .catch(function () { saveBlob(new Blob([bytes]), name + ".zz"); });
          return;
        } catch (e) {}
      }
      saveBlob(new Blob([bytes]), name + ".zz");  // 非対応環境は zlib のまま保存
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
