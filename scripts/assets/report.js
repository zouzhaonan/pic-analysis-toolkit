// plotly のレンダリングを遅延実行する。
// 多数のプロットを最初に描くと重いので、details が開いた時 / details 外なら即時に描く。
(function () {
  function renderPlot(el) {
    if (!el || el.dataset.rendered) return;
    // 非表示ブロック内は描画しない (表示された時に描く)
    if (el.closest && el.closest(".pic-select-item[hidden]")) return;
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
    root.querySelectorAll("input[type=checkbox]").forEach(function (cb) {
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

    function plotGenes(idxs) {
      plots.innerHTML = "";
      if (!idxs.length) return;
      var div = document.createElement("div");
      div.className = "pic-plot";
      div.style.height = Math.max(360, 300 + idxs.length * 8) + "px";
      plots.appendChild(div);
      var geneLabels = idxs.map(function (gi) { return data.ext[gi] || data.ens[gi]; });
      // group ごとに 1 トレース: x = 遺伝子(大カテゴリ), boxmode=group で group を横並び
      var traces = groups.map(function (g) {
        var xs = [], ys = [];
        idxs.forEach(function (gi) {
          var counts = data.counts[gi], gl = data.ext[gi] || data.ens[gi];
          for (var s = 0; s < data.samples.length; s++) {
            if (data.groups[s] === g) { xs.push(gl); ys.push(counts[s]); }
          }
        });
        var col = data.palette[g] || "#1f77b4";
        return {
          type: "box", name: g, x: xs, y: ys, boxpoints: "all", jitter: 0.5, pointpos: 0,
          marker: { color: col, size: 5 }, line: { color: col }, fillcolor: hexToRgba(col, 0.25)
        };
      });
      // 遺伝子間の縦の区切り線 (カテゴリ境界 = i-0.5)
      var seps = [];
      for (var b = 1; b < idxs.length; b++) {
        seps.push({ type: "line", x0: b - 0.5, x1: b - 0.5, xref: "x", yref: "paper", y0: 0, y1: 1,
                    line: { color: "#d7dee6", width: 1 } });
      }
      Plotly.newPlot(div, traces, {
        yaxis: { title: "Normalized count", rangemode: "tozero" },
        xaxis: { type: "category", categoryorder: "array", categoryarray: geneLabels },
        boxmode: "group", showlegend: true, shapes: seps, hovermode: false,
        legend: { orientation: "h", yanchor: "bottom", y: 1.02, x: 0 },
        margin: { l: 56, r: 16, t: 30, b: 80 }
      }, { responsive: true, displaylogo: false }).then(function () {
        attachTickTips(div, idxs);
        div.on("plotly_afterplot", function () { attachTickTips(div, idxs); });
      });
      div.dataset.rendered = "1";
      // 凡例で group を on/off した後、表示中 group だけで y を再スケール
      div.on("plotly_legendclick", function () {
        setTimeout(function () { try { Plotly.relayout(div, { "yaxis.autorange": true }); } catch (e) {} }, 60);
        return true;
      });
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
      var ticks = div.querySelectorAll(".xaxislayer-above .xtick text");
      if (!ticks.length) ticks = div.querySelectorAll(".xtick text");
      var drect = div.getBoundingClientRect();
      var fdr = (typeof data.fdr === "number") ? data.fdr : 0.1;
      ticks.forEach(function (t, k) {
        var gi = idxs[k];
        if (gi === undefined) return;
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
              return "<tr><td>" + escHtml(c) + "</td><td>" + fmtQ(qq[ci]) + "</td></tr>";
            }).join("");
            tip.innerHTML = "<div class='pic-expr-tip-h'>" + escHtml(data.ext[gi] || data.ens[gi]) +
              "</div><table><tr><th>contrast</th><th>padj</th></tr>" + rows + "</table>";
            tip.style.display = "block";
            tip.style.left = (e.clientX + 14) + "px"; tip.style.top = (e.clientY + 14) + "px";
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
    // details 外・非表示ブロック外のプロットは即描画
    document.querySelectorAll(".pic-plot").forEach(function (p) {
      if (p.closest("details")) return;
      if (p.closest(".pic-select-item[hidden]")) return;
      renderPlot(p);
    });
    // details は開いた時に描画
    document.querySelectorAll("details").forEach(function (d) {
      // 既に open のものは描画
      if (d.open) renderWithin(d);
      d.addEventListener("toggle", function () {
        if (d.open) renderWithin(d);
      });
    });
    // チェックボックス (バー or group 行列) でブロックを表示/非表示 (表示時に遅延描画)
    document.querySelectorAll("input[type=checkbox][data-target]").forEach(function (cb) {
      cb.addEventListener("change", function () {
        var t = document.getElementById(cb.dataset.target);
        if (!t) return;
        if (cb.checked) {
          t.hidden = false;
          renderWithin(t);
          // 既描画のプロットは再表示時にサイズを再計算
          t.querySelectorAll(".pic-plot").forEach(function (p) {
            if (p.dataset.rendered) { try { Plotly.Plots.resize(p); } catch (e) {} }
          });
        } else {
          t.hidden = true;
        }
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
    // source パスのコピーボタン
    document.querySelectorAll(".pic-src-copy").forEach(function (b) {
      b.addEventListener("click", function () {
        var txt = b.dataset.copy || "";
        function done() {
          var orig = "copy path";
          b.textContent = "copied"; b.classList.add("copied");
          setTimeout(function () { b.textContent = orig; b.classList.remove("copied"); }, 1200);
        }
        function fallback() {
          var ta = document.createElement("textarea");
          ta.value = txt; ta.style.position = "fixed"; ta.style.opacity = "0";
          document.body.appendChild(ta); ta.select();
          try { document.execCommand("copy"); } catch (e) {}
          ta.remove(); done();
        }
        if (navigator.clipboard && navigator.clipboard.writeText) {
          navigator.clipboard.writeText(txt).then(done, fallback);
        } else { fallback(); }
      });
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
