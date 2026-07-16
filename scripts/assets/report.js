// plotly のレンダリングを遅延実行する。
// 多数のプロットを最初に描くと重いので、details が開いた時 / details 外なら即時に描く。
(function () {
  function renderPlot(el) {
    if (!el || el.dataset.rendered) return;
    var spec = window.PIC_PLOTS && window.PIC_PLOTS[el.id];
    if (!spec) return;
    var layout = spec.layout || {};
    layout.autosize = true;
    try {
      Plotly.newPlot(el, spec.data, layout, { responsive: true, displaylogo: false });
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

  document.addEventListener("DOMContentLoaded", function () {
    // details 外のプロットは即描画
    document.querySelectorAll(".pic-plot").forEach(function (p) {
      if (!p.closest("details")) renderPlot(p);
    });
    // details は開いた時に描画
    document.querySelectorAll("details").forEach(function (d) {
      // 既に open のものは描画
      if (d.open) renderWithin(d);
      d.addEventListener("toggle", function () {
        if (d.open) renderWithin(d);
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
