/* PrimerForge front-end: fetch helpers, SVG charts and the amplicon diagram.
   Charts are hand-rolled SVG - no chart library, no build step. */
const PF = (() => {
  const SVGNS = 'http://www.w3.org/2000/svg';

  /* ---------------- utilities ---------------- */
  const esc = (s) => String(s ?? '').replace(/[&<>"']/g,
    c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

  function debounce(fn, ms) {
    let t; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); };
  }

  /* On a shared server an expired session answers 401: go to the sign-in page and
     come back here afterwards. */
  function signInAgain(r) {
    if (r.status !== 401) return false;
    const here = location.pathname + location.search;
    location.href = '/login?next=' + encodeURIComponent(here);
    return true;
  }

  async function post(url, body) {
    const r = await fetch(url, {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body)
    });
    const data = await r.json().catch(() => ({}));
    if (signInAgain(r)) throw new Error('Please sign in again.');
    if (!r.ok) throw new Error(data.error || `Request failed (${r.status})`);
    return data;
  }

  async function get(url) {
    const r = await fetch(url);
    if (signInAgain(r)) throw new Error('Please sign in again.');
    if (!r.ok) throw new Error(`Request failed (${r.status})`);
    return r.json();
  }

  let toastTimer;
  function toast(msg) {
    const el = document.getElementById('toast');
    if (!el) return;
    el.textContent = msg; el.classList.add('on');
    clearTimeout(toastTimer);
    toastTimer = setTimeout(() => el.classList.remove('on'), 1800);
  }

  function copy(text, msg = 'Copied') {
    navigator.clipboard.writeText(text).then(() => toast(msg),
      () => toast('Copy failed'));
  }

  function pollProgress(runId, onTick, onDone) {
    const tick = async () => {
      let p;
      try { p = await get('/api/progress/' + runId); }
      catch { return setTimeout(tick, 1200); }
      onTick(p);
      if (p.finished) onDone(p); else setTimeout(tick, 600);
    };
    tick();
  }

  const cssVar = (n) => getComputedStyle(document.body).getPropertyValue(n).trim();
  const fmt = (n) => n >= 1e6 ? (n / 1e6).toFixed(1) + 'M'
    : n >= 1e3 ? (n / 1e3).toFixed(n >= 1e4 ? 0 : 1) + 'k' : String(n);
  const commas = (n) => Number(n).toLocaleString();

  /* ---------------- tooltip ---------------- */
  let tipEl;
  function tip() {
    if (!tipEl) {
      tipEl = document.createElement('div');
      tipEl.className = 'tip';
      document.body.appendChild(tipEl);
    }
    return tipEl;
  }
  function showTip(evt, html) {
    const t = tip();
    t.innerHTML = html; t.classList.add('on');
    const r = t.getBoundingClientRect();
    let x = evt.clientX + 14, y = evt.clientY - r.height - 10;
    if (x + r.width > innerWidth - 8) x = evt.clientX - r.width - 14;
    if (y < 8) y = evt.clientY + 16;
    t.style.left = x + 'px'; t.style.top = y + 'px';
  }
  const hideTip = () => tipEl && tipEl.classList.remove('on');

  /* Integer y-axis ticks, deduplicated so small counts don't print "0,1,1,2". */
  function axisTicks(max, want = 3) {
    const n = Math.max(1, Math.min(want, Math.ceil(max)));
    const out = [];
    for (let i = 0; i <= n; i++) {
      const v = Math.round(max * i / n);
      if (!out.length || out[out.length - 1].v !== v) out.push({ v, f: i / n });
    }
    return out;
  }

  function el(tag, attrs = {}, parent = null) {
    const n = document.createElementNS(SVGNS, tag);
    for (const [k, v] of Object.entries(attrs)) n.setAttribute(k, v);
    if (parent) parent.appendChild(n);
    return n;
  }

  function svgRoot(host, w, h) {
    host.innerHTML = '';
    const s = el('svg', { viewBox: `0 0 ${w} ${h}`, role: 'img' });
    s.style.width = '100%'; s.style.height = 'auto';
    host.appendChild(s);
    return s;
  }

  /* ---------------- charts ---------------- */
  const charts = {};

  /* Vertical bars over time. Single series, so no legend - the title names it. */
  charts.timeBars = function (host, rows, { label = 'runs' } = {}) {
    if (!rows.length) { host.innerHTML = '<div class="empty">No activity yet.</div>'; return; }
    const W = 560, H = 190, ml = 34, mr = 8, mt = 10, mb = 26;
    const iw = W - ml - mr, ih = H - mt - mb;
    const max = Math.max(...rows.map(r => r.n), 1);
    const bw = Math.max(3, Math.min(26, iw / rows.length - 3));
    const step = iw / rows.length;
    const s = svgRoot(host, W, H);

    axisTicks(max).forEach(t => {
      const y = mt + ih - ih * t.f;
      el('line', { x1: ml, x2: W - mr, y1: y, y2: y, stroke: cssVar('--grid'), 'stroke-width': 1 }, s);
      const tx = el('text', {
        x: ml - 7, y: y + 3.5, 'text-anchor': 'end', 'font-size': 10,
        fill: cssVar('--muted'), 'font-variant-numeric': 'tabular-nums'
      }, s);
      tx.textContent = t.v;
    });
    rows.forEach((r, i) => {
      const h = Math.max(2, ih * r.n / max);
      const x = ml + i * step + (step - bw) / 2;
      const y = mt + ih - h;
      const rect = el('rect', {
        x, y, width: bw, height: h, rx: Math.min(4, bw / 2),
        fill: cssVar('--accent')
      }, s);
      rect.style.cursor = 'pointer';
      rect.addEventListener('mousemove', e => showTip(e,
        `<b>${esc(r.day)}</b><br>${r.n} ${label}`));
      rect.addEventListener('mouseleave', hideTip);
    });
    el('line', {
      x1: ml, x2: W - mr, y1: mt + ih, y2: mt + ih,
      stroke: cssVar('--axis'), 'stroke-width': 1
    }, s);
    const ticks = [0, Math.floor(rows.length / 2), rows.length - 1]
      .filter((v, i, a) => a.indexOf(v) === i && v >= 0);
    ticks.forEach(i => {
      const t = el('text', {
        x: ml + i * step + step / 2, y: H - 8, 'text-anchor': 'middle',
        'font-size': 10, fill: cssVar('--muted')
      }, s);
      t.textContent = (rows[i].day || '').slice(5);
    });
  };

  /* Horizontal bars for ranked categories. Direct value labels, no legend. */
  charts.hbars = function (host, rows, { key = 'label', val = 'n', unit = '' } = {}) {
    if (!rows.length) { host.innerHTML = '<div class="empty">Nothing to show yet.</div>'; return; }
    rows = rows.slice(0, 10);
    const rowH = 24, W = 560, H = rows.length * rowH + 8;
    const labelW = 96, valW = 42;
    const iw = W - labelW - valW - 10;
    const max = Math.max(...rows.map(r => r[val]), 1);
    const s = svgRoot(host, W, H);
    rows.forEach((r, i) => {
      const y = i * rowH + 4;
      const lbl = el('text', {
        x: labelW - 8, y: y + rowH / 2 - 1, 'text-anchor': 'end', 'font-size': 11.5,
        fill: cssVar('--ink-2'), 'dominant-baseline': 'middle'
      }, s);
      lbl.textContent = String(r[key] ?? '—').slice(0, 16);
      const w = Math.max(2, iw * r[val] / max);
      const bar = el('rect', {
        x: labelW, y: y + 3, width: w, height: rowH - 10, rx: 4,
        fill: cssVar('--accent')
      }, s);
      bar.style.cursor = 'pointer';
      bar.addEventListener('mousemove', e => showTip(e,
        `<b>${esc(r[key])}</b><br>${commas(r[val])} ${esc(unit)}`));
      bar.addEventListener('mouseleave', hideTip);
      const v = el('text', {
        x: labelW + w + 7, y: y + rowH / 2 - 1, 'font-size': 11.5,
        fill: cssVar('--ink-2'), 'dominant-baseline': 'middle',
        'font-variant-numeric': 'tabular-nums'
      }, s);
      v.textContent = commas(r[val]);
    });
  };

  /* Status breakdown. Fixed status hues, always beside a written label. */
  const STATUS_STYLE = {
    'unique': ['--good', 'Single predicted product'],
    'minor off-target': ['--warning', 'A few extra products'],
    'non-specific': ['--critical', 'Many extra products'],
    'no product': ['--critical', 'Intended product not found'],
    'error': ['--critical', 'Check failed'],
    'unavailable': ['--muted', 'Not checked'],
    'not checked': ['--muted', 'Not checked']
  };

  charts.statusBars = function (host, rows) {
    const total = rows.reduce((a, r) => a + r.n, 0);
    if (!total) { host.innerHTML = '<div class="empty">No primer pairs yet.</div>'; return; }
    host.innerHTML = '';
    rows.slice().sort((a, b) => b.n - a.n).forEach(r => {
      const [varName, desc] = STATUS_STYLE[r.status] || ['--muted', ''];
      const pct = 100 * r.n / total;
      const line = document.createElement('div');
      line.style.cssText = 'margin-bottom:11px';
      line.innerHTML =
        `<div style="display:flex;align-items:center;gap:8px;font-size:12.5px;
                     margin-bottom:5px">
           <span style="width:9px;height:9px;border-radius:2px;flex:none;
                        background:var(${varName})"></span>
           <span style="color:var(--ink)">${esc(r.status)}</span>
           <span style="color:var(--muted);font-size:11.5px">${esc(desc)}</span>
           <span style="margin-left:auto;font-variant-numeric:tabular-nums;
                        color:var(--ink-2)">${commas(r.n)}
             <span style="color:var(--muted)">(${pct.toFixed(0)}%)</span></span>
         </div>
         <div style="height:7px;background:var(--surface-2);border-radius:999px;
                     overflow:hidden">
           <div style="height:100%;width:${pct}%;background:var(${varName});
                       border-radius:999px"></div>
         </div>`;
      host.appendChild(line);
    });
  };

  /* Distribution of a continuous measure. */
  charts.histogram = function (host, values, { bins = 18, unit = '', digits = 1 } = {}) {
    values = values.filter(v => typeof v === 'number' && isFinite(v));
    if (values.length < 2) {
      host.innerHTML = '<div class="empty">Not enough data yet.</div>'; return;
    }
    const min = Math.min(...values), max = Math.max(...values);
    const span = (max - min) || 1;
    const counts = new Array(bins).fill(0);
    values.forEach(v => {
      const i = Math.min(bins - 1, Math.floor((v - min) / span * bins));
      counts[i]++;
    });
    const W = 560, H = 190, ml = 34, mr = 8, mt = 10, mb = 26;
    const iw = W - ml - mr, ih = H - mt - mb;
    const cmax = Math.max(...counts, 1);
    const bw = iw / bins;
    const s = svgRoot(host, W, H);
    axisTicks(cmax).forEach(tk => {
      const y = mt + ih - ih * tk.f;
      el('line', { x1: ml, x2: W - mr, y1: y, y2: y, stroke: cssVar('--grid'), 'stroke-width': 1 }, s);
      const t = el('text', {
        x: ml - 7, y: y + 3.5, 'text-anchor': 'end', 'font-size': 10,
        fill: cssVar('--muted'), 'font-variant-numeric': 'tabular-nums'
      }, s);
      t.textContent = tk.v;
    });
    counts.forEach((c, i) => {
      if (!c) return;
      const h = Math.max(2, ih * c / cmax);
      const r = el('rect', {
        x: ml + i * bw + 1, y: mt + ih - h, width: Math.max(1, bw - 2), height: h,
        rx: Math.min(3, (bw - 2) / 2), fill: cssVar('--accent')
      }, s);
      r.style.cursor = 'pointer';
      const lo = (min + span * i / bins).toFixed(digits);
      const hi = (min + span * (i + 1) / bins).toFixed(digits);
      r.addEventListener('mousemove', e => showTip(e,
        `<b>${lo}–${hi}${esc(unit)}</b><br>${c} primer${c === 1 ? '' : 's'}`));
      r.addEventListener('mouseleave', hideTip);
    });
    el('line', {
      x1: ml, x2: W - mr, y1: mt + ih, y2: mt + ih, stroke: cssVar('--axis'),
      'stroke-width': 1
    }, s);
    [[min, ml, 'start'], [(min + max) / 2, ml + iw / 2, 'middle'], [max, W - mr, 'end']]
      .forEach(([v, x, anchor]) => {
        const t = el('text', {
          x, y: H - 8, 'text-anchor': anchor, 'font-size': 10, fill: cssVar('--muted'),
          'font-variant-numeric': 'tabular-nums'
        }, s);
        t.textContent = v.toFixed(digits) + unit;
      });
  };

  /* ---------------- amplicon diagram ---------------- */
  /* One pair: gene track, amplicon span, both primers and the variant. */
  function amplicon(host, locus, pair, opts = {}) {
    const amp = pair.amplicon;
    const pad = Math.max(60, Math.round(pair.product_size * 0.14));
    const vStart = Math.min(locus.start, amp.start), vEnd = Math.max(locus.end, amp.end);
    const x0 = Math.min(amp.start, vStart) - pad, x1 = Math.max(amp.end, vEnd) + pad;
    const span = x1 - x0;

    const W = 900, H = 132;
    const ml = 10, mr = 10, iw = W - ml - mr;
    const X = (g) => ml + (g - x0) / span * iw;
    const s = svgRoot(host, W, H);
    s.setAttribute('class', 'ampmap');

    const yExon = 20, yAmp = 62, yAxis = 104;

    /* gene / exon track */
    const tr = locus.transcript_region;
    if (tr && tr.exons && tr.exons.length) {
      el('line', {
        x1: ml, x2: W - mr, y1: yExon + 7, y2: yExon + 7,
        stroke: cssVar('--axis'), 'stroke-width': 1.5
      }, s);
      tr.exons.forEach((ex, i) => {
        if (ex.end < x0 || ex.start > x1) return;
        const ex0 = Math.max(ex.start, x0), ex1 = Math.min(ex.end, x1);
        const r = el('rect', {
          x: X(ex0), y: yExon, width: Math.max(2, X(ex1) - X(ex0)), height: 14, rx: 2,
          fill: cssVar('--accent'), opacity: .34
        }, s);
        r.style.cursor = 'pointer';
        r.addEventListener('mousemove', e => showTip(e,
          `<b>Exon ${i + 1} of ${tr.exons.length}</b><br>${locus.chrom}:${commas(ex.start)}–${commas(ex.end)}`));
        r.addEventListener('mouseleave', hideTip);
      });
      const lab = el('text', {
        x: ml, y: yExon - 6, 'font-size': 10.5, fill: cssVar('--muted')
      }, s);
      lab.textContent = `${locus.gene || ''} ${tr.id || ''} · strand ${tr.strand > 0 ? '+' : '−'}`;
    }

    /* amplicon body */
    el('rect', {
      x: X(amp.start), y: yAmp - 8, width: X(amp.end) - X(amp.start), height: 16, rx: 3,
      fill: cssVar('--accent'), opacity: .13
    }, s);

    /* primer arrows */
    const arrow = (gs, ge, dir, o, name) => {
      const xa = X(gs), xb = X(ge);
      const w = Math.max(9, xb - xa), head = Math.min(9, w * 0.55);
      const y = yAmp, hh = 7;
      const pts = dir > 0
        ? `${xa},${y - hh} ${xa + w - head},${y - hh} ${xa + w},${y} ${xa + w - head},${y + hh} ${xa},${y + hh}`
        : `${xa + w},${y - hh} ${xa + head},${y - hh} ${xa},${y} ${xa + head},${y + hh} ${xa + w},${y + hh}`;
      const p = el('polygon', {
        points: pts, fill: cssVar('--accent'),
        stroke: cssVar('--surface'), 'stroke-width': 1.5
      }, s);
      p.style.cursor = 'pointer';
      p.addEventListener('mousemove', e => showTip(e,
        `<b>${name}</b><br><span style="font-family:monospace">${esc(o.seq)}</span>` +
        `<br>${o.len} nt · Tm ${o.tm} °C · GC ${o.gc}%` +
        `<br>${locus.chrom}:${commas(o.start)}–${commas(o.end)}`));
      p.addEventListener('mouseleave', hideTip);
    };
    arrow(pair.left.start, pair.left.end, +1, pair.left, 'Forward primer');
    arrow(pair.right.start, pair.right.end, -1, pair.right, 'Reverse primer');

    /* variant marker */
    const vx = X(locus.start);
    const vw = Math.max(2.5, X(locus.end + 1) - vx);
    el('line', {
      x1: vx + vw / 2, x2: vx + vw / 2, y1: yAmp - 26, y2: yAmp + 14,
      stroke: cssVar('--critical'), 'stroke-width': 1.2, 'stroke-dasharray': '2 2'
    }, s);
    const vm = el('rect', {
      x: vx, y: yAmp - 11, width: vw, height: 22, rx: 1.5, fill: cssVar('--critical')
    }, s);
    vm.style.cursor = 'pointer';
    vm.addEventListener('mousemove', e => showTip(e,
      `<b>${esc(locus.label || locus.input)}</b><br>${locus.chrom}:${commas(locus.start)}` +
      `${locus.end !== locus.start ? '–' + commas(locus.end) : ''}` +
      `<br>${esc(locus.ref)} → ${esc(locus.alt)}`));
    vm.addEventListener('mouseleave', hideTip);
    const vt = el('text', {
      x: vx + vw / 2, y: yAmp - 30, 'text-anchor': 'middle', 'font-size': 10.5,
      fill: cssVar('--crit-ink'), 'font-weight': 600
    }, s);
    vt.textContent = locus.hgvs_c ? locus.hgvs_c.split(':').pop()
      : (locus.hgvs_g || '').split(':').pop();

    /* distance callouts */
    const dline = (xa, xb, text) => {
      if (xb - xa < 34) return;
      const y = yAmp + 26;
      el('line', {
        x1: xa, x2: xb, y1: y, y2: y, stroke: cssVar('--muted'), 'stroke-width': 1,
        'marker-start': '', opacity: .7
      }, s);
      [xa, xb].forEach(x => el('line', {
        x1: x, x2: x, y1: y - 3, y2: y + 3, stroke: cssVar('--muted'), 'stroke-width': 1
      }, s));
      const t = el('text', {
        x: (xa + xb) / 2, y: y - 4, 'text-anchor': 'middle', 'font-size': 10,
        fill: cssVar('--muted')
      }, s);
      t.textContent = text;
    };
    dline(X(pair.left.end), vx, `${commas(pair.dist_left)} bp`);
    dline(vx + vw, X(pair.right.start), `${commas(pair.dist_right)} bp`);

    /* coordinate axis */
    el('line', {
      x1: ml, x2: W - mr, y1: yAxis, y2: yAxis, stroke: cssVar('--axis'), 'stroke-width': 1
    }, s);
    [[x0, 'start'], [Math.round((x0 + x1) / 2), 'middle'], [x1, 'end']].forEach(([g, a]) => {
      el('line', { x1: X(g), x2: X(g), y1: yAxis, y2: yAxis + 4, stroke: cssVar('--axis'), 'stroke-width': 1 }, s);
      const t = el('text', {
        x: X(g), y: yAxis + 16, 'text-anchor': a, 'font-size': 10,
        fill: cssVar('--muted'), 'font-variant-numeric': 'tabular-nums'
      }, s);
      t.textContent = commas(g);
    });
    const scale = el('text', {
      x: W - mr, y: yAxis - 7, 'text-anchor': 'end', 'font-size': 10, fill: cssVar('--muted')
    }, s);
    scale.textContent = `${commas(pair.product_size)} bp product · chr${locus.chrom}`;
  }

  /* ---------------- theme ---------------- */
  function initTheme() {
    const stored = localStorage.getItem('pf-theme');
    if (stored) document.documentElement.setAttribute('data-theme', stored);
    const btn = document.getElementById('theme-toggle');
    if (!btn) return;
    btn.addEventListener('click', () => {
      const cur = document.documentElement.getAttribute('data-theme');
      const isDark = cur ? cur === 'dark'
        : matchMedia('(prefers-color-scheme: dark)').matches;
      const next = isDark ? 'light' : 'dark';
      document.documentElement.setAttribute('data-theme', next);
      try { localStorage.setItem('pf-theme', next); } catch { /* private mode */ }
      document.dispatchEvent(new CustomEvent('pf:theme'));
    });
  }
  document.addEventListener('DOMContentLoaded', initTheme);

  return { esc, debounce, post, get, signInAgain, toast, copy, pollProgress, charts, amplicon,
           fmt, commas, showTip, hideTip, STATUS_STYLE };
})();
