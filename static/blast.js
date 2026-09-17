/* BLAST/BLAT front-end: input form, tickets list, results (karyotype, HSP plot,
   results table, hit pop-ups) and the per-hit sequence views. Mirrors the Ensembl
   BLAST/BLAT tool's behaviour; builds on PF helpers from app.js. */
const BL = (() => {
  const E = PF.esc, C = PF.commas;
  const SVGNS = 'http://www.w3.org/2000/svg';

  const svg = (tag, attrs = {}, parent = null) => {
    const n = document.createElementNS(SVGNS, tag);
    for (const [k, v] of Object.entries(attrs)) n.setAttribute(k, v);
    if (parent) parent.appendChild(n);
    return n;
  };
  const cssVar = (n) => getComputedStyle(document.body).getPropertyValue(n).trim();
  const $ = (sel, root = document) => root.querySelector(sel);
  const $$ = (sel, root = document) => Array.from(root.querySelectorAll(sel));
  const ori = (v) => (v || 1) > 0 ? 'Forward' : 'Reverse';
  const sign = (v) => (v || 1) > 0 ? '+' : '−';
  const fmtE = (e) => e === 0 ? '0.0' : (e < 0.001 || e >= 1e4 ? Number(e).toExponential(1) : String(+Number(e).toPrecision(3)));
  const since = (t) => {
    const s = Math.max(0, Math.round(Date.now() / 1000 - t));
    return s < 60 ? `${s}s` : s < 3600 ? `${Math.floor(s / 60)}m ${s % 60}s` : `${Math.floor(s / 3600)}h ${Math.floor(s % 3600 / 60)}m`;
  };
  const when = (t) => new Date(t * 1000).toLocaleString(undefined,
    { day: '2-digit', month: 'short', year: 'numeric', hour: '2-digit', minute: '2-digit' });

  /* %ID colour ramp: Ensembl's karyotype gradient (gold → dark red). */
  const RAMP = ['#ffd700', '#ffa500', '#d2691e', '#b22222', '#8b0000'];
  function rampColour(pident) {
    const t = Math.max(0, Math.min(1, (pident - 50) / 50));
    const x = t * (RAMP.length - 1), i = Math.min(RAMP.length - 2, Math.floor(x)), f = x - i;
    const a = RAMP[i], b = RAMP[i + 1];
    const mix = (k) => Math.round(parseInt(a.substr(k, 2), 16) * (1 - f) + parseInt(b.substr(k, 2), 16) * f);
    return `rgb(${mix(1)},${mix(3)},${mix(5)})`;
  }

  /* ------------------------------------------------------------------ */
  /* Sequence parsing — same rules as blastconf.parse_sequences          */
  /* ------------------------------------------------------------------ */
  function parseSequences(raw, limits) {
    const text = (raw || '').replace(/\r\n?/g, '\n').trim();
    const out = { sequences: [], invalids: 0, is_rna: false, errors: [] };
    if (!text) return out;
    const chunks = text.split(/(?=>)|\n[ \t]*\n+/).filter(c => c && c.trim());
    const seen = new Set();
    for (const chunk of chunks) {
      let lines = chunk.trim().split('\n').map(l => l.trim());
      if (lines.length === 1 && lines[0].startsWith('>')) {        // newlines lost in a copy/paste
        const words = lines[0].split(/\s+/), seqWords = [];
        let only60 = false;
        while (words.length) {
          const w = words.pop();
          if (/^[A-Za-z*-]+$/.test(w) && (!only60 || w.length === 60)) { seqWords.unshift(w); only60 = true; }
          else { words.push(w); break; }
        }
        if (seqWords.length && words.length) lines = [words.join(' '), ...seqWords];
      }
      let description = '', seq = '';
      lines.forEach((l, j) => {
        if (/^[>;]/.test(l)) { if (j === 0) description = l.slice(1).trim(); return; }
        seq += l.toUpperCase().replace(/[\d\s]+/g, '');
      });
      if (!seq) continue;
      if (!/^[A-Z*]+$/.test(seq)) {
        out.invalids++;
        out.errors.push(`${description || seq.slice(0, 20)}: invalid characters (${[...new Set(seq.replace(/[A-Z*]/g, ''))].join('')})`);
        continue;
      }
      if (seq.length > limits.max_sequence_length) {
        out.invalids++;
        out.errors.push(`${description || seq.slice(0, 20)}: longer than ${C(limits.max_sequence_length)} characters`);
        continue;
      }
      const dna = [...seq].filter(c => 'ACTUGNX'.includes(c)).length;
      const type = 100 * dna / seq.length >= limits.dna_threshold_percent ? 'dna' : 'peptide';
      if (type === 'dna' && seq.includes('U')) { out.is_rna = true; seq = seq.replace(/U/g, 'T'); }
      if (seen.has(seq)) { out.invalids++; continue; }
      seen.add(seq);
      out.sequences.push({ description, sequence: seq, type });
      if (out.sequences.length >= limits.max_num_sequences) break;
    }
    return out;
  }

  const looksLikeId = (t) => /[0-9]/.test(t) && /^[a-z][a-z0-9.\-_]{3,40}$/i.test(t.trim());

  /* ------------------------------------------------------------------ */
  /* Input form                                                          */
  /* ------------------------------------------------------------------ */
  function initForm(data) {
    const O = data.options, L = O.limits;
    const speciesAll = data.species;
    const byName = Object.fromEntries(speciesAll.map(s => [s.name, s]));
    const st = {
      sequences: [], queryTypeUser: null, editing: null,
      species: speciesAll.length ? [speciesAll[0].name] : [],
      searchToolUser: '', configs: {}, configSet: {}
    };
    const form = $('#blast-form');
    const ta = $('#seq-input'), seqList = $('#seq-list'), seqInfo = $('#seq-info');
    const seqMsg = $('#seq-msg'), inputWrap = $('#seq-input-wrap');
    const searchSel = $('#search_type'), setSel = $('#config_set');
    const errBox = $('#form-error');

    /* ----- sequences ----- */
    function message(html, cls = 'info') {
      seqMsg.innerHTML = html ? `<div class="note ${cls}" style="margin-top:8px">${html}</div>` : '';
    }

    function addText(text, replaceIndex = null) {
      const trimmed = text.trim();
      if (!trimmed) return;
      if (!/\s/.test(trimmed) && looksLikeId(trimmed) && !/^[ACGTUN]+$/i.test(trimmed)) {
        return fetchId(trimmed, replaceIndex);
      }
      const parsed = parseSequences(text, L);
      let added = 0, dup = 0;
      parsed.sequences.forEach((s, i) => {
        const clash = st.sequences.findIndex((x, k) => x.sequence === s.sequence && k !== replaceIndex);
        if (clash >= 0) { dup++; return; }
        if (replaceIndex !== null && i === 0) { st.sequences[replaceIndex] = s; added++; return; }
        if (st.sequences.length >= L.max_num_sequences) return;
        st.sequences.push(s); added++;
      });
      const bad = parsed.invalids + dup;
      const msgs = [];
      if (parsed.is_rna) msgs.push('You have added a RNA sequence. However, it has been replaced by its DNA equivalent.');
      if (parsed.errors.length) msgs.push(parsed.errors.slice(0, 3).map(E).join('<br>'));
      message(msgs.join('<br>'), parsed.errors.length ? 'bad' : 'info');
      ta.value = '';
      st.editing = null;
      renderSequences({ added, bad });
      setQueryType();
    }

    async function fetchId(id, replaceIndex) {
      message(`<span class="spin"></span> Fetching sequence for ${E(id)}…`);
      try {
        const r = await PF.post('/api/blast/fetch', { id });
        const text = r.sequences.map(s => `>${s.description}\n${s.sequence}`).join('\n\n');
        addText(text, replaceIndex);
        message(`Fetched ${r.sequences.length} sequence${r.sequences.length > 1 ? 's' : ''} for ${E(id)} from ${E(r.sequences[0].source)}.`);
      } catch (e) {
        message(E(e.message), 'bad');
      }
    }

    function renderSequences(info) {
      seqList.innerHTML = st.sequences.map((s, i) => `
        <div class="seq-card" data-i="${i}">
          <div class="row" style="gap:8px">
            <span class="badge ${s.type === 'dna' ? 'badge-muted' : 'badge-warn'}">${s.type === 'dna' ? 'DNA' : 'Protein'}</span>
            ${s.role ? `<span class="badge badge-good">${s.role === 'forward' ? 'Forward primer' : 'Reverse primer'}</span>` : ''}
            <span class="seq-desc mono">&gt;${E(s.description || `Sequence ${i + 1}`)}</span>
            <span class="small muted">${C(s.sequence.length)} ${s.type === 'dna' ? 'bp' : 'aa'}</span>
            <span class="spacer"></span>
            <button type="button" class="copy" data-act="edit" title="Edit this sequence">edit</button>
            <button type="button" class="copy" data-act="remove" title="Remove this sequence">remove</button>
          </div>
          <div class="seq-preview mono">${E(s.sequence.slice(0, 180))}${s.sequence.length > 180 ? '…' : ''}</div>
        </div>`).join('');
      const n = st.sequences.length, left = L.max_num_sequences - n;
      let msg = '';
      if (info && typeof info === 'object') {
        msg = `${info.added || 'No'} sequence${info.added === 1 ? '' : 's'} added` +
          (info.bad ? ` (${info.bad} invalid or duplicate sequence${info.bad > 1 ? 's' : ''} ignored)` : '') + ', ';
      }
      seqInfo.innerHTML = n
        ? `${msg}${left ? `${left} more sequence${left > 1 ? 's' : ''} can be added` : 'no more sequences can be added'}` +
          (left && inputWrap.hidden ? ` · <a href="#" id="add-more">Add more sequences</a>` : '')
        : '';
      inputWrap.hidden = n > 0 && st.editing === null && !(info && info.showInput);
      $('#seq-cancel').hidden = !n;
      const more = $('#add-more');
      if (more) more.addEventListener('click', (e) => { e.preventDefault(); inputWrap.hidden = false; ta.focus(); renderSequences(); });
      updateTools();
    }

    seqList.addEventListener('click', (e) => {
      const btn = e.target.closest('button[data-act]');
      if (!btn) return;
      const i = +btn.closest('.seq-card').dataset.i;
      if (btn.dataset.act === 'remove') {
        st.sequences.splice(i, 1);
        renderSequences();
        setQueryType();
      } else {
        const s = st.sequences[i];
        st.editing = i;
        ta.value = `>${s.description}\n${s.sequence.match(/.{1,60}/g).join('\n')}`;
        inputWrap.hidden = false;
        $('#seq-add').textContent = 'Save sequence';
        ta.focus();
      }
    });

    $('#seq-add').addEventListener('click', () => {
      const idx = st.editing;
      $('#seq-add').textContent = 'Add sequence(s)';
      if (idx !== null && !ta.value.trim()) { st.sequences.splice(idx, 1); st.editing = null; renderSequences(); return; }
      addText(ta.value, idx);
    });
    $('#seq-cancel').addEventListener('click', () => {
      ta.value = ''; st.editing = null; $('#seq-add').textContent = 'Add sequence(s)';
      inputWrap.hidden = st.sequences.length > 0; renderSequences();
    });
    ta.addEventListener('paste', () => setTimeout(() => {
      if (st.editing === null && ta.value.trim()) addText(ta.value);
    }, 0));
    ta.addEventListener('blur', () => {
      if (st.editing === null && ta.value.trim()) addText(ta.value);
    });
    $('#seq-file').addEventListener('change', async (e) => {
      const f = e.target.files[0];
      if (!f) return;
      const fd = new FormData(); fd.append('file', f);
      message(`<span class="spin"></span> Reading ${E(f.name)}…`);
      const r = await fetch('/api/blast/upload', { method: 'POST', body: fd });
      const j = await r.json().catch(() => ({}));
      if (PF.signInAgain(r)) return;
      e.target.value = '';
      if (!r.ok) return message(E(j.error || 'Upload failed'), 'bad');
      addText(j.text);
    });
    $$('.bf-example').forEach(b => b.addEventListener('click', () => {
      inputWrap.hidden = false;
      addText(b.dataset.seq);
    }));

    /* ----- query type ----- */
    $$('input[name=query_type]').forEach(r => r.addEventListener('change', () => {
      st.queryTypeUser = r.value; updateTools();
    }));
    function setQueryType() {
      if (!st.queryTypeUser && st.sequences.length) {
        const dna = st.sequences.filter(s => s.type === 'dna').length;
        const qt = dna >= st.sequences.length - dna ? 'dna' : 'peptide';
        $(`input[name=query_type][value=${qt}]`).checked = true;
      }
      updateTools();
    }

    /* ----- species ----- */
    function renderSpecies() {
      const box = $('#species-picked');
      box.innerHTML = st.species.length ? st.species.map(n => {
        const s = byName[n] || { display_name: n, assembly: '' };
        return `<span class="species-chip"><b>${E(s.display_name)}</b>
          <span class="muted">${E(s.scientific_name || '')} ${E(s.assembly)}</span>
          ${st.species.length > 1 ? `<button type="button" class="copy" data-rm="${E(n)}" title="Remove">×</button>` : ''}</span>`;
      }).join('') : '<span class="small muted">No species selected.</span>';
      $$('[data-rm]', box).forEach(b => b.addEventListener('click', () => {
        st.species = st.species.filter(x => x !== b.dataset.rm); renderSpecies(); updateSources();
      }));
    }
    const dlg = $('#species-dialog');
    $('#species-change').addEventListener('click', () => {
      const list = $('#species-options');
      list.innerHTML = speciesAll.map(s => `
        <label class="species-opt" data-q="${E((s.display_name + ' ' + s.scientific_name + ' ' + s.name).toLowerCase())}">
          <input type="checkbox" value="${E(s.name)}" ${st.species.includes(s.name) ? 'checked' : ''}>
          <span><b>${E(s.display_name)}</b> <span class="muted">${E(s.scientific_name)} · ${E(s.assembly)}</span>
          ${s.blat ? '<span class="badge badge-good" style="margin-left:6px"><i class="bdot"></i>BLAT</span>' : ''}</span>
        </label>`).join('') || '<div class="small muted">No species installed yet. Add one in Settings.</div>';
      $('#species-filter').value = '';
      dlg.showModal();
    });
    $('#species-filter').addEventListener('input', (e) => {
      const q = e.target.value.toLowerCase();
      $$('.species-opt', dlg).forEach(o => { o.hidden = !o.dataset.q.includes(q); });
    });
    $('#species-done').addEventListener('click', () => {
      const picked = $$('#species-options input:checked').map(i => i.value).slice(0, L.species);
      if (picked.length) st.species = picked;
      dlg.close(); renderSpecies(); updateSources();
    });

    /* ----- data sources ----- */
    const available = (src) => st.species.length && st.species.every(n => byName[n] && byName[n].sources[src]);
    function updateSources() {
      ['dna', 'peptide'].forEach(dt => {
        const sel = $(`#source_${dt}`), cur = sel.value;
        const restricted = O.restrictions[searchSel.value] || [];
        sel.innerHTML = O.sources.filter(s => s.db_type === dt && !restricted.includes(s.value)).map(s =>
          `<option value="${s.value}" ${available(s.value) ? '' : 'disabled'}>${E(s.label)}${available(s.value) ? '' : ' (not installed)'}</option>`).join('');
        const opts = [...sel.options];
        const keep = opts.find(o => o.value === cur && !o.disabled) || opts.find(o => !o.disabled) || opts[0];
        if (keep) keep.selected = true;
      });
      updateTools();
    }
    $$('input[name=db_type]').forEach(r => r.addEventListener('change', updateTools));
    ['#source_dna', '#source_peptide'].forEach(s => $(s).addEventListener('change', (e) => {
      $(`input[name=db_type][value=${s.slice(8)}]`).checked = true; updateTools();
    }));

    /* ----- search tools ----- */
    const blatAvailable = () => st.species.length && st.species.every(n => byName[n] && byName[n].blat);
    const currentQt = () => $('input[name=query_type]:checked').value;
    const currentDb = () => $('input[name=db_type]:checked').value;
    const currentSource = () => $(`#source_${currentDb()}`).value;

    function updateTools(preferred) {
      const qt = currentQt(), dt = currentDb(), src = currentSource();
      const blat = blatAvailable();
      const valid = O.search_types.filter(t =>
        t.query_type === qt && t.db_type === dt && t.sources.includes(src) &&
        !(O.restrictions[t.value] || []).includes(src) &&
        (t.value !== O.blat_value || blat) &&
        (!t.min_length || !st.sequences.length || st.sequences.every(s => s.sequence.length > t.min_length)));
      const before = searchSel.value;
      const want = preferred || st.searchToolUser || (blat ? O.blat_value : '') || before;
      searchSel.innerHTML = valid.map(t => `<option value="${t.value}">${E(t.caption)}</option>`).join('');
      const pick = valid.find(t => t.value === want) || valid.find(t => t.value === before) || valid[0];
      if (pick) searchSel.value = pick.value;
      const hint = [];
      if (!valid.length) hint.push('No search tool can search this query type against the selected database.');
      if (!blat && qt === 'dna' && dt === 'dna' && src === 'LATESTGP') {
        hint.push('BLAT is not available for ' + (st.species.length > 1 ? 'all selected species' : 'this species') + ' — install it in <a href="/settings">Settings</a>.');
      }
      const bl = O.search_types.find(t => t.value === O.blat_value);
      if (blat && st.sequences.some(s => s.sequence.length <= bl.min_length) && qt === 'dna' && src === 'LATESTGP') {
        hint.push(`BLAT needs sequences longer than ${bl.min_length} bases.`);
      }
      $('#search-hint').innerHTML = hint.join(' ');
      $('#run-btn').disabled = !valid.length || !st.sequences.length || !st.species.length || !available(src);
      if (searchSel.value !== before || !setSel.dataset.for || setSel.dataset.for !== searchSel.value) onSearchType();
    }
    searchSel.addEventListener('change', () => {
      st.searchToolUser = searchSel.value;
      const restricted = O.restrictions[searchSel.value] || [];
      if (restricted.length) updateSources(); else onSearchType();
    });

    function onSearchType() {
      const stype = searchSel.value;
      const sets = O.config_sets[stype];
      setSel.dataset.for = stype;
      $('#sens-row').hidden = !sets;
      if (sets) {
        const cur = st.configSet[stype] || 'normal';
        setSel.innerHTML = O.sensitivity_options.filter(o => sets[o.value]).map(o =>
          `<option value="${o.value}" ${o.value === cur ? 'selected' : ''}>${E(o.caption)}</option>`).join('');
      }
      if (!st.configs[stype]) st.configs[stype] = defaultsFor(stype);
      if (sets) applySet(stype, setSel.value, false);
      renderConfigs();
    }
    setSel.addEventListener('change', () => {
      st.configSet[searchSel.value] = setSel.value;
      applySet(searchSel.value, setSel.value, true);
    });

    function defaultsFor(stype) {
      return { ...O.config_defaults.all, ...(O.config_defaults[stype] || {}) };
    }
    function applySet(stype, setName, rerender) {
      const preset = (O.config_sets[stype] || {})[setName] || {};
      const cfg = st.configs[stype];
      for (const [k, v] of Object.entries(preset)) if (k in cfg) cfg[k] = String(v);
      if (rerender) renderConfigs();
    }

    const wordRange = (stype) => {
      const t = O.search_types.find(x => x.value === stype);
      return O.word_size_range[t.program] || [2, 15];
    };

    function renderConfigs() {
      const stype = searchSel.value;
      const cfg = st.configs[stype] || {};
      const host = $('#config-groups');
      if (!stype) { host.innerHTML = ''; return; }
      host.innerHTML = O.config_fields.map(group => {
        const fields = group.fields.filter(f => f.name in cfg);
        if (!fields.length) return '';
        return `<fieldset class="cfg-group"><legend>${E(group.title)}</legend><div class="cfg-grid">
          ${fields.map(f => {
            if (f.type === 'checklist') {
              return `<label class="check cfg-check"><input type="checkbox" data-cfg="${f.name}" ${cfg[f.name] === '1' ? 'checked' : ''}>
                <span>${E(f.label)}${f.help ? ` <span class="muted small">— ${E(f.help)}</span>` : ''}</span></label>`;
            }
            let values = f.values || (f.values_by || {})[cfg[f.depends_on]] || [];
            if (f.name === 'word_size') { const [lo, hi] = wordRange(stype); values = values.filter(v => +v.value >= lo && +v.value <= hi); }
            if (values.length && !values.some(v => v.value === cfg[f.name])) cfg[f.name] = values[0].value;
            return `<div><label for="cfg_${f.name}">${E(f.label)}</label>
              <select id="cfg_${f.name}" data-cfg="${f.name}">${values.map(v =>
                `<option value="${E(v.value)}" ${v.value === cfg[f.name] ? 'selected' : ''}>${E(v.caption)}</option>`).join('')}</select></div>`;
          }).join('')}</div></fieldset>`;
      }).join('');
      $$('[data-cfg]', host).forEach(el => el.addEventListener('change', () => {
        cfg[el.dataset.cfg] = el.type === 'checkbox' ? (el.checked ? '1' : '0') : el.value;
        if (['score', 'matrix'].includes(el.dataset.cfg)) renderConfigs();
      }));
    }

    /* ----- reset, prefill, submit ----- */
    function reset() {
      st.sequences = []; st.queryTypeUser = null; st.searchToolUser = ''; st.configs = {}; st.configSet = {};
      st.editing = null;
      st.species = speciesAll.length ? [speciesAll[0].name] : [];
      $('input[name=query_type][value=dna]').checked = true;
      $('input[name=db_type][value=dna]').checked = true;
      $('#description').value = ''; ta.value = ''; message(''); errBox.innerHTML = '';
      inputWrap.hidden = false;
      renderSpecies(); updateSources(); renderSequences();
    }
    $('#reset-btn').addEventListener('click', reset);

    async function prefill() {
      let pre = null;
      try { pre = JSON.parse(sessionStorage.getItem('pf_blast_prefill') || 'null'); } catch { pre = null; }
      sessionStorage.removeItem('pf_blast_prefill');
      if (data.edit_ticket) {
        try {
          const t = await PF.get('/api/blast/ticket/' + encodeURIComponent(data.edit_ticket));
          const jobs = data.edit_job ? t.jobs.filter(j => j.id === data.edit_job) : t.jobs;
          const seqs = data.edit_job && jobs.length
            ? [{ description: jobs[0].seq_desc || '', sequence: jobs[0].sequence, type: jobs[0].seq_type }]
            : t.sequences;
          const sp = data.edit_job && jobs.length ? [jobs[0].species] : t.species;
          pre = { sequences: seqs, species: sp, query_type: t.query_type, db_type: t.db_type, source: t.source,
                  search_type: t.search_type, config_set: t.config_set, configs: t.configs, description: t.description };
        } catch { message('That ticket no longer exists.', 'bad'); }
      }
      if (!pre) return;
      st.sequences = (pre.sequences || []).map(s => ({ description: s.description || '', sequence: String(s.sequence).toUpperCase(), type: s.type || parseSequences(s.sequence, L).sequences[0]?.type || 'dna', role: s.role || null }));
      if (pre.species) st.species = pre.species.filter(n => byName[n]);
      if (!st.species.length && speciesAll.length) st.species = [speciesAll[0].name];
      if (pre.query_type) { st.queryTypeUser = pre.query_type; $(`input[name=query_type][value=${pre.query_type}]`).checked = true; }
      if (pre.db_type) $(`input[name=db_type][value=${pre.db_type}]`).checked = true;
      // Record the wanted tool and preset before anything selects a tool, so the first
      // selection already uses them.
      if (pre.search_type) st.searchToolUser = pre.search_type;
      if (pre.config_set && pre.search_type) st.configSet[pre.search_type] = pre.config_set;
      renderSpecies(); updateSources();
      if (pre.source) { const o = $(`#source_${pre.db_type || 'dna'} option[value=${pre.source}]`); if (o && !o.disabled) o.selected = true; }
      renderSequences({ added: st.sequences.length });
      setQueryType();
      updateTools(pre.search_type);
      if (pre.search_type && searchSel.value === pre.search_type) {
        // Apply the preset to a fresh set of defaults, then any saved option values on top.
        st.configs[pre.search_type] = defaultsFor(pre.search_type);
        onSearchType();
        if (pre.configs) {
          Object.assign(st.configs[pre.search_type], pre.configs); renderConfigs();
          $('#adv').open = true;
        }
      }
      if (pre.description) $('#description').value = pre.description;
      if (data.edit_ticket) message('Loaded the ticket for editing. Adjust the settings and run it again as a new ticket.');
      else if (pre.message) message(E(pre.message));
    }

    form.addEventListener('submit', async (e) => {
      e.preventDefault();
      errBox.innerHTML = '';
      if (ta.value.trim() && st.editing === null) addText(ta.value);
      const btn = $('#run-btn');
      btn.disabled = true;
      try {
        const res = await PF.post('/api/blast/submit', {
          sequences: st.sequences.map(s => ({ description: s.description, sequence: s.sequence, role: s.role || undefined })),
          species: st.species, query_type: currentQt(), db_type: currentDb(), source: currentSource(),
          search_type: searchSel.value, config_set: O.config_sets[searchSel.value] ? setSel.value : '',
          configs: st.configs[searchSel.value] || {}, description: $('#description').value
        });
        PF.toast('Job submitted');
        tickets.refresh(res.ticket_id);
        const jobs = st.sequences.length * st.species.length;
        $('#run-msg').innerHTML = `Ticket submitted with ${jobs} job${jobs > 1 ? 's' : ''}. <a href="/blast/ticket/${res.ticket_id}">View results</a>`;
      } catch (err) {
        errBox.innerHTML = `<div class="note bad" style="margin-top:12px">${E(err.message)}</div>`;
      } finally {
        btn.disabled = false; updateTools();
      }
    });

    renderSpecies();
    updateSources();
    renderSequences();
    const tickets = initTickets($('#tickets-body'));
    prefill();
  }

  /* ------------------------------------------------------------------ */
  /* Recent tickets                                                      */
  /* ------------------------------------------------------------------ */
  function statusTag(j) {
    const cls = { queued: 'badge-muted', running: 'badge-muted', done: 'badge-good', failed: 'badge-crit', cancelled: 'badge-warn' }[j.status] || 'badge-muted';
    let text = { queued: 'Queued', running: 'Running', done: 'Done', failed: 'Failed', cancelled: 'Cancelled' }[j.status] || j.status;
    if (j.status === 'queued' && j.queue_position) text += ` (#${j.queue_position})`;
    if (j.status === 'running' && j.started_at) text += ` · ${since(j.started_at)}`;
    if (j.status === 'done') {
      text += `: ${j.n_hits ? C(j.n_hits) : 'No'} hit${j.n_hits === 1 ? '' : 's'} found`;
      if (!j.n_hits) return `<span class="badge badge-warn" title="This job is finished, but no hits were found. If you believe that there should be a match to your query sequence please edit the job to adjust the configuration parameters and resubmit the search."><i class="bdot"></i>${E(text)}</span>`;
    }
    const spin = j.status === 'running' ? '<span class="spin" style="width:9px;height:9px;border-width:1.5px"></span>' : '<i class="bdot"></i>';
    return `<span class="badge ${cls}" title="${E(j.message || '')}">${spin}${E(text)}</span>`;
  }

  function initTickets(body) {
    let timer = null, highlight = null;
    async function refresh(hl) {
      if (hl) highlight = hl;
      clearTimeout(timer);
      let tickets;
      try { tickets = await PF.get('/api/blast/tickets'); }
      catch { timer = setTimeout(refresh, 5000); return; }
      $('#tickets-count').textContent = `${tickets.length} ticket${tickets.length === 1 ? '' : 's'}`;
      if (!tickets.length) {
        body.innerHTML = `<tr><td colspan="4"><div class="empty">No tickets yet. Submit a search above and it will appear here.</div></td></tr>`;
      } else {
        body.innerHTML = tickets.map(t => `
          <tr class="${t.id === highlight ? 'row-new' : ''}">
            <td style="white-space:nowrap"><b>${E(t.method)}</b><div class="small muted mono">${E(t.id)}</div>
              <div class="small muted">${E(t.source_label)}</div></td>
            <td class="jobs-cell">${t.jobs.map(j => `
              <div class="job-line">
                <a href="/blast/ticket/${t.id}?job=${j.id}" class="job-name">${E(j.description || j.summary)}</a>
                <span class="small muted">${E(j.summary)}</span>
                ${statusTag(j)}
                <span class="job-acts">
                  ${j.status === 'done' ? `<a class="copy" href="/blast/ticket/${t.id}?job=${j.id}">results</a>` : ''}
                  <a class="copy" href="/blast?edit=${t.id}&job=${j.id}" title="Edit and resubmit this job">edit</a>
                  ${['queued', 'running'].includes(j.status) ? `<button class="copy" data-job="${j.id}" data-act="cancel">cancel</button>` : ''}
                  ${['failed', 'cancelled'].includes(j.status) ? `<button class="copy" data-job="${j.id}" data-act="rerun">rerun</button>` : ''}
                  ${t.jobs.length > 1 ? `<button class="copy" data-job="${j.id}" data-act="delete">delete</button>` : ''}
                </span>
              </div>`).join('')}</td>
            <td class="small muted" style="white-space:nowrap">${when(t.created_at)}</td>
            <td style="white-space:nowrap">
              <a class="btn btn-sm" href="/blast?edit=${t.id}" title="Edit and resubmit the whole ticket">Edit</a>
              <button class="btn btn-sm btn-danger" data-ticket="${t.id}">Delete</button>
            </td>
          </tr>`).join('');
      }
      const active = tickets.some(t => t.jobs.some(j => ['queued', 'running'].includes(j.status)));
      timer = setTimeout(refresh, active ? 2500 : 15000);
    }
    body.addEventListener('click', async (e) => {
      const jb = e.target.closest('[data-job]');
      const tb = e.target.closest('[data-ticket]');
      if (jb) {
        if (jb.dataset.act === 'delete' && !confirm('Delete this job and its results?')) return;
        await PF.post(`/api/blast/job/${jb.dataset.job}/${jb.dataset.act}`, {});
        PF.toast({ cancel: 'Job cancelled', rerun: 'Job queued again', delete: 'Job deleted' }[jb.dataset.act]);
        refresh();
      } else if (tb) {
        if (!confirm('Delete this ticket and all of its results?')) return;
        await PF.post(`/api/blast/ticket/${tb.dataset.ticket}/delete`, {});
        PF.toast('Ticket deleted');
        refresh();
      }
    });
    refresh();
    return { refresh };
  }

  /* ------------------------------------------------------------------ */
  /* Hit pop-up (Ensembl's ZMenu)                                        */
  /* ------------------------------------------------------------------ */
  let popup;
  function closePopup() { if (popup) { popup.remove(); popup = null; } }
  document.addEventListener('keydown', (e) => { if (e.key === 'Escape') closePopup(); });
  document.addEventListener('click', (e) => {
    if (popup && !popup.contains(e.target) && !e.target.closest('[data-hits]')) closePopup();
  });

  function showHits(evt, hits, ctx) {
    closePopup();
    hits = hits.slice().sort((a, b) => b.pident - a.pident);
    const kind = ctx.isBlat ? 'BLAT' : 'BLAST';
    popup = document.createElement('div');
    popup.className = 'zmenu';
    popup.innerHTML = `
      <div class="zmenu-head">${hits.length > 1 ? `${hits.length} ${kind} hits` : `${kind} hit`}
        <button class="copy" style="margin-left:auto" aria-label="Close">×</button></div>
      <div class="zmenu-body">${hits.slice(0, 25).map(h => hitBlock(h, ctx, hits.length > 1)).join('')}
      ${hits.length > 25 ? `<div class="small muted" style="padding:6px 0">…and ${hits.length - 25} more</div>` : ''}</div>`;
    document.body.appendChild(popup);
    popup.querySelector('.zmenu-head button').addEventListener('click', closePopup);
    const r = popup.getBoundingClientRect();
    let x = evt.clientX + 12, y = evt.clientY + 12;
    if (x + r.width > innerWidth - 10) x = Math.max(10, evt.clientX - r.width - 12);
    if (y + r.height > innerHeight - 10) y = Math.max(10, innerHeight - r.height - 10);
    popup.style.left = x + 'px'; popup.style.top = y + 'px';
    hits.forEach(h => ctx.highlight && ctx.highlight(h.idx, true));
  }

  function hitBlock(h, ctx, multiple) {
    const links = hitLinks(h, ctx);
    const rows = [
      ['Genomic bp', h.gid ? `<a href="${links.ensemblLoc}" target="_blank" rel="noopener">${E(h.gid)}:${C(h.gstart)}-${C(h.gend)}</a>` : '–'],
      ['Query bp', `${E(h.qid)}:${C(h.qstart)}-${C(h.qend)}`]
    ];
    if (!ctx.genomic) rows.push(['Target', links.target ? `<a href="${links.target}" target="_blank" rel="noopener">${E(h.tid)}</a>` : E(h.tid)]);
    const genes = geneLinks(h, ctx);
    if (genes) rows.push([ctx.genomic ? 'Overlapping Gene(s)' : 'Gene hit', genes]);
    rows.push(['Score', h.score], ['E-value', fmtE(h.evalue)], ['%ID', h.pident], ['Length', C(h.len)]);
    return `<div class="zmenu-hit">
      ${multiple ? `<div class="zmenu-cap">${h.gid ? `${E(h.gid)}:${C(h.gstart)}-${C(h.gend)}` : E(h.tid)}</div>` : ''}
      <table>${rows.map(([k, v]) => `<tr><th>${k}</th><td>${v}</td></tr>`).join('')}</table>
      <div class="zmenu-links"><a href="${links.alignment}">Alignment</a><a href="${links.query}">Query sequence</a>
        ${h.gid ? `<a href="${links.genomic}">Genomic sequence</a>` : ''}</div>
    </div>`;
  }

  function hitLinks(h, ctx) {
    const base = `/blast/job/${ctx.jobId}/hit/${h.idx}`;
    const out = { alignment: `${base}?view=alignment`, query: `${base}?view=query`, genomic: `${base}?view=genomic` };
    if (h.gid) {
      const len = h.gend - h.gstart;
      const s = Math.max(1, Math.floor(h.gstart - len * 0.05)), e = Math.floor(h.gend + len * 0.05);
      out.ensemblLoc = `${ctx.ensemblUrl}/Location/View?r=${encodeURIComponent(h.gid)}:${s}-${e}`;
    }
    if (!ctx.genomic) {
      const t = ctx.source === 'PEP_ALL' ? (h.transcript_id || h.tid) : h.tid;
      out.target = `${ctx.ensemblUrl}/Transcript/${ctx.source === 'PEP_ALL' ? 'ProteinSummary' : 'Summary'}?t=${encodeURIComponent(t.split('.')[0])}`;
    }
    return out;
  }

  function geneLinks(h, ctx) {
    const list = ctx.genomic ? (h.genes || []) : (h.gene_hit ? [h.gene_hit] : []);
    return list.map(g => `<a href="${ctx.ensemblUrl}/Gene/Summary?g=${encodeURIComponent(g.id)}" target="_blank" rel="noopener" title="${E(g.biotype || '')}">${E(g.name || g.id)}</a>`).join(', ');
  }

  /* ------------------------------------------------------------------ */
  /* Karyotype: HSP distribution on genome                               */
  /* ------------------------------------------------------------------ */
  function karyotype(host, chroms, hits, ctx) {
    const byChrom = Object.fromEntries(chroms.map(c => [c.name, c]));
    const placed = hits.filter(h => h.gid && byChrom[h.gid]);
    const unplaced = hits.filter(h => h.gid && !byChrom[h.gid]).length;
    const n = chroms.length;
    const colW = Math.max(30, Math.min(60, 1100 / n));
    const W = Math.max(640, n * colW + 20), H = 330, top = 16, maxH = 250, barW = Math.min(12, colW * 0.34);
    const maxLen = Math.max(...chroms.map(c => c.length));
    host.innerHTML = '';
    const s = svg('svg', { viewBox: `0 0 ${W} ${H}`, class: 'karyo', role: 'img', 'aria-label': 'Hits on chromosomes' });
    host.appendChild(s);
    const ink = cssVar('--ink'), surface = cssVar('--surface'), axis = cssVar('--axis'), muted = cssVar('--muted');
    const stainOpacity = { gneg: 0, gpos25: 0.22, gpos33: 0.3, gpos50: 0.42, gpos66: 0.55, gpos75: 0.62, gpos100: 0.82, gpos: 0.6, gvar: 0.14, stalk: 0.1, acen: 0 };
    const pointers = [];
    const best = hits[0];

    chroms.forEach((c, ci) => {
      const x = 10 + ci * colW + (colW - barW) / 2 - colW * 0.12;
      const h = Math.max(6, maxH * c.length / maxLen);
      const Y = (pos) => top + h * (pos - 1) / c.length;
      const g = svg('g', {}, s);
      const clipId = `clip-${ci}-${Math.random().toString(36).slice(2, 7)}`;
      const clip = svg('clipPath', { id: clipId }, g);
      svg('rect', { x, y: top, width: barW, height: h, rx: barW / 2 }, clip);
      const body = svg('g', { 'clip-path': `url(#${clipId})` }, g);
      svg('rect', { x, y: top, width: barW, height: h, fill: surface }, body);
      (c.bands || []).forEach(b => {
        const y1 = Y(b.start), y2 = Y(b.end + 1);
        if (b.stain === 'acen') {
          const mid = (y1 + y2) / 2, isP = (b.id || '').startsWith('p');
          const pts = isP ? `${x},${y1} ${x + barW},${y1} ${x + barW / 2 + 1},${y2} ${x + barW / 2 - 1},${y2}`
            : `${x + barW / 2 - 1},${y1} ${x + barW / 2 + 1},${y1} ${x + barW},${y2} ${x},${y2}`;
          svg('polygon', { points: pts, fill: muted, opacity: 0.55 }, body);
          return mid;
        }
        const op = stainOpacity[b.stain] ?? 0.3;
        if (op > 0) svg('rect', { x, y: y1, width: barW, height: Math.max(0.5, y2 - y1), fill: ink, 'fill-opacity': op }, body);
      });
      const outline = svg('rect', { x, y: top, width: barW, height: h, rx: barW / 2, fill: 'none', stroke: axis, 'stroke-width': 1 }, g);
      outline.style.cursor = 'pointer';
      const label = svg('text', { x: x + barW / 2, y: top + maxH + 18, 'text-anchor': 'middle', 'font-size': 11, fill: cssVar('--ink-2') }, s);
      label.textContent = c.name;
      [outline, label].forEach(el => {
        el.addEventListener('mousemove', (e) => PF.showTip(e, `<b>Chromosome ${E(c.name)}</b><br>${C(c.length)} bp<br>Click to open in Ensembl`));
        el.addEventListener('mouseleave', PF.hideTip);
        el.addEventListener('click', () => window.open(`${ctx.ensemblUrl}/Location/Chromosome?r=${encodeURIComponent(c.name)}`, '_blank', 'noopener'));
      });

      // Pointers, merged when they would overlap on screen.
      const mine = placed.filter(hh => hh.gid === c.name).sort((a, b) => a.gstart - b.gstart);
      const groups = [];
      mine.forEach(hh => {
        const y = Y((hh.gstart + hh.gend) / 2);
        const last = groups[groups.length - 1];
        if (last && Math.abs(last.y - y) < 5) { last.hits.push(hh); last.y = (last.y + y) / 2; }
        else groups.push({ y, hits: [hh] });
      });
      groups.forEach(gr => {
        const top1 = gr.hits.reduce((a, b) => (b.pident > a.pident ? b : a));
        const px = x + barW + 2, ay = gr.y;
        const arrow = svg('polygon', {
          points: `${px},${ay} ${px + 8},${ay - 4.5} ${px + 8},${ay - 1.5} ${px + 15},${ay - 1.5} ${px + 15},${ay + 1.5} ${px + 8},${ay + 1.5} ${px + 8},${ay + 4.5}`,
          fill: rampColour(top1.pident), stroke: surface, 'stroke-width': 0.6
        }, s);
        arrow.dataset.hits = gr.hits.map(hh => hh.idx).join(',');
        arrow.style.cursor = 'pointer';
        arrow.addEventListener('mousemove', (e) => PF.showTip(e, gr.hits.length > 1
          ? `<b>${gr.hits.length} hits</b> on ${E(c.name)}<br>best %ID ${top1.pident}`
          : `<b>${E(c.name)}:${C(top1.gstart)}-${C(top1.gend)}</b><br>%ID ${top1.pident} · E ${fmtE(top1.evalue)}`));
        arrow.addEventListener('mouseleave', PF.hideTip);
        arrow.addEventListener('click', (e) => { e.stopPropagation(); PF.hideTip(); showHits(e, gr.hits, ctx); });
        pointers.push({ el: arrow, idx: gr.hits.map(hh => hh.idx) });
        if (best && gr.hits.some(hh => hh.idx === best.idx)) {
          svg('rect', { x: px - 2, y: ay - 7, width: 19, height: 14, rx: 2, fill: 'none', stroke: ink, 'stroke-width': 1.2, 'pointer-events': 'none' }, s);
        }
      });
    });

    // %ID legend
    const lg = svg('g', { transform: `translate(${W - 230}, ${H - 26})` }, s);
    const lt = svg('text', { x: 0, y: 10, 'font-size': 10.5, fill: cssVar('--ink-2') }, lg);
    lt.textContent = '%ID';
    for (let i = 0; i < 20; i++) svg('rect', { x: 28 + i * 7, y: 2, width: 7, height: 10, fill: rampColour(50 + i * 50 / 19) }, lg);
    [['≤50', 28], ['100', 168]].forEach(([t, xx]) => {
      const tt = svg('text', { x: xx, y: 24, 'font-size': 9.5, fill: muted, 'text-anchor': 'middle' }, lg); tt.textContent = t;
    });
    const bx = svg('g', { transform: `translate(10, ${H - 24})` }, s);
    svg('rect', { x: 0, y: 0, width: 19, height: 14, rx: 2, fill: 'none', stroke: ink, 'stroke-width': 1.2 }, bx);
    const bt = svg('text', { x: 26, y: 11, 'font-size': 10.5, fill: cssVar('--ink-2') }, bx);
    bt.textContent = 'Best hit (highest score)';

    return {
      unplaced,
      highlight(idx, on) {
        pointers.forEach(p => { if (p.idx.includes(idx)) p.el.classList.toggle('ptr-on', on); });
      }
    };
  }

  /* ------------------------------------------------------------------ */
  /* HSP distribution on query sequence                                  */
  /* ------------------------------------------------------------------ */
  function queryPlot(host, hits, qlen, ctx, limit = 60) {
    const shown = hits.slice(0, limit);
    const W = 900, ml = 8, mr = 8, iw = W - ml - mr;
    const rows = [];
    shown.forEach(h => {
      let r = rows.findIndex(row => row.every(o => h.qstart > o.qend + 2 || h.qend < o.qstart - 2));
      if (r < 0) { rows.push([]); r = rows.length - 1; }
      rows[r].push(h);
    });
    const rowH = 11, top = 56, H = top + rows.length * rowH + 22;
    host.innerHTML = '';
    const s = svg('svg', { viewBox: `0 0 ${W} ${H}`, class: 'qplot', role: 'img', 'aria-label': 'HSPs along the query' });
    host.appendChild(s);
    const X = (p) => ml + iw * (p - 1) / Math.max(1, qlen);
    const ink = cssVar('--ink'), muted = cssVar('--muted'), axis = cssVar('--axis');

    // ruler
    const step = [1, 2, 5].map(m => [1, 10, 100, 1000, 10000, 100000].map(p => m * p)).flat()
      .sort((a, b) => a - b).find(v => qlen / v <= 10) || Math.ceil(qlen / 10);
    svg('line', { x1: ml, x2: W - mr, y1: 20, y2: 20, stroke: axis }, s);
    for (let p = 0; p <= qlen; p += step) {
      const xx = X(Math.max(1, p));
      svg('line', { x1: xx, x2: xx, y1: 16, y2: 20, stroke: axis }, s);
      const anchor = p === 0 ? 'start' : xx > W - mr - 24 ? 'end' : 'middle';
      const t = svg('text', { x: xx, y: 12, 'font-size': 9.5, fill: muted, 'text-anchor': anchor }, s);
      t.textContent = C(Math.max(1, p));
    }
    // query: a chain of black and white boxes, with red where the query also hits elsewhere
    const blocks = 10;
    for (let i = 0; i < blocks; i++) {
      const a = 1 + Math.floor(qlen * i / blocks), b = Math.floor(qlen * (i + 1) / blocks);
      svg('rect', { x: X(a), y: 28, width: Math.max(1, X(b) - X(a)), height: 10, fill: i % 2 ? cssVar('--surface') : ink, stroke: ink, 'stroke-width': 0.6 }, s);
    }
    const best = hits[0];
    hits.slice(1).forEach(h => {
      if (best && h.gid === best.gid && h.gstart <= best.gend && h.gend >= best.gstart) return;
      svg('rect', { x: X(h.qstart), y: 28, width: Math.max(1, X(h.qend) - X(h.qstart)), height: 10, fill: cssVar('--critical'), 'fill-opacity': 0.18, 'pointer-events': 'none' }, s);
    });
    const qlab = svg('text', { x: ml, y: 50, 'font-size': 10, fill: muted }, s);
    qlab.textContent = `Query (${C(qlen)} ${ctx.queryType === 'peptide' ? 'aa' : 'bp'}) · HSPs sorted by score`;

    rows.forEach((row, r) => row.forEach(h => {
      const y = top + r * rowH;
      const rect = svg('rect', { x: X(h.qstart), y, width: Math.max(2, X(h.qend) - X(h.qstart)), height: rowH - 3, rx: 1.5, fill: rampColour(h.pident) }, s);
      rect.dataset.hits = String(h.idx);
      rect.style.cursor = 'pointer';
      rect.addEventListener('mousemove', (e) => PF.showTip(e,
        `<b>${h.gid ? `${E(h.gid)}:${C(h.gstart)}-${C(h.gend)}` : E(h.tid)}</b><br>query ${C(h.qstart)}–${C(h.qend)} · %ID ${h.pident}<br>score ${h.score} · E ${fmtE(h.evalue)}`));
      rect.addEventListener('mouseleave', PF.hideTip);
      rect.addEventListener('click', (e) => { e.stopPropagation(); PF.hideTip(); showHits(e, [h], ctx); });
    }));
    if (hits.length > limit) {
      const t = svg('text', { x: ml, y: H - 6, 'font-size': 10, fill: muted }, s);
      t.textContent = `Showing the ${limit} best of ${C(hits.length)} HSPs`;
    }
  }

  /* ------------------------------------------------------------------ */
  /* Results table                                                       */
  /* ------------------------------------------------------------------ */
  function resultsTable(host, hits, ctx, onHover, opts = {}) {
    const genomic = ctx.genomic;
    const primer = !!opts.primer;
    const primerCol = {
      key: 'tp', title: "3′ end", sort: h => (h.three_prime_ok ? 1 : 0) + (h.full_match ? 1 : 0),
      text: h => h.full_match ? 'full match' : h.three_prime_ok ? '3′ end anneals' : '',
      html: h => h.full_match
        ? `<span class="badge badge-crit" title="The whole primer matches: it will prime here"><i class="bdot"></i>full match</span>`
        : h.three_prime_ok
          ? `<span class="badge badge-warn" title="The last 5 bases pair, so extension is possible"><i class="bdot"></i>anneals</span>`
          : '<span class="muted" title="The 3′ end does not pair, so this site is unlikely to prime">–</span>'
    };
    const links = (h) => hitLinks(h, ctx);
    const locCell = (h) => h.gid
      ? `<a href="${links(h).ensemblLoc}" target="_blank" rel="noopener" title="Region in Ensembl">${E(h.gid)}:${C(h.gstart)}-${C(h.gend)}</a>${h.gapprox ? ' <span class="muted" title="Transcript span; exact mapping unavailable">~</span>' : ''}
         <a class="small" href="${links(h).genomic}" title="View genomic sequence">[Sequence]</a>`
      : '<span class="muted">–</span>';
    const cols = [
      ...(genomic ? [
        { key: 'loc', title: 'Genomic Location', sort: h => [h.gid, h.gstart], text: h => h.gid ? `${h.gid}:${h.gstart}-${h.gend}` : '', html: locCell },
        { key: 'gene', title: 'Overlapping Gene(s)', sort: h => (h.genes || []).map(g => g.name).join(','), text: h => (h.genes || []).map(g => g.name).join(', '), html: h => h.genes ? (geneLinks(h, ctx) || '<span class="muted">–</span>') : '<span class="muted genes-pending">…</span>' },
        { key: 'gori', title: 'Orientation', sort: h => h.gori, text: h => ori(h.gori), html: h => ori(h.gori) }
      ] : [
        { key: 'tid', title: 'Subject name', sort: h => h.tid, text: h => h.tid, html: h => `<a href="${links(h).target}" target="_blank" rel="noopener">${E(h.tid)}</a>` },
        { key: 'gene', title: 'Gene hit', sort: h => (h.gene_hit || {}).name || '', text: h => (h.gene_hit || {}).name || '', html: h => geneLinks(h, ctx) || '<span class="muted">–</span>' },
        {
          key: 'ids', title: 'MANE / accession',
          sort: h => [(h.ids || {}).mane ? 0 : 1, ((h.ids || {}).mane_refseq || (h.ids || {}).refseq_mrna || 'zz')],
          text: h => [(h.ids || {}).mane, (h.ids || {}).mane_refseq || (h.ids || {}).refseq_mrna, (h.ids || {}).ensembl_transcript].filter(Boolean).join(' '),
          html: h => {
            const i = h.ids || {};
            const acc = i.refseq_mrna || i.mane_refseq || i.ensembl_transcript;
            const badge = i.mane === 'MANE Select' ? '<span class="badge badge-good"><i class="bdot"></i>MANE Select</span>'
              : i.mane ? `<span class="badge badge-warn"><i class="bdot"></i>${E(i.mane)}</span>` : '';
            const link = acc ? `<a href="/transcript/${ctx.species}/${encodeURIComponent(acc)}" title="Open the transcript">${E(acc)}</a>` : '';
            return [badge, link].filter(Boolean).join(' ') || '<span class="muted">–</span>';
          }
        },
        { key: 'tstart', title: 'Subject start', num: true, sort: h => h.tstart, text: h => h.tstart, html: h => C(h.tstart) },
        { key: 'tend', title: 'Subject end', num: true, sort: h => h.tend, text: h => h.tend, html: h => C(h.tend) },
        { key: 'tori', title: 'Subject ori', sort: h => h.tori, text: h => ori(h.tori), html: h => ori(h.tori) },
        { key: 'loc', title: 'Genomic Location', sort: h => [h.gid || '', h.gstart || 0], text: h => h.gid ? `${h.gid}:${h.gstart}-${h.gend}` : '', html: locCell },
        { key: 'gori', title: 'Orientation', sort: h => h.gori || 0, text: h => h.gori ? ori(h.gori) : '', html: h => h.gori ? ori(h.gori) : '–' }
      ]),
      ...(primer ? [primerCol] : []),
      { key: 'qid', title: 'Query name', hidden: true, sort: h => h.qid, text: h => h.qid, html: h => E(h.qid) },
      { key: 'qstart', title: 'Query start', num: true, sort: h => h.qstart, text: h => h.qstart, html: h => C(h.qstart) },
      { key: 'qend', title: 'Query end', num: true, sort: h => h.qend, text: h => h.qend, html: h => C(h.qend) },
      { key: 'qori', title: 'Query ori', hidden: true, sort: h => h.qori, text: h => ori(h.qori), html: h => ori(h.qori) },
      { key: 'len', title: 'Length', num: true, sort: h => h.len, text: h => h.len, html: h => `${C(h.len)} <a class="small" href="${links(h).query}" title="View query sequence">[Sequence]</a>` },
      { key: 'score', title: 'Score', num: true, sort: h => h.score, text: h => h.score, html: h => h.score },
      { key: 'evalue', title: 'E-val', num: true, sort: h => h.evalue, text: h => fmtE(h.evalue), html: h => fmtE(h.evalue) },
      { key: 'pident', title: '%ID', num: true, sort: h => h.pident, text: h => h.pident, html: h => `${h.pident} <a class="small" href="${links(h).alignment}" title="View alignment">[Alignment]</a>` }
    ];
    const storeKey = `pf-blast-cols-${genomic ? 'g' : 't'}`;
    let hidden;
    try { hidden = new Set(JSON.parse(localStorage.getItem(storeKey) || 'null') || cols.filter(c => c.hidden).map(c => c.key)); }
    catch { hidden = new Set(cols.filter(c => c.hidden).map(c => c.key)); }
    if (primer) hidden.delete('tp');
    let sortKey = 'score', dir = -1, page = 0, per = 25, filter = '', annealOnly = false;

    host.innerHTML = `
      <div class="toolbar rt-toolbar">
        <input type="search" class="rt-filter" placeholder="Filter results…" aria-label="Filter results">
        ${primer ? `<label class="check small"><input type="checkbox" class="rt-anneal"> Only sites where the 3′ end anneals</label>` : ''}
        <label class="small dim">Show <select class="rt-per" style="width:auto;padding:5px 8px">
          ${[10, 25, 50, 100, 0].map(v => `<option value="${v}" ${v === per ? 'selected' : ''}>${v || 'All'}</option>`).join('')}</select> entries</label>
        <details class="rt-cols"><summary class="btn btn-sm">Columns</summary><div class="rt-cols-menu">
          ${cols.map(c => `<label class="check"><input type="checkbox" value="${c.key}" ${hidden.has(c.key) ? '' : 'checked'}> ${E(c.title)}</label>`).join('')}
        </div></details>
        <span class="spacer"></span>
        <span class="small muted rt-count"></span>
      </div>
      <div class="tbl-wrap"><table class="tbl rt"><thead></thead><tbody></tbody></table></div>
      <div class="row rt-pager" style="margin-top:10px"></div>`;
    const thead = $('thead', host), tbody = $('tbody', host);

    function filtered() {
      const pool = annealOnly ? hits.filter(h => h.three_prime_ok) : hits;
      if (!filter) return pool;
      const q = filter.toLowerCase();
      return pool.filter(h => cols.some(c => String(c.text(h)).toLowerCase().includes(q)));
    }
    function sorted(list) {
      const col = cols.find(c => c.key === sortKey);
      return list.slice().sort((a, b) => {
        let x = col.sort(a), y = col.sort(b);
        if (Array.isArray(x)) { if (x[0] !== y[0]) return (x[0] > y[0] ? 1 : -1) * dir; x = x[1]; y = y[1]; }
        if (x === y) return a.idx - b.idx;
        return (x > y ? 1 : -1) * dir;
      });
    }
    function render() {
      const vis = cols.filter(c => !hidden.has(c.key));
      thead.innerHTML = `<tr>${vis.map(c => `<th class="${c.num ? 'num' : ''}" data-k="${c.key}" style="cursor:pointer">${E(c.title)}${c.key === sortKey ? (dir > 0 ? ' ▲' : ' ▼') : ''}</th>`).join('')}</tr>`;
      const list = sorted(filtered());
      const pages = per ? Math.max(1, Math.ceil(list.length / per)) : 1;
      page = Math.min(page, pages - 1);
      const slice = per ? list.slice(page * per, page * per + per) : list;
      tbody.innerHTML = slice.length ? slice.map(h => `<tr data-idx="${h.idx}">${vis.map(c => `<td class="${c.num ? 'num' : ''}">${c.html(h)}</td>`).join('')}</tr>`).join('')
        : `<tr><td colspan="${vis.length}"><div class="empty">No matching results.</div></td></tr>`;
      const from = list.length ? page * (per || list.length) + 1 : 0;
      $('.rt-count', host).textContent = `Showing ${C(from)} to ${C(Math.min(list.length, per ? (page + 1) * per : list.length))} of ${C(list.length)} entries${filter || annealOnly ? ` (filtered from ${C(hits.length)})` : ''}`;
      $('.rt-pager', host).innerHTML = pages > 1 ? `
        <button class="btn btn-sm" data-p="prev" ${page ? '' : 'disabled'}>Previous</button>
        <span class="small muted">Page ${page + 1} of ${pages}</span>
        <button class="btn btn-sm" data-p="next" ${page < pages - 1 ? '' : 'disabled'}>Next</button>` : '';
      if (genomic) loadGenes(slice);
    }
    const pendingGenes = new Set();
    async function loadGenes(slice) {
      const want = slice.filter(h => h.gid && h.genes == null && !pendingGenes.has(h.idx)).map(h => h.idx);
      if (!want.length) return;
      want.forEach(i => pendingGenes.add(i));
      try {
        const got = await PF.post(`/api/blast/job/${ctx.jobId}/genes`, { idx: want });
        want.forEach(i => { const h = hits.find(x => x.idx === i); if (h && h.genes == null) h.genes = got[String(i)] || []; });
        render();
      } catch { /* genes are decoration */ }
    }

    thead.addEventListener('click', (e) => {
      const th = e.target.closest('th[data-k]');
      if (!th) return;
      if (sortKey === th.dataset.k) dir = -dir; else { sortKey = th.dataset.k; dir = ['score', 'pident', 'len'].includes(sortKey) ? -1 : 1; }
      render();
    });
    $('.rt-filter', host).addEventListener('input', PF.debounce((e) => { filter = e.target.value.trim(); page = 0; render(); }, 150));
    $('.rt-per', host).addEventListener('change', (e) => { per = +e.target.value; page = 0; render(); });
    if (primer) $('.rt-anneal', host).addEventListener('change', (e) => { annealOnly = e.target.checked; page = 0; render(); });
    $('.rt-cols-menu', host).addEventListener('change', (e) => {
      if (e.target.checked) hidden.delete(e.target.value); else hidden.add(e.target.value);
      try { localStorage.setItem(storeKey, JSON.stringify([...hidden])); } catch { /* private mode */ }
      render();
    });
    $('.rt-pager', host).addEventListener('click', (e) => {
      const b = e.target.closest('[data-p]');
      if (!b) return;
      page += b.dataset.p === 'next' ? 1 : -1; render();
    });
    tbody.addEventListener('mouseover', (e) => { const tr = e.target.closest('tr[data-idx]'); if (tr && onHover) onHover(+tr.dataset.idx, true); });
    tbody.addEventListener('mouseout', (e) => { const tr = e.target.closest('tr[data-idx]'); if (tr && onHover) onHover(+tr.dataset.idx, false); });
    render();
    return {
      focus(idx) {
        const list = sorted(filtered());
        const pos = list.findIndex(h => h.idx === idx);
        if (pos < 0) return;
        if (per) page = Math.floor(pos / per);
        render();
        const tr = $(`tr[data-idx="${idx}"]`, tbody);
        if (tr) { tr.classList.add('row-new'); tr.scrollIntoView({ block: 'nearest', behavior: 'smooth' }); setTimeout(() => tr.classList.remove('row-new'), 1800); }
      }
    };
  }

  /* ------------------------------------------------------------------ */
  /* Results page                                                        */
  /* ------------------------------------------------------------------ */
  /* A results page shows one job, or a primer pair as one division per primer. */
  async function initResults(data) {
    $$('[data-job-act]').forEach(b => b.addEventListener('click', async () => {
      const act = b.dataset.jobAct, id = b.dataset.jobId;
      if (act === 'delete' && !confirm('Delete this job and its results?')) return;
      const r = await PF.post(`/api/blast/job/${id}/${act}`, {});
      if (act === 'delete') window.location = r.ticket_left ? `/blast/ticket/${r.ticket_left}` : '/blast';
      else window.location.reload();
    }));
    $$('.copy-command').forEach(b => b.addEventListener('click', () =>
      PF.copy(b.parentElement.querySelector('.command-text').textContent, 'Command copied')));

    data.sections.forEach(sec => {
      const root = document.getElementById(`job-${sec.job.id}`);
      if (!root) return;
      if (['queued', 'running'].includes(sec.job.status)) watchJob(sec.job, root);
      else if (sec.job.status === 'done' && sec.job.n_hits) initSection(sec, root, data);
    });
  }

  function watchJob(job, root) {
    const tick = async () => {
      let j;
      try { j = await PF.get(`/api/blast/job/${job.id}`); } catch { return setTimeout(tick, 3000); }
      if (!['queued', 'running'].includes(j.status)) return window.location.reload();
      const el = $('.job-live', root);
      if (el) el.textContent = j.status === 'running'
        ? `Running for ${since(j.started_at)}…`
        : `Queued${j.queue_position ? ` — position ${j.queue_position} in the queue` : ''}…`;
      setTimeout(tick, 2000);
    };
    tick();
  }

  async function initSection(sec, root, data) {
    const ctx = {
      jobId: sec.job.id, genomic: data.genomic, isBlat: data.is_blat, source: data.ticket.source,
      ensemblUrl: data.species.ensembl_url, queryType: data.ticket.query_type,
      species: data.species.name
    };
    const hits = await PF.get(`/api/blast/job/${sec.job.id}/hits`);
    const primer = !!sec.role;
    let kar = null;
    const highlight = (idx, on) => {
      if (kar) kar.highlight(idx, on);
      $$(`.qplot [data-hits="${idx}"]`, root).forEach(el => el.classList.toggle('ptr-on', on));
    };

    // For a primer, only sites where its 3' end pairs can prime; plot those when there are any.
    const binding = primer ? hits.filter(h => h.three_prime_ok) : hits;
    if (primer) {
      const full = hits.filter(h => h.full_match).length;
      const summary = $('.primer-summary', root);
      if (summary) summary.innerHTML =
        `<b>${C(binding.length)}</b> of ${C(hits.length)} alignment${hits.length === 1 ? '' : 's'} reach the primer's 3′ end with its last 5 bases paired` +
        ` · <b>${C(full)}</b> full-length perfect match${full === 1 ? '' : 'es'}` +
        (full > 1 ? ` <span class="badge badge-warn"><i class="bdot"></i>binds perfectly at ${C(full)} sites</span>` : '');
    }

    const karHost = $('.karyotype', root);
    if (karHost) {
      const plotted = primer && binding.length ? binding : hits;
      const draw = () => {
        if (data.species.karyotype.length && plotted.some(h => h.gid)) {
          kar = karyotype(karHost, data.species.karyotype, plotted, ctx);
          const notes = [];
          if (plotted !== hits) notes.push(`Showing the ${C(plotted.length)} sites where the 3′ end anneals.`);
          if (kar.unplaced) notes.push(`${C(kar.unplaced)} hit${kar.unplaced > 1 ? 's' : ''} on scaffolds or patches outside the chromosomes ${kar.unplaced > 1 ? 'are' : 'is'} listed in the table only.`);
          $('.karyo-note', root).textContent = notes.join(' ');
          if (plotted !== hits) $('.karyo-count', root).textContent = `${C(plotted.length)} binding site${plotted.length === 1 ? '' : 's'}`;
        } else {
          $('.karyo-section', root).hidden = true;
        }
      };
      draw();
      document.addEventListener('pf:theme', draw);
    }
    const qHost = $('.queryplot', root);
    if (qHost) {
      const qlen = sec.seq_len || hits[0].qlen || Math.max(...hits.map(h => h.qend));
      const draw = () => queryPlot(qHost, hits, qlen, ctx);
      const qSection = $('.qplot-section', root);
      qSection.addEventListener('toggle', () => { if (qSection.open) draw(); });
      document.addEventListener('pf:theme', () => { if (qSection.open) draw(); });
    }
    const table = resultsTable($('.results-table', root), hits, ctx, highlight, { primer });
    // Clicking a pointer shows its pop-up and brings the matching row into view.
    ctx.highlight = (idx) => table.focus(idx);
  }

  /* ------------------------------------------------------------------ */
  /* Hit page (Alignment / Query sequence / Genomic sequence)            */
  /* ------------------------------------------------------------------ */
  function initHit(data) {
    const form = $('#view-config');
    if (form) {
      const save = (params) => {
        document.cookie = `pf_view_${data.view}=${encodeURIComponent(JSON.stringify(params))}; path=/; max-age=31536000; SameSite=Lax`;
      };
      form.addEventListener('submit', (e) => {
        e.preventDefault();
        const params = Object.fromEntries(new FormData(form).entries());
        save(params);
        const q = new URLSearchParams({ view: data.view, ...params });
        window.location = `${window.location.pathname}?${q}`;
      });
      $('#view-reset').addEventListener('click', () => {
        // Reset the options but stay on the current tab, if the page has tabs.
        const tab = form.querySelector('input[name=tab]');
        if (tab) save({ tab: tab.value });
        else document.cookie = `pf_view_${data.view}=; path=/; max-age=0`;
        window.location = `${window.location.pathname}?view=${data.view}${tab ? `&tab=${tab.value}` : ''}`;
      });
      $$('select', form).forEach(s => s.addEventListener('change', () => form.requestSubmit()));
    }
    const blastBtn = $('#blast-this');
    if (blastBtn) blastBtn.addEventListener('click', () => {
      sessionStorage.setItem('pf_blast_prefill', JSON.stringify({
        sequences: [{ description: blastBtn.dataset.name, sequence: data.sequence }],
        species: [data.species]
      }));
      window.location = '/blast';
    });
    const copySeq = $('#copy-seq');
    if (copySeq) copySeq.addEventListener('click', () => PF.copy(data.sequence, 'Sequence copied'));

    // Annotation Ensembl REST had not served yet: fetch it in the background, then
    // re-render from the cache. Ensembl can take a minute for variant-dense regions.
    if (data.pending && data.pending.length) {
      PF.post(`/api/blast/annotation/${data.species}/fetch`, { items: data.pending })
        .then(r => {
          if (r.missing.length < data.pending.length) window.location.reload();
          else $('#pending-note').innerHTML = 'Ensembl did not return this annotation (it may be offline or busy). Reload the page to try again.';
        })
        .catch(() => { $('#pending-note').textContent = 'Could not reach Ensembl for this annotation.'; });
    }
  }

  /* Transcript page: same option panel behaviour as the hit views, plus copy/BLAST. */
  function initTranscript(data) {
    initHit(data);
    // Remember the chosen tab along with the other saved options.
    $$('a.tab[data-tab]').forEach(a => a.addEventListener('click', () => {
      let saved = {};
      try {
        const raw = (document.cookie.match(/(?:^|; )pf_view_transcript=([^;]*)/) || [])[1];
        saved = raw ? JSON.parse(decodeURIComponent(raw)) : {};
      } catch { saved = {}; }
      saved.tab = a.dataset.tab;
      document.cookie = `pf_view_transcript=${encodeURIComponent(JSON.stringify(saved))}; path=/; max-age=31536000; SameSite=Lax`;
    }));
    const blast = $('#blast-cdna');
    if (blast) blast.addEventListener('click', () => blastSequence(data.name, data.sequence, data.species));
    const copyC = $('#copy-cdna');
    if (copyC) copyC.addEventListener('click', () => PF.copy(data.sequence, `Copied ${PF.commas(data.sequence.length)} bases`));
    const copyP = $('#copy-protein');
    if (copyP) copyP.addEventListener('click', () => PF.copy(data.protein, `Copied ${PF.commas(data.protein.length)} residues`));
    const copyG = $('#copy-genomic');
    if (copyG) copyG.addEventListener('click', async () => {
      const r = await fetch(copyG.dataset.url);
      if (PF.signInAgain(r)) return;
      if (!r.ok) return PF.toast('Could not load the genomic sequence');
      const seq = (await r.text()).split(/\r?\n/).filter(l => !l.startsWith('>')).join('');
      PF.copy(seq, `Copied ${PF.commas(seq.length)} bases of genomic sequence`);
    });
  }

  /* "BLAST this sequence" from anywhere (e.g. a primer design amplicon). */
  function blastSequence(name, sequence, speciesName) {
    sessionStorage.setItem('pf_blast_prefill', JSON.stringify({
      sequences: [{ description: name, sequence }], species: speciesName ? [speciesName] : undefined
    }));
    window.location = '/blast';
  }

  /* Search both primers of a pair: one job each, shown as a forward and a reverse division.
     Primers are too short for BLAT (Ensembl requires more than 26 bases), so this uses
     BLASTN with Ensembl's "Short sequences" sensitivity. */
  function blastPrimers(pairName, forward, reverse, speciesName) {
    sessionStorage.setItem('pf_blast_prefill', JSON.stringify({
      sequences: [
        { description: `${pairName} forward primer`, sequence: forward.toUpperCase(), role: 'forward' },
        { description: `${pairName} reverse primer`, sequence: reverse.toUpperCase(), role: 'reverse' }
      ],
      species: speciesName ? [speciesName] : undefined,
      query_type: 'dna', db_type: 'dna', source: 'LATESTGP',
      search_type: 'NCBIBLAST_BLASTN', config_set: 'near_oligo',
      message: 'Primer pair loaded: each primer runs as its own BLASTN job with the "Short sequences" settings, and the results show a forward and a reverse primer division.'
    }));
    window.location = '/blast';
  }

  return { parseSequences, initForm, initResults, initHit, initTranscript, blastSequence,
           blastPrimers, rampColour };
})();
