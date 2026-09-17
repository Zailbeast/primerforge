"""PrimerForge - a local web app for designing primers from HGVS variants."""
from __future__ import annotations

import csv
import datetime
import io
import json
import os
import re
import subprocess
import sys
import threading
import webbrowser
from pathlib import Path

from flask import (Flask, Response, abort, g, jsonify, redirect, render_template,
                   request, send_file, session, url_for)

from primerforge import (annotation, auth, blastconf, blastjobs, blatserver, config, engines,
                         pipeline, presets, refseq, seqfetch, seqviews, species, store,
                         transcripts, transcriptview, variants, wsl)
from primerforge.variants import VariantError

app = Flask(__name__)
app.config["JSON_SORT_KEYS"] = False
app.secret_key = auth.secret_key()
app.config.update(SESSION_COOKIE_NAME="primerforge_session", SESSION_COOKIE_HTTPONLY=True,
                  SESSION_COOKIE_SAMESITE="Lax", SESSION_COOKIE_SECURE=config.SECURE_COOKIES,
                  PERMANENT_SESSION_LIFETIME=datetime.timedelta(days=14))
ROOT = Path(__file__).resolve().parent


@app.template_filter("ts")
def _ts(value: float | None) -> str:
    import datetime as dt
    if not value:
        return "-"
    return dt.datetime.fromtimestamp(value).strftime("%d %b %Y  %H:%M")


@app.template_filter("dur")
def _dur(run: dict) -> str:
    if not run.get("finished_at") or not run.get("created_at"):
        return "-"
    s = run["finished_at"] - run["created_at"]
    return f"{s:.1f}s" if s < 60 else f"{s / 60:.1f} min"


def _env() -> dict:
    return {
        "genome_ready": config.genome_ready(),
        "blast_ready": config.blast_ready(),
        "seq_source": refseq.source(),
        "assembly": config.ASSEMBLY,
        "blat_species": sum(1 for s in species.all_species() if s["blat"]),
    }


# --------------------------------------------------------------------------
# Pages
# --------------------------------------------------------------------------
@app.route("/")
def page_design():
    return render_template("design.html", presets=presets.PRESETS,
                           recent=store.list_runs(limit=8), env=_env())


@app.route("/runs")
def page_runs():
    return render_template("runs.html", env=_env())


@app.route("/run/<run_id>")
def page_run(run_id: str):
    run = store.get_run(run_id)
    if not run:
        abort(404)
    return render_template("run.html", run=run, env=_env())


@app.route("/settings")
def page_settings():
    return render_template("settings.html", env=_env(), native=wsl.NATIVE,
                           components=species.COMPONENTS, install_order=species.INSTALL_ORDER,
                           paths={"data": str(config.DATA_DIR),
                                  "genome": str(config.GENOME_FASTA),
                                  "blastdb": str(config.BLAST_DB)})


# --------------------------------------------------------------------------
# Design API
# --------------------------------------------------------------------------
@app.post("/api/validate")
def api_validate():
    """Offline syntax check so the input box can give instant feedback."""
    items = variants.parse_input_block(request.json.get("text", ""))
    out = []
    for item in items:
        try:
            kind, cleaned = variants.classify(item)
            out.append({"input": item, "ok": True, "kind": kind, "cleaned": cleaned})
        except VariantError as exc:
            out.append({"input": item, "ok": False, "error": str(exc)})
    return jsonify({"items": out, "n_ok": sum(1 for i in out if i["ok"]),
                    "n_bad": sum(1 for i in out if not i["ok"])})


@app.post("/api/design")
def api_design():
    data = request.json or {}
    try:
        run_id = pipeline.start_run(
            input_text=data.get("variants", ""),
            assay=data.get("assay", "sanger"),
            overrides=data.get("params") or {},
            label=data.get("label", ""),
            assembly=data.get("assembly", "GRCh38"),
            check_specificity=data.get("check_specificity", True),
        )
    except (VariantError, ValueError) as exc:
        return jsonify({"error": str(exc)}), 400
    store.set_run_owner(run_id, _username())
    return jsonify({"run_id": run_id})


@app.get("/api/progress/<run_id>")
def api_progress(run_id: str):
    prog = pipeline.get_progress(run_id)
    if prog:
        return jsonify(prog)
    run = store.get_run(run_id)
    if not run:
        abort(404)
    return jsonify({"run_id": run_id, "finished": run["status"] != "running",
                    "stage": run["status"], "done": run["n_ok"] + run["n_failed"],
                    "total": run["n_input"], "pct": 100.0})


@app.get("/api/run/<run_id>")
def api_run(run_id: str):
    run = store.get_run(run_id)
    if not run:
        abort(404)
    return jsonify(run)


@app.get("/api/runs")
def api_runs():
    return jsonify(store.list_runs(limit=int(request.args.get("limit", 500)),
                                   search=request.args.get("q", ""),
                                   assay=request.args.get("assay", "")))


@app.get("/api/stats")
def api_stats():
    return jsonify(store.stats())


@app.post("/api/run/<run_id>/delete")
def api_delete(run_id: str):
    run = store.get_run(run_id)
    if not run:
        abort(404)
    if not _may_change(run.get("owner")):
        return _not_yours("run")
    store.delete_run(run_id)
    return jsonify({"ok": True})


# --------------------------------------------------------------------------
# Exports
# --------------------------------------------------------------------------
@app.get("/run/<run_id>/export.<fmt>")
def export(run_id: str, fmt: str):
    run = store.get_run(run_id)
    if not run:
        abort(404)
    stamp = f"primerforge_{run_id}"

    if fmt == "json":
        return Response(json.dumps(run, indent=2), mimetype="application/json",
                        headers={"Content-Disposition":
                                 f"attachment; filename={stamp}.json"})

    if fmt == "fasta":
        lines = []
        for v in run["variants"]:
            if v["status"] != "ok":
                continue
            tag = (v["locus"] or {}).get("gene") or v["chrom"]
            for p in v["pairs"]:
                base = f"{tag}_{v['pos']}_pair{p['rank'] + 1}"
                lines += [f">{base}_F", p["left"]["seq"],
                          f">{base}_R", p["right"]["seq"]]
        return Response("\n".join(lines) + "\n", mimetype="text/plain",
                        headers={"Content-Disposition":
                                 f"attachment; filename={stamp}.fasta"})

    if fmt == "csv":
        buf = io.StringIO()
        w = csv.writer(buf, lineterminator="\n")
        w.writerow(["variant", "gene", "hgvs_c", "hgvs_p", "hgvs_g", "chrom", "pos",
                    "pair", "primer", "sequence", "length", "tm_c", "gc_pct",
                    "amplicon_bp", "amplicon_region", "variant_offset_in_amplicon",
                    "pair_penalty", "specificity", "warnings"])
        for v in run["variants"]:
            if v["status"] != "ok":
                w.writerow([v["input_hgvs"], "", "", "", "", "", "", "", "FAILED",
                            v["error"] or "", "", "", "", "", "", "", "", "", ""])
                continue
            loc = v["locus"] or {}
            for p in v["pairs"]:
                spec = (p.get("specificity") or {}).get("status", "")
                amp = p["amplicon"]
                region = f"{amp['chrom']}:{amp['start']}-{amp['end']}"
                for side, key in (("F", "left"), ("R", "right")):
                    o = p[key]
                    w.writerow([
                        v["input_hgvs"], loc.get("gene", ""), loc.get("hgvs_c", ""),
                        loc.get("hgvs_p", ""), loc.get("hgvs_g", ""), loc.get("chrom", ""),
                        loc.get("start", ""), p["rank"] + 1, side, o["seq"], o["len"],
                        o["tm"], o["gc"], p["product_size"], region,
                        p["variant_pos_in_amplicon"], p["penalty"], spec,
                        "; ".join(p["warnings"])])
        return Response(buf.getvalue(), mimetype="text/csv",
                        headers={"Content-Disposition":
                                 f"attachment; filename={stamp}.csv"})
    abort(404)


# --------------------------------------------------------------------------
# Local genome setup
# --------------------------------------------------------------------------
_setup_lock = threading.Lock()
_setup_proc: subprocess.Popen | None = None


def _pid_running(status: dict) -> bool:
    """Whether an unfinished install recorded in a status file is still running, e.g.
    one the Linux installer started as its own service rather than from this page."""
    pid = status.get("pid")
    if os.name == "nt" or not pid or status.get("done"):
        return False              # (os.kill on Windows would terminate the process)
    try:
        os.kill(int(pid), 0)
    except (OSError, ValueError):
        return False
    return True


@app.get("/api/setup/status")
def api_setup_status():
    sys.path.insert(0, str(ROOT))
    from tools.setup_genome import read_status
    status = read_status()
    running = (_setup_proc is not None and _setup_proc.poll() is None) or _pid_running(status)
    return jsonify({**status, "running": running, **_env()})


@app.post("/api/setup/start")
def api_setup_start():
    global _setup_proc
    sys.path.insert(0, str(ROOT))
    from tools.setup_genome import read_status
    with _setup_lock:
        if (_setup_proc is not None and _setup_proc.poll() is None) or _pid_running(read_status()):
            return jsonify({"ok": True, "already_running": True})
        log = open(config.DATA_DIR / "setup.log", "a", encoding="utf-8")
        _setup_proc = subprocess.Popen(
            [sys.executable, str(ROOT / "tools" / "setup_genome.py")],
            stdout=log, stderr=subprocess.STDOUT, cwd=str(ROOT),
            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
    return jsonify({"ok": True})


@app.post("/api/setup/reload")
def api_setup_reload():
    refseq.reset_index()
    return jsonify(_env())


# --------------------------------------------------------------------------
# BLAST/BLAT pages
# --------------------------------------------------------------------------
def _species_payload() -> list[dict]:
    out = []
    for s in species.all_species():
        m = species.load(s["name"]) or {}
        out.append({**s, "ensembl_url": f"{species.ensembl_site(m)}/{species.web_name(s['name'])}"})
    return out


@app.route("/blast")
def page_blast():
    return render_template("blast.html", env=_env(), species_list=_species_payload(),
                           form_options=blastconf.form_options(),
                           edit_ticket=request.args.get("edit", ""),
                           edit_job=request.args.get("job", ""))


def _ticket_or_404(ticket_id: str) -> dict:
    t = blastjobs.get_ticket(ticket_id)
    if not t:
        abort(404)
    return t


PRIMER_LABELS = {"forward": "Forward primer", "reverse": "Reverse primer"}


@app.route("/blast/ticket/<ticket_id>")
def page_blast_ticket(ticket_id: str):
    t = _ticket_or_404(ticket_id)
    if not t["jobs"]:
        abort(404)
    jid = request.args.get("job")
    job = next((j for j in t["jobs"] if j["id"] == jid), t["jobs"][0])

    # A forward and a reverse primer searched together are shown as one primer pair:
    # a division per primer, forward first, rather than one job at a time.
    same_species = [j for j in t["jobs"] if j["species"] == job["species"]]
    roles = {j.get("role") for j in same_species}
    pair_available = {"forward", "reverse"} <= roles
    pair = pair_available and request.args.get("layout") != "single"
    if pair:
        order = {"forward": 0, "reverse": 1}
        section_jobs = sorted(same_species,
                              key=lambda j: (order.get(j.get("role"), 2), j["job_number"]))
    else:
        section_jobs = [job]
    sections = [{
        "job": j, "label": PRIMER_LABELS.get(j.get("role")) or j.get("seq_desc") or j["summary"],
        "downloads": _download_list(t, j),
        "queue_position": blastjobs.queue_position(j["id"]) if j["status"] == "queued" else None,
    } for j in section_jobs]
    species_tabs = []
    for j in t["jobs"]:
        if j["species"] not in {s["species"] for s in species_tabs}:
            m_tab = species.load(j["species"]) or {}
            species_tabs.append({"species": j["species"], "job_id": j["id"],
                                 "label": m_tab.get("display_name") or j["species"]})

    if pair:
        forward = next(j for j in section_jobs if j.get("role") == "forward")
        base = re.sub(r"([\s_\-]+(forward|fwd|left)([\s_\-]+primer)?|[\s_\-]+(f|fw))$", "",
                      forward.get("seq_desc") or "", flags=re.I).strip()
        page_title = t["description"] or (f"Primer pair: {base}" if base else "Primer pair")
    else:
        page_title = job["description"] or job["summary"]
    js_sections = [{"job": {k: s["job"][k] for k in ("id", "status", "n_hits", "started_at")},
                    "seq_len": len(s["job"]["sequence"]), "role": s["job"].get("role")}
                   for s in sections]

    m = species.load(job["species"]) or {"name": job["species"]}
    return render_template(
        "blast_results.html", env=_env(), ticket=t, job=job, sections=sections,
        js_sections=js_sections, page_title=page_title,
        pair=pair, pair_available=pair_available, species_tabs=species_tabs,
        species_info={**species.summary(m), "karyotype": species.karyotype(m),
                      "ensembl_url": f"{species.ensembl_site(m)}/{species.web_name(m['name'])}",
                      "has_annotation": bool((m.get("gtf") or {}).get("sqlite"))},
        config_groups=blastconf.describe_configs(t["search_type"], t["configs"]),
        genomic=t["source"] in blastconf.GENOMIC_SOURCES,
        is_blat=t["search_type"] == blastconf.BLAT_VALUE)


def _download_list(ticket: dict, job: dict) -> list[dict]:
    if job["status"] != "done":
        return []
    if job.get("engine") == "blat":
        items = [{"fmt": k, "label": v[1]} for k, v in engines.BLAT_FORMATS.items()]
    else:
        items = [{"fmt": k, "label": v[1]} for k, v in engines.NCBI_FORMATS.items()
                 if k != "sam" or ticket["search_type"] == "NCBIBLAST_BLASTN"]
    items += [{"fmt": "table_csv", "label": "Results table (CSV)"},
              {"fmt": "table_tsv", "label": "Results table (TSV)"},
              {"fmt": "query", "label": "Query sequence (FASTA)"}]
    return items


VIEW_NAMES = {"alignment": "Alignment", "query": "Query sequence", "genomic": "Genomic sequence"}


def _saved_view_options(view: str) -> dict:
    """Options a sequence view page last saved (JSON, URL-encoded by the page's script)."""
    import urllib.parse
    raw = request.cookies.get(f"pf_view_{view}", "")
    for candidate in (urllib.parse.unquote(raw), raw):
        try:
            value = json.loads(candidate)
            if isinstance(value, dict):
                return value
        except ValueError:
            continue
    return {}


@app.route("/blast/job/<jid>/hit/<int:idx>")
def page_blast_hit(jid: str, idx: int):
    job = blastjobs.get_job(jid)
    if not job:
        abort(404)
    ticket = _ticket_or_404(job["ticket_id"])
    hit = blastjobs.get_hit(jid, idx)
    if not hit:
        abort(404)
    view = request.args.get("view", "alignment")
    if view not in VIEW_NAMES:
        view = "alignment"
    opts = seqviews.options(view, {**_saved_view_options(view), **request.args.to_dict()})
    m = species.load(job["species"]) or {"name": job["species"]}
    error = None
    result = None
    try:
        if view == "alignment":
            result = seqviews.alignment_view(ticket, job, hit, m, opts)
        elif view == "query":
            result = seqviews.query_view(ticket, job, hit, blastjobs.get_hits(jid), opts)
        else:
            result = seqviews.genomic_view(ticket, job, hit, blastjobs.get_hits(jid), m, opts)
    except seqviews.ViewError as exc:
        error = str(exc)
    return render_template(
        "blast_hit.html", env=_env(), ticket=ticket, job=job, hit=hit, view=view,
        view_names=VIEW_NAMES, result=result, error=error,
        option_fields=seqviews.option_form(view, opts), opts=opts,
        summary=seqviews.hit_summary(ticket, hit), n_hits=job["n_hits"],
        ensembl_url=f"{species.ensembl_site(m)}/{species.web_name(m['name'])}",
        genomic=ticket["source"] in blastconf.GENOMIC_SOURCES)


# --------------------------------------------------------------------------
# BLAST/BLAT API
# --------------------------------------------------------------------------
@app.get("/api/blast/form")
def api_blast_form():
    return jsonify({**blastconf.form_options(), "species": _species_payload()})


@app.post("/api/blast/parse")
def api_blast_parse():
    return jsonify(blastconf.parse_sequences((request.json or {}).get("text", "")))


@app.post("/api/blast/fetch")
def api_blast_fetch():
    try:
        return jsonify({"sequences": seqfetch.fetch((request.json or {}).get("id", ""))})
    except seqfetch.FetchError as exc:
        return jsonify({"error": str(exc)}), 400


@app.post("/api/blast/upload")
def api_blast_upload():
    f = request.files.get("file")
    if not f or not f.filename:
        return jsonify({"error": "Choose a file to upload."}), 400
    if not f.filename.lower().endswith((".fa", ".fasta", ".fas", ".fna", ".faa", ".txt",
                                        ".seq")):
        return jsonify({"error": "Uploaded file should be of type plain text or FASTA."}), 400
    limit = int(blastconf.MAX_SEQUENCE_LENGTH * blastconf.MAX_NUM_SEQUENCES * 1.1)
    data = f.read(limit + 1)
    if len(data) > limit:
        return jsonify({"error": f"Uploaded file should not be more than "
                                 f"{limit // (1024 * 1024)} MB."}), 400
    return jsonify({"text": data.decode("utf-8", errors="replace")})


@app.post("/api/blast/submit")
def api_blast_submit():
    try:
        ticket_id = blastjobs.submit(request.json or {})
    except blastjobs.SubmitError as exc:
        return jsonify({"error": str(exc)}), 400
    blastjobs.set_ticket_owner(ticket_id, _username())
    return jsonify({"ticket_id": ticket_id})


@app.get("/api/blast/tickets")
def api_blast_tickets():
    tickets = blastjobs.list_tickets(limit=int(request.args.get("limit", 100)))
    positions = {}
    queued = [j for t in tickets for j in t["jobs"] if j["status"] == "queued"]
    for j in sorted(queued, key=lambda j: (j["created_at"], j["job_number"])):
        positions[j["id"]] = len(positions) + 1
    for t in tickets:
        for j in t["jobs"]:
            j["queue_position"] = positions.get(j["id"])
    return jsonify(tickets)


@app.get("/api/blast/ticket/<ticket_id>")
def api_blast_ticket(ticket_id: str):
    return jsonify(_ticket_or_404(ticket_id))


@app.post("/api/blast/ticket/<ticket_id>/delete")
def api_blast_ticket_delete(ticket_id: str):
    ticket = blastjobs.get_ticket(ticket_id, with_jobs=False)
    if ticket and not _may_change(ticket.get("owner")):
        return _not_yours("ticket")
    blastjobs.delete_ticket(ticket_id)
    return jsonify({"ok": True})


@app.post("/api/blast/job/<jid>/<action>")
def api_blast_job_action(jid: str, action: str):
    job = blastjobs.get_job(jid)
    if not job:
        abort(404)
    ticket = blastjobs.get_ticket(job["ticket_id"], with_jobs=False) or {}
    if not _may_change(ticket.get("owner")):
        return _not_yours("ticket")
    if action == "cancel":
        return jsonify({"ok": blastjobs.cancel(jid)})
    if action == "rerun":
        return jsonify({"ok": blastjobs.resubmit_job(jid)})
    if action == "delete":
        return jsonify({"ok": True, "ticket_left": blastjobs.delete_job(jid)})
    abort(404)


@app.get("/api/blast/job/<jid>")
def api_blast_job(jid: str):
    job = blastjobs.get_job(jid)
    if not job:
        abort(404)
    job.pop("sequence", None)
    if job["status"] == "queued":
        job["queue_position"] = blastjobs.queue_position(jid)
    return jsonify(job)


@app.get("/api/blast/job/<jid>/hits")
def api_blast_hits(jid: str):
    if not blastjobs.get_job(jid):
        abort(404)
    return jsonify(blastjobs.get_hits(jid))


@app.post("/api/blast/job/<jid>/genes")
def api_blast_genes(jid: str):
    """Overlapping genes for hits placed without a local annotation index (Ensembl REST)."""
    job = blastjobs.get_job(jid)
    if not job:
        abort(404)
    m = species.load(job["species"]) or {"name": job["species"]}
    wanted = [int(i) for i in (request.json or {}).get("idx", [])][:60]
    found = {}
    for idx in wanted:
        hit = blastjobs.get_hit(jid, idx)
        if not hit or not hit.get("gid") or hit.get("genes") is not None:
            continue
        genes = annotation.genes_in_region(m, hit["gid"], hit["gstart"], hit["gend"])
        found[idx] = [{"id": g["id"], "name": g["name"], "biotype": g["biotype"],
                       "strand": g["strand"]} for g in genes]
    blastjobs.update_hit_genes(jid, found)
    return jsonify({str(k): v for k, v in found.items()})


@app.post("/api/blast/annotation/<name>/fetch")
def api_blast_annotation_fetch(name: str):
    """Fetch (and cache) the REST annotation a sequence view rendered without."""
    m = species.load(name)
    if not m:
        abort(404)
    fetched, failed = 0, 0
    for item in ((request.json or {}).get("items") or [])[:6]:
        fn = seqviews.lookup_function(item.get("kind"))
        try:
            fn(m, str(item["chrom"]), int(item["start"]), int(item["end"]))
            fetched += 1
        except Exception:                                         # noqa: BLE001
            failed += 1
    # A failed REST call caches nothing, so report whether anything is now available.
    still = []
    for item in ((request.json or {}).get("items") or [])[:6]:
        fn = seqviews.lookup_function(item.get("kind"))
        try:
            fn(m, str(item["chrom"]), int(item["start"]), int(item["end"]), cached_only=True)
        except annotation.NotCached:
            still.append(item)
        except Exception:                                         # noqa: BLE001
            still.append(item)
    return jsonify({"fetched": fetched, "failed": failed, "missing": still})


def _table_rows(ticket: dict, hits: list[dict]) -> tuple[list[str], list[list]]:
    genomic = ticket["source"] in blastconf.GENOMIC_SOURCES
    ori = lambda v: "Forward" if (v or 1) > 0 else "Reverse"            # noqa: E731
    loc = lambda h: f"{h['gid']}:{h['gstart']}-{h['gend']}" if h.get("gid") else ""  # noqa
    genes = lambda h: ", ".join(g["name"] for g in (h.get("genes") or []))          # noqa
    if genomic:
        head = ["Genomic Location", "Overlapping Gene(s)", "Orientation"]
        lead = lambda h: [loc(h), genes(h), ori(h.get("gori"))]                     # noqa
    else:
        head = ["Subject name", "Gene hit", "Subject start", "Subject end", "Subject ori",
                "Genomic Location", "Orientation"]
        lead = lambda h: [h["tid"], (h.get("gene_hit") or {}).get("name", ""),     # noqa
                          h["tstart"], h["tend"], ori(h["tori"]), loc(h),
                          ori(h.get("gori")) if h.get("gori") else ""]
    head += ["Query name", "Query start", "Query end", "Query ori", "Length", "Score",
             "E-val", "%ID"]
    rows = [lead(h) + [h["qid"], h["qstart"], h["qend"], ori(h["qori"]), h["len"],
                       h["score"], f"{h['evalue']:.3g}", h["pident"]] for h in hits]
    return head, rows


@app.get("/blast/job/<jid>/download/<fmt>")
def blast_download(jid: str, fmt: str):
    job = blastjobs.get_job(jid)
    if not job:
        abort(404)
    ticket = _ticket_or_404(job["ticket_id"])
    stem = f"{job['ticket_id']}_{job['job_number']}"
    workdir = blastjobs.job_dir(job["ticket_id"], jid)
    if fmt == "query":
        return Response(blastconf.fasta(job.get("seq_desc") or f"query_{job['job_number']}",
                                        job["sequence"]), mimetype="text/plain",
                        headers={"Content-Disposition": f"attachment; filename={stem}_query.fa"})
    if job["status"] != "done":
        abort(409)
    if fmt in ("table_csv", "table_tsv"):
        head, rows = _table_rows(ticket, blastjobs.get_hits(jid))
        buf = io.StringIO()
        w = csv.writer(buf, delimiter="," if fmt == "table_csv" else "\t", lineterminator="\n")
        w.writerow(head)
        w.writerows(rows)
        ext = "csv" if fmt == "table_csv" else "tsv"
        return Response(buf.getvalue(), mimetype=f"text/{ext}",
                        headers={"Content-Disposition":
                                 f"attachment; filename={stem}_results.{ext}"})
    if job.get("engine") == "blat":
        spec = engines.BLAT_FORMATS.get(fmt)
        if not spec or not (workdir / spec[0]).exists():
            abort(404)
        path, mime, ext = workdir / spec[0], spec[2], spec[3]
    else:
        spec = engines.NCBI_FORMATS.get(fmt)
        if not spec:
            abort(404)
        try:
            path = engines.render_ncbi_format(workdir, fmt)
        except engines.EngineError as exc:
            return Response(str(exc), status=500, mimetype="text/plain")
        mime, ext = spec[2], spec[3]
    return send_file(path, mimetype=mime, as_attachment=True,
                     download_name=f"{stem}_{fmt}.{ext}")


@app.get("/blast/job/<jid>/hit/<int:idx>/genomic.fa")
def blast_genomic_fasta(jid: str, idx: int):
    job = blastjobs.get_job(jid)
    hit = blastjobs.get_hit(jid, idx) if job else None
    if not hit:
        abort(404)
    ticket = _ticket_or_404(job["ticket_id"])
    m = species.load(job["species"]) or {}
    opts = seqviews.options("genomic", request.args.to_dict())
    opts.update(hsp_display="off", exon_display="off", snp_display="off")
    try:
        view = seqviews.genomic_view(ticket, job, hit, [], m, opts)
    except seqviews.ViewError as exc:
        return Response(str(exc), status=404, mimetype="text/plain")
    r = view["region"]
    return Response(view["fasta"], mimetype="text/plain", headers={
        "Content-Disposition": f"attachment; filename={r['chrom']}_{r['start']}-{r['end']}.fa"})


# --------------------------------------------------------------------------
# Transcripts: MANE Select, RefSeq accessions and the aligned cDNA view
# --------------------------------------------------------------------------
def _transcript_species(name: str | None) -> dict:
    installed = [s for s in species.all_species() if s["components"]["gtf"]]
    if name:
        m = species.load(name)
        if m:
            return m
    for s in installed:
        return species.load(s["name"])
    abort(404)


@app.route("/transcripts")
def page_transcripts():
    query = request.args.get("q", "").strip()
    name = request.args.get("species") or ""
    options = [s for s in species.all_species() if s["components"]["gtf"]]
    m = _transcript_species(name) if options else None
    results, error = [], None
    if m and query:
        try:
            results = transcripts.search(m, query)
        except Exception as exc:                                  # noqa: BLE001
            error = str(exc)
        if not results:
            error = (f"Nothing matched '{query}'. Search by gene symbol (GALE), "
                     "Ensembl ID (ENSG/ENST), RefSeq accession (NM_001008216) or HGNC ID.")
    return render_template("transcripts.html", env=_env(), query=query, results=results,
                           error=error, species_options=options,
                           species_name=m["name"] if m else "",
                           has_ids=bool(m and transcripts.available(m)),
                           ensembl_url=f"{species.ensembl_site(m)}/{species.web_name(m['name'])}"
                           if m else "")


@app.route("/transcript/<name>/<tid>")
def page_transcript(name: str, tid: str):
    m = species.load(name)
    if not m:
        abort(404)
    # RefSeq accessions are shown through the Ensembl transcript they are matched to.
    requested = tid
    resolved = tid if tid.upper().startswith("ENS") else (transcripts.ensembl_for(m, tid) or tid)
    # The last tab and options chosen on this page, overridden by the URL.
    opts = transcriptview.options({**_saved_view_options("transcript"),
                                   **request.args.to_dict()})
    view = error = None
    try:
        # Opened by RefSeq accession: show that accession's own sequence in the cDNA row.
        view = transcriptview.build(m, resolved, opts,
                                    compare_with=None if tid.upper().startswith("ENS") else tid)
    except transcriptview.ViewError as exc:
        error = str(exc)
    tx = (view or {}).get("transcript") or annotation.transcript(m, resolved)
    gene = transcripts.gene_detail(m, tx["gene_id"]) if tx else None
    return render_template(
        "transcript.html", env=_env(), species_info=species.summary(m), species_name=name,
        requested=requested, transcript=tx, gene=gene, view=view, error=error, opts=opts,
        option_fields=transcriptview.option_form(opts), tabs=transcriptview.TABS,
        csq_css=transcriptview.consequence_css(),
        ids=transcripts.annotate(m, resolved) if tx else {},
        summary=transcriptview.summary(m, tx, view) if (tx and view) else [],
        ensembl_url=f"{species.ensembl_site(m)}/{species.web_name(name)}")


@app.get("/api/transcripts/search")
def api_transcripts_search():
    m = _transcript_species(request.args.get("species"))
    return jsonify(transcripts.search(m, request.args.get("q", ""),
                                      limit=int(request.args.get("limit", 15))))


@app.get("/transcript/<name>/<tid>/<kind>.fa")
def transcript_fasta(name: str, tid: str, kind: str):
    m = species.load(name)
    if not m or kind not in ("cdna", "protein", "gene"):
        abort(404)
    resolved = tid if tid.upper().startswith("ENS") else (transcripts.ensembl_for(m, tid) or tid)
    try:
        if kind == "gene":
            # The entire genomic sequence with the requested flanks, never shortened.
            flank5 = int(transcriptview.options(request.args.to_dict())["flank5"])
            flank3 = int(transcriptview.options(request.args.to_dict())["flank3"])
            region = transcriptview.genomic_sequence(m, resolved, flank5, flank3)
            gene_name = region["transcript"].get("gene_name") or ""
            label = (f"{tid} {gene_name} {region['chrom']}:{region['start']}-{region['end']}:"
                     f"{region['strand']} genomic, 5' flank {flank5} bp, 3' flank {flank3} bp")
            seq = region["sequence"]
        else:
            view = transcriptview.build(m, resolved, transcriptview.options({"layout": "spliced"}),
                                        compare_with=None if tid.upper().startswith("ENS") else tid)
            gene_name = view["transcript"].get("gene_name") or ""
            if kind == "protein":
                seq, label = view["protein"], f"{tid} {gene_name} protein"
            else:
                seq, label = view["cdna"], f"{view.get('reference') or tid} {gene_name} cDNA"
    except transcriptview.ViewError as exc:
        return Response(str(exc), status=404, mimetype="text/plain")
    return Response(blastconf.fasta(label.strip(), seq), mimetype="text/plain", headers={
        "Content-Disposition": f"attachment; filename={tid}_{kind}.fa"})


# --------------------------------------------------------------------------
# Species databases and the BLAT server
# --------------------------------------------------------------------------
_installs: dict[str, subprocess.Popen] = {}


@app.get("/api/species")
def api_species():
    sys.path.insert(0, str(ROOT))
    from tools.setup_species import read_status
    out = []
    for s in _species_payload():
        proc = _installs.get(s["name"])
        status = read_status(s["name"])
        out.append({**s, "install": status,
                    "installing": (proc is not None and proc.poll() is None)
                    or _pid_running(status)})
    return jsonify(out)


@app.get("/api/species/ensembl")
def api_species_ensembl():
    cached = store.cache_get("ensembl:species", max_age=7 * 86400)
    if cached is None:
        try:
            from primerforge import net
            data = net.fetch_json(f"{config.ENSEMBL_REST}/info/species?"
                                  "content-type=application/json", timeout=90)
        except Exception as exc:                                  # noqa: BLE001
            return jsonify({"error": f"Could not reach Ensembl: {exc}"}), 502
        cached = sorted(({"name": s["name"], "display_name": s.get("display_name"),
                          "assembly": s.get("assembly"), "common_name": s.get("common_name"),
                          "taxon_id": s.get("taxon_id")}
                         for s in data.get("species", [])),
                        key=lambda s: (s["display_name"] or s["name"]).lower())
        store.cache_put("ensembl:species", cached)
    return jsonify(cached)


@app.post("/api/species/install")
def api_species_install():
    data = request.json or {}
    name = str(data.get("species", "")).strip().lower()
    comps = [c for c in data.get("components") or species.INSTALL_ORDER
             if c in species.COMPONENTS]
    if not name or not name.replace("_", "").isalnum():
        return jsonify({"error": "Choose a species."}), 400
    proc = _installs.get(name)
    sys.path.insert(0, str(ROOT))
    from tools.setup_species import read_status
    if (proc is not None and proc.poll() is None) or _pid_running(read_status(name)):
        return jsonify({"ok": True, "already_running": True})
    config.ensure_dirs()
    log_dir = species.species_dir(name)
    log_dir.mkdir(parents=True, exist_ok=True)
    log = open(log_dir / "install.log", "a", encoding="utf-8")
    _installs[name] = subprocess.Popen(
        [sys.executable, str(ROOT / "tools" / "setup_species.py"), name,
         "--components", ",".join(comps)],
        stdout=log, stderr=subprocess.STDOUT, cwd=str(ROOT),
        creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
    return jsonify({"ok": True})


@app.get("/api/blat/status")
def api_blat_status():
    distros = wsl.distros()
    settings = config.settings()
    try:
        blat_dir = wsl.linux_path(settings["blat_dir"]) if wsl.distro() else settings["blat_dir"]
    except (wsl.WslError, OSError, subprocess.SubprocessError):
        blat_dir = settings["blat_dir"]
    out = {"wsl": {"available": wsl.NATIVE or bool(wsl.exe()), "native": wsl.NATIVE,
                   "distros": distros, "distro": wsl.distro(), "settings": settings,
                   "blat_dir": blat_dir},
           "servers": []}
    for s in species.all_species():
        if s["components"]["blat"]:
            out["servers"].append({**blatserver.status(s["name"]),
                                   "display_name": s["display_name"]})
    return jsonify(out)


@app.post("/api/blat/<name>/<action>")
def api_blat_action(name: str, action: str):
    if action == "start":
        blatserver.start(name)
    elif action == "stop":
        threading.Thread(target=blatserver.stop, args=(name,), daemon=True).start()
    else:
        abort(404)
    return jsonify({"ok": True})


@app.post("/api/settings")
def api_settings():
    data = request.json or {}
    clean = {}
    if "wsl_distro" in data:
        clean["wsl_distro"] = str(data["wsl_distro"]).strip()
    if "blast_workers" in data:
        clean["blast_workers"] = max(1, min(8, int(data["blast_workers"])))
    if "blat_autostart" in data:
        clean["blat_autostart"] = bool(data["blat_autostart"])
    saved = config.save_settings(clean)
    wsl.reset_cache()
    return jsonify(saved)


def _autostart_blat() -> None:
    if not config.settings().get("blat_autostart"):
        return
    for s in species.all_species():
        if s["blat"]:
            blatserver.start(s["name"])


_started = False
_start_lock = threading.Lock()


def _startup() -> None:
    """Database, species registry, BLAST job workers and BLAT servers, once per process."""
    global _started
    with _start_lock:
        if _started:
            return
        _started = True
    store.init()
    auth.init()
    config.ensure_dirs()
    species.ensure_human()
    blastjobs.init()
    threading.Thread(target=_autostart_blat, daemon=True).start()


@app.before_request
def _ensure_started() -> None:
    # `flask run` imports the app without calling main(); start services on first use.
    if not _started:
        _startup()


# --------------------------------------------------------------------------
# Accounts (server installs only; the desktop app has no logins)
# --------------------------------------------------------------------------
PUBLIC_ENDPOINTS = {"static", "page_login", "api_health"}
ADMIN_ENDPOINTS = {"page_settings", "api_setup_start", "api_setup_reload", "api_species_install",
                   "api_blat_action", "api_settings", "page_users", "api_users_create",
                   "api_users_action"}


def _wants_json() -> bool:
    return request.path.startswith("/api/") or request.method != "GET"


@app.before_request
def _require_login():
    g.user = None
    if not config.AUTH_ENABLED:
        return None
    uid = session.get("uid")
    if uid:
        user = auth.get_user(uid)
        if user and user["active"] and session.get("tok") == auth.session_token(user):
            g.user = user
        else:
            session.clear()
    endpoint = request.endpoint or ""
    if endpoint in PUBLIC_ENDPOINTS:
        return None
    if g.user is None:
        if _wants_json():
            return jsonify({"error": "Your session has ended. Please sign in again.",
                            "login": url_for("page_login")}), 401
        target = request.full_path if request.query_string else request.path
        return redirect(url_for("page_login", next=target))
    if g.user["must_change"] and endpoint not in ("page_account", "page_logout"):
        if _wants_json():
            return jsonify({"error": "Choose a new password first."}), 403
        return redirect(url_for("page_account"))
    if endpoint in ADMIN_ENDPOINTS and g.user["role"] != "admin":
        if _wants_json():
            return jsonify({"error": "Only an administrator can do that."}), 403
        abort(403)
    return None


@app.context_processor
def _account_context() -> dict:
    user = g.get("user")
    return {"auth_enabled": config.AUTH_ENABLED, "current_user": user,
            "is_admin": not config.AUTH_ENABLED or bool(user and user["role"] == "admin")}


def _username() -> str | None:
    user = g.get("user")
    return user["username"] if user else None


def _may_change(owner: str | None) -> bool:
    """Admins may change anything; users only what they created."""
    user = g.get("user")
    if not config.AUTH_ENABLED:
        return True
    if not user:
        return False
    return user["role"] == "admin" or bool(owner and owner.lower() == user["username"].lower())


def _not_yours(what: str):
    return jsonify({"error": f"Only the person who created this {what} or an administrator "
                             "can change it."}), 403


def _safe_next(target: str | None) -> str:
    target = (target or "").strip()
    if not target.startswith("/") or target.startswith("//") or "\\" in target:
        return url_for("page_design")
    return target


@app.get("/api/health")
def api_health():
    return jsonify({"ok": True, "auth": config.AUTH_ENABLED})


@app.route("/login", methods=["GET", "POST"])
def page_login():
    if not config.AUTH_ENABLED:
        return redirect(url_for("page_design"))
    error = None
    next_url = _safe_next(request.values.get("next"))
    if request.method == "POST":
        try:
            user = auth.authenticate(request.form.get("username", ""),
                                     request.form.get("password", ""),
                                     request.remote_addr or "unknown")
        except auth.AuthError as exc:
            error = str(exc)
        else:
            session.clear()
            session.permanent = True
            session["uid"], session["tok"] = user["id"], auth.session_token(user)
            return redirect(next_url)
    elif g.get("user"):
        return redirect(next_url)
    return render_template("login.html", error=error, next_url=next_url,
                           username=request.form.get("username", "")), (401 if error else 200)


@app.post("/logout")
def page_logout():
    session.clear()
    return redirect(url_for("page_login") if config.AUTH_ENABLED else url_for("page_design"))


@app.route("/account", methods=["GET", "POST"])
def page_account():
    if not config.AUTH_ENABLED:
        return redirect(url_for("page_design"))
    message = error = None
    if request.method == "POST":
        user = auth.get_user(g.user["id"])
        from werkzeug.security import check_password_hash
        new = request.form.get("new_password", "")
        if not check_password_hash(user["password_hash"], request.form.get("password", "")):
            error = "Your current password is not correct."
        elif new != request.form.get("confirm_password", ""):
            error = "The new passwords do not match."
        else:
            try:
                auth.set_password(user["id"], new)
            except auth.AuthError as exc:
                error = str(exc)
            else:
                g.user = auth.get_user(user["id"])
                session["tok"] = auth.session_token(g.user)
                if user["must_change"]:
                    return redirect(url_for("page_design"))
                message = "Password changed."
    return render_template("account.html", env=_env(), message=message, error=error,
                           min_password=auth.MIN_PASSWORD)


@app.get("/admin/users")
def page_users():
    if not config.AUTH_ENABLED:
        return redirect(url_for("page_design"))
    return render_template("users.html", env=_env(), users=auth.list_users(),
                           roles=auth.ROLES, auth_enabled=config.AUTH_ENABLED)


@app.post("/api/users")
def api_users_create():
    data = request.json or {}
    password = auth.generate_password()
    try:
        user = auth.create_user(str(data.get("username", "")), password,
                                role=str(data.get("role", "user")),
                                full_name=str(data.get("full_name", "")), must_change=True)
    except auth.AuthError as exc:
        return jsonify({"error": str(exc)}), 400
    return jsonify({"user": user, "password": password})


@app.post("/api/users/<int:user_id>/<action>")
def api_users_action(user_id: int, action: str):
    data = request.json or {}
    if user_id == g.user["id"] and action in ("disable", "delete", "role", "reset-password"):
        return jsonify({"error": "Change your own password on the Account page. Your role "
                                 "and account can only be changed by another administrator."}), 400
    if not auth.get_user(user_id):
        abort(404)
    try:
        if action == "reset-password":
            password = auth.generate_password()
            auth.set_password(user_id, password, must_change=True)
            return jsonify({"ok": True, "password": password})
        if action == "role":
            auth.set_role(user_id, str(data.get("role", "")))
        elif action in ("enable", "disable"):
            auth.set_active(user_id, action == "enable")
        elif action == "delete":
            auth.delete_user(user_id)
        else:
            abort(404)
    except auth.AuthError as exc:
        return jsonify({"error": str(exc)}), 400
    return jsonify({"ok": True})


def main() -> None:
    _startup()
    port = config.PORT
    url = f"http://127.0.0.1:{port}"
    if "--no-browser" not in sys.argv:
        threading.Timer(1.2, lambda: webbrowser.open(url)).start()
    print(f"\n  PrimerForge running at {url}\n  Press Ctrl+C to stop.\n")
    app.run(host=config.HOST, port=port, debug="--debug" in sys.argv,
            use_reloader="--debug" in sys.argv, threaded=True)


if __name__ == "__main__":
    main()
