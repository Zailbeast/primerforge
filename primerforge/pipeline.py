"""Run orchestration: input text in, a stored and fully annotated run out."""
from __future__ import annotations

import threading
import time
from dataclasses import dataclass, field

from . import design, presets, refseq, specificity, store, variants
from .config import blast_ready
from .variants import VariantError


@dataclass
class Progress:
    run_id: str
    total: int
    done: int = 0
    stage: str = "starting"
    current: str = ""
    finished: bool = False
    error: str | None = None
    started_at: float = field(default_factory=time.time)

    def as_dict(self) -> dict:
        pct = 100.0 * self.done / self.total if self.total else 0.0
        return {"run_id": self.run_id, "total": self.total, "done": self.done,
                "pct": round(pct, 1), "stage": self.stage, "current": self.current,
                "finished": self.finished, "error": self.error,
                "elapsed": round(time.time() - self.started_at, 1)}


_progress: dict[str, Progress] = {}
_lock = threading.Lock()


def get_progress(run_id: str) -> dict | None:
    with _lock:
        p = _progress.get(run_id)
    return p.as_dict() if p else None


def _set(run_id: str, **kw) -> None:
    with _lock:
        p = _progress.get(run_id)
        if p:
            for k, v in kw.items():
                setattr(p, k, v)


def start_run(input_text: str, assay: str, overrides: dict, label: str = "",
              assembly: str = "GRCh38", check_specificity: bool = True) -> str:
    """Validate, create the run row, and design in a background thread."""
    items = variants.parse_input_block(input_text)
    if not items:
        raise VariantError("Enter at least one variant.")
    if len(items) > 200:
        raise VariantError(f"{len(items)} variants supplied; the limit per run is 200.")

    params = presets.build(assay, overrides)          # raises ValueError when invalid
    params["assembly"] = assembly
    params["check_specificity"] = bool(check_specificity) and blast_ready()

    run_id = store.create_run(label=label.strip() or _auto_label(items, assay),
                              assay=assay, input_raw=input_text, params=params,
                              n_input=len(items), specificity=params["check_specificity"])
    with _lock:
        _progress[run_id] = Progress(run_id=run_id, total=len(items))

    thread = threading.Thread(target=_execute, args=(run_id, items, params, assembly),
                              daemon=True)
    thread.start()
    return run_id


def _auto_label(items: list[str], assay: str) -> str:
    head = items[0][:44]
    extra = f" +{len(items) - 1} more" if len(items) > 1 else ""
    return f"{head}{extra}"


def _execute(run_id: str, items: list[str], params: dict, assembly: str) -> None:
    n_ok = n_failed = 0
    try:
        for idx, raw in enumerate(items):
            _set(run_id, done=idx, current=raw, stage="resolving")
            try:
                locus = variants.resolve(raw, assembly)
            except VariantError as exc:
                store.add_variant(run_id, idx, raw, "failed", error=str(exc))
                n_failed += 1
                continue
            except Exception as exc:                      # noqa: BLE001
                store.add_variant(run_id, idx, raw, "failed",
                                  error=f"Unexpected error while resolving: {exc}")
                n_failed += 1
                continue

            try:
                _set(run_id, stage="designing", current=locus.label())
                template, pairs, warns = design.design(locus, params)

                _set(run_id, stage="checking common variants")
                try:
                    design.annotate_primer_variants(pairs, locus, params["snp_maf"])
                except Exception:                          # noqa: BLE001
                    pass                                   # annotation is advisory

                if params["check_specificity"]:
                    _set(run_id, stage="in-silico PCR")
                    try:
                        specificity.check_pairs(
                            pairs, locus.chrom, params["product_min"],
                            params["product_max"])
                    except Exception as exc:               # noqa: BLE001
                        for p in pairs:
                            p["specificity"] = {"status": "error", "message": str(exc)}

                # Warnings accrue after the first sort, so rank once more at the end.
                pairs.sort(key=lambda p: (
                    0 if (p.get("specificity") or {}).get("status") == "unique" else 1,
                    len(p["warnings"]), p["penalty"]))
                for rank, p in enumerate(pairs):
                    p["rank"] = rank

                vid = store.add_variant(run_id, idx, raw, "ok", locus=locus.to_dict(),
                                        template=_thin(template), warnings=warns)
                store.add_pairs(vid, run_id, pairs)
                n_ok += 1
            except (design.DesignError, refseq.SequenceError) as exc:
                store.add_variant(run_id, idx, raw, "failed", error=str(exc),
                                  locus=locus.to_dict())
                n_failed += 1
            except Exception as exc:                      # noqa: BLE001
                store.add_variant(run_id, idx, raw, "failed",
                                  error=f"Unexpected error: {exc}",
                                  locus=locus.to_dict())
                n_failed += 1

        _set(run_id, done=len(items), stage="done", finished=True)
        status = "ok" if n_failed == 0 else ("partial" if n_ok else "failed")
        store.finish_run(run_id, status, n_ok, n_failed)
    except Exception as exc:                              # noqa: BLE001
        _set(run_id, finished=True, error=str(exc), stage="error")
        store.finish_run(run_id, "failed", n_ok, n_failed, error=str(exc))


def _thin(template: dict) -> dict:
    """Store template metadata without the multi-kilobase sequence."""
    return {k: v for k, v in template.items() if k != "seq"} | {
        "length": len(template["seq"])}
