"""Primer design around a resolved variant.

Pulls reference sequence around the locus, hands Primer3 a target region that
forces the variant to sit comfortably inside the product, then maps every
returned oligo back to genomic coordinates and annotates it.
"""
from __future__ import annotations

import primer3

from . import net, refseq, store
from .variants import Locus, revcomp

# Bases at the 3' end where a mismatch or SNP genuinely kills extension.
CRITICAL_3P = 5


class DesignError(RuntimeError):
    pass


# --------------------------------------------------------------------------
# Template
# --------------------------------------------------------------------------
def build_template(locus: Locus, params: dict) -> dict:
    """Reference window around the variant, plus a reference-allele sanity check."""
    flank = params["flank"]
    t_start = max(1, locus.start - flank)
    t_end = locus.end + flank
    seq = refseq.fetch(locus.chrom, t_start, t_end, locus.assembly)
    if not seq:
        raise DesignError(f"No reference sequence around {locus.chrom}:{locus.start}.")

    offset = locus.start - t_start          # 0-based index of the variant
    var_len = max(1, locus.end - locus.start + 1)

    warnings: list[str] = []
    observed = seq[offset:offset + var_len].upper()
    expected = (locus.ref or "").upper()
    if expected and expected != "-" and observed != expected:
        warnings.append(
            f"Reference mismatch: {locus.assembly} has '{observed}' at "
            f"{locus.chrom}:{locus.start}-{locus.end} but the variant states "
            f"'{expected}'. Check the transcript version and genome build.")

    return {
        "chrom": locus.chrom, "start": t_start, "end": t_start + len(seq) - 1,
        "seq": seq, "variant_offset": offset, "variant_len": var_len,
        "source": refseq.source(locus.assembly),
        "repeat_frac": round(refseq.repeat_fraction(seq), 4),
        "gc": round(refseq.gc_fraction(seq), 4),
        "warnings": warnings,
    }


# --------------------------------------------------------------------------
# Design
# --------------------------------------------------------------------------
def design(locus: Locus, params: dict) -> tuple[dict, list[dict], list[str]]:
    template = build_template(locus, params)
    seq = template["seq"]
    offset = template["variant_offset"]
    var_len = template["variant_len"]
    margin = params["margin"]
    warnings = list(template["warnings"])

    target_start = max(0, offset - margin)
    target_len = min(len(seq) - target_start, var_len + 2 * margin)
    if target_len <= 0:
        raise DesignError("The variant falls outside the retrieved template.")

    seq_args = {
        "SEQUENCE_ID": (locus.gene or locus.chrom) + f"_{locus.start}",
        "SEQUENCE_TEMPLATE": seq,
        "SEQUENCE_TARGET": [target_start, target_len],
    }

    p3 = dict(params["p3"])
    try:
        res = primer3.bindings.design_primers(seq_args=seq_args, global_args=p3)
    except (OSError, ValueError) as exc:
        raise DesignError(f"Primer3 failed: {exc}") from exc

    n = int(res.get("PRIMER_PAIR_NUM_RETURNED", 0) or 0)
    if n == 0:
        raise DesignError(_explain_failure(res, params))

    pairs = [_build_pair(i, res, template, locus, params) for i in range(n)]

    # Clean primers first, then Primer3's own penalty.
    pairs.sort(key=lambda p: (len(p["warnings"]), p["penalty"]))
    for rank, p in enumerate(pairs):
        p["rank"] = rank
    return template, pairs, warnings


def _explain_failure(res: dict, params: dict) -> str:
    """Turn Primer3's counters into something a user can act on."""
    explain = " ".join(str(res.get(k, "")) for k in
                       ("PRIMER_LEFT_EXPLAIN", "PRIMER_RIGHT_EXPLAIN",
                        "PRIMER_PAIR_EXPLAIN")).strip()
    hints = []
    low = explain.lower()
    if "high tm" in low or "low tm" in low:
        hints.append("widen the Tm window")
    if "high gc" in low or "low gc" in low or "gc clamp" in low:
        hints.append("relax the GC limits")
    if "high any compl" in low or "high end compl" in low or "hairpin" in low:
        hints.append("raise the dimer/hairpin thresholds")
    if "product size" in low or "unacceptable product size" in low:
        hints.append(f"widen the product range (currently "
                     f"{params['product_min']}-{params['product_max']} bp)")
    if "long poly" in low:
        hints.append("allow longer homopolymers")
    if not hints:
        hints.append(f"widen the product range or reduce the "
                     f"{params['margin']} bp variant-to-primer margin")
    msg = "No primer pair satisfied the constraints. Try to " + ", or ".join(hints) + "."
    return msg + (f"  [Primer3: {explain}]" if explain else "")


def _oligo(res: dict, side: str, i: int, template: dict, params: dict) -> dict:
    seq_t = template["seq"]
    pos, length = res[f"PRIMER_{side}_{i}"]
    if side == "LEFT":
        idx0 = int(pos)
    else:
        # Primer3 anchors right primers at their 5' (highest) template index.
        idx0 = int(pos) - int(length) + 1
    footprint = seq_t[idx0:idx0 + int(length)]
    g_start = template["start"] + idx0
    g_end = g_start + int(length) - 1
    seq_p = res[f"PRIMER_{side}_{i}_SEQUENCE"]
    return {
        "side": side.lower(), "seq": seq_p, "len": int(length),
        "tm": round(float(res[f"PRIMER_{side}_{i}_TM"]), 2),
        "gc": round(float(res[f"PRIMER_{side}_{i}_GC_PERCENT"]), 1),
        "self_any": round(float(res.get(f"PRIMER_{side}_{i}_SELF_ANY_TH", 0)), 1),
        "self_end": round(float(res.get(f"PRIMER_{side}_{i}_SELF_END_TH", 0)), 1),
        "hairpin": round(float(res.get(f"PRIMER_{side}_{i}_HAIRPIN_TH", 0)), 1),
        "end_stability": round(float(res.get(f"PRIMER_{side}_{i}_END_STABILITY", 0)), 2),
        "penalty": round(float(res.get(f"PRIMER_{side}_{i}_PENALTY", 0)), 3),
        "start": g_start, "end": g_end, "idx0": idx0,
        "template_footprint": footprint,
        "repeat_frac": round(refseq.repeat_fraction(footprint), 3),
        "gc_clamp_3p": footprint[-1].upper() in "GC" if side == "LEFT"
                       else footprint[0].upper() in "GC",
    }


def _build_pair(i: int, res: dict, template: dict, locus: Locus, params: dict) -> dict:
    left = _oligo(res, "LEFT", i, template, params)
    right = _oligo(res, "RIGHT", i, template, params)
    amp_start, amp_end = left["start"], right["end"]
    product = int(res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"])

    dist_left = locus.start - left["end"]        # variant 5' distance from left primer
    dist_right = right["start"] - locus.end

    warnings: list[str] = []
    for o, name in ((left, "Forward"), (right, "Reverse")):
        if o["repeat_frac"] > 0.5:
            warnings.append(f"{name} primer sits mostly in a repeat "
                            f"({o['repeat_frac'] * 100:.0f}% masked).")
        if o["hairpin"] and o["hairpin"] > params["p3"].get("PRIMER_MAX_HAIRPIN_TH", 24):
            warnings.append(f"{name} primer has a strong hairpin "
                            f"({o['hairpin']:.0f} C).")
        if not o["gc_clamp_3p"]:
            warnings.append(f"{name} primer has no G/C at its 3' end.")
    if abs(left["tm"] - right["tm"]) > 3:
        warnings.append(f"Tm mismatch of {abs(left['tm'] - right['tm']):.1f} C "
                        "between the two primers.")
    if min(dist_left, dist_right) < params["margin"]:
        warnings.append(f"Variant is only {min(dist_left, dist_right)} bp from a "
                        "primer; a Sanger read may not resolve it cleanly.")

    return {
        "rank": i, "left": left, "right": right,
        "product_size": product,
        "penalty": round(float(res[f"PRIMER_PAIR_{i}_PENALTY"]), 3),
        "compl_any": round(float(res.get(f"PRIMER_PAIR_{i}_COMPL_ANY_TH", 0)), 1),
        "compl_end": round(float(res.get(f"PRIMER_PAIR_{i}_COMPL_END_TH", 0)), 1),
        "tm_diff": round(abs(left["tm"] - right["tm"]), 2),
        "amplicon": {"chrom": locus.chrom, "start": amp_start, "end": amp_end},
        "amplicon_seq": template["seq"][left["idx0"]:right["idx0"] + right["len"]],
        "dist_left": dist_left, "dist_right": dist_right,
        "variant_pos_in_amplicon": locus.start - amp_start + 1,
        "warnings": warnings,
        "specificity": None,
        "snps": None,
    }


# --------------------------------------------------------------------------
# Common variants under the primers
# --------------------------------------------------------------------------
GNOMAD_API = "https://gnomad.broadinstitute.org/api"
GNOMAD_QUERY = """query($chrom:String!,$start:Int!,$stop:Int!){
  region(chrom:$chrom,start:$start,stop:$stop,reference_genome:GRCh38){
    variants(dataset:gnomad_r4){ variant_id pos rsids genome{af} exome{af} }
  }
}"""


def _gnomad_region(chrom: str, start: int, stop: int) -> list[dict] | None:
    """Population frequencies for a region, or None if gnomAD is unavailable.

    dbSNP alone is useless here -- it lists a variant at nearly every position.
    Only allele frequency separates a primer-killing common SNP from noise.
    """
    key = f"gnomad:{chrom}:{start}-{stop}"
    hit = store.cache_get(key, max_age=90 * 86400)
    if hit is not None:
        return hit
    try:
        payload = net.post_json(GNOMAD_API, {
            "query": GNOMAD_QUERY,
            "variables": {"chrom": chrom.replace("chr", ""),
                          "start": int(start), "stop": int(stop)}})
    except net.HttpError:
        return None                  # frequency annotation is advisory, never fatal
    region = (payload.get("data") or {}).get("region") or {}
    out = []
    for v in region.get("variants") or []:
        af = max((v.get("genome") or {}).get("af") or 0.0,
                 (v.get("exome") or {}).get("af") or 0.0)
        if af:
            out.append({"id": (v.get("rsids") or [v.get("variant_id")])[0],
                        "variant_id": v.get("variant_id"), "pos": v.get("pos"),
                        "af": round(af, 5)})
    store.cache_put(key, out)
    return out


def annotate_primer_variants(pairs: list[dict], locus: Locus, min_af: float = 0.01) -> None:
    """Flag common variants under each primer, weighted by 3'-end proximity."""
    if not pairs or locus.assembly != "GRCh38":
        return
    lo = min(p["left"]["start"] for p in pairs)
    hi = max(p["right"]["end"] for p in pairs)
    known = _gnomad_region(locus.chrom, lo, hi)
    if known is None:
        for pair in pairs:
            pair["snps"] = {"status": "unavailable", "left": [], "right": []}
        return

    common = [v for v in known if v["af"] >= min_af]
    for pair in pairs:
        hits: dict = {"status": "checked", "min_af": min_af, "left": [], "right": []}
        for side in ("left", "right"):
            o = pair[side]
            crit = set(range(o["end"] - CRITICAL_3P + 1, o["end"] + 1) if side == "left"
                       else range(o["start"], o["start"] + CRITICAL_3P))
            for v in common:
                if o["start"] <= v["pos"] <= o["end"]:
                    hits[side].append({**v, "critical": v["pos"] in crit})
        pair["snps"] = hits
        crit_hits = [h for s in ("left", "right") for h in hits[s] if h["critical"]]
        other = [h for s in ("left", "right") for h in hits[s] if not h["critical"]]
        if crit_hits:
            worst = max(crit_hits, key=lambda h: h["af"])
            pair["warnings"].append(
                f"{worst['id']} (AF {worst['af']:.1%}) sits within {CRITICAL_3P} bp of a "
                "primer 3' end - risk of allele dropout.")
        elif other:
            worst = max(other, key=lambda h: h["af"])
            pair["warnings"].append(
                f"Common variant {worst['id']} (AF {worst['af']:.1%}) lies under a primer.")
