"""Assay presets.

Each preset bundles the Primer3 thermodynamic settings with the tool-level
geometry that matters for the assay: how much reference sequence to pull, and
how far the variant must sit from either primer.
"""
from __future__ import annotations

BASE_P3 = {
    "PRIMER_TASK": "generic",
    "PRIMER_PICK_LEFT_PRIMER": 1,
    "PRIMER_PICK_RIGHT_PRIMER": 1,
    "PRIMER_PICK_INTERNAL_OLIGO": 0,
    "PRIMER_MAX_POLY_X": 4,
    "PRIMER_MAX_NS_ACCEPTED": 0,
    "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_SALT_DIVALENT": 1.5,
    "PRIMER_DNTP_CONC": 0.6,
    "PRIMER_DNA_CONC": 50.0,
    "PRIMER_TM_FORMULA": 1,          # SantaLucia 1998
    "PRIMER_SALT_CORRECTIONS": 1,
    "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1,
    "PRIMER_LOWERCASE_MASKING": 0,
}

PRESETS: dict[str, dict] = {
    "sanger": {
        "name": "Sanger confirmation",
        "blurb": "Flanking primers for dideoxy sequencing. Keeps the variant well "
                 "clear of both primers and of the noisy start of the read.",
        "flank": 1200,
        "margin": 120,               # min bases between variant and either primer
        "product": [[400, 900]],
        "opt_product": 600,
        "p3": {
            "PRIMER_OPT_SIZE": 20, "PRIMER_MIN_SIZE": 18, "PRIMER_MAX_SIZE": 27,
            "PRIMER_OPT_TM": 60.0, "PRIMER_MIN_TM": 57.0, "PRIMER_MAX_TM": 63.0,
            "PRIMER_PAIR_MAX_DIFF_TM": 3.0,
            "PRIMER_MIN_GC": 40.0, "PRIMER_OPT_GC_PERCENT": 50.0, "PRIMER_MAX_GC": 60.0,
            "PRIMER_MAX_SELF_ANY_TH": 45.0, "PRIMER_MAX_SELF_END_TH": 35.0,
            "PRIMER_PAIR_MAX_COMPL_ANY_TH": 45.0, "PRIMER_PAIR_MAX_COMPL_END_TH": 35.0,
            "PRIMER_MAX_HAIRPIN_TH": 24.0,
            "PRIMER_MAX_END_STABILITY": 9.0,
            "PRIMER_GC_CLAMP": 1,
        },
    },
    "longrange": {
        "name": "Long-range / gap-filling",
        "blurb": "Kilobase amplicons for long-range polymerases: longer, hotter "
                 "primers that survive a 68 C two-step protocol.",
        "flank": 6000,
        "margin": 300,
        "product": [[1000, 5000]],
        "opt_product": 3000,
        "p3": {
            "PRIMER_OPT_SIZE": 25, "PRIMER_MIN_SIZE": 22, "PRIMER_MAX_SIZE": 32,
            "PRIMER_OPT_TM": 66.0, "PRIMER_MIN_TM": 62.0, "PRIMER_MAX_TM": 72.0,
            "PRIMER_PAIR_MAX_DIFF_TM": 2.0,
            "PRIMER_MIN_GC": 40.0, "PRIMER_OPT_GC_PERCENT": 52.0, "PRIMER_MAX_GC": 65.0,
            "PRIMER_MAX_SELF_ANY_TH": 45.0, "PRIMER_MAX_SELF_END_TH": 35.0,
            "PRIMER_PAIR_MAX_COMPL_ANY_TH": 45.0, "PRIMER_PAIR_MAX_COMPL_END_TH": 35.0,
            "PRIMER_MAX_HAIRPIN_TH": 28.0,
            "PRIMER_MAX_END_STABILITY": 9.0,
            "PRIMER_GC_CLAMP": 2,
        },
    },
}

# Fields the UI may override, with sane bounds for validation.
TUNABLES = {
    "product_min":   ("Min product (bp)", 80, 20000),
    "product_max":   ("Max product (bp)", 120, 20000),
    "tm_min":        ("Min Tm (C)", 45.0, 80.0),
    "tm_opt":        ("Optimum Tm (C)", 45.0, 80.0),
    "tm_max":        ("Max Tm (C)", 45.0, 85.0),
    "gc_min":        ("Min GC (%)", 10.0, 80.0),
    "gc_max":        ("Max GC (%)", 20.0, 90.0),
    "size_min":      ("Min primer length", 15, 36),
    "size_opt":      ("Optimum primer length", 15, 36),
    "size_max":      ("Max primer length", 18, 36),
    "margin":        ("Min variant-to-primer distance (bp)", 0, 2000),
    "n_return":      ("Primer pairs to return", 1, 20),
}


def build(assay: str, overrides: dict | None = None) -> dict:
    """Resolve a preset plus user overrides into a concrete parameter set."""
    if assay not in PRESETS:
        raise ValueError(f"Unknown assay '{assay}'. Choose one of {sorted(PRESETS)}.")
    preset = PRESETS[assay]
    o = {k: v for k, v in (overrides or {}).items() if v not in (None, "")}

    product_min = int(o.get("product_min", preset["product"][0][0]))
    product_max = int(o.get("product_max", preset["product"][0][1]))
    if product_min >= product_max:
        raise ValueError("Minimum product size must be smaller than the maximum.")

    margin = int(o.get("margin", preset["margin"]))
    n_return = max(1, min(20, int(o.get("n_return", 5))))

    p3 = dict(BASE_P3)
    p3.update(preset["p3"])
    mapping = {
        "tm_min": "PRIMER_MIN_TM", "tm_opt": "PRIMER_OPT_TM", "tm_max": "PRIMER_MAX_TM",
        "gc_min": "PRIMER_MIN_GC", "gc_max": "PRIMER_MAX_GC",
        "size_min": "PRIMER_MIN_SIZE", "size_opt": "PRIMER_OPT_SIZE",
        "size_max": "PRIMER_MAX_SIZE",
    }
    for key, p3key in mapping.items():
        if key in o:
            p3[p3key] = float(o[key]) if "TM" in p3key or "GC" in p3key else int(o[key])
    if p3["PRIMER_MIN_TM"] > p3["PRIMER_MAX_TM"]:
        raise ValueError("Minimum Tm must not exceed maximum Tm.")
    if p3["PRIMER_MIN_SIZE"] > p3["PRIMER_MAX_SIZE"]:
        raise ValueError("Minimum primer length must not exceed maximum length.")

    p3["PRIMER_PRODUCT_SIZE_RANGE"] = [[product_min, product_max]]
    p3["PRIMER_PRODUCT_OPT_SIZE"] = int(o.get(
        "opt_product", min(max(preset["opt_product"], product_min), product_max)))
    p3["PRIMER_PAIR_WT_PRODUCT_SIZE_GT"] = 0.05
    p3["PRIMER_PAIR_WT_PRODUCT_SIZE_LT"] = 0.05
    p3["PRIMER_NUM_RETURN"] = n_return

    # Enough template to place the largest allowed product with room to spare.
    flank = max(preset["flank"], product_max + margin + 400)

    return {
        "assay": assay, "assay_name": preset["name"], "blurb": preset["blurb"],
        "flank": flank, "margin": margin, "n_return": n_return,
        "product_min": product_min, "product_max": product_max,
        "avoid_snps": bool((overrides or {}).get("avoid_snps", True)),
        "snp_maf": float((overrides or {}).get("snp_maf", 0.01)),
        "p3": p3,
    }
