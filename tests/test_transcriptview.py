"""The cDNA tab's rendering rules, copied from Ensembl's transcript sequence page.

These are the rules that were derived from the reference screenshots: which
consequence colours a base, which variant a base shows when several overlap,
and which residues the protein line turns red.
"""
from __future__ import annotations

import unittest

from primerforge import transcriptview as tv


def _transcript(seq: str = "ATGAAATTTTAA", first_pos: int = 101) -> dict:
    """A single-exon, entirely coding transcript on the forward strand."""
    cols = [{"pos": first_pos + i, "base": b, "number": 1, "cds": True, "tx": i}
            for i, b in enumerate(seq)]
    return {"exonic_columns": cols, "coding_columns": cols, "phase": 0,
            "utr5_end": None, "utr3_start": None}


def _snv(pos: int, ref: str, alt: str) -> dict:
    return {"start": pos, "end": pos, "allele_list": [ref, alt]}


class Consequences(unittest.TestCase):
    def test_severity_orders_most_severe_first(self):
        self.assertLess(tv.SEVERITY["stop_gained"], tv.SEVERITY["missense_variant"])
        self.assertLess(tv.SEVERITY["missense_variant"], tv.SEVERITY["synonymous_variant"])
        self.assertLess(tv.SEVERITY["synonymous_variant"], tv.SEVERITY["3_prime_UTR_variant"])

    def test_worst_of_several_terms(self):
        self.assertEqual(tv._worst(["synonymous_variant", "stop_gained", "missense_variant"]),
                         "stop_gained")

    def test_worst_ignores_terms_it_does_not_colour(self):
        self.assertEqual(tv._worst(["intergenic_variant", "missense_variant"]),
                         "missense_variant")
        self.assertIsNone(tv._worst(["intergenic_variant"]))

    def test_css_covers_every_consequence(self):
        css = tv.consequence_css()
        for term, colour, _label in tv.CONSEQUENCES:
            self.assertIn(tv._csq_class(term), css)
            self.assertIn(colour, css)

    def test_dark_backgrounds_get_white_ink(self):
        css = tv.consequence_css()
        self.assertIn(".csq-stop-gained{background:#ff0000;color:#fff}", css)
        # Light ones keep dark ink, and the even-exon blue when the base is in one.
        self.assertIn(".csq-missense-variant.tvc-exon-alt{color:#1044ee}", css)


class Iupac(unittest.TestCase):
    def test_two_alleles(self):
        self.assertEqual(tv._iupac({"alleles": "A/G"}, "A"), "R")
        self.assertEqual(tv._iupac({"alleles": "C/T"}, "C"), "Y")

    def test_reference_base_is_always_included(self):
        self.assertEqual(tv._iupac({"alleles": "G"}, "A"), "R")

    def test_three_alleles(self):
        self.assertEqual(tv._iupac({"alleles": "A/C/G"}, "A"), "V")

    def test_indel_has_no_single_base_code(self):
        self.assertEqual(tv._iupac({"alleles": "AAG/-"}, "A"), "*")


class Trim(unittest.TestCase):
    def test_shared_leading_base_is_dropped(self):
        self.assertEqual(tv._trim("AAAG", "A"), (1, "AAG", ""))

    def test_shared_trailing_bases_are_dropped(self):
        # GCC -> CC is the deletion of the leading G, once the shared CC is trimmed.
        self.assertEqual(tv._trim("GCC", "CC"), (0, "G", ""))

    def test_substitution_is_untouched(self):
        self.assertEqual(tv._trim("A", "G"), (0, "A", "G"))


class Translate(unittest.TestCase):
    def test_codons_to_residues(self):
        self.assertEqual(tv.translate("ATGAAATTTTAA"), "MKF*")

    def test_partial_codon_at_the_end_is_ignored(self):
        self.assertEqual(tv.translate("ATGAA"), "M")


class AnnotateVariants(unittest.TestCase):
    """meta["exonic_columns"] is annotated in place; each base keeps one variant."""

    def test_missense_colours_the_base_and_reddens_its_codon(self):
        meta = _transcript()                       # ATG AAA TTT TAA
        tv.annotate_variants(meta, [_snv(105, "A", "G")], strand=1)   # AAA -> AGA (K -> R)
        cols = meta["exonic_columns"]
        self.assertEqual(cols[4]["csq"], "missense_variant")
        self.assertEqual(cols[4]["code"], "R")
        # The whole codon goes red in the protein line, not just the changed base.
        self.assertEqual([c.get("aa_changed", False) for c in cols[3:6]], [True, True, True])
        self.assertFalse(cols[6].get("aa_changed", False))

    def test_synonymous_change_leaves_the_protein_alone(self):
        meta = _transcript()
        tv.annotate_variants(meta, [_snv(109, "T", "C")], strand=1)   # TTT -> TTC (F -> F)
        cols = meta["exonic_columns"]
        self.assertEqual(cols[8]["csq"], "synonymous_variant")
        self.assertEqual(cols[8]["code"], "Y")
        self.assertFalse(any(c.get("aa_changed") for c in cols))

    def test_stop_gained(self):
        meta = _transcript()
        tv.annotate_variants(meta, [_snv(104, "A", "T")], strand=1)   # AAA -> TAA (K -> *)
        self.assertEqual(meta["exonic_columns"][3]["csq"], "stop_gained")
        self.assertTrue(meta["exonic_columns"][3]["aa_changed"])

    def test_in_frame_deletion_marks_every_base_it_covers(self):
        meta = _transcript()
        tv.annotate_variants(meta, [{"start": 104, "end": 106,
                                     "allele_list": ["AAA", "-"]}], strand=1)
        cols = meta["exonic_columns"]
        for col in cols[3:6]:
            self.assertEqual(col["csq"], "inframe_deletion")
            self.assertEqual(col["code"], "*")       # an indel has no IUPAC code

    def test_shortest_variant_wins_on_a_shared_base(self):
        """An SNV inside a deletion keeps its own code and colour (Ensembl's draw order)."""
        meta = _transcript()
        tv.annotate_variants(meta, [
            {"start": 104, "end": 106, "allele_list": ["AAA", "-"]},   # 3 bases
            _snv(105, "A", "G"),                                        # 1 base
        ], strand=1)
        cols = meta["exonic_columns"]
        self.assertEqual(cols[4]["csq"], "missense_variant")
        self.assertEqual(cols[4]["code"], "R")
        # Its neighbours, covered only by the deletion, still show the deletion.
        self.assertEqual(cols[3]["csq"], "inframe_deletion")
        self.assertEqual(cols[5]["csq"], "inframe_deletion")

    def test_reverse_strand_alleles_are_complemented(self):
        meta = _transcript()
        # On the minus strand Ensembl reports alleles against the genome, so a
        # T/C there is an A/G on the transcript: the code is R, not Y.
        tv.annotate_variants(meta, [_snv(105, "T", "C")], strand=-1)
        self.assertEqual(meta["exonic_columns"][4]["code"], "R")

    def test_variant_outside_the_exon_is_ignored(self):
        meta = _transcript()
        tv.annotate_variants(meta, [_snv(999, "A", "G")], strand=1)
        self.assertFalse(any("csq" in c for c in meta["exonic_columns"]))


if __name__ == "__main__":
    unittest.main()
