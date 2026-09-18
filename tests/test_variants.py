"""Input parsing and HGVS classification - the first thing a user's text meets."""
from __future__ import annotations

import unittest

from primerforge import variants
from primerforge.variants import VariantError


class ParseInputBlock(unittest.TestCase):
    def test_splits_on_newlines_and_semicolons(self):
        text = "NM_000546.6:c.215C>G\nrs121913343;17-7676154-G-C"
        self.assertEqual(variants.parse_input_block(text),
                         ["NM_000546.6:c.215C>G", "rs121913343", "17-7676154-G-C"])

    def test_ignores_blank_lines_and_comments(self):
        self.assertEqual(variants.parse_input_block("# a note\n\n  rs334  \n"), ["rs334"])

    def test_comma_does_not_split_a_single_hgvs_description(self):
        # One colon and no comma: the line is one variant even though HGVS
        # ranges contain underscores and digits that look separable.
        self.assertEqual(variants.parse_input_block("NM_000546.6:c.1521_1523del"),
                         ["NM_000546.6:c.1521_1523del"])

    def test_comma_separates_plain_identifiers(self):
        self.assertEqual(variants.parse_input_block("rs334, rs1801133"), ["rs334", "rs1801133"])


class Classify(unittest.TestCase):
    def test_rsid_is_lowercased(self):
        self.assertEqual(variants.classify("RS121913343"), ("rsid", "rs121913343"))

    def test_coding_hgvs(self):
        kind, cleaned = variants.classify("NM_000546.6:c.215C>G")
        self.assertEqual((kind, cleaned), ("hgvs_c", "NM_000546.6:c.215C>G"))

    def test_genomic_hgvs(self):
        kind, cleaned = variants.classify("NC_000017.11:g.7676154G>C")
        self.assertEqual((kind, cleaned), ("hgvs_g", "NC_000017.11:g.7676154G>C"))

    def test_vcf_coordinates(self):
        self.assertEqual(variants.classify("17-7676154-G-C"), ("vcf", "17-7676154-G-C"))

    def test_unversioned_transcript_is_refused_with_advice(self):
        with self.assertRaises(VariantError) as caught:
            variants.classify("NM_000546:c.215C>G")
        self.assertIn("version", str(caught.exception))

    def test_gene_symbol_alone_is_refused(self):
        with self.assertRaises(VariantError) as caught:
            variants.classify("TP53")
        self.assertIn("variant description", str(caught.exception))

    def test_uppercase_hgvs_type_is_accepted_and_normalised(self):
        # Pasting from a paper often capitalises the type; the cleaned form is
        # what later lookups use, so it must come back lower case.
        self.assertEqual(variants.classify("NM_000546.6:C.215C>G"),
                         ("hgvs_c", "NM_000546.6:c.215C>G"))

    def test_incomplete_description_is_explained(self):
        with self.assertRaises(VariantError) as caught:
            variants.classify("NM_000546.6:c.215")
        self.assertIn("not a recognised HGVS change", str(caught.exception))

    def test_empty_input(self):
        with self.assertRaises(VariantError):
            variants.classify("   ")


class Revcomp(unittest.TestCase):
    def test_reverse_complement(self):
        self.assertEqual(variants.revcomp("ACGTN"), "NACGT")


if __name__ == "__main__":
    unittest.main()
