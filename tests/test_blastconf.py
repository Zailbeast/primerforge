"""The BLAST/BLAT form: how pasted text becomes searchable sequences."""
from __future__ import annotations

import unittest

from primerforge import blastconf


class ParseSequences(unittest.TestCase):
    def test_fasta_records_split_on_the_header(self):
        result = blastconf.parse_sequences(">one\nACGTACGT\n>two\nTTTTGGGG\n")
        self.assertEqual([s["description"] for s in result["sequences"]], ["one", "two"])
        self.assertEqual(result["sequences"][0]["sequence"], "ACGTACGT")

    def test_blank_line_starts_a_new_sequence(self):
        result = blastconf.parse_sequences("ACGTACGT\n\nTTTTGGGG")
        self.assertEqual(len(result["sequences"]), 2)

    def test_digits_and_spaces_are_dropped(self):
        # Sequence copied from a viewer arrives with position numbers in it.
        result = blastconf.parse_sequences("1 ACGT ACGT 8")
        self.assertEqual(result["sequences"][0]["sequence"], "ACGTACGT")

    def test_rna_is_converted_to_dna_and_flagged(self):
        result = blastconf.parse_sequences("ACGUACGU")
        self.assertTrue(result["is_rna"])
        self.assertEqual(result["sequences"][0]["sequence"], "ACGTACGT")

    def test_duplicates_are_dropped(self):
        result = blastconf.parse_sequences(">a\nACGTACGT\n>b\nACGTACGT\n")
        self.assertEqual(len(result["sequences"]), 1)
        self.assertEqual(result["invalids"], 1)

    def test_invalid_characters_are_reported_not_silently_dropped(self):
        result = blastconf.parse_sequences("ACGT!!!@@@")
        self.assertEqual(result["sequences"], [])
        self.assertTrue(result["errors"])
        self.assertIn("invalid characters", result["errors"][0])

    def test_peptide_is_detected(self):
        result = blastconf.parse_sequences(">p\nMKWVTFISLLFLFSSAYSRGVFRR\n")
        self.assertEqual(result["sequences"][0]["type"], "peptide")

    def test_single_line_fasta_pasted_from_a_web_page(self):
        # A copied record can arrive with its newline turned into a space.
        result = blastconf.parse_sequences(">gene one ACGTACGTACGT")
        self.assertEqual(len(result["sequences"]), 1)
        self.assertEqual(result["sequences"][0]["sequence"], "ACGTACGTACGT")

    def test_respects_the_sequence_limit(self):
        # Distinct sequences: identical ones are dropped as duplicates, and digits
        # are stripped from the input, so the numbering has to be in the bases.
        many = "\n\n".join("ACGTACGTAC" + "".join("ACGT"[(i >> s) & 3] for s in range(6))
                           for i in range(blastconf.MAX_NUM_SEQUENCES + 5))
        result = blastconf.parse_sequences(many)
        self.assertEqual(len(result["sequences"]), blastconf.MAX_NUM_SEQUENCES)

    def test_empty_input(self):
        self.assertEqual(blastconf.parse_sequences("")["sequences"], [])


class QueryType(unittest.TestCase):
    def test_mixed_input_follows_the_majority(self):
        seqs = [{"type": "dna"}, {"type": "dna"}, {"type": "peptide"}]
        self.assertEqual(blastconf.guess_query_type(seqs), "dna")
        self.assertEqual(blastconf.guess_query_type([{"type": "peptide"}] * 3), "peptide")


class Fasta(unittest.TestCase):
    def test_wraps_and_names(self):
        text = blastconf.fasta("hit", "A" * 130)
        lines = text.strip().split("\n")
        self.assertEqual(lines[0], ">hit")
        self.assertEqual([len(ln) for ln in lines[1:]], [60, 60, 10])

    def test_unnamed_sequence_gets_a_placeholder(self):
        self.assertTrue(blastconf.fasta("", "ACGT").startswith(">query"))


class BlastArgs(unittest.TestCase):
    def test_defaults_round_trip_through_normalise(self):
        defaults = blastconf.defaults_for("NCBIBLAST_BLASTN")
        clean = blastconf.normalise_configs("NCBIBLAST_BLASTN", defaults)
        args = blastconf.blast_args("NCBIBLAST_BLASTN", clean)
        self.assertEqual(args[:2], ["-task", "blastn"])
        self.assertIn("-evalue", args)

    def test_an_unknown_value_is_refused(self):
        with self.assertRaises(blastconf.ConfigError):
            blastconf.normalise_configs("NCBIBLAST_BLASTN", {"evalue": "not-a-number"})


if __name__ == "__main__":
    unittest.main()
