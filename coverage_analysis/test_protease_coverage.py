from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from coverage_analysis.protease_coverage import (
    LABELS,
    RULES,
    Protein,
    cleavage_boundaries,
    digest_mask,
    interval_mask,
    peptide_intervals,
    read_fasta,
)


class DigestionTests(unittest.TestCase):
    def test_trypsin_is_blocked_by_proline(self):
        self.assertEqual(cleavage_boundaries("AKPRA", RULES["Trypsin"]), [0, 4, 5])

    def test_n_terminal_rules(self):
        self.assertEqual(cleavage_boundaries("ADKRA", RULES["Asp-N"]), [0, 1, 5])
        self.assertEqual(cleavage_boundaries("ADKRA", RULES["Lys-N"]), [0, 2, 5])
        self.assertEqual(cleavage_boundaries("ADKRA", RULES["Arg-N"]), [0, 3, 5])

    def test_missed_cleavages_generate_joined_peptide(self):
        intervals = list(
            peptide_intervals("AAAAKAAAAKAAAA", RULES["Lys-C"], 1, 100, 1)
        )
        self.assertIn((0, 10), intervals)
        self.assertNotIn((0, 14), intervals)

    def test_label_mask_only_covers_label_bearing_peptide(self):
        mask, labelled, count = digest_mask(
            "AAAAKLLLL",
            RULES["Lys-C"],
            LABELS,
            min_length=1,
            max_length=20,
            missed_cleavages=0,
        )
        self.assertEqual(mask, interval_mask(0, 9))
        self.assertEqual(labelled["K+R"], interval_mask(0, 5))
        self.assertEqual(labelled["L"], interval_mask(5, 9))
        self.assertEqual(count, 2)

    def test_fasta_parser(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "test.fasta"
            path.write_text(
                ">sp|P1|ONE\nAAAA\n>ENSP2 description\nKK\n", encoding="utf-8"
            )
            proteins = read_fasta(path)
        self.assertEqual(proteins, [Protein("P1", "AAAA"), Protein("ENSP2", "KK")])


if __name__ == "__main__":
    unittest.main()

