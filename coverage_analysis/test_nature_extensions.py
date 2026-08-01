import csv
import gzip
import tempfile
from pathlib import Path
from unittest import TestCase

from coverage_analysis.nature_propeptide import _observed_bits_for_candidates
from coverage_analysis.nature_reference import aggregate_label_coverage


class NatureExtensionTests(TestCase):
    def test_candidate_observation_bits_are_filtered_and_combined(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "peptides.tsv.gz"
            with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(("id", "sequence", "enzyme_bits"))
                writer.writerow((1, "KEEP", 1))
                writer.writerow((2, "DROP", 2))
                writer.writerow((3, "KEEP", 4))
            self.assertEqual(
                _observed_bits_for_candidates(path, {"KEEP"}), {"KEEP": 5}
            )

    def test_label_coverage_aggregation(self):
        designs = [("Combo", ("Trypsin",))]
        labels = {"K+R": frozenset("KR"), "L": frozenset("L")}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "labels.tsv.gz"
            with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(
                    (
                        "group_id",
                        "major_protein",
                        "length",
                        "combo_k_r_covered",
                        "combo_k_r_pct",
                        "combo_l_covered",
                        "combo_l_pct",
                    )
                )
                writer.writerow((1, "P1", 10, 5, 50.0, 2, 20.0))
                writer.writerow((2, "P2", 20, 10, 50.0, 8, 40.0))
            rows = aggregate_label_coverage(path, designs, labels)
            lookup = {row["label"]: row for row in rows}
            self.assertEqual(lookup["K+R"]["covered_residues"], 15)
            self.assertEqual(lookup["K+R"]["residue_weighted_pct"], 50.0)
            self.assertEqual(lookup["L"]["covered_residues"], 10)
            self.assertAlmostEqual(lookup["L"]["residue_weighted_pct"], 100 / 3)
