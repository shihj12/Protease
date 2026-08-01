from types import SimpleNamespace
from unittest import TestCase
from unittest.mock import patch

from ase_assay.presets import PRESETS
from coverage_analysis.propeptide_quantification import (
    ProcessingJunction,
    junction_candidates,
    processing_junctions,
)


class FakeEntry(SimpleNamespace):
    def of_kind(self, kind):
        return [feature for feature in self.features if feature.kind == kind]


def feature(kind, start, end, description=""):
    return SimpleNamespace(
        kind=kind,
        start=start,
        end=end,
        description=description,
        evidence="experimental",
    )


class PropeptideQuantificationTests(TestCase):
    def test_pgrn_is_in_presets(self):
        self.assertIn("GRN", PRESETS["Lysosomal hydrolases"])

    def test_pgrn_products_add_internal_boundaries_but_not_signal_boundary(self):
        entry = FakeEntry(
            accession="P28799",
            gene="GRN",
            protein="Progranulin",
            sequence="A" * 120,
            length=120,
            features=[
                feature("SIGNAL", 1, 17),
                feature("CHAIN", 18, 120, "Progranulin"),
                feature("PEPTIDE", 18, 47, "Paragranulin"),
                feature("PEPTIDE", 58, 113, "Granulin-1"),
            ],
        )
        observed = [(j.position, j.side) for j in processing_junctions(entry)]
        self.assertEqual(observed, [(47, "C"), (57, "N"), (113, "C")])

    @patch(
        "coverage_analysis.propeptide_quantification.proteome.is_unique",
        return_value=True,
    )
    def test_exact_digest_cut_allows_intact_missed_but_not_mature(self, _unique):
        entry = FakeEntry(sequence="AAAAAAKAAAAAAA")
        junction = ProcessingJunction(
            "P1", "GENE", "Protein", 7, "N", 8, 14, "Mature", "propeptide", ""
        )
        intact, mature = junction_candidates(entry, junction, "Trypsin")
        self.assertTrue(intact)
        self.assertTrue(all(peptide.missed_cleavages >= 1 for peptide in intact))
        self.assertEqual(mature, [])

    @patch(
        "coverage_analysis.propeptide_quantification.proteome.is_unique",
        return_value=True,
    )
    def test_nonmatching_digest_cut_can_report_both_states(self, _unique):
        entry = FakeEntry(sequence="AAAAAAQAAAAAAK")
        junction = ProcessingJunction(
            "P1", "GENE", "Protein", 7, "N", 8, 14, "Mature", "propeptide", ""
        )
        intact, mature = junction_candidates(entry, junction, "Trypsin")
        self.assertTrue(intact)
        self.assertEqual([peptide.sequence for peptide in mature], ["AAAAAAK"])
