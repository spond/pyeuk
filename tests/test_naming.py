"""
Unit tests for cyclospora_pyeuk.naming:
Bidirectional naming and parsing contract, round-trip locus invariants, and edge case handling.
"""

import unittest
from cyclospora_pyeuk.naming import parse_locus_name, name_haplotype, format_de_novo_haplotype_name


class TestNamingContract(unittest.TestCase):

    def test_cdc_legacy_parsing(self):
        """Test parsing of traditional CDC MLST marker column identifiers."""
        self.assertEqual(parse_locus_name("Nu_378_PART_A_Hap_4"), "Nu_378_PART_A")
        self.assertEqual(parse_locus_name("Nu_378_PART_B_Hap_1"), "Nu_378_PART_B")
        self.assertEqual(parse_locus_name("Nu_684_PART_A_Hap_10"), "Nu_684_PART_A")
        self.assertEqual(parse_locus_name("Mt_Cmt_Hap_1"), "Mt_Cmt")
        self.assertEqual(parse_locus_name("Mt_Junction_Hap_2"), "Mt_Cmt")
        self.assertEqual(parse_locus_name("Nu_130_PART_A_NOVEL_1"), "Nu_130_PART_A")

    def test_de_novo_and_amplicon_parsing(self):
        """Test parsing of de novo, windowed, and amplicon haplotype identifiers."""
        self.assertEqual(parse_locus_name("gp60_L752bp.H01_9180"), "gp60")
        self.assertEqual(parse_locus_name("18S_L830bp.H01_32A1"), "18S")
        self.assertEqual(parse_locus_name("COWP_L500bp.H02"), "COWP")
        self.assertEqual(parse_locus_name("ITS-2.H01"), "ITS-2")
        self.assertEqual(parse_locus_name("isotype-1_beta-tubulin_L276bp.H01_A1B2"), "isotype-1_beta-tubulin")
        self.assertEqual(parse_locus_name("PvAmpSeq_11_L120bp.H03"), "PvAmpSeq_11")
        self.assertEqual(parse_locus_name("Locus_01_L150bp.H02_ABCD"), "Locus_01")

    def test_round_trip_invariants(self):
        """
        Verify the round-trip invariant:
        parse_locus_name(name_haplotype(locus, ...)) == locus
        across various loci, lengths, hashes, ranks, and styles.
        """
        test_loci = [
            "Nu_378_PART_A",
            "Nu_684_PART_B",
            "gp60",
            "18S",
            "COWP",
            "Cp_HSP70",
            "ITS-2",
            "isotype-1_beta-tubulin",
            "PvAmpSeq_11",
            "Locus_01",
            "Marker_999",
        ]

        dummy_seq = "ATGCGATCGATCGATCGATCGATCGATC"

        for locus in test_loci:
            # 1. de_novo with length and hash
            col = name_haplotype(locus, hap_id=1, sequence=dummy_seq, style="de_novo")
            parsed = parse_locus_name(col)
            self.assertEqual(parsed, locus, f"Failed round-trip for de_novo: {col} -> {parsed} != {locus}")

            # 2. de_novo without hash
            col = name_haplotype(locus, hap_id=2, sequence=dummy_seq, hash_len=0, style="de_novo")
            parsed = parse_locus_name(col)
            self.assertEqual(parsed, locus, f"Failed round-trip for de_novo no-hash: {col} -> {parsed} != {locus}")

            # 3. cdc style
            col = name_haplotype(locus, hap_id=3, style="cdc")
            parsed = parse_locus_name(col)
            self.assertEqual(parsed, locus, f"Failed round-trip for cdc style: {col} -> {parsed} != {locus}")

            # 4. novel style
            col = name_haplotype(locus, hap_id=4, style="novel")
            parsed = parse_locus_name(col)
            self.assertEqual(parsed, locus, f"Failed round-trip for novel style: {col} -> {parsed} != {locus}")

            # 5. compact style
            col = name_haplotype(locus, hap_id=5, style="compact")
            parsed = parse_locus_name(col)
            self.assertEqual(parsed, locus, f"Failed round-trip for compact style: {col} -> {parsed} != {locus}")

    def test_format_de_novo_haplotype_name_compat(self):
        """Verify backward compatibility of format_de_novo_haplotype_name."""
        seq = "ATGCGATCGATC"
        name = format_de_novo_haplotype_name("18S", sequence=seq, rank=1, include_length=True, include_hash=True)
        self.assertTrue(name.startswith("18S_L12bp.H01_"))
        self.assertEqual(parse_locus_name(name), "18S")


if __name__ == "__main__":
    unittest.main()
