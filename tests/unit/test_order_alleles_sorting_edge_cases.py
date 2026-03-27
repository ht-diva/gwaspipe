import unittest

import pandas as pd
from gwaslab.info.g_Log import Log

from gwaspipe.order_alleles.sorting import custom_alleles_sort


class TestCustomAllelesSortEdgeCases(unittest.TestCase):
    """Tests for edge cases in custom_alleles_sort function."""

    def test_empty_list(self):
        """Test sorting of empty list."""
        self.assertEqual(custom_alleles_sort([]), [])

    def test_single_allele(self):
        """Test sorting of single allele."""
        self.assertEqual(custom_alleles_sort(["A"]), ["A"])

    def test_identical_alleles(self):
        """Test sorting of identical alleles."""
        self.assertEqual(custom_alleles_sort(["A", "A"]), ["A", "A"])

    def test_mixed_length_alleles(self):
        """Test mixed-length alleles with different lengths."""
        alleles = ["A", "ATG", "AT"]
        sorted_alleles = custom_alleles_sort(alleles)
        # ATG (len 3) comes before AT (len 2) which comes before A (len 1)
        self.assertEqual(sorted_alleles, ["ATG", "AT", "A"])

    def test_multi_length_same_prefix(self):
        """Test multi-length alleles with same prefix."""
        alleles = ["AC", "AT"]
        expected = ["AC", "AT"]
        self.assertEqual(custom_alleles_sort(alleles), expected)

    def test_all_multi_length_alleles(self):
        """Test with only multi-length alleles."""
        alleles = ["AT", "CG", "TA", "GC"]
        expected = ["AT", "CG", "GC", "TA"]
        self.assertEqual(custom_alleles_sort(alleles), expected)

    def test_all_single_length_alleles(self):
        """Test with only single-length alleles."""
        alleles = ["C", "A", "T", "G"]
        expected = ["A", "C", "G", "T"]
        self.assertEqual(custom_alleles_sort(alleles), expected)

    def test_very_long_alleles(self):
        """Test with very long alleles."""
        alleles = ["A", "ATGC", "ATG", "ATGCG"]
        sorted_alleles = custom_alleles_sort(alleles)
        # Longest comes first
        self.assertEqual(sorted_alleles[0], "ATGCG")
        self.assertEqual(sorted_alleles[-1], "A")

    def test_numeric_alleles(self):
        """Test with numeric alleles (should work with string comparison)."""
        alleles = ["1", "10", "2"]
        sorted_alleles = custom_alleles_sort(alleles)
        # Should sort lexicographically (10 comes before 1, 2 comes after)
        self.assertEqual(sorted_alleles, ["10", "1", "2"])

    def test_mixed_case_alleles(self):
        """Test with mixed case alleles."""
        alleles = ["a", "A", "t", "T"]
        sorted_alleles = custom_alleles_sort(alleles)
        # Should sort based on ASCII values
        self.assertEqual(sorted_alleles, ["A", "T", "a", "t"])

    def test_special_characters(self):
        """Test with special characters in alleles."""
        alleles = ["A-", "A", "AT", "A."]
        sorted_alleles = custom_alleles_sort(alleles)
        # Multi-length come first, then single-length
        self.assertTrue(all(len(a) > 1 for a in sorted_alleles[:2]))
        # Last two should be single-length (A) and multi-length (AT)
        self.assertEqual(sorted_alleles[2:], ["AT", "A"])


if __name__ == "__main__":
    unittest.main()
