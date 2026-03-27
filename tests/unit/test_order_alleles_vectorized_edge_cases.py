import unittest

import pandas as pd
from gwaslab.info.g_Log import Log

from gwaspipe.order_alleles.vectorized import (
    _orderalleles_status_vec,
    orderalleles_status,
    parallelorderalleles_status,
    vectorizedorderalleles_status,
)


class TestOrderAllelesStatusVecEdgeCases(unittest.TestCase):
    """Tests for edge cases in _orderalleles_status_vec function."""

    def setUp(self):
        self.log = Log()

    def test_empty_dataframe(self):
        """Test with an empty DataFrame returns immediately."""
        empty_df = pd.DataFrame({"EA": [], "NEA": [], "STATUS": []})
        result = _orderalleles_status_vec(empty_df, verbose=False, log=self.log)
        self.assertTrue(result.empty)

    def test_single_row_no_swap(self):
        """Test single row where no swap is needed."""
        df = pd.DataFrame(
            {
                "EA": ["A"],
                "NEA": ["T"],
                "STATUS": [9999999],
            }
        )
        result = _orderalleles_status_vec(df, verbose=False, log=self.log)
        self.assertEqual(result.loc[0, "STATUS"], 9999999)

    def test_single_row_swap_needed(self):
        """Test single row where swap is needed."""
        df = pd.DataFrame(
            {
                "EA": ["T"],
                "NEA": ["A"],
                "STATUS": [9999999],
            }
        )
        result = _orderalleles_status_vec(df, verbose=False, log=self.log)
        self.assertEqual(result.loc[0, "STATUS"], 9999939)

    def test_very_long_alleles(self):
        """Test with very long alleles."""
        df = pd.DataFrame(
            {
                "EA": ["A"],
                "NEA": ["ATGCGCGATCGATCGATCG"],
                "STATUS": [9999999],
            }
        )
        result = _orderalleles_status_vec(df, verbose=False, log=self.log)
        # NEA is longer, so swap should be triggered
        self.assertEqual(result.loc[0, "STATUS"], 9999939)

    def test_equal_length_different_alleles(self):
        """Test with equal length alleles that are different."""
        df = pd.DataFrame(
            {
                "EA": ["AT"],
                "NEA": ["CG"],
                "STATUS": [9999999],
            }
        )
        result = _orderalleles_status_vec(df, verbose=False, log=self.log)
        # AT comes before CG alphabetically, no swap
        self.assertEqual(result.loc[0, "STATUS"], 9999999)

    def test_equal_length_same_alleles(self):
        """Test with equal length alleles that are the same."""
        df = pd.DataFrame(
            {
                "EA": ["AT"],
                "NEA": ["AT"],
                "STATUS": [9999999],
            }
        )
        result = _orderalleles_status_vec(df, verbose=False, log=self.log)
        # Alleles are identical, no swap
        self.assertEqual(result.loc[0, "STATUS"], 9999999)

    def test_mixed_lengths_in_dataframe(self):
        """Test DataFrame with mixed allele lengths."""
        df = pd.DataFrame(
            {
                "EA": ["A", "AT", "ATGC"],
                "NEA": ["T", "CG", "AT"],
                "STATUS": [9999999, 9999999, 9999999],
            }
        )
        result = _orderalleles_status_vec(df, verbose=False, log=self.log)
        # Row 0: A < T, no swap
        self.assertEqual(result.loc[0, "STATUS"], 9999999)
        # Row 1: AT < CG, no swap
        self.assertEqual(result.loc[1, "STATUS"], 9999999)
        # Row 2: ATGC > AT (length), but both start with AT, so no swap
        self.assertEqual(result.loc[2, "STATUS"], 9999999)


class TestVectorizedOrderAllelesStatusEdgeCases(unittest.TestCase):
    """Tests for edge cases in vectorizedorderalleles_status function."""

    def setUp(self):
        self.log = Log()

    def test_all_alleles_above_threshold(self):
        """Test when all alleles exceed the max_len threshold."""
        df = pd.DataFrame(
            {
                "EA": ["ATCGA", "ATCGAT"],
                "NEA": ["A", "AT"],
                "STATUS": [9900000, 9900000],
            }
        )
        result = vectorizedorderalleles_status(df, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 2)

    def test_all_alleles_below_threshold(self):
        """Test when all alleles are below the max_len threshold."""
        df = pd.DataFrame(
            {
                "EA": ["A", "T"],
                "NEA": ["T", "A"],
                "STATUS": [9900000, 9900000],
            }
        )
        result = vectorizedorderalleles_status(df, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 2)

    def test_mixed_threshold_alleles(self):
        """Test with mix of alleles above and below threshold."""
        df = pd.DataFrame(
            {
                "EA": ["A", "ATCGA"],
                "NEA": ["T", "AT"],
                "STATUS": [9900000, 9900000],
            }
        )
        result = vectorizedorderalleles_status(df, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 2)


class TestOrderAllelesStatusEdgeCases(unittest.TestCase):
    """Tests for edge cases in orderalleles_status function."""

    def setUp(self):
        self.log = Log()

    def test_empty_dataframe(self):
        """Test with empty DataFrame."""
        empty_df = pd.DataFrame({"EA": [], "NEA": [], "STATUS": []})
        result = orderalleles_status(empty_df, verbose=False, log=self.log)
        self.assertTrue(result.empty)

    def test_single_row_no_swap(self):
        """Test single row where no swap is needed."""
        df = pd.DataFrame(
            {
                "EA": ["A"],
                "NEA": ["T"],
                "STATUS": [9999999],
            }
        )
        result = orderalleles_status(df, verbose=False, log=self.log)
        self.assertEqual(result.loc[0, "STATUS"], 9999999)

    def test_single_row_swap_needed(self):
        """Test single row where swap is needed."""
        df = pd.DataFrame(
            {
                "EA": ["T"],
                "NEA": ["A"],
                "STATUS": [9999999],
            }
        )
        result = orderalleles_status(df, verbose=False, log=self.log)
        self.assertEqual(result.loc[0, "STATUS"], 9999939)


class TestParallelOrderAllelesStatusEdgeCases(unittest.TestCase):
    """Tests for edge cases in parallelorderalleles_status function."""

    def setUp(self):
        self.log = Log()
        self.test_data = pd.DataFrame(
            {
                "CHR": [1, 2, 3],
                "POS": [1000, 2000, 3000],
                "EA": ["A", "T", "C"],
                "NEA": ["T", "A", "G"],
                "STATUS": [9999999, 9999999, 9999999],
            }
        )

    def test_single_core(self):
        """Test parallel ordering with single core."""
        result = parallelorderalleles_status(self.test_data.copy(), n_cores=1, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 3)

    def test_multi_core(self):
        """Test parallel ordering with multiple cores."""
        result = parallelorderalleles_status(self.test_data.copy(), n_cores=2, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 3)


if __name__ == "__main__":
    unittest.main()
