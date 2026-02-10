import unittest
import pandas as pd
from gwaspipe.util_order_alleles import (
    custom_alleles_sort,
    vectorizedorderalleles_status,
    orderalleles_status,
    parallelorderalleles_status,
    build_snpids,
    parallelbuildsnpid,
    order_alleles,
)
from gwaslab.g_Log import Log


class TestUtilOrderAlleles(unittest.TestCase):
    def setUp(self):
        """Set up test data for all tests."""
        self.test_data = pd.DataFrame(
            {
                "CHR": [1, 2, 3, 4, 5],
                "POS": [1000, 2000, 3000, 4000, 5000],
                "EA": ["A", "T", "C", "G", "A"],
                "NEA": ["T", "A", "G", "C", "T"],
                "STATUS": ["000000", "000001", "000002", "000003", "000004"],
                "SNPID": ["1:1000:A:T", "2:2000:T:A", "3:3000:C:G", "4:4000:G:C", "5:5000:A:T"],
            }
        )
        self.log = Log()

    def test_empty_list(self):
        """Test sorting of empty list"""
        self.assertEqual(custom_alleles_sort([]), [])

    def test_single_allele(self):
        """Test sorting of single allele"""
        self.assertEqual(custom_alleles_sort(["A"]), ["A"])

    def test_single_length_alleles(self):
        """Test sorting of single-length alleles"""
        input_alleles = ["C", "A", "T", "G"]
        expected = ["A", "C", "G", "T"]
        self.assertEqual(custom_alleles_sort(input_alleles), expected)

    def test_multi_length_alleles_same_length(self):
        """Test sorting of multi-length alleles with same length"""
        input_alleles = ["AT", "CG", "TA", "GC"]
        expected = ["AT", "CG", "GC", "TA"]
        self.assertEqual(custom_alleles_sort(input_alleles), expected)

    # def test_custom_alleles_sort(self):
    #     """Test custom allele sorting function."""
    #     alleles = ['A', 'T', 'C', 'G']
    #     sorted_alleles = custom_alleles_sort(alleles)
    #     expected = ['A', 'C', 'G', 'T']
    #     self.assertEqual(sorted_alleles, expected)
    #
    #     # Test with mixed length alleles
    #     alleles = ['AT', 'A', 'T', 'AC']
    #     sorted_alleles = custom_alleles_sort(alleles)
    #     expected = ['A', 'AC', 'AT', 'T']
    #     self.assertEqual(sorted_alleles, expected)
    #
    # def test_orderalleles_status_vec(self):
    #     """Test vectorized allele ordering."""
    #     # Create a simple test case
    #     test_df = pd.DataFrame({
    #         'EA': ['A', 'T'],
    #         'NEA': ['T', 'A'],
    #         'STATUS': ['000000', '000000']
    #     })
    #     print(test_df)
    #
    #     result = _orderalleles_status_vec(test_df, verbose=False, log=self.log)
    #     print(result)
    #     # Should not swap alleles for first row since T < A
    #     self.assertEqual(result.loc[0, 'STATUS'], '000000')
    #     # Should swap for second row since A < T
    #     self.assertEqual(result.loc[1, 'STATUS'], '000003')
    #
    def test_vectorizedorderalleles_status(self):
        """Test vectorized ordering function."""
        result = vectorizedorderalleles_status(self.test_data.copy(), verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 5)

    def test_orderalleles_status(self):
        """Test basic allele ordering function."""
        result = orderalleles_status(self.test_data.copy(), verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 5)

    def test_parallelorderalleles_status(self):
        """Test parallel allele ordering function."""
        result = parallelorderalleles_status(self.test_data.copy(), n_cores=1, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 5)

    def test_build_snpids(self):
        """Test SNP ID building function."""
        result = build_snpids(self.test_data)
        expected = ["1:1000:A:T", "2:2000:T:A", "3:3000:C:G", "4:4000:G:C", "5:5000:A:T"]
        self.assertEqual(result.tolist(), expected)

    def test_parallelbuildsnpid(self):
        """Test parallel SNP ID building function."""
        result = parallelbuildsnpid(self.test_data.copy(), n_cores=1, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 5)

    def test_order_alleles(self):
        """Test main order_alleles function."""
        # Test vectorized mode
        result_v = order_alleles(self.test_data.copy(), mode="v", verbose=False, log=self.log)
        self.assertIsInstance(result_v, pd.DataFrame)
        self.assertEqual(len(result_v), 5)

        # Test parallel mode
        result_p = order_alleles(self.test_data.copy(), mode="p", n_cores=1, verbose=False, log=self.log)
        self.assertIsInstance(result_p, pd.DataFrame)
        self.assertEqual(len(result_p), 5)

    def test_order_alleles_edge_cases(self):
        """Test edge cases in order_alleles."""
        # Test with empty dataframe
        empty_df = pd.DataFrame()
        result = order_alleles(empty_df, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)

        # Test with single row
        single_row = pd.DataFrame({"CHR": [1], "POS": [1000], "EA": ["A"], "NEA": ["T"], "STATUS": ["000000"]})
        result = order_alleles(single_row, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 1)

    def test_mixed_length_alleles(self):
        """Test with mixed length alleles."""
        mixed_data = pd.DataFrame(
            {"CHR": [1, 2], "POS": [1000, 2000], "EA": ["A", "AT"], "NEA": ["T", "A"], "STATUS": ["000000", "000000"]}
        )

        result = order_alleles(mixed_data, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 2)
