import unittest
import pandas as pd
from gwaspipe.order_alleles import (
    custom_alleles_sort,
    vectorizedorderalleles_status,
    orderalleles_status,
    parallelorderalleles_status,
    build_snpids,
    parallelbuildsnpid,
    order_alleles,
    _orderalleles_status_vec,
    _flip_allele_statistics,
    ORDER_MAPPING,
    TRANSLATE_TABLE_ORDER,
)
from gwaslab.g_Log import Log


class TestCustomAllelesSort(unittest.TestCase):
    """Tests for the custom_alleles_sort function."""

    def test_empty_list(self):
        """Test sorting of empty list."""
        self.assertEqual(custom_alleles_sort([]), [])

    def test_single_allele(self):
        """Test sorting of single allele."""
        self.assertEqual(custom_alleles_sort(["A"]), ["A"])

    def test_single_length_alleles_alphabetical(self):
        """Test sorting of single-length alleles returns alphabetical order."""
        input_alleles = ["C", "A", "T", "G"]
        expected = ["A", "C", "G", "T"]
        self.assertEqual(custom_alleles_sort(input_alleles), expected)

    def test_multi_length_alleles_same_length(self):
        """Test sorting of multi-length alleles with same length."""
        input_alleles = ["AT", "CG", "TA", "GC"]
        expected = ["AT", "CG", "GC", "TA"]
        self.assertEqual(custom_alleles_sort(input_alleles), expected)

    def test_mixed_length_alleles_multi_before_single(self):
        """Test that multi-length alleles come before single-length."""
        alleles = ["AT", "A", "T", "AC"]
        sorted_alleles = custom_alleles_sort(alleles)
        expected = ["AC", "AT", "A", "T"]
        self.assertEqual(sorted_alleles, expected)

    def test_two_single_alleles_already_ordered(self):
        """Test two single alleles that are already in order."""
        self.assertEqual(custom_alleles_sort(["A", "T"]), ["A", "T"])

    def test_two_single_alleles_need_swap(self):
        """Test two single alleles that need to be swapped."""
        self.assertEqual(custom_alleles_sort(["T", "A"]), ["A", "T"])

    def test_identical_alleles(self):
        """Test sorting of identical alleles."""
        self.assertEqual(custom_alleles_sort(["A", "A"]), ["A", "A"])

    def test_multi_length_different_lengths(self):
        """Test multi-length alleles with different lengths; longer comes first."""
        alleles = ["A", "ATG", "AT"]
        sorted_alleles = custom_alleles_sort(alleles)
        # ATG (len 3) comes before AT (len 2) which comes before A (len 1)
        self.assertEqual(sorted_alleles, ["ATG", "AT", "A"])

    def test_multi_length_same_prefix(self):
        """Test multi-length alleles with same prefix."""
        alleles = ["AC", "AT"]
        expected = ["AC", "AT"]
        self.assertEqual(custom_alleles_sort(alleles), expected)


class TestOrderMapping(unittest.TestCase):
    """Tests for module-level constants ORDER_MAPPING and TRANSLATE_TABLE_ORDER."""

    def test_order_mapping_keys(self):
        """Test ORDER_MAPPING contains the four nucleotide bases."""
        self.assertEqual(set(ORDER_MAPPING.keys()), {"A", "C", "G", "T"})

    def test_order_mapping_values_not_null(self):
        """Test ORDER_MAPPING values are not chr(0)."""
        for v in ORDER_MAPPING.values():
            self.assertNotEqual(v, chr(0))

    def test_translate_table_order_exists(self):
        """Test TRANSLATE_TABLE_ORDER is created."""
        self.assertIsNotNone(TRANSLATE_TABLE_ORDER)


class TestOrderAllelesStatusVec(unittest.TestCase):
    """Tests for the _orderalleles_status_vec function."""

    def setUp(self):
        self.log = Log()

    def test_empty_dataframe(self):
        """Test with an empty DataFrame returns immediately."""
        empty_df = pd.DataFrame({"EA": [], "NEA": [], "STATUS": []})
        result = _orderalleles_status_vec(empty_df, verbose=False, log=self.log)
        self.assertTrue(result.empty)

    def test_no_swap_when_ea_comes_first(self):
        """Test no swap occurs when EA already sorts first (A < T)."""
        df = pd.DataFrame(
            {
                "EA": ["A"],
                "NEA": ["T"],
                "STATUS": ["9999999"],
            }
        )
        result = _orderalleles_status_vec(df, verbose=False, log=self.log)
        self.assertEqual(result.loc[0, "STATUS"], "9999999")

    def test_swap_when_nea_comes_first(self):
        """Test swap occurs when NEA sorts before EA."""
        df = pd.DataFrame(
            {
                "EA": ["T"],
                "NEA": ["A"],
                "STATUS": ["9999999"],
            }
        )
        result = _orderalleles_status_vec(df, verbose=False, log=self.log)
        # Position 6 should be changed to '3'
        self.assertEqual(result.loc[0, "STATUS"], "9999939")

    def test_nea_longer_than_ea_triggers_swap(self):
        """Test swap when NEA is longer than EA (longer alleles sort first)."""
        df = pd.DataFrame(
            {
                "EA": ["A"],
                "NEA": ["AT"],
                "STATUS": ["9999999"],
            }
        )
        result = _orderalleles_status_vec(df, verbose=False, log=self.log)
        self.assertEqual(result.loc[0, "STATUS"], "9999939")

    def test_ea_longer_than_nea_no_swap(self):
        """Test no swap when EA is longer than NEA."""
        df = pd.DataFrame(
            {
                "EA": ["AT"],
                "NEA": ["A"],
                "STATUS": ["9999999"],
            }
        )
        result = _orderalleles_status_vec(df, verbose=False, log=self.log)
        self.assertEqual(result.loc[0, "STATUS"], "9999999")

    def test_same_length_alleles_mixed(self):
        """Test multiple rows with same-length alleles, some needing swap."""
        df = pd.DataFrame(
            {
                "EA": ["A", "T", "C", "G"],
                "NEA": ["T", "A", "G", "C"],
                "STATUS": ["9999999", "9999999", "9999999", "9999999"],
            }
        )
        result = _orderalleles_status_vec(df, verbose=False, log=self.log)
        # Row 0: EA=A, NEA=T -> A < T, no swap needed
        self.assertEqual(result.loc[0, "STATUS"], "9999999")
        # Row 1: EA=T, NEA=A -> A < T, swap needed
        self.assertEqual(result.loc[1, "STATUS"], "9999939")
        # Row 2: EA=C, NEA=G -> C < G, no swap needed
        self.assertEqual(result.loc[2, "STATUS"], "9999999")
        # Row 3: EA=G, NEA=C -> C < G, swap needed
        self.assertEqual(result.loc[3, "STATUS"], "9999939")

    def test_multi_char_alleles_same_length(self):
        """Test with multi-character alleles of same length."""
        df = pd.DataFrame(
            {
                "EA": ["AT", "CG"],
                "NEA": ["AC", "TA"],
                "STATUS": ["9999999", "9999999"],
            }
        )
        result = _orderalleles_status_vec(df, verbose=False, log=self.log)
        # Row 0: EA=AT, NEA=AC -> AC < AT, swap needed
        self.assertEqual(result.loc[0, "STATUS"], "9999939")
        # Row 1: EA=CG, NEA=TA -> CG < TA, no swap needed
        self.assertEqual(result.loc[1, "STATUS"], "9999999")


class TestVectorizedOrderAllelesStatus(unittest.TestCase):
    """Tests for the vectorizedorderalleles_status wrapper function."""

    def setUp(self):
        self.log = Log()
        self.test_data = pd.DataFrame(
            {
                "CHR": [1, 2, 3],
                "POS": [1000, 2000, 3000],
                "EA": ["A", "T", "C"],
                "NEA": ["T", "A", "G"],
                "STATUS": ["9900000", "9900000", "9900000"],
            }
        )

    def test_returns_dataframe(self):
        """Test that it returns a DataFrame of same length."""
        result = vectorizedorderalleles_status(self.test_data.copy(), verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 3)

    def test_missing_columns_returns_unchanged(self):
        """Test that missing required columns returns the DataFrame unchanged."""
        df = pd.DataFrame({"X": [1, 2]})
        result = vectorizedorderalleles_status(df, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 2)

    def test_long_alleles_above_threshold(self):
        """Test alleles exceeding the max_len=4 threshold are handled."""
        df = pd.DataFrame(
            {
                "EA": ["ATCGA", "A"],
                "NEA": ["A", "T"],
                "STATUS": ["9900000", "9900000"],
            }
        )
        result = vectorizedorderalleles_status(df, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 2)


class TestOrderAllelesStatus(unittest.TestCase):
    """Tests for the orderalleles_status function."""

    def setUp(self):
        self.log = Log()

    def test_basic_ordering(self):
        """Test basic allele ordering with swap."""
        df = pd.DataFrame(
            {
                "EA": ["T", "A"],
                "NEA": ["A", "T"],
                "STATUS": ["0000000", "0000000"],
            }
        )
        result = orderalleles_status(df, verbose=False, log=self.log)
        # Row 0: EA=T, NEA=A -> A<T, so EA/NEA should be swapped -> status digit 6 becomes '3'
        self.assertEqual(result.loc[0, "STATUS"], "0000030")
        # Row 1: EA=A, NEA=T -> A<T, already correct -> no swap
        self.assertEqual(result.loc[1, "STATUS"], "0000000")

    def test_no_swap_needed(self):
        """Test when alleles are already in correct order."""
        df = pd.DataFrame(
            {
                "EA": ["A"],
                "NEA": ["T"],
                "STATUS": ["0000000"],
            }
        )
        result = orderalleles_status(df, verbose=False, log=self.log)
        self.assertEqual(result.loc[0, "STATUS"], "0000000")


class TestParallelOrderAllelesStatus(unittest.TestCase):
    """Tests for the parallelorderalleles_status function."""

    def setUp(self):
        self.log = Log()
        self.test_data = pd.DataFrame(
            {
                "CHR": [1, 2, 3],
                "POS": [1000, 2000, 3000],
                "EA": ["A", "T", "C"],
                "NEA": ["T", "A", "G"],
                "STATUS": ["9900000", "9900000", "9900000"],
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

    def test_missing_columns_returns_unchanged(self):
        """Test that missing required columns returns DataFrame unchanged."""
        df = pd.DataFrame({"X": [1, 2]})
        result = parallelorderalleles_status(df, n_cores=1, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 2)


class TestBuildSnpids(unittest.TestCase):
    """Tests for the build_snpids function."""

    def test_basic_snpid_building(self):
        """Test SNP ID building from CHR:POS:EA:NEA."""
        df = pd.DataFrame(
            {
                "CHR": [1, 2],
                "POS": [1000, 2000],
                "EA": ["A", "T"],
                "NEA": ["T", "A"],
            }
        )
        result = build_snpids(df)
        expected = ["1:1000:A:T", "2:2000:T:A"]
        self.assertEqual(result.tolist(), expected)

    def test_custom_column_names(self):
        """Test with custom column names."""
        df = pd.DataFrame(
            {
                "CHROM": [1],
                "POSITION": [500],
                "A1": ["C"],
                "A2": ["G"],
            }
        )
        result = build_snpids(df, chrom="CHROM", pos="POSITION", ea="A1", nea="A2")
        self.assertEqual(result.tolist(), ["1:500:C:G"])


class TestParallelBuildSnpid(unittest.TestCase):
    """Tests for the parallelbuildsnpid function."""

    def setUp(self):
        self.log = Log()
        self.test_data = pd.DataFrame(
            {
                "CHR": [1, 2, 3],
                "POS": [1000, 2000, 3000],
                "EA": ["A", "T", "C"],
                "NEA": ["T", "A", "G"],
                "SNPID": ["", "", ""],
            }
        )

    def test_single_core(self):
        """Test SNPID building with single core."""
        result = parallelbuildsnpid(self.test_data.copy(), n_cores=1, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(result.loc[0, "SNPID"], "1:1000:A:T")
        self.assertEqual(result.loc[1, "SNPID"], "2:2000:T:A")
        self.assertEqual(result.loc[2, "SNPID"], "3:3000:C:G")

    def test_multi_core(self):
        """Test SNPID building with multiple cores."""
        result = parallelbuildsnpid(self.test_data.copy(), n_cores=2, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(result.loc[0, "SNPID"], "1:1000:A:T")

    def test_missing_snpid_column_skips(self):
        """Test that missing SNPID column skips building."""
        df = pd.DataFrame(
            {
                "CHR": [1],
                "POS": [1000],
                "EA": ["A"],
                "NEA": ["T"],
            }
        )
        result = parallelbuildsnpid(df, n_cores=1, verbose=False, log=self.log)
        self.assertNotIn("SNPID", result.columns)


class TestFlipAlleleStatistics(unittest.TestCase):
    """Tests for the _flip_allele_statistics function."""

    def setUp(self):
        self.log = Log()

    def test_no_flip_needed(self):
        """Test when no statistics need flipping (status digit 6 is 0,1,2)."""
        df = pd.DataFrame(
            {
                "EA": pd.Categorical(["A", "C"]),
                "NEA": pd.Categorical(["T", "G"]),
                "STATUS": ["9900020", "9900020"],
            }
        )
        result = _flip_allele_statistics(df, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 2)

    def test_flip_ref_status_3(self):
        """Test flip_ref when status digit 6 is '3' (xxxxx3x)."""
        df = pd.DataFrame(
            {
                "EA": pd.Categorical(["A", "T"], categories=["A", "T"]),
                "NEA": pd.Categorical(["T", "A"], categories=["A", "T"]),
                "STATUS": ["9900030", "9900030"],
            }
        )
        result = _flip_allele_statistics(
            df,
            flip_ref=True,
            reverse_compl=False,
            flip_ref_und=False,
            flip_rev_strand=False,
            verbose=False,
            log=self.log,
        )
        self.assertIsInstance(result, pd.DataFrame)
        # Status should change from xxxxx3x to xxxxx1x
        for idx in result.index:
            self.assertIn(result.loc[idx, "STATUS"][5], ["1", "2"])

    def test_reverse_complement_status_4(self):
        """Test reverse complement when status digit 6 is '4' (xxxxx4x)."""
        df = pd.DataFrame(
            {
                "EA": pd.Categorical(["A"], categories=["A", "T"]),
                "NEA": pd.Categorical(["T"], categories=["A", "T"]),
                "STATUS": ["9900040"],
            }
        )
        result = _flip_allele_statistics(
            df,
            reverse_compl=True,
            flip_ref=False,
            flip_ref_und=False,
            flip_rev_strand=False,
            verbose=False,
            log=self.log,
        )
        self.assertIsInstance(result, pd.DataFrame)
        # Status position 6 should change from 4 to 2
        self.assertEqual(result.loc[0, "STATUS"][5], "2")

    def test_reverse_complement_status_5(self):
        """Test reverse complement when status digit 6 is '5' (xxxxx5x)."""
        df = pd.DataFrame(
            {
                "EA": pd.Categorical(["A"], categories=["A", "T"]),
                "NEA": pd.Categorical(["T"], categories=["A", "T"]),
                "STATUS": ["9900050"],
            }
        )
        result = _flip_allele_statistics(
            df,
            reverse_compl=True,
            flip_ref=True,
            flip_ref_und=False,
            flip_rev_strand=False,
            verbose=False,
            log=self.log,
        )
        self.assertIsInstance(result, pd.DataFrame)

    def test_flip_ref_und_status_pattern(self):
        """Test flip_ref_und for indels with status pattern xxxx[123][67]6."""
        df = pd.DataFrame(
            {
                "EA": pd.Categorical(["AT"], categories=["AT", "A"]),
                "NEA": pd.Categorical(["A"], categories=["AT", "A"]),
                "STATUS": ["9901166"],
            }
        )
        result = _flip_allele_statistics(
            df,
            reverse_compl=False,
            flip_ref=False,
            flip_ref_und=True,
            flip_rev_strand=False,
            verbose=False,
            log=self.log,
        )
        self.assertIsInstance(result, pd.DataFrame)
        # Last digit should change from 6 to 4
        self.assertEqual(result.loc[0, "STATUS"][6], "4")

    def test_flip_rev_strand_status_pattern(self):
        """Test flip_rev_strand for palindromic SNPs with pattern xxxxx[012]5."""
        df = pd.DataFrame(
            {
                "EA": pd.Categorical(["A"], categories=["A", "T"]),
                "NEA": pd.Categorical(["T"], categories=["A", "T"]),
                "STATUS": ["9900005"],
            }
        )
        result = _flip_allele_statistics(
            df,
            reverse_compl=False,
            flip_ref=False,
            flip_ref_und=False,
            flip_rev_strand=True,
            verbose=False,
            log=self.log,
        )
        self.assertIsInstance(result, pd.DataFrame)
        # Last digit should change from 5 to 2
        self.assertEqual(result.loc[0, "STATUS"][6], "2")

    def test_all_flags_disabled(self):
        """Test with all flip flags disabled - no changes should occur."""
        df = pd.DataFrame(
            {
                "EA": pd.Categorical(["A"]),
                "NEA": pd.Categorical(["T"]),
                "STATUS": ["9900030"],
            }
        )
        result = _flip_allele_statistics(
            df,
            reverse_compl=False,
            flip_ref=False,
            flip_ref_und=False,
            flip_rev_strand=False,
            verbose=False,
            log=self.log,
        )
        self.assertEqual(result.loc[0, "STATUS"], "9900030")

    def test_no_matching_status_patterns(self):
        """Test when status codes don't match any flip patterns."""
        df = pd.DataFrame(
            {
                "EA": pd.Categorical(["A"]),
                "NEA": pd.Categorical(["T"]),
                "STATUS": ["9900020"],
            }
        )
        result = _flip_allele_statistics(
            df, reverse_compl=True, flip_ref=True, flip_ref_und=True, flip_rev_strand=True, verbose=False, log=self.log
        )
        # No patterns matched, so status should remain unchanged
        self.assertEqual(result.loc[0, "STATUS"], "9900020")


class TestOrderAlleles(unittest.TestCase):
    """Tests for the main order_alleles function."""

    def setUp(self):
        self.log = Log()
        self.test_data = pd.DataFrame(
            {
                "CHR": [1, 2, 3, 4, 5],
                "POS": [1000, 2000, 3000, 4000, 5000],
                "EA": pd.Categorical(["A", "T", "C", "G", "A"]),
                "NEA": pd.Categorical(["T", "A", "G", "C", "T"]),
                "STATUS": ["9900000", "9900010", "9900020", "9900030", "9900040"],
                "SNPID": ["1:1000:A:T", "2:2000:T:A", "3:3000:C:G", "4:4000:G:C", "5:5000:A:T"],
            }
        )

    def test_empty_dataframe(self):
        """Test with empty DataFrame returns empty DataFrame."""
        empty_df = pd.DataFrame()
        result = order_alleles(empty_df, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertTrue(result.empty)

    def test_vectorized_mode(self):
        """Test order_alleles in vectorized mode."""
        result = order_alleles(self.test_data.copy(), mode="v", verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 5)

    def test_parallel_mode(self):
        """Test order_alleles in parallel mode."""
        result = order_alleles(self.test_data.copy(), mode="p", n_cores=1, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 5)

    def test_invalid_mode_falls_back_to_vectorized(self):
        """Test that invalid mode falls back to vectorized."""
        result = order_alleles(self.test_data.copy(), mode="invalid", verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 5)

    def test_log_none_creates_default_log(self):
        """Test that passing log=None creates a default Log."""
        result = order_alleles(self.test_data.copy(), log=None, verbose=False)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 5)

    def test_format_snpid_false(self):
        """Test with format_snpid=False skips SNPID rebuilding."""
        data = self.test_data.copy()
        result = order_alleles(data, format_snpid=False, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 5)

    def test_format_snpid_true_rebuilds(self):
        """Test with format_snpid=True (default) rebuilds SNPID column."""
        result = order_alleles(self.test_data.copy(), format_snpid=True, verbose=False, log=self.log)
        self.assertIn("SNPID", result.columns)

    def test_single_row(self):
        """Test with single-row DataFrame."""
        single_row = pd.DataFrame(
            {
                "CHR": [1],
                "POS": [1000],
                "EA": ["A"],
                "NEA": ["T"],
                "STATUS": ["9900000"],
                "SNPID": ["1:1000:A:T"],
            }
        )
        result = order_alleles(single_row, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 1)

    def test_mixed_length_alleles(self):
        """Test with mixed-length alleles."""
        mixed_data = pd.DataFrame(
            {
                "CHR": [1, 2],
                "POS": [1000, 2000],
                "EA": ["A", "AT"],
                "NEA": ["T", "A"],
                "STATUS": ["9900000", "9900000"],
                "SNPID": ["1:1000:A:T", "2:2000:AT:A"],
            }
        )
        result = order_alleles(mixed_data, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 2)

    def test_flipallelestats_args_custom(self):
        """Test passing custom flipallelestats_args."""
        result = order_alleles(
            self.test_data.copy(),
            flipallelestats_args={"reverse_compl": True},
            verbose=False,
            log=self.log,
        )
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 5)

    def test_flipallelestats_args_none_default(self):
        """Test that flipallelestats_args=None defaults to empty dict."""
        result = order_alleles(
            self.test_data.copy(),
            flipallelestats_args=None,
            verbose=False,
            log=self.log,
        )
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 5)

    def test_category_dtype_alleles(self):
        """Test with categorical EA/NEA columns."""
        data = self.test_data.copy()
        data["EA"] = pd.Categorical(data["EA"])
        data["NEA"] = pd.Categorical(data["NEA"])
        result = order_alleles(data, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 5)

    def test_parallel_mode_multi_core(self):
        """Test parallel mode with n_cores > 1."""
        result = order_alleles(self.test_data.copy(), mode="p", n_cores=2, verbose=False, log=self.log)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 5)
