import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import gwaslab as gl
import pandas as pd

from gwaspipe.gwaspipe import SumstatsManager


class TestSumstatsManager(unittest.TestCase):
    """Tests for the SumstatsManager class."""

    def setUp(self):
        self.test_data_path = Path("tests/data/test_sumstats.tsv")
        self.formatbook_path = Path("data/formatbook.json")
        self.pid = False
        self.bcfliftover = False

    def test_init_with_plink_pvar_format(self):
        """Test initialization with plink_pvar format."""
        sm = SumstatsManager(
            input_path=str(self.test_data_path),
            input_format="plink_pvar",
            input_separator="\t",
            input_study=None,
            formatbook_path=self.formatbook_path,
            pid=self.pid,
            bcfliftover=self.bcfliftover,
        )
        self.assertIsNotNone(sm.mysumstats)
        self.assertIsInstance(sm.mysumstats, gl.Sumstats)

    def test_init_with_pickle_format(self):
        """Test initialization with pickle format."""
        # Create a minimal pickle file for testing
        test_df = pd.DataFrame(
            {
                "CHR": [1],
                "POS": [1000],
                "EA": ["A"],
                "NEA": ["T"],
                "BETA": [0.5],
                "SE": [0.1],
                "P": [1e-5],
            }
        )
        sumstats = gl.Sumstats(test_df)
        pickle_path = Path("tests/data/test.pkl")
        gl.dump_pickle(sumstats, str(pickle_path))

        sm = SumstatsManager(
            input_path=str(pickle_path),
            input_format="pickle",
            input_separator="\t",
            input_study=None,
            formatbook_path=self.formatbook_path,
            pid=self.pid,
            bcfliftover=self.bcfliftover,
        )
        self.assertIsNotNone(sm.mysumstats)
        pickle_path.unlink()

    def test_make_gwaslab_snpid(self):
        """Test _make_gwaslab_snpid method."""
        # Use a test file that has CHR, POS, EA, NEA columns
        test_data_path = Path("tests/data/test_with_chr_pos.pkl")
        sm = SumstatsManager(
            input_path=str(test_data_path),
            input_format="pickle",
            input_separator="\t",
            input_study=None,
            formatbook_path=self.formatbook_path,
            pid=self.pid,
            bcfliftover=self.bcfliftover,
        )
        snpid = sm._make_gwaslab_snpid()
        self.assertIsInstance(snpid, pd.Series)
        self.assertTrue(all(snpid.str.contains(":")))

    def test_float_dict_custom(self):
        """Test float_dict_custom method."""
        # Use a test file that has BETA column
        test_data_path = Path("tests/data/test_with_beta.pkl")
        sm = SumstatsManager(
            input_path=str(test_data_path),
            input_format="pickle",
            input_separator="\t",
            input_study=None,
            formatbook_path=self.formatbook_path,
            pid=self.pid,
            bcfliftover=self.bcfliftover,
        )
        gp = {"float_formats": {"BETA": ".3f"}}
        float_dict = sm.float_dict_custom(gp)
        self.assertIsInstance(float_dict, dict)
        self.assertIn("BETA", float_dict)

    def test_order_alleles(self):
        """Test order_alleles method."""
        sm = SumstatsManager(
            input_path=str(self.test_data_path),
            input_format="plink_pvar",
            input_separator="\t",
            input_study=None,
            formatbook_path=self.formatbook_path,
            pid=self.pid,
            bcfliftover=self.bcfliftover,
        )
        initial_shape = sm.mysumstats.data.shape
        sm.order_alleles(n_cores=1, mode="v")
        self.assertEqual(sm.mysumstats.data.shape, initial_shape)


class TestCLIOptions(unittest.TestCase):
    """Tests for CLI option parsing."""

    @patch("gwaspipe.gwaspipe.ConfigurationManager")
    @patch("gwaspipe.gwaspipe.SumstatsManager")
    def test_main_with_minimal_options(self, mock_sm_class, mock_cm_class):
        """Test main function with minimal required options."""
        mock_cm = MagicMock()
        mock_cm.log_file_path = Path("test.log")
        mock_cm.formatbook_path = Path("data/formatbook.json")
        mock_cm.run_sequence = ()
        mock_cm.filename_settings = (None, None)
        mock_cm_class.return_value = mock_cm

        mock_sm = MagicMock()
        mock_sm_class.return_value = mock_sm

        import sys
        from io import StringIO

        from gwaspipe.gwaspipe import main

        # Test with minimal options
        testargs = [
            "gwaspipe",
            "-c",
            "tests/data/test_config.yaml",
            "-i",
            "tests/data/test_sumstats.tsv",
            "-f",
            "plink_pvar",
            "-o",
            "results",
        ]

        with patch.object(sys, "argv", testargs):
            try:
                main()
            except SystemExit:
                pass

        mock_cm_class.assert_called_once()


if __name__ == "__main__":
    unittest.main()
