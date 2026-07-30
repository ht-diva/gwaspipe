import json
import unittest
from datetime import UTC, datetime
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import MagicMock, patch

import gwaslab as gl
import pandas as pd

from gwaspipe import __version__
from gwaspipe.gwaspipe import (
    AssemblyValidationError,
    SumstatsManager,
    _require_validated_assembly,
    _write_run_provenance,
    validate_declared_assembly,
)


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

    @patch("gwaspipe.gwaspipe.ConfigurationManager")
    @patch("gwaspipe.gwaspipe.SumstatsManager")
    def test_main_forwards_canonicalization_parameters_for_both_step_names(self, mock_sm_class, mock_cm_class):
        """The preferred canonicalization step and its alias forward identical parameters."""
        sort_params = {"mode": "p", "n_cores": 2, "format_snpid": False, "verbose": False}
        mock_cm = MagicMock()
        mock_cm.log_file_path = Path("test.log")
        mock_cm.formatbook_path = Path("data/formatbook.json")
        mock_cm.filename_settings = (None, None)
        mock_cm.step.return_value = ({"run": True, "workspace": "default"}, sort_params)
        mock_cm_class.return_value = mock_cm

        from click.testing import CliRunner

        from gwaspipe.gwaspipe import main

        with TemporaryDirectory() as temporary_directory:
            mock_cm.root_path = temporary_directory
            for step_name in ("canonicalize_effect_alleles", "sort_alphabetically"):
                with self.subTest(step_name=step_name):
                    mock_cm.run_sequence = (step_name,)
                    mock_sm = MagicMock()
                    mock_sm.mysumstats.data.columns = []
                    mock_sm_class.return_value = mock_sm
                    result = CliRunner().invoke(
                        main,
                        [
                            "-c",
                            "tests/data/test_config.yaml",
                            "-i",
                            "tests/data/test_sumstats.tsv",
                            "-f",
                            "plink_pvar",
                            "-o",
                            temporary_directory,
                        ],
                    )

                    self.assertEqual(result.exit_code, 0, result.output)
                    mock_sm.order_alleles.assert_called_once_with(**sort_params)

    @patch("gwaspipe.gwaspipe.ConfigurationManager")
    @patch("gwaspipe.gwaspipe.SumstatsManager")
    def test_main_accepts_both_conflicting_snpid_step_names(self, mock_sm_class, mock_cm_class):
        """The preferred step name and its deprecated alias have identical behaviour."""
        from click.testing import CliRunner

        from gwaspipe.gwaspipe import main

        mock_cm = MagicMock()
        mock_cm.log_file_path = Path("test.log")
        mock_cm.formatbook_path = Path("data/formatbook.json")
        mock_cm.filename_settings = (None, None)
        mock_cm.step.return_value = ({"run": True, "workspace": "default"}, {})
        mock_cm_class.return_value = mock_cm

        with TemporaryDirectory() as temporary_directory:
            mock_cm.root_path = temporary_directory
            for step_name in ("filter_conflicting_snpids", "check_ambiguous_snps"):
                with self.subTest(step_name=step_name):
                    mock_cm.run_sequence = (step_name,)
                    mock_sm = MagicMock()
                    mock_sm.mysumstats.data = pd.DataFrame(
                        {
                            "SNPID": ["1:100:A:G", "1:100:A:G"],
                            "EAF": [0.2, 0.2],
                            "BETA": [0.1, 0.1],
                            "SE": [0.01, 0.01],
                            "CHR": [1, 1],
                            "POS": [100, 100],
                        }
                    )
                    mock_sm_class.return_value = mock_sm

                    result = CliRunner().invoke(
                        main,
                        [
                            "-c",
                            "tests/data/test_config.yaml",
                            "-i",
                            "tests/data/test_sumstats.tsv",
                            "-f",
                            "plink_pvar",
                            "-o",
                            temporary_directory,
                        ],
                    )

                    self.assertEqual(result.exit_code, 0, result.output)
                    self.assertEqual(len(mock_sm.mysumstats.data), 1)


class TestAssemblyValidation(unittest.TestCase):
    def setUp(self):
        self.sumstats = MagicMock()
        self.sumstats.data = pd.DataFrame({"CHR": [1], "POS": [100]})
        self.sumstats.meta = {"gwaslab": {}}
        self.config = {
            "genome_assembly": "GRCh38",
            "reference_resources": {"reference_fasta": "GRCh38, release 109"},
            "assembly_validation": {"min_hapmap3_matches": 10000, "allow_override": False},
        }

    @patch("gwaspipe.gwaspipe._hapmap3_match_counts", return_value={"19": 2, "38": 10000})
    def test_records_successful_assembly_audit(self, mock_match_counts):
        audit = validate_declared_assembly(self.sumstats, self.config)

        self.assertEqual(audit["decision"], "passed")
        self.assertEqual(audit["inferred_assembly"], "GRCh38")
        self.assertEqual(audit["hapmap3_match_counts"], {"GRCh37": 2, "GRCh38": 10000})
        self.assertEqual(self.sumstats.meta["gwaspipe"]["assembly_validation"], audit)
        self.sumstats.infer_build.assert_called_once_with()
        mock_match_counts.assert_called_once_with(self.sumstats.data)

    @patch("gwaspipe.gwaspipe._hapmap3_match_counts", return_value={"19": 10000, "38": 2})
    def test_rejects_discordant_assembly_without_override(self, _):
        with self.assertRaisesRegex(AssemblyValidationError, "disagree"):
            validate_declared_assembly(self.sumstats, self.config)

        self.assertEqual(self.sumstats.meta["gwaspipe"]["assembly_validation"]["decision"], "failed")

    @patch("gwaspipe.gwaspipe._hapmap3_match_counts", return_value={"19": 2, "38": 10000})
    def test_rejects_missing_assembly_declaration_without_override(self, _):
        self.config.pop("genome_assembly")

        with self.assertRaisesRegex(AssemblyValidationError, "missing or unsupported"):
            validate_declared_assembly(self.sumstats, self.config)

    @patch("gwaspipe.gwaspipe._hapmap3_match_counts", return_value={"19": 10000, "38": 10000})
    def test_rejects_unknown_inference_without_override(self, _):
        with self.assertRaisesRegex(AssemblyValidationError, "ambiguous inferred"):
            validate_declared_assembly(self.sumstats, self.config)

    @patch("gwaspipe.gwaspipe._hapmap3_match_counts", return_value={"19": 0, "38": 4})
    def test_rejects_low_evidence_without_override(self, _):
        with self.assertRaisesRegex(AssemblyValidationError, "fewer than 10000"):
            validate_declared_assembly(self.sumstats, self.config)

    @patch("gwaspipe.gwaspipe._hapmap3_match_counts", return_value={"19": 0, "38": 4})
    def test_records_explicit_override_for_low_evidence(self, _):
        self.config["assembly_validation"] = {
            "min_hapmap3_matches": 10000,
            "allow_override": True,
            "override_reason": "Validated against the source study manifest.",
        }

        audit = validate_declared_assembly(self.sumstats, self.config)

        self.assertEqual(audit["decision"], "overridden")
        self.assertEqual(audit["override_reason"], "Validated against the source study manifest.")

    def test_requires_a_completed_validation(self):
        with self.assertRaisesRegex(AssemblyValidationError, "infer_build validation"):
            _require_validated_assembly(self.sumstats)

    def test_writes_json_provenance_sidecar(self):
        self.sumstats.meta["gwaspipe"] = {"assembly_validation": {"decision": "passed"}}
        with TemporaryDirectory() as temporary_directory:
            output_path = Path(temporary_directory, "summary_statistics")
            source_path = Path("input.tsv")
            _write_run_provenance(output_path, self.sumstats, source_path)
            provenance_path = Path(f"{output_path}.provenance.json")

            self.assertTrue(provenance_path.exists())
            provenance = json.loads(provenance_path.read_text())
            self.assertEqual(provenance["gwaspipe_version"], __version__)
            timestamp = datetime.fromisoformat(provenance["timestamp_utc"])
            self.assertEqual(timestamp.tzinfo, UTC)
            self.assertEqual(provenance["source_path"], str(source_path))
            self.assertEqual(provenance["output_path"], str(output_path))
            self.assertEqual(provenance["assembly_validation"]["decision"], "passed")


if __name__ == "__main__":
    unittest.main()
