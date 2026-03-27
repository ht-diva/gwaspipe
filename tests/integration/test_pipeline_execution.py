import shutil
import tempfile
import unittest
from pathlib import Path

import gwaslab as gl
import pandas as pd

from gwaspipe.configuring import ConfigurationManager
from gwaspipe.gwaspipe import SumstatsManager


class TestPipelineExecution(unittest.TestCase):
    """Integration tests for pipeline execution."""

    def setUp(self):
        # Create a temporary directory for test outputs
        self.test_dir = tempfile.mkdtemp()
        self.config_file = Path("tests/data/test_config.yaml")
        self.test_data_path = Path("tests/data/test_sumstats.tsv")

    def tearDown(self):
        # Clean up temporary directory
        if Path(self.test_dir).exists():
            shutil.rmtree(self.test_dir)

    def test_basic_pipeline_execution(self):
        """Test basic pipeline execution with minimal configuration."""
        # Create a minimal config for testing
        test_config = {
            "root_path": self.test_dir,
            "log_filename": "test.log",
            "filename_mask": None,
            "filename_sep": None,
            "run_sequence": {"basic_check": "basic_check"},
            "steps": {"basic_check": {"params": {"run": True, "workspace": "default"}, "gl_params": {}}},
        }

        # Write test config
        import yaml

        config_path = Path(self.test_dir, "test_config.yaml")
        with open(config_path, "w") as f:
            yaml.dump(test_config, f)

        # Test pipeline execution
        cm = ConfigurationManager(config_file=config_path, root_path=self.test_dir)
        sm = SumstatsManager(
            input_path=str(self.test_data_path),
            input_format="plink_pvar",
            input_separator="\t",
            input_study=None,
            formatbook_path=cm.formatbook_path,
            pid=False,
            bcfliftover=False,
        )

        # Execute basic_check step
        params, gl_params = cm.step("basic_check")
        if params.get("run", False):
            sm.mysumstats.basic_check(**gl_params)

        # Verify results
        self.assertIsNotNone(sm.mysumstats)
        self.assertTrue(sm.mysumstats.data is not None)

    def test_multiple_step_pipeline(self):
        """Test pipeline with multiple steps."""
        test_config = {
            "root_path": self.test_dir,
            "log_filename": "test.log",
            "filename_mask": None,
            "filename_sep": None,
            "run_sequence": {"basic_check": "basic_check", "sort_alphabetically": "sort_alphabetically"},
            "steps": {
                "basic_check": {"params": {"run": True, "workspace": "default"}, "gl_params": {}},
                "sort_alphabetically": {"params": {"run": True, "workspace": "default"}, "gl_params": {"n_cores": 1}},
            },
        }

        import yaml

        config_path = Path(self.test_dir, "test_config.yaml")
        with open(config_path, "w") as f:
            yaml.dump(test_config, f)

        cm = ConfigurationManager(config_file=config_path, root_path=self.test_dir)
        sm = SumstatsManager(
            input_path=str(self.test_data_path),
            input_format="plink_pvar",
            input_separator="\t",
            input_study=None,
            formatbook_path=cm.formatbook_path,
            pid=False,
            bcfliftover=False,
        )

        # Execute multiple steps
        for step in cm.run_sequence:
            params, gl_params = cm.step(step)
            if params.get("run", False):
                if step == "basic_check":
                    sm.mysumstats.basic_check(**gl_params)
                elif step == "sort_alphabetically":
                    sm.order_alleles(n_cores=gl_params.get("n_cores", 1))

        # Verify results
        self.assertIsNotNone(sm.mysumstats)
        self.assertTrue(sm.mysumstats.data is not None)

    def test_pipeline_with_different_formats(self):
        """Test the pipeline with different input formats."""
        formats_to_test = ["plink_pvar"]  # Add more formats as needed

        for fmt in formats_to_test:
            with self.subTest(format=fmt):
                test_config = {
                    "root_path": self.test_dir,
                    "log_filename": "test.log",
                    "filename_mask": None,
                    "filename_sep": None,
                    "run_sequence": {"basic_check": "basic_check"},
                    "steps": {"basic_check": {"params": {"run": True, "workspace": "default"}, "gl_params": {}}},
                }

                import yaml

                config_path = Path(self.test_dir, f"test_config_{fmt}.yaml")
                with open(config_path, "w") as f:
                    yaml.dump(test_config, f)

                cm = ConfigurationManager(config_file=config_path, root_path=self.test_dir)
                sm = SumstatsManager(
                    input_path=str(self.test_data_path),
                    input_format=fmt,
                    input_separator="\t",
                    input_study=None,
                    formatbook_path=cm.formatbook_path,
                    pid=False,
                    bcfliftover=False,
                )

                params, gl_params = cm.step("basic_check")
                if params.get("run", False):
                    sm.mysumstats.basic_check(**gl_params)

                self.assertIsNotNone(sm.mysumstats)


if __name__ == "__main__":
    unittest.main()
