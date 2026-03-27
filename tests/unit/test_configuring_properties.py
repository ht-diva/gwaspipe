import unittest
from pathlib import Path

from gwaspipe.configuring import ConfigurationManager, SingletonConfigurationManager, load_config


class TestLoadConfig(unittest.TestCase):
    """Tests for the load_config function."""

    def test_load_valid_config(self):
        """Test loading a valid configuration file."""
        config_path = Path("tests/data/test_config.yaml")
        config = load_config(config_path)
        self.assertIsInstance(config, dict)
        self.assertIn("root_path", config)
        self.assertIn("log_filename", config)

    def test_load_nonexistent_config(self):
        """Test loading a non-existent configuration file."""
        config_path = Path("nonexistent.yaml")
        with self.assertRaises(FileNotFoundError):
            load_config(config_path)

    def test_load_none_config(self):
        """Test loading with None as a config file."""
        with self.assertRaises(FileNotFoundError):
            load_config(None)


class TestConfigurationManagerProperties(unittest.TestCase):
    """Tests for ConfigurationManager properties."""

    def setUp(self):
        # Clear the singleton instances before each test
        SingletonConfigurationManager._instances = {}
        self.config_file = Path("tests/data/test_config.yaml")

    def test_filename_settings_property(self):
        """Test filename_settings property."""
        cm = ConfigurationManager(config_file=self.config_file)
        mask, sep = cm.filename_settings
        self.assertIsNone(mask)
        self.assertIsNone(sep)

    def test_formatbook_path_property_with_custom_path(self):
        """Test formatbook_path property with a custom path."""
        # Create a temporary formatbook file
        import tempfile

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            custom_path = Path(f.name)
            custom_path.write_text('{"test": "data"}')

        try:
            cm = ConfigurationManager(config_file=self.config_file, formatbook_file=custom_path)
            self.assertEqual(cm.formatbook_path, custom_path)
        finally:
            custom_path.unlink()

    def test_formatbook_path_property_with_default_path(self):
        """Test formatbook_path property with default path."""
        cm = ConfigurationManager(config_file=self.config_file)
        default_path = Path(__file__).parent.parent.parent / "src" / "gwaspipe" / "data" / "formatbook.json"
        self.assertEqual(cm.formatbook_path, default_path)

    def test_log_file_path_property(self):
        """Test log_file_path property."""
        cm = ConfigurationManager(config_file=self.config_file)
        log_path = cm.log_file_path
        self.assertIsInstance(log_path, Path)
        self.assertTrue(str(log_path).endswith("file.log"))

    def test_root_path_property_with_custom_path(self):
        """Test root_path property with custom path."""
        custom_path = Path("custom_root")
        cm = ConfigurationManager(config_file=self.config_file, root_path=custom_path)
        self.assertEqual(cm.root_path, custom_path)

    def test_root_path_property_with_config_path(self):
        """Test root_path property with path from config."""
        cm = ConfigurationManager(config_file=self.config_file)
        self.assertIsInstance(cm.root_path, Path)

    def test_run_sequence_property(self):
        """Test run_sequence property."""
        cm = ConfigurationManager(config_file=self.config_file)
        run_sequence = cm.run_sequence
        self.assertIsInstance(run_sequence, tuple)

    def test_steps_property(self):
        """Test steps property."""
        cm = ConfigurationManager(config_file=self.config_file)
        steps = cm.steps
        self.assertIsInstance(steps, dict)

    def test_step_method_with_existing_step(self):
        """Test step method with existing step."""
        cm = ConfigurationManager(config_file=self.config_file)
        params, gl_params = cm.step("basic_check")
        self.assertIsInstance(params, dict)
        self.assertIn("run", params)
        self.assertIsInstance(gl_params, dict)

    def test_step_method_with_nonexistent_step(self):
        """Test step method with non-existent step."""
        cm = ConfigurationManager(config_file=self.config_file)
        params, gl_params = cm.step("nonexistent_step")
        self.assertIsInstance(params, dict)
        self.assertEqual(params, {"run": False, "workspace": ""})
        self.assertEqual(gl_params, {})


class TestSingletonPattern(unittest.TestCase):
    """Tests for the singleton pattern implementation."""

    def setUp(self):
        # Clear the singleton instances before each test
        SingletonConfigurationManager._instances = {}

    def test_singleton_creates_single_instance(self):
        """Test that only one instance is created."""
        cm1 = ConfigurationManager(config_file="tests/data/test_config.yaml")
        cm2 = ConfigurationManager(config_file="tests/data/test_config.yaml")
        self.assertIs(cm1, cm2)

    def test_singleton_with_different_configs(self):
        """Test singleton behavior with different configs."""
        cm1 = ConfigurationManager(config_file="tests/data/test_config.yaml")
        cm2 = ConfigurationManager(config_file=None)
        self.assertIs(cm1, cm2)


if __name__ == "__main__":
    unittest.main()
