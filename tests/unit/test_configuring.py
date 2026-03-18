import unittest
from pathlib import Path

from gwaspipe.configuring import ConfigurationManager, SingletonConfigurationManager


class TestConfigurationManagerMissing(unittest.TestCase):
    def setUp(self):
        # Clear the singleton instances before each test
        SingletonConfigurationManager._instances = {}

    def test_missing_config_file(self):
        # Test case where configuration file is missing
        config_file = Path("nonexistent.config")
        with self.assertRaises(FileNotFoundError):
            ConfigurationManager(config_file=config_file)

        # Verify the error message
        msg = str(config_file)
        self.assertEqual(msg, "nonexistent.config")

    def test_no_config_file_provided(self):
        # Test case where no configuration file is provided
        config_file = None
        with self.assertRaises(FileNotFoundError):
            ConfigurationManager(config_file=config_file)

        # Verify the error message
        msg = "No configuration file provided"
        self.assertEqual(msg, msg)


class TestConfigurationManager(unittest.TestCase):
    def setUp(self) -> None:
        # Clear the singleton instances before each test
        SingletonConfigurationManager._instances = {}

        self.rootpath = "tests/results"
        self.rootpath_name = self.rootpath.split("/")[1]

        self.test_configfile = "tests/data/test_config.yaml"
        self.log_filename = "file.log"

        self.formatbook_path = "data/formatbook.json"

        self.params = ["run", "workspace"]

    def test_singleton_pattern(self):
        cm1 = ConfigurationManager(config_file=self.test_configfile)
        cm2 = ConfigurationManager(config_file=None)
        self.assertEqual(cm1, cm2)

    def test_root_path_from_configfile(self):
        config = ConfigurationManager(config_file=self.test_configfile)
        assert isinstance(config, ConfigurationManager)
        assert isinstance(config.root_path, Path)
        assert config.root_path.name == self.rootpath_name
        assert str(config.root_path) == self.rootpath

        config = ConfigurationManager(config_file=Path(self.test_configfile))
        assert isinstance(config, ConfigurationManager)
        assert isinstance(config.root_path, Path)
        assert config.root_path.name == self.rootpath_name
        assert str(config.root_path) == self.rootpath

    def test_root_path_from_cli(self):
        config = ConfigurationManager(config_file=Path(self.test_configfile), root_path=self.rootpath)
        assert isinstance(config, ConfigurationManager)
        assert isinstance(config.root_path, Path)
        assert config.root_path.name == self.rootpath_name
        assert str(config.root_path) == self.rootpath

        config = ConfigurationManager(config_file=Path(self.test_configfile), root_path=Path(self.rootpath))
        assert isinstance(config, ConfigurationManager)
        assert isinstance(config.root_path, Path)
        assert config.root_path.name == self.rootpath_name
        assert str(config.root_path) == self.rootpath

    def test_log_file_path(self):
        config = ConfigurationManager(config_file=Path(self.test_configfile))
        assert config.log_file_path == Path(self.rootpath, self.log_filename)

    def test_formatbook_path(self):
        config = ConfigurationManager(config_file=Path(self.test_configfile))
        assert str(self.formatbook_path) in str(config.formatbook_path)

    def test_step_default(self):
        config = ConfigurationManager(config_file=Path(self.test_configfile))
        params, gl_params = config.step("name")
        assert isinstance(params, dict) and len(params) == 2
        self.assertIn(self.params[0], params.keys())
        self.assertFalse(params[self.params[0]])
        self.assertIn(self.params[1], params.keys())
        assert isinstance(params[self.params[1]], str)
        assert len(params[self.params[1]]) == 0
        assert isinstance(gl_params, dict) and len(gl_params) == 0

    def test_run_sequence_default(self):
        config = ConfigurationManager(config_file=Path(self.test_configfile))
        assert isinstance(config.run_sequence, tuple)
        assert len(config.run_sequence) == 0

    def test_filename_settings_default(self):
        config = ConfigurationManager(config_file=Path(self.test_configfile))
        mask, sep = config.filename_settings
        assert mask is None
        assert sep is None
