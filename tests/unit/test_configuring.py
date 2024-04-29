import unittest
from pathlib import Path

from src.configuring import ConfigurationManager


class TestConfigurationManager(unittest.TestCase):
    def setUp(self) -> None:
        self.rootpath = "tests/results"
        self.rootpath_name = self.rootpath.split("/")[1]

        self.test_configfile = "tests/data/test_config.yaml"
        self.log_filename = "file.log"

        self.formatbook_path = "data/formatbook.json"

        self.params = ["run", "workspace"]

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
        config = ConfigurationManager(root_path=self.rootpath_name)
        assert isinstance(config, ConfigurationManager)
        assert isinstance(config.root_path, Path)
        assert config.root_path.name == self.rootpath_name
        assert str(config.root_path) == self.rootpath

        config = ConfigurationManager(root_path=Path(self.rootpath_name))
        assert isinstance(config, ConfigurationManager)
        assert isinstance(config.root_path, Path)
        assert config.root_path.name == self.rootpath_name
        assert str(config.root_path) == self.rootpath

    def test_log_file_path(self):
        config = ConfigurationManager(config_file=Path(self.test_configfile))
        assert config.log_file_path == Path(self.rootpath, self.log_filename)

    def test_formatbook_path(self):
        config = ConfigurationManager(config_file=Path(self.test_configfile))
        assert config.formatbook_path == self.formatbook_path

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
        assert tuple(mask) == tuple([True, False])
        assert sep == "."
