from pathlib import Path

from ruamel.yaml import YAML


class SingletonConfigurationManager(type):
    """Metaclass."""

    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(SingletonConfigurationManager, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class ConfigurationManager(metaclass=SingletonConfigurationManager):
    def __init__(self, config_file: Path | str = None, root_path: Path | str = None):
        self.config = None
        config_path = Path(config_file) if config_file else None
        if config_path:
            if config_path.exists():
                yaml = YAML(typ="safe")
                with open(config_file, "r") as file:
                    self.config = yaml.load(file)
            else:
                msg = (
                    f"Configuration file {config_path.name} not found"
                    if config_path
                    else "No configuration file provided"
                )
                raise FileNotFoundError(msg)
        self._root_path = Path(root_path) if root_path else None

    @property
    def root_path(self):
        root_path = self._root_path or self.config["root_path"]
        path = Path(root_path)
        path.parent.mkdir(exist_ok=True)
        return path

    @property
    def log_file_path(self):
        return Path(self.root_path, self.config["log_filename"])

    @property
    def formatbook_path(self):
        return self.config["formatbook_path"]

    @property
    def steps(self):
        return self.config.get("steps", {})

    def step(self, step_name):
        # Get the step configuration with default parameters
        step_config = self.steps.get(step_name, {})
        default_params = {"run": False, "workspace": ""}
        params = {**default_params, **step_config.get("params", {})}

        # Get the GL parameters if available
        gl_params = step_config.get("gl_params", {})

        return params, gl_params

    @property
    def run_sequence(self):
        run_sequence = tuple(self.config.get("run_sequence", {}).values())
        return run_sequence

    @property
    def filename_settings(self):
        mask = self.config.get("filename_mask", [True, False])
        sep = self.config.get("filename_sep", ".")
        return mask, sep
