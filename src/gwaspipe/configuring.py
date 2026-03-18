from pathlib import Path

from ruamel.yaml import YAML


def load_config(config_file: Path | str | None) -> dict:
    """Load and validate configuration from YAML file."""
    if not config_file:
        raise FileNotFoundError("No configuration file provided")

    config_path = Path(config_file)
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file {config_path.name} not found")

    with config_path.open("r") as file:
        return YAML(typ="safe").load(file)


class SingletonConfigurationManager(type):
    """Metaclass."""

    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(SingletonConfigurationManager, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class ConfigurationManager(metaclass=SingletonConfigurationManager):
    def __init__(
        self,
        config_file: Path | str | None = None,
        formatbook_file: Path | str | None = None,
        root_path: Path | str | None = None,
    ):
        self.config = load_config(config_file)
        self._root_path = Path(root_path) if root_path else None
        self._formatbook_path = Path(formatbook_file) if formatbook_file else None

    @property
    def filename_settings(self):
        mask = self.config.get("filename_mask")
        sep = self.config.get("filename_sep")
        return mask, sep

    @property
    def formatbook_path(self) -> Path:
        """Return the path to the formatbook file, using default if custom path doesn't exist."""
        if self._formatbook_path and self._formatbook_path.exists():
            return self._formatbook_path
        return Path(__file__).parent / "data" / "formatbook.json"

    @property
    def log_file_path(self) -> Path:
        return Path(self.root_path, self.config["log_filename"])

    @property
    def root_path(self) -> Path:
        root_path = self._root_path or self.config.get("root_path")
        if not root_path:
            raise ValueError("Root path not specified in config or constructor")

        path = Path(root_path)
        path.parent.mkdir(exist_ok=True)
        return path

    @property
    def run_sequence(self):
        run_sequence = tuple(self.config.get("run_sequence", {}).values())
        return run_sequence

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
