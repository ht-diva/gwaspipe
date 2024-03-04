from pathlib import Path

from ruamel.yaml import YAML


class SingletonConfigurationManager(type):
    """Metaclass."""

    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(SingletonConfigurationManager, cls).__call__(
                *args, **kwargs
            )
        return cls._instances[cls]


class ConfigurationManager(metaclass=SingletonConfigurationManager):
    def __init__(self, config_file=None, root_path=None):
        self.c = None
        if not isinstance(config_file, Path):
            config_file = Path(config_file)
        if config_file.exists():
            yaml = YAML(typ="safe")
            with open(config_file, "r") as file:
                self.c = yaml.load(file)
        else:
            msg = f"configuration file {config_file.name} not found"
            exit(msg)
        self._root_path = root_path

    @property
    def root_path(self):
        root_path = self._root_path if self._root_path else self.c["root_path"]
        path = Path(root_path)
        path.mkdir(exist_ok=True)
        return path

    @property
    def log_file_path(self):
        return Path(self.root_path, self.c["log_filename"])

    @property
    def report_if_filename(self):
        return self.c["report_if_filename"]

    @property
    def formatbook_path(self):
        return self.c["formatbook_path"]

    @property
    def steps(self):
        return self.c["steps"].keys()

    def step(self, name):
        # use get method with default values
        params = self.c["steps"].get(name, {}).get("params", {"run": False,
                                                              "workspace": ""})
        for k,v in [('run', False), ('workspace', '')]:
            if k not in params.keys():
                params[k] = v

        gl_params = self.c["steps"].get(name, {}).get("gl_params", {})

        return params, gl_params

    @property
    def run_sequence(self):
        run_sequence = tuple(self.c["run_sequence"].values())
        return run_sequence

    @property
    def filename_settings(self):
        mask = self.c.get('filename_mask', False)
        sep = self.c.get('filename_sep', '.')
        return mask, sep
