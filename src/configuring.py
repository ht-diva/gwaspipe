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
    def __init__(self, config_file=None):
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

    @property
    def root_path(self):
        path = Path(self.c["root_path"])
        path.mkdir(exist_ok=True)
        return path

    @property
    def log_file_path(self):
        return Path(self.root_path, self.c["log_filename"])

    @property
    def report_filename(self):
        return self.c["report_filename"]

    @property
    def steps(self):
        return self.c["steps"].keys()

    def step(self, name):
        params = {'run': False,
                  'workspace': ''}
        gl_params = {}
        if name in self.steps:
            for k, v in self.c["steps"][name].items():
                if k in params.keys():
                    params[k] = v
                else:
                    gl_params[k] = v
        return params, gl_params

    @property
    def run_sequence(self):
        run_sequence = tuple(self.c["run_sequence"].values())
        return run_sequence
