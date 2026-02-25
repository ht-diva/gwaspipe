import importlib.metadata


from cloup import Context, HelpFormatter, HelpTheme, Style
from gwaslab.info.g_Log import Log
from loguru import logger as a_logger


__all__ = ["__appname__", "__version__", "formatter_settings", "context_settings", "logger", "Log"]

__appname__ = __name__

try:
    __version__ = importlib.metadata.version(__appname__)
except importlib.metadata.PackageNotFoundError:
    __version__ = "unknown"

formatter_settings = HelpFormatter.settings(
    theme=HelpTheme(
        invoked_command=Style(fg="red"),
        heading=Style(fg="bright_white", bold=True),
        constraint=Style(fg="magenta"),
        col1=Style(fg="yellow"),
    )
)
context_settings = Context.settings(formatter_settings=formatter_settings)

logger = a_logger
logger.remove()
