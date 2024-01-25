from appdirs import user_log_dir
from cloup import Context, HelpFormatter, HelpTheme, Style
from loguru import logger as a_logger

__all__ = ["__appname__", "context_settings", "log_file", "logger"]

__appname__ = "gwaspipe"
log_file = user_log_dir(__appname__)

# Check the docs for all available arguments of HelpFormatter and HelpTheme.
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
