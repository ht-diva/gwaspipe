import time

from cloup import Context, HelpFormatter, HelpTheme, Style
from loguru import logger as a_logger

__all__ = ["__appname__", "context_settings", "logger"]

__appname__ = "gwaspipe"

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


# GWASLab compatible logger
class Log:
    def __init__(self, path):
        self.path = path

    def write(self, *message, end="\n", show_time=True, verbose=True):
        timestamp = str(time.strftime('%Y/%m/%d %H:%M:%S')) if show_time else ""
        line = (timestamp + " " if show_time else "") + " ".join(map(str, message)) + end

        if verbose:
            print(line, end="")

        with open(self.path, "a") as f:
            f.write(line)
