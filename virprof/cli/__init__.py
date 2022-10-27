"""
CLI module

The commands are each in files named as the command is named.
"""

import os
import logging
import csv

from importlib import import_module

import click
import tqdm  # type: ignore


#: Increase CSV field size limit to 2GB
csv.field_size_limit(2**31)


class TqdmHandler(logging.Handler):
    """Logging handler passing writes through tqdm

    This is used so progress bar and log messages co-habitate stdout without
    clobbering lines.

    """

    def __init__(self) -> None:
        logging.Handler.__init__(self)

    def emit(self, record: logging.LogRecord) -> None:
        try:
            msg = self.format(record)
            tqdm.tqdm.write(msg)
        except (KeyboardInterrupt, SystemExit, RecursionError):
            raise
        except Exception:  # pylint: disable=broad-except
            self.handleError(record)


def setup_logging() -> None:
    """Sets up python logging facility"""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S",
        handlers=[TqdmHandler()],
    )


def setup_profiling(ctx, _attr, value) -> None:
    """Start yappi profiler

    Profiling ends and results are printed when the program ends.

    Requires "yappi" to have been installed.
    """
    if not value or ctx.resilient_parsing:
        return value
    # pylint: disable=import-outside-toplevel,import-error
    import yappi  # type: ignore
    import atexit

    def dump_profile() -> None:
        """Print the profile at exit"""
        stats = yappi.get_func_stats()
        stats.sort("ttot")
        stats.print_all(
            columns={
                0: ("name", 120),
                1: ("ncall", 10),
                2: ("tsub", 8),
                3: ("ttot", 8),
                4: ("tavg", 8),
            }
        )
        yappi.stop()

    atexit.register(dump_profile)
    yappi.start()
    return value


def setup_debug(ctx, _attr, value) -> None:
    """Installs handler to print traceback on receiving SIGUSR1

    This helps figure out where the code is at without interrupting
    long running calculations.
    """

    if not value or ctx.resilient_parsing:
        return value
    # pylint: disable=import-outside-toplevel
    import signal
    import traceback

    def handlesigusr1(sig, frame):
        print("------------------")
        print(f"Dumping stack (caught signal {sig})")
        traceback.print_stack(f=frame, limit=5)
        print("------------------")

    signal.signal(signal.SIGUSR1, handlesigusr1)
    return value


class AutoLoadCommand(click.MultiCommand):
    """Click command type to auto-generate commands from modules

    Scans for modules in the same directory as this one and creates
    commands for each module defining a function cli().
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.commands_folder = os.path.dirname(__file__)

    def list_commands(self, ctx):
        return sorted(
            [
                filename[:-3].replace("_", "-")
                for filename in os.listdir(self.commands_folder)
                if filename.endswith(".py")
                and not filename == "__init__.py"
                and not filename.startswith(".")
                and not filename.startswith("_")
            ]
        )

    def get_command(self, ctx, cmd_name):
        mod = import_module("." + cmd_name.replace("-", "_"), __name__)
        try:
            return mod.cli
        except AttributeError:  # no .cli
            return None


@click.command(cls=AutoLoadCommand)
@click.option(
    "--profile",
    is_flag=True,
    expose_value=False,
    callback=setup_profiling,
    help="Dump profiling info; requires yappi installed",
)
@click.option(
    "--traceback-on-usr1",
    is_flag=True,
    expose_value=False,
    callback=setup_debug,
    help="Dump stack trace on receiving signal SIGUSR1",
)
def main() -> None:
    """Use any of the subcommands"""
    setup_logging()
