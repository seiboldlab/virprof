import click
import pathlib
from .. import get_resource
from subprocess import run, DEVNULL
from shutil import copy2
import logging

LOG = logging.getLogger(__name__)


@click.command()
@click.argument("path", type=click.Path(path_type=pathlib.Path))
@click.option("--force", is_flag=True)
def cli(path, force) -> bool:
    """Set up a demo/test directory"""
    testdata = get_resource("test_data")
    ymp_yml = get_resource("example_ymp.yml")
    config_yml = get_resource("config.yml")
    whitelist_txt = get_resource("respiratory_virus_whitelist.txt")
    logout = DEVNULL if LOG.getEffectiveLevel() >= logging.INFO else None

    if not path.exists():
        path.mkdir(parents=True)
    elif not path.is_dir():
        raise click.ClickException("Not a directory: {path}")
    elif not force and any(path.iterdir()):
        raise click.ClickException("Directory not empty: {path}")

    run(testdata / "install_testdata.sh", cwd=path, check=True, stdout=logout)

    click.echo("Installing ymp.yml")
    copy2(ymp_yml, path / "ymp.yml")
    click.echo("Installing config.yml")
    copy2(config_yml, path)
    click.echo("Installing 'respiratory_virus_whitelist.txt'")
    copy2(whitelist_txt, path)

    return True
