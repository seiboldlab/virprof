import click
import pathlib
from .. import get_resource
from subprocess import run

@click.command()
@click.argument("path", type=click.Path(path_type=pathlib.Path))
def cli(path) -> bool:
    """Set up a demo/test directory"""
    testdata = get_resource("test_data")

    if not path.exists():
        path.mkdir(parents=True)
    elif not path.is_dir():
        raise click.ClickException("Not a directory: {path}")
    elif any(path.iterdir()):
        raise click.ClickException("Directory not empty: {path}")

    run(testdata / "install_testdata.sh", cwd = path, check = True)
    
    return True
