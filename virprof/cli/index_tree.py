import click

from ..taxonomy import load_taxonomy


@click.command()
@click.option("--out", "-o", type=click.File("w"), required=True, help="Output binary")
@click.option(
    "--library",
    "-b",
    type=click.Choice(["graph_tool", "networkx"]),
    default="graph_tool",
    help="Tree library to use",
)
@click.option(
    "--ncbi-taxonomy",
    type=click.Path(),
    required=True,
    help="Path to the NCBI taxonomy dump directory",
)
def cli(out: click.utils.LazyFile, library: str, ncbi_taxonomy: str) -> bool:
    """Parse NCBI taxonomy from dump files and write tree to binary"""
    taxonomy = load_taxonomy(ncbi_taxonomy, library=library)
    taxonomy.save_tree_binary(out.name)
    return True
