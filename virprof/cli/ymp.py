import click
import os

@click.command(context_settings=dict(
    ignore_unknown_options=True,
    allow_extra_args=True
))
@click.pass_context
def cli(ctx) -> None:
    """Run ymp command (for container passthrough)"""
    os.execvp("ymp", ["ymp"] + ctx.args)
    
