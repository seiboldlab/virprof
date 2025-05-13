"""Virprof pipeline custom functions for pathogen recovery"""

from pathlib import Path

def get_resource(path):
    roots = [Path(__file__).parent.parent, Path(__file__).parent / "pipeline"]
    for root in roots:
        loc = root / path
        if loc.exists():
            return loc
    raise Exception(f"Unable to find {path}. Please check your installation")

def get_ymp_yml():
    return str(get_resource("virprof.yml"))
