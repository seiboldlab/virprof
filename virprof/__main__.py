"""Allow calling CLI via python -m virprof"""

from virprof.cli import main

if __name__ == "__main__":
    # pylint: disable=unexpected-keyword-arg, no-value-for-parameter
    main(prog_name="python -m virprof")
