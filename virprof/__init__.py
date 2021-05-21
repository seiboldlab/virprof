"""Contains the custom code for blast based binning"""

import sys
import os

if sys.version_info.major < 3 or sys.version_info.minor < 6:
    print("ERROR: Python version must be at least 3.6!")
    print("Found ", sys.version)
    sys.exit(1)


def get_ymp_yml():
    try:
        import virprof.pipeline

        paths = virprof.pipeline.__path__
    except ModuleNotFoundError:
        import virprof

        paths = [os.path.dirname(modpath) for modpath in virprof.__path__]
    for path in paths:
        ymp_yml = os.path.join(path, "virprof.yml")
        if os.path.exists(ymp_yml):
            return ymp_yml
    return None
