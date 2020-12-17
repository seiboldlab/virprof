"""Contains the custom code for blast based binning"""

import sys

if sys.version_info.major < 3 or sys.version_info.minor < 6:
    print("ERROR: Python version must be at least 3.6!")
    print("Found ", sys.version)
    sys.exit(1)
