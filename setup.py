#!/usr/bin/env python
from setuptools import setup

import os
import shutil
import itertools
from glob import glob

extra_package_data = {
    "virprof/pipeline": [
        "R/*.R",
        "bin/*",
        "conda_envs/*/*.txt",
        "data/*",
        "rules/*",
        "test_data/*",
        "*.yml",
        "README.md",
    ]
}
extra_exclude_package_data = {"virprof.pipeline": ["environment.yml"]}
package_data = {"virprof.pipeline": ["*", "**/*"]}
for submodule in extra_package_data:
    if os.path.exists(submodule):
        raise RuntimeError(
            "Cannot add data submodule '%s': submodule already exists" % submodule
        )
try:
    for submodule, patterns in extra_package_data.items():
        expanded = map(glob, patterns)
        matches = itertools.chain.from_iterable(expanded)
        files = filter(os.path.isfile, matches)
        for src in files:
            dst = os.path.join(submodule, src)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            print("%s -> %s" % (src, dst))
            shutil.copy(src, dst)
        with open(os.path.join(submodule, "__init__.py"), "w"):
            pass
    setup(package_data=package_data)
finally:
    for submodule in extra_package_data:
        shutil.rmtree(submodule)
