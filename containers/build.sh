#!/bin/bash
VERS=${1-master}
YMP_VERS=master
echo "Building Virprof $VERS"
apptainer build \
  --force \
  --build-arg version=$VERS \
  --build-arg ymp_version=$YMP_VERS \
  virprof-$VERS.sif virprof.def \
  |& tee virprof-$VERS.build.log
