#!/usr/bin/env bash

# exit when any command fails
set -o nounset
set -o errexit
set -x
set -euox pipefail
set -e

# force clean it up
function cleanup() {
    echo "REMOVING OLD FILES IF EXISTS";
    rm -rf build/temp*
    rm -rf build/lib.linux*
    rm -rf dist/*
    rm -rf __pycache__/
    rm -rf *pyc
    rm -rf *so
    rm -rf pydbretina/DBRetina.egg-info/
    rm -rf build/bdist.linux-x86_64
}

trap cleanup EXIT
cleanup

$(which python) -m pip uninstall DBRetina -y || true
$(which python) -m pip install .

echo "BUILD COMPLETE"
