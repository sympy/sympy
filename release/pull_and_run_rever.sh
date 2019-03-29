#!/bin/bash

# This is a separate script from rever so that it can pull any updates to the
# rever file. The first argument should be the branch name, and the second
# argument should be the version (can be omitted if it is the same as the
# branch name).

set -e
set -x

if [[ -z $2 ]]; then
    $2=$1
fi
git pull
git checkout $1
git pull

shift
/opt/conda/bin/rever "$@"
