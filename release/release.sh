#!/bin/bash

# This is the main entrypoint script to do the release. Run it like

#     ./release.sh <branch> <version>

# If the version is the same as the branch you can omit it.

# You may need to run the script with sudo on Linux. The only requirement for
# the script to work is Docker.

set -e
set -x

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

if [[ -z $2 ]]; then
    $2=$1
fi

docker run -t -v "$parent_path/release-$2":/root/release sympy/sympy-release "$@"
