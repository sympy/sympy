#!/bin/bash
set -e
set -x

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

docker build -f Dockerfile . -t sympy-release
docker run -v "$parent_path/release-$1":/home/release sympy-release "$@"
