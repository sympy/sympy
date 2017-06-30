#!/bin/bash
set -e
set -x

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

docker build -f Dockerfile-base . -t sympy-release-base
docker build -f Dockerfile --no-cache . -t sympy-release
docker run -v "$parent_path/release-$1" sympy-release "$@"
