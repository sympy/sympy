#! /usr/bin/env bash

set -o errexit

echo "Testing SPHINX"
cd doc
make html
make man
make latexpdf LATEXMKOPTS="-halt-on-error -xelatex -silent" || {
    echo "An error had occurred during the LaTeX build";
    tail -n 1000 _build/latex/*.log;
    sleep 1; # A guard against the CI running tail concurrently.
    exit 1;
}
