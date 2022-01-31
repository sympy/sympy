#! /usr/bin/env bash

set -o errexit

echo "Testing SPHINX"
cd doc
make html
make man
make latexpdf LATEXMKOPTS="-halt-on-error -xelatex -silent" || {
    echo "An error had occured during the LaTeX build";
    tail -n 1000 _build/latex/*.log;
    sleep 1; # A guard against travis running tail concurrently.
    exit 1;
}
