#! /usr/bin/env bash

# Exit on error
set -e
# Echo each command
set -x

if [[ "${TEST_SPHINX}" == "true" ]]; then
    cd doc
    make html-errors
    make clean
    make latex
    cd _build/latex
    export LATEXOPTIONS="-interaction=nonstopmode"
    make all
else
    # We change directories to make sure that we test the installed version of
    # sympy.
    mkdir empty
    cd empty

    if [[ "${TEST_DOCTESTS}" == "true" ]]; then
        cat << EOF | python
import sympy
if not sympy.doctest():
    raise Exception('Tests failed')
EOF
        cd ..
        bin/doctest doc/
    else
        cat << EOF | python
import sympy
if not sympy.test():
    raise Exception('Tests failed')
EOF
        fi
fi
