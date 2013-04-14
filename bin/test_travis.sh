#! /usr/bin/env bash

# Exit on error
set -e
# Echo each command
set -x

if [[ "${TEST_SPHINX}" == "true" ]]; then
    cd doc
    make html-errors
else
    # We change directories to make sure that we test the installed version of
    # sympy.
    mkdir empty
    cd empty
    cat << EOF | python
import sympy
t1=sympy.test()
t2=sympy.doctest()
if not (t1 and t2):
    raise Exception('Tests failed')
EOF

    cd ..
    bin/doctest doc/
fi
