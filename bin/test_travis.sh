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

    if [[ "${TEST_DOCTESTS}" == "true" ]]; then
        cat << EOF | python
import sympy
if not sympy.doctest():
    raise Exception('Tests failed')
EOF
        cd ..
        bin/doctest doc/
    elif [[ "${TEST_SLOW}" == "true" ]]; then
        cat << EOF | python
import sympy
if not sympy.test(slow=True, timeout=550):
    # Travis times out if no activity is seen for 10 minutes
    raise Exception('Tests failed')
EOF
    else
        cat << EOF | python
import sympy
if not sympy.test():
    raise Exception('Tests failed')
EOF
        fi
fi
