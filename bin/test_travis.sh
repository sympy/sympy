#! /usr/bin/env bash

# Exit on error
set -e
# Echo each command
set -x

if [[ "${TEST_SPHINX}" == "true" ]]; then
    cd doc
    make html-errors
    make man
    make latex
    cd _build/latex
    export LATEXOPTIONS="-interaction=nonstopmode"
    make all
fi

if [[ "${TEST_SAGE}" == "true" ]]; then
    sage -v
    sage -python bin/test sympy/external/tests/test_sage.py
fi

# We change directories to make sure that we test the installed version of
# sympy.
mkdir empty
cd empty

if [[ "${TEST_ASCII}" == "true" ]]; then
    export OLD_LANG=$LANG
    export LANG=c
    cat <<EOF | python
import sympy
sympy.test('print')
EOF
    cd ..
    bin/doctest
    export LANG=$OLD_LANG
fi

if [[ "${TEST_DOCTESTS}" == "true" ]]; then
    cat << EOF | python
import sympy
if not sympy.doctest():
    raise Exception('Tests failed')
EOF
    cd ..
    bin/doctest doc/
fi

if [[ "${TEST_SLOW}" == "true" ]]; then
    cat << EOF | python
import sympy
if not sympy.test(split='${SPLIT}', slow=True):
    # Travis times out if no activity is seen for 10 minutes. It also times
    # out if the whole tests run for more than 50 minutes.
    raise Exception('Tests failed')
EOF
fi

if [[ "${TEST_THEANO}" == "true" ]]; then
    cat << EOF | python
import sympy
if not sympy.test('*theano*'):
    raise Exception('Tests failed')
EOF
fi

if [[ "${TEST_GMPY}" == "true" ]]; then
    cat << EOF | python
import sympy
if not sympy.test('sympy/polys/') and sympy.doctest('sympy/polys/'):
    raise Exception('Tests failed')
EOF
fi

if [[ "${TEST_MATPLOTLIB}" == "true" ]]; then
    cat << EOF | python
import sympy
if not sympy.test('sympy/plotting') and sympy.doctest('sympy/plotting'):
    raise Exception('Tests failed')
EOF
fi

if [[ "${TEST_AUTOWRAP}" == "true" ]]; then
    cat << EOF | python
import sympy
if not sympy.test('sympy/external/tests/test_autowrap.py'):
    raise Exception('Tests failed')
EOF
fi

if [[ "${TEST_SYMPY}" == "true" ]]; then
    cat << EOF | python
import sympy
if not sympy.test(split='${SPLIT}'):
   raise Exception('Tests failed')
EOF
fi
