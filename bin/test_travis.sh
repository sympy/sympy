#! /usr/bin/env bash

# Exit on error
set -e
# Echo each command
set -x

if [[ "${TEST_SPHINX}" == "true" ]]; then
    echo "Testing SPHINX"
    cd doc
    make html-errors
    make man
    make latex
    cd _build/latex
    export LATEXOPTIONS="-interaction=nonstopmode"
    make all
fi

if [[ "${TEST_SAGE}" == "true" ]]; then
    echo "Testing SAGE"
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
print('Testing ASCII')
import sympy
sympy.test('print')
EOF
    cd ..
    bin/doctest
    export LANG=$OLD_LANG
fi

if [[ "${TEST_DOCTESTS}" == "true" ]]; then
    cat << EOF | python
print('Testing DOCTESTS')
import sympy
if not sympy.doctest():
    raise Exception('Tests failed')
EOF
    cd ..
    bin/doctest doc/
fi

if [[ "${TEST_SLOW}" == "true" ]]; then
    cat << EOF | python
print('Testing SLOW')
import sympy
if not sympy.test(split='${SPLIT}', slow=True):
    # Travis times out if no activity is seen for 10 minutes. It also times
    # out if the whole tests run for more than 50 minutes.
    raise Exception('Tests failed')
EOF
fi

if [[ "${TEST_THEANO}" == "true" ]]; then
    cat << EOF | python
print('Testing THEANO')
import sympy
if not sympy.test('*theano*'):
    raise Exception('Tests failed')
EOF
fi

if [[ "${TEST_GMPY}" == "true" ]]; then
    cat << EOF | python
print('Testing GMPY')
import sympy
if not (sympy.test('sympy/polys/') and sympy.doctest('sympy/polys/')):
    raise Exception('Tests failed')
EOF
fi

if [[ "${TEST_MATPLOTLIB}" == "true" ]]; then
    cat << EOF | python
print('Testing MATPLOTLIB')
# Set matplotlib so that it works correctly in headless Travis. We have to do
# this here because it doesn't work after the sympy plotting module is
# imported.
import matplotlib
matplotlib.use("Agg")
import sympy
# Unfortunately, we have to use subprocess=False so that the above will be
# applied, so no hash randomization here.
if not (sympy.test('sympy/plotting', subprocess=False) and
    sympy.doctest('sympy/plotting', subprocess=False)):
    raise Exception('Tests failed')
EOF
fi

if [[ "${TEST_AUTOWRAP}" == "true" ]]; then
    cat << EOF | python
print('Testing AUTOWRAP')
import sympy
if not sympy.test('sympy/external/tests/test_autowrap.py'):
    raise Exception('Tests failed')
EOF
fi

if [[ "${TEST_SYMPY}" == "true" ]]; then
    cat << EOF | python
print('Testing SYMPY, split ${SPLIT}')
import sympy
if not sympy.test(split='${SPLIT}'):
   raise Exception('Tests failed')
EOF
fi
