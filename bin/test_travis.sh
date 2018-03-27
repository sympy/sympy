#! /usr/bin/env bash

# Exit on error
set -e
# Echo each command
set -x

if [[ "${TEST_SETUP}" == "true" ]]; then
    python bin/test_setup.py
fi

if [[ "${TEST_SPHINX}" == "true" ]]; then
    echo "Testing SPHINX"
    cd doc
    make html
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
    sage -t sympy/external/tests/test_sage.py
fi

# We change directories to make sure that we test the installed version of
# sympy.
mkdir empty
cd empty

if [[ "${TEST_ASCII}" == "true" ]]; then
    export OLD_LC_ALL=$LC_ALL
    export LC_ALL=C
    cat <<EOF | python
print('Testing ASCII')
try:
    print(u'\u2713')
except UnicodeEncodeError:
    pass
else:
    raise Exception('Not an ASCII-only environment')
import sympy
if not (sympy.test('print') and sympy.doctest()):
    raise Exception('Tests failed')
EOF
    export LC_ALL=$OLD_LC_ALL
fi

if [[ "${TEST_DOCTESTS}" == "true" ]]; then
    # -We:invalid makes invalid escape sequences error in Python 3.6. See
    # -#12028.
    cat << EOF | python -We:invalid
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
if not sympy.test(split='${SPLIT}', slow=True, verbose=True):
    raise Exception('Tests failed')
EOF
fi

# lambdify with tensorflow and numexpr is tested here
if [[ "${TEST_OPT_DEPENDENCY}" == *"numpy"* ]]; then
    cat << EOF | python
print('Testing NUMPY')
import sympy
if not (sympy.test('*numpy*', 'sympy/core/tests/test_numbers.py',
                   'sympy/matrices/', 'sympy/physics/quantum/',
                   'sympy/core/tests/test_sympify.py',
                   'sympy/utilities/tests/test_lambdify.py',
                   blacklist=['sympy/physics/quantum/tests/test_circuitplot.py'])
        and sympy.doctest('sympy/matrices/', 'sympy/utilities/lambdify.py')):
    raise Exception('Tests failed')
EOF
fi

if [[ "${TEST_OPT_DEPENDENCY}" == *"scipy"* ]]; then
    cat << EOF | python
print('Testing SCIPY')
import sympy
# scipy matrices are tested in numpy testing
if not sympy.test('sympy/external/tests/test_scipy.py'):
    raise Exception('Tests failed')
EOF
fi

if [[ "${TEST_OPT_DEPENDENCY}" == *"llvmlite"* ]]; then
    cat << EOF | python
print('Testing LLVMJIT')
import sympy
if not (sympy.test('sympy/printing/tests/test_llvmjit.py')
        and sympy.doctest('sympy/printing/llvmjitcode.py')):
    raise Exception('Tests failed')
EOF
fi

if [[ "${TEST_OPT_DEPENDENCY}" == *"theano"* ]]; then
    cat << EOF | python
print('Testing THEANO')
import sympy
if not sympy.test('*theano*'):
    raise Exception('Tests failed')
EOF
fi

if [[ "${TEST_OPT_DEPENDENCY}" == *"gmpy"* ]]; then
    cat << EOF | python
print('Testing GMPY')
import sympy
if not (sympy.test('sympy/polys/') and sympy.doctest('sympy/polys/')):
    raise Exception('Tests failed')
EOF
fi

if [[ "${TEST_OPT_DEPENDENCY}" == *"matplotlib"* ]]; then
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
if not (sympy.test('sympy/plotting', 'sympy/physics/quantum/tests/test_circuitplot.py',
    subprocess=False) and sympy.doctest('sympy/plotting', subprocess=False)):
    raise Exception('Tests failed')
EOF
fi

if [[ "${TEST_OPT_DEPENDENCY}" == *"autowrap"* ]]; then
    cat << EOF | python
print('Testing AUTOWRAP')
import sympy
if not (sympy.test('sympy/external/tests/test_autowrap.py')
        and sympy.doctest('sympy/utilities/autowrap.py')):
    raise Exception('Tests failed')
EOF
fi

if [[ "${TEST_OPT_DEPENDENCY}" == *"ipython"* ]]; then
    cat << EOF | python
print('Testing IPYTHON')
import sympy
if not sympy.test('*ipython*'):
    raise Exception('Tests failed')
EOF
fi

if [[ "${TEST_SYMPY}" == "true" ]]; then
    # -We:invalid makes invalid escape sequences error in Python 3.6. See
    # -#12028.
    cat << EOF | python -We:invalid
print('Testing SYMPY, split ${SPLIT}')
import sympy
if not sympy.test(split='${SPLIT}'):
   raise Exception('Tests failed')
EOF
fi


if [[ "${TEST_OPT_DEPENDENCY}" == *"symengine"* ]]; then
    export USE_SYMENGINE=1
    cat << EOF | python
print('Testing SYMENGINE')
import sympy
if not sympy.test('sympy/physics/mechanics'):
    raise Exception('Tests failed')
if not sympy.test('sympy/liealgebras'):
    raise Exception('Tests failed')
EOF
    unset USE_SYMENGINE
fi

if [[ "${TEST_OPT_DEPENDENCY}" == *"antlr"* ]]; then
    cat << EOF | python -We:invalid
print('Testing ANTLR')
import sympy
if not sympy.test('sympy/parsing/tests/test_latex'):
    raise Exception('Tests failed')
EOF
fi
