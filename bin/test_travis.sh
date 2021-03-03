#! /usr/bin/env bash

# Exit on error
set -e
# Echo each command
set -x

if [[ "${TEST_FLAKE8}" == "true" ]]; then
    flake8 sympy;
fi

if [[ "${TEST_SETUP}" == "true" ]]; then
    python bin/test_setup.py
fi

if [[ "${TEST_PY2_IMPORT}" == "true" ]]; then
    python bin/test_py2_import.py
fi

if [[ "${TEST_SPHINX}" == "true" ]]; then
    echo "Testing SPHINX"
    cd doc
    make html
    make man
    make latex
    cd _build/latex
    export LATEXMKOPTS="-halt-on-error -xelatex -silent"
    make all || {
        echo "An error had occured during the LaTeX build";
        tail -n 1000 *.log;
        sleep 1; # A guard against travis running tail concurrently.
        exit -1;
    }
fi

if [[ "${TEST_SAGE}" == "true" ]]; then
    echo "Testing SAGE"
    source deactivate
    source activate sage
    sage -v
    sage -python bin/test sympy/external/tests/test_sage.py
    PYTHONPATH=. sage -t sympy/external/tests/test_sage.py
    export MPMATH_NOSAGE=1
    source deactivate
    source activate test-environment
fi

if [[ -n "${TEST_OPT_DEPENDENCY}" ]]; then
    python bin/test_external_imports.py
    python bin/test_submodule_imports.py
    python bin/test_executable.py
fi

# We change directories to make sure that we test the installed version of
# sympy.
mkdir empty
cd empty

if [[ "${TEST_COVERAGE}" == "true" ]]; then
    rm -f $TRAVIS_BUILD_DIR/.coverage.* $TRAVIS_BUILD_DIR/.coverage
    cat << EOF | python
import distutils.sysconfig
import os

with open(os.path.join(distutils.sysconfig.get_python_lib(), 'coverage.pth'), 'w') as pth:
    pth.write('import sys; exec(%r)\n' % '''\
try:
    import coverage
except ImportError:
    pass
else:
    coverage.process_startup()
''')
EOF
    export COVERAGE_PROCESS_START=$TRAVIS_BUILD_DIR/coveragerc_travis
fi

if [[ "${TEST_ASCII}" == "true" ]]; then
    # Force Python to act like pre-3.7 where LC_ALL=C causes
    # UnicodeEncodeErrors. Once the lowest Python version we support is 3.7,
    # we can consider dropping this test entirely. See PEP 538.
    export PYTHONIOENCODING=ascii:strict
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


# Test tensorflow 1 support (which is a separate conda environment)
if [[ "${TEST_TENSORFLOW_1}" == "true" ]]; then
    cat << EOF | python
print('Testing optional dependencies')

import sympy
test_list = [
    "tensorflow",
]

doctest_list = [
    "tensorflow",
]

if not (sympy.test(*test_list) and sympy.doctest(*doctest_list)):
    raise Exception('Tests failed')
EOF

    # lambdify with tensorflow and numexpr is tested here

    # TODO: Generate these tests automatically
elif [[ -n "${TEST_OPT_DEPENDENCY}" ]]; then
    cat << EOF | python
print('Testing optional dependencies')

import sympy
test_list = [
    # numpy
    '*numpy*',
    'sympy/core/',
    'sympy/matrices/',
    'sympy/physics/quantum/',
    'sympy/utilities/tests/test_lambdify.py',

    # scipy
    '*scipy*',

    # llvmlite
    '*llvm*',

    # theano
    '*theano*',

    # gmpy
    'polys',

    # autowrap
    '*autowrap*',

    # ipython
    '*ipython*',

    # antlr, lfortran, clang
    'sympy/parsing/',

    # matchpy
    '*rubi*',

    # codegen
    'sympy/codegen/',
    'sympy/utilities/tests/test_codegen',
    'sympy/utilities/_compilation/tests/test_compilation',

    # cloudpickle
    'pickling',

    # pycosat
    'sympy/logic',
    'sympy/assumptions',

    #stats
    'sympy/stats',

]

blacklist = [
    'sympy/physics/quantum/tests/test_circuitplot.py',
]

doctest_list = [
    # numpy
    'sympy/matrices/',
    'sympy/utilities/lambdify.py',

    # scipy
    '*scipy*',

    # llvmlite
    '*llvm*',

    # theano
    '*theano*',

    # gmpy
    'polys',

    # autowrap
    '*autowrap*',

    # ipython
    '*ipython*',

    # antlr, lfortran, clang
    'sympy/parsing/',

    # matchpy
    '*rubi*',

    # codegen
    'sympy/codegen/',

    # pycosat
    'sympy/logic',
    'sympy/assumptions',

    #stats
    'sympy/stats',

]

if not (sympy.test(*test_list, blacklist=blacklist) and sympy.doctest(*doctest_list)):
    raise Exception('Tests failed')
EOF
    cd ..
    bin/doctest doc/src/modules/numeric-computation.rst
fi

# This is separate because it needs to be run with subprocess=False
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
if [[ "${TEST_COVERAGE}" == "true" ]]; then
    unset COVERAGE_PROCESS_START
fi

if [[ "${TEST_EXAMPLES}" == "true" ]]; then
    # No need to change directory if executed after the rst doctest
    examples/all.py -q
fi
