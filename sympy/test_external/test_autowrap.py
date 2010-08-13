try:
    import numpy
except:
    numpy = None

from sympy import symbols, Eq
from sympy.tensor import IndexedBase, Idx
from sympy.utilities.autowrap import autowrap
from sympy.utilities.pytest import XFAIL, skip

#
# test runners used by several language-backend combinations
#

def runtest_autowrap_twice(language, backend):
    a, b, c = symbols('a b c')
    f = autowrap((((a + b)/c)**5).expand(), language, backend)
    g = autowrap((((a + b)/c)**4).expand(), language, backend)

    # check that autowrap updates the module name.  Else, g gives the same as f
    assert f(1, -2, 1) == -1.0
    assert g(1, -2, 1) ==  1.0

def runtest_autowrap_trace(language, backend):
    if not numpy: skip()
    A = IndexedBase('A')
    n = symbols('n', integer=True)
    i = Idx('i', n)
    trace = autowrap(A[i, i], language, backend)
    assert trace(numpy.eye(100)) == 100

def runtest_autowrap_matrix_vector(language, backend):
    if not numpy: skip()
    A, x, y = map(IndexedBase, ['A', 'x', 'y'])
    n, m = symbols('n m', integer=True)
    i = Idx('i', m)
    j = Idx('j', n)
    expr = Eq(y[i], A[i, j]*x[j])
    mv = autowrap(expr, language, backend)

    # compare with numpy's dot product
    M = numpy.random.rand(10, 20)
    x = numpy.random.rand(20)
    y = numpy.dot(M, x)
    assert numpy.sum(numpy.abs(y - mv(M, x))) < 1e-13

def runtest_autowrap_matrix_matrix(language, backend):
    if not numpy: skip()
    A, B, C = map(IndexedBase, ['A', 'B', 'C'])
    n, m, d = symbols('n m d', integer=True)
    i = Idx('i', m)
    j = Idx('j', n)
    k = Idx('k', d)
    expr = Eq(C[i, j], A[i, k]*B[k, j])
    matmat = autowrap(expr, language, backend)

    # compare with numpy's dot product
    M1 = numpy.random.rand(10, 20)
    M2 = numpy.random.rand(20, 15)
    M3 = numpy.dot(M1, M2)
    assert numpy.sum(numpy.abs(M3 - matmat(M1, M2))) < 1e-13


#
# tests of language-backend combinations
#
try:
    import Cython
    has_Cython = True
except ImportError:
    has_Cython = False

try:
    import numpy.f2py
    has_f2py = True
except ImportError:
    has_f2py = False

# f2py

def test_wrap_twice_f95_f2py():
    if not has_f2py: skip()
    runtest_autowrap_twice('f95', 'f2py')

def test_autowrap_trace_f95_f2py():
    if not has_f2py: skip()
    runtest_autowrap_trace('f95', 'f2py')

def test_autowrap_matrix_vector_f95_f2py():
    if not has_f2py: skip()
    runtest_autowrap_matrix_vector('f95', 'f2py')

def test_autowrap_matrix_matrix_f95_f2py():
    if not has_f2py: skip()
    runtest_autowrap_matrix_matrix('f95', 'f2py')

# Cython

def test_wrap_twice_c_cython():
    if not has_Cython: skip()
    runtest_autowrap_twice('C', 'cython')

@XFAIL
def test_autowrap_trace_C_Cython():
    if not has_Cython: skip()
    runtest_autowrap_trace('C', 'cython')

@XFAIL
def test_autowrap_matrix_vector_C_cython():
    if not has_Cython: skip()
    runtest_autowrap_matrix_vector('C', 'cython')

@XFAIL
def test_autowrap_matrix_matrix_C_cython():
    if not has_Cython: skip()
    runtest_autowrap_matrix_matrix('C', 'cython')
