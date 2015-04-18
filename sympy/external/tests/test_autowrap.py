from sympy import symbols, Eq
from sympy.external import import_module
from sympy.tensor import IndexedBase, Idx, EinsteinSum
from sympy.utilities.autowrap import autowrap, ufuncify, CodeWrapError
from sympy.utilities.pytest import skip

numpy = import_module('numpy', min_module_version='1.6.1')
Cython = import_module('Cython', min_module_version='0.15.1')
f2py = import_module('numpy.f2py', __import__kwargs={'fromlist': ['f2py']})

f2pyworks = False
if f2py:
    try:
        autowrap(symbols('x'), 'f95', 'f2py')
    except (CodeWrapError, ImportError, OSError):
        f2pyworks = False
    else:
        f2pyworks = True

a, b, c = symbols('a b c')
n, m, d = symbols('n m d', integer=True)
A, B, C = symbols('A B C', cls=IndexedBase)
i = Idx('i', m)
j = Idx('j', n)
k = Idx('k', d)


def has_module(module):
    """
    Return True if module exists, otherwise run skip().

    module should be a string.
    """
    # To give a string of the module name to skip(), this function takes a
    # string.  So we don't waste time running import_module() more than once,
    # just map the three modules tested here in this dict.
    modnames = {'numpy': numpy, 'Cython': Cython, 'f2py': f2py}

    if modnames[module]:
        if module == 'f2py' and not f2pyworks:
            skip("Couldn't run f2py.")
        return True
    skip("Couldn't import %s." % module)

#
# test runners used by several language-backend combinations
#


def runtest_autowrap_twice(language, backend):
    f = autowrap((((a + b)/c)**5).expand(), language, backend)
    g = autowrap((((a + b)/c)**4).expand(), language, backend)

    # check that autowrap updates the module name.  Else, g gives the same as f
    assert f(1, -2, 1) == -1.0
    assert g(1, -2, 1) == 1.0


def runtest_autowrap_trace(language, backend):
    has_module('numpy')
    trace = autowrap(EinsteinSum(A[i, i]), language, backend)
    assert trace(numpy.eye(100)) == 100


def runtest_autowrap_matrix_vector(language, backend):
    has_module('numpy')
    x, y = symbols('x y', cls=IndexedBase)
    expr = Eq(y[i], EinsteinSum(A[i, j]*x[j]))
    mv = autowrap(expr, language, backend)

    # compare with numpy's dot product
    M = numpy.random.rand(10, 20)
    x = numpy.random.rand(20)
    y = numpy.dot(M, x)
    assert numpy.sum(numpy.abs(y - mv(M, x))) < 1e-13


def runtest_autowrap_matrix_matrix(language, backend):
    has_module('numpy')
    expr = Eq(C[i, j], EinsteinSum(A[i, k]*B[k, j]))
    matmat = autowrap(expr, language, backend)

    # compare with numpy's dot product
    M1 = numpy.random.rand(10, 20)
    M2 = numpy.random.rand(20, 15)
    M3 = numpy.dot(M1, M2)
    assert numpy.sum(numpy.abs(M3 - matmat(M1, M2))) < 1e-13


def runtest_ufuncify(language, backend):
    has_module('numpy')
    a, b, c = symbols('a b c')
    fabc = ufuncify([a, b, c], a*b + c, backend=backend)
    facb = ufuncify([a, c, b], a*b + c, backend=backend)
    grid = numpy.linspace(-2, 2, 50)
    b = numpy.linspace(-5, 4, 50)
    c = numpy.linspace(-1, 1, 50)
    expected = grid*b + c
    numpy.testing.assert_allclose(fabc(grid, b, c), expected)
    numpy.testing.assert_allclose(facb(grid, c, b), expected)

#
# tests of language-backend combinations
#

# f2py


def test_wrap_twice_f95_f2py():
    has_module('f2py')
    runtest_autowrap_twice('f95', 'f2py')


def test_autowrap_trace_f95_f2py():
    has_module('f2py')
    runtest_autowrap_trace('f95', 'f2py')


def test_autowrap_matrix_vector_f95_f2py():
    has_module('f2py')
    runtest_autowrap_matrix_vector('f95', 'f2py')


def test_autowrap_matrix_matrix_f95_f2py():
    has_module('f2py')
    runtest_autowrap_matrix_matrix('f95', 'f2py')


def test_ufuncify_f95_f2py():
    has_module('f2py')
    runtest_ufuncify('f95', 'f2py')


# Cython

def test_wrap_twice_c_cython():
    has_module('Cython')
    runtest_autowrap_twice('C', 'cython')


def test_autowrap_trace_C_Cython():
    has_module('Cython')
    runtest_autowrap_trace('C', 'cython')


def test_autowrap_matrix_vector_C_cython():
    has_module('Cython')
    runtest_autowrap_matrix_vector('C', 'cython')


def test_autowrap_matrix_matrix_C_cython():
    has_module('Cython')
    runtest_autowrap_matrix_matrix('C', 'cython')


def test_ufuncify_C_Cython():
    has_module('Cython')
    runtest_ufuncify('C', 'cython')


# Numpy

def test_ufuncify_numpy():
    # This test doesn't use Cython, but if Cython works, then there is a valid
    # C compiler, which is needed.
    has_module('Cython')
    runtest_ufuncify('C', 'numpy')
