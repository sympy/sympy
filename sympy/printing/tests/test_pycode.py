# -*- coding: utf-8 -*-
from __future__ import absolute_import

from sympy.codegen import Assignment
from sympy.codegen.ast import none
from sympy.codegen.matrix_nodes import MatrixSolve
from sympy.core import Expr, Mod, symbols, Eq, Le, Gt, zoo, oo, Rational
from sympy.core.numbers import pi
from sympy.functions import acos, erf, gamma,Piecewise , jacobi, loggamma, legendre, hermite, laguerre, digamma, besselk
from sympy.functions import gegenbauer, chebyshevt, chebyshevu, RisingFactorial, besselj, factorial, bessely, besseli
from sympy.functions import acosh, asin, asinh, atan, atanh, assoc_laguerre, sign, sqrt, beta, fresnels, fresnelc, erfc
from sympy.logic import And, Or
from sympy.codegen.cfunctions import exp2
from sympy.matrices import SparseMatrix, MatrixSymbol, Identity
from sympy.printing.pycode import (
    MpmathPrinter, NumPyPrinter, PythonCodePrinter, pycode, SciPyPrinter,
    SymPyPrinter
)
from sympy.utilities.pytest import raises
from sympy.tensor import IndexedBase

x, y, z = symbols('x y z')
p = IndexedBase("p")

def test_PythonCodePrinter():
    prntr = PythonCodePrinter()

    assert not prntr.module_imports

    assert prntr.doprint(x**y) == 'x**y'
    assert prntr.doprint(Mod(x, 2)) == 'x % 2'
    assert prntr.doprint(And(x, y)) == 'x and y'
    assert prntr.doprint(Or(x, y)) == 'x or y'
    assert not prntr.module_imports

    assert prntr.doprint(pi) == 'math.pi'
    assert prntr.module_imports == {'math': {'pi'}}

    assert prntr.doprint(x**Rational(1, 2)) == 'math.sqrt(x)'
    assert prntr.doprint(sqrt(x)) == 'math.sqrt(x)'
    assert prntr.module_imports == {'math': {'pi', 'sqrt'}}

    assert prntr.doprint(acos(x)) == 'math.acos(x)'
    assert prntr.doprint(Assignment(x, 2)) == 'x = 2'
    assert prntr.doprint(Piecewise((1, Eq(x, 0)),
                        (2, x>6))) == '((1) if (x == 0) else (2) if (x > 6) else None)'
    assert prntr.doprint(Piecewise((2, Le(x, 0)),
                        (3, Gt(x, 0)), evaluate=False)) == '((2) if (x <= 0) else'\
                                                        ' (3) if (x > 0) else None)'
    assert prntr.doprint(sign(x)) == '(0.0 if x == 0 else math.copysign(1, x))'
    assert prntr.doprint(p[0, 1]) == 'p[0, 1]'


def test_PythonCodePrinter_standard():
    import sys
    prntr = PythonCodePrinter({'standard':None})

    python_version = sys.version_info.major
    if python_version == 2:
        assert prntr.standard == 'python2'
    if python_version == 3:
        assert prntr.standard == 'python3'

    raises(ValueError, lambda: PythonCodePrinter({'standard':'python4'}))

def test_MpmathPrinter():
    p = MpmathPrinter()
    assert p.doprint(sign(x)) == 'mpmath.sign(x)'
    assert p.doprint(Rational(1, 2)) == 'mpmath.mpf(1)/mpmath.mpf(2)'
    assert p.doprint(loggamma(x)) == 'mpmath.loggamma(x)'
    assert p.doprint(fresnels(x)) == 'mpmath.fresnels(x)'
    assert p.doprint(fresnelc(x)) == 'mpmath.fresnelc(x)'

def test_NumPyPrinter():
    p = NumPyPrinter()
    assert p.doprint(sign(x)) == 'numpy.sign(x)'
    assert p.doprint(acos(x)) == 'numpy.arccos(x)'
    assert p.doprint(acosh(x)) == 'numpy.arccosh(x)'
    assert p.doprint(asin(x)) == 'numpy.arcsin(x)'
    assert p.doprint(asinh(x)) == 'numpy.arcshinh(x)'
    assert p.doprint(atan(x)) == 'numpy.arctan(x)'
    assert p.doprint(atanh(x)) == 'numpy.arctanh(x)'
    assert p.doprint(exp2(x)) == 'numpy.exp2(x)'
    A = MatrixSymbol("A", 2, 2)
    assert p.doprint(A**(-1)) == "numpy.linalg.inv(A)"
    assert p.doprint(A**5) == "numpy.linalg.matrix_power(A, 5)"
    assert p.doprint(Identity(3)) == "numpy.eye(3)"

    u = MatrixSymbol('x', 2, 1)
    v = MatrixSymbol('y', 2, 1)
    assert p.doprint(MatrixSolve(A, u)) == 'numpy.linalg.solve(A, x)'
    assert p.doprint(MatrixSolve(A, u) + v) == 'numpy.linalg.solve(A, x) + y'
    # Workaround for numpy negative integer power errors
    assert p.doprint(x**-1) == 'x**(-1.0)'
    assert p.doprint(x**-2) == 'x**(-2.0)'


def test_SciPyPrinter():
    p = SciPyPrinter()
    expr = acos(x)
    assert 'numpy' not in p.module_imports
    assert p.doprint(expr) == 'numpy.arccos(x)'
    assert 'numpy' in p.module_imports
    assert not any(m.startswith('scipy') for m in p.module_imports)
    smat = SparseMatrix(2, 5, {(0, 1): 3})
    assert p.doprint(smat) == 'scipy.sparse.coo_matrix([3], ([0], [1]), shape=(2, 5))'
    assert 'scipy.sparse' in p.module_imports
    assert p.doprint(erf(x)) == 'scipy.special.erf(x)'
    assert p.doprint(erfc(x)) == 'scipy.special.erfc(x)'
    assert p.doprint(besselj(x, y)) == 'scipy.special.jv(x, y)'
    assert p.doprint(bessely(x, y)) == 'scipy.special.yv(x, y)'
    assert p.doprint(besseli(x, y)) == 'scipy.special.iv(x, y)'
    assert p.doprint(besselk(x, y)) == 'scipy.special.kv(x, y)'
    assert p.doprint(factorial(x)) == 'scipy.special.factorial(x, exact=true)'
    assert p.doprint(gamma(x)) == 'scipy.special.gamma(x)'
    assert p.doprint(loggamma(x)) == 'scipy.special.gammaln(x)'
    assert p.doprint(digamma(x)) == 'scipy.special.psi(x)'
    assert p.doprint(RisingFactorial(x, y)) == 'scipy.special.poch(x, y)'
    assert p.doprint(jacobi(1, x, y, z)) == 'scipy.special.eval_jacobi(1, x, y, z)'
    assert p.doprint(gegenbauer(x, y, z)) == 'scipy.special.eval_gegenbauer(x, y, z)'
    assert p.doprint(chebyshevt(x, y)) == 'scipy.special.eval_chebyt(x, y)'
    assert p.doprint(chebyshevu(x, y)) == 'scipy.special.eval_chebyu(x, y)'
    assert p.doprint(legendre(x, y)) == 'scipy.special.eval_legendre(x, y)'
    assert p.doprint(hermite(x, y)) == 'scipy.special.eval_hermite(x, y)'
    assert p.doprint(laguerre(x, 0, y)) == 'scipy.special.eval_laguerre(x, y)'
    assert p.doprint(assoc_laguerre(x,y, z)) == 'scipy.special.eval_genlaguerre(x, y, z)'
    assert p.doprint(beta(x, y)) == 'scipy.special.beta(x, y)'


def test_pycode_reserved_words():
    s1, s2 = symbols('if else')
    raises(ValueError, lambda: pycode(s1 + s2, error_on_reserved=True))
    py_str = pycode(s1 + s2)
    assert py_str in ('else_ + if_', 'if_ + else_')


def test_sqrt():
    prntr = PythonCodePrinter()
    assert prntr._print_Pow(sqrt(x), rational=False) == 'math.sqrt(x)'
    assert prntr._print_Pow(1/sqrt(x), rational=False) == '1/math.sqrt(x)'

    prntr = PythonCodePrinter({'standard' : 'python2'})
    assert prntr._print_Pow(sqrt(x), rational=True) == 'x**(1./2.)'
    assert prntr._print_Pow(1/sqrt(x), rational=True) == 'x**(-1./2.)'

    prntr = PythonCodePrinter({'standard' : 'python3'})
    assert prntr._print_Pow(sqrt(x), rational=True) == 'x**(1/2)'
    assert prntr._print_Pow(1/sqrt(x), rational=True) == 'x**(-1/2)'

    prntr = MpmathPrinter()
    assert prntr._print_Pow(sqrt(x), rational=False) == 'mpmath.sqrt(x)'
    assert prntr._print_Pow(sqrt(x), rational=True) == \
        "x**(mpmath.mpf(1)/mpmath.mpf(2))"

    prntr = NumPyPrinter()
    assert prntr._print_Pow(sqrt(x), rational=False) == 'numpy.sqrt(x)'
    assert prntr._print_Pow(sqrt(x), rational=True) == 'x**(1/2)'

    prntr = SciPyPrinter()
    assert prntr._print_Pow(sqrt(x), rational=False) == 'numpy.sqrt(x)'
    assert prntr._print_Pow(sqrt(x), rational=True) == 'x**(1/2)'

    prntr = SymPyPrinter()
    assert prntr._print_Pow(sqrt(x), rational=False) == 'sympy.sqrt(x)'
    assert prntr._print_Pow(sqrt(x), rational=True) == 'x**(1/2)'


class CustomPrintedObject(Expr):
    def _numpycode(self, printer):
        return 'numpy'

    def _mpmathcode(self, printer):
        return 'mpmath'


def test_printmethod():
    obj = CustomPrintedObject()
    assert NumPyPrinter().doprint(obj) == 'numpy'
    assert MpmathPrinter().doprint(obj) == 'mpmath'


def test_codegen_ast_nodes():
    assert pycode(none) == 'None'


def test_issue_14283():
    prntr = PythonCodePrinter()

    assert prntr.doprint(zoo) == "float('nan')"
    assert prntr.doprint(-oo) == "float('-inf')"

def test_NumPyPrinter_print_seq():
    n = NumPyPrinter()

    assert n._print_seq(range(2)) == '(0, 1,)'


def test_issue_16535_16536():
    from sympy import lowergamma, uppergamma

    a = symbols('a')
    expr1 = lowergamma(a, x)
    expr2 = uppergamma(a, x)

    prntr = SciPyPrinter()
    assert prntr.doprint(expr1) == 'scipy.special.gamma(a)*scipy.special.gammainc(a, x)'
    assert prntr.doprint(expr2) == 'scipy.special.gamma(a)*scipy.special.gammaincc(a, x)'

    prntr = NumPyPrinter()
    assert prntr.doprint(expr1) == '  # Not supported in Python with NumPy:\n  # lowergamma\nlowergamma(a, x)'
    assert prntr.doprint(expr2) == '  # Not supported in Python with NumPy:\n  # uppergamma\nuppergamma(a, x)'

    prntr = PythonCodePrinter()
    assert prntr.doprint(expr1) == '  # Not supported in Python:\n  # lowergamma\nlowergamma(a, x)'
    assert prntr.doprint(expr2) == '  # Not supported in Python:\n  # uppergamma\nuppergamma(a, x)'


def test_fresnel_integrals():
    from sympy import fresnelc, fresnels

    expr1 = fresnelc(x)
    expr2 = fresnels(x)

    prntr = SciPyPrinter()
    assert prntr.doprint(expr1) == 'scipy.special.fresnel(x)[1]'
    assert prntr.doprint(expr2) == 'scipy.special.fresnel(x)[0]'

    prntr = NumPyPrinter()
    assert prntr.doprint(expr1) == '  # Not supported in Python with NumPy:\n  # fresnelc\nfresnelc(x)'
    assert prntr.doprint(expr2) == '  # Not supported in Python with NumPy:\n  # fresnels\nfresnels(x)'

    prntr = PythonCodePrinter()
    assert prntr.doprint(expr1) == '  # Not supported in Python:\n  # fresnelc\nfresnelc(x)'
    assert prntr.doprint(expr2) == '  # Not supported in Python:\n  # fresnels\nfresnels(x)'

    prntr = MpmathPrinter()
    assert prntr.doprint(expr1) == 'mpmath.fresnelc(x)'
    assert prntr.doprint(expr2) == 'mpmath.fresnels(x)'


def test_beta():
    from sympy import beta

    expr = beta(x, y)

    prntr = SciPyPrinter()
    assert prntr.doprint(expr) == 'scipy.special.beta(x, y)'

    prntr = NumPyPrinter()
    assert prntr.doprint(expr) == 'math.gamma(x)*math.gamma(y)/math.gamma(x + y)'

    prntr = PythonCodePrinter()
    assert prntr.doprint(expr) == 'math.gamma(x)*math.gamma(y)/math.gamma(x + y)'

    prntr = MpmathPrinter()
    assert prntr.doprint(expr) ==  'mpmath.beta(x, y)'
