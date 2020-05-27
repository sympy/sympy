from __future__ import absolute_import

from sympy.codegen import Assignment
from sympy.codegen.ast import none
from sympy.codegen.matrix_nodes import MatrixSolve
from sympy.core import Expr, Mod, symbols, Eq, Le, Gt, zoo, oo, Rational, Pow
from sympy.core.numbers import pi
from sympy.core.singleton import S
from sympy.functions import acos, KroneckerDelta, Piecewise, sign, sqrt
from sympy.logic import And, Or
from sympy.matrices import SparseMatrix, MatrixSymbol, Identity
from sympy.printing.pycode import (
    MpmathPrinter, NumPyPrinter, PythonCodePrinter, pycode, SciPyPrinter,
    SymPyPrinter
)
from sympy.testing.pytest import raises
from sympy.tensor import IndexedBase
from sympy.testing.pytest import skip
from sympy.external import import_module

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
    assert prntr.doprint(KroneckerDelta(x,y)) == '(1 if x == y else 0)'


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

    assert p.doprint(S.Exp1) == 'mpmath.e'
    assert p.doprint(S.Pi) == 'mpmath.pi'
    assert p.doprint(S.GoldenRatio) == 'mpmath.phi'
    assert p.doprint(S.EulerGamma) == 'mpmath.euler'
    assert p.doprint(S.NaN) == 'mpmath.nan'
    assert p.doprint(S.Infinity) == 'mpmath.inf'
    assert p.doprint(S.NegativeInfinity) == 'mpmath.ninf'


def test_NumPyPrinter():
    from sympy import (Lambda, ZeroMatrix, OneMatrix, FunctionMatrix,
        HadamardProduct, KroneckerProduct, Adjoint, DiagonalOf,
        DiagMatrix, DiagonalMatrix)
    from sympy.abc import a, b
    p = NumPyPrinter()
    assert p.doprint(sign(x)) == 'numpy.sign(x)'
    A = MatrixSymbol("A", 2, 2)
    B = MatrixSymbol("B", 2, 2)
    C = MatrixSymbol("C", 1, 5)
    D = MatrixSymbol("D", 3, 4)
    assert p.doprint(A**(-1)) == "numpy.linalg.inv(A)"
    assert p.doprint(A**5) == "numpy.linalg.matrix_power(A, 5)"
    assert p.doprint(Identity(3)) == "numpy.eye(3)"

    u = MatrixSymbol('x', 2, 1)
    v = MatrixSymbol('y', 2, 1)
    assert p.doprint(MatrixSolve(A, u)) == 'numpy.linalg.solve(A, x)'
    assert p.doprint(MatrixSolve(A, u) + v) == 'numpy.linalg.solve(A, x) + y'

    assert p.doprint(ZeroMatrix(2, 3)) == "numpy.zeros((2, 3))"
    assert p.doprint(OneMatrix(2, 3)) == "numpy.ones((2, 3))"
    assert p.doprint(FunctionMatrix(4, 5, Lambda((a, b), a + b))) == \
        "numpy.fromfunction(lambda a, b: a + b, (4, 5))"
    assert p.doprint(HadamardProduct(A, B)) == "numpy.multiply(A, B)"
    assert p.doprint(KroneckerProduct(A, B)) == "numpy.kron(A, B)"
    assert p.doprint(Adjoint(A)) == "numpy.conjugate(numpy.transpose(A))"
    assert p.doprint(DiagonalOf(A)) == "numpy.reshape(numpy.diag(A), (-1, 1))"
    assert p.doprint(DiagMatrix(C)) == "numpy.diagflat(C)"
    assert p.doprint(DiagonalMatrix(D)) == "numpy.multiply(D, numpy.eye(3, 4))"

    # Workaround for numpy negative integer power errors
    assert p.doprint(x**-1) == 'x**(-1.0)'
    assert p.doprint(x**-2) == 'x**(-2.0)'

    expr = Pow(2, -1, evaluate=False)
    assert p.doprint(expr) == "2**(-1.0)"

    assert p.doprint(S.Exp1) == 'numpy.e'
    assert p.doprint(S.Pi) == 'numpy.pi'
    assert p.doprint(S.EulerGamma) == 'numpy.euler_gamma'
    assert p.doprint(S.NaN) == 'numpy.nan'
    assert p.doprint(S.Infinity) == 'numpy.PINF'
    assert p.doprint(S.NegativeInfinity) == 'numpy.NINF'


def test_issue_18770():
    numpy = import_module('numpy')
    if not numpy:
        skip("numpy not installed.")

    from sympy import lambdify, Min, Max

    expr1 = Min(0.1*x + 3, x + 1, 0.5*x + 1)
    func = lambdify(x, expr1, "numpy")
    assert (func(numpy.linspace(0, 3, 3)) == [1.0 , 1.75, 2.5 ]).all()
    assert  func(4) == 3

    expr1 = Max(x**2 , x**3)
    func = lambdify(x,expr1, "numpy")
    assert (func(numpy.linspace(-1 , 2, 4)) == [1, 0, 1, 8] ).all()
    assert func(4) == 64


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

    assert p.doprint(S.GoldenRatio) == 'scipy.constants.golden_ratio'
    assert p.doprint(S.Pi) == 'scipy.constants.pi'
    assert p.doprint(S.Exp1) == 'numpy.e'


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

    prntr = PythonCodePrinter({'allow_unknown_functions': True})
    assert prntr.doprint(expr) == 'math.gamma(x)*math.gamma(y)/math.gamma(x + y)'

    prntr = MpmathPrinter()
    assert prntr.doprint(expr) ==  'mpmath.beta(x, y)'

def test_airy():
    from sympy import airyai, airybi

    expr1 = airyai(x)
    expr2 = airybi(x)

    prntr = SciPyPrinter()
    assert prntr.doprint(expr1) == 'scipy.special.airy(x)[0]'
    assert prntr.doprint(expr2) == 'scipy.special.airy(x)[2]'

    prntr = NumPyPrinter()
    assert prntr.doprint(expr1) == '  # Not supported in Python with NumPy:\n  # airyai\nairyai(x)'
    assert prntr.doprint(expr2) == '  # Not supported in Python with NumPy:\n  # airybi\nairybi(x)'

    prntr = PythonCodePrinter()
    assert prntr.doprint(expr1) == '  # Not supported in Python:\n  # airyai\nairyai(x)'
    assert prntr.doprint(expr2) == '  # Not supported in Python:\n  # airybi\nairybi(x)'

def test_airy_prime():
    from sympy import airyaiprime, airybiprime

    expr1 = airyaiprime(x)
    expr2 = airybiprime(x)

    prntr = SciPyPrinter()
    assert prntr.doprint(expr1) == 'scipy.special.airy(x)[1]'
    assert prntr.doprint(expr2) == 'scipy.special.airy(x)[3]'

    prntr = NumPyPrinter()
    assert prntr.doprint(expr1) == '  # Not supported in Python with NumPy:\n  # airyaiprime\nairyaiprime(x)'
    assert prntr.doprint(expr2) == '  # Not supported in Python with NumPy:\n  # airybiprime\nairybiprime(x)'

    prntr = PythonCodePrinter()
    assert prntr.doprint(expr1) == '  # Not supported in Python:\n  # airyaiprime\nairyaiprime(x)'
    assert prntr.doprint(expr2) == '  # Not supported in Python:\n  # airybiprime\nairybiprime(x)'
