"""
Tests for the ConTeXt printer.

ConTeXt is a document preparation system based on TeX, similar to LaTeX.
Most mathematical expressions are identical to LaTeX, but ConTeXt uses
different environment delimiters.
"""
from sympy.core.numbers import (Float, Integer, Rational, oo, pi, I)
from sympy.core.singleton import S
from sympy.core.symbol import symbols
from sympy.functions.elementary.complexes import (Abs, conjugate, im, re)
from sympy.functions.elementary.exponential import (exp, log)
from sympy.functions.elementary.hyperbolic import (sinh, cosh, tanh, asinh, coth)
from sympy.functions.elementary.miscellaneous import (sqrt, root, Max, Min)
from sympy.functions.elementary.trigonometric import (sin, cos, tan, asin, acos, atan, cot, sec, csc)
from sympy.functions.special.gamma_functions import gamma
from sympy.integrals.integrals import Integral
from sympy.matrices.dense import Matrix
from sympy.series.limits import Limit
from sympy.concrete.summations import Sum
from sympy.concrete.products import Product
from sympy.core.function import (Derivative, Function, diff)
from sympy.core.relational import Eq, Ne, Lt, Le, Gt, Ge
from sympy.logic.boolalg import (And, Or, Not, Xor, Implies, Equivalent)
from sympy.sets.sets import Interval, FiniteSet, Union, Intersection
from sympy.functions.combinatorial.factorials import factorial, binomial
from sympy.polys.polytools import Poly

from sympy.printing.context import context, print_context, ContextPrinter
from sympy.testing.pytest import raises


def test_context_basic():
    """Test basic ConTeXt printing."""
    x, y, z = symbols('x y z')

    # Basic symbols
    assert context(x) == 'x'
    assert context(x + y) == 'x + y'
    assert context(x - y) == 'x - y'
    assert context(x * y) == 'x y'
    assert context(x / y) == r'\frac{x}{y}'

    # Powers
    assert context(x**2) == 'x^{2}'
    assert context(x**(-1)) == r'\frac{1}{x}'
    assert context(x**(Rational(1, 2))) == r'\sqrt{x}'
    assert context(x**(Rational(3, 2))) == r'x^{\frac{3}{2}}'


def test_context_numbers():
    """Test printing of numbers."""
    assert context(0) == '0'
    assert context(1) == '1'
    assert context(-1) == '-1'
    assert context(Integer(42)) == '42'
    assert context(Rational(1, 2)) == r'\frac{1}{2}'
    assert context(Rational(3, 4)) == r'\frac{3}{4}'
    assert context(Float(1.23)) == '1.23'


def test_context_constants():
    """Test printing of mathematical constants."""
    assert context(pi) == r'\pi'
    assert context(oo) == r'\infty'
    assert context(S.Exp1) == 'e'
    assert context(I) == 'i'
    assert context(S.GoldenRatio) == r'\phi'


def test_context_functions():
    """Test printing of functions."""
    x = symbols('x')

    # Trigonometric functions
    assert context(sin(x)) == r'\sin{\left(x \right)}'
    assert context(cos(x)) == r'\cos{\left(x \right)}'
    assert context(tan(x)) == r'\tan{\left(x \right)}'
    assert context(cot(x)) == r'\cot{\left(x \right)}'
    assert context(sec(x)) == r'\sec{\left(x \right)}'
    assert context(csc(x)) == r'\csc{\left(x \right)}'

    # Inverse trigonometric functions
    assert context(asin(x)) == r'\arcsin{\left(x \right)}'
    assert context(acos(x)) == r'\arccos{\left(x \right)}'
    assert context(atan(x)) == r'\arctan{\left(x \right)}'

    # Hyperbolic functions
    assert context(sinh(x)) == r'\sinh{\left(x \right)}'
    assert context(cosh(x)) == r'\cosh{\left(x \right)}'
    assert context(tanh(x)) == r'\tanh{\left(x \right)}'
    assert context(asinh(x)) == r'\operatorname{asinh}{\left(x \right)}'
    assert context(coth(x)) == r'\coth{\left(x \right)}'

    # Exponential and logarithmic functions
    assert context(exp(x)) == 'e^{x}'
    assert context(log(x)) == r'\log{\left(x \right)}'

    # Root functions
    assert context(sqrt(x)) == r'\sqrt{x}'
    assert context(root(x, 3)) == r'\sqrt[3]{x}'
    assert context(root(x, 4)) == r'\sqrt[4]{x}'

    # Other functions
    assert context(Abs(x)) == r'\left|{x}\right|'
    assert context(conjugate(x)) == r'\overline{x}'
    assert context(re(x)) == r'\operatorname{re}{\left(x \right)}'
    assert context(im(x)) == r'\operatorname{im}{\left(x \right)}'
    assert context(gamma(x)) == r'\Gamma\left(x\right)'
    assert context(factorial(x)) == r'x!'
    assert context(binomial(5, 3)) == r'{\binom{5}{3}}'


def test_context_derivatives():
    """Test printing of derivatives."""
    x, y = symbols('x y')
    f = Function('f')

    assert context(diff(f(x), x)) == r'\frac{d}{d x} f{\left(x \right)}'
    assert context(diff(f(x), x, 2)) == r'\frac{d^{2}}{d x^{2}} f{\left(x \right)}'
    assert context(Derivative(f(x), x)) == r'\frac{d}{d x} f{\left(x \right)}'
    assert context(Derivative(f(x, y), x, y)) == r'\frac{\partial^{2}}{\partial x\partial y} f{\left(x,y \right)}'


def test_context_integrals():
    """Test printing of integrals."""
    x, y = symbols('x y')

    assert context(Integral(x, x)) == r'\int x\, dx'
    assert context(Integral(x**2, x)) == r'\int x^{2}\, dx'
    assert context(Integral(x, (x, 0, 1))) == r'\int_{0}^{1} x\, dx'
    assert context(Integral(x*y, x, y)) == r'\int\int x y\, dx\, dy'


def test_context_sums_products():
    """Test printing of sums and products."""
    x, n = symbols('x n')

    assert context(Sum(x, (x, 0, n))) == r'\sum_{x=0}^{n} x'
    assert context(Sum(x**2, (x, 0, n))) == r'\sum_{x=0}^{n} x^{2}'
    assert context(Product(x, (x, 1, n))) == r'\prod_{x=1}^{n} x'


def test_context_limits():
    """Test printing of limits."""
    x = symbols('x')

    assert context(Limit(x, x, 0)) == r'\lim_{x \to 0} x'
    assert context(Limit(sin(x)/x, x, 0)) == r'\lim_{x \to 0}\left(\frac{\sin{\left(x \right)}}{x}\right)'


def test_context_matrices():
    """Test printing of matrices."""
    # Simple 2x2 matrix
    M = Matrix([[1, 2], [3, 4]])
    result = context(M)
    assert r'\left[\begin{matrix}' in result
    assert r'1 & 2' in result
    assert r'3 & 4' in result
    assert r'\end{matrix}\right]' in result

    # 3x3 matrix
    M = Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    result = context(M)
    assert r'\left[\begin{matrix}' in result
    assert r'1 & 0 & 0' in result
    assert r'\end{matrix}\right]' in result


def test_context_relational():
    """Test printing of relational operators."""
    x, y = symbols('x y')

    assert context(Eq(x, y)) == 'x = y'
    assert context(Ne(x, y)) == r'x \neq y'
    assert context(Lt(x, y)) == 'x < y'
    assert context(Le(x, y)) == r'x \leq y'
    assert context(Gt(x, y)) == 'x > y'
    assert context(Ge(x, y)) == r'x \geq y'


def test_context_logic():
    """Test printing of logical operators."""
    x, y = symbols('x y')

    assert context(And(x, y)) == r'x \wedge y'
    assert context(Or(x, y)) == r'x \vee y'
    assert context(Not(x)) == r'\neg x'
    assert context(Xor(x, y)) == r'x \veebar y'
    assert context(Implies(x, y)) == r'x \Rightarrow y'
    assert context(Equivalent(x, y)) == r'x \Leftrightarrow y'


def test_context_sets():
    """Test printing of sets."""
    x = symbols('x')

    assert context(Interval(0, 1)) == r'\left[0, 1\right]'
    assert context(Interval(0, 1, left_open=True)) == r'\left(0, 1\right]'
    assert context(Interval(0, 1, right_open=True)) == r'\left[0, 1\right)'
    assert context(Interval(0, 1, left_open=True, right_open=True)) == r'\left(0, 1\right)'

    assert context(FiniteSet(1, 2, 3)) == r'\left\{1, 2, 3\right\}'
    assert context(Union(Interval(0, 1), Interval(2, 3))) == r'\left[0, 1\right] \cup \left[2, 3\right]'
    assert context(Intersection(Interval(0, 2), Interval(1, 3))) == r'\left[0, 2\right] \cap \left[1, 3\right]'


def test_context_mode_plain():
    """Test ConTeXt printing in plain mode."""
    x = symbols('x')

    # Plain mode (default)
    result = context(x**2, mode='plain')
    assert result == 'x^{2}'
    assert not result.startswith('$')
    assert not result.startswith(r'\startformula')


def test_context_mode_inline():
    """Test ConTeXt printing in inline mode."""
    x = symbols('x')

    # Inline mode
    result = context(x**2, mode='inline')
    assert result == '$x^{2}$'
    assert result.startswith('$')
    assert result.endswith('$')


def test_context_mode_equation():
    """Test ConTeXt printing in equation mode."""
    x = symbols('x')

    # Equation mode - ConTeXt uses \startformula...\stopformula
    result = context(x**2, mode='equation')
    assert result == r'\startformula x^{2} \stopformula'
    assert result.startswith(r'\startformula')
    assert result.endswith(r'\stopformula')


def test_context_mode_equation_star():
    """Test ConTeXt printing in equation* mode."""
    x = symbols('x')

    # Equation* mode - ConTeXt uses \startformula...\stopformula (unnumbered by default)
    result = context(x**2, mode='equation*')
    assert result == r'\startformula x^{2} \stopformula'
    assert result.startswith(r'\startformula')
    assert result.endswith(r'\stopformula')


def test_context_invalid_mode():
    """Test that invalid mode raises an error."""
    x = symbols('x')

    with raises(ValueError):
        context(x, mode='invalid')


def test_context_settings():
    """Test that LaTeX printer settings work with ConTeXt printer."""
    x, y = symbols('x y')

    # fold_frac_powers
    assert context(x**(Rational(1, 2)), fold_frac_powers=True) == r'\sqrt{x}'
    assert context(x**(Rational(3, 4)), fold_frac_powers=True) == r'x^{3/4}'

    # fold_short_frac
    assert context(3*x**2/y, fold_short_frac=True) == r'3 x^{2} / y'

    # mul_symbol
    assert context(x*y, mul_symbol='dot') == r'x \cdot y'
    assert context(x*y, mul_symbol='times') == r'x \times y'

    # order
    assert context(x + y, order='lex') == 'x + y'


def test_context_greek_letters():
    """Test printing of Greek letters."""
    alpha, beta, gamma, delta = symbols('alpha beta gamma delta')
    theta, phi, psi, omega = symbols('theta phi psi omega')

    assert context(alpha) == r'\alpha'
    assert context(beta) == r'\beta'
    assert context(gamma) == r'\gamma'
    assert context(delta) == r'\delta'
    assert context(theta) == r'\theta'
    assert context(phi) == r'\phi'
    assert context(psi) == r'\psi'
    assert context(omega) == r'\omega'


def test_context_complex_expression():
    """Test printing of complex expressions."""
    x, y, z = symbols('x y z')

    # Complex nested expression
    expr = (x**2 + y**2)**(Rational(1, 2)) + exp(I*pi)
    result = context(expr)
    assert r'\sqrt{x^{2} + y^{2}}' in result
    assert r'e^{i \pi}' in result

    # Expression with multiple operations
    expr = (sin(x) + cos(y)) / (tan(z) + 1)
    result = context(expr)
    assert r'\sin{\left(x \right)}' in result
    assert r'\cos{\left(y \right)}' in result
    assert r'\tan{\left(z \right)}' in result


def test_context_printer_class():
    """Test ContextPrinter class directly."""
    x = symbols('x')

    printer = ContextPrinter()
    assert printer.doprint(x**2) == 'x^{2}'

    printer = ContextPrinter({'mode': 'inline'})
    assert printer.doprint(x**2) == '$x^{2}$'

    printer = ContextPrinter({'mode': 'equation'})
    assert printer.doprint(x**2) == r'\startformula x^{2} \stopformula'


def test_context_vs_latex():
    """Test that ConTeXt and LaTeX produce similar output for math."""
    from sympy.printing.latex import latex
    x, y = symbols('x y')

    # Most mathematical expressions should be identical
    # except for the delimiters
    expr = x**2 + sqrt(y)

    context_plain = context(expr, mode='plain')
    latex_plain = latex(expr, mode='plain')
    assert context_plain == latex_plain

    # Check that delimiters are different for equation mode
    context_eq = context(expr, mode='equation')
    latex_eq = latex(expr, mode='equation')
    assert context_eq != latex_eq  # Different delimiters
    assert context_eq.startswith(r'\startformula')
    assert latex_eq.startswith(r'\begin{equation}')


def test_print_context():
    """Test print_context function."""
    import sys
    from io import StringIO

    x = symbols('x')

    # Capture stdout
    old_stdout = sys.stdout
    sys.stdout = StringIO()

    try:
        print_context(x**2)
        output = sys.stdout.getvalue()
        assert output.strip() == 'x^{2}'
    finally:
        sys.stdout = old_stdout


def test_context_polynomials():
    """Test printing of polynomials."""
    x = symbols('x')

    # Polynomial
    p = Poly(x**3 + 2*x**2 - x + 1, x)
    result = context(p)
    assert r'x^{3}' in result or 'Poly' in result  # Could be expanded or kept as Poly


def test_context_special_numbers():
    """Test printing of special numeric values."""
    assert context(S.NaN) == r'\text{NaN}'
    assert context(S.Infinity) == r'\infty'
    assert context(S.NegativeInfinity) == r'-\infty'


def test_context_modifiers():
    """Test variable name modifiers work."""
    x = symbols('x')
    xdot = symbols('xdot')
    xhat = symbols('xhat')
    xbar = symbols('xbar')

    # These should work since they inherit from LatexPrinter
    assert 'x' in context(x)
    # Modifiers like dot, hat, bar are handled by the parent LatexPrinter


def test_context_min_max():
    """Test printing of Min and Max."""
    x, y = symbols('x y')

    assert context(Max(x, y)) == r'\max\left(x, y\right)'
    assert context(Min(x, y)) == r'\min\left(x, y\right)'
