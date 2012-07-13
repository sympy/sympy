from sympy import (symbols, Rational, Symbol, Integral, log, diff, sin, exp,
    Function, factorial, factorial2, floor, ceiling, Abs, re, im, conjugate,
    Order, Piecewise, Matrix, asin, Interval, EmptySet, Union, S, Sum, Product,
    Limit, oo, Poly, Float, lowergamma, uppergamma, hyper, meijerg, polar_lift,
    Lambda, Poly, RootOf, RootSum, sqrt, Dict, catalan, Min, Max, cot, coth,
    re, im, root, arg, zeta, dirichlet_eta, binomial, RisingFactorial,
    FallingFactorial, polylog, lerchphi, Ei, expint, Si, Ci, Shi, Chi, gamma,
    Tuple, MellinTransform, InverseMellinTransform, LaplaceTransform,
    InverseLaplaceTransform, FourierTransform, InverseFourierTransform,
    SineTransform, InverseSineTransform, CosineTransform,
    InverseCosineTransform, FiniteSet, TransformationSet, Range, Subs)

from sympy.abc import mu, tau
from sympy.printing.latex import latex
from sympy.utilities.pytest import XFAIL, raises
from sympy.functions import DiracDelta
from sympy.logic import Implies

x, y, z, t = symbols('x y z t')
k, n = symbols('k n', integer=True)

def test_printmethod():
    class R(Abs):
        def _latex(self, printer):
            return "foo(%s)" % printer._print(self.args[0])
    assert latex(R(x)) == "foo(x)"
    class R(Abs):
        def _latex(self, printer):
            return "foo"
    assert latex(R(x)) == "foo"

def test_latex_basic():
    assert latex(1+x) == "x + 1"
    assert latex(x**2) == "x^{2}"
    assert latex(x**(1+x)) == "x^{x + 1}"
    assert latex(x**3+x+1+x**2) == "x^{3} + x^{2} + x + 1"

    assert latex(2*x*y) == "2 x y"
    assert latex(2*x*y, mul_symbol='dot') == r"2 \cdot x \cdot y"

    assert latex(sqrt(x)) == r"\sqrt{x}"
    assert latex(x**Rational(1,3)) == r"\sqrt[3]{x}"
    assert latex(sqrt(x)**3) == r"x^{\frac{3}{2}}"
    assert latex(sqrt(x),itex=True) == r"\sqrt{x}"
    assert latex(x**Rational(1,3),itex=True) == r"\root{3}{x}"
    assert latex(sqrt(x)**3,itex=True) == r"x^{\frac{3}{2}}"
    assert latex(x**Rational(3,4)) == r"x^{\frac{3}{4}}"
    assert latex(x**Rational(3,4), fold_frac_powers=True) == "x^{3/4}"

    assert latex(1.5e20*x) == r"1.5 \times 10^{20} x"
    assert latex(1.5e20*x, mul_symbol='dot') == r"1.5 \cdot 10^{20} \cdot x"

    assert latex(1/sin(x)) == r"\frac{1}{\sin{\left (x \right )}}"
    assert latex(sin(x)**-1) == r"\frac{1}{\sin{\left (x \right )}}"

    assert latex(~x) == r"\neg x"
    assert latex(x & y) == r"x \wedge y"
    assert latex(x & y & z) == r"x \wedge y \wedge z"
    assert latex(x | y) == r"x \vee y"
    assert latex(x | y | z) == r"x \vee y \vee z"
    assert latex((x & y) | z) == r"z \vee \left(x \wedge y\right)"
    assert latex(Implies(x,y)) == r"x \Rightarrow y"
    assert latex(~(x >> ~y)) == r"\neg (x \Rightarrow \neg y)"

    assert latex(~x, symbol_names={x: "x_i"}) == r"\neg x_i"
    assert latex(x & y, symbol_names={x: "x_i", y: "y_i"}) == \
        r"x_i \wedge y_i"
    assert latex(x & y & z, symbol_names={x: "x_i", y: "y_i", z: "z_i"}) == \
        r"x_i \wedge y_i \wedge z_i"
    assert latex(x | y, symbol_names={x: "x_i", y: "y_i"}) == r"x_i \vee y_i"
    assert latex(x | y | z, symbol_names={x: "x_i", y: "y_i", z: "z_i"}) == \
        r"x_i \vee y_i \vee z_i"
    assert latex((x & y) | z, symbol_names={x: "x_i", y: "y_i", z: "z_i"}) ==\
        r"z_i \vee \left(x_i \wedge y_i\right)"
    assert latex(Implies(x,y), symbol_names={x: "x_i", y: "y_i"}) == \
        r"x_i \Rightarrow y_i"

def test_latex_Float():
    assert latex(Float(1.0e100)) == r"1.0 \times 10^{100}"
    assert latex(Float(1.0e-100)) == r"1.0 \times 10^{-100}"
    assert latex(Float(1.0e-100), mul_symbol="dot") == r"1.0 \cdot 10^{-100}"
    assert latex(1.0*oo) == r"\infty"
    assert latex(-1.0*oo) == r"- \infty"

def test_latex_symbols():
    Gamma, lmbda, rho = map(Symbol, ('Gamma', 'lambda', 'rho'))
    mass, volume = map(Symbol, ('mass', 'volume'))
    assert latex(Gamma + lmbda) == r"\Gamma + \lambda"
    assert latex(Gamma * lmbda) == r"\Gamma \lambda"
    assert latex(Symbol('q21')) == r"q_{21}"
    assert latex(Symbol('epsilon0')) == r"\epsilon_{0}"
    assert latex(Symbol('91')) == r"91"
    assert latex(Symbol('alpha_new')) == r"\alpha_{new}"
    assert latex(Symbol('C^orig')) == r"C^{orig}"

@XFAIL
def test_latex_symbols_failing():
    assert latex(volume * rho == mass) == r"\rho \mathrm{volume} = \mathrm{mass}"
    assert latex(volume / mass * rho == 1) == r"\rho \mathrm{volume} {\mathrm{mass}}^{(-1)} = 1"
    assert latex(mass**3 * volume**3) == r"{\mathrm{mass}}^{3} \cdot {\mathrm{volume}}^{3}"

def test_latex_functions():
    assert latex(exp(x)) == "e^{x}"
    assert latex(exp(1) + exp(2)) == "e + e^{2}"

    f = Function('f')
    assert latex(f(x)) == '\\operatorname{f}{\\left (x \\right )}'

    beta = Function('beta')

    assert latex(beta(x)) == r"\beta{\left (x \right )}"
    assert latex(sin(x)) == r"\sin{\left (x \right )}"
    assert latex(sin(x), fold_func_brackets=True) == r"\sin {x}"
    assert latex(sin(2*x**2), fold_func_brackets=True) == \
        r"\sin {2 x^{2}}"
    assert latex(sin(x**2), fold_func_brackets=True) == \
        r"\sin {x^{2}}"

    assert latex(asin(x)**2) == r"\operatorname{asin}^{2}{\left (x \right )}"
    assert latex(asin(x)**2, inv_trig_style="full") == \
        r"\arcsin^{2}{\left (x \right )}"
    assert latex(asin(x)**2, inv_trig_style="power") == \
        r"\sin^{-1}{\left (x \right )}^{2}"
    assert latex(asin(x**2), inv_trig_style="power",
                 fold_func_brackets=True) == \
        r"\sin^{-1} {x^{2}}"

    assert latex(factorial(k)) == r"k!"
    assert latex(factorial(-k)) == r"\left(- k\right)!"

    assert latex(factorial2(k)) == r"k!!"
    assert latex(factorial2(-k)) == r"\left(- k\right)!!"

    assert latex(binomial(2, k)) == r"{\binom{2}{k}}"

    assert latex(FallingFactorial(3, k)) == r"{\left(3\right)}_{\left(k\right)}"
    assert latex(RisingFactorial(3, k)) == r"{\left(3\right)}^{\left(k\right)}"

    assert latex(floor(x)) == r"\lfloor{x}\rfloor"
    assert latex(ceiling(x)) == r"\lceil{x}\rceil"
    assert latex(Min(x, 2, x**3)) == r"\min\left(2, x, x^{3}\right)"
    assert latex(Max(x, 2, x**3)) == r"\max\left(2, x, x^{3}\right)"
    assert latex(Abs(x)) == r"\lvert{x}\rvert"
    assert latex(re(x)) == r"\Re{x}"
    assert latex(re(x + y)) == r"\Re{x} + \Re{y}"
    assert latex(im(x)) == r"\Im{x}"
    assert latex(conjugate(x)) == r"\overline{x}"
    assert latex(gamma(x)) == r"\Gamma\left(x\right)"
    assert latex(Order(x)) == r"\mathcal{O}\left(x\right)"
    assert latex(lowergamma(x, y)) == r'\gamma\left(x, y\right)'
    assert latex(uppergamma(x, y)) == r'\Gamma\left(x, y\right)'

    assert latex(cot(x)) == r'\cot{\left (x \right )}'
    assert latex(coth(x)) == r'\coth{\left (x \right )}'
    assert latex(re(x)) == r'\Re{x}'
    assert latex(im(x)) == r'\Im{x}'
    assert latex(root(x, y)) == r'x^{\frac{1}{y}}'
    assert latex(arg(x)) == r'\arg{\left (x \right )}'
    assert latex(zeta(x)) == r'\zeta\left(x\right)'

    assert latex(zeta(x)) == r"\zeta\left(x\right)"
    assert latex(zeta(x)**2) == r"\zeta^{2}\left(x\right)"
    assert latex(zeta(x, y)) == r"\zeta\left(x, y\right)"
    assert latex(zeta(x, y)**2) == r"\zeta^{2}\left(x, y\right)"
    assert latex(dirichlet_eta(x)) == r"\eta\left(x\right)"
    assert latex(dirichlet_eta(x)**2) == r"\eta^{2}\left(x\right)"
    assert latex(polylog(x, y)) == r"\operatorname{Li}_{x}\left(y\right)"
    assert latex(polylog(x, y)**2) == r"\operatorname{Li}_{x}^{2}\left(y\right)"
    assert latex(lerchphi(x, y, n)) == r"\Phi\left(x, y, n\right)"
    assert latex(lerchphi(x, y, n)**2) == r"\Phi^{2}\left(x, y, n\right)"

    assert latex(Ei(x)) == r'\operatorname{Ei}{\left (x \right )}'
    assert latex(Ei(x)**2) == r'\operatorname{Ei}^{2}{\left (x \right )}'
    assert latex(expint(x, y)**2) == r'\operatorname{E}_{x}^{2}\left(y\right)'
    assert latex(Shi(x)**2) == r'\operatorname{Shi}^{2}{\left (x \right )}'
    assert latex(Si(x)**2) == r'\operatorname{Si}^{2}{\left (x \right )}'
    assert latex(Ci(x)**2) == r'\operatorname{Ci}^{2}{\left (x \right )}'
    assert latex(Chi(x)**2) == r'\operatorname{Chi}^{2}{\left (x \right )}'

    # Test latex printing of function names with "_"
    assert latex(polar_lift(0)) == r"\operatorname{polar\_lift}{\left (0 \right )}"
    assert latex(polar_lift(0)**3) == r"\operatorname{polar\_lift}^{3}{\left (0 \right )}"


def test_hyper_printing():
    from sympy import pi, Tuple
    from sympy.abc import x, z

    assert latex(meijerg(Tuple(pi, pi, x), Tuple(1), \
                         (0,1), Tuple(1, 2, 3/pi),z)) == \
             r'{G_{4, 5}^{2, 3}\left(\begin{matrix} \pi, \pi, x & 1 \\0, 1 & 1, 2, \frac{3}{\pi} \end{matrix} \middle| {z} \right)}'
    assert latex(meijerg(Tuple(), Tuple(1), (0,), Tuple(),z)) == \
             r'{G_{1, 1}^{1, 0}\left(\begin{matrix}  & 1 \\0 &  \end{matrix} \middle| {z} \right)}'
    assert latex(hyper((x, 2), (3,), z)) == \
               r'{{}_{2}F_{1}\left(\begin{matrix} x, 2 ' \
               r'\\ 3 \end{matrix}\middle| {z} \right)}'
    assert latex(hyper(Tuple(), Tuple(1), z)) == \
               r'{{}_{0}F_{1}\left(\begin{matrix}  ' \
               r'\\ 1 \end{matrix}\middle| {z} \right)}'

def test_latex_bessel():
    from sympy.functions.special.bessel import (besselj, bessely, besseli,
            besselk, hankel1, hankel2, jn, yn)
    from sympy.abc import z
    assert latex(besselj(n, z**2)**k) == r'J^{k}_{n}\left(z^{2}\right)'
    assert latex(bessely(n, z)) == r'Y_{n}\left(z\right)'
    assert latex(besseli(n, z)) == r'I_{n}\left(z\right)'
    assert latex(besselk(n, z)) == r'K_{n}\left(z\right)'
    assert latex(hankel1(n, z**2)**2) == \
              r'\left(H^{(1)}_{n}\left(z^{2}\right)\right)^{2}'
    assert latex(hankel2(n, z)) == r'H^{(2)}_{n}\left(z\right)'
    assert latex(jn(n, z)) == r'j_{n}\left(z\right)'
    assert latex(yn(n, z)) == r'y_{n}\left(z\right)'

def test_latex_fresnel():
    from sympy.functions.special.error_functions import (fresnels, fresnelc)
    from sympy.abc import z
    assert latex(fresnels(z)) == r'S\left(z\right)'
    assert latex(fresnelc(z)) == r'C\left(z\right)'
    assert latex(fresnels(z)**2) == r'S^{2}\left(z\right)'
    assert latex(fresnelc(z)**2) == r'C^{2}\left(z\right)'

def test_latex_brackets():
    assert latex((-1)**x) == r"\left(-1\right)^{x}"

def test_latex_derivatives():
    assert latex(diff(x**3, x, evaluate=False)) == \
    r"\frac{\partial}{\partial x} x^{3}"
    assert latex(diff(sin(x)+x**2, x, evaluate=False)) == \
    r"\frac{\partial}{\partial x}\left(x^{2} + \sin{\left (x \right )}\right)"

def test_latex_subs():
    assert latex(Subs(x*y, (x, y), (1, 2))) == r'\left. x y \right|_{\substack{ x=1\\ y=2 }}'

def test_latex_integrals():
    assert latex(Integral(log(x), x)) == r"\int \log{\left (x \right )}\, dx"
    assert latex(Integral(x**2, (x,0,1))) == r"\int_{0}^{1} x^{2}\, dx"
    assert latex(Integral(x**2, (x,10,20))) == r"\int_{10}^{20} x^{2}\, dx"
    assert latex(Integral(y*x**2, (x,0,1), y)) == r"\int\int_{0}^{1} x^{2} y\, dx\, dy"
    assert latex(Integral(y*x**2, (x,0,1), y), mode='equation*') \
        == r"\begin{equation*}\int\int\limits_{0}^{1} x^{2} y\, dx\, dy\end{equation*}"
    assert latex(Integral(y*x**2, (x,0,1), y), mode='equation*', itex=True) \
        == r"$$\int\int_{0}^{1} x^{2} y\, dx\, dy$$"
    assert latex(Integral(x, (x, 0))) == r"\int^{0} x\, dx"
    assert latex(Integral(x*y, x, y)) == r"\iint x y\, dx\, dy"
    assert latex(Integral(x*y*z, x, y, z)) == r"\iiint x y z\, dx\, dy\, dz"
    assert latex(Integral(x*y*z*t, x, y, z, t)) == \
        r"\iiiint t x y z\, dx\, dy\, dz\, dt"
    assert latex(Integral(x, x, x, x, x, x, x)) == \
        r"\int\int\int\int\int\int x\, dx\, dx\, dx\, dx\, dx\, dx"
    assert latex(Integral(x, x, y, (z, 0, 1))) == \
        r"\int_{0}^{1}\int\int x\, dx\, dy\, dz"

def test_latex_sets():
    for s in (FiniteSet, frozenset, set):
        assert latex(s([x*y, x**2])) == r"\left\{x^{2}, x y\right\}"
        assert latex(s(range(1, 6))) == r"\left\{1, 2, 3, 4, 5\right\}"
        assert latex(s(range(1, 13))) == \
            r"\left\{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12\right\}"

def test_latex_Range():
    assert latex(Range(1, 51)) ==\
            r'\left\{1, 2, \ldots, 50\right\}'
    assert latex(Range(1, 4)) == r'\left\{1, 2, 3\right\}'

def test_latex_intervals():
    a = Symbol('a', real=True)
    assert latex(Interval(0, 0)) == r"\left\{0\right\}"
    assert latex(Interval(0, a)) == r"\left[0, a\right]"
    assert latex(Interval(0, a, False, False)) == r"\left[0, a\right]"
    assert latex(Interval(0, a, True, False)) == r"\left(0, a\right]"
    assert latex(Interval(0, a, False, True)) == r"\left[0, a\right)"
    assert latex(Interval(0, a, True, True)) == r"\left(0, a\right)"

def test_latex_emptyset():
    assert latex(S.EmptySet) == r"\emptyset"

def test_latex_union():
    assert latex(Union(Interval(0, 1), Interval(2, 3))) == \
        r"\left[0, 1\right] \cup \left[2, 3\right]"
    assert latex(Union(Interval(1, 1), Interval(2, 2), Interval(3, 4))) == \
        r"\left\{1, 2\right\} \cup \left[3, 4\right]"

def test_latex_productset():
    line = Interval(0,1)
    bigline = Interval(0, 10)
    fset = FiniteSet(1, 2, 3)
    assert latex(line**2) == r"%s^2"%latex(line)
    assert latex(line * bigline * fset) == r"%s \times %s \times %s"%(
            latex(line), latex(bigline), latex(fset))

def test_latex_Naturals():
    assert latex(S.Naturals) == r"\mathbb{N}"
    assert latex(S.Integers) == r"\mathbb{Z}"

def test_latex_TransformationSet():
    x = Symbol('x')
    assert latex(TransformationSet(Lambda(x, x**2), S.Naturals)) == \
            r"\left\{x^{2}\; |\; x \in \mathbb{N}\right\}"

def test_latex_sum():
    assert latex(Sum(x*y**2, (x, -2, 2), (y, -5, 5))) == \
        r"\sum_{\substack{-2 \leq x \leq 2\\-5 \leq y \leq 5}} x y^{2}"
    assert latex(Sum(x**2, (x, -2, 2))) == \
        r"\sum_{x=-2}^{2} x^{2}"
    assert latex(Sum(x**2 + y, (x, -2, 2))) == \
        r"\sum_{x=-2}^{2} \left(x^{2} + y\right)"

def test_latex_product():
    assert latex(Product(x*y**2, (x, -2, 2), (y, -5, 5))) == \
        r"\prod_{\substack{-2 \leq x \leq 2\\-5 \leq y \leq 5}} x y^{2}"
    assert latex(Product(x**2, (x, -2, 2))) == \
        r"\prod_{x=-2}^{2} x^{2}"
    assert latex(Product(x**2 + y, (x, -2, 2))) == \
        r"\prod_{x=-2}^{2} \left(x^{2} + y\right)"

def test_latex_limits():
    assert latex(Limit(x, x, oo)) == r"\lim_{x \to \infty} x"

def test_issue469():
    beta = Symbol(r'\beta')
    y = beta+x
    assert latex(y) in [r'\beta + x', r'x + \beta']

    beta = Symbol(r'beta')
    y = beta+x
    assert latex(y) in [r'\beta + x', r'x + \beta']

def test_latex():
    assert latex((2*tau)**Rational(7,2)) == "8 \\sqrt{2} \\tau^{\\frac{7}{2}}"
    assert latex((2*mu)**Rational(7,2), mode='equation*') == \
            "\\begin{equation*}8 \\sqrt{2} \\mu^{\\frac{7}{2}}\\end{equation*}"
    assert latex((2*mu)**Rational(7,2), mode='equation', itex=True) == \
            "$$8 \\sqrt{2} \\mu^{\\frac{7}{2}}$$"
    assert latex([2/x, y]) =="\\begin{bmatrix}\\frac{2}{x}, & y\\end{bmatrix}"


def test_latex_dict():
    d = {Rational(1): 1, x**2: 2, x: 3, x**3: 4}
    assert latex(d) == '\\begin{Bmatrix}1 : 1, & x : 3, & x^{2} : 2, & x^{3} : 4\\end{Bmatrix}'
    D = Dict(d)
    assert latex(D) == '\\begin{Bmatrix}1 : 1, & x : 3, & x^{2} : 2, & x^{3} : 4\\end{Bmatrix}'

def test_latex_rational():
    #tests issue 874
    assert latex(-Rational(1,2)) == "- \\frac{1}{2}"
    assert latex(Rational(-1,2)) == "- \\frac{1}{2}"
    assert latex(Rational(1,-2)) == "- \\frac{1}{2}"
    assert latex(-Rational(-1,2)) == "\\frac{1}{2}"
    assert latex(-Rational(1,2)*x) == "- \\frac{1}{2} x"
    assert latex(-Rational(1,2)*x+Rational(-2,3)*y) in [
            "- \\frac{1}{2} x - \\frac{2}{3} y",
            "- \\frac{2}{3} y - \\frac{1}{2} x",
            ]

def test_latex_inverse():
    #tests issue 1030
    assert latex(1/x) == "\\frac{1}{x}"
    assert latex(1/(x+y)) in ["\\frac{1}{x + y}", "\\frac{1}{y + x}"]

def test_latex_DiracDelta():
    assert latex(DiracDelta(x)) == "\\delta\\left(x\\right)"
    assert latex(DiracDelta(x,0)) == "\\delta\\left(x\\right)"
    assert latex(DiracDelta(x,5)) == "\\delta^{\\left( 5 \\right)}\\left( x \\right)"

def test_mode():
    expr = x+y
    assert latex(expr) == 'x + y'
    assert latex(expr, mode='plain') == 'x + y'
    assert latex(expr, mode='inline') == '$x + y$'
    assert latex(expr, mode='equation*')== '\\begin{equation*}x + y\\end{equation*}'
    assert latex(expr, mode='equation')== '\\begin{equation}x + y\\end{equation}'

def test_latex_Piecewise():
    p = Piecewise((x,x<1),(x**2,True))
    assert latex(p) == "\\begin{cases} x & \\text{for}\: x < 1 \\\\x^{2} &" \
                       " \\text{otherwise} \\end{cases}"
    assert latex(p, itex=True) == "\\begin{cases} x & \\text{for}\: x \\lt 1 \\\\x^{2} &" \
                                  " \\text{otherwise} \\end{cases}"
    p = Piecewise((x, x < 0), (0, x >= 0))
    assert latex(p) == "\\begin{cases} x & \\text{for}\\: x < 0 \\\\0 &" \
                       " \\text{for}\\: x \\geq 0 \\end{cases}"

def test_latex_Matrix():
    M = Matrix([[1+x, y],[y, x-1]])
    assert latex(M) == '\\left[\\begin{smallmatrix}x + 1 & y\\\\y & x -'\
                       '1\\end{smallmatrix}\\right]'
    settings = {'mat_str' : 'bmatrix'}
    assert latex(M, **settings) == '\\left[\\begin{bmatrix}x + 1 & y\\\\y &'\
           ' x -1\\end{bmatrix}\\right]'
    settings['mat_delim'] = None
    assert latex(M, **settings) == '\\begin{bmatrix}x + 1 & y\\\\y & x -1'\
                       '\\end{bmatrix}'
    assert latex(M) == '\\left[\\begin{smallmatrix}x + 1 & y\\\\y & x -1'\
                       '\\end{smallmatrix}\\right]'

def test_latex_mul_symbol():
    assert latex(4*4**x, mul_symbol='times') == "4 \\times 4^{x}"
    assert latex(4*4**x, mul_symbol='dot') == "4 \\cdot 4^{x}"
    assert latex(4*4**x, mul_symbol='ldot') == "4 \,.\, 4^{x}"

    assert latex(4*x, mul_symbol='times') == "4 \\times x"
    assert latex(4*x, mul_symbol='dot') == "4 \\cdot x"
    assert latex(4*x, mul_symbol='ldot') == "4 \,.\, x"

def test_latex_Poly():
    assert latex(Poly(x**2 + 2 * x, x)) == r"x^{2} + 2 x"

def test_latex_issue1282():
    y = 4*4**log(2)
    assert latex(y) == '4 \\times 4^{\\log{\\left (2 \\right )}}'
    assert latex(1/y) == '\\frac{1}{4 \\times 4^{\\log{\\left (2 \\right )}}}'

def test_latex_issue1477():
    assert latex(Symbol("beta_13_2")) == r"\beta_{13 2}"
    assert latex(Symbol("beta_132_20")) == r"\beta_{132 20}"
    assert latex(Symbol("beta_13")) == r"\beta_{13}"
    assert latex(Symbol("x_a_b")) == r"x_{a b}"
    assert latex(Symbol("x_1_2_3")) == r"x_{1 2 3}"
    assert latex(Symbol("x_a_b1")) == r"x_{a b1}"
    assert latex(Symbol("x_a_1")) == r"x_{a 1}"
    assert latex(Symbol("x_1_a")) == r"x_{1 a}"
    assert latex(Symbol("x_1^aa")) == r"x^{aa}_{1}"
    assert latex(Symbol("x_1__aa")) == r"x^{aa}_{1}"
    assert latex(Symbol("x_11^a")) == r"x^{a}_{11}"
    assert latex(Symbol("x_11__a")) == r"x^{a}_{11}"
    assert latex(Symbol("x_a_a_a_a")) == r"x_{a a a a}"
    assert latex(Symbol("x_a_a^a^a")) == r"x^{a a}_{a a}"
    assert latex(Symbol("x_a_a__a__a")) == r"x^{a a}_{a a}"
    assert latex(Symbol("alpha_11")) == r"\alpha_{11}"
    assert latex(Symbol("alpha_11_11")) == r"\alpha_{11 11}"
    assert latex(Symbol("alpha_alpha")) == r"\alpha_{\alpha}"
    assert latex(Symbol("alpha^aleph")) == r"\alpha^{\aleph}"
    assert latex(Symbol("alpha__aleph")) == r"\alpha^{\aleph}"

def test_latex_pow_fraction():
    x = Symbol('x')
    # Testing exp
    assert 'e^{-x}' in latex(exp(-x)/2).replace(' ', '') # Remove Whitespace

    # Testing just e^{-x} in case future changes alter behavior of muls or fracs
    # In particular current output is \frac{1}{2}e^{- x} but perhaps this will
    # change to \frac{e^{-x}}{2}

    # Testing general, non-exp, power
    assert '3^{-x}' in latex(3**-x/2).replace(' ', '')

def test_noncommutative():
    A, B, C = symbols('A,B,C', commutative=False)

    assert latex(A*B*C**-1) == "A B C^{-1}"
    assert latex(C**-1*A*B) == "C^{-1} A B"
    assert latex(A*C**-1*B) == "A C^{-1} B"

def test_latex_order():
    expr = x**3 + x**2*y + 3*x*y**3 + y**4

    assert latex(expr, order='lex') == "x^{3} + x^{2} y + 3 x y^{3} + y^{4}"
    assert latex(expr, order='rev-lex') == "y^{4} + 3 x y^{3} + x^{2} y + x^{3}"

def test_latex_Lambda():
    assert latex(Lambda(x, x + 1)) == \
        r"\Lambda {\left (x, x + 1 \right )}"
    assert latex(Lambda((x, y), x + 1)) == \
        r"\Lambda {\left (\begin{pmatrix}x, & y\end{pmatrix}, x + 1 \right )}"

def test_latex_Poly():
    assert latex(Poly(x/y, x)) == \
        r"\operatorname{Poly}{\left( \frac{x}{y}, x, domain=\mathbb{Z}\left(y\right) \right)}"
    assert latex(Poly(2.0*x + y)) == \
        r"\operatorname{Poly}{\left( 2.0 x + 1.0 y, x, y, domain=\mathbb{R} \right)}"

def test_latex_RootOf():
    assert latex(RootOf(x**5 + x + 3, 0)) == \
        r"\operatorname{RootOf} {\left(x^{5} + x + 3, 0\right)}"

def test_latex_RootSum():
    assert latex(RootSum(x**5 + x + 3, sin)) == \
        r"\operatorname{RootSum} {\left(x^{5} + x + 3, \Lambda {\left (x, \sin{\left (x \right )} \right )}\right)}"

def test_settings():
    raises(TypeError, lambda: latex(x*y, method="garbage"))

def test_latex_numbers():
    assert latex(catalan(n)) == r"C_{n}"

def test_lamda():
    assert latex(Symbol('lamda')) == r"\lambda"
    assert latex(Symbol('Lamda')) == r"\Lambda"

def test_custom_symbol_names():
    x = Symbol('x')
    y = Symbol('y')
    assert latex(x) == "x"
    assert latex(x, symbol_names={x:"x_i"}) == "x_i"
    assert latex(x + y, symbol_names={x:"x_i"}) == "x_i + y"
    assert latex(x**2, symbol_names={x:"x_i"}) == "x_i^{2}"
    assert latex(x + y, symbol_names={x:"x_i", y:"y_j"}) == "x_i + y_j"

def test_matAdd():
    from sympy import MatrixSymbol
    from sympy.printing.latex import LatexPrinter
    C = MatrixSymbol('C', 5, 5)
    B = MatrixSymbol('B', 5, 5)
    l = LatexPrinter()
    assert l._print_MatAdd(C - 2*B) in ['- 2 B + C', 'C - 2 B']
    assert l._print_MatAdd(C + 2*B) in ['2 B + C', 'C + 2 B']
    assert l._print_MatAdd(B - 2*C) in ['B - 2 C', '- 2 C + B']
    assert l._print_MatAdd(B + 2*C) in ['B + 2 C', '2 C + B']

def test_matMul():
    from sympy import MatrixSymbol
    from sympy.printing.latex import LatexPrinter
    A = MatrixSymbol('A', 5, 5)
    B = MatrixSymbol('B', 5, 5)
    x = Symbol('x')
    l = LatexPrinter()
    assert l._print_MatMul(2*A) == '2 A'
    assert l._print_MatMul(2*x*A) == '2 x A'
    assert l._print_MatMul(-2*A) == '- 2 A'
    assert l._print_MatMul(1.0*A) == '1.0 A'
    assert l._print_MatMul(sqrt(2)*A) == r'\sqrt{2} A'
    assert l._print_MatMul(-sqrt(2)*A) == r'- \sqrt{2} A'
    assert l._print_MatMul(2*sqrt(2)*x*A) == r'2 \sqrt{2} x A'
    assert l._print_MatMul(-2*A*(A+2*B)) in [r'- 2 A \left(A + 2 B\right)',
        r'- 2 A \left(2 B + A\right)']

def test_latex_RandomDomain():
    from sympy.stats import Normal, Die, Exponential, pspace, where
    X = Normal('x1', 0, 1)
    assert latex(where(X>0)) == "Domain: 0 < x_{1}"

    D = Die('d1', 6)
    assert latex(where(D>4)) == r"Domain: d_{1} = 5 \vee d_{1} = 6"

    A = Exponential('a', 1)
    B = Exponential('b', 1)
    assert latex(pspace(Tuple(A,B)).domain) =="Domain: 0 \leq a \wedge 0 \leq b"

def test_PrettyPoly():
    from sympy.polys.domains import QQ
    F = QQ.frac_field(x, y)
    R = QQ[x, y]

    assert latex(F.convert(x/(x + y))) == latex(x/(x + y))
    assert latex(R.convert(x + y)) == latex(x + y)

def test_integral_transforms():
    x = Symbol("x")
    k = Symbol("k")
    f = Function("f")
    a = Symbol("a")
    b = Symbol("b")

    assert latex(MellinTransform(f(x), x, k)) == r"\mathcal{M}_{x}\left[\operatorname{f}{\left (x \right )}\right]\left(k\right)"
    assert latex(InverseMellinTransform(f(k), k, x, a,b)) == r"\mathcal{M}^{-1}_{k}\left[\operatorname{f}{\left (k \right )}\right]\left(x\right)"

    assert latex(LaplaceTransform(f(x), x, k)) == r"\mathcal{L}_{x}\left[\operatorname{f}{\left (x \right )}\right]\left(k\right)"
    assert latex(InverseLaplaceTransform(f(k), k, x, (a,b))) == r"\mathcal{L}^{-1}_{k}\left[\operatorname{f}{\left (k \right )}\right]\left(x\right)"

    assert latex(FourierTransform(f(x), x, k)) == r"\mathcal{F}_{x}\left[\operatorname{f}{\left (x \right )}\right]\left(k\right)"
    assert latex(InverseFourierTransform(f(k), k, x)) == r"\mathcal{F}^{-1}_{k}\left[\operatorname{f}{\left (k \right )}\right]\left(x\right)"

    assert latex(CosineTransform(f(x), x, k)) == r"\mathcal{COS}_{x}\left[\operatorname{f}{\left (x \right )}\right]\left(k\right)"
    assert latex(InverseCosineTransform(f(k), k, x)) == r"\mathcal{COS}^{-1}_{k}\left[\operatorname{f}{\left (k \right )}\right]\left(x\right)"

    assert latex(SineTransform(f(x), x, k)) == r"\mathcal{SIN}_{x}\left[\operatorname{f}{\left (x \right )}\right]\left(k\right)"
    assert latex(InverseSineTransform(f(k), k, x)) == r"\mathcal{SIN}^{-1}_{k}\left[\operatorname{f}{\left (k \right )}\right]\left(x\right)"

def test_PolynomialRing():
    from sympy.polys.domains import QQ
    assert latex(QQ[x, y]) == r"\mathbb{Q}\left[x, y\right]"
    assert latex(QQ.poly_ring(x, y, order="ilex")) == \
            r"S_<^{-1}\mathbb{Q}\left[x, y\right]"

def test_categories():
    from sympy.categories import (Object, Morphism, IdentityMorphism,
                                  NamedMorphism, CompositeMorphism,
                                  Category, Diagram)

    A1 = Object("A1")
    A2 = Object("A2")
    A3 = Object("A3")

    f1 = NamedMorphism(A1, A2, "f1")
    f2 = NamedMorphism(A2, A3, "f2")
    id_A1 = IdentityMorphism(A1)

    K1 = Category("K1")

    assert latex(A1) == "A_{1}"
    assert latex(f1) == "f_{1}:A_{1}\\rightarrow A_{2}"
    assert latex(id_A1) == "id:A_{1}\\rightarrow A_{1}"
    assert latex(f2*f1) == "f_{2}\\circ f_{1}:A_{1}\\rightarrow A_{3}"

    assert latex(K1) == "\mathbf{K_{1}}"

    d = Diagram()
    assert latex(d) == "\emptyset"

    d = Diagram({f1:"unique", f2:S.EmptySet})
    assert latex(d) == "\\begin{Bmatrix}f_{2}\\circ f_{1}:A_{1}" \
           "\\rightarrow A_{3} : \\emptyset, & id:A_{1}\\rightarrow " \
           "A_{1} : \\emptyset, & id:A_{2}\\rightarrow A_{2} : " \
           "\\emptyset, & id:A_{3}\\rightarrow A_{3} : \\emptyset, " \
           "& f_{1}:A_{1}\\rightarrow A_{2} : \\left\\{unique\\right\\}, " \
           "& f_{2}:A_{2}\\rightarrow A_{3} : \\emptyset\\end{Bmatrix}"

    d = Diagram({f1:"unique", f2:S.EmptySet}, {f2 * f1: "unique"})
    assert latex(d) == "\\begin{Bmatrix}f_{2}\\circ f_{1}:A_{1}" \
           "\\rightarrow A_{3} : \\emptyset, & id:A_{1}\\rightarrow " \
           "A_{1} : \\emptyset, & id:A_{2}\\rightarrow A_{2} : " \
           "\\emptyset, & id:A_{3}\\rightarrow A_{3} : \\emptyset, " \
           "& f_{1}:A_{1}\\rightarrow A_{2} : \\left\\{unique\\right\\}," \
           " & f_{2}:A_{2}\\rightarrow A_{3} : \\emptyset\\end{Bmatrix}" \
           "\\Longrightarrow \\begin{Bmatrix}f_{2}\\circ f_{1}:A_{1}" \
           "\\rightarrow A_{3} : \\left\\{unique\\right\\}\\end{Bmatrix}"

def test_Modules():
    from sympy.polys.domains import QQ
    from sympy import homomorphism
    R = QQ[x, y]
    F = R.free_module(2)
    M = F.submodule([x, y], [1, x**2])

    assert latex(F) == r"{\mathbb{Q}\left[x, y\right]}^{2}"
    assert latex(M) == \
        r"\left< {\left[ {x},{y} \right]},{\left[ {1},{x^{2}} \right]} \right>"

    I = R.ideal(x**2, y)
    assert latex(I) == r"\left< {x^{2}},{y} \right>"

    Q = F / M
    assert latex(Q) == r"\frac{{\mathbb{Q}\left[x, y\right]}^{2}}{\left< {\left[ {x},{y} \right]},{\left[ {1},{x^{2}} \right]} \right>}"
    assert latex(Q.submodule([1, x**3/2], [2, y])) == \
        r"\left< {{\left[ {1},{\frac{1}{2} x^{3}} \right]} + {\left< {\left[ {x},{y} \right]},{\left[ {1},{x^{2}} \right]} \right>}},{{\left[ {2},{y} \right]} + {\left< {\left[ {x},{y} \right]},{\left[ {1},{x^{2}} \right]} \right>}} \right>"

    h = homomorphism(QQ[x].free_module(2), QQ[x].free_module(2), [0, 0])

    assert latex(h) == r"{\left[\begin{smallmatrix}0 & 0\\0 & 0\end{smallmatrix}\right]} : {{\mathbb{Q}\left[x\right]}^{2}} \to {{\mathbb{Q}\left[x\right]}^{2}}"

def test_QuotientRing():
    from sympy.polys.domains import QQ
    R = QQ[x]/[x**2 + 1]

    assert latex(R) == r"\frac{\mathbb{Q}\left[x\right]}{\left< {x^{2} + 1} \right>}"
    assert latex(R.one) == r"{1} + {\left< {x^{2} + 1} \right>}"
