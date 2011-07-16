from sympy import (symbols, Rational, Symbol, Integral, log, diff, sin, exp,
        Function, factorial, floor, ceiling, Abs, re, im, conjugate, gamma,
        Order, Piecewise, Matrix, asin, Interval, EmptySet, Union, S, Sum,
        Limit, oo, Poly, Float, lowergamma, uppergamma, hyper, meijerg)
from sympy.abc import mu, tau
from sympy.printing.latex import latex
from sympy.utilities.pytest import XFAIL, raises
from sympy.functions import DiracDelta

x,y = symbols('x,y')
k,n = symbols('k,n', integer=True)

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

    assert latex(x**(Rational(1,2))) == r"\sqrt{x}"
    assert latex(x**(Rational(1,3))) == r"\sqrt[3]{x}"
    assert latex(x**(Rational(3,2))) == r"x^{\frac{3}{2}}"
    assert latex(x**(Rational(1,2)),itex=True) == r"\sqrt{x}"
    assert latex(x**(Rational(1,3)),itex=True) == r"\root{3}{x}"
    assert latex(x**(Rational(3,2)),itex=True) == r"x^{\frac{3}{2}}"
    assert latex(x**(Rational(3,4))) == r"x^{\frac{3}{4}}"
    assert latex(x**(Rational(3,4)), fold_frac_powers=True) == "x^{3/4}"

    assert latex(1.5e20*x) == r"1.5 \times 10^{20} x"
    assert latex(1.5e20*x, mul_symbol='dot') == r"1.5 \cdot 10^{20} \cdot x"

def test_latex_Float():
    assert latex(Float(1.0e100)) == r"1.0 \times 10^{100}"
    assert latex(Float(1.0e-100)) == r"1.0 \times 10^{-100}"
    latex(Float(1.0e-100), mul_symbol="dot") == r"1.0 \cdot 10^{-100}"
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

    #assert latex(volume * rho == mass) == r"\rho \mathrm{volume} = \mathrm{mass}"
    #assert latex(volume / mass * rho == 1) == r"\rho \mathrm{volume} {\mathrm{mass}}^{(-1)} = 1"
    #assert latex(mass**3 * volume**3) == r"{\mathrm{mass}}^{3} \cdot {\mathrm{volume}}^{3}"

def test_latex_functions():
    assert latex(exp(x)) == "e^{x}"
    assert latex(exp(1)+exp(2)) == "e + e^{2}"

    f = Function('f')
    assert latex(f(x)) == '\\operatorname{f}\\left(x\\right)'

    beta = Function('beta')

    assert latex(beta(x)) == r"\operatorname{beta}\left(x\right)"
    assert latex(sin(x)) == r"\operatorname{sin}\left(x\right)"
    assert latex(sin(x), fold_func_brackets=True) == r"\operatorname{sin}x"
    assert latex(sin(2*x**2), fold_func_brackets=True) == \
    r"\operatorname{sin}2 x^{2}"
    assert latex(sin(x**2), fold_func_brackets=True) == \
    r"\operatorname{sin}x^{2}"

    assert latex(asin(x)**2) == r"\operatorname{asin}^{2}\left(x\right)"
    assert latex(asin(x)**2,inv_trig_style="full") == \
        r"\operatorname{arcsin}^{2}\left(x\right)"
    assert latex(asin(x)**2,inv_trig_style="power") == \
        r"\operatorname{sin}^{-1}\left(x\right)^{2}"
    assert latex(asin(x**2),inv_trig_style="power",fold_func_brackets=True) == \
        r"\operatorname{sin}^{-1}x^{2}"

    assert latex(factorial(k)) == r"k!"
    assert latex(factorial(-k)) == r"\left(- k\right)!"

    assert latex(floor(x)) == r"\lfloor{x}\rfloor"
    assert latex(ceiling(x)) == r"\lceil{x}\rceil"
    assert latex(Abs(x)) == r"\lvert{x}\rvert"
    assert latex(re(x)) == r"\Re{x}"
    assert latex(im(x)) == r"\Im{x}"
    assert latex(conjugate(x)) == r"\overline{x}"
    assert latex(gamma(x)) == r"\operatorname{\Gamma}\left(x\right)"
    assert latex(Order(x)) == r"\operatorname{\mathcal{O}}\left(x\right)"
    assert latex(lowergamma(x, y)) == r'\operatorname{\gamma}\left(x, y\right)'
    assert latex(uppergamma(x, y)) == r'\operatorname{\Gamma}\left(x, y\right)'

def test_hyper_printing():
    from sympy import pi, Tuple
    from sympy.abc import x, z

    assert latex(meijerg(Tuple(pi, pi, x), Tuple(1), \
                         (0,1), Tuple(1, 2, 3/pi),z)) == \
             r'{G_{4, 5}^{2, 3}\left.\left(\begin{matrix} \pi, \pi, x & 1 \\0, 1 & 1, 2, \frac{3}{\pi} \end{matrix} \right| {z} \right)}'
    assert latex(meijerg(Tuple(), Tuple(1), (0,), Tuple(),z)) == \
             r'{G_{1, 1}^{1, 0}\left.\left(\begin{matrix}  & 1 \\0 &  \end{matrix} \right| {z} \right)}'
    assert latex(hyper((x, 2), (3,), z)) == \
               r'{{}_{2}F_{1}\left.\left(\begin{matrix} x, 2 ' \
               r'\\ 3 \end{matrix}\right| {z} \right)}'
    assert latex(hyper(Tuple(), Tuple(1), z)) == \
               r'{{}_{0}F_{1}\left.\left(\begin{matrix}  ' \
               r'\\ 1 \end{matrix}\right| {z} \right)}'

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

def test_latex_brackets():
    assert latex((-1)**x) == r"\left(-1\right)^{x}"

def test_latex_derivatives():
    assert latex(diff(x**3, x, evaluate=False)) == \
    r"\frac{\partial}{\partial x} x^{3}"
    assert latex(diff(sin(x)+x**2, x, evaluate=False)) == \
    r"\frac{\partial}{\partial x}\left(x^{2} + \operatorname{sin}\left(x\right)\right)"

def test_latex_integrals():
    assert latex(Integral(log(x), x)) == r"\int \operatorname{log}\left(x\right)\,dx"
    assert latex(Integral(x**2, (x,0,1))) == r"\int_{0}^{1} x^{2}\,dx"
    assert latex(Integral(x**2, (x,10,20))) == r"\int_{10}^{20} x^{2}\,dx"
    assert latex(Integral(y*x**2, (x,0,1), y)) == r"\int\int_{0}^{1} x^{2} y\,dx dy"
    assert latex(Integral(y*x**2, (x,0,1), y), mode='equation*') \
        == r"\begin{equation*}\int\int\limits_{0}^{1} x^{2} y\,dx dy\end{equation*}"
    assert latex(Integral(y*x**2, (x,0,1), y), mode='equation*', itex=True) \
        == r"$$\int\int_{0}^{1} x^{2} y\,dx dy$$"

def test_latex_intervals():
    a = Symbol('a', real=True)
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
        r"\left[3, 4\right] \cup \left\{1, 2\right\}"

def test_latex_sum():
    assert latex(Sum(x*y**2, (x, -2, 2), (y, -5, 5))) == \
        r"\sum_{\substack{-2 \leq x \leq 2\\-5 \leq y \leq 5}} x y^{2}"
    assert latex(Sum(x**2, (x, -2, 2))) == \
        r"\sum_{x=-2}^{2} x^{2}"
    assert latex(Sum(x**2 + y, (x, -2, 2))) == \
        r"\sum_{x=-2}^{2} \left(x^{2} + y\right)"

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

def test_latex_Matrix():
    M = Matrix([[1+x, y],[y, x-1]])
    assert latex(M) == '\\left(\\begin{smallmatrix}x + 1 & y\\\\y & x -'\
                       '1\\end{smallmatrix}\\right)'
    settings = {'mat_str' : 'bmatrix'}
    assert latex(M, **settings) == '\\left(\\begin{bmatrix}x + 1 & y\\\\y &'\
           ' x -1\\end{bmatrix}\\right)'
    settings['mat_delim'] = None
    assert latex(M, **settings) == '\\begin{bmatrix}x + 1 & y\\\\y & x -1'\
                       '\\end{bmatrix}'
    assert latex(M) == '\\left(\\begin{smallmatrix}x + 1 & y\\\\y & x -1'\
                       '\\end{smallmatrix}\\right)'

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
    assert latex(y) == '4 \\times 4^{\\operatorname{log}\\left(2\\right)}'
    assert latex(1/y) == '\\frac{1}{4 \\times 4^{\\operatorname{log}\\left(2\\right)}}'

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

def test_latex_order():
    expr = x**3 + x**2*y + 3*x*y**3 + y**4

    assert latex(expr, order='lex') == "x^{3} + x^{2} y + 3 x y^{3} + y^{4}"
    assert latex(expr, order='rev-lex') == "y^{4} + 3 x y^{3} + x^{2} y + x^{3}"

def test_settings():
    raises(TypeError, 'latex(x*y, method="garbage")')
