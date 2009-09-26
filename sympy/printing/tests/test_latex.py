from sympy import symbols, Rational, Symbol, Integral, log, diff, sin, exp, \
        Function, factorial, floor, ceiling, abs, re, im, conjugate, gamma, \
        Order, Piecewise, Matrix, asin, Interval, EmptySet, Union, S
from sympy.abc import mu, tau
from sympy.printing.latex import latex
from sympy.utilities.pytest import XFAIL
from sympy.functions import DiracDelta

x,y = symbols('xy')
k,n = symbols('kn', integer=True)

def test_printmethod():
    class R(abs):
        def _latex_(self, printer):
            return "foo(%s)" % printer._print(self.args[0])
    assert latex(R(x)) == "$foo(x)$"
    class R(abs):
        def _latex_(self, printer):
            return "foo"
    assert latex(R(x)) == "$foo$"

def test_latex_basic():
    assert latex(1+x) == "$1 + x$"
    assert latex(x**2) == "$x^{2}$"
    assert latex(x**(1+x)) == "$x^{1 + x}$"
    assert latex(x**3+x+1+x**2) == "$1 + x + x^{2} + x^{3}$"

    assert latex(2*x*y) == "$2 x y$"
    assert latex(2*x*y, mul_symbol='dot') == r"$2 \cdot x \cdot y$"

    assert latex(x**(Rational(1,2))) == r"$\sqrt{x}$"
    assert latex(x**(Rational(3,2))) == r"$\sqrt[3]{x}$"
    assert latex(x**(Rational(3,4))) == r"$x^{\frac{3}{4}}$"
    assert latex(x**(Rational(3,4)), fold_frac_powers=True) == "$x^{3/4}$"

    assert latex(1.5e20*x) == r"$1.5 \times 10^{20} x$"
    assert latex(1.5e20*x, mul_symbol='dot') == r"$1.5 \cdot 10^{20} \cdot x$"


def test_latex_symbols():
    Gamma, lmbda, rho = map(Symbol, ('Gamma', 'lambda', 'rho'))
    mass, volume = map(Symbol, ('mass', 'volume'))
    assert latex(Gamma + lmbda) == r"$\Gamma + \lambda$"
    assert latex(Gamma * lmbda) == r"$\Gamma \lambda$"
    assert latex(Symbol('q21')) == r"$q_{21}$"
    assert latex(Symbol('epsilon0')) == r"$\epsilon_{0}$"
    assert latex(Symbol('91')) == r"$91$"
    assert latex(Symbol('alpha_new')) == r"$\alpha_{new}$"
    assert latex(Symbol('C^orig')) == r"$C^{orig}$"

    #assert latex(volume * rho == mass) == r"$\rho \mathrm{volume} = \mathrm{mass}$"
    #assert latex(volume / mass * rho == 1) == r"$\rho \mathrm{volume} {\mathrm{mass}}^{(-1)} = 1$"
    #assert latex(mass**3 * volume**3) == r"${\mathrm{mass}}^{3} \cdot {\mathrm{volume}}^{3}$"

def test_latex_functions():
    assert latex(exp(x)) == "$e^{x}$"
    assert latex(exp(1)+exp(2)) == "$e + e^{2}$"

    f = Function('f')
    assert latex(f(x)) == '$\\operatorname{f}\\left(x\\right)$'

    beta = Function('beta')

    assert latex(beta(x)) == r"$\operatorname{beta}\left(x\right)$"
    assert latex(sin(x)) == r"$\operatorname{sin}\left(x\right)$"
    assert latex(sin(x), fold_func_brackets=True) == r"$\operatorname{sin}x$"
    assert latex(sin(2*x**2), fold_func_brackets=True) == \
    r"$\operatorname{sin}2 x^{2}$"
    assert latex(sin(x**2), fold_func_brackets=True) == \
    r"$\operatorname{sin}x^{2}$"

    assert latex(asin(x)**2) == r"$\operatorname{asin}^{2}\left(x\right)$"
    assert latex(asin(x)**2,inv_trig_style="full") == \
        r"$\operatorname{arcsin}^{2}\left(x\right)$"
    assert latex(asin(x)**2,inv_trig_style="power") == \
        r"$\operatorname{sin}^{-1}\left(x\right)^{2}$"
    assert latex(asin(x**2),inv_trig_style="power",fold_func_brackets=True) == \
        r"$\operatorname{sin}^{-1}x^{2}$"

    assert latex(factorial(k)) == r"$k!$"
    assert latex(factorial(-k)) == r"$\left(- k\right)!$"

    assert latex(floor(x)) == r"$\lfloor{x}\rfloor$"
    assert latex(ceiling(x)) == r"$\lceil{x}\rceil$"
    assert latex(abs(x)) == r"$\lvert{x}\rvert$"
    assert latex(re(x)) == r"$\Re{x}$"
    assert latex(im(x)) == r"$\Im{x}$"
    assert latex(conjugate(x)) == r"$\overline{x}$"
    assert latex(gamma(x)) == r"$\operatorname{\Gamma}\left(x\right)$"
    assert latex(Order(x)) == r"$\operatorname{\mathcal{O}}\left(x\right)$"

def test_latex_derivatives():
    assert latex(diff(x**3, x, evaluate=False)) == \
    r"$\frac{\partial}{\partial x} x^{3}$"
    assert latex(diff(sin(x)+x**2, x, evaluate=False)) == \
    r"$\frac{\partial}{\partial x}\left(\operatorname{sin}\left(x\right) + x^{2}\right)$"

def test_latex_integrals():
    assert latex(Integral(log(x), x)) == r"$\int \operatorname{log}\left(x\right)\,dx$"
    assert latex(Integral(x**2, (x,0,1))) == r"$\int_{0}^{1} x^{2}\,dx$"
    assert latex(Integral(x**2, (x,10,20))) == r"$\int_{10}^{20} x^{2}\,dx$"
    assert latex(Integral(y*x**2, (x,0,1), y)) == r"$\int\int_{0}^{1} y x^{2}\,dx dy$"

def test_latex_intervals():
    a = Symbol('a', real=True)
    assert latex(Interval(0, a)) == r"$\left[0, a\right]$"
    assert latex(Interval(0, a, False, False)) == r"$\left[0, a\right]$"
    assert latex(Interval(0, a, True, False)) == r"$\left(0, a\right]$"
    assert latex(Interval(0, a, False, True)) == r"$\left[0, a\right)$"
    assert latex(Interval(0, a, True, True)) == r"$\left(0, a\right)$"

def test_latex_emptyset():
    assert latex(S.EmptySet) == r"$\emptyset$"

def test_latex_union():
    assert latex(Union(Interval(0, 1), Interval(2, 3))) == r"$\left[0, 1\right] \cup \left[2, 3\right]$"

@XFAIL
def test_latex_limits():
    assert latex(limit(x, x, oo, evaluate=False)) == r"$\lim_{x \to \infty}x$"

def test_issue469():
    beta = Symbol(r'\beta')
    y = beta+x
    assert latex(y) in [r'$\beta + x$', r'$x + \beta$']

    beta = Symbol(r'beta')
    y = beta+x
    assert latex(y) in [r'$\beta + x$', r'$x + \beta$']

def test_latex():
    assert latex((2*tau)**Rational(7,2)) == "$8 \\sqrt{2} \\sqrt[7]{\\tau}$"
    assert latex((2*mu)**Rational(7,2), inline=False) == \
            "\\begin{equation*}8 \\sqrt{2} \\sqrt[7]{\\mu}\\end{equation*}"
    assert latex([2/x, y]) =="$\\begin{bmatrix}\\frac{2}{x}, & y\\end{bmatrix}$"


def test_latex_dict():
    d = {Rational(1): 1, x**2: 2, x: 3, x**3: 4}
    assert latex(d) == '$\\begin{Bmatrix}1 : 1, & x : 3, & x^{2} : 2, & x^{3} : 4\\end{Bmatrix}$'

def test_latex_rational():
    #tests issue 874
    assert latex(-Rational(1,2)) == "$- \\frac{1}{2}$"
    assert latex(Rational(-1,2)) == "$- \\frac{1}{2}$"
    assert latex(Rational(1,-2)) == "$- \\frac{1}{2}$"
    assert latex(-Rational(-1,2)) == "$\\frac{1}{2}$"
    assert latex(-Rational(1,2)*x) == "$- \\frac{1}{2} x$"
    assert latex(-Rational(1,2)*x+Rational(-2,3)*y) in [
            "$- \\frac{1}{2} x - \\frac{2}{3} y$",
            "$- \\frac{2}{3} y - \\frac{1}{2} x$",
            ]

def test_latex_inverse():
    #tests issue 1030
    assert latex(1/x) == "$\\frac{1}{x}$"
    assert latex(1/(x+y)) in ["$\\frac{1}{x + y}$", "$\\frac{1}{y + x}$"]

def test_latex_DiracDelta():
    assert latex(DiracDelta(x)) == "$\\delta\\left(x\\right)$"
    assert latex(DiracDelta(x,0)) == "$\\delta\\left(x\\right)$"
    assert latex(DiracDelta(x,5)) == "$\\delta^{\\left( 5 \\right)}\\left( x \\right)$"

def test_latex_Piecewise():
    p = Piecewise((x,x<1),(x**2,True))
    assert latex(p) == "$\\left\\{\\begin{array}{cl} x & for x < 1 \\\\x^{2} &" \
                       " \\textrm{otherwise} \\end{array}\\right.$"

def test_latex_Matrix():
    M = Matrix([[1+x, y],[y, x-1]])
    assert latex(M) == '$\\left(\\begin{smallmatrix}1 + x & y\\\\y & -1 + '\
                       'x\\end{smallmatrix}\\right)$'
    profile = {'mat_str' : 'bmatrix'}
    assert latex(M, profile) == '$\\left(\\begin{bmatrix}1 + x & y\\\\y & -1 + '+ \
           'x\\end{bmatrix}\\right)$'
    profile['mat_delim'] = None
    assert latex(M, profile) == '$\\begin{bmatrix}1 + x & y\\\\y & -1 + '\
                       'x\\end{bmatrix}$'

def test_latex_mul_symbol():
    assert latex(4*4**x, mul_symbol='times') == "$4 \\times 4^{x}$"
    assert latex(4*4**x, mul_symbol='dot') == "$4 \\cdot 4^{x}$"
    assert latex(4*4**x, mul_symbol='ldot') == "$4 \,.\, 4^{x}$"

    assert latex(4*x, mul_symbol='times') == "$4 \\times x$"
    assert latex(4*x, mul_symbol='dot') == "$4 \\cdot x$"
    assert latex(4*x, mul_symbol='ldot') == "$4 \,.\, x$"

def test_latex_issue1282():
    y = 4*4**log(2)
    assert latex(y) == '$4 \\times 4^{\\operatorname{log}\\left(2\\right)}$'
    assert latex(1/y) == '$\\frac{1}{4 \\times 4^{\\operatorname{log}\\left(2\\right)}}$'

def test_latex_issue1477():
    assert latex(Symbol("beta_13_2")) == r"$\beta_{13 2}$"
    assert latex(Symbol("beta_132_20")) == r"$\beta_{132 20}$"
    assert latex(Symbol("beta_13")) == r"$\beta_{13}$"
    assert latex(Symbol("x_a_b")) == r"$x_{a b}$"
    assert latex(Symbol("x_1_2_3")) == r"$x_{1 2 3}$"
    assert latex(Symbol("x_a_b1")) == r"$x_{a b1}$"
    assert latex(Symbol("x_a_1")) == r"$x_{a 1}$"
    assert latex(Symbol("x_1_a")) == r"$x_{1 a}$"
    assert latex(Symbol("x_1^aa")) == r"$x^{aa}_{1}$"
    assert latex(Symbol("x_1__aa")) == r"$x^{aa}_{1}$"
    assert latex(Symbol("x_11^a")) == r"$x^{a}_{11}$"
    assert latex(Symbol("x_11__a")) == r"$x^{a}_{11}$"
    assert latex(Symbol("x_a_a_a_a")) == r"$x_{a a a a}$"
    assert latex(Symbol("x_a_a^a^a")) == r"$x^{a a}_{a a}$"
    assert latex(Symbol("x_a_a__a__a")) == r"$x^{a a}_{a a}$"
    assert latex(Symbol("alpha_11")) == r"$\alpha_{11}$"
    assert latex(Symbol("alpha_11_11")) == r"$\alpha_{11 11}$"
    assert latex(Symbol("alpha_alpha")) == r"$\alpha_{\alpha}$"
    assert latex(Symbol("alpha^aleph")) == r"$\alpha^{\aleph}$"
    assert latex(Symbol("alpha__aleph")) == r"$\alpha^{\aleph}$"


def test_mainvar():
    expr = 3*x*y**3+x**2*y+x**3+y**4
    profile_y = {'mainvar' : y}
    assert latex(expr, profile_y) == '$x^{3} + y x^{2} + 3 x y^{3} + y^{4}$'
    profile_x = {'mainvar' : x}
    assert latex(expr, profile_x) == '$y^{4} + 3 x y^{3} + y x^{2} + x^{3}$'
    profile_y['descending'] = True
    assert latex(expr, profile_y) == '$y^{4} + 3 x y^{3} + y x^{2} + x^{3}$'
    profile_x['descending'] = True
    assert latex(expr, profile_x) == '$x^{3} + y x^{2} + 3 x y^{3} + y^{4}$'
