import py

from sympy import symbols, Rational, Symbol, Integral, log, diff, sin, exp, \
        Function, factorial
from sympy.abc import mu, tau
from sympy.printing.latex import latex
from sympy.utilities.pytest import XFAIL

x,y = symbols('xy')
k,n = symbols('kn', integer=True)

def test_latex_basic():
    assert latex(1+x) in ["$1 + x$", '$x + 1$']
    assert latex(x**2) == "${x}^{2}$"
    assert latex(x**(1+x)) in ["${x}^{1 + x}$", '${x}^{x + 1}$']

def test_latex_symbols():
    Gamma, lmbda, rho = map(Symbol, ('Gamma', 'lambda', 'rho'))
    mass, volume = map(Symbol, ('mass', 'volume'))
    assert latex(Gamma + lmbda) in [r"$\Gamma + \lambda$", '$\lambda + \Gamma$']
    assert latex(Gamma * lmbda) == r"$\Gamma \lambda$"
    #assert latex(volume * rho == mass) == r"$\rho \mathrm{volume} = \mathrm{mass}$"
    #assert latex(volume / mass * rho == 1) == r"$\rho \mathrm{volume} {\mathrm{mass}}^{(-1)} = 1$"
    #assert latex(mass**3 * volume**3) == r"${\mathrm{mass}}^{3} \cdot {\mathrm{volume}}^{3}$"

def test_latex_functions():
    assert latex(exp(x)) == "${e}^{x}$"

    f = Function('f')
    assert latex(f(x)) == '$\\operatorname{f}\\left(x\\right)$'

    beta = Function('beta')

    assert latex(beta(x)) == r"$\operatorname{beta}\left(x\right)$"
    assert latex(sin(x)) == r"$\operatorname{sin}\left(x\right)$"

    assert latex(factorial(k)) == r"$k!$"
    assert latex(factorial(-k)) == r"$\left(- k\right)!$"

def test_latex_derivatives():
    assert latex(diff(x**3, x, evaluate=False)) == \
    r"$\frac{\partial}{\partial x} {x}^{3}$"
    assert latex(diff(sin(x)+x**2, x, evaluate=False)) in [
        r"$\frac{\partial}{\partial x}\left({x}^{2} + \operatorname{sin}\left(x\right)\right)$",
        r'$\frac{\partial}{\partial x}\left(\operatorname{sin}\left(x\right) + {x}^{2}\right)$']

def test_latex_integrals():
    assert latex(Integral(log(x), x)) == r"$\int \operatorname{log}\left(x\right)\,dx$"
    assert latex(Integral(x**2, (x,0,1))) == r"$\int_{0}^{1} {x}^{2}\,dx$"
    assert latex(Integral(x**2, (x,10,20))) == r"$\int_{10}^{20} {x}^{2}\,dx$"
    assert latex(Integral(y*x**2, (x,0,1), y)) == r"$\int\int_{0}^{1} y {x}^{2}\,dx dy$"

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
