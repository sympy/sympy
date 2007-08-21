import py

from sympy import *
from sympy.printing.latex import latex

x = Symbol('x')

def test_latex_basic():
    assert latex(1+x) == "$1+x$"
    assert latex(x**2) == "${x}^{2}$"
    assert latex(x**(1+x)) == "${x}^{(1+x)}$"

def test_latex_symbols():
    Gamma, lmbda, rho = map(Symbol, ('Gamma', 'lambda', 'rho'))
    mass, volume = map(Symbol, ('mass', 'volume'))
    assert latex(Gamma + lmbda) == "$\Gamma+\lambda$"
    assert latex(Gamma * lmbda) == "$\Gamma \lambda$"
    assert latex(volume * rho == mass) == r"$\rho \cdot \mathrm{volume} = \mathrm{mass}$"
    assert latex(volume / mass * rho == 1) == r"$\rho \cdot \mathrm{volume} \cdot {\mathrm{mass}}^{(-1)} = 1$"
    assert latex(mass**3 * volume**3) == r"${\mathrm{mass}}^{3} \cdot {\mathrm{volume}}^{3}$"

def test_latex_functions():
    assert latex(exp(x)) == "${e}^{x}$"
    assert latex(exp(I*pi)) == "${e}^{\pi \mathrm{i}}$"

    f = Function('f')
    beta = Function('beta')
    assert latex(f(x)) == r"$f\left(x\right)$"
    assert latex(beta(x)) == r"$\beta\left(x\right)$"
    assert latex(sin(x)) == r"$\mathrm{sin}\left(x\right)$"

    #assert latex(factorial(x, evaluate=False)) == "$x!$"
    #assert latex(factorial(-4, evaluate=False)) == "$(-4)!$"
    #assert latex(factorial(-x, evaluate=False)) == "$(- x)!$"

def test_latex_derivatives():
    assert latex(diff(x**3, x, evaluate=False)) == \
        r"$\frac{\partial}{\partial x} {x}^{3}$"
    assert latex(diff(sin(x)+x**2, x, evaluate=False)) == \
        r"$\frac{\partial}{\partial x} \left({x}^{2}+\mathrm{sin}\left(x\right)\right)$"

def test_latex_integrals():
    assert latex(integrate(log(x), x, evaluate=False)) == r"$\int \mathrm{log}\left(x\right)\,dx$"
    assert latex(integrate(x**2, x==(0,1), evaluate=False)) == r"$\int^0_1 {x}^{2}\,dx$"

#def test_latex_limits():
#    assert latex(limit(x, x, oo, evaluate=False)) == r"$\lim_{x \to \infty}x$"
