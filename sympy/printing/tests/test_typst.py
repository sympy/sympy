from sympy.concrete.summations import Sum
from sympy.core.function import Function
from sympy.core.mul import Mul
from sympy.core.numbers import (I, Rational, oo)
from sympy.core.power import Pow
from sympy.core.singleton import S
from sympy.core.symbol import (Symbol, symbols)
from sympy.functions.elementary.complexes import Abs
from sympy.functions.elementary.exponential import (exp, log)
from sympy.functions.elementary.integers import (ceiling, floor)
from sympy.printing.typst import (typst, translate, greek_letters_set,
                                  typst_greek_dictionary)
from sympy.series.limits import Limit

x, y = symbols('x, y')

def test_printmethod():
    class R(Abs):
        def _typst(self, printer):
            return "foo(%s)" % printer._print(self.args[0])
    assert typst(R(x)) == r"foo(x)"

    class R(Abs):
        def _typst(self, printer):
            return "foo"
    assert typst(R(x)) == r"foo"


def test_typst_basic():
    assert typst(1 + x) == r"x + 1"
    assert typst(x**2) == r"x^(2)"
    assert typst(x**(1 + x)) == r"x^(x + 1)"
    assert typst(x**3 + x + 1 + x**2) == r"x^(3) + x^(2) + x + 1"

    assert typst(2*x*y) == r"2 x y"
    assert typst(2*x*y, mul_symbol='dot') == r"2 dot x dot y"
    assert typst(3*x**2*y, mul_symbol='#h(1cm)') == r"3#h(1cm)x^(2)#h(1cm)y"
    assert typst(1.5*3**x, mul_symbol='#h(1cm)') == r"1.5#h(1cm)3^(x)"

    assert typst(x**S.Half**5) == r"root(32, x)"
    assert typst(Mul(S.Half, x**2, -5, evaluate=False)) == r"1/2 x^(2) (-5)"
    assert typst(Mul(S.Half, x**2, 5, evaluate=False)) == r"1/2 x^(2) 5"
    assert typst(Mul(-5, -5, evaluate=False)) == r"(-5) (-5)"
    assert typst(Mul(5, -5, evaluate=False)) == r"5 (-5)"
    assert typst(Mul(S.Half, -5, S.Half, evaluate=False)) == r"1/2 (-5) 1/2"
    assert typst(Mul(5, I, 5, evaluate=False)) == r"5 i 5"
    assert typst(Mul(5, I, -5, evaluate=False)) == r"5 i (-5)"
    assert typst(Mul(Pow(x, 2), S.Half*x + 1)) == r"x^(2) (x/2 + 1)"
    assert typst(Mul(Pow(x, 3), Rational(2, 3)*x + 1)) == r"x^(3) ((2 x)/3 + 1)"
    assert typst(Mul(Pow(x, 11), 2*x + 1)) == r"x^(11) (2 x + 1)"


def test_typst_limits():
    assert typst(Limit(x, x, oo)) == r"lim_(x -> infinity) x"

    # issue 8175
    f = Function('f')
    assert typst(Limit(f(x), x, 0)) == r"lim_(x -> 0^+) f(x)"
    assert typst(Limit(f(x), x, 0, "-")) == \
        r"lim_(x -> 0^-) f(x)"

    # issue #10806
    assert typst(Limit(f(x), x, 0)**2) == \
        r"(lim_(x -> 0^+) f(x))^(2)"
    # bi-directional limit
    assert typst(Limit(f(x), x, 0, dir='+-')) == \
        r"lim_(x -> 0) f(x)"


def test_typst_log():
    assert typst(log(x)) == r"log(x)"
    assert typst(log(x), ln_notation=True) == r"ln(x)"
    assert typst(log(x) + log(y)) == \
        r"log(x) + log(y)"
    assert typst(log(x) + log(y), ln_notation=True) == \
        r"ln(x) + ln(y)"
    assert typst(pow(log(x), x)) == r"log(x)^(x)"
    assert typst(pow(log(x), x), ln_notation=True) == \
        r"ln(x)^(x)"


def test_typst_sum():
    assert typst(Sum(x*y**2, (x, -2, 2), (y, -5, 5))) == \
        r"sum_(-2 <= x <= 2 \ -5 <= y <= 5) x y^(2)"
    assert typst(Sum(x**2, (x, -2, 2))) == \
        r"sum_(x=-2)^(2) x^(2)"
    assert typst(Sum(x**2 + y, (x, -2, 2))) == \
        r"sum_(x=-2)^(2) (x^(2) + y)"
    assert typst(Sum(x**2 + y, (x, -2, 2))**2) == \
        r"(sum_(x=-2)^(2) (x^(2) + y))^(2)"


def test_typst_symbols():
    Gamma, lmbda, rho = symbols('Gamma, lambda, rho')
    tau, Tau, TAU, taU = symbols('tau, Tau, TAU, taU')
    assert typst(tau) == r"tau"
    assert typst(Tau) == r"Tau"
    assert typst(TAU) == r"tau"
    assert typst(taU) == r"tau"
    # Check that all capitalized greek letters are handled explicitly
    capitalized_letters = {l.capitalize() for l in greek_letters_set}
    assert len(capitalized_letters - set(typst_greek_dictionary.keys())) == 0
    assert typst(Gamma + lmbda) == r"Gamma + lambda"
    assert typst(Gamma * lmbda) == r"Gamma lambda"
    assert typst(Symbol('q1')) == r"q_(1)"
    assert typst(Symbol('q21')) == r"q_(21)"
    assert typst(Symbol('epsilon0')) == r"epsilon.alt_(0)"
    assert typst(Symbol('omega1')) == r"omega_(1)"
    assert typst(Symbol('91')) == r"91"
    assert typst(Symbol('alpha_new')) == r"alpha_(new)"
    assert typst(Symbol('C^orig')) == r"C^(orig)"
    assert typst(Symbol('x^alpha')) == r"x^(alpha)"
    assert typst(Symbol('beta^alpha')) == r"beta^(alpha)"
    assert typst(Symbol('e^Alpha')) == r"e^(Alpha)"
    assert typst(Symbol('omega_alpha^beta')) == r"omega^(beta)_(alpha)"
    assert typst(Symbol('omega') ** Symbol('beta')) == r"omega^(beta)"


def test_typst_functions():
    assert typst(exp(x)) == r"e^(x)"
    assert typst(exp(1) + exp(2)) == r"e + e^(2)"


    assert typst(Abs(x)) == r"abs(x)"
    assert typst(Abs(x)**2) == r"abs(x)^(2)"

    assert typst(floor(x)) == r"floor(x)"
    assert typst(ceiling(x)) == r"ceil(x)"
    assert typst(floor(x)**2) == r"floor(x)^(2)"
    assert typst(ceiling(x)**2) == r"ceil(x)^(2)"


def test_translate():
    s = 'Alpha'
    assert translate(s) == r'Alpha'
    s = 'Beta'
    assert translate(s) == r'Beta'
    s = 'Eta'
    assert translate(s) == r'Eta'
    s = 'omicron'
    assert translate(s) == r'omicron'
    s = 'Pi'
    assert translate(s) == r'Pi'
    s = 'pi'
    assert translate(s) == r'pi'
    s = 'LamdaHatDOT'
    assert translate(s) == r'accent(accent(Lambda, hat), dot)'


def test_other_symbols():
    from sympy.printing.typst import other_symbols
    for s in other_symbols:
        assert typst(symbols(s)) == r"" + s


def test_modifiers():
    # Test each modifier individually in the simplest case
    # (with funny capitalizations)
    assert typst(symbols("xMathring")) == r"accent(x, circle)"
    assert typst(symbols("xCheck")) == r"accent(x, caron)"
    assert typst(symbols("xBreve")) == r"accent(x, breve)"
    assert typst(symbols("xAcute")) == r"accent(x, acute)"
    assert typst(symbols("xGrave")) == r"accent(x, grave)"
    assert typst(symbols("xTilde")) == r"accent(x, tilde)"
    assert typst(symbols("xPrime")) == r"(x)'"
    assert typst(symbols("xddDDot")) == r"accent(x, dot.quad)"
    assert typst(symbols("xDdDot")) == r"accent(x, dot.triple)"
    assert typst(symbols("xDDot")) == r"accent(x, dot.double)"
    assert typst(symbols("xBold")) == r"bold(x)"
    assert typst(symbols("xnOrM")) == r"norm(x)"
    assert typst(symbols("xAVG")) == r"lr(angle.l x angle.r)"
    assert typst(symbols("xHat")) == r"accent(x, hat)"
    assert typst(symbols("xDot")) == r"accent(x, dot)"
    assert typst(symbols("xBar")) == r"accent(x, macron)"
    assert typst(symbols("xVec")) == r"accent(x, arrow)"
    assert typst(symbols("xAbs")) == r"abs(x)"
    assert typst(symbols("xMag")) == r"abs(x)"
    assert typst(symbols("xPrM")) == r"(x)'"
    assert typst(symbols("xBM")) == r"bold(x)"
    # Test strings that are *only* the names of modifiers
    assert typst(symbols("Mathring")) == r"Mathring"
    assert typst(symbols("Check")) == r"Check"
    assert typst(symbols("Breve")) == r"Breve"
    assert typst(symbols("Acute")) == r"Acute"
    assert typst(symbols("Grave")) == r"Grave"
    assert typst(symbols("Tilde")) == r"Tilde"
    assert typst(symbols("Prime")) == r"Prime"
    assert typst(symbols("DDot")) == r"accent(D, dot)"
    assert typst(symbols("Bold")) == r"Bold"
    assert typst(symbols("NORm")) == r"NORm"
    assert typst(symbols("AVG")) == r"AVG"
    assert typst(symbols("Hat")) == r"Hat"
    assert typst(symbols("Dot")) == r"Dot"
    assert typst(symbols("Bar")) == r"Bar"
    assert typst(symbols("Vec")) == r"Vec"
    assert typst(symbols("Abs")) == r"Abs"
    assert typst(symbols("Mag")) == r"Mag"
    assert typst(symbols("PrM")) == r"PrM"
    assert typst(symbols("BM")) == r"BM"
    assert typst(symbols("hbar")) == r"plank.reduce"
    # Check a few combinations
    assert typst(symbols("xvecdot")) == r"accent(accent(x, arrow), dot)"
    assert typst(symbols("xDotVec")) == r"accent(accent(x, dot), arrow)"
    assert typst(symbols("xHATNorm")) == r"norm(accent(x, hat))"
    # Check a couple big, ugly combinations
    assert typst(symbols('xMathringBm_yCheckPRM__zbreveAbs')) == \
        r"bold(accent(x, circle))^(abs(accent(z, breve)))_((accent(y, caron))')"
    assert typst(symbols('alphadothat_nVECDOT__tTildePrime')) == \
        r"accent(accent(alpha, dot), hat)^((accent(t, tilde))')_(accent(accent(n, arrow), dot))"


def test_greek_symbols():
    assert typst(Symbol('alpha'))   == 'alpha'
    assert typst(Symbol('beta'))    == 'beta'
    assert typst(Symbol('gamma'))   == 'gamma'
    assert typst(Symbol('delta'))   == 'delta'
    assert typst(Symbol('epsilon')) == 'epsilon.alt'
    assert typst(Symbol('zeta'))    == 'zeta'
    assert typst(Symbol('eta'))     == 'eta'
    assert typst(Symbol('theta'))   == 'theta'
    assert typst(Symbol('iota'))    == 'iota'
    assert typst(Symbol('kappa'))   == 'kappa'
    assert typst(Symbol('lambda'))  == 'lambda'
    assert typst(Symbol('mu'))      == 'mu'
    assert typst(Symbol('nu'))      == 'nu'
    assert typst(Symbol('xi'))      == 'xi'
    assert typst(Symbol('omicron')) == 'omicron'
    assert typst(Symbol('pi'))      == 'pi'
    assert typst(Symbol('rho'))     == 'rho'
    assert typst(Symbol('sigma'))   == 'sigma'
    assert typst(Symbol('tau'))     == 'tau'
    assert typst(Symbol('upsilon')) == 'upsilon'
    assert typst(Symbol('phi'))     == 'phi.alt'
    assert typst(Symbol('chi'))     == 'chi'
    assert typst(Symbol('psi'))     == 'psi'
    assert typst(Symbol('omega'))   == 'omega'

    assert typst(Symbol('Alpha'))   == 'Alpha'
    assert typst(Symbol('Beta'))    == 'Beta'
    assert typst(Symbol('Gamma'))   == 'Gamma'
    assert typst(Symbol('Delta'))   == 'Delta'
    assert typst(Symbol('Epsilon')) == 'Epsilon'
    assert typst(Symbol('Zeta'))    == 'Zeta'
    assert typst(Symbol('Eta'))     == 'Eta'
    assert typst(Symbol('Theta'))   == 'Theta'
    assert typst(Symbol('Iota'))    == 'Iota'
    assert typst(Symbol('Kappa'))   == 'Kappa'
    assert typst(Symbol('Lambda'))  == 'Lambda'
    assert typst(Symbol('Mu'))      == 'Mu'
    assert typst(Symbol('Nu'))      == 'Nu'
    assert typst(Symbol('Xi'))      == 'Xi'
    assert typst(Symbol('Omicron')) == 'Omicron'
    assert typst(Symbol('Pi'))      == 'Pi'
    assert typst(Symbol('Rho'))     == 'Rho'
    assert typst(Symbol('Sigma'))   == 'Sigma'
    assert typst(Symbol('Tau'))     == 'Tau'
    assert typst(Symbol('Upsilon')) == 'Upsilon'
    assert typst(Symbol('Phi'))     == 'Phi'
    assert typst(Symbol('Chi'))     == 'Chi'
    assert typst(Symbol('Psi'))     == 'Psi'
    assert typst(Symbol('Omega'))   == 'Omega'

    assert typst(Symbol('varepsilon')) == 'epsilon'
    assert typst(Symbol('varkappa')) == 'kappa'
    assert typst(Symbol('varphi')) == 'phi'
    assert typst(Symbol('varpi')) == 'pi.alt'
    assert typst(Symbol('varrho')) == 'rho.alt'
    assert typst(Symbol('varsigma')) == 'sigma.alt'
    assert typst(Symbol('vartheta')) == 'theta.alt'
