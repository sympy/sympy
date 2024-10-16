from sympy.combinatorics.permutations import Cycle, Permutation, AppliedPermutation
from sympy.concrete.summations import Sum
from sympy.core.function import (Derivative, Function, diff)
from sympy.core.mul import Mul
from sympy.core.numbers import (Float, I, Rational, oo)
from sympy.core.power import Pow
from sympy.core.singleton import S
from sympy.core.symbol import (Symbol, symbols)
from sympy.functions.elementary.complexes import (Abs, arg, conjugate, im, polar_lift, re)
from sympy.functions.elementary.exponential import (exp, log)
from sympy.functions.elementary.integers import (ceiling, floor)
from sympy.functions.elementary.miscellaneous import Max
from sympy.functions.elementary.trigonometric import sin
from sympy.logic.boolalg import (And, Or, Xor, Equivalent, false, Not, true)
from sympy.printing.typst import (typst, translate, greek_letters_set,
                                  typst_greek_dictionary)
from sympy.series.limits import Limit
from sympy.vector import CoordSys3D, Cross, Curl, Dot, Divergence, Gradient, Laplacian


from sympy.testing.pytest import warns_deprecated_sympy

x, y = symbols('x, y')
n = symbols('n', integer=True)

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


def test_latex_builtins():
    assert typst(True) == r'upright("True")'
    assert typst(False) == r'upright("False")'
    assert typst(None) == r'upright("None")'
    assert typst(true) == r'upright("True")'
    assert typst(false) == r'upright("False")'


def test_typst_cycle():
    assert typst(Cycle(1, 2, 4)) == r"(1 space 2 space 4)"
    assert typst(Cycle(1, 2)(4, 5, 6)) == \
        r"(1 space 2)(4 space 5 space 6)"
    assert typst(Cycle()) == r"()"


def test_typst_permutation():
    assert typst(Permutation(1, 2, 4)) == r"(1 space 2 space 4)"
    assert typst(Permutation(1, 2)(4, 5, 6)) == \
        r"(1 space 2)(4 space 5 space 6)"
    assert typst(Permutation()) == r"()"
    assert typst(Permutation(2, 4)*Permutation(5)) == \
        r"(2 space 4)(5)"
    assert typst(Permutation(5)) == r"(5)"

    assert typst(Permutation(0, 1), perm_cyclic=False) == \
        r"mat(0, 1; 1, 0)"
    assert typst(Permutation(0, 1)(2, 3), perm_cyclic=False) == \
        r"mat(0, 1, 2, 3; 1, 0, 3, 2)"
    assert typst(Permutation(), perm_cyclic=False) == \
        r"()"

    with warns_deprecated_sympy():
        old_print_cyclic = Permutation.print_cyclic
        Permutation.print_cyclic = False
        assert typst(Permutation(0, 1)(2, 3)) == \
            r"mat(0, 1, 2, 3; 1, 0, 3, 2)"
        Permutation.print_cyclic = old_print_cyclic

def test_latex_Float():
    assert typst(Float(1.0e100)) == r"1.0 dot 10^(100)"
    assert typst(Float(1.0e-100)) == r"1.0 dot 10^(-100)"
    assert typst(Float(1.0e-100), mul_symbol="times") == \
        r"1.0 times 10^(-100)"
    assert typst(Float('10000.0'), full_prec=False, min=-2, max=2) == \
        r"1.0 dot 10^(4)"
    assert typst(Float('10000.0'), full_prec=False, min=-2, max=4) == \
        r"1.0 dot 10^(4)"
    assert typst(Float('10000.0'), full_prec=False, min=-2, max=5) == \
        r"10000.0"
    assert typst(Float('0.099999'), full_prec=True,  min=-2, max=5) == \
        r"9.99990000000000 dot 10^(-2)"

def test_latex_vector_expressions():
    A = CoordSys3D('A')

    assert typst(Cross(A.i, A.j*A.x*3+A.k)) == \
        r"bold(hat(i)_(A)) times ((3 bold(x_(A)))bold(hat(j)_(A)) + bold(hat(k)_(A)))"
    assert typst(Cross(A.i, A.j)) == \
        r"bold(hat(i)_(A)) times bold(hat(j)_(A))"
    assert typst(x*Cross(A.i, A.j)) == \
        r"x (bold(hat(i)_(A)) times bold(hat(j)_(A)))"
    assert typst(Cross(x*A.i, A.j)) == \
        r'- bold(hat(j)_(A)) times ((x)bold(hat(i)_(A)))'

    assert typst(Curl(3*A.x*A.j)) == \
        r"nabla times ((3 bold(x_(A)))bold(hat(j)_(A)))"
    assert typst(Curl(3*A.x*A.j+A.i)) == \
        r"nabla times (bold(hat(i)_(A)) + (3 bold(x_(A)))bold(hat(j)_(A)))"
    assert typst(Curl(3*x*A.x*A.j)) == \
        r"nabla times ((3 bold(x_(A)) x)bold(hat(j)_(A)))"
    assert typst(x*Curl(3*A.x*A.j)) == \
        r"x (nabla times ((3 bold(x_(A)))bold(hat(j)_(A))))"

    assert typst(Divergence(3*A.x*A.j+A.i)) == \
        r"nabla dot (bold(hat(i)_(A)) + (3 bold(x_(A)))bold(hat(j)_(A)))"
    assert typst(Divergence(3*A.x*A.j)) == \
        r"nabla dot ((3 bold(x_(A)))bold(hat(j)_(A)))"
    assert typst(x*Divergence(3*A.x*A.j)) == \
        r"x (nabla dot ((3 bold(x_(A)))bold(hat(j)_(A))))"

    assert typst(Dot(A.i, A.j*A.x*3+A.k)) == \
        r"bold(hat(i)_(A)) dot ((3 bold(x_(A)))bold(hat(j)_(A)) + bold(hat(k)_(A)))"
    assert typst(Dot(A.i, A.j)) == \
        r"bold(hat(i)_(A)) dot bold(hat(j)_(A))"
    assert typst(Dot(x*A.i, A.j)) == \
        r"bold(hat(j)_(A)) dot ((x)bold(hat(i)_(A)))"
    assert typst(x*Dot(A.i, A.j)) == \
        r"x (bold(hat(i)_(A)) dot bold(hat(j)_(A)))"

    assert typst(Gradient(A.x)) == r"nabla bold(x_(A))"
    assert typst(Gradient(A.x + 3*A.y)) == \
        r"nabla (bold(x_(A)) + 3 bold(y_(A)))"
    assert typst(x*Gradient(A.x)) == r"x (nabla bold(x_(A)))"
    assert typst(Gradient(x*A.x)) == r"nabla (bold(x_(A)) x)"

    assert typst(Laplacian(A.x)) == r"Delta bold(x_(A))"
    assert typst(Laplacian(A.x + 3*A.y)) == \
        r"Delta (bold(x_(A)) + 3 bold(y_(A)))"
    assert typst(x*Laplacian(A.x)) == r"x (Delta bold(x_(A)))"
    assert typst(Laplacian(x*A.x)) == r"Delta (bold(x_(A)) x)"

def test_AppliedPermutation():
    p = Permutation(0, 1, 2)
    x = Symbol('x')
    assert typst(AppliedPermutation(p, x)) == \
        r'sigma_((0 space 1 space 2))(x)'

def test_imaginary_unit():
    assert typst(1 + I) == r'1 + i'
    assert typst(1 + I, imaginary_unit='i') == r'1 + i'
    assert typst(1 + I, imaginary_unit='j') == r'1 + j'
    assert typst(1 + I, imaginary_unit='foo') == r'1 + foo'
    assert typst(I, imaginary_unit="ti") == r'text(i)'
    assert typst(I, imaginary_unit="tj") == r'text(j)'


def test_text_re_im():
    assert typst(im(x), gothic_re_im=True) == r'Im(x)'
    assert typst(im(x), gothic_re_im=False) == r'upright("im")(x)'
    assert typst(re(x), gothic_re_im=True) == r'Re(x)'
    assert typst(re(x), gothic_re_im=False) == r'upright("re")(x)'

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

    assert typst(conjugate(x)) == r"overline(x)"
    assert typst(conjugate(x)**2) == r"overline(x)^(2)"
    assert typst(conjugate(x**2)) == r"overline(x)^(2)"

    assert typst(arg(x)) == r'arg(x)'

    # Test latex printing of function names with "_"
    assert typst(polar_lift(0)) == \
        r'upright("polar_lift")(0)'
    assert typst(polar_lift(0)**3) == \
        r'upright("polar_lift")^(3)(0)'



def test_boolean_args_order():
    syms = symbols('a:f')

    expr = And(*syms)
    assert typst(expr) == r'a and b and c and d and e and f'

    expr = Or(*syms)
    assert typst(expr) == r'a or b or c or d or e or f'

    expr = Equivalent(*syms)
    assert typst(expr) == \
        r'a <=> b <=> c <=> d <=> e <=> f'

    expr = Xor(*syms)
    assert typst(expr) == \
        r'a \u{22BB} b \u{22BB} c \u{22BB} d \u{22BB} e \u{22BB} f'

def test_issue_7180():
    assert typst(Equivalent(x, y)) == r"x <=> y"
    assert typst(Not(Equivalent(x, y))) == r"x arrow.l.r.double.not y"


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

def test_typst_derivatives():
    # regular "d" for ordinary derivatives
    assert typst(diff(x**3, x, evaluate=False)) == \
        r"d / (d x) x^(3)"
    assert typst(diff(sin(x) + x**2, x, evaluate=False)) == \
        r"d / (d x) (x^(2) + sin(x))"
    assert typst(diff(diff(sin(x) + x**2, x, evaluate=False), evaluate=False))\
        == \
        r"d^(2) / (d x^(2)) (x^(2) + sin(x))"
    assert typst(diff(diff(diff(sin(x) + x**2, x, evaluate=False), evaluate=False), evaluate=False)) == \
        r"d^(3) / (d x^(3)) (x^(2) + sin(x))"

    # partial for partial derivatives
    assert typst(diff(sin(x * y), x, evaluate=False)) == \
        r"partial / (partial x) sin(x y)"
    assert typst(diff(sin(x * y) + x**2, x, evaluate=False)) == \
        r"partial / (partial x) (x^(2) + sin(x y))"
    assert typst(diff(diff(sin(x*y) + x**2, x, evaluate=False), x, evaluate=False)) == \
        r"partial^(2) / (partial x^(2)) (x^(2) + sin(x y))"
    assert typst(diff(diff(diff(sin(x*y) + x**2, x, evaluate=False), x, evaluate=False), x, evaluate=False)) == \
        r"partial^(3) / (partial x^(3)) (x^(2) + sin(x y))"

    # mixed partial derivatives
    f = Function("f")
    assert typst(diff(diff(f(x, y), x, evaluate=False), y, evaluate=False)) == \
        r"partial^(2) / (partial y partial x) " + typst(f(x, y))

    assert typst(diff(diff(diff(f(x, y), x, evaluate=False), x, evaluate=False), y, evaluate=False)) == \
        r"partial^(3) / (partial y partial x^(2)) " + typst(f(x, y))

    # for negative nested Derivative
    assert typst(diff(-diff(y**2,x,evaluate=False),x,evaluate=False)) == r'd / (d x) (- d / (d x) y^(2))'
    assert typst(diff(diff(-diff(diff(y,x,evaluate=False),x,evaluate=False),x,evaluate=False),x,evaluate=False)) == \
        r'd^(2) / (d x^(2)) (- d^(2) / (d x^(2)) y)'

    # # use ordinary d when one of the variables has been integrated out
    # assert typst(diff(Integral(exp(-x*y), (x, 0, oo)), y, evaluate=False)) == \
    #     r"(d / d y) int(limits_(0)^(infty) e^(- x y), dx"

    # Derivative wrapped in power:
    assert typst(diff(x, x, evaluate=False)**2) == \
        r"(d / (d x) x)^(2)"

    assert typst(diff(f(x), x)**2) == \
        r"(d / (d x) f(x))^(2)"

    assert typst(diff(f(x), (x, n))) == \
        r"d^(n) / (d x^(n)) f(x)"

    x1 = Symbol('x1')
    x2 = Symbol('x2')
    assert typst(diff(f(x1, x2), x1)) == r'partial / (partial x_(1)) f(x_(1),x_(2))'

    n1 = Symbol('n1')
    assert typst(diff(f(x), (x, n1))) == r'd^(n_(1)) / (d x^(n_(1))) f(x)'

    n2 = Symbol('n2')
    assert typst(diff(f(x), (x, Max(n1, n2)))) == \
        r'd^(max(n_(1), n_(2))) / (d x^(max(n_(1), n_(2)))) f(x)'

    # set diff operator
    assert typst(diff(f(x), x), diff_operator="rd") == r'bold(d) / (bold(d) x) f(x)'

def test_issue_17092():
    x_star = Symbol('x^*')
    assert typst(Derivative(x_star, x_star,2)) == r'd^(2) / (d (x^(*))^(2)) x^(*)'
