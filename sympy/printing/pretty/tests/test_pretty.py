# -*- coding: utf-8 -*-

from sympy import (
    Add, And, Basic, Derivative, Dict, Eq, Equivalent, FF,
    FiniteSet, Function, Ge, Gt, I, Implies, Integral,
    Lambda, Le, Limit, Lt, Matrix, Mul, Nand, Ne, Nor, Not, O, Or,
    Pow, Product, QQ, RR, Rational, Ray, RootOf, RootSum, S,
    Segment, Subs, Sum, Symbol, Tuple, Xor, ZZ, conjugate,
    groebner, oo, pi, symbols, ilex, grlex, Range, Contains)
from sympy.functions import (Abs, Chi, Ci, Ei, KroneckerDelta,
    Piecewise, Shi, Si, atan2, binomial, catalan, ceiling, cos,
    euler, exp, expint, factorial, factorial2, floor, gamma, hyper, log,
    lowergamma, meijerg, sin, sqrt, subfactorial, tan, uppergamma,
    elliptic_k, elliptic_f, elliptic_e, elliptic_pi)

from sympy.printing.pretty import pretty as xpretty
from sympy.printing.pretty import pprint

from sympy.physics.units import joule

from sympy.utilities.pytest import raises, XFAIL
from sympy.core.trace import Tr

from sympy.core.compatibility import u_decode as u
from sympy.core.compatibility import range

a, b, x, y, z, k = symbols('a,b,x,y,z,k')
th = Symbol('theta')
ph = Symbol('phi')

"""
Expressions whose pretty-printing is tested here:
(A '#' to the right of an expression indicates that its various acceptable
orderings are accounted for by the tests.)


BASIC EXPRESSIONS:

oo
(x**2)
1/x
y*x**-2
x**Rational(-5,2)
(-2)**x
Pow(3, 1, evaluate=False)
(x**2 + x + 1)  #
1-x  #
1-2*x  #
x/y
-x/y
(x+2)/y  #
(1+x)*y  #3
-5*x/(x+10)  # correct placement of negative sign
1 - Rational(3,2)*(x+1)
-(-x + 5)*(-x - 2*sqrt(2) + 5) - (-y + 5)*(-y + 5) # issue 5524


ORDERING:

x**2 + x + 1
1 - x
1 - 2*x
2*x**4 + y**2 - x**2 + y**3


RELATIONAL:

Eq(x, y)
Lt(x, y)
Gt(x, y)
Le(x, y)
Ge(x, y)
Ne(x/(y+1), y**2)  #


RATIONAL NUMBERS:

y*x**-2
y**Rational(3,2) * x**Rational(-5,2)
sin(x)**3/tan(x)**2


FUNCTIONS (ABS, CONJ, EXP, FUNCTION BRACES, FACTORIAL, FLOOR, CEILING):

(2*x + exp(x))  #
Abs(x)
Abs(x/(x**2+1)) #
Abs(1 / (y - Abs(x)))
factorial(n)
factorial(2*n)
subfactorial(n)
subfactorial(2*n)
factorial(factorial(factorial(n)))
factorial(n+1) #
conjugate(x)
conjugate(f(x+1)) #
f(x)
f(x, y)
f(x/(y+1), y) #
sin(x)**2
conjugate(a+b*I)
conjugate(exp(a+b*I))
conjugate( f(1 + conjugate(f(x))) ) #
f(x/(y+1), y)  # denom of first arg
floor(1 / (y - floor(x)))
ceiling(1 / (y - ceiling(x)))


SQRT:

sqrt(2)
2**Rational(1,3)
2**Rational(1,1000)
sqrt(x**2 + 1)
(1 + sqrt(5))**Rational(1,3)
2**(1/x)
sqrt(2+pi)
(2+(1+x**2)/(2+x))**Rational(1,4)+(1+x**Rational(1,1000))/sqrt(3+x**2)


DERIVATIVES:

Derivative(log(x), x, evaluate=False)
Derivative(log(x), x, evaluate=False) + x  #
Derivative(log(x) + x**2, x, y, evaluate=False)
Derivative(2*x*y, y, x, evaluate=False) + x**2  #
beta(alpha).diff(alpha)


INTEGRALS:

Integral(log(x), x)
Integral(x**2, x)
Integral((sin(x))**2 / (tan(x))**2)
Integral(x**(2**x), x)
Integral(x**2, (x,1,2))
Integral(x**2, (x,Rational(1,2),10))
Integral(x**2*y**2, x,y)
Integral(x**2, (x, None, 1))
Integral(x**2, (x, 1, None))
Integral(sin(th)/cos(ph), (th,0,pi), (ph, 0, 2*pi))


MATRICES:

Matrix([[x**2+1, 1], [y, x+y]])  #
Matrix([[x/y, y, th], [0, exp(I*k*ph), 1]])


PIECEWISE:

Piecewise((x,x<1),(x**2,True))


SEQUENCES (TUPLES, LISTS, DICTIONARIES):

()
[]
{}
(1/x,)
[x**2, 1/x, x, y, sin(th)**2/cos(ph)**2]
(x**2, 1/x, x, y, sin(th)**2/cos(ph)**2)
{x: sin(x)}
{1/x: 1/y, x: sin(x)**2}  #
[x**2]
(x**2,)
{x**2: 1}


LIMITS:

Limit(x, x, oo)
Limit(x**2, x, 0)
Limit(1/x, x, 0)
Limit(sin(x)/x, x, 0)


UNITS:

joule => kg*m**2/s


SUBS:

Subs(f(x), x, ph**2)
Subs(f(x).diff(x), x, 0)
Subs(f(x).diff(x)/y, (x, y), (0, Rational(1, 2)))


ORDER:

O(1)
O(1/x)
O(x**2 + y**2)

"""


def pretty(expr, order=None):
    """ASCII pretty-printing"""
    return xpretty(expr, order=order, use_unicode=False, wrap_line=False)


def upretty(expr, order=None):
    """Unicode pretty-printing"""
    return xpretty(expr, order=order, use_unicode=True, wrap_line=False)


def test_pretty_ascii_str():
    assert pretty( 'xxx' ) == 'xxx'
    assert pretty( "xxx" ) == 'xxx'
    assert pretty( 'xxx\'xxx' ) == 'xxx\'xxx'
    assert pretty( 'xxx"xxx' ) == 'xxx\"xxx'
    assert pretty( 'xxx\"xxx' ) == 'xxx\"xxx'
    assert pretty( "xxx'xxx" ) == 'xxx\'xxx'
    assert pretty( "xxx\'xxx" ) == 'xxx\'xxx'
    assert pretty( "xxx\"xxx" ) == 'xxx\"xxx'
    assert pretty( "xxx\"xxx\'xxx" ) == 'xxx"xxx\'xxx'
    assert pretty( "xxx\nxxx" ) == 'xxx\nxxx'


def test_pretty_unicode_str():
    assert pretty( u('xxx') ) == u('xxx')
    assert pretty( u('xxx') ) == u('xxx')
    assert pretty( u('xxx\'xxx') ) == u('xxx\'xxx')
    assert pretty( u('xxx"xxx') ) == u('xxx\"xxx')
    assert pretty( u('xxx\"xxx') ) == u('xxx\"xxx')
    assert pretty( u("xxx'xxx") ) == u('xxx\'xxx')
    assert pretty( u("xxx\'xxx") ) == u('xxx\'xxx')
    assert pretty( u("xxx\"xxx") ) == u('xxx\"xxx')
    assert pretty( u("xxx\"xxx\'xxx") ) == u('xxx"xxx\'xxx')
    assert pretty( u("xxx\nxxx") ) == u('xxx\nxxx')


def test_upretty_greek():
    assert upretty( oo ) == u('∞')
    assert upretty( Symbol('alpha^+_1') ) == u('α⁺₁')
    assert upretty( Symbol('beta') ) == u('β')
    assert upretty(Symbol('lambda')) == u('λ')


def test_upretty_multiindex():
    assert upretty( Symbol('beta12') ) == u('β₁₂')
    assert upretty( Symbol('Y00') ) == u('Y₀₀')
    assert upretty( Symbol('Y_00') ) == u('Y₀₀')
    assert upretty( Symbol('F^+-') ) == u('F⁺⁻')


def test_upretty_sub_super():
    assert upretty( Symbol('beta_1_2') ) == u('β₁ ₂')
    assert upretty( Symbol('beta^1^2') ) == u('β¹ ²')
    assert upretty( Symbol('beta_1^2') ) == u('β²₁')
    assert upretty( Symbol('beta_10_20') ) == u('β₁₀ ₂₀')
    assert upretty( Symbol('beta_ax_gamma^i') ) == u('βⁱₐₓ ᵧ')
    assert upretty( Symbol("F^1^2_3_4") ) == u('F¹ ²₃ ₄')
    assert upretty( Symbol("F_1_2^3^4") ) == u('F³ ⁴₁ ₂')
    assert upretty( Symbol("F_1_2_3_4") ) == u('F₁ ₂ ₃ ₄')
    assert upretty( Symbol("F^1^2^3^4") ) == u('F¹ ² ³ ⁴')


def test_upretty_subs_missingin_24():
    assert upretty( Symbol('F_beta') ) == u('Fᵦ')
    assert upretty( Symbol('F_gamma') ) == u('Fᵧ')
    assert upretty( Symbol('F_rho') ) == u('Fᵨ')
    assert upretty( Symbol('F_phi') ) == u('Fᵩ')
    assert upretty( Symbol('F_chi') ) == u('Fᵪ')

    assert upretty( Symbol('F_a') ) == u('Fₐ')
    assert upretty( Symbol('F_e') ) == u('Fₑ')
    assert upretty( Symbol('F_i') ) == u('Fᵢ')
    assert upretty( Symbol('F_o') ) == u('Fₒ')
    assert upretty( Symbol('F_u') ) == u('Fᵤ')
    assert upretty( Symbol('F_r') ) == u('Fᵣ')
    assert upretty( Symbol('F_v') ) == u('Fᵥ')
    assert upretty( Symbol('F_x') ) == u('Fₓ')

def test_upretty_modifiers():
    # Accents
    assert upretty( Symbol('Fmathring') ) == u('F̊')
    assert upretty( Symbol('Fddddot') ) == u('F̈̈')
    assert upretty( Symbol('Fdddot') ) == u('F̈̇')
    assert upretty( Symbol('Fddot') ) == u('F̈')
    assert upretty( Symbol('Fdot') ) == u('Ḟ')
    assert upretty( Symbol('Fcheck') ) == u('F̌')
    assert upretty( Symbol('Fbreve') ) == u('F̆')
    assert upretty( Symbol('Facute') ) == u('F́')
    assert upretty( Symbol('Fgrave') ) == u('F̀')
    assert upretty( Symbol('Ftilde') ) == u('F̃')
    assert upretty( Symbol('Fhat') ) == u('F̂')
    assert upretty( Symbol('Fbar') ) == u('F̅')
    assert upretty( Symbol('Fvec') ) == u('F⃗')
    assert upretty( Symbol('Fprime') ) == u('F ̍')
    assert upretty( Symbol('Fprm') ) == u('F ̍')
    # No faces are actually implemented, but test to make sure the modifiers are stripped
    assert upretty( Symbol('Fbold') ) == u('Fbold')
    assert upretty( Symbol('Fbm') ) == u('Fbm')
    assert upretty( Symbol('Fcal') ) == u('Fcal')
    assert upretty( Symbol('Fscr') ) == u('Fscr')
    assert upretty( Symbol('Ffrak') ) == u('Ffrak')
    # Brackets
    assert upretty( Symbol('Fnorm') ) == u('‖F‖')
    assert upretty( Symbol('Favg') ) == u('⟨F⟩')
    assert upretty( Symbol('Fabs') ) == u('|F|')
    assert upretty( Symbol('Fmag') ) == u('|F|')
    # Combinations
    assert upretty( Symbol('xvecdot') ) == u('x⃗̇')
    assert upretty( Symbol('xDotVec') ) == u('ẋ⃗')
    assert upretty( Symbol('xHATNorm') ) == u('‖x̂‖')
    assert upretty( Symbol('xMathring_yCheckPRM__zbreveAbs') ) == u('x̊_y̌ ̍__|z̆|')
    assert upretty( Symbol('alphadothat_nVECDOT__tTildePrime') ) == u('α̇̂_n⃗̇__t̃ ̍')
    assert upretty( Symbol('x_dot') ) == u('x_dot')
    assert upretty( Symbol('x__dot') ) == u('x__dot')


def test_pretty_basic():
    assert pretty( -Rational(1)/2 ) == '-1/2'
    assert pretty( -Rational(13)/22 ) == \
"""\
-13 \n\
----\n\
 22 \
"""
    expr = oo
    ascii_str = \
"""\
oo\
"""
    ucode_str = \
u("""\
∞\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = (x**2)
    ascii_str = \
"""\
 2\n\
x \
"""
    ucode_str = \
u("""\
 2\n\
x \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = 1/x
    ascii_str = \
"""\
1\n\
-\n\
x\
"""
    ucode_str = \
u("""\
1\n\
─\n\
x\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    # not the same as 1/x
    expr = x**-1.0
    ascii_str = \
"""\
 -1.0\n\
x    \
"""
    ucode_str = \
("""\
 -1.0\n\
x    \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    # see issue #2860
    expr = S(2)**-1.0
    ascii_str = \
"""\
 -1.0\n\
2    \
"""
    ucode_str = \
("""\
 -1.0\n\
2    \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = y*x**-2
    ascii_str = \
"""\
y \n\
--\n\
 2\n\
x \
"""
    ucode_str = \
u("""\
y \n\
──\n\
 2\n\
x \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = x**Rational(-5, 2)
    ascii_str = \
"""\
 1  \n\
----\n\
 5/2\n\
x   \
"""
    ucode_str = \
u("""\
 1  \n\
────\n\
 5/2\n\
x   \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = (-2)**x
    ascii_str = \
"""\
    x\n\
(-2) \
"""
    ucode_str = \
u("""\
    x\n\
(-2) \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    # See issue 4923
    expr = Pow(3, 1, evaluate=False)
    ascii_str = \
"""\
 1\n\
3 \
"""
    ucode_str = \
u("""\
 1\n\
3 \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = (x**2 + x + 1)
    ascii_str_1 = \
"""\
         2\n\
1 + x + x \
"""
    ascii_str_2 = \
"""\
 2        \n\
x  + x + 1\
"""
    ascii_str_3 = \
"""\
 2        \n\
x  + 1 + x\
"""
    ucode_str_1 = \
u("""\
         2\n\
1 + x + x \
""")
    ucode_str_2 = \
u("""\
 2        \n\
x  + x + 1\
""")
    ucode_str_3 = \
u("""\
 2        \n\
x  + 1 + x\
""")
    assert pretty(expr) in [ascii_str_1, ascii_str_2, ascii_str_3]
    assert upretty(expr) in [ucode_str_1, ucode_str_2, ucode_str_3]

    expr = 1 - x
    ascii_str_1 = \
"""\
1 - x\
"""
    ascii_str_2 = \
"""\
-x + 1\
"""
    ucode_str_1 = \
u("""\
1 - x\
""")
    ucode_str_2 = \
u("""\
-x + 1\
""")
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = 1 - 2*x
    ascii_str_1 = \
"""\
1 - 2*x\
"""
    ascii_str_2 = \
"""\
-2*x + 1\
"""
    ucode_str_1 = \
u("""\
1 - 2⋅x\
""")
    ucode_str_2 = \
u("""\
-2⋅x + 1\
""")
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = x/y
    ascii_str = \
"""\
x\n\
-\n\
y\
"""
    ucode_str = \
u("""\
x\n\
─\n\
y\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = -x/y
    ascii_str = \
"""\
-x \n\
---\n\
 y \
"""
    ucode_str = \
u("""\
-x \n\
───\n\
 y \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = (x + 2)/y
    ascii_str_1 = \
"""\
2 + x\n\
-----\n\
  y  \
"""
    ascii_str_2 = \
"""\
x + 2\n\
-----\n\
  y  \
"""
    ucode_str_1 = \
u("""\
2 + x\n\
─────\n\
  y  \
""")
    ucode_str_2 = \
u("""\
x + 2\n\
─────\n\
  y  \
""")
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = (1 + x)*y
    ascii_str_1 = \
"""\
y*(1 + x)\
"""
    ascii_str_2 = \
"""\
(1 + x)*y\
"""
    ascii_str_3 = \
"""\
y*(x + 1)\
"""
    ucode_str_1 = \
u("""\
y⋅(1 + x)\
""")
    ucode_str_2 = \
u("""\
(1 + x)⋅y\
""")
    ucode_str_3 = \
u("""\
y⋅(x + 1)\
""")
    assert pretty(expr) in [ascii_str_1, ascii_str_2, ascii_str_3]
    assert upretty(expr) in [ucode_str_1, ucode_str_2, ucode_str_3]

    # Test for correct placement of the negative sign
    expr = -5*x/(x + 10)
    ascii_str_1 = \
"""\
-5*x  \n\
------\n\
10 + x\
"""
    ascii_str_2 = \
"""\
-5*x  \n\
------\n\
x + 10\
"""
    ucode_str_1 = \
u("""\
-5⋅x  \n\
──────\n\
10 + x\
""")
    ucode_str_2 = \
u("""\
-5⋅x  \n\
──────\n\
x + 10\
""")
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = -S(1)/2 - 3*x
    ascii_str = \
"""\
-3*x - 1/2\
"""
    ucode_str = \
u("""\
-3⋅x - 1/2\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = S(1)/2 - 3*x
    ascii_str = \
"""\
-3*x + 1/2\
"""
    ucode_str = \
u("""\
-3⋅x + 1/2\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = -S(1)/2 - 3*x/2
    ascii_str = \
"""\
  3*x   1\n\
- --- - -\n\
   2    2\
"""
    ucode_str = \
u("""\
  3⋅x   1\n\
- ─── - ─\n\
   2    2\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = S(1)/2 - 3*x/2
    ascii_str = \
"""\
  3*x   1\n\
- --- + -\n\
   2    2\
"""
    ucode_str = \
u("""\
  3⋅x   1\n\
- ─── + ─\n\
   2    2\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_negative_fractions():
    expr = -x/y
    ascii_str =\
"""\
-x \n\
---\n\
 y \
"""
    ucode_str =\
u("""\
-x \n\
───\n\
 y \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str
    expr = -x*z/y
    ascii_str =\
"""\
-x*z \n\
-----\n\
  y  \
"""
    ucode_str =\
u("""\
-x⋅z \n\
─────\n\
  y  \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str
    expr = x**2/y
    ascii_str =\
"""\
 2\n\
x \n\
--\n\
y \
"""
    ucode_str =\
u("""\
 2\n\
x \n\
──\n\
y \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str
    expr = -x**2/y
    ascii_str =\
"""\
  2 \n\
-x  \n\
----\n\
 y  \
"""
    ucode_str =\
u("""\
  2 \n\
-x  \n\
────\n\
 y  \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str
    expr = -x/(y*z)
    ascii_str =\
"""\
-x \n\
---\n\
y*z\
"""
    ucode_str =\
u("""\
-x \n\
───\n\
y⋅z\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str
    expr = -a/y**2
    ascii_str =\
"""\
-a \n\
---\n\
  2\n\
 y \
"""
    ucode_str =\
u("""\
-a \n\
───\n\
  2\n\
 y \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str
    expr = y**(-a/b)
    ascii_str =\
"""\
 -a \n\
 ---\n\
  b \n\
y   \
"""
    ucode_str =\
u("""\
 -a \n\
 ───\n\
  b \n\
y   \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str
    expr = -1/y**2
    ascii_str =\
"""\
-1 \n\
---\n\
  2\n\
 y \
"""
    ucode_str =\
u("""\
-1 \n\
───\n\
  2\n\
 y \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str
    expr = -10/b**2
    ascii_str =\
"""\
-10 \n\
----\n\
  2 \n\
 b  \
"""
    ucode_str =\
u("""\
-10 \n\
────\n\
  2 \n\
 b  \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str
    expr = Rational(-200, 37)
    ascii_str =\
"""\
-200 \n\
-----\n\
  37 \
"""
    ucode_str =\
u("""\
-200 \n\
─────\n\
  37 \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

def test_issue_5524():
    assert pretty(-(-x + 5)*(-x - 2*sqrt(2) + 5) - (-y + 5)*(-y + 5)) == \
"""\
        /         ___    \\           2\n\
(x - 5)*\\-x - 2*\\/ 2  + 5/ - (-y + 5) \
"""

    assert upretty(-(-x + 5)*(-x - 2*sqrt(2) + 5) - (-y + 5)*(-y + 5)) == \
u("""\
        ⎛         ___    ⎞           2\n\
(x - 5)⋅⎝-x - 2⋅╲╱ 2  + 5⎠ - (-y + 5) \
""")


def test_pretty_ordering():
    assert pretty(x**2 + x + 1, order='lex') == \
"""\
 2        \n\
x  + x + 1\
"""
    assert pretty(x**2 + x + 1, order='rev-lex') == \
"""\
         2\n\
1 + x + x \
"""
    assert pretty(1 - x, order='lex') == '-x + 1'
    assert pretty(1 - x, order='rev-lex') == '1 - x'

    assert pretty(1 - 2*x, order='lex') == '-2*x + 1'
    assert pretty(1 - 2*x, order='rev-lex') == '1 - 2*x'

    f = 2*x**4 + y**2 - x**2 + y**3
    assert pretty(f, order=None) == \
"""\
   4    2    3    2\n\
2*x  - x  + y  + y \
"""
    assert pretty(f, order='lex') == \
"""\
   4    2    3    2\n\
2*x  - x  + y  + y \
"""
    assert pretty(f, order='rev-lex') == \
"""\
 2    3    2      4\n\
y  + y  - x  + 2*x \
"""

    expr = x - x**3/6 + x**5/120 + O(x**6)
    ascii_str = \
"""\
     3     5        \n\
    x     x     / 6\\\n\
x - -- + --- + O\\x /\n\
    6    120        \
"""
    ucode_str = \
u("""\
     3     5        \n\
    x     x     ⎛ 6⎞\n\
x - ── + ─── + O⎝x ⎠\n\
    6    120        \
""")
    assert pretty(expr, order=None) == ascii_str
    assert upretty(expr, order=None) == ucode_str

    assert pretty(expr, order='lex') == ascii_str
    assert upretty(expr, order='lex') == ucode_str

    assert pretty(expr, order='rev-lex') == ascii_str
    assert upretty(expr, order='rev-lex') == ucode_str


def test_pretty_relational():
    expr = Eq(x, y)
    ascii_str = \
"""\
x = y\
"""
    ucode_str = \
u("""\
x = y\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Lt(x, y)
    ascii_str = \
"""\
x < y\
"""
    ucode_str = \
u("""\
x < y\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Gt(x, y)
    ascii_str = \
"""\
x > y\
"""
    ucode_str = \
u("""\
x > y\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Le(x, y)
    ascii_str = \
"""\
x <= y\
"""
    ucode_str = \
u("""\
x ≤ y\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Ge(x, y)
    ascii_str = \
"""\
x >= y\
"""
    ucode_str = \
u("""\
x ≥ y\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Ne(x/(y + 1), y**2)
    ascii_str_1 = \
"""\
  x       2\n\
----- != y \n\
1 + y      \
"""
    ascii_str_2 = \
"""\
  x       2\n\
----- != y \n\
y + 1      \
"""
    ucode_str_1 = \
u("""\
  x      2\n\
───── ≠ y \n\
1 + y     \
""")
    ucode_str_2 = \
u("""\
  x      2\n\
───── ≠ y \n\
y + 1     \
""")
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]


def test_issue_7117():
    # See also issue #5031 (hence the evaluate=False in these).
    e = Eq(x + 1, x/2)
    q = Mul(2, e, evaluate=False)
    assert upretty(q) == u("""\
  ⎛        x⎞\n\
2⋅⎜x + 1 = ─⎟\n\
  ⎝        2⎠\
""")
    q = Add(e, 6, evaluate=False)
    assert upretty(q) == u("""\
    ⎛        x⎞\n\
6 + ⎜x + 1 = ─⎟\n\
    ⎝        2⎠\
""")
    q = Pow(e, 2, evaluate=False)
    assert upretty(q) == u("""\
           2\n\
⎛        x⎞ \n\
⎜x + 1 = ─⎟ \n\
⎝        2⎠ \
""")
    e2 = Eq(x, 2)
    q = Mul(e, e2, evaluate=False)
    assert upretty(q) == u("""\
⎛        x⎞        \n\
⎜x + 1 = ─⎟⋅(x = 2)\n\
⎝        2⎠        \
""")


def test_pretty_rational():
    expr = y*x**-2
    ascii_str = \
"""\
y \n\
--\n\
 2\n\
x \
"""
    ucode_str = \
u("""\
y \n\
──\n\
 2\n\
x \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = y**Rational(3, 2) * x**Rational(-5, 2)
    ascii_str = \
"""\
 3/2\n\
y   \n\
----\n\
 5/2\n\
x   \
"""
    ucode_str = \
u("""\
 3/2\n\
y   \n\
────\n\
 5/2\n\
x   \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = sin(x)**3/tan(x)**2
    ascii_str = \
"""\
   3   \n\
sin (x)\n\
-------\n\
   2   \n\
tan (x)\
"""
    ucode_str = \
u("""\
   3   \n\
sin (x)\n\
───────\n\
   2   \n\
tan (x)\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_functions():
    """Tests for Abs, conjugate, exp, function braces, and factorial."""
    expr = (2*x + exp(x))
    ascii_str_1 = \
"""\
       x\n\
2*x + e \
"""
    ascii_str_2 = \
"""\
 x      \n\
e  + 2*x\
"""
    ucode_str_1 = \
u("""\
       x\n\
2⋅x + ℯ \
""")
    ucode_str_2 = \
u("""\
 x     \n\
ℯ + 2⋅x\
""")
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = Abs(x)
    ascii_str = \
"""\
|x|\
"""
    ucode_str = \
u("""\
│x│\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Abs(x/(x**2 + 1))
    ascii_str_1 = \
"""\
|  x   |\n\
|------|\n\
|     2|\n\
|1 + x |\
"""
    ascii_str_2 = \
"""\
|  x   |\n\
|------|\n\
| 2    |\n\
|x  + 1|\
"""
    ucode_str_1 = \
u("""\
│  x   │\n\
│──────│\n\
│     2│\n\
│1 + x │\
""")
    ucode_str_2 = \
u("""\
│  x   │\n\
│──────│\n\
│ 2    │\n\
│x  + 1│\
""")
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = Abs(1 / (y - Abs(x)))
    ascii_str = \
"""\
|   1   |\n\
|-------|\n\
|y - |x||\
"""
    ucode_str = \
u("""\
│   1   │\n\
│───────│\n\
│y - │x││\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    n = Symbol('n', integer=True)
    expr = factorial(n)
    ascii_str = \
"""\
n!\
"""
    ucode_str = \
u("""\
n!\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = factorial(2*n)
    ascii_str = \
"""\
(2*n)!\
"""
    ucode_str = \
u("""\
(2⋅n)!\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = factorial(factorial(factorial(n)))
    ascii_str = \
"""\
((n!)!)!\
"""
    ucode_str = \
u("""\
((n!)!)!\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = factorial(n + 1)
    ascii_str_1 = \
"""\
(1 + n)!\
"""
    ascii_str_2 = \
"""\
(n + 1)!\
"""
    ucode_str_1 = \
u("""\
(1 + n)!\
""")
    ucode_str_2 = \
u("""\
(n + 1)!\
""")

    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = subfactorial(n)
    ascii_str = \
"""\
!n\
"""
    ucode_str = \
u("""\
!n\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = subfactorial(2*n)
    ascii_str = \
"""\
!(2*n)\
"""
    ucode_str = \
u("""\
!(2⋅n)\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    n = Symbol('n', integer=True)
    expr = factorial2(n)
    ascii_str = \
"""\
n!!\
"""
    ucode_str = \
u("""\
n!!\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = factorial2(2*n)
    ascii_str = \
"""\
(2*n)!!\
"""
    ucode_str = \
u("""\
(2⋅n)!!\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = factorial2(factorial2(factorial2(n)))
    ascii_str = \
"""\
((n!!)!!)!!\
"""
    ucode_str = \
u("""\
((n!!)!!)!!\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = factorial2(n + 1)
    ascii_str_1 = \
"""\
(1 + n)!!\
"""
    ascii_str_2 = \
"""\
(n + 1)!!\
"""
    ucode_str_1 = \
u("""\
(1 + n)!!\
""")
    ucode_str_2 = \
u("""\
(n + 1)!!\
""")

    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = 2*binomial(n, k)
    ascii_str = \
"""\
  /n\\\n\
2*| |\n\
  \k/\
"""
    ucode_str = \
u("""\
  ⎛n⎞\n\
2⋅⎜ ⎟\n\
  ⎝k⎠\
""")

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = 2*binomial(2*n, k)
    ascii_str = \
"""\
  /2*n\\\n\
2*|   |\n\
  \ k /\
"""
    ucode_str = \
u("""\
  ⎛2⋅n⎞\n\
2⋅⎜   ⎟\n\
  ⎝ k ⎠\
""")

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = 2*binomial(n**2, k)
    ascii_str = \
"""\
  / 2\\\n\
  |n |\n\
2*|  |\n\
  \k /\
"""
    ucode_str = \
u("""\
  ⎛ 2⎞\n\
  ⎜n ⎟\n\
2⋅⎜  ⎟\n\
  ⎝k ⎠\
""")

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = catalan(n)
    ascii_str = \
"""\
C \n\
 n\
"""
    ucode_str = \
u("""\
C \n\
 n\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = conjugate(x)
    ascii_str = \
"""\
_\n\
x\
"""
    ucode_str = \
u("""\
_\n\
x\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    f = Function('f')
    expr = conjugate(f(x + 1))
    ascii_str_1 = \
"""\
________\n\
f(1 + x)\
"""
    ascii_str_2 = \
"""\
________\n\
f(x + 1)\
"""
    ucode_str_1 = \
u("""\
________\n\
f(1 + x)\
""")
    ucode_str_2 = \
u("""\
________\n\
f(x + 1)\
""")
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = f(x)
    ascii_str = \
"""\
f(x)\
"""
    ucode_str = \
u("""\
f(x)\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = f(x, y)
    ascii_str = \
"""\
f(x, y)\
"""
    ucode_str = \
u("""\
f(x, y)\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = f(x/(y + 1), y)
    ascii_str_1 = \
"""\
 /  x     \\\n\
f|-----, y|\n\
 \\1 + y   /\
"""
    ascii_str_2 = \
"""\
 /  x     \\\n\
f|-----, y|\n\
 \\y + 1   /\
"""
    ucode_str_1 = \
u("""\
 ⎛  x     ⎞\n\
f⎜─────, y⎟\n\
 ⎝1 + y   ⎠\
""")
    ucode_str_2 = \
u("""\
 ⎛  x     ⎞\n\
f⎜─────, y⎟\n\
 ⎝y + 1   ⎠\
""")
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = sin(x)**2
    ascii_str = \
"""\
   2   \n\
sin (x)\
"""
    ucode_str = \
u("""\
   2   \n\
sin (x)\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = conjugate(a + b*I)
    ascii_str = \
"""\
_     _\n\
a - I*b\
"""
    ucode_str = \
u("""\
_     _\n\
a - ⅈ⋅b\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = conjugate(exp(a + b*I))
    ascii_str = \
"""\
 _     _\n\
 a - I*b\n\
e       \
"""
    ucode_str = \
u("""\
 _     _\n\
 a - ⅈ⋅b\n\
ℯ       \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = conjugate( f(1 + conjugate(f(x))) )
    ascii_str_1 = \
"""\
___________\n\
 /    ____\\\n\
f\\1 + f(x)/\
"""
    ascii_str_2 = \
"""\
___________\n\
 /____    \\\n\
f\\f(x) + 1/\
"""
    ucode_str_1 = \
u("""\
___________\n\
 ⎛    ____⎞\n\
f⎝1 + f(x)⎠\
""")
    ucode_str_2 = \
u("""\
___________\n\
 ⎛____    ⎞\n\
f⎝f(x) + 1⎠\
""")
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = f(x/(y + 1), y)
    ascii_str_1 = \
"""\
 /  x     \\\n\
f|-----, y|\n\
 \\1 + y   /\
"""
    ascii_str_2 = \
"""\
 /  x     \\\n\
f|-----, y|\n\
 \\y + 1   /\
"""
    ucode_str_1 = \
u("""\
 ⎛  x     ⎞\n\
f⎜─────, y⎟\n\
 ⎝1 + y   ⎠\
""")
    ucode_str_2 = \
u("""\
 ⎛  x     ⎞\n\
f⎜─────, y⎟\n\
 ⎝y + 1   ⎠\
""")
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = floor(1 / (y - floor(x)))
    ascii_str = \
"""\
     /     1      \\\n\
floor|------------|\n\
     \y - floor(x)/\
"""
    ucode_str = \
u("""\
⎢   1   ⎥\n\
⎢───────⎥\n\
⎣y - ⌊x⌋⎦\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = ceiling(1 / (y - ceiling(x)))
    ascii_str = \
"""\
       /      1       \\\n\
ceiling|--------------|\n\
       \y - ceiling(x)/\
"""
    ucode_str = \
u("""\
⎡   1   ⎤\n\
⎢───────⎥\n\
⎢y - ⌈x⌉⎥\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = euler(n)
    ascii_str = \
"""\
E \n\
 n\
"""
    ucode_str = \
u("""\
E \n\
 n\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = euler(1/(1 + 1/(1 + 1/n)))
    ascii_str = \
"""\
E         \n\
     1    \n\
 ---------\n\
       1  \n\
 1 + -----\n\
         1\n\
     1 + -\n\
         n\
"""

    ucode_str = \
u("""\
E         \n\
     1    \n\
 ─────────\n\
       1  \n\
 1 + ─────\n\
         1\n\
     1 + ─\n\
         n\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_sqrt():
    expr = sqrt(2)
    ascii_str = \
"""\
  ___\n\
\/ 2 \
"""
    ucode_str = \
u("""\
  ___\n\
╲╱ 2 \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = 2**Rational(1, 3)
    ascii_str = \
"""\
3 ___\n\
\/ 2 \
"""
    ucode_str = \
u("""\
3 ___\n\
╲╱ 2 \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = 2**Rational(1, 1000)
    ascii_str = \
"""\
1000___\n\
  \/ 2 \
"""
    ucode_str = \
u("""\
1000___\n\
  ╲╱ 2 \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = sqrt(x**2 + 1)
    ascii_str = \
"""\
   ________\n\
  /  2     \n\
\/  x  + 1 \
"""
    ucode_str = \
u("""\
   ________\n\
  ╱  2     \n\
╲╱  x  + 1 \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = (1 + sqrt(5))**Rational(1, 3)
    ascii_str = \
"""\
   ___________\n\
3 /       ___ \n\
\/  1 + \/ 5  \
"""
    ucode_str = \
u("""\
   ___________\n\
3 ╱       ___ \n\
╲╱  1 + ╲╱ 5  \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = 2**(1/x)
    ascii_str = \
"""\
x ___\n\
\/ 2 \
"""
    ucode_str = \
u("""\
x ___\n\
╲╱ 2 \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = sqrt(2 + pi)
    ascii_str = \
"""\
  ________\n\
\/ 2 + pi \
"""
    ucode_str = \
u("""\
  _______\n\
╲╱ 2 + π \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = (2 + (
        1 + x**2)/(2 + x))**Rational(1, 4) + (1 + x**Rational(1, 1000))/sqrt(3 + x**2)
    ascii_str = \
"""\
     ____________              \n\
    /      2        1000___    \n\
   /      x  + 1      \/ x  + 1\n\
4 /   2 + ------  + -----------\n\
\/        x + 2        ________\n\
                      /  2     \n\
                    \/  x  + 3 \
"""
    ucode_str = \
u("""\
     ____________              \n\
    ╱      2        1000___    \n\
   ╱      x  + 1      ╲╱ x  + 1\n\
4 ╱   2 + ──────  + ───────────\n\
╲╱        x + 2        ________\n\
                      ╱  2     \n\
                    ╲╱  x  + 3 \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_KroneckerDelta():
    x, y = symbols("x, y")
    expr = KroneckerDelta(x, y)
    ascii_str = \
"""\
d   \n\
 x,y\
"""
    ucode_str = \
u("""\
δ   \n\
 x,y\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_product():
    n, m, k, l = symbols('n m k l')
    f = symbols('f', cls=Function)
    expr = Product(f((n/3)**2), (n, k**2, l))

    unicode_str = \
u("""\
    l           \n\
┬────────┬      \n\
│        │  ⎛ 2⎞\n\
│        │  ⎜n ⎟\n\
│        │ f⎜──⎟\n\
│        │  ⎝9 ⎠\n\
│        │      \n\
       2        \n\
  n = k         """)
    ascii_str = \
"""\
    l           \n\
__________      \n\
|        |  / 2\\\n\
|        |  |n |\n\
|        | f|--|\n\
|        |  \9 /\n\
|        |      \n\
       2        \n\
  n = k         """

    assert pretty(expr) == ascii_str
    assert upretty(expr) == unicode_str

    expr = Product(f((n/3)**2), (n, k**2, l), (l, 1, m))

    unicode_str = \
u("""\
    m          l           \n\
┬────────┬ ┬────────┬      \n\
│        │ │        │  ⎛ 2⎞\n\
│        │ │        │  ⎜n ⎟\n\
│        │ │        │ f⎜──⎟\n\
│        │ │        │  ⎝9 ⎠\n\
│        │ │        │      \n\
  l = 1           2        \n\
             n = k         """)
    ascii_str = \
"""\
    m          l           \n\
__________ __________      \n\
|        | |        |  / 2\\\n\
|        | |        |  |n |\n\
|        | |        | f|--|\n\
|        | |        |  \9 /\n\
|        | |        |      \n\
  l = 1           2        \n\
             n = k         """

    assert pretty(expr) == ascii_str
    assert upretty(expr) == unicode_str


def test_pretty_lambda():
    # S.IdentityFunction is a special case
    expr = Lambda(y, y)
    assert pretty(expr) == "x -> x"
    assert upretty(expr) == u("x ↦ x")

    expr = Lambda(x, x+1)
    assert pretty(expr) == "x -> x + 1"
    assert upretty(expr) == u("x ↦ x + 1")

    expr = Lambda(x, x**2)
    ascii_str = \
"""\
      2\n\
x -> x \
"""
    ucode_str = \
u("""\
     2\n\
x ↦ x \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Lambda(x, x**2)**2
    ascii_str = \
"""\
         2
/      2\\ \n\
\\x -> x / \
"""
    ucode_str = \
u("""\
        2
⎛     2⎞ \n\
⎝x ↦ x ⎠ \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Lambda((x, y), x)
    ascii_str = "(x, y) -> x"
    ucode_str = u("(x, y) ↦ x")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Lambda((x, y), x**2)
    ascii_str = \
"""\
           2\n\
(x, y) -> x \
"""
    ucode_str = \
u("""\
          2\n\
(x, y) ↦ x \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_order():
    expr = O(1)
    ascii_str = \
"""\
O(1)\
"""
    ucode_str = \
u("""\
O(1)\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = O(1/x)
    ascii_str = \
"""\
 /1\\\n\
O|-|\n\
 \\x/\
"""
    ucode_str = \
u("""\
 ⎛1⎞\n\
O⎜─⎟\n\
 ⎝x⎠\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = O(x**2 + y**2)
    ascii_str = \
"""\
 / 2    2                  \\\n\
O\\x  + y ; (x, y) -> (0, 0)/\
"""
    ucode_str = \
u("""\
 ⎛ 2    2                 ⎞\n\
O⎝x  + y ; (x, y) → (0, 0)⎠\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = O(1, (x, oo))
    ascii_str = \
"""\
O(1; x -> oo)\
"""
    ucode_str = \
u("""\
O(1; x → ∞)\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = O(1/x, (x, oo))
    ascii_str = \
"""\
 /1         \\\n\
O|-; x -> oo|\n\
 \\x         /\
"""
    ucode_str = \
u("""\
 ⎛1       ⎞\n\
O⎜─; x → ∞⎟\n\
 ⎝x       ⎠\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = O(x**2 + y**2, (x, oo), (y, oo))
    ascii_str = \
"""\
 / 2    2                    \\\n\
O\\x  + y ; (x, y) -> (oo, oo)/\
"""
    ucode_str = \
u("""\
 ⎛ 2    2                 ⎞\n\
O⎝x  + y ; (x, y) → (∞, ∞)⎠\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_derivatives():
    # Simple
    expr = Derivative(log(x), x, evaluate=False)
    ascii_str = \
"""\
d         \n\
--(log(x))\n\
dx        \
"""
    ucode_str = \
u("""\
d         \n\
──(log(x))\n\
dx        \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Derivative(log(x), x, evaluate=False) + x
    ascii_str_1 = \
"""\
    d         \n\
x + --(log(x))\n\
    dx        \
"""
    ascii_str_2 = \
"""\
d             \n\
--(log(x)) + x\n\
dx            \
"""
    ucode_str_1 = \
u("""\
    d         \n\
x + ──(log(x))\n\
    dx        \
""")
    ucode_str_2 = \
u("""\
d             \n\
──(log(x)) + x\n\
dx            \
""")
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    # basic partial derivatives
    expr = Derivative(log(x + y) + x, x)
    ascii_str_1 = \
"""\
d                 \n\
--(log(x + y) + x)\n\
dx                \
"""
    ascii_str_2 = \
"""\
d                 \n\
--(x + log(x + y))\n\
dx                \
"""
    ucode_str_1 = \
u("""\
∂                 \n\
──(log(x + y) + x)\n\
∂x                \
""")
    ucode_str_2 = \
u("""\
∂                 \n\
──(x + log(x + y))\n\
∂x                \
""")
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2], upretty(expr)

    # Multiple symbols
    expr = Derivative(log(x) + x**2, x, y)
    ascii_str_1 = \
"""\
   2              \n\
  d  /          2\\\n\
-----\log(x) + x /\n\
dy dx             \
"""
    ascii_str_2 = \
"""\
   2              \n\
  d  / 2         \\\n\
-----\\x  + log(x)/\n\
dy dx             \
"""
    ucode_str_1 = \
u("""\
   2              \n\
  d  ⎛          2⎞\n\
─────⎝log(x) + x ⎠\n\
dy dx             \
""")
    ucode_str_2 = \
u("""\
   2              \n\
  d  ⎛ 2         ⎞\n\
─────⎝x  + log(x)⎠\n\
dy dx             \
""")
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = Derivative(2*x*y, y, x) + x**2
    ascii_str_1 = \
"""\
   2             \n\
  d             2\n\
-----(2*x*y) + x \n\
dx dy            \
"""
    ascii_str_2 = \
"""\
        2        \n\
 2     d         \n\
x  + -----(2*x*y)\n\
     dx dy       \
"""
    ucode_str_1 = \
u("""\
   2             \n\
  ∂             2\n\
─────(2⋅x⋅y) + x \n\
∂x ∂y            \
""")
    ucode_str_2 = \
u("""\
        2        \n\
 2     ∂         \n\
x  + ─────(2⋅x⋅y)\n\
     ∂x ∂y       \
""")
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = Derivative(2*x*y, x, x)
    ascii_str = \
"""\
  2       \n\
 d        \n\
---(2*x*y)\n\
  2       \n\
dx        \
"""
    ucode_str = \
u("""\
  2       \n\
 ∂        \n\
───(2⋅x⋅y)\n\
  2       \n\
∂x        \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Derivative(2*x*y, x, 17)
    ascii_str = \
"""\
 17        \n\
d          \n\
----(2*x*y)\n\
  17       \n\
dx         \
"""
    ucode_str = \
u("""\
 17        \n\
∂          \n\
────(2⋅x⋅y)\n\
  17       \n\
∂x         \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Derivative(2*x*y, x, x, y)
    ascii_str = \
"""\
   3         \n\
  d          \n\
------(2*x*y)\n\
     2       \n\
dy dx        \
"""
    ucode_str = \
u("""\
   3         \n\
  ∂          \n\
──────(2⋅x⋅y)\n\
     2       \n\
∂y ∂x        \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    # Greek letters
    alpha = Symbol('alpha')
    beta = Function('beta')
    expr = beta(alpha).diff(alpha)
    ascii_str = \
"""\
  d                \n\
------(beta(alpha))\n\
dalpha             \
"""
    ucode_str = \
u("""\
d       \n\
──(β(α))\n\
dα      \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_integrals():
    expr = Integral(log(x), x)
    ascii_str = \
"""\
  /         \n\
 |          \n\
 | log(x) dx\n\
 |          \n\
/           \
"""
    ucode_str = \
u("""\
⌠          \n\
⎮ log(x) dx\n\
⌡          \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Integral(x**2, x)
    ascii_str = \
"""\
  /     \n\
 |      \n\
 |  2   \n\
 | x  dx\n\
 |      \n\
/       \
"""
    ucode_str = \
u("""\
⌠      \n\
⎮  2   \n\
⎮ x  dx\n\
⌡      \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Integral((sin(x))**2 / (tan(x))**2)
    ascii_str = \
"""\
  /          \n\
 |           \n\
 |    2      \n\
 | sin (x)   \n\
 | ------- dx\n\
 |    2      \n\
 | tan (x)   \n\
 |           \n\
/            \
"""
    ucode_str = \
u("""\
⌠           \n\
⎮    2      \n\
⎮ sin (x)   \n\
⎮ ─────── dx\n\
⎮    2      \n\
⎮ tan (x)   \n\
⌡           \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Integral(x**(2**x), x)
    ascii_str = \
"""\
  /        \n\
 |         \n\
 |  / x\   \n\
 |  \\2 /   \n\
 | x     dx\n\
 |         \n\
/          \
"""
    ucode_str = \
u("""\
⌠         \n\
⎮  ⎛ x⎞   \n\
⎮  ⎝2 ⎠   \n\
⎮ x     dx\n\
⌡         \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Integral(x**2, (x, 1, 2))
    ascii_str = \
"""\
  2      \n\
  /      \n\
 |       \n\
 |   2   \n\
 |  x  dx\n\
 |       \n\
/        \n\
1        \
"""
    ucode_str = \
u("""\
2      \n\
⌠      \n\
⎮  2   \n\
⎮ x  dx\n\
⌡      \n\
1      \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Integral(x**2, (x, Rational(1, 2), 10))
    ascii_str = \
"""\
 10      \n\
  /      \n\
 |       \n\
 |   2   \n\
 |  x  dx\n\
 |       \n\
/        \n\
1/2      \
"""
    ucode_str = \
u("""\
 10      \n\
 ⌠       \n\
 ⎮   2   \n\
 ⎮  x  dx\n\
 ⌡       \n\
1/2      \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Integral(x**2*y**2, x, y)
    ascii_str = \
"""\
  /  /           \n\
 |  |            \n\
 |  |  2  2      \n\
 |  | x *y  dx dy\n\
 |  |            \n\
/  /             \
"""
    ucode_str = \
u("""\
⌠ ⌠            \n\
⎮ ⎮  2  2      \n\
⎮ ⎮ x ⋅y  dx dy\n\
⌡ ⌡            \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Integral(sin(th)/cos(ph), (th, 0, pi), (ph, 0, 2*pi))
    ascii_str = \
"""\
 2*pi pi                           \n\
   /   /                           \n\
  |   |                            \n\
  |   |  sin(theta)                \n\
  |   |  ---------- d(theta) d(phi)\n\
  |   |   cos(phi)                 \n\
  |   |                            \n\
 /   /                             \n\
 0   0                             \
"""
    ucode_str = \
u("""\
2⋅π π             \n\
 ⌠  ⌠             \n\
 ⎮  ⎮ sin(θ)      \n\
 ⎮  ⎮ ────── dθ dφ\n\
 ⎮  ⎮ cos(φ)      \n\
 ⌡  ⌡             \n\
 0  0             \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_matrix():
    # Empty Matrix
    expr = Matrix()
    ascii_str = "[]"
    unicode_str = "[]"
    assert pretty(expr) == ascii_str
    assert upretty(expr) == unicode_str
    expr = Matrix(2, 0, lambda i, j: 0)
    ascii_str = "[]"
    unicode_str = "[]"
    assert pretty(expr) == ascii_str
    assert upretty(expr) == unicode_str
    expr = Matrix(0, 2, lambda i, j: 0)
    ascii_str = "[]"
    unicode_str = "[]"
    assert pretty(expr) == ascii_str
    assert upretty(expr) == unicode_str
    expr = Matrix([[x**2 + 1, 1], [y, x + y]])
    ascii_str_1 = \
"""\
[     2       ]
[1 + x     1  ]
[             ]
[  y     x + y]\
"""
    ascii_str_2 = \
"""\
[ 2           ]
[x  + 1    1  ]
[             ]
[  y     x + y]\
"""
    ucode_str_1 = \
u("""\
⎡     2       ⎤
⎢1 + x     1  ⎥
⎢             ⎥
⎣  y     x + y⎦\
""")
    ucode_str_2 = \
u("""\
⎡ 2           ⎤
⎢x  + 1    1  ⎥
⎢             ⎥
⎣  y     x + y⎦\
""")
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = Matrix([[x/y, y, th], [0, exp(I*k*ph), 1]])
    ascii_str = \
"""\
[x                 ]
[-     y      theta]
[y                 ]
[                  ]
[    I*k*phi       ]
[0  e           1  ]\
"""
    ucode_str = \
u("""\
⎡x           ⎤
⎢─    y     θ⎥
⎢y           ⎥
⎢            ⎥
⎢    ⅈ⋅k⋅φ   ⎥
⎣0  ℯ       1⎦\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_Adjoint():
    from sympy.matrices import Adjoint, Inverse, MatrixSymbol, Transpose
    X = MatrixSymbol('X', 2, 2)
    Y = MatrixSymbol('Y', 2, 2)
    assert pretty(Adjoint(X)) == " +\nX "
    assert pretty(Adjoint(X + Y)) == "       +\n(X + Y) "
    assert pretty(Adjoint(X) + Adjoint(Y)) == " +    +\nX  + Y "
    assert pretty(Adjoint(X*Y)) == "     +\n(X*Y) "
    assert pretty(Adjoint(Y)*Adjoint(X)) == " +  +\nY *X "
    assert pretty(Adjoint(X**2)) == "    +\n/ 2\\ \n\\X / "
    assert pretty(Adjoint(X)**2) == "    2\n/ +\\ \n\\X / "
    assert pretty(Adjoint(Inverse(X))) == "     +\n/ -1\\ \n\\X  / "
    assert pretty(Inverse(Adjoint(X))) == "    -1\n/ +\\  \n\\X /  "
    assert pretty(Adjoint(Transpose(X))) == "    +\n/ T\\ \n\\X / "
    assert pretty(Transpose(Adjoint(X))) == "    T\n/ +\\ \n\\X / "
    assert upretty(Adjoint(X)) == u(" †\nX ")
    assert upretty(Adjoint(X + Y)) == u("       †\n(X + Y) ")
    assert upretty(Adjoint(X) + Adjoint(Y)) == u(" †    †\nX  + Y ")
    assert upretty(Adjoint(X*Y)) == u("     †\n(X⋅Y) ")
    assert upretty(Adjoint(Y)*Adjoint(X)) == u(" †  †\nY ⋅X ")
    assert upretty(Adjoint(X**2)) == \
        u("    †\n⎛ 2⎞ \n⎝X ⎠ ")
    assert upretty(Adjoint(X)**2) == \
        u("    2\n⎛ †⎞ \n⎝X ⎠ ")
    assert upretty(Adjoint(Inverse(X))) == \
        u("     †\n⎛ -1⎞ \n⎝X  ⎠ ")
    assert upretty(Inverse(Adjoint(X))) == \
        u("    -1\n⎛ †⎞  \n⎝X ⎠  ")
    assert upretty(Adjoint(Transpose(X))) == \
        u("    †\n⎛ T⎞ \n⎝X ⎠ ")
    assert upretty(Transpose(Adjoint(X))) == \
        u("    T\n⎛ †⎞ \n⎝X ⎠ ")


def test_pretty_piecewise():
    expr = Piecewise((x, x < 1), (x**2, True))
    ascii_str = \
"""\
/x   for x < 1\n\
|             \n\
< 2           \n\
|x   otherwise\n\
\             \
"""
    ucode_str = \
u("""\
⎧x   for x < 1\n\
⎪             \n\
⎨ 2           \n\
⎪x   otherwise\n\
⎩             \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = -Piecewise((x, x < 1), (x**2, True))
    ascii_str = \
"""\
 //x   for x < 1\\\n\
 ||             |\n\
-|< 2           |\n\
 ||x   otherwise|\n\
 \\\\             /\
"""
    ucode_str = \
u("""\
 ⎛⎧x   for x < 1⎞\n\
 ⎜⎪             ⎟\n\
-⎜⎨ 2           ⎟\n\
 ⎜⎪x   otherwise⎟\n\
 ⎝⎩             ⎠\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = x + Piecewise((x, x > 0), (y, True)) + Piecewise((x/y, x < 2),
    (y**2, x > 2), (1, True)) + 1
    ascii_str = \
"""\
                      //x            \    \n\
                      ||-   for x < 2|    \n\
                      ||y            |    \n\
    //x  for x > 0\   ||             |    \n\
x + |<            | + |< 2           | + 1\n\
    \\\\y  otherwise/   ||y   for x > 2|    \n\
                      ||             |    \n\
                      ||1   otherwise|    \n\
                      \\\\             /    \
"""
    ucode_str = \
u("""\
                      ⎛⎧x            ⎞    \n\
                      ⎜⎪─   for x < 2⎟    \n\
                      ⎜⎪y            ⎟    \n\
    ⎛⎧x  for x > 0⎞   ⎜⎪             ⎟    \n\
x + ⎜⎨            ⎟ + ⎜⎨ 2           ⎟ + 1\n\
    ⎝⎩y  otherwise⎠   ⎜⎪y   for x > 2⎟    \n\
                      ⎜⎪             ⎟    \n\
                      ⎜⎪1   otherwise⎟    \n\
                      ⎝⎩             ⎠    \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = x - Piecewise((x, x > 0), (y, True)) + Piecewise((x/y, x < 2),
    (y**2, x > 2), (1, True)) + 1
    ascii_str = \
"""\
                      //x            \    \n\
                      ||-   for x < 2|    \n\
                      ||y            |    \n\
    //x  for x > 0\   ||             |    \n\
x - |<            | + |< 2           | + 1\n\
    \\\\y  otherwise/   ||y   for x > 2|    \n\
                      ||             |    \n\
                      ||1   otherwise|    \n\
                      \\\\             /    \
"""
    ucode_str = \
u("""\
                      ⎛⎧x            ⎞    \n\
                      ⎜⎪─   for x < 2⎟    \n\
                      ⎜⎪y            ⎟    \n\
    ⎛⎧x  for x > 0⎞   ⎜⎪             ⎟    \n\
x - ⎜⎨            ⎟ + ⎜⎨ 2           ⎟ + 1\n\
    ⎝⎩y  otherwise⎠   ⎜⎪y   for x > 2⎟    \n\
                      ⎜⎪             ⎟    \n\
                      ⎜⎪1   otherwise⎟    \n\
                      ⎝⎩             ⎠    \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = x*Piecewise((x, x > 0), (y, True))
    ascii_str = \
"""\
  //x  for x > 0\\\n\
x*|<            |\n\
  \\\\y  otherwise/\
"""
    ucode_str = \
u("""\
  ⎛⎧x  for x > 0⎞\n\
x⋅⎜⎨            ⎟\n\
  ⎝⎩y  otherwise⎠\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Piecewise((x, x > 0), (y, True))*Piecewise((x/y, x < 2), (y**2, x >
    2), (1, True))
    ascii_str = \
"""\
                //x            \\\n\
                ||-   for x < 2|\n\
                ||y            |\n\
//x  for x > 0\ ||             |\n\
|<            |*|< 2           |\n\
\\\\y  otherwise/ ||y   for x > 2|\n\
                ||             |\n\
                ||1   otherwise|\n\
                \\\\             /\
"""
    ucode_str = \
u("""\
                ⎛⎧x            ⎞\n\
                ⎜⎪─   for x < 2⎟\n\
                ⎜⎪y            ⎟\n\
⎛⎧x  for x > 0⎞ ⎜⎪             ⎟\n\
⎜⎨            ⎟⋅⎜⎨ 2           ⎟\n\
⎝⎩y  otherwise⎠ ⎜⎪y   for x > 2⎟\n\
                ⎜⎪             ⎟\n\
                ⎜⎪1   otherwise⎟\n\
                ⎝⎩             ⎠\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = -Piecewise((x, x > 0), (y, True))*Piecewise((x/y, x < 2), (y**2, x
        > 2), (1, True))
    ascii_str = \
"""\
                 //x            \\\n\
                 ||-   for x < 2|\n\
                 ||y            |\n\
 //x  for x > 0\ ||             |\n\
-|<            |*|< 2           |\n\
 \\\\y  otherwise/ ||y   for x > 2|\n\
                 ||             |\n\
                 ||1   otherwise|\n\
                 \\\\             /\
"""
    ucode_str = \
u("""\
                 ⎛⎧x            ⎞\n\
                 ⎜⎪─   for x < 2⎟\n\
                 ⎜⎪y            ⎟\n\
 ⎛⎧x  for x > 0⎞ ⎜⎪             ⎟\n\
-⎜⎨            ⎟⋅⎜⎨ 2           ⎟\n\
 ⎝⎩y  otherwise⎠ ⎜⎪y   for x > 2⎟\n\
                 ⎜⎪             ⎟\n\
                 ⎜⎪1   otherwise⎟\n\
                 ⎝⎩             ⎠\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Piecewise((0, Abs(1/y) < 1), (1, Abs(y) < 1), (y*meijerg(((2, 1),
        ()), ((), (1, 0)), 1/y), True))
    ascii_str = \
"""\
/                                |1|    \n\
|            0               for |-| < 1\n\
|                                |y|    \n\
|                                       \n\
<            1               for |y| < 1\n\
|                                       \n\
|   __0, 2 /2, 1       | 1\             \n\
|y*/__     |           | -|   otherwise \n\
\  \\_|2, 2 \      1, 0 | y/             \
"""
    ucode_str = \
u("""\
⎧                                │1│    \n\
⎪            0               for │─│ < 1\n\
⎪                                │y│    \n\
⎪                                       \n\
⎨            1               for │y│ < 1\n\
⎪                                       \n\
⎪  ╭─╮0, 2 ⎛2, 1       │ 1⎞             \n\
⎪y⋅│╶┐     ⎜           │ ─⎟   otherwise \n\
⎩  ╰─╯2, 2 ⎝      1, 0 │ y⎠             \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    # XXX: We have to use evaluate=False here because Piecewise._eval_power
    # denests the power.
    expr = Pow(Piecewise((x, x > 0), (y, True)), 2, evaluate=False)
    ascii_str = \
"""\
               2\n\
//x  for x > 0\ \n\
|<            | \n\
\\\\y  otherwise/ \
"""
    ucode_str = \
u("""\
               2\n\
⎛⎧x  for x > 0⎞ \n\
⎜⎨            ⎟ \n\
⎝⎩y  otherwise⎠ \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

def test_pretty_seq():
    expr = ()
    ascii_str = \
"""\
()\
"""
    ucode_str = \
u("""\
()\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = []
    ascii_str = \
"""\
[]\
"""
    ucode_str = \
u("""\
[]\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = {}
    expr_2 = {}
    ascii_str = \
"""\
{}\
"""
    ucode_str = \
u("""\
{}\
""")
    assert pretty(expr) == ascii_str
    assert pretty(expr_2) == ascii_str
    assert upretty(expr) == ucode_str
    assert upretty(expr_2) == ucode_str

    expr = (1/x,)
    ascii_str = \
"""\
 1  \n\
(-,)\n\
 x  \
"""
    ucode_str = \
u("""\
⎛1 ⎞\n\
⎜─,⎟\n\
⎝x ⎠\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = [x**2, 1/x, x, y, sin(th)**2/cos(ph)**2]
    ascii_str = \
"""\
                 2        \n\
  2  1        sin (theta) \n\
[x , -, x, y, -----------]\n\
     x            2       \n\
               cos (phi)  \
"""
    ucode_str = \
u("""\
⎡                2   ⎤\n\
⎢ 2  1        sin (θ)⎥\n\
⎢x , ─, x, y, ───────⎥\n\
⎢    x           2   ⎥\n\
⎣             cos (φ)⎦\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = (x**2, 1/x, x, y, sin(th)**2/cos(ph)**2)
    ascii_str = \
"""\
                 2        \n\
  2  1        sin (theta) \n\
(x , -, x, y, -----------)\n\
     x            2       \n\
               cos (phi)  \
"""
    ucode_str = \
u("""\
⎛                2   ⎞\n\
⎜ 2  1        sin (θ)⎟\n\
⎜x , ─, x, y, ───────⎟\n\
⎜    x           2   ⎟\n\
⎝             cos (φ)⎠\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Tuple(x**2, 1/x, x, y, sin(th)**2/cos(ph)**2)
    ascii_str = \
"""\
                 2        \n\
  2  1        sin (theta) \n\
(x , -, x, y, -----------)\n\
     x            2       \n\
               cos (phi)  \
"""
    ucode_str = \
u("""\
⎛                2   ⎞\n\
⎜ 2  1        sin (θ)⎟\n\
⎜x , ─, x, y, ───────⎟\n\
⎜    x           2   ⎟\n\
⎝             cos (φ)⎠\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = {x: sin(x)}
    expr_2 = Dict({x: sin(x)})
    ascii_str = \
"""\
{x: sin(x)}\
"""
    ucode_str = \
u("""\
{x: sin(x)}\
""")
    assert pretty(expr) == ascii_str
    assert pretty(expr_2) == ascii_str
    assert upretty(expr) == ucode_str
    assert upretty(expr_2) == ucode_str

    expr = {1/x: 1/y, x: sin(x)**2}
    expr_2 = Dict({1/x: 1/y, x: sin(x)**2})
    ascii_str = \
"""\
 1  1        2    \n\
{-: -, x: sin (x)}\n\
 x  y             \
"""
    ucode_str = \
u("""\
⎧1  1        2   ⎫\n\
⎨─: ─, x: sin (x)⎬\n\
⎩x  y            ⎭\
""")
    assert pretty(expr) == ascii_str
    assert pretty(expr_2) == ascii_str
    assert upretty(expr) == ucode_str
    assert upretty(expr_2) == ucode_str

    # There used to be a bug with pretty-printing sequences of even height.
    expr = [x**2]
    ascii_str = \
"""\
  2 \n\
[x ]\
"""
    ucode_str = \
u("""\
⎡ 2⎤\n\
⎣x ⎦\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = (x**2,)
    ascii_str = \
"""\
  2  \n\
(x ,)\
"""
    ucode_str = \
u("""\
⎛ 2 ⎞\n\
⎝x ,⎠\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Tuple(x**2)
    ascii_str = \
"""\
  2  \n\
(x ,)\
"""
    ucode_str = \
u("""\
⎛ 2 ⎞\n\
⎝x ,⎠\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = {x**2: 1}
    expr_2 = Dict({x**2: 1})
    ascii_str = \
"""\
  2    \n\
{x : 1}\
"""
    ucode_str = \
u("""\
⎧ 2   ⎫\n\
⎨x : 1⎬\n\
⎩     ⎭\
""")
    assert pretty(expr) == ascii_str
    assert pretty(expr_2) == ascii_str
    assert upretty(expr) == ucode_str
    assert upretty(expr_2) == ucode_str


def test_any_object_in_sequence():
    # Cf. issue 5306
    b1 = Basic()
    b2 = Basic(Basic())

    expr = [b2, b1]
    assert pretty(expr) == "[Basic(Basic()), Basic()]"
    assert upretty(expr) == u("[Basic(Basic()), Basic()]")

    expr = set([b2, b1])
    assert pretty(expr) == "set([Basic(), Basic(Basic())])"
    assert upretty(expr) == u("set([Basic(), Basic(Basic())])")

    expr = {b2: b1, b1: b2}
    expr2 = Dict({b2: b1, b1: b2})
    assert pretty(expr) == "{Basic(): Basic(Basic()), Basic(Basic()): Basic()}"
    assert pretty(
        expr2) == "{Basic(): Basic(Basic()), Basic(Basic()): Basic()}"
    assert upretty(
        expr) == u("{Basic(): Basic(Basic()), Basic(Basic()): Basic()}")
    assert upretty(
        expr2) == u("{Basic(): Basic(Basic()), Basic(Basic()): Basic()}")


def test_pretty_sets():
    s = FiniteSet
    assert pretty(s(*[x*y, x**2])) == \
"""\
  2      \n\
{x , x*y}\
"""
    assert pretty(s(*range(1, 6))) == "{1, 2, 3, 4, 5}"
    assert pretty(s(*range(1, 13))) == "{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}"
    for s in (frozenset, set):
        assert pretty(s([x*y, x**2])) == \
"""\
%s   2       \n\
%s([x , x*y])\
""" % (" " * len(s.__name__), s.__name__)
        assert pretty(s(range(1, 6))) == "%s([1, 2, 3, 4, 5])" % s.__name__
        assert pretty(s(range(1, 13))) == \
            "%s([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])" % s.__name__

    assert pretty(Range(0, 3, 1)) == '{0, 1, 2}'

    ascii_str = '{0, 1, ..., 29}'
    ucode_str = u('{0, 1, …, 29}')
    assert pretty(Range(0, 30, 1)) == ascii_str
    assert upretty(Range(0, 30, 1)) == ucode_str

    ascii_str = '{0, 2, ..., oo}'
    ucode_str = u('{0, 2, …, ∞}')
    assert pretty(Range(0, oo, 2)) == ascii_str
    assert upretty(Range(0, oo, 2)) == ucode_str

    ascii_str = '{-oo, ..., -3, -2}'
    ucode_str = u('{-∞, …, -3, -2}')
    assert pretty(Range(-2, -oo, -1)) == ascii_str
    assert upretty(Range(-2, -oo, -1)) == ucode_str


def test_pretty_limits():
    expr = Limit(x, x, oo)
    ascii_str = \
"""\
 lim x\n\
x->oo \
"""
    ucode_str = \
u("""\
lim x\n\
x─→∞ \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Limit(x**2, x, 0)
    ascii_str = \
"""\
      2\n\
 lim x \n\
x->0+  \
"""
    ucode_str = \
u("""\
      2\n\
 lim x \n\
x─→0⁺  \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Limit(1/x, x, 0)
    ascii_str = \
"""\
     1\n\
 lim -\n\
x->0+x\
"""
    ucode_str = \
u("""\
     1\n\
 lim ─\n\
x─→0⁺x\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Limit(sin(x)/x, x, 0)
    ascii_str = \
"""\
     sin(x)\n\
 lim ------\n\
x->0+  x   \
"""
    ucode_str = \
u("""\
     sin(x)\n\
 lim ──────\n\
x─→0⁺  x   \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Limit(sin(x)/x, x, 0, "-")
    ascii_str = \
"""\
     sin(x)\n\
 lim ------\n\
x->0-  x   \
"""
    ucode_str = \
u("""\
     sin(x)\n\
 lim ──────\n\
x─→0⁻  x   \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_RootOf():
    expr = RootOf(x**5 + 11*x - 2, 0)
    ascii_str = \
"""\
      / 5              \\\n\
RootOf\\x  + 11*x - 2, 0/\
"""
    ucode_str = \
u("""\
      ⎛ 5              ⎞\n\
RootOf⎝x  + 11⋅x - 2, 0⎠\
""")

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_RootSum():
    expr = RootSum(x**5 + 11*x - 2, auto=False)
    ascii_str = \
"""\
       / 5           \\\n\
RootSum\\x  + 11*x - 2/\
"""
    ucode_str = \
u("""\
       ⎛ 5           ⎞\n\
RootSum⎝x  + 11⋅x - 2⎠\
""")

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = RootSum(x**5 + 11*x - 2, Lambda(z, exp(z)))
    ascii_str = \
"""\
       / 5                   z\\\n\
RootSum\\x  + 11*x - 2, z -> e /\
"""
    ucode_str = \
u("""\
       ⎛ 5                  z⎞\n\
RootSum⎝x  + 11⋅x - 2, z ↦ ℯ ⎠\
""")

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_GroebnerBasis():
    expr = groebner([], x, y)

    ascii_str = \
"""\
GroebnerBasis([], x, y, domain=ZZ, order=lex)\
"""
    ucode_str = \
u("""\
GroebnerBasis([], x, y, domain=ℤ, order=lex)\
""")

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    F = [x**2 - 3*y - x + 1, y**2 - 2*x + y - 1]
    expr = groebner(F, x, y, order='grlex')

    ascii_str = \
"""\
             /[ 2                 2              ]                              \\\n\
GroebnerBasis\\[x  - x - 3*y + 1, y  - 2*x + y - 1], x, y, domain=ZZ, order=grlex/\
"""
    ucode_str = \
u("""\
             ⎛⎡ 2                 2              ⎤                             ⎞\n\
GroebnerBasis⎝⎣x  - x - 3⋅y + 1, y  - 2⋅x + y - 1⎦, x, y, domain=ℤ, order=grlex⎠\
""")

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = expr.fglm('lex')

    ascii_str = \
"""\
             /[       2           4      3      2           ]                            \\\n\
GroebnerBasis\\[2*x - y  - y + 1, y  + 2*y  - 3*y  - 16*y + 7], x, y, domain=ZZ, order=lex/\
"""
    ucode_str = \
u("""\
             ⎛⎡       2           4      3      2           ⎤                           ⎞\n\
GroebnerBasis⎝⎣2⋅x - y  - y + 1, y  + 2⋅y  - 3⋅y  - 16⋅y + 7⎦, x, y, domain=ℤ, order=lex⎠\
""")

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_Boolean():
    expr = Not(x, evaluate=False)

    assert pretty(expr) == "Not(x)"
    assert upretty(expr) == u("¬x")

    expr = And(x, y)

    assert pretty(expr) == "And(x, y)"
    assert upretty(expr) == u("x ∧ y")

    expr = Or(x, y)

    assert pretty(expr) == "Or(x, y)"
    assert upretty(expr) == u("x ∨ y")

    syms = symbols('a:f')
    expr = And(*syms)

    assert pretty(expr) == "And(a, b, c, d, e, f)"
    assert upretty(expr) == u("a ∧ b ∧ c ∧ d ∧ e ∧ f")

    expr = Or(*syms)

    assert pretty(expr) == "Or(a, b, c, d, e, f)"
    assert upretty(expr) == u("a ∨ b ∨ c ∨ d ∨ e ∨ f")

    expr = Xor(x, y, evaluate=False)

    assert pretty(expr) == "Xor(x, y)"
    assert upretty(expr) == u("x ⊻ y")

    expr = Nand(x, y, evaluate=False)

    assert pretty(expr) == "Nand(x, y)"
    assert upretty(expr) == u("x ⊼ y")

    expr = Nor(x, y, evaluate=False)

    assert pretty(expr) == "Nor(x, y)"
    assert upretty(expr) == u("x ⊽ y")

    expr = Implies(x, y, evaluate=False)

    assert pretty(expr) == "Implies(x, y)"
    assert upretty(expr) == u("x → y")

    # don't sort args
    expr = Implies(y, x, evaluate=False)

    assert pretty(expr) == "Implies(y, x)"
    assert upretty(expr) == u("y → x")

    expr = Equivalent(x, y, evaluate=False)

    assert pretty(expr) == "Equivalent(x, y)"
    assert upretty(expr) == u("x ≡ y")

    expr = Equivalent(y, x, evaluate=False)

    assert pretty(expr) == "Equivalent(x, y)"
    assert upretty(expr) == u("x ≡ y")


def test_pretty_Domain():
    expr = FF(23)

    assert pretty(expr) == "GF(23)"
    assert upretty(expr) == u("ℤ₂₃")

    expr = ZZ

    assert pretty(expr) == "ZZ"
    assert upretty(expr) == u("ℤ")

    expr = QQ

    assert pretty(expr) == "QQ"
    assert upretty(expr) == u("ℚ")

    expr = RR

    assert pretty(expr) == "RR"
    assert upretty(expr) == u("ℝ")

    expr = QQ[x]

    assert pretty(expr) == "QQ[x]"
    assert upretty(expr) == u("ℚ[x]")

    expr = QQ[x, y]

    assert pretty(expr) == "QQ[x, y]"
    assert upretty(expr) == u("ℚ[x, y]")

    expr = ZZ.frac_field(x)

    assert pretty(expr) == "ZZ(x)"
    assert upretty(expr) == u("ℤ(x)")

    expr = ZZ.frac_field(x, y)

    assert pretty(expr) == "ZZ(x, y)"
    assert upretty(expr) == u("ℤ(x, y)")

    expr = QQ.poly_ring(x, y, order=grlex)

    assert pretty(expr) == "QQ[x, y, order=grlex]"
    assert upretty(expr) == u("ℚ[x, y, order=grlex]")

    expr = QQ.poly_ring(x, y, order=ilex)

    assert pretty(expr) == "QQ[x, y, order=ilex]"
    assert upretty(expr) == u("ℚ[x, y, order=ilex]")


def test_pretty_prec():
    assert xpretty(S("0.3"), full_prec=True, wrap_line=False) == "0.300000000000000"
    assert xpretty(S("0.3"), full_prec="auto", wrap_line=False) == "0.300000000000000"
    assert xpretty(S("0.3"), full_prec=False, wrap_line=False) == "0.3"
    assert xpretty(S("0.3")*x, full_prec=True, use_unicode=False, wrap_line=False) in [
        "0.300000000000000*x",
        "x*0.300000000000000"
    ]
    assert xpretty(S("0.3")*x, full_prec="auto", use_unicode=False, wrap_line=False) in [
        "0.3*x",
        "x*0.3"
    ]
    assert xpretty(S("0.3")*x, full_prec=False, use_unicode=False, wrap_line=False) in [
        "0.3*x",
        "x*0.3"
    ]


def test_pprint():
    import sys
    from sympy.core.compatibility import StringIO
    fd = StringIO()
    sso = sys.stdout
    sys.stdout = fd
    try:
        pprint(pi, use_unicode=False, wrap_line=False)
    finally:
        sys.stdout = sso
    assert fd.getvalue() == 'pi\n'


def test_pretty_class():
    """Test that the printer dispatcher correctly handles classes."""
    class C:
        pass   # C has no .__class__ and this was causing problems

    class D(object):
        pass

    assert pretty( C ) == str( C )
    assert pretty( D ) == str( D )


def test_pretty_no_wrap_line():
    huge_expr = 0
    for i in range(20):
        huge_expr += i*sin(i + x)
    assert xpretty(huge_expr            ).find('\n') != -1
    assert xpretty(huge_expr, wrap_line=False).find('\n') == -1


def test_settings():
    raises(TypeError, lambda: pretty(S(4), method="garbage"))


def test_pretty_sum():
    from sympy.abc import x, a, b, k, m, n

    expr = Sum(k**k, (k, 0, n))
    ascii_str = \
"""\
  n     \n\
 ___    \n\
 \\  `   \n\
  \\    k\n\
  /   k \n\
 /__,   \n\
k = 0   \
"""
    ucode_str = \
u("""\
  n     \n\
 ___    \n\
 ╲      \n\
  ╲    k\n\
  ╱   k \n\
 ╱      \n\
 ‾‾‾    \n\
k = 0   \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum(k**(Integral(x**n, (x, -oo, oo))), (k, 0, n**n))
    ascii_str = \
"""\
    n             \n\
   n              \n\
______            \n\
\\     `           \n\
 \\        oo      \n\
  \\        /      \n\
   \\      |       \n\
    \\     |   n   \n\
     )    |  x  dx\n\
    /     |       \n\
   /     /        \n\
  /      -oo      \n\
 /      k         \n\
/_____,           \n\
 k = 0            \
"""
    ucode_str = \
u("""\
   n            \n\
  n             \n\
______          \n\
╲               \n\
 ╲      ∞       \n\
  ╲     ⌠       \n\
   ╲    ⎮   n   \n\
    ╲   ⎮  x  dx\n\
    ╱   ⌡       \n\
   ╱    -∞      \n\
  ╱    k        \n\
 ╱              \n\
╱               \n\
‾‾‾‾‾‾          \n\
k = 0           \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum(k**(
        Integral(x**n, (x, -oo, oo))), (k, 0, Integral(x**x, (x, -oo, oo))))
    ascii_str = \
"""\
 oo                 \n\
  /                 \n\
 |                  \n\
 |   x              \n\
 |  x  dx           \n\
 |                  \n\
/                   \n\
-oo                 \n\
 ______             \n\
 \\     `            \n\
  \\         oo      \n\
   \\         /      \n\
    \\       |       \n\
     \\      |   n   \n\
      )     |  x  dx\n\
     /      |       \n\
    /      /        \n\
   /       -oo      \n\
  /       k         \n\
 /_____,            \n\
  k = 0             \
"""
    ucode_str = \
u("""\
∞                 \n\
⌠                 \n\
⎮   x             \n\
⎮  x  dx          \n\
⌡                 \n\
-∞                \n\
 ______           \n\
 ╲                \n\
  ╲       ∞       \n\
   ╲      ⌠       \n\
    ╲     ⎮   n   \n\
     ╲    ⎮  x  dx\n\
     ╱    ⌡       \n\
    ╱     -∞      \n\
   ╱     k        \n\
  ╱               \n\
 ╱                \n\
 ‾‾‾‾‾‾           \n\
 k = 0            \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum(k**(Integral(x**n, (x, -oo, oo))), (
        k, x + n + x**2 + n**2 + (x/n) + (1/x), Integral(x**x, (x, -oo, oo))))
    ascii_str = \
"""\
          oo                          \n\
           /                          \n\
          |                           \n\
          |   x                       \n\
          |  x  dx                    \n\
          |                           \n\
         /                            \n\
         -oo                          \n\
          ______                      \n\
          \     `                     \n\
           \                  oo      \n\
            \                  /      \n\
             \                |       \n\
              \               |   n   \n\
               )              |  x  dx\n\
              /               |       \n\
             /               /        \n\
            /                -oo      \n\
           /                k         \n\
          /_____,                     \n\
     2        2       1   x           \n\
k = n  + n + x  + x + - + -           \n\
                      x   n           \
"""
    ucode_str = \
u("""\
          ∞                          \n\
          ⌠                          \n\
          ⎮   x                      \n\
          ⎮  x  dx                   \n\
          ⌡                          \n\
          -∞                         \n\
           ______                    \n\
           ╲                         \n\
            ╲                ∞       \n\
             ╲               ⌠       \n\
              ╲              ⎮   n   \n\
               ╲             ⎮  x  dx\n\
               ╱             ⌡       \n\
              ╱              -∞      \n\
             ╱              k        \n\
            ╱                        \n\
           ╱                         \n\
           ‾‾‾‾‾‾                    \n\
     2        2       1   x          \n\
k = n  + n + x  + x + ─ + ─          \n\
                      x   n          \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum(k**(
        Integral(x**n, (x, -oo, oo))), (k, 0, x + n + x**2 + n**2 + (x/n) + (1/x)))
    ascii_str = \
"""\
 2        2       1   x           \n\
n  + n + x  + x + - + -           \n\
                  x   n           \n\
        ______                    \n\
        \     `                   \n\
         \                oo      \n\
          \                /      \n\
           \              |       \n\
            \             |   n   \n\
             )            |  x  dx\n\
            /             |       \n\
           /             /        \n\
          /              -oo      \n\
         /              k         \n\
        /_____,                   \n\
         k = 0                    \
"""
    ucode_str = \
u("""\
 2        2       1   x          \n\
n  + n + x  + x + ─ + ─          \n\
                  x   n          \n\
         ______                  \n\
         ╲                       \n\
          ╲              ∞       \n\
           ╲             ⌠       \n\
            ╲            ⎮   n   \n\
             ╲           ⎮  x  dx\n\
             ╱           ⌡       \n\
            ╱            -∞      \n\
           ╱            k        \n\
          ╱                      \n\
         ╱                       \n\
         ‾‾‾‾‾‾                  \n\
         k = 0                   \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum(x, (x, 0, oo))
    ascii_str = \
"""\
  oo   \n\
 __    \n\
 \\ `   \n\
  )   x\n\
 /_,   \n\
x = 0  \
"""
    ucode_str = \
u("""\
  ∞    \n\
 ___   \n\
 ╲     \n\
  ╲   x\n\
  ╱    \n\
 ╱     \n\
 ‾‾‾   \n\
x = 0  \
""")

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum(x**2, (x, 0, oo))
    ascii_str = \
u("""\
  oo    \n\
 ___    \n\
 \\  `   \n\
  \\    2\n\
  /   x \n\
 /__,   \n\
x = 0   \
""")
    ucode_str = \
u("""\
  ∞     \n\
 ___    \n\
 ╲      \n\
  ╲    2\n\
  ╱   x \n\
 ╱      \n\
 ‾‾‾    \n\
x = 0   \
""")

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum(x/2, (x, 0, oo))
    ascii_str = \
"""\
  oo   \n\
 ___   \n\
 \\  `  \n\
  \\   x\n\
   )  -\n\
  /   2\n\
 /__,  \n\
x = 0  \
"""
    ucode_str = \
u("""\
  ∞    \n\
 ____  \n\
 ╲     \n\
  ╲   x\n\
   ╲  ─\n\
   ╱  2\n\
  ╱    \n\
 ╱     \n\
 ‾‾‾‾  \n\
x = 0  \
""")

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum(x**3/2, (x, 0, oo))
    ascii_str = \
"""\
  oo    \n\
____    \n\
\\   `   \n\
 \\     3\n\
  \\   x \n\
  /   --\n\
 /    2 \n\
/___,   \n\
x = 0   \
"""
    ucode_str = \
u("""\
  ∞     \n\
 ____   \n\
 ╲      \n\
  ╲    3\n\
   ╲  x \n\
   ╱  ──\n\
  ╱   2 \n\
 ╱      \n\
 ‾‾‾‾   \n\
x = 0   \
""")

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum((x**3*y**(x/2))**n, (x, 0, oo))
    ascii_str = \
"""\
  oo          \n\
____          \n\
\\   `         \n\
 \\           n\n\
  \\   /    x\\ \n\
   )  |    -| \n\
  /   | 3  2| \n\
 /    \\x *y / \n\
/___,         \n\
x = 0         \
"""
    ucode_str = \
u("""\
  ∞           \n\
_____         \n\
╲             \n\
 ╲           n\n\
  ╲   ⎛    x⎞ \n\
   ╲  ⎜    ─⎟ \n\
   ╱  ⎜ 3  2⎟ \n\
  ╱   ⎝x ⋅y ⎠ \n\
 ╱            \n\
╱             \n\
‾‾‾‾‾         \n\
x = 0         \
""")

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum(1/x**2, (x, 0, oo))
    ascii_str = \
"""\
  oo    \n\
____    \n\
\\   `   \n\
 \\    1 \n\
  \\   --\n\
  /    2\n\
 /    x \n\
/___,   \n\
x = 0   \
"""
    ucode_str = \
u("""\
  ∞     \n\
 ____   \n\
 ╲      \n\
  ╲   1 \n\
   ╲  ──\n\
   ╱   2\n\
  ╱   x \n\
 ╱      \n\
 ‾‾‾‾   \n\
x = 0   \
""")

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum(1/y**(a/b), (x, 0, oo))
    ascii_str = \
"""\
  oo      \n\
____      \n\
\\   `     \n\
 \\     -a \n\
  \\    ---\n\
  /     b \n\
 /    y   \n\
/___,     \n\
x = 0     \
"""
    ucode_str = \
u("""\
  ∞       \n\
 ____     \n\
 ╲        \n\
  ╲    -a \n\
   ╲   ───\n\
   ╱    b \n\
  ╱   y   \n\
 ╱        \n\
 ‾‾‾‾     \n\
x = 0     \
""")

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum(1/y**(a/b), (x, 0, oo), (y, 1, 2))
    ascii_str = \
"""\
  2     oo     \n\
____  ____     \n\
\\   ` \\   `    \n\
 \\     \\     -a\n\
  \\     \\    --\n\
  /     /    b \n\
 /     /    y  \n\
/___, /___,    \n\
y = 1 x = 0    \
"""
    ucode_str = \
u("""\
  2     ∞      \n\
____  ____     \n\
╲     ╲        \n\
 ╲     ╲     -a\n\
  ╲     ╲    ──\n\
  ╱     ╱    b \n\
 ╱     ╱    y  \n\
╱     ╱        \n\
‾‾‾‾  ‾‾‾‾     \n\
y = 1 x = 0    \
""")
    expr = Sum(1/(1 + 1/(
        1 + 1/k)) + 1, (k, 111, 1 + 1/n), (k, 1/(1 + m), oo)) + 1/(1 + 1/k)
    ascii_str = \
"""\
               1                         \n\
           1 + -                         \n\
    oo         n                         \n\
  _____    _____                         \n\
  \\    `   \\    `                        \n\
   \\        \\     /        1    \\        \n\
    \\        \\    |1 + ---------|        \n\
     \\        \\   |          1  |     1  \n\
      )        )  |    1 + -----| + -----\n\
     /        /   |            1|       1\n\
    /        /    |        1 + -|   1 + -\n\
   /        /     \\            k/       k\n\
  /____,   /____,                        \n\
      1   k = 111                        \n\
k = -----                                \n\
    m + 1                                \
"""
    ucode_str = \
u("""\
               1                         \n\
           1 + ─                         \n\
    ∞          n                         \n\
  ______   ______                        \n\
  ╲        ╲                             \n\
   ╲        ╲     ⎛        1    ⎞        \n\
    ╲        ╲    ⎜1 + ─────────⎟        \n\
     ╲        ╲   ⎜          1  ⎟        \n\
      ╲        ╲  ⎜    1 + ─────⎟     1  \n\
      ╱        ╱  ⎜            1⎟ + ─────\n\
     ╱        ╱   ⎜        1 + ─⎟       1\n\
    ╱        ╱    ⎝            k⎠   1 + ─\n\
   ╱        ╱                           k\n\
  ╱        ╱                             \n\
  ‾‾‾‾‾‾   ‾‾‾‾‾‾                        \n\
      1   k = 111                        \n\
k = ─────                                \n\
    m + 1                                \
""")

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_units():
    expr = joule
    ascii_str = \
"""\
    2\n\
kg*m \n\
-----\n\
   2 \n\
  s  \
"""

    unicode_str = \
u("""\
    2\n\
kg⋅m \n\
─────\n\
   2 \n\
  s  \
""")
    assert upretty(expr) == unicode_str
    assert pretty(expr) == ascii_str


def test_pretty_Subs():
    f = Function('f')
    expr = Subs(f(x), x, ph**2)
    ascii_str = \
"""\
(f(x))|     2\n\
      |x=phi \
"""
    unicode_str = \
u("""\
(f(x))│   2\n\
      │x=φ \
""")

    assert pretty(expr) == ascii_str
    assert upretty(expr) == unicode_str

    expr = Subs(f(x).diff(x), x, 0)
    ascii_str = \
"""\
/d       \\|   \n\
|--(f(x))||   \n\
\\dx      /|x=0\
"""
    unicode_str = \
u("""\
⎛d       ⎞│   \n\
⎜──(f(x))⎟│   \n\
⎝dx      ⎠│x=0\
""")

    assert pretty(expr) == ascii_str
    assert upretty(expr) == unicode_str

    expr = Subs(f(x).diff(x)/y, (x, y), (0, Rational(1, 2)))
    ascii_str = \
"""\
/d       \\|          \n\
|--(f(x))||          \n\
|dx      ||          \n\
|--------||          \n\
\\   y    /|x=0, y=1/2\
"""
    unicode_str = \
u("""\
⎛d       ⎞│          \n\
⎜──(f(x))⎟│          \n\
⎜dx      ⎟│          \n\
⎜────────⎟│          \n\
⎝   y    ⎠│x=0, y=1/2\
""")

    assert pretty(expr) == ascii_str
    assert upretty(expr) == unicode_str


def test_gammas():
    assert upretty(lowergamma(x, y)) == u("γ(x, y)")
    assert upretty(uppergamma(x, y)) == u("Γ(x, y)")
    assert xpretty(gamma(x), use_unicode=True) == u('Γ(x)')


def test_hyper():
    expr = hyper((), (), z)
    ucode_str = \
u("""\
 ┌─  ⎛  │  ⎞\n\
 ├─  ⎜  │ z⎟\n\
0╵ 0 ⎝  │  ⎠\
""")
    ascii_str = \
"""\
  _         \n\
 |_  /  |  \\\n\
 |   |  | z|\n\
0  0 \\  |  /\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = hyper((), (1,), x)
    ucode_str = \
u("""\
 ┌─  ⎛  │  ⎞\n\
 ├─  ⎜  │ x⎟\n\
0╵ 1 ⎝1 │  ⎠\
""")
    ascii_str = \
"""\
  _         \n\
 |_  /  |  \\\n\
 |   |  | x|\n\
0  1 \\1 |  /\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = hyper([2], [1], x)
    ucode_str = \
u("""\
 ┌─  ⎛2 │  ⎞\n\
 ├─  ⎜  │ x⎟\n\
1╵ 1 ⎝1 │  ⎠\
""")
    ascii_str = \
"""\
  _         \n\
 |_  /2 |  \\\n\
 |   |  | x|\n\
1  1 \\1 |  /\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = hyper((pi/3, -2*k), (3, 4, 5, -3), x)
    ucode_str = \
u("""\
     ⎛  π         │  ⎞\n\
 ┌─  ⎜  ─, -2⋅k   │  ⎟\n\
 ├─  ⎜  3         │ x⎟\n\
2╵ 4 ⎜            │  ⎟\n\
     ⎝3, 4, 5, -3 │  ⎠\
""")
    ascii_str = \
"""\
                      \n\
  _  /  pi        |  \\\n\
 |_  |  --, -2*k  |  |\n\
 |   |  3         | x|\n\
2  4 |            |  |\n\
     \\3, 4, 5, -3 |  /\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = hyper((pi, S('2/3'), -2*k), (3, 4, 5, -3), x**2)
    ucode_str = \
u("""\
 ┌─  ⎛π, 2/3, -2⋅k │  2⎞\n\
 ├─  ⎜             │ x ⎟\n\
3╵ 4 ⎝3, 4, 5, -3  │   ⎠\
""")
    ascii_str = \
"""\
  _                      \n\
 |_  /pi, 2/3, -2*k |  2\\\n\
 |   |              | x |\n\
3  4 \\ 3, 4, 5, -3  |   /\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = hyper([1, 2], [3, 4], 1/(1/(1/(1/x + 1) + 1) + 1))
    ucode_str = \
u("""\
     ⎛     │       1      ⎞\n\
     ⎜     │ ─────────────⎟\n\
     ⎜     │         1    ⎟\n\
 ┌─  ⎜1, 2 │ 1 + ─────────⎟\n\
 ├─  ⎜     │           1  ⎟\n\
2╵ 2 ⎜3, 4 │     1 + ─────⎟\n\
     ⎜     │             1⎟\n\
     ⎜     │         1 + ─⎟\n\
     ⎝     │             x⎠\
""")

    ascii_str = \
"""\
                           \n\
     /     |       1      \\\n\
     |     | -------------|\n\
  _  |     |         1    |\n\
 |_  |1, 2 | 1 + ---------|\n\
 |   |     |           1  |\n\
2  2 |3, 4 |     1 + -----|\n\
     |     |             1|\n\
     |     |         1 + -|\n\
     \\     |             x/\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_meijerg():
    expr = meijerg([pi, pi, x], [1], [0, 1], [1, 2, 3], z)
    ucode_str = \
u("""\
╭─╮2, 3 ⎛π, π, x     1    │  ⎞\n\
│╶┐     ⎜                 │ z⎟\n\
╰─╯4, 5 ⎝ 0, 1    1, 2, 3 │  ⎠\
""")
    ascii_str = \
"""\
 __2, 3 /pi, pi, x     1    |  \\\n\
/__     |                   | z|\n\
\_|4, 5 \\  0, 1     1, 2, 3 |  /\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = meijerg([1, pi/7], [2, pi, 5], [], [], z**2)
    ucode_str = \
u("""\
        ⎛   π          │   ⎞\n\
╭─╮0, 2 ⎜1, ─  2, π, 5 │  2⎟\n\
│╶┐     ⎜   7          │ z ⎟\n\
╰─╯5, 0 ⎜              │   ⎟\n\
        ⎝              │   ⎠\
""")
    ascii_str = \
"""\
        /   pi           |   \\\n\
 __0, 2 |1, --  2, pi, 5 |  2|\n\
/__     |   7            | z |\n\
\_|5, 0 |                |   |\n\
        \\                |   /\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    ucode_str = \
u("""\
╭─╮ 1, 10 ⎛1, 1, 1, 1, 1, 1, 1, 1, 1, 1  1 │  ⎞\n\
│╶┐       ⎜                                │ z⎟\n\
╰─╯11,  2 ⎝             1                1 │  ⎠\
""")
    ascii_str = \
"""\
 __ 1, 10 /1, 1, 1, 1, 1, 1, 1, 1, 1, 1  1 |  \\\n\
/__       |                                | z|\n\
\_|11,  2 \\             1                1 |  /\
"""

    expr = meijerg([1]*10, [1], [1], [1], z)
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = meijerg([1, 2, ], [4, 3], [3], [4, 5], 1/(1/(1/(1/x + 1) + 1) + 1))

    ucode_str = \
u("""\
        ⎛           │       1      ⎞\n\
        ⎜           │ ─────────────⎟\n\
        ⎜           │         1    ⎟\n\
╭─╮1, 2 ⎜1, 2  4, 3 │ 1 + ─────────⎟\n\
│╶┐     ⎜           │           1  ⎟\n\
╰─╯4, 3 ⎜ 3    4, 5 │     1 + ─────⎟\n\
        ⎜           │             1⎟\n\
        ⎜           │         1 + ─⎟\n\
        ⎝           │             x⎠\
""")

    ascii_str = \
"""\
        /           |       1      \\\n\
        |           | -------------|\n\
        |           |         1    |\n\
 __1, 2 |1, 2  4, 3 | 1 + ---------|\n\
/__     |           |           1  |\n\
\_|4, 3 | 3    4, 5 |     1 + -----|\n\
        |           |             1|\n\
        |           |         1 + -|\n\
        \\           |             x/\
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Integral(expr, x)

    ucode_str = \
u("""\
⌠                                        \n\
⎮         ⎛           │       1      ⎞   \n\
⎮         ⎜           │ ─────────────⎟   \n\
⎮         ⎜           │         1    ⎟   \n\
⎮ ╭─╮1, 2 ⎜1, 2  4, 3 │ 1 + ─────────⎟   \n\
⎮ │╶┐     ⎜           │           1  ⎟ dx\n\
⎮ ╰─╯4, 3 ⎜ 3    4, 5 │     1 + ─────⎟   \n\
⎮         ⎜           │             1⎟   \n\
⎮         ⎜           │         1 + ─⎟   \n\
⎮         ⎝           │             x⎠   \n\
⌡                                        \
""")

    ascii_str = \
"""\
  /                                       \n\
 |                                        \n\
 |         /           |       1      \\   \n\
 |         |           | -------------|   \n\
 |         |           |         1    |   \n\
 |  __1, 2 |1, 2  4, 3 | 1 + ---------|   \n\
 | /__     |           |           1  | dx\n\
 | \_|4, 3 | 3    4, 5 |     1 + -----|   \n\
 |         |           |             1|   \n\
 |         |           |         1 + -|   \n\
 |         \\           |             x/   \n\
 |                                        \n\
/                                         \
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_noncommutative():
    A, B, C = symbols('A,B,C', commutative=False)

    expr = A*B*C**-1
    ascii_str = \
"""\
     -1\n\
A*B*C  \
"""
    ucode_str = \
u("""\
     -1\n\
A⋅B⋅C  \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = C**-1*A*B
    ascii_str = \
"""\
 -1    \n\
C  *A*B\
"""
    ucode_str = \
u("""\
 -1    \n\
C  ⋅A⋅B\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = A*C**-1*B
    ascii_str = \
"""\
   -1  \n\
A*C  *B\
"""
    ucode_str = \
u("""\
   -1  \n\
A⋅C  ⋅B\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = A*C**-1*B/x
    ascii_str = \
"""\
   -1  \n\
A*C  *B\n\
-------\n\
   x   \
"""
    ucode_str = \
u("""\
   -1  \n\
A⋅C  ⋅B\n\
───────\n\
   x   \
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_special_functions():
    x, y = symbols("x y")

    # atan2
    expr = atan2(y/sqrt(200), sqrt(x))
    ascii_str = \
"""\
     /  ___         \\\n\
     |\\/ 2 *y    ___|\n\
atan2|-------, \\/ x |\n\
     \\   20         /\
"""
    ucode_str = \
u("""\
     ⎛  ___         ⎞\n\
     ⎜╲╱ 2 ⋅y    ___⎟\n\
atan2⎜───────, ╲╱ x ⎟\n\
     ⎝   20         ⎠\
""")
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_geometry():
    e = Segment((0, 1), (0, 2))
    assert pretty(e) == 'Segment(Point(0, 1), Point(0, 2))'
    e = Ray((1, 1), angle=4.02*pi)
    assert pretty(e) == 'Ray(Point(1, 1), Point(2, tan(pi/50) + 1))'


def test_expint():
    expr = Ei(x)
    string = 'Ei(x)'
    assert pretty(expr) == string
    assert upretty(expr) == string

    expr = expint(1, z)
    ucode_str = u("E₁(z)")
    ascii_str = "expint(1, z)"
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    assert pretty(Shi(x)) == 'Shi(x)'
    assert pretty(Si(x)) == 'Si(x)'
    assert pretty(Ci(x)) == 'Ci(x)'
    assert pretty(Chi(x)) == 'Chi(x)'
    assert upretty(Shi(x)) == 'Shi(x)'
    assert upretty(Si(x)) == 'Si(x)'
    assert upretty(Ci(x)) == 'Ci(x)'
    assert upretty(Chi(x)) == 'Chi(x)'


def test_elliptic_functions():
    ascii_str = \
"""\
 /  1  \\\n\
K|-----|\n\
 \z + 1/\
"""
    ucode_str = \
u("""\
 ⎛  1  ⎞\n\
K⎜─────⎟\n\
 ⎝z + 1⎠\
""")
    expr = elliptic_k(1/(z + 1))
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    ascii_str = \
"""\
 / |  1  \\\n\
F|1|-----|\n\
 \ |z + 1/\
"""
    ucode_str = \
u("""\
 ⎛ │  1  ⎞\n\
F⎜1│─────⎟\n\
 ⎝ │z + 1⎠\
""")
    expr = elliptic_f(1, 1/(1 + z))
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    ascii_str = \
"""\
 /  1  \\\n\
E|-----|\n\
 \z + 1/\
"""
    ucode_str = \
u("""\
 ⎛  1  ⎞\n\
E⎜─────⎟\n\
 ⎝z + 1⎠\
""")
    expr = elliptic_e(1/(z + 1))
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    ascii_str = \
"""\
 / |  1  \\\n\
E|1|-----|\n\
 \ |z + 1/\
"""
    ucode_str = \
u("""\
 ⎛ │  1  ⎞\n\
E⎜1│─────⎟\n\
 ⎝ │z + 1⎠\
""")
    expr = elliptic_e(1, 1/(1 + z))
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    ascii_str = \
"""\
  / |4\\\n\
Pi|3|-|\n\
  \ |x/\
"""
    ucode_str = \
u("""\
 ⎛ │4⎞\n\
Π⎜3│─⎟\n\
 ⎝ │x⎠\
""")
    expr = elliptic_pi(3, 4/x)
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    ascii_str = \
"""\
  /   4| \\\n\
Pi|3; -|6|\n\
  \   x| /\
"""
    ucode_str = \
u("""\
 ⎛   4│ ⎞\n\
Π⎜3; ─│6⎟\n\
 ⎝   x│ ⎠\
""")
    expr = elliptic_pi(3, 4/x, 6)
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_RandomDomain():
    from sympy.stats import Normal, Die, Exponential, pspace, where
    X = Normal('x1', 0, 1)
    assert upretty(where(X > 0)) == u("Domain: 0 < x₁ ∧ x₁ < ∞")

    D = Die('d1', 6)
    assert upretty(where(D > 4)) == u('Domain: d₁ = 5 ∨ d₁ = 6')

    A = Exponential('a', 1)
    B = Exponential('b', 1)
    assert upretty(pspace(Tuple(A, B)).domain) == \
        u('Domain: 0 ≤ a ∧ 0 ≤ b ∧ a < ∞ ∧ b < ∞')


def test_PrettyPoly():
    F = QQ.frac_field(x, y)
    R = QQ.poly_ring(x, y)

    expr = F.convert(x/(x + y))
    assert pretty(expr) == "x/(x + y)"
    assert upretty(expr) == u("x/(x + y)")

    expr = R.convert(x + y)
    assert pretty(expr) == "x + y"
    assert upretty(expr) == u("x + y")


def test_issue_6285():
    assert pretty(Pow(2, -5, evaluate=False)) == '1 \n--\n 5\n2 '
    assert pretty(Pow(x, (1/pi))) == 'pi___\n\\/ x '


def test_issue_6359():
    assert pretty(Integral(x**2, x)**2) == \
"""\
          2
/  /     \ \n\
| |      | \n\
| |  2   | \n\
| | x  dx| \n\
| |      | \n\
\/       / \
"""
    assert upretty(Integral(x**2, x)**2) == \
u("""\
         2
⎛⌠      ⎞ \n\
⎜⎮  2   ⎟ \n\
⎜⎮ x  dx⎟ \n\
⎝⌡      ⎠ \
""")

    assert pretty(Sum(x**2, (x, 0, 1))**2) == \
"""\
          2
/  1     \\ \n\
| ___    | \n\
| \\  `   | \n\
|  \\    2| \n\
|  /   x | \n\
| /__,   | \n\
\\x = 0   / \
"""
    assert upretty(Sum(x**2, (x, 0, 1))**2) == \
u("""\
          2
⎛  1     ⎞ \n\
⎜ ___    ⎟ \n\
⎜ ╲      ⎟ \n\
⎜  ╲    2⎟ \n\
⎜  ╱   x ⎟ \n\
⎜ ╱      ⎟ \n\
⎜ ‾‾‾    ⎟ \n\
⎝x = 0   ⎠ \
""")

    assert pretty(Product(x**2, (x, 1, 2))**2) == \
"""\
           2
/  2      \\ \n\
|______   | \n\
||    |  2| \n\
||    | x | \n\
||    |   | \n\
\\x = 1    / \
"""
    assert upretty(Product(x**2, (x, 1, 2))**2) == \
u("""\
           2
⎛  2      ⎞ \n\
⎜┬────┬   ⎟ \n\
⎜│    │  2⎟ \n\
⎜│    │ x ⎟ \n\
⎜│    │   ⎟ \n\
⎝x = 1    ⎠ \
""")

    f = Function('f')
    assert pretty(Derivative(f(x), x)**2) == \
"""\
          2
/d       \\ \n\
|--(f(x))| \n\
\\dx      / \
"""
    assert upretty(Derivative(f(x), x)**2) == \
u("""\
          2
⎛d       ⎞ \n\
⎜──(f(x))⎟ \n\
⎝dx      ⎠ \
""")

def test_issue_6739():
    ascii_str = \
"""\
  1  \n\
-----\n\
  ___\n\
\/ x \
"""
    ucode_str = \
u("""\
  1  \n\
─────\n\
  ___\n\
╲╱ x \
""")
    assert pretty(1/sqrt(x)) == ascii_str
    assert upretty(1/sqrt(x)) == ucode_str


def test_complicated_symbol_unchanged():
    for symb_name in ["dexpr2_d1tau", "dexpr2^d1tau"]:
        assert pretty(Symbol(symb_name)) == symb_name


def test_categories():
    from sympy.categories import (Object, IdentityMorphism,
        NamedMorphism, Category, Diagram, DiagramGrid)

    A1 = Object("A1")
    A2 = Object("A2")
    A3 = Object("A3")

    f1 = NamedMorphism(A1, A2, "f1")
    f2 = NamedMorphism(A2, A3, "f2")
    id_A1 = IdentityMorphism(A1)

    K1 = Category("K1")

    assert pretty(A1) == "A1"
    assert upretty(A1) == u("A₁")

    assert pretty(f1) == "f1:A1-->A2"
    assert upretty(f1) == u("f₁:A₁——▶A₂")
    assert pretty(id_A1) == "id:A1-->A1"
    assert upretty(id_A1) == u("id:A₁——▶A₁")

    assert pretty(f2*f1) == "f2*f1:A1-->A3"
    assert upretty(f2*f1) == u("f₂∘f₁:A₁——▶A₃")

    assert pretty(K1) == "K1"
    assert upretty(K1) == u("K₁")

    # Test how diagrams are printed.
    d = Diagram()
    assert pretty(d) == "EmptySet()"
    assert upretty(d) == u("∅")

    d = Diagram({f1: "unique", f2: S.EmptySet})
    assert pretty(d) == "{f2*f1:A1-->A3: EmptySet(), id:A1-->A1: " \
        "EmptySet(), id:A2-->A2: EmptySet(), id:A3-->A3: " \
        "EmptySet(), f1:A1-->A2: {unique}, f2:A2-->A3: EmptySet()}"

    assert upretty(d) == u("{f₂∘f₁:A₁——▶A₃: ∅, id:A₁——▶A₁: ∅, " \
        "id:A₂——▶A₂: ∅, id:A₃——▶A₃: ∅, f₁:A₁——▶A₂: {unique}, f₂:A₂——▶A₃: ∅}")

    d = Diagram({f1: "unique", f2: S.EmptySet}, {f2 * f1: "unique"})
    assert pretty(d) == "{f2*f1:A1-->A3: EmptySet(), id:A1-->A1: " \
        "EmptySet(), id:A2-->A2: EmptySet(), id:A3-->A3: " \
        "EmptySet(), f1:A1-->A2: {unique}, f2:A2-->A3: EmptySet()}" \
        " ==> {f2*f1:A1-->A3: {unique}}"
    assert upretty(d) == u("{f₂∘f₁:A₁——▶A₃: ∅, id:A₁——▶A₁: ∅, id:A₂——▶A₂: " \
        "∅, id:A₃——▶A₃: ∅, f₁:A₁——▶A₂: {unique}, f₂:A₂——▶A₃: ∅}" \
        " ══▶ {f₂∘f₁:A₁——▶A₃: {unique}}")

    grid = DiagramGrid(d)
    assert pretty(grid) == "A1  A2\n      \nA3    "
    assert upretty(grid) == u("A₁  A₂\n      \nA₃    ")


def test_PrettyModules():
    R = QQ.old_poly_ring(x, y)
    F = R.free_module(2)
    M = F.submodule([x, y], [1, x**2])

    ucode_str = \
u("""\
       2\n\
ℚ[x, y] \
""")
    ascii_str = \
"""\
        2\n\
QQ[x, y] \
"""

    assert upretty(F) == ucode_str
    assert pretty(F) == ascii_str

    ucode_str = \
u("""\
╱        ⎡    2⎤╲\n\
╲[x, y], ⎣1, x ⎦╱\
""")
    ascii_str = \
"""\
              2  \n\
<[x, y], [1, x ]>\
"""

    assert upretty(M) == ucode_str
    assert pretty(M) == ascii_str

    I = R.ideal(x**2, y)

    ucode_str = \
u("""\
╱ 2   ╲\n\
╲x , y╱\
""")

    ascii_str = \
"""\
  2    \n\
<x , y>\
"""

    assert upretty(I) == ucode_str
    assert pretty(I) == ascii_str

    Q = F / M

    ucode_str = \
u("""\
            2    \n\
     ℚ[x, y]     \n\
─────────────────\n\
╱        ⎡    2⎤╲\n\
╲[x, y], ⎣1, x ⎦╱\
""")

    ascii_str = \
"""\
            2    \n\
    QQ[x, y]     \n\
-----------------\n\
              2  \n\
<[x, y], [1, x ]>\
"""

    assert upretty(Q) == ucode_str
    assert pretty(Q) == ascii_str

    ucode_str = \
u("""\
╱⎡    3⎤                                                ╲\n\
│⎢   x ⎥   ╱        ⎡    2⎤╲           ╱        ⎡    2⎤╲│\n\
│⎢1, ──⎥ + ╲[x, y], ⎣1, x ⎦╱, [2, y] + ╲[x, y], ⎣1, x ⎦╱│\n\
╲⎣   2 ⎦                                                ╱\
""")

    ascii_str = \
"""\
      3                                                  \n\
     x                   2                           2   \n\
<[1, --] + <[x, y], [1, x ]>, [2, y] + <[x, y], [1, x ]>>\n\
     2                                                   \
"""


def test_QuotientRing():
    R = QQ.old_poly_ring(x)/[x**2 + 1]

    ucode_str = \
u("""\
  ℚ[x]  \n\
────────\n\
╱ 2    ╲\n\
╲x  + 1╱\
""")

    ascii_str = \
"""\
 QQ[x]  \n\
--------\n\
  2     \n\
<x  + 1>\
"""

    assert upretty(R) == ucode_str
    assert pretty(R) == ascii_str

    ucode_str = \
u("""\
    ╱ 2    ╲\n\
1 + ╲x  + 1╱\
""")

    ascii_str = \
"""\
      2     \n\
1 + <x  + 1>\
"""

    assert upretty(R.one) == ucode_str
    assert pretty(R.one) == ascii_str


def test_Homomorphism():
    from sympy.polys.agca import homomorphism

    R = QQ.old_poly_ring(x)

    expr = homomorphism(R.free_module(1), R.free_module(1), [0])

    ucode_str = \
u("""\
          1         1\n\
[0] : ℚ[x]  ──> ℚ[x] \
""")

    ascii_str = \
"""\
           1          1\n\
[0] : QQ[x]  --> QQ[x] \
"""

    assert upretty(expr) == ucode_str
    assert pretty(expr) == ascii_str

    expr = homomorphism(R.free_module(2), R.free_module(2), [0, 0])

    ucode_str = \
u("""\
⎡0  0⎤       2         2\n\
⎢    ⎥ : ℚ[x]  ──> ℚ[x] \n\
⎣0  0⎦                  \
""")

    ascii_str = \
"""\
[0  0]        2          2\n\
[    ] : QQ[x]  --> QQ[x] \n\
[0  0]                    \
"""

    assert upretty(expr) == ucode_str
    assert pretty(expr) == ascii_str

    expr = homomorphism(R.free_module(1), R.free_module(1) / [[x]], [0])

    ucode_str = \
u("""\
                    1\n\
          1     ℚ[x] \n\
[0] : ℚ[x]  ──> ─────\n\
                <[x]>\
""")

    ascii_str = \
"""\
                      1\n\
           1     QQ[x] \n\
[0] : QQ[x]  --> ------\n\
                 <[x]> \
"""

    assert upretty(expr) == ucode_str
    assert pretty(expr) == ascii_str


def test_Tr():
    A, B = symbols('A B', commutative=False)
    t = Tr(A*B)
    assert pretty(t) == r'Tr(A*B)'
    assert upretty(t) == u('Tr(A⋅B)')


def test_pretty_Add():
    eq = Mul(-2, x - 2, evaluate=False) + 5
    assert pretty(eq) == '-2*(x - 2) + 5'


def test_issue_7179():
    assert upretty(Not(Equivalent(x, y))) == u('x ≢ y')
    assert upretty(Not(Implies(x, y))) == u('x ↛ y')


def test_issue_7180():
    assert upretty(Equivalent(x, y)) == u('x ≡ y')


def test_pretty_Complement():
    assert pretty(S.Reals - S.Naturals) == '(-oo, oo) \ Naturals()'
    assert upretty(S.Reals - S.Naturals) == u('ℝ \ ℕ')


def test_pretty_SymmetricDifference():
    from sympy import SymmetricDifference, Interval
    from sympy.utilities.pytest import raises
    assert upretty(SymmetricDifference(Interval(2,3), Interval(3,5), \
           evaluate = False)) == u('[2, 3] ∆ [3, 5]')
    with raises(NotImplementedError):
        pretty(SymmetricDifference(Interval(2,3), Interval(3,5), evaluate = False))


def test_pretty_Contains():
    assert pretty(Contains(x, S.Integers)) == 'Contains(x, Integers())'
    assert upretty(Contains(x, S.Integers)) == u('x ∈ ℤ')


def test_issue_8292():
    from sympy.core import sympify
    e = sympify('((x+x**4)/(x-1))-(2*(x-1)**4/(x-1)**4)', evaluate=False)
    ucode_str = \
u("""\
           4    4    \n\
  2⋅(x - 1)    x  + x\n\
- ────────── + ──────\n\
          4    x - 1 \n\
   (x - 1)           \
""")
    ascii_str = \
"""\
           4    4    \n\
  2*(x - 1)    x  + x\n\
- ---------- + ------\n\
          4    x - 1 \n\
   (x - 1)           \
"""
    assert pretty(e) == ascii_str
    assert upretty(e) == ucode_str


def test_issue_4335():
    expr = -y(x).diff(x)
    ucode_str = \
u("""\
 d       \n\
-──(y(x))\n\
 dx      \
""")
    ascii_str = \
"""\
  d       \n\
- --(y(x))\n\
  dx      \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_issue_8344():
    from sympy.core import sympify
    e = sympify('2*x*y**2/1**2 + 1', evaluate=False)
    ucode_str = \
u("""\
     2    \n\
2⋅x⋅y     \n\
────── + 1\n\
   2      \n\
  1       \
""")
    assert upretty(e) == ucode_str


def test_issue_6324():
    x = Pow(2, 3, evaluate=False)
    y = Pow(10, -2, evaluate=False)
    e = Mul(x, y, evaluate=False)
    ucode_str = \
u("""\
  3\n\
 2 \n\
───\n\
  2\n\
10 \
""")
    assert upretty(e) == ucode_str


@XFAIL
def test_issue_9057():
    from sympy.functions import beta
    # beta function should take 2 arguments
    raises(TypeError, lambda: beta(2))

    # now beta is a pure Symbol
    from sympy.abc import alpha, phi, beta, t

    # as expected, no errors raised by these
    alpha(2)
    alpha(2.5)
    alpha(t)
    beta(2)
    beta(t)
    phi(2)
    phi(t)

    # but this raises TypeError (became the beta function )
    beta(2.5)


def test_issue_9057_phi():
    from sympy.abc import phi
    # similar to above, but segfaults on my system
    phi(2.5)


def test_issue_6134():
    # FIXME: this segfaults for me---as in test_issue_9057_phi---probably
    # because somewhere it tries numerical quadrature on the integral).
    from sympy.abc import lamda, phi, t

    e = lamda*x*Integral(phi(t)*pi*sin(pi*t), (t, 0, 1)) + lamda*x**2*Integral(phi(t)*2*pi*sin(2*pi*t), (t, 0, 1))
    ucode_str = \
u("""\
     1                              1                   \n\
   2 ⌠                              ⌠                   \n\
λ⋅x ⋅⎮ 2⋅π⋅φ(t)⋅sin(2⋅π⋅t) dt + λ⋅x⋅⎮ π⋅φ(t)⋅sin(π⋅t) dt\n\
     ⌡                              ⌡                   \n\
     0                              0                   \
""")
    assert upretty(e) == ucode_str
