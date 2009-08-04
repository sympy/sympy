# -*- coding: utf-8 -*-

from sympy import symbols, Symbol, sin, cos, Matrix, Integral, pi, sqrt,\
        Function, Rational, tan, oo, Limit, ceiling, floor, conjugate, exp, I
from sympy.printing.pretty import pretty

x,y,k = symbols('xyk')
th  = Symbol('theta')
ph  = Symbol('phi')

def upretty(expr):
    return pretty(expr, use_unicode=True)

def test_upretty_greek():
    assert upretty( oo ) == u'∞'
    assert upretty( Symbol('alpha^+_1') )   ==  u'α⁺₁'
    assert upretty( Symbol('beta') )    == u'β'

def test_upretty_multiindex():
    assert upretty( Symbol('beta12') )  == u'β₁₂'
    assert upretty( Symbol('Y00') )     == u'Y₀₀'
    assert upretty( Symbol('Y_00') )    == u'Y₀₀'
    assert upretty( Symbol('F^+-') )    == u'F⁺⁻'

def test_upretty_sub_super():
    assert upretty( Symbol('beta_1_2') )  == u'β₁ ₂'
    assert upretty( Symbol('beta^1^2') )  == u'β¹ ²'
    assert upretty( Symbol('beta_1^2') )  == u'β²₁'
    assert upretty( Symbol('beta_10_20') )  == u'β₁₀ ₂₀'
    assert upretty( Symbol('beta_ax_gamma^i') )  == u'βⁱₐₓ ᵧ'
    assert upretty( Symbol("F^1^2_3_4") )  == u'F¹ ²₃ ₄'
    assert upretty( Symbol("F_1_2^3^4") )  == u'F³ ⁴₁ ₂'
    assert upretty( Symbol("F_1_2_3_4") )  == u'F₁ ₂ ₃ ₄'
    assert upretty( Symbol("F^1^2^3^4") )  == u'F¹ ² ³ ⁴'

def test_upretty_subs_missingin_24():
    assert upretty( Symbol('F_beta') )  == u'Fᵦ'
    assert upretty( Symbol('F_gamma') ) == u'Fᵧ'
    assert upretty( Symbol('F_rho') )   == u'Fᵨ'
    assert upretty( Symbol('F_phi') )   == u'Fᵩ'
    assert upretty( Symbol('F_chi') )   == u'Fᵪ'

    assert upretty( Symbol('F_a') )     == u'Fₐ'
    assert upretty( Symbol('F_e') )     == u'Fₑ'
    assert upretty( Symbol('F_i') )     == u'Fᵢ'
    assert upretty( Symbol('F_o') )     == u'Fₒ'
    assert upretty( Symbol('F_r') )     == u'Fᵣ'
    assert upretty( Symbol('F_u') )     == u'Fᵤ'
    assert upretty( Symbol('F_v') )     == u'Fᵥ'
    assert upretty( Symbol('F_x') )     == u'Fₓ'


def test_upretty_nicemul():
    assert upretty(x*y) ==  u'x⋅y'


def test_upretty_nicerat():
    u = upretty(y*x**-2)
    s = \
u"""\
y \n\
──\n\
 2\n\
x \
"""
    assert u == s

    u = upretty(y**Rational(3,2) * x**Rational(-5,2))
    s = \
u"""\
 3/2
y   \n\
────\n\
 5/2\n\
x   \
"""
    assert u == s

    u = upretty(sin(x)**3/tan(x)**2)
    s = \
u"""\
   3   \n\
sin (x)\n\
───────\n\
   2   \n\
tan (x)\
"""

def test_upretty_funcbraces():
    f = Function('f')
    u = upretty(f(x/(y+1), y))
    s1 = \
u"""\
 ⎛  x     ⎞
f⎜─────, y⎟
 ⎝1 + y   ⎠\
"""
    s2 = \
u"""\
 ⎛  x     ⎞
f⎜─────, y⎟
 ⎝y + 1   ⎠\
"""
    assert u in [s1, s2]

def test_upretty_sqrt():
    u = upretty( sqrt((sqrt(x+1))+1) )
    s1 = \
u"""\
   ⎽⎽⎽⎽⎽⎽⎽⎽⎽⎽⎽⎽⎽⎽⎽
  ╱       ⎽⎽⎽⎽⎽⎽⎽ \n\
╲╱  1 + ╲╱ 1 + x  \
"""
    s2 = \
u"""\
   ⎽⎽⎽⎽⎽⎽⎽⎽⎽⎽⎽⎽⎽⎽⎽
  ╱   ⎽⎽⎽⎽⎽⎽⎽     \n\
╲╱  ╲╱ x + 1  + 1 \
"""
    assert u in [s1, s2]

def test_upretty_conjugate():
    f = Function('f')
    u = upretty( conjugate( f(1 + conjugate(f(x))) ) )
    s1 = \
u"""\
⎽⎽⎽⎽⎽⎽⎽⎽⎽⎽⎽
 ⎛    ⎽⎽⎽⎽⎞
f⎝1 + f(x)⎠\
"""
    s2 = \
u"""\
⎽⎽⎽⎽⎽⎽⎽⎽⎽⎽⎽
 ⎛⎽⎽⎽⎽    ⎞
f⎝f(x) + 1⎠\
"""
    assert u in [s1, s2]

def test_upretty_abs():
    assert upretty(abs(x)) == u'│x│'

    u = upretty( abs(1 / (y - abs(x))) )
    s = \
u"""\
│   1   │
│───────│
│y - │x││\
"""
    assert u == s

def test_upretty_floor():
    assert upretty(floor(x)) == u'⌊x⌋'

    u = upretty( floor(1 / (y - floor(x))) )
    s = \
u"""\
⎢   1   ⎥
⎢───────⎥
⎣y - ⌊x⌋⎦\
"""
    assert u == s

def test_upretty_ceiling():
    assert upretty(ceiling(x)) == u'⌈x⌉'

    u = upretty( ceiling(1 / (y - ceiling(x))) )
    s = \
u"""\
⎡   1   ⎤
⎢───────⎥
⎢y - ⌈x⌉⎥\
"""
    assert u == s


def test_upretty_diff():
    alpha = Symbol('alpha')
    beta  = Function('beta')

    u = upretty( beta(alpha).diff(alpha) )
    s = \
u"""\
d       \n\
──(β(α))\n\
dα      \
"""
    assert u == s

def test_upretty_integral():
    u = upretty( Integral(sin(th)/cos(ph), (th,0,pi), (ph, 0, 2*pi)) )
    s = \
u"""\
2⋅π π             \n\
 ⌠  ⌠             \n\
 ⎮  ⎮ sin(θ)      \n\
 ⎮  ⎮ ────── dθ dφ\n\
 ⎮  ⎮ cos(φ)      \n\
 ⌡  ⌡             \n\
 0  0             \
"""
    assert u == s

    u = upretty( Integral(x**2*sin(y), (x,0,1), (y,0,pi)) )
    s = \
u"""\
π 1                \n\
⌠ ⌠                \n\
⎮ ⎮  2             \n\
⎮ ⎮ x ⋅sin(y) dx dy\n\
⌡ ⌡                \n\
0 0                \
"""
    assert u == s

    u = upretty( Integral(sin(x), (x, None, 1)) )
    s = \
u"""1          \n\
⌠          \n\
⎮ sin(x) dx\n\
⌡          \n\
           \
"""
    assert u == s

def test_upretty_limit():
    u = upretty( Limit(x, x, oo) )
    s = \
u"""\
lim x
x->∞ \
"""
    assert u == s

    u = upretty( Limit(x**2, x, 0) )
    s = \
u"""\
     2
lim x \n\
x->0  \
"""
    assert u == s

    u = upretty( Limit(1/x, x, 0) )
    s = \
u"""\
    1
lim ─
x->0x\
"""
    assert u == s

    u = upretty( Limit(sin(x)/x, x, 0) )
    s = \
u"""\
    sin(x)
lim ──────
x->0  x   \
"""
    assert u == s


def test_upretty_matrix():
    u = upretty( Matrix([[x**2+1, 1], [y, x+y]]) )
    s1 = \
u"""\
⎡     2       ⎤
⎢1 + x     1  ⎥
⎢             ⎥
⎣  y     x + y⎦\
"""
    s2 = \
u"""\
⎡ 2           ⎤
⎢x  + 1    1  ⎥
⎢             ⎥
⎣  y     y + x⎦\
"""
    assert u in [s1, s2]


def test_upretty_matrix2():
    m = Matrix([[x/y, y, th], [0, exp(I*k*ph), 1]])
    u = upretty(m)
    s = \
u"""\
⎡x           ⎤
⎢─    y     θ⎥
⎢y           ⎥
⎢            ⎥
⎢    ⅈ⋅k⋅φ   ⎥
⎣0  ℯ       1⎦\
"""
    assert u == s



def test_upretty_seq():
    assert upretty([]) == '[]'
    assert upretty(()) == '()'
    assert upretty({}) == '{}'
    assert upretty(set()) == 'set()'
    assert upretty(frozenset()) == 'frozenset()'

    e = [x**2, 1/x, x, y, sin(th)**2/cos(ph)**2]
    u = upretty(e)
    s = \
u"""\
⎡                2   ⎤
⎢ 2  1        sin (θ)⎥
⎢x , ─, x, y, ───────⎥
⎢    x           2   ⎥
⎣             cos (φ)⎦\
"""
    assert u == s

    assert upretty((x,)) == '(x,)'

    p = upretty((1/x,))
    s = \
u"""\
⎛1 ⎞
⎜─,⎟
⎝x ⎠\
"""
    assert p == s


    e = tuple(e)
    u = upretty(e)
    s = \
u"""\
⎛                2   ⎞
⎜ 2  1        sin (θ)⎟
⎜x , ─, x, y, ───────⎟
⎜    x           2   ⎟
⎝             cos (φ)⎠\
"""
    assert u == s

    e_= dict(enumerate(e))
    u = upretty(e_)
    s = \
u"""\
⎧                               2   ⎫
⎪   1      2                 sin (θ)⎪
⎨1: ─, 0: x , 2: x, 3: y, 4: ───────⎬
⎪   x                           2   ⎪
⎩                            cos (φ)⎭\
"""
    assert u == s

    e_= set(e)
    u = upretty(e_)
    s = \
u"""\
   ⎛                2   ⎞
   ⎜      1   2  sin (θ)⎟
set⎜x, y, ─, x , ───────⎟
   ⎜      x         2   ⎟
   ⎝             cos (φ)⎠\
"""

    e = {x: sin(x)}
    u = upretty(e)
    s = \
u"""\
{x: sin(x)}\
"""
    assert u == s

    e = {1/x: 1/y, x: sin(x)**2}
    u = upretty(e)
    s = \
u"""\
⎧      2     1  1⎫
⎨x: sin (x), ─: ─⎬
⎩            x  y⎭\
"""
    assert u == s


def test_upretty_seq_even():
    """there used to be a bug when pprinting sequences with even height"""
    u = upretty([x**2])
    s = \
u"""\
⎡ 2⎤
⎣x ⎦\
"""
    assert u == s

    u = upretty((x**2,))
    s = \
u"""\
⎛ 2 ⎞
⎝x ,⎠\
"""
    assert u == s

    u = upretty({x**2: 1})
    s = \
u"""\
⎧ 2   ⎫
⎨x : 1⎬
⎩     ⎭\
"""
    assert u == s

