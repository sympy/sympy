# -*- coding: utf-8 -*-

from sympy import symbols, Symbol, sin, cos, Matrix, Integral, pi, sqrt, \
        Function, Rational, tan, oo
from sympy.printing.pretty import pretty

x,y = symbols('xy')
th  = Symbol('theta')
ph  = Symbol('phi')

def upretty(expr):
    return pretty(expr, True)

def test_upretty_greek():
    assert upretty( oo ) == u'∞'
    assert upretty( Symbol('alpha^+_1') )   ==  u'α⁺₁'
    assert upretty( Symbol('beta') )    == u'β'

def test_upretty_multiindex():
    assert upretty( Symbol('beta12') )  == u'β₁₂'
    assert upretty( Symbol('Y00') )     == u'Y₀₀'
    assert upretty( Symbol('Y_00') )    == u'Y₀₀'
    assert upretty( Symbol('F^+-') )    == u'F⁺⁻'

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


def test_upretty_nicerat():
    u = upretty(y*x**-2)
    s = \
u"""\
y 
──
 2
x \
"""
    assert u == s

    u = upretty(y**Rational(3,2) * x**Rational(-5,2))
    s = \
u"""\
 3/2
y   
────
 5/2
x   \
"""
    assert u == s

    u = upretty(sin(x)**3/tan(x)**2)
    s = \
u"""\
   3   
sin (x)
───────
   2   
tan (x)\
"""

def test_upretty_funcbraces():
    f = Function('f')
    u = upretty(f(x/(y+1), y))
    s = \
u"""\
 ⎛  x     ⎞
f⎜─────, y⎟
 ⎝1 + y   ⎠\
"""
    assert u == s

def test_upretty_sqrt():
    u = upretty( sqrt((sqrt(x+1))+1) )
    s = \
u"""\
   ⎽⎽⎽⎽⎽⎽⎽⎽⎽⎽⎽⎽⎽⎽⎽
  ╱       ⎽⎽⎽⎽⎽⎽⎽ 
╲╱  1 + ╲╱ 1 + x  \
"""
    assert u == s

def test_upretty_integral():
    u = upretty( Integral(sin(th)/cos(ph), (th,0,pi), (ph, 0, 2*pi)) )
    s = \
u"""\
2*π π             
 ⌠  ⌠             
 ⎮  ⎮ sin(θ)      
 ⎮  ⎮ ────── dθ dφ
 ⎮  ⎮ cos(φ)      
 ⌡  ⌡             
 0  0             \
"""
    assert u == s

    u = upretty( Integral(x**2*sin(y), (x,0,1), (y,0,pi)) )
    s = \
u"""\
π 1                
⌠ ⌠                
⎮ ⎮  2             
⎮ ⎮ x *sin(y) dx dy
⌡ ⌡                
0 0                \
"""
    assert u == s


def test_upretty_matrix():
    u = upretty( Matrix([[x**2+1, 1], [y, x+y]]) )
    s = \
u"""\
⎡     2       ⎤
⎢1 + x       1⎥
⎢             ⎥
⎣     y  x + y⎦\
"""
    assert u == s


def test_upretty_seq():
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
