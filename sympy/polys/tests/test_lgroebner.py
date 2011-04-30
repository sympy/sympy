"""Tests for groebner basis using lgroebner"""

from sympy.polys.lpoly import *
from sympy.polys.lgroebner import *
from sympy import *

def test_S_poly():
    lp = LPoly('x,y,z', QQ,O_lex)
    f1 = lp('x*y^2*z-x*y*z')
    f2 = lp('x^2*y*z - z^2')
    p3 = S_poly((f1.leading_expv(),f1),(f2.leading_expv(),f2))
    assert p3 == lp('-x^2*y*z + y*z^2')

def test_c1():
  # [CLO] p92
  lp,x,y = lgens('x,y',QQ,O_grlex)
  f3 = x**2
  f4 = x*y
  f5 = lp('y^2 - 1/2*x')
  gr = groebner_basis((f3,f4,f5))
  f3a = x**2 + x*y
  assert gr == [f3,f4,f5]

def test_c2():
    # [CLO] p94
    lp,x,y,z,w = lgens('x,y,z,w',QQ,O_lex)
    f0 = 3*x - 6*y - 2*z
    f1 = 2*x - 4*y + 4*w
    f2 = x - 2*y - z - w
    gr = groebner_basis((f0,f1,f2))
    assert gr == [x-2*y+2*w,z+3*w]

def test_c3():
    """from [CLO] p96"""
    lp = LPoly('x,y,z', QQ,O_lex)
    f0 = lp('x^2+y^2+z^2-1')
    f1 = lp('x^2+z^2-y')
    f2 = lp('x-z')
    gr = groebner_basis((f0,f1,f2))
    assert gr == [lp('+x -z'),lp('  +y -2*z^2'),lp('  +z^4 +1/2*z^2 -1/4')]
    " gr[2] = +z^4 +1/2*z^2 -1/4; solve gr[2] = 0 and back substitute"

def test_c4():
    """from [CLO] p97"""
    lp = LPoly('la,x,y,z', QQ,O_lex)
    f0 = lp('3*x^2 + 2*y*z - 2*x*la')
    f1 = lp('2*x*z-2*y*la')
    f2 = lp('2*x*y - 2*z - 2*z*la')
    f3 = lp('x^2+y^2+z^2-1')
    gr = groebner_basis((f0,f1,f2,f3))
    assert gr[-1] == lp('+z^7 -1763/1152*z^5 +655/1152*z^3 -11/288*z')
  
def test_c5():
    """from [CLO] p99"""
    lp,t,x,y,z = lgens('t,x,y,z',QQ,O_lex)
    f0 = t**4 - x
    f1 = t**3 - y
    f2 = t**2 - z
    gr = groebner_basis((f0,f1,f2))
    assert gr == [t**2 -z,t*y -z**2,t*z -y,x -z**2,y**2 -z**3]
  
def test_c6():
    """from [CLO] p100"""
    lp = LPoly('t,u,x,y,z',QQ,O_lex)
    f0 = lp('t + u - x')
    f1 = lp('t^2 + 2*t*u - y')
    f2 = lp('t^3 + 3*t^2*u - z')
    gr = groebner_basis((f0,f1,f2))
    assert gr == [lp('t +u -x'),lp('  +u^2 -x^2 +y'),lp('  +u*x^2 -u*y -x^3 +3/2*x*y -1/2*z'),lp('  +u*x*y -u*z -x^2*y -x*z +2*y^2'),lp('  +u*x*z -u*y^2 +x^2*z -1/2*x*y^2 -1/2*y*z'),lp('  +u*y^3 -u*z^2 -2*x^2*y*z +1/2*x*y^3 -x*z^2 +5/2*y^2*z'),lp('  +x^3*z -3/4*x^2*y^2 -3/2*x*y*z +y^3 +1/4*z^2')]
  
  
def test_s1():
    # http://www.sagemath.org/doc/constructions/polynomials.html
    lp,a,b,c,d = lgens('a,b,c,d',QQ,O_lex)
    I = (a + b + c + d, a*b + a*d + b*c + c*d, a*b*c + a*b*d + a*c*d + b*c*d, a*b*c*d - 1)
    gr = groebner_basis(I)
    assert gr == [lp('a +b +c +d'),lp('  +b^2 +2*b*d +d^2'),lp('  +b*c -b*d +c^2*d^4 +c*d -2*d^2'),lp('  +b*d^4 -b +d^5 -d'),lp('  +c^3*d^2 +c^2*d^3 -c -d'),lp('  +c^2*d^6 -c^2*d^2 -d^4 +1')]

def test_s2():
    # http://docs.sympy.org/modules/polys/wester.html
    lp,c,s = lgens('s,c', QQ,O_lex)
    f0 = (1 - c**2)**5 * (1 - s**2)**5 * (c**2 + s**2)**10
    f1 = lp('c^2 + s^2 - 1')
    gr = groebner_basis((f0,f1))
    assert gr == [lp('s^2 +c^2 -1'),lp('  +c^20 -5*c^18 +10*c^16 -10*c^14 +5*c^12 -c^10')]

def test_s2b():
    # http://docs.sympy.org/modules/polys/wester.html
    lp,c,s = lgens('s,c', QQ,O_lex)
    f0 = lp('1/7 - c^2')**5 * lp('1/7 - s^2')**5 * lp('c^2 + s^2')**10
    f1 = lp('c^2 + s^2 - 1/7')
    gr = groebner_basis((f0,f1))
    assert gr == [ lp('s^2 +c^2 -1/7'), lp('c^20 -5/7*c^18 +10/49*c^16 -10/343*c^14 +5/2401*c^12 -1/16807*c^10')]

def test_s3():
  """from talk 'Groebner basis and applications', by Hillar
  find the defining equation over the rationals for the algebraic number z
  defined as the solution to
  p(z)=z^5 + sqrt(2)*z^3 - a^2*z + a = 0
  where a is a solution to a^3 + a - 1 = 0
  z^30 -2*z^26 +2*z^22 +6*z^21 +2*z^20 -6*z^18 +14*z^17 +4*z^16 +2*z^15 -
  15*z^14 -26*z^13 -5*z^12 +22*z^11 -11*z^10 -8*z^9 -18*z^8 +5*z^6 +z^4 -
  2*z^3 +2*z^2+1
  """
  lp = LPoly('a,x,z',QQ,O_lex)
  f0 = lp('x^2 - 2')
  f1 = lp('a^3+a-1')
  f2=lp('z^5+x*z^3-a^2*z+a')
  gr = groebner_basis((f0,f1,f2))
  "the element of gr involving only z is the solution"

  assert gr[-1] == lp('z^30 -2*z^26 +2*z^22 +6*z^21 +2*z^20 -6*z^18 +14*z^17 +4*z^16 +2*z^15 -15*z^14 -26*z^13 -5*z^12 +22*z^11 -11*z^10 -8*z^9 -18*z^8 +5*z^6 +z^4 -2*z^3 +2*z^2 +1')

def mora(n, exp_bits):
    " [CLO] p112 "
    gens = ['x%d' % i for i in range(4)]
    #assert(n**2 + 1 < 2**exp_bits)
    lp = LPoly(gens,QQ, O_grlex)
    x0,x1,x2,x3 = lp.gens()
    I = [x0**(n+1) - x1*x2**(n-1)*x3, x0*x1**(n-1) - x2**n, x0**n*x2-x1**n*x3]
    gr = groebner_basis(I)
    return lp, gr

def test_mora():
    """ if n**2 + 1 > 2**exp_bits there is RPolyOverflowError
    """
    lp, gr = mora(3,6)
    x0,x1,x2,x3 = lp.gens()
    assert gr == [x0**4 -x1*x2**2*x3,x0**3*x2 -x1**3*x3,x0**2*x2**4 -x1**5*x3,
        x0*x1**2 -x2**3,x0*x2**7 -x1**7*x3,x1**9*x3 -x2**10]

def test_sy1():
    # from sympy mailing list
    # http://groups.google.com/group/sympy/browse_thread/thread/8b8dd87817bfae51/c5764132f0b99d6e
    lp,x,y,z = lgens('x,y,z', QQ,O_lex)
    I = [x**3+x+1, y**2+y+1, (x+y)*z-(x**2+y)]
    gr = groebner_basis(I)
    assert gr == [ lp('+x +155/2067*z^5 -355/689*z^4 +6062/2067*z^3 -3687/689*z^2 +6878/2067*z -25/53'),  lp('+y +4/53*z^5 -91/159*z^4 +523/159*z^3 -387/53*z^2 +1043/159*z -308/159'),  lp('+z^6 -7*z^5 +41*z^4 -82*z^3 +89*z^2 -46*z +13')]

