
import sys
sys.path.append(".")

import py

from sympy import *

## sympy/modules/polynomials/base.py

def test_Polynomial():
    x = Symbol("x")
    y = Symbol("y")
    z = Symbol('z')
    f = Polynomial(x+2)
    g = Polynomial(y**2-1)
    h = f + g
    assert f.var == [x]
    assert f.coeffs == ((1, 1), (2, 0))
    assert str(f) == "2 + x"
    assert repr(f) == "Polynomial(2 + x, ((1, 1), (2, 0)), [x], 'grevlex')"
    assert f == 2 + x
    assert f == Polynomial(None, f.coeffs, f.var, f.order)
    assert f.nth_coeff(0) == 2
    assert f.nth_coeff(2) == 0
    assert h.var == [x, y]
    assert h.coeffs == ((1, 0, 2), (1, 1, 0), (1, 0, 0))
    h = f*Polynomial(y, var=x)
    assert h.var == [x]
    assert h.coeffs == ((y, 1), (2*y, 0))
    h = f*y
    assert h == (x + 2)*y
    assert not isinstance(h, Polynomial)
    h = Polynomial(h, var=y)
    assert h.var == [y]
    assert h.coeffs == ((2+x, 1),)
    assert Polynomial(1, var=x).diff(x) == Polynomial(0, var=x)
    assert Polynomial(x**3*y).diff(x) == Polynomial(3*x**2*y)
    assert Polynomial(coeffs=((Integer(1), Integer(1)),), var=x).diff(x) \
           == Polynomial(1, var=x)
    assert Polynomial(x**2 + y**2)(3,-4) == 25
    assert Polynomial(y*x)(-z, z) == -z**2

    #TODO better test that differs between all orders ?
    from sympy import sin
    assert Polynomial(1).coeffs == ((1,),)
    assert Polynomial(x).coeffs == ((1, 1),)
    assert Polynomial(x**2+y**3, order='lex').coeffs \
                   == ((1, 2, 0), (1, 0, 3))
    assert Polynomial(x**2+y**3, var=[y, x]).coeffs == ((1,3,0), (1,0,2))
    assert Polynomial(x*y).coeffs == ((1, 1, 1),)
    assert Polynomial(x**2*y**4 + sin(z)*x**3 + x*y**5,
                      var=[x, y], order='lex').coeffs \
           == ((sin(z), 3, 0), (1, 2, 4), (1, 1, 5))
    assert Polynomial(x**2*y**4 + sin(z)*x**3 + x*y**5,
                      var=[x, y], order='grlex').coeffs \
           == ((1, 2, 4), (1, 1, 5), (sin(z), 3, 0))
    assert Polynomial(x**2*y**4 + sin(z)*x**3 + x**5*y,
                      var=[x, y], order='grevlex').coeffs \
               == ((1, 5, 1), (1, 2, 4), (sin(z), 3, 0))
    assert Polynomial(z*x + x**2*y**2 + x**3*y,
                      var=[z, x, y], order='1-el').coeffs \
           == ((1, 1, 1, 0), (1, 0, 3, 1), (1, 0, 2, 2))

    py.test.raises(PolynomialException, "Polynomial(sqrt(x),var=x)")
    py.test.raises(PolynomialException, "Polynomial(sin(x),var=x)")

    assert 3*x**2 == Polynomial(coeffs=((Integer(3), Integer(2)),),
                                var=x).sympy_expr
    assert 2*x + 3*x**2 - 5 \
           == Polynomial(coeffs=((Integer(-5), Integer(0)),
                                 (Integer(2), Integer(1)),
                                 (Integer(3), Integer(2))),
                         var=[x]).sympy_expr
    assert 2*x**100 + 3*x**2 - 5 \
           == Polynomial(coeffs=((Integer(-5), Integer(0)),
                                 (Integer(3), Integer(2)),
                                 (Integer(2), Integer(100))),
                         var=[x]).sympy_expr

    assert sqrt(y)*x == Polynomial(coeffs=((sqrt(y), Integer(1)),),
                                   var=[x]).sympy_expr
    p = Polynomial(x/3 + 12*y + x**2/8)
    assert p.as_integer() == (24, Polynomial(3*x**2 + 8*x + 288*y))
    assert p.as_monic() == (Rational(1,8), Polynomial(x**2 + 96*y + 8*x/3))
    assert p.as_primitive() == (0, p)
    p = Polynomial(100*x + 12*y + 8*x**2)
    assert p.as_primitive() == (4, Polynomial(2*x**2 + 3*y + 25*x))
    assert p.leading_coeff() == Rational(8)
    assert p.leading_term() == Polynomial(8*x**2)

## sympy/modules/polynomials/common.py

def test_coeff_ring():
    from sympy.polynomials.common import coeff_ring
    x = Symbol("x")
    assert coeff_ring([Rational(2)]) == 'int'
    assert coeff_ring([Rational(2), Rational(1,2)]) == 'rat'
    assert coeff_ring([Rational(2)**Rational(1,2)]) == 'real'
    assert coeff_ring([Pi]) == 'real'
    assert coeff_ring([Real(2.1), Rational(-1)**Rational(1,2)]) == 'cplx'
    assert coeff_ring([I, x]) == 'sym'

## sympy/modules/polynomials/wrapper.py

def test_div():
    x = Symbol("x")
    y = Symbol('y')

    assert div(x**3-12*x**2-42, x-3, x) == (x**2-9*x-27, -123)
    assert div(x**3-12*x**2-42, x**2+x-3, x) == (x-13, 16*x-81)

    assert div(2+2*x+x**2, 1, x) == (2+2*x+x**2, 0)
    assert div(2+2*x+x**2, 2, x, coeff='int') == (1+x, x**2)

    assert div(3*x**3, x**2, x) == (3*x, 0)

    assert div(1,1) == (1, 0)
    assert div(1,x,x) == (0, 1)
    assert div(x*y+2*x+y,x,x) == (2+y, y)
    assert div(x*y+2*x+y,x,y) == (2+(1+1/x)*y, 0)

    assert div(x*y**2 + 1, [x*y+1, y+1], [x,y]) == ([y, -1], 2)
    assert div(x**2*y+x*y**2+y**2, [x*y-1, y**2-1], [x, y]) \
           == ([x+y, 1], 1+x+y)
    assert div(x**2*y+x*y**2+y**2, [y**2-1, x*y-1], [x, y]) \
           == ([1+x, x], 1+2*x)

def test_factor():
    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")
    assert factor(Rational(3, 8)*x**2 - Rational(3, 2)) \
           == Rational(3, 8)*((x + 2)*(x - 2))
    assert factor(x**3-1) == (x-1)*(x**2+x+1)
    assert factor(x**2+2*x+1) == (x+1)**2
    assert factor(x**3-3*x**2+3*x-1) == (x-1)**3
    assert factor(x**2+x-2) == (x-1)*(x+2)
    assert factor(x**3-x) == x*(x-1)*(x+1)
    assert factor(x**6-1) == (1+x**2-x)*(1+x)*(1+x+x**2)*(-1+x)
    assert factor(2*x**2+5*x+2) == (2+x)*(1+2*x)

    assert factor(x**2 + y**2) == x**2 + y**2
    assert factor(x*y + x*z + y*z) == x*y + x*z + y*z
    assert factor(x*(y+1) + x*z) == x*(z + y + 1)
    
def test_gcd():
    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")
    
    assert gcd(x**2, x, x) == x
    assert gcd(3*x**2, x, x) == x
    assert gcd(3*x**2, 6*x, x, coeff='rat') == x
    assert gcd(3*x**2, 6*x, x) == 3*x
    assert gcd(x**2+2*x+1, x+1, x) == x+1
    assert gcd(x**2+2*x+2, x+1, x) == 1

    assert gcd(x**2+2*x+1, 2+2*x, x) == 1+x
    assert gcd(x**2+2*x+2, 2+2*x, x) == 1

    assert gcd(4, 6) == Rational(2)
    assert gcd(6, 4, coeff='rat') == Rational(1)
    assert gcd(x, y) == Rational(1)
    assert gcd(sin(z)*(x+y), x**2+2*x*y+y**2, [x, y]) == x+y

def test_groebner():
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')

    assert groebner(y*x, x) == [x]
    assert groebner(y*x, x, reduced=False) == [x*y]
    assert groebner(x*y, z) == [1]
    
    # This one already is a Groebner base.
    assert groebner([y-x**2, z-x**3], [y,z,x], 'lex', reduced=False) \
           == [-x**2+y, z-x**3]

    assert groebner([x**3-2*x*y, x**2*y-2*y**2+x], [x,y], 'grlex',
                    reduced=False) \
           == [x**3-2*x*y, x+x**2*y-2*y**2, -x**2, -2*x*y, -2*y**2+x]
    assert groebner([x**3-2*x*y, x**2*y-2*y**2+x], [x,y], 'grlex') \
           == [x**2, x*y, Rational(-1,2)*x+y**2]

def test_lcm():
    x = Symbol('x')
    y = Symbol('y')

    assert lcm(6, 4) == Rational(12)
    assert lcm(6, 4, coeff='rat') == Rational(1)
    assert lcm(4, y) == 4*y
    assert lcm(x, y) == x*y
    assert lcm(y*(x+1), x, x) ==x+x**2
    assert lcm(2*x, x**2) == 2*x**2 

def test_real_roots():
    x = Symbol('x')

    f = x-1
    assert real_roots(f) == 1
    assert real_roots(f, None, Rational(0)) == 0
    assert real_roots(f, Rational(0), Rational(1)) == 1
    assert real_roots(f, Rational(1), None) == 0
    f = x**2 - 4
    assert real_roots(f) == 2
    assert real_roots(f, None, Rational(0)) == 1
    assert real_roots(f, Rational(-1), Rational(1)) == 0

def test_resultant():
     x, a, b, c, = [Symbol(y) for y in ['x', 'a', 'b', 'c']]

     s_res = resultant(x**2-1, x**3-x**2+2, x, method='sylvester').expand()
     b_res = resultant(x**2-1, x**3-x**2+2, x, method='bezout').expand()

     assert b_res == s_res == 0

     s_res = resultant(3*x**3-x, 5*x**2+1, x, method='sylvester').expand()
     b_res = resultant(3*x**3-x, 5*x**2+1, x, method='bezout').expand()

     assert b_res == s_res == 64

     s_res = resultant(x**2-2*x+7, x**3-x+5, x, method='sylvester').expand()
     b_res = resultant(x**2-2*x+7, x**3-x+5, x, method='bezout').expand()

     assert b_res == s_res == 265

     s_res = resultant((x-a)**2-2, a**2-3, a, method='sylvester').expand()
     b_res = resultant((x-a)**2-2, a**2-3, a, method='bezout').expand()

     assert b_res == s_res == 1 - 10*x**2 + x**4

     s_res = resultant((x-1)*(x-2)*(x-3), (x-4)*(x-5)*(x-6), x, method='sylvester').expand()
     b_res = resultant((x-1)*(x-2)*(x-3), (x-4)*(x-5)*(x-6), x, method='bezout').expand()

     assert b_res == s_res == -8640

     s_res = resultant((x-1)*(x-2)*(x-3), (x-4)*(x-5)*(x-1), x, method='sylvester').expand()
     b_res = resultant((x-1)*(x-2)*(x-3), (x-4)*(x-5)*(x-1), x, method='bezout').expand()

     assert b_res == s_res == 0

     s_res = resultant(x**3-1, x**3+2*x**2+2*x-1, x, method='sylvester').expand()
     b_res = resultant(x**3-1, x**3+2*x**2+2*x-1, x, method='bezout').expand()

     assert b_res == s_res == 16

     s_res = resultant(x**8-2, x-1, x, method='sylvester').expand()
     b_res = resultant(x**8-2, x-1, x, method='bezout').expand()

     assert b_res == s_res == -1

     s_res = resultant(3*x**2+2*a*x+3*a**2-2, 3*x**2-2*a*x+3*a**2-2, x, method='sylvester').expand()
     b_res = resultant(3*x**2+2*a*x+3*a**2-2, 3*x**2-2*a*x+3*a**2-2, x, method='bezout').expand()

     assert b_res == s_res == 144*a**4 - 96*a**2

     s_res = resultant((x-a)*(x-b), x-c, x, method='sylvester').expand()
     b_res = resultant((x-a)*(x-b), x-c, x, method='bezout').expand()

     assert b_res == s_res == ((a-c)*(b-c)).expand()

def test_roots():
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    x = Symbol("x")
    assert roots(x**2-3*x+2) == [1,2]  
    assert roots(x**2-3*x/2+Rational(1,2)) == [1,Rational(1,2)]
    assert roots(2*x**2-3*x+1) == [1,Rational(1,2)]
    assert roots(x**2-1) == [-1, 1]
    assert roots(x**2+1) == [-I, I]
    assert roots(x**3-1) == [-(-1)**Rational(1,3),
                             (-1)** Rational(1,3)/2
                             -(-1)**Rational(5,6)/2 *3**Rational(1,2),
                             (-1)**Rational(1,3)/2
                             +(-1)**Rational(5,6)/2*3**Rational(1,2)]
    assert roots(x**3) == [0]
    assert roots(x**3-x) == [0,-1,1]
    assert roots(Rational(2),x) == []
    assert roots(a*x**2 + b*x + c, var=[x]) == \
           [-b/(a*2)+(((b/a)**2-4*c/a)**Rational(1,2))/2,
            -b/(a*2)-(((b/a)**2-4*c/a)**Rational(1,2))/2]
    assert roots(x**3 + x**2 + x + 1) == [-1, -I, I]
    assert roots(x**4 - 1) == [1, exp(Pi*I/2), exp(Pi*I), exp(3*Pi*I/2)]
    assert roots(x**4 + 1) == [(-1)**Rational(1,4),
                               (-1)**Rational(1,4)*exp(Pi*I/2),
                               (-1)**Rational(1,4)*exp(Pi*I),
                               (-1)**Rational(1,4)*exp(3*Pi*I/2)]
    assert roots(x**5 - Rational(3,2)) == \
           [Rational(1,2)**Rational(1,5)*3**Rational(1,5),
            Rational(1,2)**Rational(1,5)*3**Rational(1,5)*exp(2*Pi*I/5),
            Rational(1,2)**Rational(1,5)*3**Rational(1,5)*exp(4*Pi*I/5),
            Rational(1,2)**Rational(1,5)*3**Rational(1,5)*exp(6*Pi*I/5),
            Rational(1,2)**Rational(1,5)*3**Rational(1,5)*exp(8*Pi*I/5)]

def test_solve_system():
    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")
    assert solve_system(x-1) == [[1]]
    assert solve_system([2*x - 3, 3*y/2 - 2*x, z - 5*y]) \
           == [[Rational(3, 2), 2, 10]]
    assert solve_system([y - x, y - x - 1]) == []
    assert solve_system([y - x**2, y + x**2]) == [[0, 0]]
    assert solve_system([y - x**2, y + x**2 + 1]) == \
           [[-I*2**Rational(1,2)/2, Rational(-1,2)],
            [I*2**Rational(1,2)/2, Rational(-1,2)]]
           
def test_sqf():
    x = Symbol("x")
    assert sqf(3*x**2, x) == 3*x**2
    assert sqf(x**2+2*x+1, x) == (x+1)**2
    assert sqf(x**5 - x**4 - x + 1) == (x-1)**2*(x**3 + x**2 + x + 1)

def test_sqf_part():
    x = Symbol('x')
    assert sqf_part(3*x**2, x) == 3*x
    assert sqf_part(x**2 + 2*x + 1, x) == x+1
    assert sqf_part(x**5 - x**4 - x + 1) == x**4 - 1

## sympy/modules/polynomials/ideals.py
    
## def test_Ideal():
##     x = Symbol('x')
##     y = Symbol('y')
##     z = Symbol('z')

##     # TODO: more complete tests?
##     assert len(Ideal()) == 1
##     assert not x in Ideal()
##     I = Ideal([x,y], [x,y])
##     assert x*y**2 in I
##     assert z*x in I
##     assert not z in I
##     assert I == I + Ideal()
##     assert Ideal() == Ideal()*I
##     assert I + I == I
##     assert I == I % Ideal()
##     assert Ideal() == Ideal(x*y, [x,y]) % I
##     assert Ideal(z, [x,y,z]) == Ideal([x,z], [x,y,z]) % Ideal([x,y], [x,y,z])

## sympy/modules/polynomials/roots_.py

def test_sturm():
    from sympy.polynomials import roots_
    x = Symbol('x')

    f = Polynomial(Rational(5), var=x)
    assert roots_.sturm(f) == [f]
    f = Polynomial(2*x)
    assert roots_.sturm(f) == [f, Polynomial(2, var=x)]
    f = Polynomial(x**3 - 2*x**2 + 3*x -5)
    assert roots_.sturm(f) == \
           [Polynomial(-5-2*x**2+x**3+3*x), Polynomial(3+3*x**2-4*x),
            Polynomial(Rational(13,3)-Rational(10,9)*x),
            Polynomial(Rational(-3303,100), var=x)]
