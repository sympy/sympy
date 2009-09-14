from sympy import Basic, S, Symbol, Real, Integer, Rational,  \
    sin, cos, exp, log, oo, sqrt, symbols, Integral, sympify, \
    WildFunction, Poly, Function, Derivative, Number, pi, var

from sympy.utilities.pytest import XFAIL, raises

class DummyNumber(object):
    """
    Minimal implementation of a number that works with SymPy.

    If one has a Number class (e.g. Sage Integer, or some other custom class)
    that one wants to work well with SymPy, one has to implement at least the
    methods of this class DummyNumber, resp. it's subclasses I5 and F1_1.

    Basically, one just needs to implement either __int__() or __float__() and
    then one needs to make sure that the class works with Python integers and
    with itself.
    """

    def __radd__(self, a):
        if isinstance(a, (int, float)):
            return a + self.number
        return NotImplemented

    def __truediv__(a, b):
        return a.__div__(b)

    def __rtruediv__(a, b):
        return a.__rdiv__(b)

    def __add__(self, a):
        if isinstance(a, (int, float, DummyNumber)):
            return self.number + a
        return NotImplemented

    def __rsub__(self, a):
        if isinstance(a, (int, float)):
            return a - self.number
        return NotImplemented

    def __sub__(self, a):
        if isinstance(a, (int, float, DummyNumber)):
            return self.number - a
        return NotImplemented

    def __rmul__(self, a):
        if isinstance(a, (int, float)):
            return a * self.number
        return NotImplemented

    def __mul__(self, a):
        if isinstance(a, (int, float, DummyNumber)):
            return self.number * a
        return NotImplemented

    def __rdiv__(self, a):
        if isinstance(a, (int, float)):
            return a / self.number
        return NotImplemented

    def __div__(self, a):
        if isinstance(a, (int, float, DummyNumber)):
            return self.number / a
        return NotImplemented

    def __rpow__(self, a):
        if isinstance(a, (int, float)):
            return a ** self.number
        return NotImplemented

    def __pow__(self, a):
        if isinstance(a, (int, float, DummyNumber)):
            return self.number ** a
        return NotImplemented

    def __pos__(self):
        return self.number

    def __neg__(self):
        return - self.number

class I5(DummyNumber):
    number = 5
    def __int__(self):
        return self.number

class F1_1(DummyNumber):
    number = 1.1
    def __float__(self):
        return self.number

x,y,z,t = symbols('xyzt')

i5 = I5()
f1_1 = F1_1()


# basic sympy objects
basic_objs = [
    Rational(2),
    Real("1.3"),
    x,
    y,
    pow(x,y)*y,
]

# all supported objects
all_objs = basic_objs + [
    5,
    5.5,
    i5,
    f1_1
]

def dotest(s):
    for x in all_objs:
        for y in all_objs:
            s(x,y)

def test_basic():
    def s(a,b):
        x = a
        x = +a
        x = -a
        x = a+b
        x = a-b
        x = a*b
        x = a/b
        x = a**b
    dotest(s)

def test_ibasic():
    def s(a,b):
        x = a
        x += b
        x = a
        x -= b
        x = a
        x *= b
        x = a
        x /= b
    dotest(s)

def test_basic_nostr():
    for obj in basic_objs:
        for op in ['+','-','*','/','**']:
            raises(TypeError, "obj %s '1'" % op)

def test_leadterm():
    assert (3+2*x**(log(3)/log(2)-1)).leadterm(x) == (3,0)

    assert (1/x**2+1+x+x**2).leadterm(x)[1] == -2
    assert (1/x+1+x+x**2).leadterm(x)[1] == -1
    assert (x**2+1/x).leadterm(x)[1] == -1
    assert (1+x**2).leadterm(x)[1] == 0
    assert (x+1).leadterm(x)[1] == 0
    assert (x+x**2).leadterm(x)[1] == 1
    assert (x**2).leadterm(x)[1] == 2

def test_as_leading_term():
    assert (3+2*x**(log(3)/log(2)-1)).as_leading_term(x) == 3
    assert (1/x**2+1+x+x**2).as_leading_term(x) == 1/x**2
    assert (1/x+1+x+x**2).as_leading_term(x) == 1/x
    assert (x**2+1/x).as_leading_term(x) == 1/x
    assert (1+x**2).as_leading_term(x) == 1
    assert (x+1).as_leading_term(x) == 1
    assert (x+x**2).as_leading_term(x) == x
    assert (x**2).as_leading_term(x) == x**2

def test_leadterm2():
    assert (x*cos(1)*cos(1 + sin(1)) + sin(1 + sin(1))).leadterm(x) == \
            (sin(1 + sin(1)), 0)

def test_leadterm3():
    assert (y+z+x).leadterm(x) == (y+z, 0)

def test_as_leading_term2():
    assert (x*cos(1)*cos(1 + sin(1)) + sin(1 + sin(1))).as_leading_term(x) == \
            sin(1 + sin(1))

def test_as_leading_term3():
    assert (2+pi+x).as_leading_term(x) == 2 + pi
    assert (2*x+pi*x+x**2).as_leading_term(x) == 2*x + pi*x

def test_atoms():
   assert sorted(list(x.atoms())) == [x]
   assert sorted(list((1+x).atoms())) == sorted([1, x])

   assert sorted(list((1+2*cos(x)).atoms(Symbol))) == [x]
   assert sorted(list((1+2*cos(x)).atoms(Symbol,Number))) == sorted([1, 2, x])

   assert sorted(list((2*(x**(y**x))).atoms())) == sorted([2, x, y])

   assert sorted(list(Rational(1,2).atoms())) == [S.Half]
   assert sorted(list(Rational(1,2).atoms(Symbol))) == []

   assert sorted(list(sin(oo).atoms(oo))) == [oo]

   assert sorted(list(Poly(0, x).atoms())) == [S.Zero]
   assert sorted(list(Poly(1, x).atoms())) == [S.One]
   assert sorted(list(Poly(x, x).atoms())) == sorted([S.One, x])
   assert sorted(list(Poly(x, x, y).atoms())) == sorted([S.One, x])
   assert sorted(list(Poly(x + y, x, y).atoms())) == sorted([S.One, x, y])
   assert sorted(list(Poly(x + y, x, y, z).atoms())) == sorted([S.One, x, y])
   assert sorted(list(Poly(x + y*t, x, y, z).atoms())) == \
           sorted([S.One, t, x, y])

def test_is_polynomial():
    z = Symbol('z')

    k = Symbol('k', nonnegative=True, integer=True)

    assert Rational(2).is_polynomial(x, y, z) == True
    assert (S.Pi).is_polynomial(x, y, z) == True

    assert x.is_polynomial(x) == True
    assert x.is_polynomial(y) == True

    assert (x**2).is_polynomial(x) == True
    assert (x**2).is_polynomial(y) == True

    assert (x**(-2)).is_polynomial(x) == False
    assert (x**(-2)).is_polynomial(y) == True

    assert (2**x).is_polynomial(x) == False
    assert (2**x).is_polynomial(y) == True

    assert (x**k).is_polynomial(x) == True
    assert (x**k).is_polynomial(k) == False
    assert (x**x).is_polynomial(x) == False
    assert (k**k).is_polynomial(k) == False
    assert (k**x).is_polynomial(k) == None

    assert (x**(-k)).is_polynomial(x) == None
    assert ((2*x)**k).is_polynomial(x) == True

    assert (x**2 + 3*x - 8).is_polynomial(x) == True
    assert (x**2 + 3*x - 8).is_polynomial(y) == True

    assert (x**2 + 3*x - 8).is_polynomial() == True

    assert sqrt(x).is_polynomial(x) == False
    assert (x**S.Half).is_polynomial(x) == False
    assert (x**Rational(3,2)).is_polynomial(x) == False

    assert (x**2 + 3*x*sqrt(y) - 8).is_polynomial(x) == True
    assert (x**2 + 3*x*sqrt(y) - 8).is_polynomial(y) == False

    assert ((x**2)*(y**2) + x*(y**2) + y*x + exp(2)).is_polynomial() == True
    assert ((x**2)*(y**2) + x*(y**2) + y*x + exp(x)).is_polynomial() == False

    assert ((x**2)*(y**2) + x*(y**2) + y*x + exp(2)).is_polynomial(x, y) == True
    assert ((x**2)*(y**2) + x*(y**2) + y*x + exp(x)).is_polynomial(x, y) == False

def test_is_rational_function():
    x,y = symbols('xy')

    assert Integer(1).is_rational_function() == True
    assert Integer(1).is_rational_function(x) == True

    assert Rational(17,54).is_rational_function() == True
    assert Rational(17,54).is_rational_function(x) == True

    assert (12/x).is_rational_function() == True
    assert (12/x).is_rational_function(x) == True

    assert (x/y).is_rational_function() == True
    assert (x/y).is_rational_function(x) == True
    assert (x/y).is_rational_function(x, y) == True

    assert (x**2+1/x/y).is_rational_function() == True
    assert (x**2+1/x/y).is_rational_function(x) == True
    assert (x**2+1/x/y).is_rational_function(x, y) == True

    assert (sin(y)/x).is_rational_function() == False
    assert (sin(y)/x).is_rational_function(y) == False
    assert (sin(y)/x).is_rational_function(x) == True
    assert (sin(y)/x).is_rational_function(x, y) == False

def test_SAGE1():
    #see http://code.google.com/p/sympy/issues/detail?id=247
    class MyInt:
        def _sympy_(self):
            return Integer(5)
    m = MyInt()
    e = Rational(2)*m
    assert e == 10

    raises(TypeError, "Rational(2)*MyInt")

def test_SAGE2():
    class MyInt(object):
        def __int__(self):
            return 5
    assert sympify(MyInt()) == 5
    e = Rational(2)*MyInt()
    assert e == 10

    raises(TypeError, "Rational(2)*MyInt")

def test_SAGE3():
    class MySymbol:
        def __rmul__(self, other):
            return ('mys', other, self)

    o = MySymbol()
    e = x*o

    assert e == ('mys', x, o)

def test_len():
    x, y, z = symbols("xyz")
    e = x*y
    assert len(e.args) == 2
    e = x+y+z
    assert len(e.args) == 3

def test_doit():
    a = Integral(x**2, x)

    assert isinstance(a.doit(), Integral) == False

    assert isinstance(a.doit(integrals=True), Integral) == False
    assert isinstance(a.doit(integrals=False), Integral) == True

    assert (2*Integral(x, x)).doit() == x**2

def test_attribute_error():
    raises(AttributeError, "x.cos()")
    raises(AttributeError, "x.sin()")
    raises(AttributeError, "x.exp()")

def test_args():
    assert (x*y).args[:] in ((x, y), (y, x))
    assert (x+y).args[:] in ((x, y), (y, x))
    assert (x*y+1).args[:] in ((x*y, 1), (1, x*y))
    assert sin(x*y).args[:] == (x*y,)
    assert sin(x*y).args[0] == x*y
    assert (x**y).args[:] == (x,y)
    assert (x**y).args[0] == x
    assert (x**y).args[1] == y

def test_iter_basic_args():
    assert list(sin(x*y).iter_basic_args()) == [x*y]
    assert list((x**y).iter_basic_args()) == [x, y]

    assert list(Poly(0, x).iter_basic_args()) == [S.Zero]
    assert list(Poly(1, x).iter_basic_args()) == [S.One]
    assert list(Poly(x, x).iter_basic_args()) == [S.One, x]

def test_noncommutative_expand_issue658():
    A, B, C = symbols('ABC', commutative=False)
    assert A*B - B*A != 0
    assert (A*(A+B)*B).expand() == A**2*B + A*B**2
    assert (A*(A+B+C)*B).expand() == A**2*B + A*B**2 + A*C*B

def test_as_independent():
    assert (2*x*sin(x)+y+x).as_independent(x) == (y, x + 2*x*sin(x))
    assert (2*x*sin(x)+y+x).as_independent(y) == (x + 2*x*sin(x), y)

    assert (2*x*sin(x)+y+x).as_independent(x, y) == (0, y + x + 2*x*sin(x))

    assert (x*sin(x)*cos(y)).as_independent(x) == (cos(y), x*sin(x))
    assert (x*sin(x)*cos(y)).as_independent(y) == (x*sin(x), cos(y))

    assert (x*sin(x)*cos(y)).as_independent(x, y) == (1, x*sin(x)*cos(y))

    assert (sin(x)).as_independent(x) == (1, sin(x))
    assert (sin(x)).as_independent(y) == (sin(x), 1)

    assert (2*sin(x)).as_independent(x) == (2, sin(x))
    assert (2*sin(x)).as_independent(y) == (2*sin(x), 1)

def test_subs_dict():
    a,b,c,d,e = symbols('abcde')

    assert (sin(x))._subs_dict({ x : 1, sin(x) : 2}) == 2
    assert (sin(x))._subs_dict([(x, 1), (sin(x), 2)]) == 2

    expr = sqrt(sin(2*x))*sin(exp(x)*x)*cos(2*x) + sin(2*x)

    seq = [ (sqrt(sin(2*x)),a), (cos(2*x),b), (sin(2*x),c), (x,d), (exp(x),e) ]
    assert expr._subs_dict(seq) == c + a*b*sin(d*e)

    seq = [ (sqrt(sin(2*x)),a), (sin(2*x),c), (cos(2*x),b), (x,d), (exp(x),e) ]
    assert expr._subs_dict(seq) == c + a*b*sin(d*e)

def test_subs_list():
    x,y = symbols('xy')

    assert (sin(x))._subs_list([(sin(x), 2), (x, 1)]) == 2
    assert (sin(x))._subs_list([(x, 1), (sin(x), 2)]) == sin(1)

    assert (x+y)._subs_list([(x, 3), (y, x**2)]) == 3 + x**2
    assert (x+y)._subs_list([(y, x**2), (x, 3)]) == 12

def test_call():
    a,b,c,d,e = symbols('abcde')

    assert sin(x)({ x : 1, sin(x) : 2}) == 2

    expr = sqrt(sin(2*x))*sin(exp(x)*x)*cos(2*x) + sin(2*x)

    assert expr({ sqrt(sin(2*x)) : a, cos(2*x) : b, sin(2*x) : c, x : d, exp(x) : e}) == c + a*b*sin(d*e)

def test_has():
    x, y = symbols("xy")
    f = Function("f")
    g = Function("g")
    assert sin(x).has(x)
    assert sin(x).has(sin)
    assert not sin(x).has(y)
    assert not sin(x).has(cos)
    assert f(x).has(x)
    assert f(x).has(f)
    assert not f(x).has(y)
    assert not f(x).has(g)

    assert f(x).diff(x).has(x)
    assert f(x).diff(x).has(f)
    assert f(x).diff(x).has(Derivative)
    assert not f(x).diff(x).has(y)
    assert not f(x).diff(x).has(g)
    assert not f(x).diff(x).has(sin)

def test_has_any_symbols():
    x,y,z,t,u = symbols('xyztu')

    i = Integer(4400)

    assert i.has_any_symbols(x) == False

    assert (i*x**i).has_any_symbols(x) == True
    assert (i*y**i).has_any_symbols(x) == False
    assert (i*y**i).has_any_symbols(x, y) == True

    expr = x**2*y + sin(2**t + log(z))

    assert expr.has_any_symbols(u) == False

    assert expr.has_any_symbols(x) == True
    assert expr.has_any_symbols(y) == True
    assert expr.has_any_symbols(z) == True
    assert expr.has_any_symbols(t) == True

    assert expr.has_any_symbols(x, y, z, t) == True
    assert expr.has_any_symbols(x, y, z, t, u)  == True

    from sympy.physics.units import m, s

    assert (x*m/s).has_any_symbols(x) == True
    assert (x*m/s).has_all_symbols(x) == True

    assert (x*m/s).has_any_symbols(y, z) == False
    assert (x*m/s).has_all_symbols(x, y) == False

    poly = Poly(x**2 + x*y*sin(z), x, y, t)

    assert poly.has_any_symbols(x) == True
    assert poly.has_any_symbols(x, y, z) == True
    assert poly.has_any_symbols(x, y, z, t) == True

    assert poly.has_all_symbols(x, y, z) == True
    assert poly.has_all_symbols(x, y, z, t) == False

def test_has_all_symbols():
    x,y,z,t,u = symbols('xyztu')

    i = Integer(4400)

    assert i.has_all_symbols(x) == False

    assert (i*x**i).has_all_symbols(x) == True
    assert (i*y**i).has_all_symbols(x) == False

    expr = x**2*y + sin(2**t + log(z))

    assert expr.has_all_symbols(y, z, t) == True
    assert expr.has_all_symbols(x, z, t) == True
    assert expr.has_all_symbols(x, y, t) == True
    assert expr.has_all_symbols(x, y, z) == True

    assert expr.has_all_symbols(y, u, t) == False
    assert expr.has_all_symbols(x, z, u) == False
    assert expr.has_all_symbols(u, y, z) == False

    assert expr.has_all_symbols(x, y, z, t) == True
    assert expr.has_all_symbols(x, y, z, t, u) == False

def test_as_poly_basic():
    x, y = symbols('xy')

    f = x**2 + 2*x*y

    assert f.as_poly().as_basic() == f
    assert f.as_poly(x, y).as_basic() == f

    assert (f + sin(x)).as_poly(x, y) is None

    p = Poly(f, x, y)

    assert p.as_poly() == p

def test_nonzero():
    assert bool(S.Zero) == False
    assert bool(S.One)  == True
    assert bool(x)      == True
    assert bool(x+y)    == True
    assert bool(x-x)    == False
    assert bool(x*y)    == True
    assert bool(x*1)    == True
    assert bool(x*0)    == False

def test_is_number():
    x, y = symbols('xy')
    g = WildFunction('g')

    assert Real(3.14).is_number == True
    assert Integer(737).is_number == True
    assert Rational(3, 2).is_number == True
    assert Rational(8).is_number == True
    assert x.is_number == False
    assert (2*x).is_number == False
    assert (x + y).is_number == False
    assert log(2).is_number == True
    assert log(x).is_number == False
    assert (2 + log(2)).is_number == True
    assert (8+log(2)).is_number == True
    assert (2 + log(x)).is_number == False
    assert (8+log(2)+x).is_number == False
    assert (2*g).is_number == False
    assert (1+x**2/x-x).is_number == True

    # test extensibility of .is_number
    # on subinstances of Basic
    class A(Basic):
        pass
    a = A()
    assert a.is_number == False


# TODO write more tests for as_coeff_factors
def test_as_coeff_factors():
    x = Symbol('x')

    assert     x .as_coeff_factors() == ( 0, (x,))
    assert (-1+x).as_coeff_factors() == (-1, (x,))
    assert ( 2+x).as_coeff_factors() == ( 2, (x,))
    assert ( 1+x).as_coeff_factors() == ( 1, (x,))


def test_as_coeff_exponent():
    x, y = symbols("xy")
    assert (3*x**4).as_coeff_exponent(x) == (3, 4)
    assert (2*x**3).as_coeff_exponent(x) == (2, 3)
    assert (4*x**2).as_coeff_exponent(x) == (4, 2)
    assert (6*x**1).as_coeff_exponent(x) == (6, 1)
    assert (3*x**0).as_coeff_exponent(x) == (3, 0)
    assert (2*x**0).as_coeff_exponent(x) == (2, 0)
    assert (1*x**0).as_coeff_exponent(x) == (1, 0)
    assert (0*x**0).as_coeff_exponent(x) == (0, 0)
    assert (-1*x**0).as_coeff_exponent(x) == (-1, 0)
    assert (-2*x**0).as_coeff_exponent(x) == (-2, 0)
    assert (2*x**3+pi*x**3).as_coeff_exponent(x) == (2+pi, 3)
    assert (x*log(2)/(2*x + pi*x)).as_coeff_exponent(x) == \
            (log(2)/(2+pi), 0)

def test_extractions():
    x, y = symbols("xy")
    n = Symbol("n", integer=True)
    assert ((x*y)**3).extract_multiplicatively(x**2 * y) == x*y**2
    assert ((x*y)**3).extract_multiplicatively(x**4 * y) == None
    assert (2*x).extract_multiplicatively(2) == x
    assert (2*x).extract_multiplicatively(3) == None
    assert (2*x).extract_multiplicatively(-1) == None
    assert (Rational(1,2)*x).extract_multiplicatively(3) == x/6
    assert (x**(Rational(1,2))).extract_multiplicatively(x) == None
    assert (x**(Rational(1,2))).extract_multiplicatively(1/x) == x**(Rational(3,2))

    assert ((x*y)**3).extract_additively(1) == None
    assert (x+1).extract_additively(x) == 1
    assert (x+1).extract_additively(2*x) == None
    assert (x+1).extract_additively(-x) == 1+2*x
    assert (-x+1).extract_additively(2*x) == 1-3*x

    assert (Integer(-3)).could_extract_minus_sign() == True
    assert (-n*x+x).could_extract_minus_sign() != (n*x-x).could_extract_minus_sign()
    assert (x-y).could_extract_minus_sign() != (-x+y).could_extract_minus_sign()
    assert (1-x-y).could_extract_minus_sign() == True
    assert (1-x+y).could_extract_minus_sign() == False
    assert ((-x-x*y)/y).could_extract_minus_sign() == True
    assert (-(x+x*y)/y).could_extract_minus_sign() ==  True
    assert ((x+x*y)/(-y)).could_extract_minus_sign() == True
    assert ((x+x*y)/y).could_extract_minus_sign() == False

def test_coeff():
    assert (3+2*x+4*x**2).coeff(1) == None
    assert (3+2*x+4*x**2).coeff(-1) == None
    assert (3+2*x+4*x**2).coeff(x) == 2
    assert (3+2*x+4*x**2).coeff(x**2) == 4
    assert (3+2*x+4*x**2).coeff(x**3) == None

    assert (-x/8 + x*y).coeff(x) == -S(1)/8 + y
    assert (-x/8 + x*y).coeff(-x) == S(1)/8 - y
    assert (-x/8 + x*y).coeff(2*x) == -S(1)/16 + y/2
    assert (x/8 + x*y).coeff(2*y*x) == S(1)/2
    assert (x/8 + x*y).coeff(y*x/2) == 2

    f = Function('f')
    assert (2*f(x) + 3*f(x).diff(x)).coeff(f(x)) == 2

def test_coeff2():
    var('r, kappa')
    psi = Function("psi")
    g = 1/r**2 * (2*r*psi(r).diff(r, 1) + r**2 * psi(r).diff(r, 2))
    g = g.expand()
    assert g.coeff((psi(r).diff(r))) == 2/r

def test_coeff2_0():
    var('r, kappa')
    psi = Function("psi")
    g = 1/r**2 * (2*r*psi(r).diff(r, 1) + r**2 * psi(r).diff(r, 2))
    g = g.expand()

    assert g.coeff(psi(r).diff(r, 2)) == 1

def test_coeff_expand():
    x, y, z = symbols('x y z')
    expr = z*(x+y)**2
    expr2 = z*(x+y)**2 + z*(2*x + 2*y)**2
    assert expr.coeff(z) == 2*x*y + x**2 + y**2
    assert expr.coeff(z, expand=False) == (x+y)**2
    assert expr2.coeff(z) == 10*x*y + 5*x**2 + 5*y**2
    assert expr2.coeff(z, expand=False) == (x+y)**2 + (2*x + 2*y)**2

def test_integrate():
    assert (log(x)).integrate((x, 0, 1)) == -1
    assert sin(x).integrate(x) == -cos(x)
    assert sin(x).integrate(('x',0,1)) == 1 - cos(1)

def test_count_ops():
    f = (x*y + 3/y)**(3 + 2)
    assert f.count_ops() == Symbol('ADD') + 2*Symbol('MUL') + 2*Symbol('POW')
    assert f.count_ops(symbolic=False) == 5
