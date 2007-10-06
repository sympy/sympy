from sympy import Basic, Symbol, Real, Rational, cos, exp, log, oo, sqrt, \
    symbols, Integral


def dotest(s):
    x = Symbol("x")
    y = Symbol("y")
    l = [
    Rational(2),
    Real("1.3"),
    x,
    y,
    pow(x,y)*y,
    5,
    5.5
    ]
    for x in l:
        for y in l:
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

def test_ldegree():
    x = Symbol("x")
    assert (1/x**2+1+x+x**2).ldegree(x)==-2
    assert (1/x+1+x+x**2).ldegree(x)==-1
    assert (x**2+1/x).ldegree(x)==-1
    assert (1+x**2).ldegree(x)==0
    assert (x+1).ldegree(x)==0
    assert (x+x**2).ldegree(x)==1
    assert (x**2).ldegree(x)==2

def test_leadterm():
    x = Symbol("x")
    assert (3+2*x**(log(3)/log(2)-1)).leadterm(x)==(3,0)

def test_atoms():
   x = Symbol('x')
   y = Symbol('y')
   assert list((1+x).atoms()) == [1,x]
   assert list(x.atoms()) == [x]
   assert list((1+2*cos(x)).atoms(type=Symbol)) == [x]
   assert list((2*(x**(y**x))).atoms()) == [2,x,y]
   assert list(Rational(1,2).atoms()) == [Rational(1,2)]
   assert list(Rational(1,2).atoms(type=type(oo))) == []

def test_is_polynomial():
    x, y, z = map(Symbol, 'xyz')

    k = Symbol('k', nonnegative=True, integer=True)

    assert Rational(2).is_polynomial(x, y, z) == True
    assert (Basic.Pi()).is_polynomial(x, y, z) == True

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
    assert (x**Basic.Half()).is_polynomial(x) == False
    assert (x**Rational(3,2)).is_polynomial(x) == False

    assert (x**2 + 3*x*sqrt(y) - 8).is_polynomial(x) == True
    assert (x**2 + 3*x*sqrt(y) - 8).is_polynomial(y) == False

    assert ((x**2)*(y**2) + x*(y**2) + y*x + exp(2)).is_polynomial() == True
    assert ((x**2)*(y**2) + x*(y**2) + y*x + exp(x)).is_polynomial() == False

    assert ((x**2)*(y**2) + x*(y**2) + y*x + exp(2)).is_polynomial(x, y) == True
    assert ((x**2)*(y**2) + x*(y**2) + y*x + exp(x)).is_polynomial(x, y) == False

def test_SAGE():
    #see http://code.google.com/p/sympy/issues/detail?id=247
    class MyInt:
        def __str__(self):
            return "5"
    m = MyInt()
    e = Rational(2)*m
    assert e == 10

def test_ordering():
    from sympy.core.methods import ArithMeths
    class CustomClass1(Basic, ArithMeths): pass
    class CustomClass2(Basic, ArithMeths): pass
    cc1 = CustomClass1(commutative=True)
    cc2 = CustomClass2(commutative=True)

    assert str(Rational(2)*cc1) == '2*CustomClass1()'
    assert str(cc1*Rational(2)) == '2*CustomClass1()'
    assert str(cc1*Real("1.5")) == '1.5*CustomClass1()'
    assert str(cc2*Rational(2)) == '2*CustomClass2()'
    assert str(cc2*Rational(2)*cc1) == '2*CustomClass1()*CustomClass2()'
    assert str(cc1*Rational(2)*cc2) == '2*CustomClass1()*CustomClass2()'

def test_len():
    x, y, z = symbols("xyz")
    e = x*y
    assert len(e) == 2
    e = x+y+z
    assert len(e) == 3

def test_doit():
    x = Symbol('x')
    a = Integral(x**2, x)

    assert isinstance(a.doit(), Integral) == False

    assert isinstance(a.doit(integrals=True), Integral) == False
    assert isinstance(a.doit(integrals=False), Integral) == True

    assert (2*Integral(x, x)).doit() == x**2

def test_is_number():
    x = Symbol("x")
    assert Rational(8).is_number
    assert not x.is_number
    assert (8+log(2)).is_number
    assert not (8+log(2)+x).is_number
    assert (1+x**2/x-x).is_number
