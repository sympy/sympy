from sympy import Basic, Symbol, Real, Rational, cos, exp, log, oo


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

def _test_print_tree():
    # XXX print_tree no longer exists
    x = Symbol("x")
    y = Symbol("y")
    e = (2*x-(7*x**2 - 2) + 3*y)
    e.print_tree()

def test_atoms():
   x = Symbol('x')
   y = Symbol('y')
   assert list((1+x).atoms()) == [1,x]
   assert list(x.atoms()) == [x]
   assert list((1+2*cos(x)).atoms(type=Symbol)) == [x]
   assert list((2*(x**(y**x))).atoms()) == [2,x,y]
   assert list(Rational(1,2).atoms()) == [Rational(1,2)]
   assert list(Rational(1,2).atoms(type=type(oo))) == []
