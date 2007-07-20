from sympy.core import *

x = Symbol('x')
y = Symbol('y')
z = Symbol('z')
zero = Basic.Zero()
o = Basic.One()


def test_simple_1():
    assert Order(2*x) == Order(x)
    assert Order(x)*3 == Order(x)
    assert -28*Order(x) == Order(x)
    assert Order(-23) == Order(1)
    assert Order(exp(x)),Order(1 == x)
    assert Order(exp(1/x)).expr == exp(1/x)
    assert Order(x*exp(1/x)).expr == x*exp(1/x)
    assert Order(x**(o/3)).expr == x**(o/3)
    assert Order(x**(5*o/3)).expr == x**(5*o/3)

def test_simple_2():
    assert Order(2*x)*x == Order(x**2)
    assert Order(2*x)/x,Order(1 == x)
    assert Order(2*x)*x*exp(1/x) == Order(x**2*exp(1/x))
    assert (Order(2*x)*x*exp(1/x)/ln(x)**3).expr == x**2*exp(1/x)*ln(x)**-3

def test_simple_3():
    assert Order(x)+x == Order(x)
    assert Order(x)+2 == 2+Order(x)
    assert Order(x)+x**2 == Order(x)
    assert Order(x)+1/x == 1/x+Order(x)
    assert Order(1/x)+1/x**2 == 1/x**2+Order(1/x)
    assert Order(x)+exp(1/x) == Order(x)+exp(1/x)
        
def test_simple_4():
    assert Order(x)**2 == Order(x**2)
    assert Order(x**3)**-2 == Order(x**-6)

def test_simple_5():
    assert Order(x)+Order(x**2) == Order(x)
    assert Order(x)+Order(x**-2) == Order(x**-2)
    assert Order(x)+Order(1/x) == Order(1/x)

def test_simple_6():
    assert Order(x)-Order(x) == Order(x)
    assert Order(x)+Order(1) == Order(1)
    assert Order(x)+Order(x**2) == Order(x)
    assert Order(1/x)+Order(1) == Order(1/x)
    assert Order(x)+Order(exp(1/x)) == Order(exp(1/x))
    assert Order(x**3)+Order(exp(2/x)) == Order(exp(2/x))
    assert Order(x**-3)+Order(exp(2/x)) == Order(exp(2/x))

def test_contains_0():
    assert Order(1,x).contains(Order(1,x))
    assert Order(1,x).contains(Order(1))
    assert Order(1).contains(Order(1,x))

def test_contains_1():
    assert Order(x).contains(Order(x))
    assert Order(x).contains(Order(x**2))
    assert not Order(x**2).contains(Order(x))
    assert not Order(x).contains(Order(1/x))
    assert not Order(1/x).contains(Order(exp(1/x)))
    assert not Order(x).contains(Order(exp(1/x)))
    assert Order(1/x).contains(Order(x))
    assert Order(exp(1/x)).contains(Order(x))
    assert Order(exp(1/x)).contains(Order(1/x))
    assert Order(exp(1/x)).contains(Order(exp(1/x)))
    assert Order(exp(2/x)).contains(Order(exp(1/x)))
    assert not Order(exp(1/x)).contains(Order(exp(2/x)))

def test_contains_2():
    assert Order(x).contains(Order(y)) is None
    assert Order(x).contains(Order(y*x))
    assert Order(y*x).contains(Order(x))
    assert Order(y).contains(Order(x*y))
    assert Order(x).contains(Order(y**2*x))

def test_contains_3():
    assert Order(x*y**2).contains(Order(x**2*y)) is None
    assert Order(x**2*y).contains(Order(x*y**2)) is None

def test_add_1():
    assert Order(x+x) == Order(x)
    assert Order(3*x-2*x**2) == Order(x)
    assert Order(1+x), Order(1 == x)
    assert Order(1+1/x) == Order(1/x)
    assert Order(ln(x)+1/ln(x)) == Order(ln(x))
    assert Order(exp(1/x)+x) == Order(exp(1/x))
    assert Order(exp(1/x)+1/x**20) == Order(exp(1/x))

def test_ln_args():
    assert Order(ln(2*x)).expr == ln(x) # ln(2*x) -> ln(2)+ln(x
    assert Order(ln(y*x)).expr == ln(x)+ln(y) # ln(x*y) -> ln(x)+ln(y
    assert Order(ln(x**3)).expr == ln(x) # ln(x**3) -> 3*ln(x
    assert Order(ln(2*x)).expr == ln(x) # ln(2*x) -> ln(2)+ln(x

def test_multivar_0():
    assert Order(x*y).expr == x*y
    assert Order(x*y**2).expr == x*y**2
    assert Order(x*y,x).expr == x
    assert Order(x*y**2,y).expr == y**2
    assert Order(x*y*z).expr == x*y*z
    assert Order(x/y).expr == x/y
    assert Order(x*exp(1/y)).expr == x*exp(1/y)
    assert Order(exp(1/x)*exp(1/y)).expr == exp(1/x)*exp(1/y)
    assert Order(exp(x)*exp(1/y)).expr == exp(1/y)

def test_multivar_1():
    assert Order(x+y).expr == x+y
    assert Order(x+2*y).expr == x+y
    assert (Order(x+y)+x).expr == (x+y)
    assert (Order(x+y)+x**2) == Order(x+y)
    assert (Order(x+y)+1/x) == 1/x+Order(x+y)
    assert Order(x**2+y*x).expr == x**2+y*x
        
def test_multivar_2():
    assert Order(x**2*y+y**2*x,x,y).expr == x**2*y+y**2*x

def test_multivar_mul_1():
    assert Order(x+y)*x,Order(x**2+y*x,x == y)

def test_multivar_3():
    assert (Order(x)+Order(y))[:],(Order(x) == Order(y))
    assert Order(x)+Order(y)+Order(x+y) == Order(x+y)
    assert (Order(x**2*y)+Order(y**2*x))[:],(Order(x*y**2) == Order(y*x**2))
    assert (Order(x**2*y)+Order(y*x)) == Order(x*y)
    
def test_w():
    print
    for k,v in Order._cache.items():
        if isinstance(k, Basic.Symbol):
            print k,v
