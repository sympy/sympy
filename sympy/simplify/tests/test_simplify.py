from sympy import Symbol, symbols, together, hypersimp, factorial, binomial, \
        collect, Function, powsimp, separate, sin, exp, Rational, fraction, \
        simplify, trigsimp, cos, tan, cot, log, ratsimp, Matrix, pi, integrate, \
        solve, nsimplify, GoldenRatio, sqrt, I, sympify, atan, Derivative

from sympy.utilities.pytest import XFAIL

def test_ratsimp():
    x = Symbol("x")
    y = Symbol("y")
    e = 1/x+1/y
    assert e != (x+y)/(x*y)
    assert ratsimp(e) == (x+y)/(x*y)

    e = 1/(1+1/x)
    assert ratsimp(e) == x/(x+1)
    assert ratsimp(exp(e)) == exp(x/(x+1))

def test_ratsimp2():
    x = Symbol("x")
    e = 1/(1+1/x)
    assert (x+1)*ratsimp(e)/x == 1

@XFAIL
def test_ratsimp_X1():
    e = -x-y-(x+y)**(-1)*y**2+(x+y)**(-1)*x**2
    assert e != -2*y
    assert ratsimp(e) == -2*y

@XFAIL
def test_ratsimp_X2():
    e = x/(x+y)+y/(x+y)
    assert e != 1
    assert ratsimp(e) == 1

def test_trigsimp():
    x,y = map(Symbol, 'xy')

    assert trigsimp(1 - sin(x)**2) == cos(x)**2
    assert trigsimp(1 - cos(x)**2) == sin(x)**2
    assert trigsimp(sin(x)**2 + cos(x)**2) == 1
    assert trigsimp(1 + tan(x)**2) == 1/cos(x)**2
    assert trigsimp(1/cos(x)**2 - 1) == tan(x)**2
    assert trigsimp(1/cos(x)**2 - tan(x)**2) == 1
    assert trigsimp(1 + cot(x)**2) == 1/sin(x)**2
    assert trigsimp(1/sin(x)**2 - 1) == cot(x)**2
    assert trigsimp(1/sin(x)**2 - cot(x)**2) == 1

    assert trigsimp(5*cos(x)**2 + 5*sin(x)**2) == 5
    assert trigsimp(5*cos(x/2)**2 + 2*sin(x/2)**2) in \
                [2 + 3*cos(x/2)**2, 5 - 3*sin(x/2)**2]

    assert trigsimp(cos(0.12345)**2 + sin(0.12345)**2) == 1
    e = 2*sin(x)**2 + 2*cos(x)**2
    assert trigsimp(log(e), deep=True) == log(2)

@XFAIL
def test_factorial_simplify():
    # There are more tests in test_factorials.py. These are just to
    # ensure that simplify() calls factorial_simplify correctly
    from sympy.specfun.factorials import factorial
    x = Symbol('x')
    assert simplify(factorial(x)/x) == factorial(x-1)
    assert simplify(factorial(factorial(x))) == factorial(factorial(x))

def test_simplify():
    x,y,z,k,n,m,w,f,s,A = symbols('xyzknmwfsA')

    e = 1/x + 1/y
    assert e != (x+y)/(x*y)
    assert simplify(e) == (x+y)/(x*y)

    e = A**2*s**4/(4*pi*k*m**3)
    assert simplify(e) == e

    e = (4+4*x-2*(2+2*x))/(2+2*x)
    assert simplify(e) == 0

    e = (-4*x*y**2-2*y**3-2*x**2*y)/(x+y)**2
    assert simplify(e) == -2*y

    e = -x-y-(x+y)**(-1)*y**2+(x+y)**(-1)*x**2
    assert simplify(e) == -2*y

    e = (x+x*y)/x
    assert simplify(e) == 1 + y

    e = (f(x)+y*f(x))/f(x)
    assert simplify(e) == 1 + y

    e = (2 * (1/n - cos(n * pi)/n))/pi
    assert simplify(e) == (2 - 2*cos(pi*n))/(pi*n)

    e = integrate(1/(x**3+1), x).diff(x)
    assert simplify(e) == 1/(x**3+1)

    e = integrate(x/(x**2+3*x+1), x).diff(x)
    assert simplify(e) == x/(x**2+3*x+1)

    A = Matrix([[2*k-m*w**2, -k],[-k,k-m*w**2]]).inv()

    assert simplify((A*Matrix([0,f]))[1]) == \
        (2*f*k - f*m*w**2)/(k**2 - 3*k*m*w**2 + m**2*w**4)

    a,b,c,d,e,f,g,h,i = symbols('abcdefghi')

    f_1 = x*a + y*b + z*c - 1
    f_2 = x*d + y*e + z*f - 1
    f_3 = x*g + y*h + z*i - 1

    solutions = solve([f_1,f_2,f_3], x,y,z, simplified=False)

    assert simplify(solutions[y]) == \
        (a*i+c*d+f*g-a*f-c*g-d*i)/(a*e*i+b*f*g+c*d*h-a*f*h-b*d*i-c*e*g)

def test_simplify_fail1():
    x = Symbol('x')
    y = Symbol('y')
    e = (x+y)**2/(-4*x*y**2-2*y**3-2*x**2*y)
    assert simplify(e) == 1 / (-2*y)

def test_fraction():
    x, y, z = map(Symbol, 'xyz')

    assert fraction(Rational(1, 2)) == (1, 2)

    assert fraction(x) == (x, 1)
    assert fraction(1/x) == (1, x)
    assert fraction(x/y) == (x, y)
    assert fraction(x/2) == (x, 2)

    assert fraction(x*y/z) == (x*y, z)
    assert fraction(x/(y*z)) == (x, y*z)

    assert fraction(1/y**2) == (1, y**2)
    assert fraction(x/y**2) == (x, y**2)

    assert fraction((x**2+1)/y) == (x**2+1, y)
    assert fraction(x*(y+1)/y**7) == (x*(y+1), y**7)

    assert fraction(exp(-x), exact=True) == (exp(-x), 1)

def test_together():
    x, y, z = map(Symbol, 'xyz')

    assert together(1/x) == 1/x

    assert together(1/x + 1) == (x+1)/x
    assert together(1/x + x) == (x**2+1)/x

    assert together(1/x + Rational(1, 2)) == (x+2)/(2*x)

    assert together(1/x + 2/y) == (2*x+y)/(y*x)
    assert together(1/(1 + 1/x)) == x/(1+x)
    assert together(x/(1 + 1/x)) == x**2/(1+x)

    assert together(1/x + 1/y + 1/z) == (x*y + x*z + y*z)/(x*y*z)

    assert together(1/(x*y) + 1/(x*y)**2) == y**(-2)*x**(-2)*(1+x*y)
    assert together(1/(x*y) + 1/(x*y)**4) == y**(-4)*x**(-4)*(1+x**3*y**3)
    assert together(1/(x**7*y) + 1/(x*y)**4) == y**(-4)*x**(-7)*(x**3+y**3)

    assert together(sin(1/x+1/y)) == sin(1/x+1/y)
    assert together(sin(1/x+1/y), deep=True) == sin((x+y)/(x*y))

    assert together(Rational(1,2) + x/2) == (x+1)/2

    assert together(1/x**y + 1/x**(y-1)) == x**(-y)*(1 + x)

def test_separate():
    x, y, z = map(Symbol, 'xyz')

    assert separate((x*y*z)**4) == x**4*y**4*z**4
    assert separate((x*y*z)**x) == x**x*y**x*z**x
    assert separate((x*(y*z)**2)**3) == x**3*y**6*z**6

    assert separate((sin((x*y)**2)*y)**z) == sin((x*y)**2)**z*y**z
    assert separate((sin((x*y)**2)*y)**z, deep=True) == sin(x**2*y**2)**z*y**z

    assert separate(exp(x)**2) == exp(2*x)
    assert separate((exp(x)*exp(y))**2) == exp(2*x)*exp(2*y)

    assert separate((exp((x*y)**z)*exp(y))**2) == exp(2*(x*y)**z)*exp(2*y)
    assert separate((exp((x*y)**z)*exp(y))**2, deep=True) == exp(2*x**z*y**z)*exp(2*y)

@XFAIL
def test_separate_X1():
    assert separate((exp(x)*exp(y))**z) == exp(x*z)*exp(y*z)

def test_powsimp():
    x,y,n = symbols('xyn')
    f = Function('f')
    assert powsimp( 4**x * 2**(-x) * 2**(-x) ) == 1
    assert powsimp( (-4)**x * (-2)**(-x) * 2**(-x) ) == 1

    assert powsimp( f(4**x * 2**(-x) * 2**(-x)) )   == f(4**x * 2**(-x) * 2**(-x))
    assert powsimp( f(4**x * 2**(-x) * 2**(-x)), deep = True )  == f(1)
    x,y = symbols('xy', nonnegative=True)
    n = Symbol('n', real=True)
    assert powsimp( y**n * (y/x)**(-n) ) == x**n

def test_collect():
    x,y,n = symbols('xyn')
    assert collect(2*x**2 + y*x**2 + 3*x*y, [x]) == x**2*(2+y) + 3*x*y
    assert collect(2*x**2 + y*x**2 + 3*x*y, [y]) == 2*x**2 + y*(x**2+3*x)
    assert collect((1+ (x+y) + (x+y)**2).expand(), [x,y]) == 1 + y + x*(1 + 2*y) + x**2  + y**2

def test_collect_D():
    D = Derivative
    f = Function('f')
    x,a,b = symbols('xab')
    fx  = D(f(x), x)
    fxx = D(f(x), x,x)

    assert collect(a*fx + b*fx, fx) == (a + b)*fx
    assert collect(a*D(fx,x) + b*D(fx,x), fx)   == (a + b)*D(fx, x)
    assert collect(a*fxx     + b*fxx    , fx)   == (a + b)*D(fx, x)

def test_collect_D_0():
    D = Derivative
    f = Function('f')
    x,a,b = symbols('xab')
    fxx = D(f(x), x,x)

    # collect does not distinguish nested derivatives, so it returns
    #                                           -- (a + b)*D(D(f,x), x)
    assert collect(a*fxx     + b*fxx    , fxx)  == (a + b)*fxx

def test_hypersimp():
    n, k = symbols('nk', integer=True)

    assert hypersimp(factorial(k), k) == k + 1
    assert hypersimp(factorial(k**2), k) is None

    assert hypersimp(1/factorial(k), k) == 1/(k + 1)

    assert hypersimp(2**k/factorial(k)**2, k) == 2/(k**2+2*k+1)

    assert hypersimp(binomial(n, k), k) == (n-k)/(k+1)
    assert hypersimp(binomial(n+1, k), k) == (n-k+1)/(k+1)

    term = (4*k+1)*factorial(k)/factorial(2*k+1)
    assert hypersimp(term, k) == (4*k + 5)/(6 + 16*k**2 + 28*k)

    term = 1/((2*k-1)*factorial(2*k+1))
    assert hypersimp(term, k) == (2*k-1)/(6 + 22*k + 24*k**2 + 8*k**3)

    term = binomial(n, k)*(-1)**k/factorial(k)
    assert hypersimp(term, k) == (k - n)/(k**2+2*k+1)

def test_together2():
    x, y, z = symbols("xyz")
    assert together(1/(x*y) + 1/y**2) == 1/x*y**(-2)*(x + y)
    assert together(1/(1 + 1/x)) == x/(1 + x)
    x = symbols("x", nonnegative=True)
    y = symbols("y", real=True)
    assert together(1/x**y + 1/x**(y-1)) == x**(-y)*(1 + x)

def test_nsimplify():
    x = Symbol("x")
    assert nsimplify(0) == 0
    assert nsimplify(-1) == -1
    assert nsimplify(1) == 1
    assert nsimplify(1+x) == 1+x
    assert nsimplify(2.7) == Rational(27,10)
    assert nsimplify(1-GoldenRatio) == (1-sqrt(5))/2
    assert nsimplify((1+sqrt(5))/4, [GoldenRatio]) == GoldenRatio/2
    assert nsimplify(2/GoldenRatio, [GoldenRatio]) == 2*GoldenRatio - 2
    assert nsimplify(exp(5*pi*I/3, evaluate=False)) == sympify('1/2 - I*3**(1/2)/2')
    assert nsimplify(sin(3*pi/5, evaluate=False)) == sympify('(5/8 + 1/8*5**(1/2))**(1/2)')
    assert nsimplify(sqrt(atan('1', evaluate=False))*(2+I), [pi]) == sqrt(pi) + sqrt(pi)/2*I
    assert nsimplify(2 + exp(2*atan('1/4')*I)) == sympify('49/17 + 8*I/17')
    assert nsimplify(pi, tolerance=0.01) == Rational(22,7)
    assert nsimplify(pi, tolerance=0.001) == Rational(355,113)
    assert nsimplify(0.33333, tolerance=1e-4) == Rational(1,3)
    assert nsimplify(2.0**(1/3.), tolerance=0.001) == Rational(635,504)
    assert nsimplify(2.0**(1/3.), tolerance=0.001, full=True) == 2**Rational(1,3)
