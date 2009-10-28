from sympy import Symbol, symbols, together, hypersimp, factorial, binomial, \
        collect, Function, powsimp, separate, sin, exp, Rational, fraction, \
        simplify, trigsimp, cos, tan, cot, log, ratsimp, Matrix, pi, integrate, \
        solve, nsimplify, GoldenRatio, sqrt, E, I, sympify, atan, Derivative, \
        S, diff, oo, Eq, Integer, gamma, acos, Integral, logcombine, separatevars
from sympy.utilities import all
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

def test_trigsimp1():
    x, y = symbols('x y')

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

    assert trigsimp(sin(x)/cos(x)) == tan(x)
    assert trigsimp(2*tan(x)*cos(x)) == 2*sin(x)
    assert trigsimp(cot(x)**3*sin(x)**3) == cos(x)**3
    assert trigsimp(y*tan(x)**2/sin(x)**2) == y/cos(x)**2
    assert trigsimp(cot(x)/cos(x)) == 1/sin(x)

    assert trigsimp(cos(0.12345)**2 + sin(0.12345)**2) == 1
    e = 2*sin(x)**2 + 2*cos(x)**2
    assert trigsimp(log(e), deep=True) == log(2)

def test_trigsimp2():
    x, y = symbols('x y')
    assert trigsimp(cos(x)**2*sin(y)**2 + cos(x)**2*cos(y)**2 + sin(x)**2,
            recursive=True) == 1
    assert trigsimp(sin(x)**2*sin(y)**2 + sin(x)**2*cos(y)**2 + cos(x)**2,
            recursive=True) == 1

def test_issue1274():
    x = Symbol("x")
    assert abs(trigsimp(2.0*sin(x)**2+2.0*cos(x)**2)-2.0) < 1e-10

def test_trigsimp3():
    x, y = symbols('x y')
    assert trigsimp(sin(x)/cos(x)) == tan(x)
    assert trigsimp(sin(x)**2/cos(x)**2) == tan(x)**2
    assert trigsimp(sin(x)**3/cos(x)**3) == tan(x)**3
    assert trigsimp(sin(x)**10/cos(x)**10) == tan(x)**10

    assert trigsimp(cos(x)/sin(x)) == 1/tan(x)
    assert trigsimp(cos(x)**2/sin(x)**2) == 1/tan(x)**2
    assert trigsimp(cos(x)**10/sin(x)**10) == 1/tan(x)**10

    assert trigsimp(tan(x)) == trigsimp(sin(x)/cos(x))

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

    assert all(simplify(tmp)==tmp for tmp in [I,E,oo,x,-x,-oo,-E,-I])

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

def test_simplify_issue_1308():
    assert simplify(exp(-Rational(1,2)) + exp(-Rational(3,2))) == \
        (1 + E)*exp(-Rational(3,2))
    assert simplify(exp(1)+exp(-exp(1))) == (1 + exp(1 + E))*exp(-E)


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

def test_separate_X1():
    x, y, z = map(Symbol, 'xyz')
    assert separate((exp(x)*exp(y))**z) == exp(x*z)*exp(y*z)

def test_powsimp():
    x,y,z,n = symbols('xyzn')
    f = Function('f')
    assert powsimp( 4**x * 2**(-x) * 2**(-x) ) == 1
    assert powsimp( (-4)**x * (-2)**(-x) * 2**(-x) ) == 1

    assert powsimp( f(4**x * 2**(-x) * 2**(-x)) )   == f(4**x * 2**(-x) * 2**(-x))
    assert powsimp( f(4**x * 2**(-x) * 2**(-x)), deep = True )  == f(1)
    assert exp(x)*exp(y) == exp(x)*exp(y)
    assert powsimp(exp(x)*exp(y)) == exp(x+y)
    assert powsimp(exp(x)*exp(y)*2**x*2**y) == (2*E)**(x + y)
    assert powsimp(exp(x)*exp(y)*2**x*2**y, combine='exp') == exp(x+y)*2**(x+y)
    assert powsimp(exp(x)*exp(y)*exp(2)*sin(x)+sin(y)+2**x*2**y) == exp(2+x+y)*sin(x)+sin(y)+2**(x+y)
    assert powsimp(sin(exp(x)*exp(y))) == sin(exp(x)*exp(y))
    assert powsimp(sin(exp(x)*exp(y)), deep=True) == sin(exp(x+y))
    assert powsimp(x**2*x**y) == x**(2+y)
    # This should remain factored, because 'exp' with deep=True is supposed
    # to act like old automatic exponent combining.
    assert powsimp((1 + E*exp(E))*exp(-E), combine='exp', deep=True) == (1 + exp(1 + E))*exp(-E)
    assert powsimp((1 + E*exp(E))*exp(-E), deep=True) == exp(1) + exp(-E)
    # This should not change without deep.  Otherwise, simplify() will fail.
    assert powsimp((1 + E*exp(E))*exp(-E)) == (1 + E*exp(E))*exp(-E)
    assert powsimp((1 + E*exp(E))*exp(-E), combine='exp') == (1 + E*exp(E))*exp(-E)
    assert powsimp((1 + E*exp(E))*exp(-E), combine='base') == (1 + E*exp(E))*exp(-E)
    x,y = symbols('xy', nonnegative=True)
    n = Symbol('n', real=True)
    assert powsimp( y**n * (y/x)**(-n) ) == x**n
    assert powsimp(x**(x**(x*y)*y**(x*y))*y**(x**(x*y)*y**(x*y)),deep=True) == (x*y)**(x*y)**(x*y)
    assert powsimp(2**(2**(2*x)*x), deep=False) == 2**(2**(2*x)*x)
    assert powsimp(2**(2**(2*x)*x), deep=True) == 2**(x*4**x)
    assert powsimp(exp(-x + exp(-x)*exp(-x*log(x))), deep=False, combine='exp') == exp(-x + exp(-x)*exp(-x*log(x)))
    assert powsimp(exp(-x + exp(-x)*exp(-x*log(x))), deep=False, combine='exp') == exp(-x + exp(-x)*exp(-x*log(x)))
    assert powsimp((x+y)/(3*z), deep=False, combine='exp') == (x+y)/(3*z)
    assert powsimp((x/3+y/3)/z, deep=True, combine='exp') == (x/3+y/3)/z
    assert powsimp(exp(x)/(1 + exp(x)*exp(y)), deep=True) == exp(x)/(1 + exp(x + y))
    assert powsimp(x*y**(z**x*z**y), deep=True) == x*y**(z**(x + y))
    assert powsimp((z**x*z**y)**x, deep=True) == (z**(x + y))**x
    assert powsimp(x*(z**x*z**y)**x, deep=True) == x*(z**(x + y))**x


def test_collect_1():
    """Collect with respect to a Symbol"""
    x, y, z, n = symbols('xyzn')
    assert collect( x + y*x, x ) == x * (1 + y)
    assert collect( x + x**2, x ) == x + x**2
    assert collect( x**2 + y*x**2, x ) == (x**2)*(1+y)
    assert collect( x**2 + y*x, x ) == x*y + x**2
    assert collect( 2*x**2 + y*x**2 + 3*x*y, [x] ) == x**2*(2+y) + 3*x*y
    assert collect( 2*x**2 + y*x**2 + 3*x*y, [y] ) == 2*x**2 + y*(x**2+3*x)

    assert collect( ((1 + y + x)**4).expand(), x) == ((1 + y)**4).expand() + \
                x*(4*(1 + y)**3).expand() + x**2*(6*(1 + y)**2).expand() + \
                x**3*(4*(1 + y)).expand() + x**4

def test_collect_2():
    """Collect with respect to a sum"""
    a, b, x = symbols('abx')
    assert collect(a*(cos(x)+sin(x)) + b*(cos(x)+sin(x)), sin(x)+cos(x)) == (a + b)*(cos(x) + sin(x))

def test_collect_3():
    """Collect with respect to a product"""
    a, b, c = symbols('abc')
    f = Function('f')
    x,y,z, n = symbols('xyzn')

    assert collect(-x/8 + x*y, -x) == -x*(S.One/8 - y)

    assert collect( 1 + x*(y**2), x*y ) == 1 + x*(y**2)
    assert collect( x*y + a*x*y, x*y) == x*y*(1 + a)
    assert collect( 1 + x*y + a*x*y, x*y) == 1 + x*y*(1 + a)
    assert collect(a*x*f(x) + b*(x*f(x)), x*f(x)) == x*(a + b)*f(x)

    assert collect(a*x*log(x) + b*(x*log(x)), x*log(x)) == x*(a + b)*log(x)
    assert collect(a*x**2*log(x)**2 + b*(x*log(x))**2, x*log(x)) == x**2*log(x)**2*(a + b)

    # with respect to a product of three symbols
    assert collect(y*x*z+a*x*y*z, x*y*z) == (1 + a)*x*y*z

def test_collect_4():
    """Collect with respect to a power"""
    a, b, c, x = symbols('abcx')

    assert collect(a*x**c + b*x**c, x**c) == x**c*(a + b)
    assert collect(a*x**(2*c) + b*x**(2*c), x**c) == (x**2)**c*(a + b)

def test_collect_5():
    """Collect with respect to a tuple"""
    a, x, y, z, n = symbols('axyzn')
    assert collect(x**2*y**4 + z*(x*y**2)**2 + z + a*z, [x*y**2, z]) in [
                z*(1 + a + x**2*y**4) + x**2*y**4,
                z*(1 + a) + x**2*y**4*(1 + z) ]
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
    # 1685
    assert collect(5*f(x)+3*fx, fx) == 5*f(x) + 3*fx
    assert collect(f(x) + f(x)*diff(f(x), x) + x*diff(f(x), x)*f(x), f(x).diff(x)) ==\
    (x*f(x) + f(x))*D(f(x), x) + f(x)
    assert collect(f(x) + f(x)*diff(f(x), x) + x*diff(f(x), x)*f(x), f(x).diff(x), exact=True) ==\
    (x*f(x) + f(x))*D(f(x), x) + f(x)
    assert collect(1/f(x) + 1/f(x)*diff(f(x), x) + x*diff(f(x), x)/f(x), f(x).diff(x), exact=True) ==\
    (1/f(x) + x/f(x))*D(f(x), x) + 1/f(x)

@XFAIL
def collect_issues():
    assert collect(1/f(x) + 1/f(x)*diff(f(x), x) + x*diff(f(x), x)/f(x), f(x).diff(x)) !=\
    (1 + x*D(f(x), x) + D(f(x), x))/f(x)

def test_collect_D_0():
    D = Derivative
    f = Function('f')
    x,a,b = symbols('xab')
    fxx = D(f(x), x,x)

    # collect does not distinguish nested derivatives, so it returns
    #                                           -- (a + b)*D(D(f,x), x)
    assert collect(a*fxx     + b*fxx    , fxx)  == (a + b)*fxx

def test_separatevars():
    x,y,z,n = symbols('xyzn')
    assert separatevars(2*n*x*z+2*x*y*z) == 2*x*z*(n+y)
    assert separatevars(x*z+x*y*z) == x*z*(1+y)
    assert separatevars(pi*x*z+pi*x*y*z) == pi*x*z*(1+y)
    assert separatevars(x*y**2*sin(x) + x*sin(x)*sin(y)) == x*(sin(y) + y**2)*sin(x)
    assert separatevars(x*exp(x+y)+x*exp(x)) == x*(1 + exp(y))*exp(x)
    assert separatevars((x*(y+1))**z) == x**z*(1 + y)**z
    assert separatevars(1+x+y+x*y) == (x+1)*(y+1)

@XFAIL
def test_separatevars_advanced_factor():
    # If factor() is ever improved to factor non-symbolic expressions, this
    # should XPASS
    assert separatevars(1 + log(x)*log(y) + log(x) + log(y)) == (log(x) + 1)*(log(y) + 1)
    assert separatevars(1 + x - log(z) - x*log(z) - exp(y)*log(z) - \
        x*exp(y)*log(z) + x*exp(y) + exp(y)) == \
        (1 + x)*(1 - log(z))*(1 + exp(y))
    x, y = symbols('xy', positive=True)
    assert separatevars(1 + log(x**log(y)) + log(x*y)) == (log(x) + 1)*(log(y) + 1)

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

def test_extract_minus_sign():
    x = Symbol("x")
    y = Symbol("y")
    a = Symbol("a")
    b = Symbol("b")
    assert simplify(-x/-y) == x/y
    assert simplify(-x/y) == -x/y
    assert simplify(x/y) == x/y
    assert simplify(x/-y) == -x/y
    assert simplify(-x/0) == -oo*x
    assert simplify(S(-5)/0) == -oo
    assert simplify(-a*x/(-y-b)) == a*x/(b + y)

def test_diff():
    x = Symbol("x")
    y = Symbol("y")
    f = Function("f")
    g = Function("g")
    assert simplify(g(x).diff(x)*f(x).diff(x)-f(x).diff(x)*g(x).diff(x)) == 0
    assert simplify(2*f(x)*f(x).diff(x)-diff(f(x)**2,x)) == 0
    assert simplify(diff(1/f(x),x)+f(x).diff(x)/f(x)**2) == 0
    assert simplify(f(x).diff(x,y)-f(x).diff(y,x)) == 0

def test_logcombine_1():
    x, y = symbols("xy")
    a = Symbol("a")
    z, w = symbols("zw", positive=True)
    b = Symbol("b", real=True)
    assert logcombine(log(x)+2*log(y)) == log(x) + 2*log(y)
    assert logcombine(log(x)+2*log(y), assume_pos_real=True) == log(x*y**2)
    assert logcombine(a*log(w)+log(z)) == a*log(w) + log(z)
    assert logcombine(b*log(z)+b*log(x)) == log(z**b) + b*log(x)
    assert logcombine(b*log(z)-log(w)) == log(z**b/w)
    assert logcombine(log(x)*log(z)) == log(x)*log(z)
    assert logcombine(log(w)*log(x)) == log(w)*log(x)
    assert logcombine(cos(-2*log(z)+b*log(w))) == cos(log(w**b/z**2))
    assert logcombine(log(log(x)-log(y))-log(z), assume_pos_real=True) == \
        log(log((x/y)**(1/z)))
    assert logcombine((2+I)*log(x), assume_pos_real=True) == I*log(x)+log(x**2)
    assert logcombine((x**2+log(x)-log(y))/(x*y), assume_pos_real=True) == \
        log(x**(1/(x*y))*y**(-1/(x*y)))+x/y
    assert logcombine(log(x)*2*log(y)+log(z), assume_pos_real=True) == \
        log(z*y**log(x**2))
    assert logcombine((x*y+sqrt(x**4+y**4)+log(x)-log(y))/(pi*x**Rational(2,3)*\
        y**Rational(3,2)), assume_pos_real=True) == \
        log(x**(1/(pi*x**Rational(2,3)*y**Rational(3,2)))*y**(-1/(pi*\
        x**Rational(2,3)*y**Rational(3,2)))) + (x**4 + y**4)**Rational(1,2)/(pi*\
        x**Rational(2,3)*y**Rational(3,2)) + x**Rational(1,3)/(pi*y**Rational(1,2))
    assert logcombine(Eq(log(x), -2*log(y)), assume_pos_real=True) == \
        Eq(log(x*y**2), Integer(0))
    assert logcombine(Eq(y, x*acos(-log(x/y))), assume_pos_real=True) == \
        Eq(y, x*acos(log(y/x)))
    assert logcombine(gamma(-log(x/y))*acos(-log(x/y)), assume_pos_real=True) == \
        acos(log(y/x))*gamma(log(y/x))
    assert logcombine((2+3*I)*log(x), assume_pos_real=True) == \
        log(x**2)+3*I*log(x)
    assert logcombine(Eq(y, -log(x)), assume_pos_real=True) == Eq(y, log(1/x))
    assert logcombine(Integral((sin(x**2)+cos(x**3))/x,x), assume_pos_real=True) == \
        Integral((sin(x**2)+cos(x**3))/x,x)
    assert logcombine(Integral((sin(x**2)+cos(x**3))/x,x)+ (2+3*I)*log(x), \
        assume_pos_real=True) == log(x**2)+3*I*log(x) + \
        Integral((sin(x**2)+cos(x**3))/x,x)

@XFAIL
def test_logcombine_2():
    # The same as one of the tests above, but with Rational(a,b) replaced with a/b.
    # This fails because of a bug in matches.  See issue 1274.
    x, y = symbols("xy")
    assert logcombine((x*y+sqrt(x**4+y**4)+log(x)-log(y))/(pi*x**(2/3)*y**(3/2)), \
        assume_pos_real=True) == log(x**(1/(pi*x**(2/3)*y**(3/2)))*y**(-1/\
        (pi*x**(2/3)*y**(3/2)))) + (x**4 + y**4)**(1/2)/(pi*x**(2/3)*y**(3/2)) + \
        x**(1/3)/(pi*y**(1/2))
