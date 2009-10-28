from sympy.core.evalf import PrecisionExhausted, complex_accuracy, from_float

from sympy import pi, I, Symbol, Add, Rational, exp, sqrt, sin, cos, \
    fibonacci, Integral, oo, E, atan, log, integrate, floor, ceiling, \
    factorial, binomial, Sum, zeta, Catalan, Pow, GoldenRatio, sympify, sstr

from sympy.utilities.pytest import raises

x = Symbol('x')
y = Symbol('y')
n = Symbol('n')

def NS(e, n=15, **options):
    return sstr(sympify(e).evalf(n, **options), full_prec=True)

def test_evalf_helpers():
    assert complex_accuracy((from_float(2.0),None,35,None)) == 35
    assert complex_accuracy((from_float(2.0),from_float(10.0),35,100)) == 37
    assert complex_accuracy((from_float(2.0),from_float(1000.0),35,100)) == 43
    assert complex_accuracy((from_float(2.0),from_float(10.0),100,35)) == 35
    assert complex_accuracy((from_float(2.0),from_float(1000.0),100,35)) == 35

def test_evalf_basic():
    assert NS('pi',15) == '3.14159265358979'
    assert NS('2/3',10) == '0.6666666667'
    assert NS('355/113-pi',6) == '2.66764e-7'
    assert NS('16*atan(1/5)-4*atan(1/239)', 15) == '3.14159265358979'

def test_cancellation():
    assert NS(Add(pi,Rational(1,10**1000),-pi,evaluate=False),15,maxprec=1200) == '1.00000000000000e-1000'

def test_evalf_powers():
    assert NS('pi**(10**20)',10) == '1.339148777e+49714987269413385435'
    assert NS(pi**(10**100),10) == ('4.946362032e+4971498726941338543512682882'
          '9089887365167832438044244613405349992494711208'
          '95526746555473864642912223')
    assert NS('2**(1/10**50)',15) == '1.00000000000000'
    assert NS('2**(1/10**50)-1',15) == '6.93147180559945e-51'

# Evaluation of Rump's ill-conditioned polynomial
def test_evalf_rump():
    a = 1335*y**6/4+x**2*(11*x**2*y**2-y**6-121*y**4-2)+11*y**8/2+x/(2*y)
    assert NS(a, 15, subs={x:77617, y:33096}) == '-0.827396059946821'

def test_evalf_complex():
    assert NS('2*sqrt(pi)*I',10) == '3.544907702*I'
    assert NS('3+3*I',15) == '3.00000000000000 + 3.00000000000000*I'
    assert NS('E+pi*I',15) == '2.71828182845905 + 3.14159265358979*I'
    assert NS('pi * (3+4*I)',15) == '9.42477796076938 + 12.5663706143592*I'
    assert NS('I*(2+I)',15) == '-1.00000000000000 + 2.00000000000000*I'
    #assert NS('(pi+E*I)*(E+pi*I)',15) in ('.0e-15 + 17.25866050002*I', '.0e-17 + 17.25866050002*I', '-.0e-17 + 17.25866050002*I')
    assert NS('(pi+E*I)*(E+pi*I)',15,chop=True) == '17.2586605000200*I'

def test_evalf_complex_powers():
    assert NS('(E+pi*I)**100000000000000000') == \
        '-3.58896782867793e+61850354284995199 + 4.58581754997159e+61850354284995199*I'
    # XXX: rewrite if a+a*I simplification introduced in sympy
    #assert NS('(pi + pi*I)**2') in ('.0e-15 + 19.7392088021787*I', '.0e-16 + 19.7392088021787*I')
    assert NS('(pi + pi*I)**2', chop=True) == '19.7392088021787*I'
    assert NS('(pi + 1/10**8 + pi*I)**2') == '6.2831853e-8 + 19.7392088650106*I'
    assert NS('(pi + 1/10**12 + pi*I)**2') == '6.283e-12 + 19.7392088021850*I'
    #assert NS('(pi + pi*I)**4') == '-389.63636413601 + .0e-14*I'
    assert NS('(pi + pi*I)**4', chop=True) == '-389.636364136010'
    assert NS('(pi + 1/10**8 + pi*I)**4') == '-389.636366616512 + 2.4805021e-6*I'
    assert NS('(pi + 1/10**12 + pi*I)**4') == '-389.636364136258 + 2.481e-10*I'
    assert NS('(10000*pi + 10000*pi*I)**4', chop=True) == '-3.89636364136010e+18'

def test_evalf_exponentiation():
    assert NS(sqrt(-pi)) == '1.77245385090552*I'
    assert NS(Pow(pi*I, Rational(1,2), evaluate=False)) == '1.25331413731550 + 1.25331413731550*I'
    assert NS(pi**I) == '0.413292116101594 + 0.910598499212615*I'
    assert NS(pi**(E+I/3)) == '20.8438653991931 + 8.36343473930031*I'
    assert NS((pi+I/3)**(E+I/3)) == '17.2442906093590 + 13.6839376767037*I'
    assert NS(exp(pi)) == '23.1406926327793'
    assert NS(exp(pi+E*I)) == '-21.0981542849657 + 9.50576358282422*I'
    assert NS(pi**pi) == '36.4621596072079'
    assert NS((-pi)**pi) == '-32.9138577418939 - 15.6897116534332*I'
    assert NS((-pi)**(-pi)) == '-0.0247567717232697 + 0.0118013091280262*I'

# An example from Smith, "Multiple Precision Complex Arithmetic and Functions"
def test_evalf_complex_cancellation():
    A = Rational('63287/100000')
    B = Rational('52498/100000')
    C = Rational('69301/100000')
    D = Rational('83542/100000')
    F = Rational('2231321613/2500000000')
    # XXX: the number of returned mantissa digits in the real part could
    # change with the implementation. What matters is that the returned digits are
    # correct.
    assert NS((A+B*I)*(C+D*I),6) == '6.44862e-6 + 0.892529*I'
    assert NS((A+B*I)*(C+D*I),10) == '6.447099821e-6 + 0.8925286452*I'
    assert NS((A+B*I)*(C+D*I) - F*I, 5) in ('6.4471e-6 - .0e-15*I', '6.4471e-6 + .0e-15*I')

def test_evalf_logs():
    assert NS("log(3+pi*I)", 15) == '1.46877619736226 + 0.808448792630022*I'
    assert NS("log(pi*I)", 15) == '1.14472988584940 + 1.57079632679490*I'

def test_evalf_trig():
    assert NS('sin(1)',15) == '0.841470984807897'
    assert NS('cos(1)',15) == '0.540302305868140'
    assert NS('sin(10**-6)',15) == '9.99999999999833e-7'
    assert NS('cos(10**-6)',15) == '0.999999999999500'
    assert NS('sin(E*10**100)',15) == '0.409160531722613'
    # Some input near roots
    assert NS(sin(exp(pi*sqrt(163))*pi), 15) == '-2.35596641936785e-12'
    assert NS(sin(pi*10**100 + Rational(7,10**5), evaluate=False), 15, maxprec=120) == \
        '6.99999999428333e-5'
    assert NS(sin(Rational(7,10**5), evaluate=False), 15) == \
        '6.99999999428333e-5'

# Check detection of various false identities
def test_evalf_near_integers():
    # Binet's formula
    f = lambda n: ((1+sqrt(5))**n)/(2**n * sqrt(5))
    assert NS(f(5000) - fibonacci(5000), 10, maxprec=1500) == '5.156009964e-1046'
    # Some near-integer identities from
    # http://mathworld.wolfram.com/AlmostInteger.html
    assert NS('sin(2017*2**(1/5))',15) == '-1.00000000000000'
    assert NS('sin(2017*2**(1/5))',20) == '-0.99999999999999997857'
    assert NS('1+sin(2017*2**(1/5))',15) == '2.14322287389390e-17'
    assert NS('45 - 613*E/37 + 35/991', 15) == '6.03764498766326e-11'

def test_evalf_ramanujan():
    assert NS(exp(pi*sqrt(163)) - 640320**3 - 744, 10) == '-7.499274028e-13'
    # A related identity
    A = 262537412640768744*exp(-pi*sqrt(163))
    B = 196884*exp(-2*pi*sqrt(163))
    C = 103378831900730205293632*exp(-3*pi*sqrt(163))
    assert NS(1-A-B+C,10) == '1.613679005e-59'

# Input that for various reasons have failed at some point
def test_evalf_bugs():
    assert NS(sin(1)+exp(-10**10),10) == NS(sin(1),10)
    assert NS(exp(10**10)+sin(1),10) == NS(exp(10**10),10)
    assert NS('log(1+1/10**50)',20) == '1.0000000000000000000e-50'
    assert NS('log(10**100,10)',10) == '100.0000000'
    assert NS('log(2)',10) == '0.6931471806'
    assert NS('(sin(x)-x)/x**3', 15, subs={x:'1/10**50'}) == '-0.166666666666667'
    assert NS(sin(1)+Rational(1,10**100)*I,15) == '0.841470984807897 + 1.00000000000000e-100*I'
    assert x.evalf() == x
    assert NS((1+I)**2*I,6) == '-2.00000 + 2.32831e-10*I'
    d={n: (-1)**Rational(6,7), y: (-1)**Rational(4,7), x: (-1)**Rational(2,7)}
    assert NS((x*(1+y*(1 + n))).subs(d).evalf(),6) == '0.346011 + 0.433884*I'
    assert NS(((-I-sqrt(2)*I)**2).evalf()) == '-5.82842712474619'
    assert NS((1+I)**2*I,15) == '-2.00000000000000 + 2.16840434497101e-19*I'
    #1659 (1/2):
    assert NS(pi.evalf(69) - pi) == '-4.43863937855894e-71'
    #1659 (2/2): With the bug present, this still only fails if the
    # terms are in the order given here. This is not generally the case,
    # because the order depends on the hashes of the terms.
    assert NS(20 - 5008329267844*n**25 - 477638700*n**37 - 19*n,
              subs={n:.01}) == '19.8100000000000'

def test_evalf_integer_parts():
    a = floor(log(8)/log(2) - exp(-1000), evaluate=False)
    b = floor(log(8)/log(2), evaluate=False)
    raises(PrecisionExhausted, "a.evalf()")
    assert a.evalf(chop=True) == 3
    assert a.evalf(maxprec=500) == 2
    raises(PrecisionExhausted, "b.evalf()")
    raises(PrecisionExhausted, "b.evalf(maxprec=500)")
    assert b.evalf(chop=True) == 3
    assert int(floor(factorial(50)/E,evaluate=False).evalf()) == \
        11188719610782480504630258070757734324011354208865721592720336800L
    assert int(ceiling(factorial(50)/E,evaluate=False).evalf()) == \
        11188719610782480504630258070757734324011354208865721592720336801L
    assert int(floor((GoldenRatio**999 / sqrt(5) + Rational(1,2))).evalf(1000)) == fibonacci(999)
    assert int(floor((GoldenRatio**1000 / sqrt(5) + Rational(1,2))).evalf(1000)) == fibonacci(1000)

def test_evalf_trig_zero_detection():
    a = sin(160*pi, evaluate=False)
    t = a.evalf(maxprec=100)
    assert abs(t) < 1e-100
    assert t._prec < 2
    assert a.evalf(chop=True) == 0
    raises(PrecisionExhausted, "a.evalf(strict=True)")

def test_evalf_divergent_series():
    n = Symbol('n', integer=True)
    raises(ValueError, 'Sum(1/n, (n, 1, oo)).evalf()')
    raises(ValueError, 'Sum(n/(n**2+1), (n, 1, oo)).evalf()')
    raises(ValueError, 'Sum((-1)**n, (n, 1, oo)).evalf()')
    raises(ValueError, 'Sum((-1)**n, (n, 1, oo)).evalf()')
    raises(ValueError, 'Sum(n**2, (n, 1, oo)).evalf()')
    raises(ValueError, 'Sum(2**n, (n, 1, oo)).evalf()')
    raises(ValueError, 'Sum((-2)**n, (n, 1, oo)).evalf()')

def test_evalf_py_methods():
    assert abs(float(pi+1) - 4.1415926535897932) < 1e-10
    assert abs(complex(pi+1) - 4.1415926535897932) < 1e-10
    assert abs(complex(pi+E*I) - (3.1415926535897931+2.7182818284590451j)) < 1e-10
    raises(ValueError, "float(pi+x)")
    raises(ValueError, "complex(pi+x)")

def test_evalf_power_subs_bugs():
    assert (x**2).evalf(subs={x:0}) == 0
    assert sqrt(x).evalf(subs={x:0}) == 0
    assert (x**Rational(2,3)).evalf(subs={x:0}) == 0
    assert (x**x).evalf(subs={x:0}) == 1
    assert (3**x).evalf(subs={x:0}) == 1
    assert exp(x).evalf(subs={x:0}) == 1
    assert ((2+I)**x).evalf(subs={x:0}) == 1
    assert (0**x).evalf(subs={x:0}) == 1
