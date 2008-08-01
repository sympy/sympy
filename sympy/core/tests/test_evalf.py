from sympy.core.evalf import *

from sympy import pi, I, Symbol, Add, Rational, exp, sqrt, sin, cos, fibonacci, \
    Integral, oo, E, atan, log, integrate, floor, ceiling, factorial

import py

x = Symbol('x')
y = Symbol('y')

def NS(e, n=15, **options):
    return str(sympify(e).evalf(n, **options))

def test_helpers():
    assert complex_accuracy((from_float(2.0),None,35,None)) == 35
    assert complex_accuracy((from_float(2.0),from_float(10.0),35,100)) == 37
    assert complex_accuracy((from_float(2.0),from_float(1000.0),35,100)) == 43
    assert complex_accuracy((from_float(2.0),from_float(10.0),100,35)) == 35
    assert complex_accuracy((from_float(2.0),from_float(1000.0),100,35)) == 35

def test_basic():
    assert NS('pi',15) == '3.14159265358979'
    assert NS('2/3',10) == '0.6666666667'
    assert NS('355/113-pi',6) == '2.66764e-7'
    assert NS('16*atan(1/5)-4*atan(1/239)', 15) == '3.14159265358979'

def test_cancellation():
    assert NS(Add(pi,Rational(1,10**1000),-pi,evaluate=False),15,maxprec=1200)=='1.0e-1000'

def test_powers():
    assert NS('pi**(10**20)',10) == '1.339148777e+49714987269413385435'
    assert NS(pi**(10**100),10) == ('4.946362032e+4971498726941338543512682882'
          '9089887365167832438044244613405349992494711208'
          '95526746555473864642912223')
    assert NS('2**(1/10**50)',15) == '1.0'
    assert NS('2**(1/10**50)-1',15) == '6.93147180559945e-51'

# Evaluation of Rump's ill-conditioned polynomial
def test_rump():
    a = 1335*y**6/4+x**2*(11*x**2*y**2-y**6-121*y**4-2)+11*y**8/2+x/(2*y)
    assert NS(a, 15, subs={x:77617, y:33096}) == '-0.827396059946821'

def test_complex():
    assert NS('2*sqrt(pi)*I',10) == '3.544907702*I'
    assert NS('3+3*I',15) == '3.0 + 3.0*I'
    assert NS('E+pi*I',15) == '2.71828182845905 + 3.14159265358979*I'
    assert NS('pi * (3+4*I)',15) == '9.42477796076938 + 12.5663706143592*I'
    assert NS('I*(2+I)',15) == '-1.0 + 2.0*I'
    #assert NS('(pi+E*I)*(E+pi*I)',15) in ('.0e-15 + 17.25866050002*I', '.0e-17 + 17.25866050002*I', '-.0e-17 + 17.25866050002*I')
    assert NS('(pi+E*I)*(E+pi*I)',15,chop=True) == '17.25866050002*I'

def test_complex_powers():
    assert NS('(E+pi*I)**100000000000000000') == \
        '-3.58896782867793e+61850354284995199 + 4.58581754997159e+61850354284995199*I'
    # XXX: rewrite if a+a*I simplification introduced in sympy
    #assert NS('(pi + pi*I)**2') in ('.0e-15 + 19.7392088021787*I', '.0e-16 + 19.7392088021787*I')
    assert NS('(pi + pi*I)**2', chop=True) == '19.7392088021787*I'
    assert NS('(pi + 1/10**8 + pi*I)**2') == '6.2831853e-8 + 19.7392088650106*I'
    assert NS('(pi + 1/10**12 + pi*I)**2') == '6.283e-12 + 19.739208802185*I'
    #assert NS('(pi + pi*I)**4') == '-389.63636413601 + .0e-14*I'
    assert NS('(pi + pi*I)**4', chop=True) == '-389.63636413601'
    assert NS('(pi + 1/10**8 + pi*I)**4') == '-389.636366616512 + 2.4805021e-6*I'
    assert NS('(pi + 1/10**12 + pi*I)**4') == '-389.636364136258 + 2.481e-10*I'
    assert NS('(10000*pi + 10000*pi*I)**4', chop=True) == '-3.8963636413601e+18'

# An example from Smith, "Multiple Precision Complex Arithmetic and Functions"
def test_complex_cancellation():
    A = Rational('63287/100000')
    B = Rational('52498/100000')
    C = Rational('69301/100000')
    D = Rational('83542/100000')
    F = Rational('2231321613/2500000000')
    # XXX: the number of returned mantissa digits in the real part could
    # change with the implementation. What matters is that the returned digits are
    # correct.
    assert NS((A+B*I)*(C+D*I),6) in ('6.45e-6 + 0.892529*I', '6.4e-6 + 0.892529*I')
    assert NS((A+B*I)*(C+D*I),10) == '6.4471e-6 + 0.8925286452*I'
    assert NS((A+B*I)*(C+D*I) - F*I, 5) in ('6.4471e-6 - .0e-13*I', '6.4471e-6 + 2.0e-15*I')

def test_logs():
    assert NS("log(3+pi*I)", 15) == '1.46877619736226 + 0.808448792630022*I'

def test_trig():
    assert NS('sin(1)',15) == '0.841470984807897'
    assert NS('cos(1)',15) == '0.54030230586814'
    assert NS('sin(10**-6)',15) == '9.99999999999833e-7'
    assert NS('cos(10**-6)',15) == '0.9999999999995'
    assert NS('sin(E*10**100)',15) == '0.409160531722613'
    # Some input near roots
    assert NS(sin(exp(pi*sqrt(163))*pi), 15) == '-2.35596641936785e-12'
    assert NS(sin(pi*10**100 + Rational(7,10**5), evaluate=False), 15, maxprec=120) == \
        '6.99999999428333e-5'
    assert NS(sin(Rational(7,10**5), evaluate=False), 15) == \
        '6.99999999428333e-5'

# Check detection of various false identities
def test_near_integers():
    # Binet's formula
    f = lambda n: ((1+sqrt(5))**n)/(2**n * sqrt(5))
    assert NS(f(5000) - fibonacci(5000), 10, maxprec=1500) == '5.156009964e-1046'
    # Some near-integer identities from
    # http://mathworld.wolfram.com/AlmostInteger.html
    assert NS('sin(2017*2**(1/5))',15) == '-1.0'
    assert NS('sin(2017*2**(1/5))',20) == '-0.99999999999999997857'
    assert NS('1+sin(2017*2**(1/5))',15) == '2.1432228738939e-17'
    assert NS('45 - 613*E/37 + 35/991', 15) == '6.03764498766326e-11'

def test_ramanujan():
    assert NS(exp(pi*sqrt(163)) - 640320**3 - 744, 10) == '-7.499274028e-13'
    # A related identity
    A = 262537412640768744*exp(-pi*sqrt(163))
    B = 196884*exp(-2*pi*sqrt(163))
    C = 103378831900730205293632*exp(-3*pi*sqrt(163))
    assert NS(1-A-B+C,10) == '1.613679005e-59'

def test_z_integrals():
    assert NS(Integral(x, (x, 2, 5)), 15) == '10.5'
    gauss = Integral(exp(-x**2), (x, -oo, oo))
    assert NS(gauss, 15) == '1.77245385090552'
    assert NS(gauss**2 - pi + E*Rational(1,10**20), 15) in ('2.71828182845904e-20', '2.71828182845905e-20')
    # A monster of an integral from http://mathworld.wolfram.com/DefiniteIntegral.html
    t = Symbol('t')
    a = 8*sqrt(3)/(1+3*t**2)
    b = 16*sqrt(2)*(3*t+1)*(4*t**2+t+1)**Rational(3,2)
    c = (3*t**2+1)*(11*t**2+2*t+3)**2
    d = sqrt(2)*(249*t**2+54*t+65)/(11*t**2+2*t+3)**2
    f = a - b/c - d
    assert NS(Integral(f, (t, 0, 1)), 50) == NS((3*sqrt(2)-49*pi+162*atan(sqrt(2)))/12,50)
    # http://mathworld.wolfram.com/VardisIntegral.html
    assert NS(Integral(log(log(1/x))/(1+x+x**2), (x, 0, 1)), 15) == NS('pi/sqrt(3) * log(2*pi**(5/6) / gamma(1/6))', 15)
    # http://mathworld.wolfram.com/AhmedsIntegral.html
    assert NS(Integral(atan(sqrt(x**2+2))/(sqrt(x**2+2)*(x**2+1)), (x, 0, 1)), 15) == NS(5*pi**2/96, 15)
    # http://mathworld.wolfram.com/AbelsIntegral.html
    assert NS(Integral(x/((exp(pi*x)-exp(-pi*x))*(x**2+1)), (x, 0, oo)), 15) == NS('log(2)/2-1/4',15)
    # Complex part trimming
    # http://mathworld.wolfram.com/VardisIntegral.html
    assert NS(Integral(log(log(sin(x)/cos(x))), (x, pi/4, pi/2)), 15, chop=True) == \
        NS('pi/4*log(4*pi**3/gamma(1/4)**4)', 15)
    #
    # Endpoints causing trouble (rounding error in integration points -> complex log)
    assert NS(2+Integral(log(2*cos(x/2)), (x, -pi, pi)), 17, chop=True) == '2.0'
    assert NS(2+Integral(log(2*cos(x/2)), (x, -pi, pi)), 20, chop=True) == '2.0'
    assert NS(2+Integral(log(2*cos(x/2)), (x, -pi, pi)), 22, chop=True) == '2.0'
    # Needs zero handling
    assert NS(pi - 4*Integral('sqrt(1-x**2)', (x, 0, 1)), 15, maxprec=30, chop=True) in ('0.0', '0')

def test_z_issue_939():
    # http://code.google.com/p/sympy/issues/detail?id=939
    assert NS(integrate(1/(x**5+1), x).subs(x, 4), chop=True) == '-0.000976138910649103'
    assert NS(Integral(1/(x**5+1), (x, 2, 4))) == '0.014436108888674'
    assert NS(integrate(1/(x**5+1), (x, 2, 4)), chop=True) == '0.014436108888674'

def xtest_failing_integrals():
    #---
    # Double integrals not implemented
    assert NS(Integral(sqrt(x)+x*y, (x, 1, 2), (y, -1, 1)), 15) == '2.43790283299492'
    # double integral + zero detection
    assert NS(Integral(sin(x+x*y), (x, -1, 1), (y, -1, 1)), 15) == '0.0'

# Input that for various reasons have failed at some point
def test_bugs():
    assert NS(sin(1)+exp(-10**10),10) == NS(sin(1),10)
    assert NS(exp(10**10)+sin(1),10) == NS(exp(10**10),10)
    assert NS('log(1+1/10**50)',20) == '1.0e-50'
    assert NS('log(10**100,10)',10) == '100.0'
    assert NS('log(2)',10) == '0.6931471806'
    assert NS('(sin(x)-x)/x**3', 15, subs={x:'1/10**50'}) == '-0.166666666666667'
    assert NS(sin(1)+Rational(1,10**100)*I,15) == '0.841470984807897 + 1.0e-100*I'
    assert x.evalf() == x

def test_integer_parts():
    a = floor(log(8)/log(2) - exp(-1000), evaluate=False)
    b = floor(log(8)/log(2), evaluate=False)
    py.test.raises(PrecisionExhausted, "a.evalf()")
    assert a.evalf(chop=True) == 3
    assert a.evalf(maxprec=500) == 2
    py.test.raises(PrecisionExhausted, "b.evalf()")
    py.test.raises(PrecisionExhausted, "b.evalf(maxprec=500)")
    assert b.evalf(chop=True) == 3
    assert int(floor(factorial(50)/E,evaluate=False).evalf()) == \
        11188719610782480504630258070757734324011354208865721592720336800L
    assert int(ceiling(factorial(50)/E,evaluate=False).evalf()) == \
        11188719610782480504630258070757734324011354208865721592720336801L

def test_trig_zero_detection():
    a = sin(160*pi, evaluate=False)
    t = a.evalf(maxprec=100)
    assert abs(t) < 1e-100
    assert t._prec < 2
    assert a.evalf(chop=True) == 0
    assert py.test.raises(PrecisionExhausted, "a.evalf(strict=True)")
