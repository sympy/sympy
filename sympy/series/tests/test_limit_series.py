from sympy import *
from sympy.series.limits_series import mrv_compare, mrv2
from sympy.utilities.pytest import XFAIL

x = Symbol('x', real=True, unbounded=True)
z = Basic.Zero()
o = Basic.One()


###
def mrv(expr, var):
    d,md = {},{}
    mrv2(expr, var, d, md)
    return set(md.keys())
###

"""
def test_MrvCompareLeadTerm_simple_1():
    assert mrv_compare(x,x,x) == '='

def test_MrvCompareLeadTerm_simple_2():
    assert mrv_compare(x,ln(x),x) == '>'
    assert mrv_compare(ln(x),x,x) == '<'

def test_MrvTestCase_simple_x():
    expr = x + 1/x
    assert mrv(expr,x) == set([x])
    d,md = {},{}
    r = mrv2(expr,x,d,md)
    assert r == expr
    assert md == {x:x}

def test_MrvTestCase_simple_poly():
    expr = x**2
    assert mrv(expr,x) == set([x])
    d,md = {},{}
    r = mrv2(expr,x,d,md)
    assert r == expr
    assert md == {x:x}

def test_MrvTestCase_simple_log():
    expr = log(x)
    assert mrv(expr,x) == set([x])
    d,md = {},{}
    r = mrv2(expr,x,d,md)
    assert r == expr
    assert md == {x:x}

def test_MrvTestCase_simple_exp():
    expr = exp(x)
    assert mrv(expr,x) == set([exp(x)])
    d,md = {},{}
    r = mrv2(expr,x,d,md)
    assert md == {expr:r}

def test_MrvTestCase_simple_exp2():
    expr = exp(x**2)
    assert mrv(expr,x) == set([exp(x**2)])
    d,md = {},{}
    r = mrv2(expr,x,d,md)
    assert md == {expr:r}

def test_MrvTestCase_page41_1():
    expr = exp(x+1/x)
    assert mrv(expr,x) == set([expr])
    d,md = {},{}
    r = mrv2(expr,x,d,md)
    assert md == {expr:r}

def test_MrvTestCase_page41_2():
    expr = exp(x+exp(-exp(x)))
    assert mrv(expr,x) == set([exp(exp(x))])
    d,md = {},{}
    r = mrv2(expr,x,d,md)
    assert md.keys()[0] == exp(exp(x))

def test_MrvTestCase_page41_3():
    expr = exp(x+exp(-x))
    assert mrv(expr,x) == set([exp(x+exp(-x)), exp(x)])
    d,md = {},{}
    r = mrv2(expr,x,d,md)
    assert set(md.keys()) == set([exp(x+exp(-x)), exp(x)])

def test_MrvTestCase_page41_ex3_13():
    expr = exp(x+exp(-x**2))
    assert mrv(expr,x) == set([exp(x**2)])
    d,md = {},{}
    r = mrv2(expr,x,d,md)
    assert set(md.keys()) == set([exp(x**2)])

def test_MrvTestCase_page41_4():
    expr = exp(1/x+exp(-x))
    assert mrv(expr,x) == set([exp(x)])
    d,md = {},{}
    r = mrv2(expr,x,d,md)
    assert set(md.keys()) == set([exp(x)])

def test_MrvTestCase_page41_5():
    expr = exp(x**2)+x*exp(x)+ln(x)**x/x
    assert mrv(expr,x) == set([exp(x**2)])
    d,md = {},{}
    r = mrv2(expr,x,d,md)
    assert set(md.keys()) == set([exp(x**2)])

def test_MrvTestCase_page41_6():
    expr = exp(x)*(exp(1/x+exp(-x))-exp(1/x))
    assert mrv(expr,x) == set([exp(x)])
    d,md = {},{}
    r = mrv2(expr,x,d,md)
    assert set(md.keys()) == set([exp(x)])

def test_MrvTestCase_page41_7():
    expr = log(x**2+2*exp(exp(3*x**3*log(x))))
    assert mrv(expr,x) == set([exp(exp(3*x**3*log(x)))])

def test_MrvTestCase_page41_8():
    expr = log(x-log(x))/log(x)
    assert mrv(expr,x) == set([x])
    d,md = {},{}
    r = mrv2(expr,x,d,md)
    assert set(md.keys()) == set([x])

def test_MrvTestCase_page43_ex3_15():
    expr = (exp(1/x-exp(-x))-exp(1/x))/exp(-x)
    assert mrv(expr,x) == set([exp(x)])
    d,md = {},{}
    r = mrv2(expr,x,d,md)
    assert set(md.keys()) == set([exp(x)])

def test_MrvTestCase_page44_ex3_17():
    expr = 1/exp(-x+exp(-x))-exp(x)
    assert mrv(expr,x) == set([exp(x), exp(x-exp(-x))])
    d,md = {},{}
    r = mrv2(expr,x,d,md)
    assert set(md.keys()) == set([exp(x), exp(x-exp(-x))])

@XFAIL
def test_MrvTestCase_page47_ex3_21():
    h = exp(-x/(1+exp(-x)))
    expr = exp(h)*exp(-x/(1+h))*exp(exp(-x+h))/h**2-exp(x)+x
    expected = set([1/h,exp(x),exp(x-h),exp(x/(1+h))])
    # XXX Incorrect result
    assert mrv(expr,x).difference(expected) == set()

def test_MrvTestCase_page51_ex3_25():
    expr = (ln(ln(x)+ln(ln(x)))-ln(ln(x)))/ln(ln(x)+ln(ln(ln(x))))*ln(x)
    assert mrv(expr,x) == set([x])
    d,md = {},{}
    r = mrv2(expr,x,d,md)
    assert set(md.keys()) == set([x])

@XFAIL
def test_MrvTestCase_page56_ex3_27():
    # XXX Fails due to infinite recursion
    expr = exp(-x+exp(-x)*exp(-x*ln(x)))
    assert mrv(expr,x) == set([exp(x*log(x))])
    d,md = {},{}
    r = mrv2(expr,x,d,md)
    assert set(md.keys()) == set([exp(x*log(x))])

def test_MrvTestCase_page60_sec3_5_1():
    expr1 = log(log(x*exp(x*exp(x))+1))
    assert mrv(expr1,x) == set([exp(x*exp(x))])
    d,md = {},{}
    r = mrv2(expr1,x,d,md)
    assert set(md.keys()) == set([exp(x*exp(x))])

def test_MrvTestCase_page60_sec3_5_2():
    expr2 = exp(exp(log(log(x)+1/x)))
    c = mrv_compare(expr2,x,x)
    assert c == '=',`c`
    assert mrv(expr2,x) == set([x])
    d,md = {},{}
    r = mrv2(expr2,x,d,md)
    assert set(md.keys()) == set([x])

def test_MrvTestCase_page60_sec3_5():
    expr1 = log(log(x*exp(x*exp(x))+1))
    expr2 = exp(exp(log(log(x)+1/x)))
    expr = expr1 - expr2
    assert mrv(expr,x) == set([exp(x*exp(x))])
    d,md = {},{}
    r = mrv2(expr1,x,d,md)
    assert set(md.keys()) == set([exp(x*exp(x))])

def test_MrvLimitTestCase_simple():
    assert x.limit(x,2) == 2

def test_MrvLimitTestCase_simple_inf():
    assert x.limit(x,oo) == oo

def test_MrvLimitTestCase_page2():
    assert (x**7/exp(x)).limit(x,oo) == 0

@XFAIL
def test_MrvLimitTestCase_page4():
    expr = 1/(x**(log(log(log(log(1/x))))-1))
    assert expr.limit(x,0) == oo

def test_MrvLimitTestCase_page4_1():
    expr = (log(log(log(log(x))))-1)*log(x)
    assert expr.limit(x,oo) == oo

def test_MrvLimitTestCase_page4_2():
    expr = (log(log(log(x)))-1)*x
    assert expr.limit(x,oo) == oo

def test_MrvLimitTestCase_page12_ex2_5():
    expr = sqrt(ln(x+1))-sqrt(ln(x))
    assert expr.limit(x,oo) == 0

def test_MrvLimitTestCase_page12_ex2_6():
    s = Symbol('s')
    expr = ((1+x)**s-1)/x
    assert expr.limit(x,0) == s

def test_MrvLimitTestCase_page13_ex2_7():
    n = Symbol('n', positive=True)
    m = Symbol('m', positive=True)
    expr = (x**(1/n)-1)/(x**(1/m)-1)
    assert expr.limit(x,1) == m/n

def test_MrvLimitTestCase_page14_ex2_9():
    n = Symbol('n')
    expr = x**n/exp(x)
    assert expr.limit(x,oo) == 0

def test_MrvLimitTestCase_page15_ex2_10():
    expr = x/(x-1)-1/ln(x)
    assert expr.limit(x,1) == Rational(1,2)

def test_MrvLimitTestCase_page15_ex2_11():
    expr = x*ln(x)
    assert expr.limit(x,0) == 0

    expr = x**x
    assert expr.limit(x,0) == 1

def test_MrvLimitTestCase_page16_ex2_13():
    expr = sin(x)/x
    assert expr.limit(x,0) == 1

@XFAIL
def test_MrvLimitTestCase_page16_ex2_14():  # enable it after defining phi
    expr = exp(exp(phi(phi(x))))/x
    assert expr.limit(x,oo) == exp(-Rational(1,2))

def test_MrvLimitTestCase_page18_ex2_15():
    expr = exp(exp(exp(exp(x-1/exp(x)))))/exp(exp(exp(exp(x))))
    assert expr.limit(x,oo) == 0

def test_MrvLimitTestCase_page27_ex2_17():
    expr = exp(x)*(sin(1/x+exp(-x))-sin(1/x))
    assert expr.limit(x,oo) == 1

def test_MrvLimitTestCase_page43_ex3_15():
    expr = (exp(1/x-exp(-x))-exp(1/x))/exp(-x)
    assert expr.limit(x,oo) == -1

def test_MrvLimitTestCase_page44_ex3_17():
    expr = 1/exp(-x+exp(-x))-exp(x)
    assert expr.limit(x,oo) == -1

def test_MrvLimitTestCase_page47_ex3_21():
    h = exp(-x/(1+exp(-x)))
    expr = exp(h)*exp(-x/(1+h))*exp(exp(-x+h))/h**2-exp(x)+x
    assert expr.limit(x,oo) == 2

def test_MrvLimitTestCase_page47_ex3_21_1():
    expr = exp((-x) / (1 + exp(-x)))
    assert expr.limit(x,oo) == 0

def test_MrvLimitTestCase_page51_ex3_25():
    expr = (ln(ln(x)+ln(ln(x)))-ln(ln(x)))/ln(ln(x)+ln(ln(ln(x))))*ln(x)
    assert expr.limit(x,oo) == 1

def test_MrvLimitTestCase_bug_0():
    expr = ln(x+x**2)/ln(x) # (ln(x)+ln(1+x))/ln(x)
    assert expr.limit(x,0) == 1

def test_MrvLimitTestCase_page77_ex5_2():
    expr = exp(sin(1/x+exp(-x))-sin(1/x))
    assert expr.limit(x,oo) == 1

@XFAIL
def test_MrvLimitTestCaseWorkInProgress_page7(): # enable it after defining erf
    expr = (erf(x-exp(x-exp(-exp(x))))-erf(x))*exp(exp(x))*exp(x**2)
    assert expr.limit(x,0) == -2/sqrt(pi)

def test_MrvLimitTestCaseWorkInProgress_page14_complicated(): # incorrect result
    r = Symbol('r',positive=True, bounded=True)
    R = sqrt(sqrt(x**4+2*x**2*(r**2+1)+(r**2-1)**2)+x**2+r**2-1)
    expr = x/R
    assert expr.limit(x,0) == sqrt((1-r**2)/2)

def test_MrvLimitTestCaseWorkInProgress_page60_sec3_5():
    expr = ln(ln(x*exp(x*exp(x))+1))-exp(exp(ln(ln(x))+1/x))
    assert expr.limit(x,oo) == 0

def test_MrvLimitTestCaseWorkInProgress_page78_ex5_3():
    expr = exp(1/sin(x))/exp(cos(x) / sin(x))
    assert expr.limit(x,0) == 1

def test_MrvLimitTestCaseWorkInProgress_page86(): # need interval calculus
    expr = exp(-x)/cos(x)
    assert expr.limit(x,oo) == nan

def test_MrvLimitTestCaseWorkInProgress_page107(): # returns incorrect result 0
    w = Symbol('w')
    expr = w/(sqrt(1+w)*sin(x)**2+sqrt(1-w)*cos(x)**2-1)
    assert expr.limit(w,0) == 2/(1-2*cos(x)**2)

def test_MrvLimitTestCaseWorkInProgress_page108(): # returns incorrect result 0
    w = exp(-x)
    expr = w/(sqrt(1+w)*sin(1/x)**2+sqrt(1-w)*cos(1/x)**2-1)
    assert expr.limit(x,oo) == -2

def test_MrvLimitTestCaseWorkInProgress_page111(): # returns incorrect result 0
    expr = ln(1-(ln(exp(x)/x-1)+ln(x))/x)/x
    assert expr.limit(x,oo) == -1

def test_MrvLimitTestCaseWorkInProgress_page112(): # need branch cut support
    expr = ln(-1+x*I)
    assert expr.limit(x,0) == pi*I
    expr = ln(-1-x*I)
    assert expr.limit(x,0) == -pi*I

def test_MrvLimitTestCaseWorkInProgress_page114(): # returns incorrect result -1, need csgn
    expr = 1/(x*(1+(1/x-1)**(1/x-1)))
    assert expr.limit(x,oo) == -1/(I*pi+1)

def test_MrvLimitTestCaseWorkInProgress_page116(): # returns incorrect?? result -oo
    expr = ln(x*(x+1)/ln(exp(x)+exp(ln(x)**2)*exp(x**2))+1/ln(x))
    assert expr.limit(x,oo) == 0


    # page 122-..
def test_MrvLimitTestCaseComparison_8_1():
    expr = exp(x)*(exp(1/x-exp(-x))-exp(1/x))
    assert expr.limit(x,oo) == -1

def test_MrvLimitTestCaseComparison_8_2(): # returns incorrect result -oo
    expr = exp(x)*(exp(1/x+exp(-x)+exp(-x**2))-exp(1/x-exp(-exp(x))))
    assert expr.limit(x,oo) == 1

def test_MrvLimitTestCaseComparison_8_3(): # fails with assertion error
    expr = exp(exp(x-exp(-x))/(1-1/x))-exp(exp(x))
    assert expr.limit(x,oo) == oo

def test_MrvLimitTestCaseComparison_8_4(): # runs forever
    expr = exp(exp(exp(x)/(1-1/x)))-exp(exp(exp(x)/(1-1/x-ln(x)**(-ln(x)))))
    assert expr.limit(x,oo) == -oo

def test_MrvLimitTestCaseComparison_8_5(): # fails with value error
    expr = exp(exp(exp(x+exp(-x))))/exp(exp(exp(x)))
    assert expr.limit(x,oo) == oo

def test_MrvLimitTestCaseComparison_8_6(): # returns incorrect result 1
    expr = exp(exp(exp(x)))/exp(exp(exp(x-exp(-exp(x)))))
    assert expr.limit(x,oo) == oo

def test_MrvLimitTestCaseComparison_8_7(): # fails with infinite recursion
    expr = exp(exp(exp(x)))/exp(exp(exp(x-exp(-exp(exp(x))))))
    assert expr.limit(x,oo) == 1

def test_MrvLimitTestCaseComparison_8_8(): # fails with infinite recursion
    expr = exp(exp(x))/exp(exp(x-exp(-exp(exp(x)))))
    assert expr.limit(x,oo) == 1

def test_MrvLimitTestCaseComparison_8_9(): # returns incorrect result oo
    expr = ln(x)**2 * exp(sqrt(ln(x))*(ln(ln(x)))**2*exp(sqrt(ln(ln(x)))*ln(ln(ln(x)))**3)) / sqrt(x)
    assert expr.limit(x,oo) == 0

def test_MrvLimitTestCaseComparison_8_10(): # returns incorrect result 0
    expr = (x*ln(x)*(ln(x*exp(x)-x**2))**2)/ln(ln(x**2+2*exp(exp(3*x**3*ln(x)))))
    assert expr.limit(x,oo) == Rational(1,3)

def test_MrvLimitTestCaseComparison_8_11(): # returns incorrect result -oo
    expr = (exp(x*exp(-x)/(exp(-x)+exp(-2*x**2/(1+x))))-exp(x))/x
    assert expr.limit(x,oo) == -exp(2)

def test_MrvLimitTestCaseComparison_8_12(): # fails with infinite recursion
    expr = (3**x+5**x)**(1/x)
    assert expr.limit(x,oo) == 5

def test_MrvLimitTestCaseComparison_8_13():
    expr = x/ln(x**(ln(x**(ln(2)/ln(x)))))
    assert expr.limit(x,oo) == oo

def test_MrvLimitTestCaseComparison_8_14(): # fails with infinite recursion
    expr = exp(exp(2*ln(x**5+x)*ln(ln(x))))/exp(exp(10*ln(x)*ln(ln(x))))
    assert expr.limit(x,oo) == oo

def test_MrvLimitTestCaseComparison_8_15():
    expr = 4 * exp(exp(5*x**(-Rational(5,7))/2+21*x**(Rational(6,11))/8+2*x**-8+54*x**(Rational(49,45))/17))**8 \
               / ln(ln(-ln(4*x**(-Rational(5,14))/3))) ** Rational(7,6) / 9
    assert expr.limit(x,oo) == oo

def test_MrvLimitTestCaseComparison_8_16(): # fails with assertion error
    expr = (exp(4*x*exp(-x)/(1/exp(x)+1/exp(2*x**2/(x+1))))-exp(x))/exp(x)**4
    assert expr.limit(x,oo) == 1

def test_MrvLimitTestCaseComparison_8_17(): # fails with infinite recursion
    expr = exp(x*exp(-x)/(exp(-x)+exp(-2*x**2/(x+1))))/exp(x)
    assert expr.limit(x,oo) == 1

def test_MrvLimitTestCaseComparison_8_18():
    expr = exp(exp(-x/(1+exp(-x))))*exp(-x/(1+exp(-x/(1+exp(-x))))) * exp(exp(-x+exp(-x/(1+exp(-x))))) / exp(-x/(1+exp(-x)))**2 - exp(x) + x
    assert expr.limit(x,oo) == 2

def test_MrvLimitTestCaseComparison_8_19():
    expr = (ln(ln(x)+ln(ln(x)))-ln(ln(x)))/ln(ln(x)+ln(ln(ln(x))))*ln(x)
    assert expr.limit(x,oo) == 1

def test_MrvLimitTestCaseComparison_8_20(): # returns incorrect result 1
    expr = exp(ln(ln(x+exp(ln(x)*ln(ln(x)))))/ln(ln(ln(exp(x)+x+ln(x)))))
    assert expr.limit(x,oo) == exp(1)

def test_MrvLimitTestCaseComparison_8_21(): # fails with assertion error
    expr = exp(x)*(sin(1/x+exp(-x))-sin(1/x+exp(-x**2)))
    assert expr.limit(x,oo) == 1

def test_MrvLimitTestCaseComparison_8_22():
    expr = exp(exp(x))*(exp(sin(1/x+exp(-exp(x))))-exp(sin(1/x)))
    assert expr.limit(x,oo) == 1

@XFAIL
def test_MrvLimitTestCaseComparison_8_23(): # need erf
    expr = (erf(x-exp(-exp(x)))-erf(x))*exp(exp(x))*exp(x**2)
    assert expr.limit(x,oo) == -2/sqrt(pi)

    # ...

@XFAIL
def test_MrvLimitTestCaseComparison_8_37(): # need max_/min_ implementations to handle unbounded args
    expr = max_(x, exp(x))/ln(min_(exp(-x),exp(-exp(x))))
    assert expr.limit(x,oo) == 1

def test_interface():
    assert limit(sin(x)/x, x, 0) == 1

def test_loglogloglimit():
    #from Gruntz somewhere
    e=(log(log(x)+log(log(x)))-log(log(x)))/log(log(x)+log(log(log(x))))*log(x)
    assert limit(e,x,oo) == 1


def test_logexplimit():
    a = Symbol('a')
    assert limit(x-ln(1+exp(x)), x, oo) == 0
    assert limit(x-ln(a+exp(x)), x, oo) == 0
    assert limit(exp(x)/(1+exp(x)), x, oo) == 1
    assert limit(exp(x)/(a+exp(x)), x, oo) == 1
"""
