from sympy import S, Symbol, sin, cos
from sympy.integrals.rubi.SineIntegrationRules import intsin5, intsin12

# The various functions in SineIntegrationRules return the antiderivatives of
# the following expressions:
#   intsin5(a,b,c,n,x)
#       (c*sin(a+b*x))**n
#   intsin6(a,b,c,d,n,x)
#       (a+b*sin(c+d*x))**n
#   intsin9(a,b,c,d,e,f,m,n,x)
#       (a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**n
#   intsin(a,b,c,d,e,f,A,B,m,n,x)
#       x*(a+b*sin(e+f*x))**m*(A+B*sin(e+f*x))*(c+d*sin(e+f*x))**n
#   intsin12(a,b,c,d,e,f,A,B,C,m,n,x)
#       (a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**n*(A+B*sin(e+f*x)+C*sin(e+f*x)**2)
# where the parameters (a,b,c,d,e,f,A,B,C) and exponents (m,n) are arbitrary
# numeric or symbolic expressions, including zeros and ones.

def test_intsin5():
    x = Symbol("x")
    assert intsin5(S(1), S(2), S(3), S(2), x) == 9*x/2 - 9*sin(2*x + 1)*cos(2*x + 1)/4
    assert intsin5(S(1), S(2), S(3), -S(2), x) == -cos(2*x + 1)/(18*sin(2*x + 1))

def test_intsin12():
    x = Symbol("x")
    assert intsin12(S(0), S(3), S(1), S(0), S(1), S(2), S(1), S(0), S(0), S(2), S(0), x) == 9*x/2 - 9*sin(2*x + 1)*cos(2*x + 1)/4
