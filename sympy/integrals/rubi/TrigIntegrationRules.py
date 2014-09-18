from sympy import pi
from .SineIntegrationRules import intsin12

def intcos12(a,b,c,d,e,f,A,B,C,m,n,x):
    # intsin12 integrates:
    # (a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**n*(A+B*sin(e+f*x)+C*sin(e+f*x)**2)
    # To change all sin() into cos(), we substitute:
    #     e -> e + pi/2
    return intsin12(a,b,c,d,e+pi/2,f,A,B,C,m,n,x)
