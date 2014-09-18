from sympy import pi, I
from .SineIntegrationRules import intsin12

def intcos12(a,b,c,d,e,f,A,B,C,m,n,x):
    # intsin12 integrates:
    # (a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**n*(A+B*sin(e+f*x)+C*sin(e+f*x)**2)
    # To change all sin() into cos(), we substitute:
    #     e -> e + pi/2
    return intsin12(a,b,c,d,e+pi/2,f,A,B,C,m,n,x)

def intsinh12(a,b,c,d,e,f,A,B,C,m,n,x):
    # intsin12 integrates:
    # (a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**n*(A+B*sin(e+f*x)+C*sin(e+f*x)**2)
    # To change all sin() into sinh(), we substitute:
    #     e -> I*e
    #     f -> I*f
    #     b -> -I*b
    #     d -> -I*d
    #     B -> -I*B
    #     C -> -C
    return intsin12(a,-I*b,c,-I*d,I*e,I*f,A,-I*B,-C,m,n,x)
