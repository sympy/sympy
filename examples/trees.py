import sys
sys.path.append("..")
from sympy import *
from sympy.modules.polynomials import coeff

def T(x):
    return x+x**2+2*x**3 + 4*x**4 + 9*x**5 + 20*x**6 + 48 * x**7 + \
            115* x**8 + 286*x**9+719*x**10

def A(x):
    return 1 + T(x) - T(x)**2/2 + T(x**2)/2


x=Symbol("x")
s = A(x).expand()
num = [coeff(s, x, n) for n in range(11)]

print s
print num
