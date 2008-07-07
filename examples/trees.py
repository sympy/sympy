#!/usr/bin/env python
import iam_sympy_example

from sympy import *

def T(x):
    return x+x**2+2*x**3 + 4*x**4 + 9*x**5 + 20*x**6 + 48 * x**7 + \
            115* x**8 + 286*x**9+719*x**10

def A(x):
    return 1 + T(x) - T(x)**2/2 + T(x**2)/2


x=Symbol("x")
s = Poly(A(x), x)
num = [s.coeff(n) for n in range(11)]

print s.as_basic()
print num
