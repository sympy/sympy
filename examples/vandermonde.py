from sympy import sqrt, Symbol
from sympy.modules.matrices import *
from sympy.modules.polynomials import *
from sympy.core.power import *
from sympy.core.numbers import *

w, x, y, z = Symbol('w'), Symbol('x'), Symbol('y'), Symbol('z')
L = [x,y,z]
V = eye(len(L))
for i in range(len(L)):
    for j in range(len(L)):
        V[i,j] = L[i]**j
print V.det()
