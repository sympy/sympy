#!/usr/bin/env python
import iam_sympy_example

"""
Go to the sympy root directory and execute:

$ python examples/fem_test.py
[  1/60,     0, -1/360,     0, -1/90, -1/360]
[     0,  4/45,      0,  2/45,  2/45,  -1/90]
[-1/360,     0,   1/60, -1/90,     0, -1/360]
[     0,  2/45,  -1/90,  4/45,  2/45,      0]
[ -1/90,  2/45,      0,  2/45,  4/45,      0]
[-1/360, -1/90, -1/360,     0,     0,   1/60]
"""

from fem import *

t = ReferenceSimplex(2)
fe = Lagrange(2,2)

u = 0
#compute u = sum_i u_i N_i
us = []
for i in range(0, fe.nbf()):
   ui = Symbol("u_%d" % i)
   us.append(ui)
   u += ui*fe.N[i]


J = zeronm(fe.nbf(), fe.nbf())
for i in range(0, fe.nbf()):
   Fi = u*fe.N[i]
   for j in range(0, fe.nbf()):
       uj = us[j]
       integrands = diff(Fi, uj)
       J[j,i] = t.integrate(integrands)


pprint(J)
