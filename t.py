#!/usr/bin/env python

from sympy import *

# x = Symbol ('x')
# m = Matrix(8, 8, [x+i for i in range (64)]) [:4, :4]
# m = Matrix(8, 8, ([1+x, 1-x]*4 + [1-x, 1+x]*4)*4) [:4, :4]

# r = m**4

# print (r)


M = Matrix(S('''[
    [             -3/4,       45/32 - 37*I/16,                   0,                     0],
    [-149/64 + 49*I/32, -177/128 - 1369*I/128,                   0, -2063/256 + 541*I/128],
    [                0,         9/4 + 55*I/16, 2473/256 + 137*I/64,                     0],
    [                0,                     0,                   0, -177/128 - 1369*I/128]]'''))
r = M.inv('LDL')

print (r)
