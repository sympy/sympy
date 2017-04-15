from __future__ import division, print_function

from sympy import *

n = Symbol('n')

expr = cos(1/n) - 1

expr2 = expr.subs(n, 1/n)

t = expr2.series(n, 0, 4)
t1 = expr2.series(n, 0, 6)
t2 = expr2.series(n, 0, 2)


O(expr, (n, S.Infinity))

s = Sum(expr, (n, 1, oo))

s1 = s.is_convergent()

print(Sum(expr, (n, 1, oo)).is_convergent())

expr2 = 1 - cos(1/n)

print(Sum(expr2, (n, 1, oo)).is_convergent())







