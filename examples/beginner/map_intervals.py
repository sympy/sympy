#!/usr/bin/env python

from sympy import Interval
from sympy.abc import a, b, x
from sympy import factor

e = Interval.map(Interval(10, 12), Interval(-1, 1), x)
print(e.subs(x, 10))
print(e.subs(x, 11))
print(e.subs(x, 12))

e = Interval.map(Interval(a, b), Interval(-1, 1), x)
print(factor(e))
