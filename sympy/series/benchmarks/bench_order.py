from __future__ import division, print_function

from sympy import Add, O, Symbol

x = Symbol('x')
l = list(x**i for i in range(1000))
l.append(O(x**1001))

def timeit_order_1x():
    _ = Add(*l)
