from sympy.core.add import Add
from sympy.core.symbol import Symbol
from sympy.series.order import O

x = Symbol('x')
l = list(x**i for i in range(1000))
l.append(O(x**1001))

def timeit_order_1x():
    _ = Add(*l)
