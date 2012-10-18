from sympy import Symbol, limit, oo

x = Symbol('x')


def timeit_limit_1x():
    limit(1/x, x, oo)
