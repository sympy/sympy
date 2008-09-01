from sympy import factor, Symbol

x = Symbol('x')


def timeit_factor_x4_x_1():
    factor(x**4+x+1)


def bench_factor_x56_1():
    factor(x**56-1)
