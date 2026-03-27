from sympy.core import symbols, I

x, y, z = symbols('x,y,z')

p = 3*x**2*y*z**7 + 7*x*y*z**2 + 4*x + x*y**4
e = (x + y + z + 1)**32


def timeit_expand_nothing_todo():
    p.expand()


def bench_expand_32():
    """(x+y+z+1)**32  -> expand"""
    e.expand()


def timeit_expand_complex_number_1():
    ((2 + 3*I)**1000).expand(complex=True)


def timeit_expand_complex_number_2():
    ((2 + 3*I/4)**1000).expand(complex=True)

# Additional benchmarks to extend expand() coverage:
# - large polynomial expansion
# - multivariable expressions
# - nested expressions
# - trigonometric cases

def timeit_expand_large_polynomial():
    expr = (x + 1)**100
    expr.expand()


def timeit_expand_multivariable():
    expr = (x + y + z)**20
    expr.expand()


def timeit_expand_nested_expression():
    expr = ((x + 1)**10 + (x + 2)**10)**5
    expr.expand()


def timeit_expand_repeated():
    expr = (x + 1)**20
    expr.expand()


def timeit_expand_trig():
    from sympy import sin
    expr = (sin(x + 1) + sin(x + 2))**5
    expr.expand()
