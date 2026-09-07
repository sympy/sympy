from sympy import symbols, expand

x = symbols('x')

def timeit_expand_small():
    expand((x + 1)**10)

def timeit_expand_medium():
    expand((x + 1)**50)

def timeit_expand_large():
    expand((x + 1)**100)
