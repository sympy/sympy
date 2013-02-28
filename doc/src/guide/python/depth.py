from sympy.core import Atom

def depth(expr):
    if isinstance(expr, Atom):
        return 1
    else:
        return 1 + max([ depth(arg) for arg in expr.args ])
