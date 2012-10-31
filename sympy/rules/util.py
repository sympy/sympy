from sympy import Basic

new = Basic.__new__

def is_leaf(x):
    return not isinstance(x, Basic) or x.is_Atom

def children(x):
    return x.args
