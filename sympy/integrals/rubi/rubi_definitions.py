def eq(a, b):
    return (a-b) == 0

def gt(a, b):
    return (a-b).is_positive

def lt(a, b):
    return (a-b).is_negative

def ge(a, b):
    return eq(a, b) or gt(a, b)

def le(a, b):
    return eq(a, b) or lt(a, b)

def integer(a):
    return a.is_integer

def subst(a, x, y):
    return a.subs(x, y)

def integrate(a, x):
    return a.integrate(x)
