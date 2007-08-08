from __future__ import division

from sympy import Basic, Symbol
from math import sin, cos, tan, asin, acos, atan, log, pi
Pi = pi

def lambdify(*args):
    """
    >>> from sympy import symbols, sqrt
    >>> x,y,z = symbols('xyz')
    >>> f = lambdify(x**2, [x])
    >>> f(2)
    4
    >>> f = lambdify([z,y,x], [x,y,z])
    >>> f(1,2,3)
    (3, 2, 1)
    >>> f = lambdify(sqrt(x), [x])
    >>> f(4)
    2.0
    """
    lstr = lambdastr(*args)
    #print lstr
    return eval(lstr)

def lambdastr(*args):
    """
    >>> from sympy import symbols
    >>> x,y,z = symbols('xyz')
    >>> lambdastr(x**2, [x])
    'lambda x: (x**2)'
    >>> lambdastr([z,y,x], [x,y,z])
    'lambda x,y,z: (z,y,x)'
    """
    assert len(args) == 2

    if isinstance(args[0], str) and isinstance(args[1], str):
        exprs, vargs = args

    elif isinstance(args[1], (list, tuple)):
        vargs = args[1]

        for v in vargs:
            assert isinstance(v, Symbol)

        if isinstance(args[0], (list, tuple)):
            exprs = list(args[0])
        else: exprs = [args[0]]

        for e in xrange(len(exprs)):
            exprs[e] = Basic.sympify(exprs[e])
            for a in exprs[e].atoms(type=Symbol):
                assert a in vargs

        vargs = ','.join(str(v) for v in vargs)
        exprs = ','.join(str(e) for e in exprs)

    else: raise ValueError("Lambdification requires arguments "
                           "of the form expr(s), vars. Examples: "
                           "(x**2, [x]) or ([x**2, 1/y], [x,y]).")

    return "lambda %s: (%s)" % (vargs, exprs)

if __name__ == '__main__':
    from sympy import symbols
    x,y,z = symbols('xyz')
    print lambdastr(x**2, [x])
    print lambdastr([z,y,x], [x,y,z])
