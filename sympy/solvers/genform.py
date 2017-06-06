from sympy.solvers.solveset import linsolve
from sympy import S, Dummy, Symbol, pi
from sympy.simplify import (simplify, collect, powsimp, posify, powdenest,
                            nsimplify, denom, logcombine)
from sympy.polys import factor
from sympy.solvers.solvers import solve
from sympy.sets import ConditionSet, imageset
from sympy.core.function import Lambda


def genform(*args):
    """
    To compute General Form of any system.
    Right now able to generate linear general form
    Examples
    ========

    >>> from sympy.solvers.genform import genform
    >>> from sympy.abc import x
    >>> from sympy import sqrt, sin, cos
    >>> from sympy.solvers.solvers import solve
    >>> genform(pi/6,pi/2,pi)
    ImageSet(Lambda(_n, pi*(2*_n - 1)/6), Naturals())
    >>> genform(pi/2 , 3*pi/2)
    ImageSet(Lambda(_n, pi*(2*_n - 1)/2), Naturals())
    >>> genform(solve((sin(x)-1)*(sin(x)+1),x)[1],solve((sin(x)-1)*(sin(x)+1),x)[2])
    ImageSet(Lambda(_n, pi*(2*_n - 1)/2), Naturals())


    TODO implement for quadratic general form.
    """
    if len(args) == 1:
        return args[0]
    # TODO check linear general form or quadratic general form.
    # https://github.com/sympy/sympy/wiki/
    #          Solveset-and-Solver-Discussion#simplified-general-form-for-trigonometric-solution
    else:
        a = Dummy('a', real=True)
        b = Dummy('b', real=True)
        n = Dummy('n', real =True)
        equations = []
        for i in xrange(0,len(args)):
            eq = a*(i+1) + b - args[i]
            equations.append(eq)

        a, b= list(linsolve(equations,(a,b)))[0]
        res = simplify(a*n + b)
        return imageset(Lambda(n, factor(res)), S.Naturals)
