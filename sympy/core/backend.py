import os
USE_SYMENGINE = os.getenv('USE_SYMENGINE', '0')
USE_SYMENGINE = USE_SYMENGINE.lower() in ('1', 't', 'true')

if USE_SYMENGINE:
    from symengine import (Symbol, Integer, sympify, S,
        SympifyError, exp, log, gamma, sqrt, I, E, pi, Matrix,
        sin, cos, tan, cot, csc, sec, asin, acos, atan, acot, acsc, asec,
        sinh, cosh, tanh, coth, asinh, acosh, atanh, acoth,
        lambdify, symarray, diff, zeros, eye, diag, ones, zeros,
        expand, Function, symbols, var, Add, Mul, Derivative,
        ImmutableMatrix, MatrixBase, Rational, Basic, atan2, Dummy, Eq,
        oo, sign)
    from symengine.lib.symengine_wrapper import gcd as igcd, And, Piecewise, FiniteSet
    from symengine import AppliedUndef
else:
    from sympy import (Symbol, Integer, sympify, S,
        SympifyError, exp, log, gamma, sqrt, I, E, pi, Matrix,
        sin, cos, tan, cot, csc, sec, asin, acos, atan, acot, acsc, asec,
        sinh, cosh, tanh, coth, asinh, acosh, atanh, acoth,
        symarray, diff, zeros, eye, diag, ones, zeros, expand,
        Function, symbols, var, Add, Mul, Derivative, Set,
        ImmutableMatrix, MatrixBase, Rational, Basic, igcd, And,
        Piecewise, atan2, Dummy, Eq, oo, sign)
    from sympy.core.function import AppliedUndef
    from sympy.utilities import lambdify
    from sympy.sets import FiniteSet
