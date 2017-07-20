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
        ImmutableMatrix, MatrixBase, Rational, Basic, Pow, oo, zoo, nan,
        Float, perfect_power, sign, ceiling, EulerGamma)
    from symengine.lib.symengine_wrapper import gcd, gcd as igcd, Piecewise
    from symengine import AppliedUndef
else:
    from sympy import (Symbol, Integer, sympify, S,
        SympifyError, exp, log, gamma, sqrt, I, E, pi, Matrix,
        sin, cos, tan, cot, csc, sec, asin, acos, atan, acot, acsc, asec,
        sinh, cosh, tanh, coth, asinh, acosh, atanh, acoth,
        lambdify, symarray, diff, zeros, eye, diag, ones, zeros,
        expand, Function, symbols, var, Add, Mul, Derivative,
        ImmutableMatrix, MatrixBase, Rational, Basic, igcd, Pow, oo, zoo, nan,
        gcd, Float, sign, ceiling, EulerGamma)
    from sympy.core.function import AppliedUndef
    from sympy.ntheory.factor_ import perfect_power
