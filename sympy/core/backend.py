import os
USE_SYMENGINE = os.getenv('USE_SYMENGINE')

if USE_SYMENGINE:
    from symengine.sympy_compat import (Symbol, Integer, sympify, S,
        SympifyError, exp, log, gamma, sqrt, I, E, pi, Matrix,
        sin, cos, tan, cot, csc, sec, asin, acos, atan, acot, acsc, asec,
        sinh, cosh, tanh, coth, asinh, acosh, atanh, acoth,
        lambdify, symarray, diff, zeros, eye, diag, ones, zeros,
        expand, Function, symbols, var, Add, Mul, Derivative)
    from symengine.sympy_compat import AppliedUndef
    #TODO: Fix this
    from symengine.lib.symengine_wrapper import (Matrix as ImmutableMatrix,
        MatrixBase)
else:
    from sympy import (Symbol, Integer, sympify, S,
        SympifyError, exp, log, gamma, sqrt, I, E, pi, Matrix,
        sin, cos, tan, cot, csc, sec, asin, acos, atan, acot, acsc, asec,
        sinh, cosh, tanh, coth, asinh, acosh, atanh, acoth,
        lambdify, symarray, diff, zeros, eye, diag, ones, zeros,
        expand, Function, symbols, var, Add, Mul, Derivative)
    from sympy.core.function import AppliedUndef
    from sympy import ImmutableMatrix, MatrixBase
