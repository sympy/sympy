import os
USE_SYMENGINE = os.getenv('USE_SYMENGINE')

if USE_SYMENGINE:
    import symengine as sm
    from symengine import (Symbol, Integer, sympify, sympify as S,
        SympifyError, exp, log, gamma, sqrt, I, E, pi, Matrix,
        sin, cos, tan, cot, csc, sec, asin, acos, atan, acot, acsc, asec,
        sinh, cosh, tanh, coth, asinh, acosh, atanh, acoth,
        Lambdify as lambdify, symarray, diff, zeros, eye, diag, ones, zeros,
        expand, UndefFunction as Function, symbols, var)
    from symengine.lib.symengine_wrapper import (Function as FunctionClass,
        Add as AddClass, Mul as MulClass, add as Add, mul as Mul)
else:
    import sympy as sm
    from sympy import (Symbol, Integer, sympify, S, SympifyError,
        Add, Mul, exp, log, gamma, sqrt, I, E, pi, Matrix,
        sin, cos, tan, cot, csc, sec, asin, acos, atan, acot, acsc, asec,
        sinh, cosh, tanh, coth, asinh, acosh, atanh, acoth, lambdify,
        symarray, diff, zeros, eye, diag, ones, zeros, expand,
        Function, symbols, var)
    from sympy import (Function as FunctionClass, Add as AddClass,
        Mul as MulClass)
