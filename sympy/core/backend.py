import os
USE_SYMENGINE = os.getenv('USE_SYMENGINE', '0')
USE_SYMENGINE = USE_SYMENGINE.lower() in ('1', 't', 'true')  # type: ignore

if USE_SYMENGINE:
    from symengine import (Symbol, Integer, sympify, S,
        SympifyError, exp, log, gamma, sqrt, I, E, pi, Matrix,
        sin, cos, tan, cot, csc, sec, asin, acos, atan, acot, acsc, asec,
        sinh, cosh, tanh, coth, asinh, acosh, atanh, acoth,
        lambdify, symarray, diff, zeros, eye, diag, ones,
        expand, Function, symbols, var, Add, Mul, Derivative,
        ImmutableMatrix, MatrixBase, Rational, Basic)
    from symengine.lib.symengine_wrapper import gcd as igcd
    from symengine import AppliedUndef

    class Matrix(Matrix):
        def simplify(self):
            super().simplify()
            return self
    class ImmutableMatrix(ImmutableMatrix):
        def simplify(self):
            super().simplify()
            return self
    old_eye, eye = eye, lambda *args: Matrix(old_eye(*args))
    old_zeros, zeros = zeros, lambda *args: Matrix(old_zeros(*args))
else:
    from sympy import (Symbol, Integer, sympify, S,
        SympifyError, exp, log, gamma, sqrt, I, E, pi, Matrix,
        sin, cos, tan, cot, csc, sec, asin, acos, atan, acot, acsc, asec,
        sinh, cosh, tanh, coth, asinh, acosh, atanh, acoth,
        lambdify, symarray, diff, zeros, eye, diag, ones,
        expand, Function, symbols, var, Add, Mul, Derivative,
        ImmutableMatrix, MatrixBase, Rational, Basic, igcd)
    from sympy.core.function import AppliedUndef

__all__ = [
    'Symbol', 'Integer', 'sympify', 'S', 'SympifyError', 'exp', 'log',
    'gamma', 'sqrt', 'I', 'E', 'pi', 'Matrix', 'sin', 'cos', 'tan', 'cot',
    'csc', 'sec', 'asin', 'acos', 'atan', 'acot', 'acsc', 'asec', 'sinh',
    'cosh', 'tanh', 'coth', 'asinh', 'acosh', 'atanh', 'acoth', 'lambdify',
    'symarray', 'diff', 'zeros', 'eye', 'diag', 'ones', 'expand', 'Function',
    'symbols', 'var', 'Add', 'Mul', 'Derivative', 'ImmutableMatrix',
    'MatrixBase', 'Rational', 'Basic', 'igcd', 'AppliedUndef',
]
