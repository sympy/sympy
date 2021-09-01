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
else:
    from sympy import (Symbol, Integer, sympify, S,
        SympifyError, exp, log, gamma, sqrt, I, E, pi, Matrix,
        sin, cos, tan, cot, csc, sec, asin, acos, atan, acot, acsc, asec,
        sinh, cosh, tanh, coth, asinh, acosh, atanh, acoth,
        lambdify, symarray, diff, zeros, eye, diag, ones,
        expand, Function, symbols, var, Add, Mul, Derivative,
        ImmutableMatrix, MatrixBase, Rational, Basic, igcd)
    from sympy.core.function import AppliedUndef


#
# XXX: Handling of immutable and mutable matrices in SymEngine is inconsistent
# with SymPy's matrix classes in at least SymEngine version 0.7.0. Until that
# is fixed the function below is needed for consistent behaviour when
# attempting to simplify a matrix.
#
# Expected behaviour of a SymPy mutable/immutable matrix .simplify() method:
#
#   Matrix.simplify() : works in place, returns None
#   ImmutableMatrix.simplify() : returns a simplified copy
#
# In SymEngine both mutable and immutable matrices simplify in place and return
# None. This is inconsistent with the matrix being "immutable" and also the
# returned None leads to problems in the mechanics module.
#
# The simplify function should not be used because simplify(M) sympifies the
# matrix M and the SymEngine matrices all sympify to SymPy matrices. If we want
# to work with SymEngine matrices then we need to use their .simplify() method
# but that method does not work correctly with immutable matrices.
#
# The _simplify_matrix function can be removed when the SymEngine bug is fixed.
# Since this should be a temporary problem we do not make this function part of
# the public API.
#
#   SymEngine issue: https://github.com/symengine/symengine.py/issues/363
#

def _simplify_matrix(M):
    """Return a simplified copy of the matrix M"""
    assert isinstance(M, (Matrix, ImmutableMatrix))
    Mnew = M.as_mutable() # makes a copy if mutable
    Mnew.simplify()
    if isinstance(M, ImmutableMatrix):
        Mnew = Mnew.as_immutable()
    return Mnew


__all__ = [
    'Symbol', 'Integer', 'sympify', 'S', 'SympifyError', 'exp', 'log',
    'gamma', 'sqrt', 'I', 'E', 'pi', 'Matrix', 'sin', 'cos', 'tan', 'cot',
    'csc', 'sec', 'asin', 'acos', 'atan', 'acot', 'acsc', 'asec', 'sinh',
    'cosh', 'tanh', 'coth', 'asinh', 'acosh', 'atanh', 'acoth', 'lambdify',
    'symarray', 'diff', 'zeros', 'eye', 'diag', 'ones', 'expand', 'Function',
    'symbols', 'var', 'Add', 'Mul', 'Derivative', 'ImmutableMatrix',
    'MatrixBase', 'Rational', 'Basic', 'igcd', 'AppliedUndef',
]
