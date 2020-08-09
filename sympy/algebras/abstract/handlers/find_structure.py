"""
Finds the most general algebraic structure for the expressions

Examples
========

>>> from sympy import S
>>> from sympy.algebras.abstract.handlers import find_structure

>>> find_structure(S.One) == S.ComplexesField
True

"""

from sympy import S, Expr
from sympy.core.sympify import _sympify
from sympy.map import (
    Map, FunctionSet, ConstantMap,
)
from sympy.algebras import (
    Ring,
    FunctionAdditionOperator, AbelianGroup,
    FunctionScalarMultiplicationOperator,
    VectorSpace, FunctionVectorMultiplicationOperator,
    Algebra,
)
from sympy.multipledispatch import dispatch

__all__ = [
    'find_structure'
]

@dispatch(object)  # type: ignore # noqa:F811
def find_structure(a): # noqa:F811
    return find_structure(_sympify(a))

@dispatch(Expr)  # type: ignore # noqa:F811
def find_structure(a): # noqa:F811
    return S.ComplexesField

@dispatch(Map)  # type: ignore # noqa:F811
def find_structure(a): # noqa:F811
    X, F = a.domain, a.codomain

    if not isinstance(F, Ring):
        raise NotImplementedError("Only operation between scalar functions is supported")

    fs = FunctionSet(domain=X, codomain=F)
    zerofunc = ConstantMap(F.add_op.identity, domain=X)
    fadd = FunctionAdditionOperator(fs**2, fs, zerofunc)
    G = AbelianGroup('', (fs,), (fadd,))
    fsmul = FunctionScalarMultiplicationOperator(F*G, G)
    FS = VectorSpace('', (F, G), (fsmul,))
    ff_mul = FunctionVectorMultiplicationOperator(FS*FS, G)
    A = Algebra('', (FS,), (ff_mul,))
    return A
