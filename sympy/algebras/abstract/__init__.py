__all__ = [
    "AlgebraicStructure",
    "Magma",
    "Semigroup",
    "LeftQuasigroup", "RightQuasigroup", "Quasigroup",
    "Monoid",
    "Loop",
    "Group", 'AbelianGroup',
    "Ring", "CommutativeRing",
    "Field",
    "Module",
    "VectorSpace",
    "Algebra",
    'Module',
    'VectorSpace',
    'FunctionAdditionOperator', 'FunctionAddition',
    'FunctionScalarMultiplicationOperator',
    'FunctionMultiplication',
    'FunctionVectorMultiplicationOperator',
    'FunctionExponent',
]

from .structure import AlgebraicStructure

from .group import (
    Magma, Semigroup,
    LeftQuasigroup, RightQuasigroup, Quasigroup, Monoid, Loop,
    Group, AbelianGroup,
)

from .ring import (
    Ring, CommutativeRing, Field,
)

from .module import (
    Module,
    VectorSpace,
    FunctionAdditionOperator, FunctionAddition,
    FunctionScalarMultiplicationOperator,
    FunctionMultiplication,
    FunctionVectorMultiplicationOperator,
    FunctionExponent,
)

from .algebra import (
    Algebra,
)

### Defined structures

from sympy.core.singleton import S, Singleton
from sympy.map.add import NumericAdditionOperator
from sympy.map.mul import NumericMultiplicationOperator

class IntegersRing(CommutativeRing, metaclass=Singleton):
    """
    The ring of integers.

    Examples
    ========

    >>> from sympy import S
    >>> S.IntegersRing
    Z

    >>> x = S.IntegersRing.element('x')

    >>> x.is_integer
    True
    >>> S.IntegersRing.add(x, x, evaluate=True)
    2*x
    >>> S.IntegersRing.sub(x, x, evaluate=True)
    0
    >>> S.IntegersRing.mul(x, -1, evaluate=True)
    (-1)*x
    >>> S.IntegersRing.pow(x, 4, evaluate=True)
    x**4

    """

    def __new__(cls, *args, **kwargs):
        name = 'Z'
        sets = (S.Integers,)
        operators = (
            NumericAdditionOperator(S.Integers**2, S.Integers),
            NumericMultiplicationOperator(S.Integers**2, S.Integers)
        )
        return super().__new__(cls, name, sets, operators)

class RealsField(Field, metaclass=Singleton):
    """
    The field of real numbers.

    Examples
    ========

    >>> from sympy import S
    >>> S.RealsField
    R

    >>> x = S.RealsField.element('x')

    >>> x.is_real
    True
    >>> S.RealsField.add(x, x, evaluate=True)
    2*x
    >>> S.RealsField.sub(x, x, evaluate=True)
    0
    >>> S.RealsField.mul(x, -1, evaluate=True)
    (-1)*x
    >>> S.RealsField.pow(x, 4, evaluate=True)
    x**4
    >>> S.RealsField.div(x, 2, evaluate=True)
    x*(2**(-1))

    """

    def __new__(cls, *args, **kwargs):
        name = 'R'
        sets = (S.Reals,)
        operators = (
            NumericAdditionOperator(S.Reals**2, S.Reals),
            NumericMultiplicationOperator(S.Reals**2, S.Reals)
        )
        return super().__new__(cls, name, sets, operators)

class ComplexesField(Field, metaclass=Singleton):
    """
    The field of complex numbers.

    Examples
    ========

    >>> from sympy import S
    >>> S.ComplexesField
    C

    >>> x = S.ComplexesField.element('x')

    >>> x.is_complex
    True
    >>> S.ComplexesField.add(x, x, evaluate=True)
    2*x
    >>> S.ComplexesField.sub(x, x, evaluate=True)
    0
    >>> S.ComplexesField.mul(x, -1, evaluate=True)
    (-1)*x
    >>> S.ComplexesField.pow(x, 4, evaluate=True)
    x**4
    >>> S.ComplexesField.div(x, 2, evaluate=True)
    x*(2**(-1))

    """

    def __new__(cls, *args, **kwargs):
        name = 'C'
        sets = (S.Complexes,)
        operators = (
            NumericAdditionOperator(S.Complexes**2, S.Complexes),
            NumericMultiplicationOperator(S.Complexes**2, S.Complexes)
        )
        return super().__new__(cls, name, sets, operators)
