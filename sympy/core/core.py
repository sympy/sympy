""" The core's core. """
from __future__ import annotations


# used for canonical ordering of symbolic sequences
# via __cmp__ method:
# FIXME this is *so* irrelevant and outdated!
ordering_of_classes = [
    # singleton numbers
    'Zero', 'One', 'Half', 'Infinity', 'NaN', 'NegativeOne', 'NegativeInfinity',
    # numbers
    'Integer', 'Rational', 'Float',
    # singleton symbols
    'Exp1', 'Pi', 'ImaginaryUnit',
    # symbols
    'Symbol', 'Wild', 'Temporary',
    # arithmetic operations
    'Pow', 'Mul', 'Add',
    # function values
    'Derivative', 'Integral',
    # defined singleton functions
    'Abs', 'Sign', 'Sqrt',
    'Floor', 'Ceiling',
    'Re', 'Im', 'Arg',
    'Conjugate',
    'Exp', 'Log',
    'Sin', 'Cos', 'Tan', 'Cot', 'ASin', 'ACos', 'ATan', 'ACot',
    'Sinh', 'Cosh', 'Tanh', 'Coth', 'ASinh', 'ACosh', 'ATanh', 'ACoth',
    'RisingFactorial', 'FallingFactorial',
    'factorial', 'binomial',
    'Gamma', 'LowerGamma', 'UpperGamma', 'PolyGamma',
    'Erf',
    # special polynomials
    'Chebyshev', 'Chebyshev2',
    # undefined functions
    'Function', 'WildFunction',
    # anonymous functions
    'Lambda',
    # Landau O symbol
    'Order',
    # relational operations
    'Equality', 'Unequality', 'StrictGreaterThan', 'StrictLessThan',
    'GreaterThan', 'LessThan',
]


class Registry:
    """
    Base class for registry objects.

    Registries map a name to an object using attribute notation. Registry
    classes behave singletonically: all their instances share the same state,
    which is stored in the class object.

    All subclasses should set `__slots__ = ()`.
    """
    __slots__ = ()

    def __setattr__(self, name, obj):
        setattr(self.__class__, name, obj)

    def __delattr__(self, name):
        delattr(self.__class__, name)
