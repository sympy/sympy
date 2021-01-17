"""
Module to efficiently partition SymPy objects.

This system is introduced because class of SymPy object does not always
represent the mathematical classification of the entity. For example,
``Integral(1, x)`` and ``Integral(Matrix([1,2]), x)`` are both instance
of ``Integral`` class. However the former is number and the latter is
matrix.

One way to resolve this is defining subclass for each mathematical type,
such as ``MatAdd`` for the addition between matrices. Basic algebraic
operation such as addition or multiplication take this approach, but
defining every class for every mathematical object is not scalable.

Therefore, we define the "kind" of the object and let the expression
infer the kind of itself from its arguments. Function and class can
filter the arguments by their kind, and behave differently according to
the type of itself.

This module defines basic kinds for core objects. Other kinds such as
``ArrayKind`` or ``MatrixKind`` can be found in corresponding modules.

.. notes::
       This approach is experimental, and can be replaced or deleted in the future.
       See https://github.com/sympy/sympy/pull/20549.
"""


class KindMeta(type):
    """
    Metaclass for ``Kind``.

    Assigns empty ``dict`` as class attribute ``_inst`` for every class,
    in order to endow singleton-like behavior.
    """
    def __new__(cls, clsname, bases, dct):
        dct['_inst'] = {}
        return super().__new__(cls, clsname, bases, dct)


class Kind(object, metaclass=KindMeta):
    """
    Base class for kinds.

    Kind of the object represents the mathematical classification that
    the entity falls into. It is expected that functions and classes
    recognize and filter the argument by its kind.

    Kind of every object must be carefully selected so that it shows the
    intention of design. Expressions may have different kind according
    to the kind of its arguements. For example, arguements of ``Add``
    must have common kind since addition is group operator, and the
    resulting ``Add()`` has the same kind.

    For the performance, each kind is as broad as possible and is not
    based on set theory. For example, ``NumberKind`` includes not only
    complex number but expression containing ``S.Infinity`` or ``S.NaN``
    which are not strictly number.

    Kind may have arguments as parameter. For example, ``MatrixKind()``
    may be constructed with one element which represents the kind of its
    elements.

    ``Kind`` behaves in singleton-like fashion. Same signature will
    return the same object.

    """
    def __new__(cls, *args):
        if args in cls._inst:
            inst = cls._inst[args]
        else:
            inst = super().__new__(cls)
            cls._inst[args] = inst
        return inst


class UndefinedKind(Kind):
    """
    Default kind for all SymPy object. If the kind is not defined for
    the object, or if the object cannot infer the kind from its
    arguments, this will be returned.

    Examples
    ========

    >>> from sympy import Expr
    >>> Expr().kind
    UndefinedKind
    """
    def __new__(cls):
        return super().__new__(cls)

    def __repr__(self):
        return "UndefinedKind"

UndefinedKind = UndefinedKind()


class NumberKind(Kind):
    """
    Kind for all numeric object.

    This kind represents every number, including complex numbers,
    infinity and ``S.NaN``. Other objects such as quaternions do not
    have this kind.

    Most ``Expr`` are initially designed to represent the number, so
    this will be the most common kind in SymPy core. For example
    ``Symbol()``, which represents a scalar, has this kind as long as it
    is commutative.

    Numbers form a field. Any operation between number-kind objects will
    result this kind as well.

    Examples
    ========

    >>> from sympy import S, oo, Symbol
    >>> S.One.kind
    NumberKind
    >>> (-oo).kind
    NumberKind
    >>> S.NaN.kind
    NumberKind

    Commutative symbol are treated as number.

    >>> x = Symbol('x')
    >>> x.kind
    NumberKind
    >>> Symbol('y', commutative=False).kind
    UndefinedKind

    Operation between numbers results number.

    >>> (x+1).kind
    NumberKind

    See Also
    ========

    sympy.core.expr.Expr.is_Number : check if the object is strictly
    subclass of ``Number`` class.

    sympy.core.expr.Expr.is_number : check if the object is number
    without any free symbol.

    """
    def __new__(cls):
        return super().__new__(cls)

    def __repr__(self):
        return "NumberKind"

NumberKind = NumberKind()


class BooleanKind(Kind):
    """
    Kind for boolean objects.

    SymPy's ``S.true``, ``S.false``, and built-in ``True`` and ``False``
    have this kind. Boolean number ``1`` and ``0`` are not relevent.

    Examples
    ========

    >>> from sympy import S, Q
    >>> S.true.kind
    BooleanKind
    >>> Q.even(3).kind
    BooleanKind
    """
    def __new__(cls):
        return super().__new__(cls)

    def __repr__(self):
        return "BooleanKind"

BooleanKind = BooleanKind()
