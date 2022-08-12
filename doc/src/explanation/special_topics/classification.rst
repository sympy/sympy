===============================
Classification of SymPy objects
===============================

There are several ways of how SymPy object is classified.

class
=====

Like any other object in Python, SymPy expression is an instance of class. You can
get the class of the object with built-in `type()` function, and check it with
`isinstance()` function.

    >>> from sympy import Add
    >>> from sympy.abc import x,y
    >>> type(x + y)
    <class 'sympy.core.add.Add'>
    >>> isinstance(x + y, Add)
    True

Classes represent only the programmatic structures of the objects, and does not
distinguish the mathematical difference between them. For example, the integral
of number and the integral of matrix both have the class `Integral`, although the
former is number and the latter is matrix.

    >>> from sympy import MatrixSymbol, Integral
    >>> A = MatrixSymbol('A', 2, 2)
    >>> type(Integral(1, x))
    <class 'sympy.integrals.integrals.Integral'>
    >>> type(Integral(A, x))
    <class 'sympy.integrals.integrals.Integral'>

.. _kind_classification:

kind
====

Kind indicates what mathematical object does the expression represent.
You can retrieve the kind of expression with `.kind` property.

    >>> Integral(1, x).kind
    NumberKind
    >>> Integral(A, x).kind
    MatrixKind(NumberKind)

This result shows that `Integral(1, x)` is number, and `Integral(A, x)` is matrix with number element.

Since the class cannot guarantee to catch this difference, kind of the object is very important.
For example, if you are building a function or class that is designed to work only for
numbers, you should consider filtering the arguments with `NumberKind` so that the user
does not naively pass unsupported objects such as `Integral(A, x)`.

For the performance, set theory is not implemented in kind system. For example,

    `NumberKind` does not distinguish the real number and complex number.

    >>> from sympy import pi, I
    >>> pi.kind
    NumberKind
    >>> I.kind
    NumberKind

    SymPy's `Set` and kind are not compatible.

    >>> from sympy import S
    >>> from sympy.core.kind import NumberKind
    >>> S.Reals.is_subset(S.Complexes)
    True
    >>> S.Reals.is_subset(NumberKind)
    Traceback (most recent call last):
    ...
    ValueError: Unknown argument 'NumberKind'

sets and assumptions
====================

If you want to classify the object in strictly mathematical way, you may need
SymPy's sets and assumptions.

    >>> from sympy import ask, Q
    >>> S.One in S.Reals
    True
    >>> ask(Q.even(2*x), Q.odd(x))
    True

See `assumptions` module and `sets` module for more information.

func
====

`func` is the head of the object, and it is used to recurse over the expression tree.

    >>> Add(x + y).func
    <class 'sympy.core.add.Add'>
    >>> Add(x + x).func
    <class 'sympy.core.mul.Mul'>
    >>> Q.even(x).func
    <class 'sympy.assumptions.assume.AppliedPredicate'>

As you can see, resulting head may be a class or another SymPy object.
Keep this in mind when you classify the object with this attribute.
See :ref:`tutorial-manipulation` for detailed information.
