.. _booleans-guide:

=============================
 Symbolic and fuzzy booleans
=============================


This page describes what a symbolic :class:`~.Boolean` in SymPy is and also
how that relates to three-valued fuzzy-bools that are used in many parts of
SymPy. It also discusses some common problems that arise when writing code
that uses three-valued logic and how to handle them correctly.


Symbolic Boolean vs three valued bool
=====================================

Assumptions queries like ``x.ispositive`` give fuzzy-bool ``True``,
``False`` or ``None`` results [#fuzzy]_. These are low-level Python objects
rather than SymPy's symbolic :class:`~.Boolean` expressions.

    >>> from sympy import Symbol, symbols
    >>> xpos = Symbol('xpos', positive=True)
    >>> xneg = Symbol('xneg', negative=True)
    >>> x = Symbol('x')
    >>> print(xpos.is_positive)
    True
    >>> print(xneg.is_positive)
    False
    >>> print(x.is_positive)
    None

A ``None`` result as a fuzzy-bool should be interpreted as meaning "maybe" or
"unknown".

An example of a symbolic :class:`~.Boolean` class in SymPy can be found when
using inequalities. When an inequality is not known to be true or false a
:class:`~.Boolean` can represent indeterminate results symbolically:

    >>> xpos > 0
    True
    >>> xneg > 0
    False
    >>> x > 0
    x > 0
    >>> type(x > 0)
    <class 'sympy.core.relational.StrictGreaterThan'>

The last example shows what happens when an inequality is indeterminate: we
get an instance of :class:`~.StrictGreaterThan` which represents the
inequality as a symbolic expression. Internally when attempting to evaluate an
inequality like ``a > b`` SymPy will compute ``(a - b).is_extended_positive``.
If the result is ``True`` or ``False`` then SymPy's symbolic ``S.true`` or
``S.false`` will be returned. If the result is ``None`` then an unevaluated
:class:`~.StrictGreaterThan` is returned as shown for ``x > 0`` above.

It is not obvious that queries like ``xpos > 0`` return ``S.true`` rather than
``True`` because both objects display in the same way but we can check this
using the Python ``is`` operator:

    >>> from sympy import S
    >>> xpos.is_positive is True
    True
    >>> xpos.is_positive is S.true
    False
    >>> (xpos > 0) is True
    False
    >>> (xpos > 0) is S.true
    True

There is no general symbolic analogue of ``None`` in SymPy. In the cases where
a low-level assumptions query gives ``None`` the symbolic query will result in
an unevaluated symbolic :class:`~.Boolean` (e.g, ``x > 0``).  We can use a
symbolic :class:`~.Boolean` as part of a symbolic expression such as a
:class:`~.Piecewise`:

    >>> from sympy import Piecewise
    >>> p = Piecewise((1, x > 0), (2, True))
    >>> p
    Piecewise((1, x > 0), (2, True))
    >>> p.subs(x, 3)
    1

Here ``p`` represents an expression that will be equal to ``1`` if ``x > 0``
or otherwise it will be equal to ``2``. The unevaluated :class:`~.Boolean` inequality
``x > 0`` represents the condition for deciding the value of the expression
symbolically. When we substitute a value for ``x`` the inequality will resolve
to ``S.true`` and then the :class:`~.Piecewise` can evaluate to ``1`` or ``2``.

The same will not work when using a fuzzy-bool instead of a symbolic
:class:`~.Boolean`:

    >>> p2 = Piecewise((1, x.is_positive), (2, True))
    Traceback (most recent call last):
    ...
    TypeError: Second argument must be a Boolean, not `NoneType`

The :class:`~.Piecewise` can not use ``None`` as the condition because unlike the
inequality ``x > 0`` it gives no information. With the inequality it is
possible to decide in future if the condition might ``True`` or ``False``
once a value for ``x`` is known. A value of ``None`` can not be used in that
way so it is rejected.

.. note:: We can use ``True`` in the :class:`~.Piecewise` because ``True`` sympifies
          to ``S.true``. Sympifying ``None`` just gives ``None`` again which
          is not a valid symbolic SymPy object.

There are many other symbolic :class:`~.Boolean` types in SymPy. The same
considerations about the differences between fuzzy bool and symbolic
:class:`~.Boolean` apply to all other SymPy :class:`~.Boolean` types. To give
a different example there is :class:`~.Contains` which represents the
statement that an object is contained in a set:

    >>> from sympy import Reals, Contains
    >>> x = Symbol('x', real=True)
    >>> y = Symbol('y')
    >>> Contains(x, Reals)
    True
    >>> Contains(y, Reals)
    Contains(y, Reals)
    >>> Contains(y, Reals).subs(y, 1)
    True

The Python operator corresponding to :class:`~.Contains` is ``in``. A quirk of
``in`` is that it can only evaluate to a ``bool`` (``True`` or ``False``) so
if the result is indeterminate then an exception will be raised:

    >>> from sympy import I
    >>> 2 in Reals
    True
    >>> I in Reals
    False
    >>> x in Reals
    True
    >>> y in Reals
    Traceback (most recent call last):
    ...
    TypeError: did not evaluate to a bool: (-oo < y) & (y < oo)

The exception can be avoided by using ``Contains(x, Reals)`` or
``Reals.contains(x)`` rather than ``x in Reals``.


Three-valued logic with fuzzy bools
===================================

Whether we use the fuzzy-bool or symbolic :class:`~.Boolean` we always need to be
aware of the possibility that a query might be indeterminate. How to write
code that handles this is different in the two cases though. We will look at
fuzzy-bools first.

Consider the following function:

    >>> def both_positive(a, b):
    ...     """ask whether a and b are both positive"""
    ...     if a.is_positive and b.is_positive:
    ...         return True
    ...     else:
    ...         return False

The ``both_positive`` function is supposed to tell us whether or not ``a`` and
``b`` are both positive. However the ``both_positive`` function will fail if
either of the ``is_positive`` queries gives ``None``:

    >>> print(both_positive(S(1), S(1)))
    True
    >>> print(both_positive(S(1), S(-1)))
    False
    >>> print(both_positive(S(-1), S(-1)))
    False
    >>> x = Symbol('x') # may or may not be positive
    >>> print(both_positive(S(1), x))
    False

.. note:: We need to sympify the arguments to this function using ``S``
          because the assumptions are only defined on SymPy objects and not
          regular Python :class:`int` objects.

Here ``False`` is incorrect because it is *possible* that ``x`` is positive in
which case both arguments would be positive. We get ``False`` here because
``x.is_positive`` gives ``None`` and Python will treat ``None`` as "falsey".

In order to handle all possible cases correctly we need to separate the logic
for identifying the ``True`` and ``False`` cases. An improved function might
be:

    >>> def both_positive_better(a, b):
    ...     """ask whether a and b are both positive"""
    ...     if a.is_positive is False or b.is_positive is False:
    ...         return False
    ...     elif a.is_positive is True and b.is_positive is True:
    ...         return True
    ...     else:
    ...         return None

This function now can handle all cases of ``True``, ``False`` or ``None`` for
both ``a`` and ``b`` and will always return a fuzzy bool representing whether
the statement "``a`` and ``b`` are both positive" is true, false or unknown:

    >>> print(both_positive_better(S(1), S(1)))
    True
    >>> print(both_positive_better(S(1), S(-1)))
    False
    >>> x = Symbol('x')
    >>> y = Symbol('y', positive=True)
    >>> print(both_positive_better(S(1), x))
    None
    >>> print(both_positive_better(S(-1), x))
    False
    >>> print(both_positive_better(S(1), y))
    True

Another case that we need to be careful of when using fuzzy-bools is negation
with Python's ``not`` operator e.g.:

    >>> x = Symbol('x')
    >>> print(x.is_positive)
    None
    >>> not x.is_positive
    True

The correct negation of a fuzzy bool ``None`` is ``None`` again. If we do not
know whether the statement "``x`` is positive" is ``True`` or ``False`` then
we also do not know whether its negation "``x`` is not positive" is ``True``
or ``False``. The reason we get ``True`` instead is again because ``None`` is
considered "falsey". When ``None`` is used with a logical operator such as
``not`` it will first be converted to a :class:`bool` and then negated:

    >>> bool(None)
    False
    >>> not bool(None)
    True
    >>> not None
    True

The fact that ``None`` is treated as falsey can be useful if used correctly.
For example we may want to do something only if ``x`` is known to positive in
which case we can do

    >>> x = Symbol('x', positive=True)
    >>> if x.is_positive:
    ...     print("x is definitely positive")
    ... else:
    ...     print("x may or may not be positive")
    x is definitely positive

Provided we understand that an alternate condition branch refers to two cases
(``False`` and ``None``) then this can be a useful way of writing
conditionals.  When we really do need to distinguish all cases then we need to
use things like ``x.is_positive is False``.  What we need to be careful of
though is using Python's binary logic operators like ``not`` or ``and`` with
fuzzy bools as they will not handle the indeterminate case correctly.

In fact SymPy has internal functions that are designed to handle fuzzy-bools
correctly:

    >>> from sympy.core.logic import fuzzy_not, fuzzy_and
    >>> print(fuzzy_not(True))
    False
    >>> print(fuzzy_not(False))
    True
    >>> print(fuzzy_not(None))
    None
    >>> print(fuzzy_and([True, True]))
    True
    >>> print(fuzzy_and([True, None]))
    None
    >>> print(fuzzy_and([False, None]))
    False

Using the ``fuzzy_and`` function we can write the ``both_positive`` function
more simply:

    >>> def both_positive_best(a, b):
    ...     """ask whether a and b are both positive"""
    ...     return fuzzy_and([a.is_positive, b.is_positive])

Making use of ``fuzzy_and``, ``fuzzy_or`` and ``fuzzy_not`` leads to simpler
code and can also reduce the chance of introducing a logic error because the
code can look more like it would in the case of ordinary binary logic.


Three-valued logic with symbolic Booleans
=========================================

When working with symbolic :class:`~.Boolean` rather than fuzzy-bool the issue of
``None`` silently being treated as falsey does not arise so it is easier not
to end up with a logic error. However instead the indeterminate case will
often lead to an exception being raised if not handled carefully.

We will try to implement the ``both_positive`` function this time using
symbolic :class:`~.Boolean`:

    >>> def both_positive(a, b):
    ...     """ask whether a and b are both positive"""
    ...     if a > 0 and b > 0:
    ...         return S.true
    ...     else:
    ...         return S.false

The first difference is that we return the symbolic :class:`~.Boolean` objects
``S.true`` and ``S.false`` rather than ``True`` and ``False``. The second
difference is that we test e.g. ``a > 0`` rather than ``a.is_positive``.
Trying this out we get

    >>> both_positive(1, 2)
    True
    >>> both_positive(-1, 1)
    False
    >>> x = Symbol('x')  # may or may not be positive
    >>> both_positive(x, 1)
    Traceback (most recent call last):
    ...
    TypeError: cannot determine truth value of Relational

What happens now is that testing ``x > 0`` gives an exception when ``x`` is
not known to be positive or not positive. More precisely ``x > 0`` does not
give an exception but ``if x > 0`` does and that is because the ``if``
statement implicitly calls ``bool(x > 0)`` which raises.

    >>> x > 0
    x > 0
    >>> bool(x > 0)
    Traceback (most recent call last):
    ...
    TypeError: cannot determine truth value of Relational
    >>> if x > 0:
    ...     print("x is positive")
    Traceback (most recent call last):
    ...
    TypeError: cannot determine truth value of Relational

The Python expression ``x > 0`` creates a SymPy :class:`~.Boolean`. Since in this
case the :class:`~.Boolean` can not evaluate to ``True`` or ``False`` we get an
unevaluated :class:`~.StrictGreaterThan`. Attempting to force that into a
``bool`` with ``bool(x > 0)`` raises an exception. That is because a regular
Python ``bool`` must be either ``True`` or ``False`` and neither of those
are known to be correct in this case.

The same kind of issue arises when using ``and``, ``or`` or ``not`` with
symbolic :class:`~.Boolean`. The solution is to use SymPy's symbolic
:class:`~.And`, :class:`~.Or` and :class:`~.Not` or equivalently Python's
bitwise logical operators ``&``, ``|`` and ``~``:

    >>> from sympy import And, Or, Not
    >>> x > 0
    x > 0
    >>> x > 0 and x < 1
    Traceback (most recent call last):
    ...
    TypeError: cannot determine truth value of Relational
    >>> And(x > 0, x < 1)
    (x > 0) & (x < 1)
    >>> (x > 0) & (x < 1)
    (x > 0) & (x < 1)
    >>> Or(x < 0, x > 1)
    (x > 1) | (x < 0)
    >>> Not(x < 0)
    x >= 0
    >>> ~(x < 0)
    x >= 0

As before we can make a better version of ``both_positive`` if we avoid
directly using a SymPy :class:`~.Boolean` in an ``if``, ``and``, ``or``, or ``not``.
Instead we can test whether or not the :class:`~.Boolean` has evaluated to ``S.true``
or ``S.false``:

    >>> def both_positive_better(a, b):
    ...     """ask whether a and b are both positive"""
    ...     if (a > 0) is S.false or (b > 0) is S.false:
    ...         return S.false
    ...     elif (a > 0) is S.true and (b > 0) is S.true:
    ...         return S.true
    ...     else:
    ...         return And(a > 0, b > 0)

Now with this version we don't get any exceptions and if the result is
indeterminate we will get a symbolic :class:`~.Boolean` representing the conditions
under which the statement "``a`` and ``b`` are both positive" would be true:

    >>> both_positive_better(S(1), S(2))
    True
    >>> both_positive_better(S(1), S(-1))
    False
    >>> x, y = symbols("x, y")
    >>> both_positive_better(x, y + 1)
    (x > 0) & (y + 1 > 0)
    >>> both_positive_better(x, S(3))
    x > 0

The last case shows that actually using the :class:`~.And` with a condition that is
known to be true simplifies the :class:`~.And`. In fact we have

    >>> And(x > 0, 3 > 0)
    x > 0
    >>> And(4 > 0, 3 > 0)
    True
    >>> And(-1 > 0, 3 > 0)
    False

What this means is that we can improve ``both_positive_better``. The
different cases are not needed at all. Instead we can simply return the
:class:`~.And` and let it simplify if possible:

    >>> def both_positive_best(a, b):
    ...     """ask whether a and b are both positive"""
    ...     return And(a > 0, b > 0)

Now this will work with any symbolic real objects and produce a symbolic
result. We can also substitute into the result to see how it would work for
particular values:

    >>> both_positive_best(2, 1)
    True
    >>> both_positive_best(-1, 2)
    False
    >>> both_positive_best(x, 3)
    x > 0
    >>> condition = both_positive_best(x/y, x + y)
    >>> condition
    (x + y > 0) & (x/y > 0)
    >>> condition.subs(x, 1)
    (1/y > 0) & (y + 1 > 0)
    >>> condition.subs(x, 1).subs(y, 2)
    True

The idea when working with symbolic :class:`~.Boolean` objects is as much as possible
to avoid trying to branch on them with ``if/else`` and other logical operators
like ``and`` etc. Instead think of computing a condition and passing it around
as a variable. The elementary symbolic operations like :class:`~.And`,
:class:`~.Or` and :class:`~.Not` can then take care of the logic for you.

.. rubric:: Footnotes

.. [#fuzzy] Note that what is referred to in SymPy as a "fuzzy bool" is really
   about using three-valued logic. In normal usage "fuzzy logic" refers to a
   system where logical values are continuous in between zero and one which is
   something different from three-valued logic.
