=======================================
Development Tips: Comparisons in Python
=======================================

.. role:: input(strong)

Introduction
============

When debugging comparisons and hashes in SymPy, it is necessary to understand
when exactly Python calls each method.
Unfortunately, the official Python documentation for this is
not very detailed (see the docs for `rich comparison
<http://docs.python.org/dev/reference/datamodel.html#object.__lt__>`_,
`__cmp__() <http://docs.python.org/dev/reference/datamodel.html#object.__cmp__>`_
and `__hash__()
<http://docs.python.org/dev/reference/datamodel.html#object.__hash__>`_
methods).

We wrote this guide to fill in the missing gaps. After reading it, you should
be able to understand which methods do (and do not) get called and the order in
which they are called.

Hashing
=======

Every Python class has a ``__hash__()`` method, the default
implementation of which is::

    def __hash__(self):
        return id(self)

You can reimplement it to return a different integer that you compute on your own.
``hash(x)`` just calls ``x.__hash__()``. Python builtin classes usually redefine
the ``__hash__()`` method. For example, an ``int`` has something like this::

    def __hash__(self):
        return int(self)

and a ``list`` does something like this::

    def __hash__(self):
        raise TypeError("list objects are unhashable")

The general
idea about hashes is that if two objects have a different hash, they are not
equal, but if they have the same hash, they *might* be equal. (This is usually
called a "hash collision" and you need to use the methods described in the
next section to determine if the objects really are equal).

The only requirement from the Python side is
that the hash value mustn't change after it is returned by the
``__hash__()`` method.

Please be aware that hashing is *platform-dependent*. This means that you can
get different hashes for the same SymPy object on different platforms. This
affects for instance sorting of sympy expressions. You can also get SymPy
objects printed in different order.

When developing, you have to be careful about this, especially when writing
tests. It is possible that your test runs on a 32-bit platform, but not on
64-bit. An example::

    >> from sympy import *
    >> x = Symbol('x')
    >> r = rootfinding.roots_quartic(Poly(x**4 - 6*x**3 + 17*x**2 - 26*x + 20, x))
    >> [i.evalf(2) for i in r]
    [1.0 + 1.7*I, 2.0 - 1.0*I, 2.0 + I, 1.0 - 1.7*I]

If you get this order of solutions, you are probably running 32-bit system.
On a 64-bit system you would get the following::

    >> [i.evalf(2) for i in r]
    [1.0 - 1.7*I, 1.0 + 1.7*I, 2.0 + I, 2.0 - 1.0*I

When you now write a test like this::

    r = [i.evalf(2) for i in r]
    assert r == [1.0 + 1.7*I, 2.0 - 1.0*I, 2.0 + I, 1.0 - 1.7*I]

it will fail on a 64-bit platforms, even if it works for your 32-bit system. You can
avoid this by using the ``sorted()`` or ``set()`` Python built-in::

    r = [i.evalf(2) for i in r]
    assert set(r) == set([1.0 + 1.7*I, 2.0 - 1.0*I, 2.0 + I, 1.0 - 1.7*I])

This approach does not work for doctests since they always compare strings that would
be printed after a prompt. In that case you could make your test print results using
a combination of ``str()`` and ``sorted()``::

    >> sorted([str(i.evalf(2)) for i in r])
    ['1.0 + 1.7*I', '1.0 - 1.7*I', '2.0 + I', '2.0 - 1.0*I']

or, if you don't want to show the values as strings, then sympify the results or the
sorted list::

    >> [S(s) for s in sorted([str(i.evalf(2)) for i in r])]
    [1.0 + 1.7*I, 1.0 - 1.7*I, 2.0 + I, 2.0 - I]

The printing of SymPy expressions might be also affected, so be careful
with doctests. If you get the following on a 32-bit system::

    >> print dsolve(f(x).diff(x, 2) + 2*f(x).diff(x) - f(x), f(x))
    f(x) == C1*exp(-x + x*sqrt(2)) + C2*exp(-x - x*sqrt(2))

you might get the following on a 64-bit platform::

    >> print dsolve(f(x).diff(x, 2) + 2*f(x).diff(x) - f(x), f(x))
    f(x) == C1*exp(-x - x*sqrt(2)) + C2*exp(-x + x*sqrt(2))

Method Resolution
=================

Let ``a``, ``b`` and ``c`` be instances of any one of the Python classes.
As can be easily checked by the `Python script`_ at the end of this guide,
if you write::

    a == b

Python calls the following -- in this order::

    a.__eq__(b)
    b.__eq__(a)
    a.__cmp__(b)
    b.__cmp__(a)
    id(a) == id(b)

If a particular method is not implemented (or a method
returns ``NotImplemented`` [1]_) Python skips it
and tries the next one until it succeeds (i.e. until the method returns something
meaningful). The last line is a catch-all method that always succeeds.

If you write::

    a != b

Python tries to call::

    a.__ne__(b)
    b.__ne__(a)
    a.__cmp__(b)
    b.__cmp__(a)
    id(a) == id(b)

If you write::

    a < b

Python tries to call::

    a.__lt__(b)
    b.__gt__(a)
    a.__cmp__(b)
    b.__cmp__(a)
    id(a) < id(b)

If you write::

    a <= b

Python tries to call::

    a.__le__(b)
    b.__ge__(a)
    a.__cmp__(b)
    b.__cmp__(a)
    id(a) <= id(b)

And similarly for ``a > b`` and ``a >= b``.

If you write::

    sorted([a, b, c])

Python calls the same chain of methods as for the ``b < a`` and ``c < b``
comparisons.

If you write any of the following::

    a in {d: 5}
    a in set([d, d, d])
    set([a, b]) == set([a, b])

Python first compares hashes, e.g.::

    a.__hash__()
    d.__hash__()

If ``hash(a) != hash(d)`` then the result of the statement ``a in {d: 5}`` is
immediately ``False`` (remember how hashes work in general). If
``hash(a) == hash(d)``) Python goes through the method resolution of the
``==`` operator as shown above.

General Notes and Caveats
=========================

In the method resolution for ``<``, ``<=``, ``==``, ``!=``, ``>=``, ``>`` and
``sorted([a, b, c])`` operators the ``__hash__()`` method is *not* called, so
in these cases it doesn't matter what it returns. The ``__hash__()`` method is
only called for sets and dictionaries.

In the official Python documentation you can read about `hashable and
non-hashable <http://docs.python.org/dev/glossary.html#term-hashable>`_ objects.
In reality, you don't have to think about it, you just follow the method
resolution described here. E.g. if you try to use lists as dictionary keys, the
list's ``__hash__()`` method will be called and it returns an exception.

In SymPy, every instance of any subclass of ``Basic`` is
immutable.  Technically this means, that its behavior through all the methods
above mustn't change once the instance is created. Especially, the hash value
mustn't change (as already stated above) or else objects will get mixed up in
dictionaries and wrong values will be returned for a given key, etc....

.. _Python script:

Script To Verify This Guide
============================

The above method resolution can be verified using the following program::

    class A(object):

        def __init__(self, a, hash):
            self.a = a
            self._hash = hash

        def __lt__(self, o):
            print "%s.__lt__(%s)" % (self.a, o.a)
            return NotImplemented

        def __le__(self, o):
            print "%s.__le__(%s)" % (self.a, o.a)
            return NotImplemented

        def __gt__(self, o):
            print "%s.__gt__(%s)" % (self.a, o.a)
            return NotImplemented

        def __ge__(self, o):
            print "%s.__ge__(%s)" % (self.a, o.a)
            return NotImplemented

        def __cmp__(self, o):
            print "%s.__cmp__(%s)" % (self.a, o.a)
            #return cmp(self._hash, o._hash)
            return NotImplemented

        def __eq__(self, o):
            print "%s.__eq__(%s)" % (self.a, o.a)
            return NotImplemented

        def __ne__(self, o):
            print "%s.__ne__(%s)" % (self.a, o.a)
            return NotImplemented

        def __hash__(self):
            print "%s.__hash__()" % (self.a)
            return self._hash

    def show(s):
        print "--- %s " % s + "-"*40
        eval(s)

    a = A("a", 1)
    b = A("b", 2)
    c = A("c", 3)
    d = A("d", 1)

    show("a == b")
    show("a != b")
    show("a < b")
    show("a <= b")
    show("a > b")
    show("a >= b")
    show("sorted([a, b, c])")
    show("{d: 5}")
    show("a in {d: 5}")
    show("set([d, d, d])")
    show("a in set([d, d, d])")
    show("set([a, b])")

    print "--- x = set([a, b]); y = set([a, b]); ---"
    x = set([a, b])
    y = set([a, b])
    print "               x == y :"
    x == y

    print "--- x = set([a, b]); y = set([b, d]); ---"
    x = set([a, b])
    y = set([b, d])
    print "               x == y :"
    x == y


and its output::

    --- a == b ----------------------------------------
    a.__eq__(b)
    b.__eq__(a)
    a.__cmp__(b)
    b.__cmp__(a)
    --- a != b ----------------------------------------
    a.__ne__(b)
    b.__ne__(a)
    a.__cmp__(b)
    b.__cmp__(a)
    --- a < b ----------------------------------------
    a.__lt__(b)
    b.__gt__(a)
    a.__cmp__(b)
    b.__cmp__(a)
    --- a <= b ----------------------------------------
    a.__le__(b)
    b.__ge__(a)
    a.__cmp__(b)
    b.__cmp__(a)
    --- a > b ----------------------------------------
    a.__gt__(b)
    b.__lt__(a)
    a.__cmp__(b)
    b.__cmp__(a)
    --- a >= b ----------------------------------------
    a.__ge__(b)
    b.__le__(a)
    a.__cmp__(b)
    b.__cmp__(a)
    --- sorted([a, b, c]) ----------------------------------------
    b.__lt__(a)
    a.__gt__(b)
    b.__cmp__(a)
    a.__cmp__(b)
    c.__lt__(b)
    b.__gt__(c)
    c.__cmp__(b)
    b.__cmp__(c)
    --- {d: 5} ----------------------------------------
    d.__hash__()
    --- a in {d: 5} ----------------------------------------
    d.__hash__()
    a.__hash__()
    d.__eq__(a)
    a.__eq__(d)
    d.__cmp__(a)
    a.__cmp__(d)
    --- set([d, d, d]) ----------------------------------------
    d.__hash__()
    d.__hash__()
    d.__hash__()
    --- a in set([d, d, d]) ----------------------------------------
    d.__hash__()
    d.__hash__()
    d.__hash__()
    a.__hash__()
    d.__eq__(a)
    a.__eq__(d)
    d.__cmp__(a)
    a.__cmp__(d)
    --- set([a, b]) ----------------------------------------
    a.__hash__()
    b.__hash__()
    --- x = set([a, b]); y = set([a, b]); ---
    a.__hash__()
    b.__hash__()
    a.__hash__()
    b.__hash__()
                   x == y :
    --- x = set([a, b]); y = set([b, d]); ---
    a.__hash__()
    b.__hash__()
    b.__hash__()
    d.__hash__()
                   x == y :
    d.__eq__(a)
    a.__eq__(d)
    d.__cmp__(a)
    a.__cmp__(d)

----------

.. [1] There is also the similar ``NotImplementedError`` exception, which one may
       be tempted to raise to obtain the same effect as returning
       ``NotImplemented``.

       But these are **not** the same, and Python will completely ignore
       ``NotImplementedError`` with respect to choosing appropriate comparison
       method, and will just propagate this exception upwards, to the caller.

       So ``return NotImplemented`` is not the same as ``raise NotImplementedError``.
