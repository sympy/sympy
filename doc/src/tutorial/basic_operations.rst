.. _tutorial-basic:

==================
 Basic Operations
==================

Here we discuss some of the most basic operations needed for expression
manipulation in SymPy.  Some more advanced operations will be discussed later
in the :ref:`advanced expression manipulation <tutorial-manipulation>` section.

    >>> from sympy import *
    >>> x, y, z = symbols("x y z")

Expressions can be compared using a regular python syntax::

    >>> from sympy.abc import x, y
    >>> x + y == y + x
    True

    >>> x + y == y - x
    False

Sometimes, you need to have a unique symbol, for example as a temporary one in
some calculation, which is going to be substituted for something else at the
end anyway. This is achieved using ``Dummy("x")``. So, to sum it
up::

    >>> from sympy import Symbol, Dummy
    >>> Symbol("x") == Symbol("x")
    True

    >>> Dummy("x") == Dummy("x")
    False

For computation, all expressions need to be in a canonical form, this is done
during the creation of the particular instance and only inexpensive operations
are performed, necessary to put the expression in the canonical form.
Whenever you construct an expression, for example ``Add(x, x)``, the
``Add.__new__()`` is called and it determines what to return. In this case::

    >>> from sympy import Add
    >>> from sympy.abc import x
    >>> e = Add(x, x)
    >>> e
    2*x

    >>> type(e)
    <class 'sympy.core.mul.Mul'>

``e`` is actually an instance of ``Mul(2, x)``, because ``Add.__new__()``
returned ``Mul``.

Substitution
============

One of the most common things you might want to do with a mathematical
expression is substitution.  Substitution replaces all instances of something
in an expression with something else.  It is done using the ``subs`` method.
For example

    >>> expr = cos(x) + 1
    >>> expr.subs(x, y)
    cos(y) + 1

Substitution is usually done for one of two reasons:

1. Evaluating an expression at a point. For example, if our expression is
   ``cos(x) + 1`` and we want to evaluate it at the point ``x = 0``, so that
   we get ``cos(0) + 1``, which is 2.

   >>> expr.subs(x, 0)
   2

2. Replacing a subexpression with another subexpression.  There are two
   reasons we might want to do this.  The first is if we are trying to build
   an expression that has some symmetry, such as `x^{x^{x^x}}`.  To build
   this, we might start with ``x**y``, and replace ``y`` with ``x**y``.  We
   would then get ``x**(x**y)``.  If we replaced ``y`` in this new expression
   with ``x**x``, we would get ``x**(x**(x**x))``, the desired expression.

   >>> expr = x**y
   >>> expr
   x**y
   >>> expr = expr.subs(y, x**y)
   >>> expr
   x**(x**y)
   >>> expr = expr.subs(y, x**x)
   >>> expr
   x**(x**(x**x))

   The second is if we want to perform a very controlled simplification, or
   perhaps a simplification that SymPy is otherwise unable to do.  For
   example, say we have `\sin(2x) + \cos(2x)`, and we want to replace
   `\sin(2x)` with `2\sin(x)\cos(x)`.  As we will learn later, the function
   ``expand_trig`` does this.  However, this function will also expand
   `\cos(2x)`, which we may not want.  While there are ways to perform such
   precise simplification, and we will learn some of them in the
   :ref:`advanced expression manipulation <tutorial-manipulation>` section, an
   easy way is to just replace `\sin(2x)` with `2\sin(x)\cos(x)`.

   >>> expr = sin(2*x) + cos(2*x)
   >>> expand_trig(expr)
   2*sin(x)*cos(x) + 2*cos(x)**2 - 1
   >>> expr.subs(sin(2*x), 2*sin(x)*cos(x))
   2*sin(x)*cos(x) + cos(2*x)

There are two important things to note about ``subs``.  First, it returns a
new expression.  SymPy objects are immutable.  That means that ``subs`` does
not modify it in-place.  For example

   >>> expr = cos(x)
   >>> expr.subs(x, 0)
   1
   >>> expr
   cos(x)
   >>> x
   x

.. sidebar:: Quick Tip

   SymPy expressions are immutable.  No function will change them in-place.

Here, we see that performing ``expr.subs(x, 0)`` leaves ``expr`` unchanged.
In fact, since SymPy expressions are immutable, no function will change them
in-place.  All functions will return new expressions.

To perform multiple substitutions at once, pass a list of ``(old, new)`` pairs
to ``subs``.

    >>> expr = x**3 + 4*x*y - z
    >>> expr.subs([(x, 2), (y, 4), (z, 0)])
    40

It is often useful to combine this with a list comprehension to do a large set
of similar replacements all at once.  For example, say we had `x^4 - 4x^3 + 4x^2 -
2x + 3` and we wanted to replace all instances of `x` that have an even power
with `y`, to get `y^4 - 4x^3 + 4y^2 - 2x + 3`.

    >>> expr = x**4 - 4*x**3 + 4*x**2 - 2*x + 3
    >>> replacements = [(x**i, y**i) for i in range(5) if i % 2 == 0]
    >>> expr.subs(replacements)
    -4*x**3 - 2*x + y**4 + 4*y**2 + 3

Converting Strings to SymPy Expressions
=======================================

The ``sympify`` function (that's ``sympify``, not to be confused with
``simplify``) can be used to convert strings into SymPy expressions.

For example

    >>> str_expr = "x**2 + 3*x - 1/2"
    >>> expr = sympify(str_expr)
    >>> expr
    x**2 + 3*x - 1/2
    >>> expr.subs(x, 2)
    19/2

.. warning:: ``sympify`` uses ``eval``.  Don't use it on unsanitized input.

``evalf``
=========

To evaluate a numerical expression into a floating point number, use
``evalf``.

    >>> expr = sqrt(8)
    >>> expr.evalf()
    2.82842712474619

SymPy can evaluate floating point expressions to arbitrary precision.  By
default, 15 digits of precision are used, but you can pass any number as the
argument to ``evalf``.  Let's compute the first 100 digits of `\pi`.

    >>> pi.evalf(100)
    3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068

To numerically evaluate an expression with a Symbol at a point, we might use
``subs`` followed by ``evalf``, but it is more efficient and numerically
stable to pass the substitution to ``evalf`` using the ``subs`` flag, which
takes a dictionary of ``Symbol: point`` pairs.

    >>> expr = cos(2*x)
    >>> expr.evalf(subs={x: 2.4})
    0.0874989834394464

Sometimes there are roundoff errors smaller than the desired precision that
remain after an expression is evaluated. Such numbers can be removed at the
user's discretion by setting the ``chop`` flag to True.

    >>> one = cos(1)**2 + sin(1)**2
    >>> (one - 1).evalf()
    -0.e-124
    >>> (one - 1).evalf(chop=True)
    0

``lambdify``
============

``subs`` and ``evalf`` are good if you want to do simple evaluation, but if
you intend to evaluate an expression at many points, there are more efficient
ways.  For example, if you wanted to evaluate an expression at a thousand
points, using SymPy would be far slower than it needs to be, especially if you
only care about machine precision.  Instead, you should use libraries like
`NumPy <http://www.numpy.org/>`_ and `SciPy <http://www.scipy.org/>`_.

The easiest way to convert a SymPy expression to an expression that can be
numerically evaluated is to use the ``lambdify`` function.  ``lambdify`` acts
like a ``lambda`` function, except it converts the SymPy names to the names of
the given numerical library, usually NumPy.  For example

    >>> import numpy # doctest:+SKIP
    >>> a = numpy.arange(10) # doctest:+SKIP
    >>> expr = sin(x)
    >>> f = lambdify(x, expr, "numpy") # doctest:+SKIP
    >>> f(a) # doctest:+SKIP
    [ 0.          0.84147098  0.90929743  0.14112001 -0.7568025  -0.95892427
     -0.2794155   0.6569866   0.98935825  0.41211849]

You can use other libraries than NumPy. For example, to use the standard
library math module, use ``"math"``.

    >>> f = lambdify(x, expr, "math")
    >>> f(0.1)
    0.0998334166468

To use lambdify with numerical libraries that it does not know about, pass a
dictionary of ``sympy_name:numerical_function`` pairs.  For example

    >>> def mysin(x):
    ...     """
    ...     My sine. Note that this is only accurate for small x.
    ...     """
    ...     return x
    >>> f = lambdify(x, expr, {"sin":mysin})
    >>> f(0.1)
    0.1


Functions
=========

How to create a new function with one variable::

    class sign(Function):

        nargs = 1

        @classmethod
        def eval(cls, arg):
            if isinstance(arg, Basic.NaN):
                return S.NaN
            if isinstance(arg, Basic.Zero):
                return S.Zero
            if arg.is_positive:
                return S.One
            if arg.is_negative:
                return S.NegativeOne
            if isinstance(arg, Basic.Mul):
                coeff, terms = arg.as_coeff_mul()
                if not isinstance(coeff, Basic.One):
                    return cls(coeff) * cls(Basic.Mul(*terms))

        is_finite = True

        def _eval_conjugate(self):
            return self

        def _eval_is_zero(self):
            return isinstance(self[0], Basic.Zero)

and that's it. The ``_eval_*`` functions are called when something is needed.
The ``eval`` is called when the class is about to be instantiated and it
should return either some simplified instance of some other class or if the
class should be unmodified, return ``None`` (see ``core/function.py`` in
``Function.__new__`` for implementation details). See also tests in
`sympy/functions/elementary/tests/test_interface.py
<https://github.com/sympy/sympy/blob/master/sympy/functions/elementary/tests/test_interface.py>`_ that test this interface. You can use them to create your own new functions.

The applied function ``sign(x)`` is constructed using
::

    sign(x)

both inside and outside of SymPy. Unapplied functions ``sign`` is just the class
itself::

    sign

both inside and outside of SymPy. This is the current structure of classes in
SymPy::

    class BasicType(type):
        pass
    class MetaBasicMeths(BasicType):
        ...
    class BasicMeths(AssumeMeths):
        __metaclass__ = MetaBasicMeths
        ...
    class Basic(BasicMeths):
        ...
    class FunctionClass(MetaBasicMeths):
        ...
    class Function(Basic, RelMeths, ArithMeths):
        __metaclass__ = FunctionClass
        ...

The exact names of the classes and the names of the methods and how they work
can be changed in the future.

This is how to create a function with two variables::

    class chebyshevt_root(Function):
        nargs = 2

        @classmethod
        def eval(cls, n, k):
            if not 0 <= k < n:
                raise ValueError("must have 0 <= k < n")
            return cos(S.Pi*(2*k + 1)/(2*n))


.. note:: the first argument of a @classmethod should be ``cls`` (i.e. not
          ``self``).

Here it's how to define a derivative of the function::

    >>> from sympy import Function, sympify, cos
    >>> class my_function(Function):
    ...     nargs = 1
    ...
    ...     def fdiff(self, argindex = 1):
    ...         return cos(self.args[0])
    ...
    ...     @classmethod
    ...     def eval(cls, arg):
    ...         arg = sympify(arg)
    ...         if arg == 0:
    ...             return sympify(0)

So guess what this ``my_function`` is going to be? Well, it's derivative is
``cos`` and the function value at 0 is 0, but let's pretend we don't know::

    >>> from sympy import pprint
    >>> pprint(my_function(x).series(x, 0, 10))
         3     5     7       9
        x     x     x       x       / 10\
    x - -- + --- - ---- + ------ + O\x  /
        6    120   5040   362880

Looks familiar indeed::

    >>> from sympy import sin
    >>> pprint(sin(x).series(x, 0, 10))
         3     5     7       9
        x     x     x       x       / 10\
    x - -- + --- - ---- + ------ + O\x  /
        6    120   5040   362880

Let's try a more complicated example. Let's define the derivative in terms of
the function itself::

    >>> class what_am_i(Function):
    ...     nargs = 1
    ...
    ...     def fdiff(self, argindex = 1):
    ...         return 1 - what_am_i(self.args[0])**2
    ...
    ...     @classmethod
    ...     def eval(cls, arg):
    ...         arg = sympify(arg)
    ...         if arg == 0:
    ...             return sympify(0)

So what is ``what_am_i``?  Let's try it::

    >>> pprint(what_am_i(x).series(x, 0, 10))
         3      5       7       9
        x    2*x    17*x    62*x     / 10\
    x - -- + ---- - ----- + ----- + O\x  /
        3     15     315     2835

Well, it's ``tanh``::

    >>> from sympy import tanh
    >>> pprint(tanh(x).series(x, 0, 10))
         3      5       7       9
        x    2*x    17*x    62*x     / 10\
    x - -- + ---- - ----- + ----- + O\x  /
        3     15     315     2835

The new functions we just defined are regular SymPy objects, you
can use them all over SymPy, e.g.::

    >>> from sympy import limit
    >>> limit(what_am_i(x)/x, x, 0)
    1


Common tasks
------------

Please use the same way as is shown below all across SymPy.

**accessing parameters**::

    >>> from sympy import sign, sin
    >>> from sympy.abc import x, y, z

    >>> e = sign(x**2)
    >>> e.args
    (x**2,)

    >>> e.args[0]
    x**2

    Number arguments (in Adds and Muls) will always be the first argument;
    other arguments might be in arbitrary order:
    >>> (1 + x + y*z).args[0]
    1
    >>> (1 + x + y*z).args[1] in (x, y*z)
    True

    >>> (y*z).args
    (y, z)

    >>> sin(y*z).args
    (y*z,)

Never use internal methods or variables, prefixed with "``_``" (example: don't
use ``_args``, use ``.args`` instead).

**testing the structure of a SymPy expression**

Applied functions::

    >>> from sympy import sign, exp, Function
    >>> e = sign(x**2)

    >>> isinstance(e, sign)
    True

    >>> isinstance(e, exp)
    False

    >>> isinstance(e, Function)
    True

So ``e`` is a ``sign(z)`` function, but not ``exp(z)`` function.

Unapplied functions::

    >>> from sympy import sign, exp, FunctionClass
    >>> e = sign

    >>> f = exp

    >>> g = Add

    >>> isinstance(e, FunctionClass)
    True

    >>> isinstance(f, FunctionClass)
    True

    >>> isinstance(g, FunctionClass)
    False

    >>> g is Add
    True

So ``e`` and ``f`` are functions, ``g`` is not a function.

.. TODO: Write an advanced numerics section
