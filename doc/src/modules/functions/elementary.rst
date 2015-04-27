Elementary
==========

This module implements elementary functions, as well as functions like ``Abs``,
``Max``, etc.


Abs
---

Returns the absolute value of the argument.

Examples::

    >>> from sympy.functions import Abs
    >>> Abs(-1)
    1

.. autoclass:: sympy.functions.elementary.complexes.Abs
   :members:

acos
----

.. autoclass:: sympy.functions.elementary.trigonometric.acos
   :members:

acosh
-----

.. autoclass:: sympy.functions.elementary.hyperbolic.acosh
   :members:

acot
----

.. autoclass:: sympy.functions.elementary.trigonometric.acot
   :members:

acoth
-----

.. autoclass:: sympy.functions.elementary.hyperbolic.acoth
   :members:

arg
---

Returns the argument (in radians) of a complex number. For a real
number, the argument is always 0.

Examples::

    >>> from sympy.functions import arg
    >>> from sympy import I, sqrt
    >>> arg(2.0)
    0
    >>> arg(I)
    pi/2
    >>> arg(sqrt(2) + I*sqrt(2))
    pi/4

.. autoclass:: sympy.functions.elementary.complexes.arg
   :members:

asin
----

.. autoclass:: sympy.functions.elementary.trigonometric.asin
   :members:

asinh
-----

.. autoclass:: sympy.functions.elementary.hyperbolic.asinh
   :members:

atan
----

.. autoclass:: sympy.functions.elementary.trigonometric.atan
   :members:

atan2
-----

This function is like `\operatorname{atan}`, but considers the sign of both
arguments in order to correctly determine the quadrant of its result.

.. autoclass:: sympy.functions.elementary.trigonometric.atan2
   :members:

atanh
-----

.. autoclass:: sympy.functions.elementary.hyperbolic.atanh
   :members:

ceiling
-------

.. autoclass:: sympy.functions.elementary.integers.ceiling
   :members:

conjugate
---------

Returns the `complex conjugate <http://en.wikipedia.org/wiki/Complex_conjugation>`_
of an argument. In mathematics, the complex conjugate of a complex number is given
by changing the sign of the imaginary part. Thus, the conjugate of the complex number

    :math:`a + ib`

(where a and b are real numbers) is

    :math:`a - ib`

Examples::

    >>> from sympy.functions import conjugate
    >>> from sympy import I
    >>> conjugate(2)
    2
    >>> conjugate(I)
    -I

.. autoclass:: sympy.functions.elementary.complexes.conjugate
   :members:

cos
---

Returns the trigonmetric cosine value of the argument. Given that the argument is an angle value in radians, the cosine of the argument is the ratio of the length of the adjacent side to the length of the hypotenuse of a triangle.

Examples::

    >>> from sympy.functions import cos
    >>> cos(pi)
    -1
    >>> cos(pi/3)
    1/2
    >>> cos(-5*pi/3)
    -sqrt(3)/2
    
.. autoclass:: sympy.functions.elementary.trigonometric.cos
   :members:

cosh
----
.. autoclass:: sympy.functions.elementary.hyperbolic.cosh
   :members:

cot
---

.. autoclass:: sympy.functions.elementary.trigonometric.cot
   :members:

coth
----

.. autoclass:: sympy.functions.elementary.hyperbolic.coth
   :members:

exp
---

.. autoclass:: sympy.functions.elementary.exponential.exp
   :members:

.. seealso:: classes :py:class:`sympy.functions.elementary.exponential.log`

ExprCondPair
------------

.. autoclass:: sympy.functions.elementary.piecewise.ExprCondPair
   :members:

floor
-----

.. autoclass:: sympy.functions.elementary.integers.floor
   :members:

HyperbolicFunction
------------------

.. autoclass:: sympy.functions.elementary.hyperbolic.HyperbolicFunction
   :members:

IdentityFunction
----------------

.. autoclass:: sympy.functions.elementary.miscellaneous.IdentityFunction
   :members:

im
--

Returns the imaginary part of an expression.

Examples::


    >>> from sympy.functions import im
    >>> from sympy import I
    >>> im(2+3*I)
    3

.. autoclass: sympy.functions.elementary.im
   :members:

.. seealso::
   :py:class:`sympy.functions.elementary.complexes.re`

LambertW
--------

.. autoclass:: sympy.functions.elementary.exponential.LambertW
   :members:

log
---

.. autoclass:: sympy.functions.elementary.exponential.log
   :members:

.. seealso:: classes :py:class:`sympy.functions.elementary.exponential.exp`

Min
---

Returns the minimum of two (comparable) expressions.

Examples::

    >>> from sympy.functions import Min
    >>> Min(1,2)
    1
    >>> from sympy.abc import x
    >>> Min(1, x)
    Min(1, x)

It is named ``Min`` and not ``min`` to avoid conflicts with the built-in function ``min``.

.. autoclass:: sympy.functions.elementary.miscellaneous.Min
   :members:


Max
---

Returns the maximum of two (comparable) expressions

It is named ``Max`` and not ``max`` to avoid conflicts with the built-in function ``max``.

.. autoclass:: sympy.functions.elementary.miscellaneous.Max
   :members:

Piecewise
---------

.. autoclass:: sympy.functions.elementary.piecewise.Piecewise
   :members:

.. autofunction:: sympy.functions.elementary.piecewise.piecewise_fold

re
--

Return the real part of an expression.

Examples::

    >>> from sympy.functions import re
    >>> from sympy import I
    >>> re(2+3*I)
    2

.. autoclass:: sympy.functions.elementary.complexes.re
   :members:

.. seealso::
   :py:class:`sympy.functions.elementary.complexes.im`

root
----

.. autofunction:: sympy.functions.elementary.miscellaneous.root

RoundFunction
-------------

.. autoclass:: sympy.functions.elementary.integers.RoundFunction

sin
---

.. autoclass:: sympy.functions.elementary.trigonometric.sin
   :members:

sinh
----

.. autoclass:: sympy.functions.elementary.hyperbolic.sinh
   :members:

sqrt
----

Returns the square root of an expression. It is equivalent to raise to ``Rational(1,2)``.

    >>> from sympy.functions import sqrt
    >>> from sympy import Rational
    >>> sqrt(2) == 2**Rational(1,2)
    True

.. autofunction:: sympy.functions.elementary.miscellaneous.sqrt

sign
----

.. autoclass:: sympy.functions.elementary.complexes.sign
   :members:

tan
---

.. autoclass:: sympy.functions.elementary.trigonometric.tan
   :members:

tanh
----

.. autoclass:: sympy.functions.elementary.hyperbolic.tanh
   :members:
