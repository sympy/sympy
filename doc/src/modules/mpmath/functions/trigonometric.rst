Trigonometric functions
-----------------------

Except where otherwise noted, the trigonometric functions
take a radian angle as input and the inverse trigonometric
functions return radian angles.

The ordinary trigonometric functions are single-valued
functions defined everywhere in the complex plane
(except at the poles of tan, sec, csc, and cot).
They are defined generally via the exponential function,
e.g.

.. math ::

    \cos(x) = \frac{e^{ix} + e^{-ix}}{2}.

The inverse trigonometric functions are multivalued,
thus requiring branch cuts, and are generally real-valued
only on a part of the real line. Definitions and branch cuts
are given in the documentation of each function.
The branch cut conventions used by mpmath are essentially
the same as those found in most standard mathematical software,
such as Mathematica and Python's own ``cmath`` libary (as of Python 2.6;
earlier Python versions implement some functions
erroneously).

Degree-radian conversion
...........................................................

:func:`degrees`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.degrees(x)

:func:`radians`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.radians(x)

Trigonometric functions
.......................

:func:`cos`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.cos(x, **kwargs)

:func:`sin`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.sin(x, **kwargs)

:func:`tan`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.tan(x, **kwargs)

:func:`sec`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.sec(x)

:func:`csc`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.csc(x)

:func:`cot`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.cot(x)

Trigonometric functions with modified argument
........................................................

:func:`cospi`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.cospi(x, **kwargs)

:func:`sinpi`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.sinpi(x, **kwargs)

Inverse trigonometric functions
................................................

:func:`acos`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.acos(x, **kwargs)

:func:`asin`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.asin(x, **kwargs)

:func:`atan`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.atan(x, **kwargs)

:func:`atan2`
^^^^^^^^^^^^^
.. autofunction:: mpmath.atan2(y, x)

:func:`asec`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.asec(x)

:func:`acsc`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.acsc(x)

:func:`acot`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.acot(x)

Sinc function
.............

:func:`sinc`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.sinc(x)

:func:`sincpi`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.sincpi(x)
