Sums, products, limits and extrapolation
----------------------------------------

The functions listed here permit approximation of infinite
sums, products, and other sequence limits.
Use :func:`mpmath.fsum` and :func:`mpmath.fprod`
for summation and multiplication of finite sequences.

Summation
..........................................

:func:`nsum`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.nsum

:func:`sumem`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.sumem

:func:`sumap`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.sumap

Products
...............................

:func:`nprod`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.nprod

Limits (``limit``)
..................

:func:`limit`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.limit

Extrapolation
..........................................

The following functions provide a direct interface to
extrapolation algorithms. :func:`nsum` and :func:`limit` essentially
work by calling the following functions with an increasing
number of terms until the extrapolated limit is accurate enough.

The following functions may be useful to call directly if the
precise number of terms needed to achieve a desired accuracy is
known in advance, or if one wishes to study the convergence
properties of the algorithms.


:func:`richardson`
^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.richardson

:func:`shanks`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.shanks

:func:`levin`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.levin

:func:`cohen_alt`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.cohen_alt

