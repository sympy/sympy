Utility functions
===============================================

This page lists functions that perform basic operations
on numbers or aid general programming.

Conversion and printing
-----------------------

:func:`mpmathify` / :func:`convert`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.mpmathify(x, strings=True)

:func:`nstr`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.nstr(x, n=6, **kwargs)

:func:`nprint`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.nprint(x, n=6, **kwargs)

Arithmetic operations
---------------------

See also :func:`mpmath.sqrt`, :func:`mpmath.exp` etc., listed
in :doc:`functions/powers`

:func:`fadd`
^^^^^^^^^^^^^
.. autofunction:: mpmath.fadd

:func:`fsub`
^^^^^^^^^^^^^
.. autofunction:: mpmath.fsub

:func:`fneg`
^^^^^^^^^^^^^
.. autofunction:: mpmath.fneg

:func:`fmul`
^^^^^^^^^^^^^
.. autofunction:: mpmath.fmul

:func:`fdiv`
^^^^^^^^^^^^^
.. autofunction:: mpmath.fdiv

:func:`fmod`
^^^^^^^^^^^^^
.. autofunction:: mpmath.fmod(x, y)

:func:`fsum`
^^^^^^^^^^^^^
.. autofunction:: mpmath.fsum(terms, absolute=False, squared=False)

:func:`fprod`
^^^^^^^^^^^^^
.. autofunction:: mpmath.fprod(factors)

:func:`fdot`
^^^^^^^^^^^^^
.. autofunction:: mpmath.fdot(A, B=None, conjugate=False)

Complex components
------------------

:func:`fabs`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.fabs(x)

:func:`sign`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.sign(x)

:func:`re`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.re(x)

:func:`im`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.im(x)

:func:`arg`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.arg(x)

:func:`conj`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.conj(x)

:func:`polar`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.polar(x)

:func:`rect`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.rect(x)

Integer and fractional parts
-----------------------------

:func:`floor`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.floor(x)

:func:`ceil`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.ceil(x)

:func:`nint`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.nint(x)

:func:`frac`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.frac(x)

Tolerances and approximate comparisons
--------------------------------------

:func:`chop`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.chop(x, tol=None)

:func:`almosteq`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.almosteq(s, t, rel_eps=None, abs_eps=None)

Properties of numbers
-------------------------------------

:func:`isinf`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.isinf(x)

:func:`isnan`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.isnan(x)

:func:`isnormal`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.isnormal(x)

:func:`isfinite`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.isfinite(x)

:func:`isint`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.isint(x, gaussian=False)

:func:`ldexp`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.ldexp(x, n)

:func:`frexp`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.frexp(x, n)

:func:`mag`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.mag(x)

:func:`nint_distance`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.nint_distance(x)

.. :func:`absmin`
.. ^^^^^^^^^^^^^^^^^^^^
.. .. autofunction:: mpmath.absmin(x)
.. .. autofunction:: mpmath.absmax(x)

Number generation
-----------------

:func:`fraction`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.fraction(p,q)

:func:`rand`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.rand()

:func:`arange`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.arange(*args)

:func:`linspace`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.linspace(*args, **kwargs)

Precision management
--------------------

:func:`autoprec`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.autoprec

:func:`workprec`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.workprec

:func:`workdps`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.workdps

:func:`extraprec`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.extraprec

:func:`extradps`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.extradps

Performance and debugging
------------------------------------

:func:`memoize`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.memoize

:func:`maxcalls`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.maxcalls

:func:`monitor`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.monitor

:func:`timing`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.timing
