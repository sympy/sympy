Exponential integrals and error functions
-----------------------------------------

Exponential integrals give closed-form solutions to a large class of commonly occurring transcendental integrals that cannot be evaluated using elementary functions. Integrals of this type include those with an integrand of the form `t^a e^{t}` or `e^{-x^2}`, the latter giving rise to the Gaussian (or normal) probability distribution.

The most general function in this section is the incomplete gamma function, to which all others can be reduced. The incomplete gamma function, in turn, can be expressed using hypergeometric functions (see :doc:`hypergeometric`).

Incomplete gamma functions
..........................

:func:`gammainc`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.gammainc(z, a=0, b=inf, regularized=False)


Exponential integrals
.....................

:func:`ei`
^^^^^^^^^^
.. autofunction:: mpmath.ei(x, **kwargs)

:func:`e1`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.e1(x, **kwargs)

:func:`expint`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.expint(*args)


Logarithmic integral
....................

:func:`li`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.li(x, **kwargs)


Trigonometric integrals
.......................

:func:`ci`
^^^^^^^^^^
.. autofunction:: mpmath.ci(x, **kwargs)

:func:`si`
^^^^^^^^^^
.. autofunction:: mpmath.si(x, **kwargs)


Hyperbolic integrals
....................

:func:`chi`
^^^^^^^^^^^
.. autofunction:: mpmath.chi(x, **kwargs)

:func:`shi`
^^^^^^^^^^^
.. autofunction:: mpmath.shi(x, **kwargs)


Error functions
...............

:func:`erf`
^^^^^^^^^^^
.. autofunction:: mpmath.erf(x, **kwargs)

:func:`erfc`
^^^^^^^^^^^^
.. autofunction:: mpmath.erfc(x, **kwargs)

:func:`erfi`
^^^^^^^^^^^^
.. autofunction:: mpmath.erfi(x)

:func:`erfinv`
^^^^^^^^^^^^^^
.. autofunction:: mpmath.erfinv(x)

The normal distribution
....................................................

:func:`npdf`
^^^^^^^^^^^^
.. autofunction:: mpmath.npdf(x, mu=0, sigma=1)

:func:`ncdf`
^^^^^^^^^^^^
.. autofunction:: mpmath.ncdf(x, mu=0, sigma=1)


Fresnel integrals
......................................................

:func:`fresnels`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.fresnels(x)

:func:`fresnelc`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.fresnelc(x)
