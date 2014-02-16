Factorials and gamma functions
------------------------------

Factorials and factorial-like sums and products are basic tools of combinatorics and number theory. Much like the exponential function is fundamental to differential equations and analysis in general, the factorial function (and its extension to complex numbers, the gamma function) is fundamental to difference equations and functional equations.

A large selection of factorial-like functions is implemented in mpmath. All functions support complex arguments, and arguments may be arbitrarily large. Results are numerical approximations, so to compute *exact* values a high enough precision must be set manually::

    >>> mp.dps = 15; mp.pretty = True
    >>> fac(100)
    9.33262154439442e+157
    >>> print int(_)    # most digits are wrong
    93326215443944150965646704795953882578400970373184098831012889540582227238570431
    295066113089288327277825849664006524270554535976289719382852181865895959724032
    >>> mp.dps = 160
    >>> fac(100)
    93326215443944152681699238856266700490715968264381621468592963895217599993229915
    608941463976156518286253697920827223758251185210916864000000000000000000000000.0

The gamma and polygamma functions are closely related to :doc:`zeta`. See also :doc:`qfunctions` for q-analogs of factorial-like functions.


Factorials
..........

:func:`factorial`/:func:`fac`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: mpmath.factorial(x, **kwargs)

:func:`fac2`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.fac2(x)

Binomial coefficients 
....................................................

:func:`binomial`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.binomial(n,k)


Gamma function
..............

:func:`gamma`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.gamma(x, **kwargs)

:func:`rgamma`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.rgamma(x, **kwargs)

:func:`gammaprod`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.gammaprod(a, b)

:func:`loggamma`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.loggamma(x)


Rising and falling factorials
.............................

:func:`rf`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.rf(x,n)

:func:`ff`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.ff(x,n)

Beta function
.............

:func:`beta`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.beta(x,y)

:func:`betainc`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.betainc(a,b,x1=0,x2=1,regularized=False)


Super- and hyperfactorials
..........................

:func:`superfac`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.superfac(z)

:func:`hyperfac`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.hyperfac(z)

:func:`barnesg`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.barnesg(z)


Polygamma functions and harmonic numbers
........................................

:func:`psi`/:func:`digamma`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.psi(m, z)

.. autofunction:: mpmath.digamma(z)

:func:`harmonic`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.harmonic(z)
