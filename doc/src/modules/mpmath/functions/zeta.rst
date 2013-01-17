Zeta functions, L-series and polylogarithms
-------------------------------------------

This section includes the Riemann zeta functions
and associated functions pertaining to analytic number theory.


Riemann and Hurwitz zeta functions
..................................................

:func:`zeta`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.zeta(s,a=1,derivative=0)


Dirichlet L-series
..................................................

:func:`altzeta`
^^^^^^^^^^^^^^^
.. autofunction:: mpmath.altzeta(s)

:func:`dirichlet`
^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.dirichlet(s,chi,derivative=0)


Stieltjes constants
...................

:func:`stieltjes`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.stieltjes(n,a=1)


Zeta function zeros
......................................

These functions are used for the study of the Riemann zeta function
in the critical strip.

:func:`zetazero`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.zetazero(n, verbose=False)

:func:`nzeros`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.nzeros(t)

:func:`siegelz`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.siegelz(t)

:func:`siegeltheta`
^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.siegeltheta(t)

:func:`grampoint`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.grampoint(n)

:func:`backlunds`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.backlunds(t)


Lerch transcendent
................................

:func:`lerchphi`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.lerchphi(z,s,a)


Polylogarithms and Clausen functions
.......................................

:func:`polylog`
^^^^^^^^^^^^^^^
.. autofunction:: mpmath.polylog(s,z)

:func:`clsin`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.clsin(s, z)

:func:`clcos`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.clcos(s, z)

:func:`polyexp`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.polyexp(s,z)


Zeta function variants
..........................

:func:`primezeta`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.primezeta(s)

:func:`secondzeta`
^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.secondzeta(s, a=0.015, **kwargs)
