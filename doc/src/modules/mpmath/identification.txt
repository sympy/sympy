Number identification
=====================

Most function in mpmath are concerned with producing approximations from exact mathematical formulas. It is also useful to consider the inverse problem: given only a decimal approximation for a number, such as 0.7320508075688772935274463, is it possible to find an exact formula?

Subject to certain restrictions, such "reverse engineering" is indeed possible thanks to the existence of *integer relation algorithms*. Mpmath implements the PSLQ algorithm (developed by H. Ferguson), which is one such algorithm.

Automated number recognition based on PSLQ is not a silver bullet. Any occurring transcendental constants (`\pi`, `e`, etc) must be guessed by the user, and the relation between those constants in the formula must be linear (such as `x = 3 \pi + 4 e`). More complex formulas can be found by combining PSLQ with functional transformations; however, this is only feasible to a limited extent since the computation time grows exponentially with the number of operations that need to be combined.

The number identification facilities in mpmath are inspired by the `Inverse Symbolic Calculator <http://oldweb.cecm.sfu.ca/projects/ISC/ISCmain.html>`_ (ISC). The ISC is more powerful than mpmath, as it uses a lookup table of millions of precomputed constants (thereby mitigating the problem with exponential complexity).

Constant recognition
-----------------------------------

:func:`identify`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.identify

Algebraic identification
---------------------------------------

:func:`findpoly`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.findpoly

Integer relations (PSLQ)
----------------------------

:func:`pslq`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.pslq
