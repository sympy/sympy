Orthogonal polynomials
----------------------

An orthogonal polynomial sequence is a sequence of polynomials `P_0(x), P_1(x), \ldots` of degree `0, 1, \ldots`, which are mutually orthogonal in the sense that 

.. math ::

    \int_S P_n(x) P_m(x) w(x) dx = 
    \begin{cases}
    c_n \ne 0 & \text{if $m = n$} \\
    0         & \text{if $m \ne n$}
    \end{cases}

where `S` is some domain (e.g. an interval `[a,b] \in \mathbb{R}`) and `w(x)` is a fixed *weight function*. A sequence of orthogonal polynomials is determined completely by `w`, `S`, and a normalization convention (e.g. `c_n = 1`). Applications of orthogonal polynomials include function approximation and solution of differential equations.

Orthogonal polynomials are sometimes defined using the differential equations they satisfy (as functions of `x`) or the recurrence relations they satisfy with respect to the order `n`. Other ways of defining orthogonal polynomials include differentiation formulas and generating functions. The standard orthogonal polynomials can also be represented as hypergeometric series (see :doc:`hypergeometric`), more specifically using the Gauss hypergeometric function `\,_2F_1` in most cases. The following functions are generally implemented using hypergeometric functions since this is computationally efficient and easily generalizes.

For more information, see the `Wikipedia article on orthogonal polynomials <http://en.wikipedia.org/wiki/Orthogonal_polynomials>`_.

Legendre functions
.......................................

:func:`legendre`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.legendre(n, x)

:func:`legenp`
^^^^^^^^^^^^^^^
.. autofunction:: mpmath.legenp(n, m, z, type=2)

:func:`legenq`
^^^^^^^^^^^^^^^
.. autofunction:: mpmath.legenq(n, m, z, type=2)

Chebyshev polynomials
.....................

:func:`chebyt`
^^^^^^^^^^^^^^^
.. autofunction:: mpmath.chebyt(n, x)

:func:`chebyu`
^^^^^^^^^^^^^^^
.. autofunction:: mpmath.chebyu(n, x)

Jacobi polynomials
..................

:func:`jacobi`
^^^^^^^^^^^^^^
.. autofunction:: mpmath.jacobi(n, a, b, z)

Gegenbauer polynomials
.....................................

:func:`gegenbauer`
^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.gegenbauer(n, a, z)

Hermite polynomials
.....................................

:func:`hermite`
^^^^^^^^^^^^^^^
.. autofunction:: mpmath.hermite(n, z)

Laguerre polynomials
.......................................

:func:`laguerre`
^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.laguerre(n, a, z)

Spherical harmonics
.....................................

:func:`spherharm`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.spherharm(l, m, theta, phi)
