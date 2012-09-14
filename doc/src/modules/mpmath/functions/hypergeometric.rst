Hypergeometric functions
------------------------

The functions listed in :doc:`expintegrals`, :doc:`bessel` and
:doc:`orthogonal`, and many other functions as well, are merely
particular instances of the generalized hypergeometric function `\,_pF_q`.
The functions listed in the following section enable efficient
direct evaluation of the underlying hypergeometric series, as
well as linear combinations, limits with respect to parameters,
and analytic continuations thereof. Extensions to twodimensional
series are also provided. See also the basic or q-analog of
the hypergeometric series in :doc:`qfunctions`.

For convenience, most of the hypergeometric series of low order are
provided as standalone functions. They can equivalently be evaluated using
:func:`~mpmath.hyper`. As will be demonstrated in the respective docstrings,
all the ``hyp#f#`` functions implement analytic continuations and/or asymptotic
expansions with respect to the argument `z`, thereby permitting evaluation
for `z` anywhere in the complex plane. Functions of higher degree can be
computed via :func:`~mpmath.hyper`, but generally only in rapidly convergent
instances.

Most hypergeometric and hypergeometric-derived functions accept optional
keyword arguments to specify options for :func:`hypercomb` or
:func:`hyper`. Some useful options are *maxprec*, *maxterms*,
*zeroprec*, *accurate_small*, *hmag*, *force_series*,
*asymp_tol* and *eliminate*. These options give control over what to
do in case of slow convergence, extreme loss of accuracy or
evaluation at zeros (these two cases cannot generally be
distinguished from each other automatically),
and singular parameter combinations.

Common hypergeometric series
............................

:func:`hyp0f1`
^^^^^^^^^^^^^^
.. autofunction:: mpmath.hyp0f1(a, z)

:func:`hyp1f1`
^^^^^^^^^^^^^^
.. autofunction:: mpmath.hyp1f1(a, b, z)

:func:`hyp1f2`
^^^^^^^^^^^^^^
.. autofunction:: mpmath.hyp1f2(a1, b1, b2, z)

:func:`hyp2f0`
^^^^^^^^^^^^^^
.. autofunction:: mpmath.hyp2f0(a, b, z)

:func:`hyp2f1`
^^^^^^^^^^^^^^
.. autofunction:: mpmath.hyp2f1(a, b, c, z)

:func:`hyp2f2`
^^^^^^^^^^^^^^
.. autofunction:: mpmath.hyp2f2(a1, a2, b1, b2, z)

:func:`hyp2f3`
^^^^^^^^^^^^^^
.. autofunction:: mpmath.hyp2f3(a1, a2, b1, b2, b3, z)

:func:`hyp3f2`
^^^^^^^^^^^^^^
.. autofunction:: mpmath.hyp3f2(a1, a2, a3, b1, b2, z)

Generalized hypergeometric functions
....................................

:func:`hyper`
^^^^^^^^^^^^^^
.. autofunction:: mpmath.hyper(a_s, b_s, z)

:func:`hypercomb`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.hypercomb

Meijer G-function
...................................

:func:`meijerg`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.meijerg(a_s,b_s,z,r=1,**kwargs)

Bilateral hypergeometric series
...............................

:func:`bihyper`
^^^^^^^^^^^^^^^
.. autofunction:: mpmath.bihyper(a_s,b_s,z,**kwargs)

Hypergeometric functions of two variables
...............................................

:func:`hyper2d`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.hyper2d(a,b,x,y,**kwargs)

:func:`appellf1`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.appellf1(a,b1,b2,c,x,y,**kwargs)

:func:`appellf2`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.appellf2(a,b1,b2,c1,c2,x,y,**kwargs)

:func:`appellf3`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.appellf3(a1,a2,b1,b2,c,x,y,**kwargs)

:func:`appellf4`
^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.appellf4(a,b,c1,c2,x,y,**kwargs)

