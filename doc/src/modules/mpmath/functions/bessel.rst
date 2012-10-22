Bessel functions and related functions
--------------------------------------

The functions in this section arise as solutions to various differential equations in physics, typically describing wavelike oscillatory behavior or a combination of oscillation and exponential decay or growth. Mathematically, they are special cases of the confluent hypergeometric functions `\,_0F_1`, `\,_1F_1` and `\,_1F_2` (see :doc:`hypergeometric`).


Bessel functions
...................................................

:func:`besselj`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: mpmath.besselj(n,x,derivative=0)
.. autofunction:: mpmath.j0(x)
.. autofunction:: mpmath.j1(x)

:func:`bessely`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.bessely(n,x,derivative=0)

:func:`besseli`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.besseli(n,x,derivative=0)

:func:`besselk`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.besselk(n,x)


Bessel function zeros
...............................

:func:`besseljzero`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.besseljzero(v,m,derivative=0)

:func:`besselyzero`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.besselyzero(v,m,derivative=0)


Hankel functions
................

:func:`hankel1`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.hankel1(n,x)

:func:`hankel2`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.hankel2(n,x)


Kelvin functions
................

:func:`ber`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.ber

:func:`bei`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: mpmath.bei

:func:`ker`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.ker

:func:`kei`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.kei


Struve functions
...................................................

:func:`struveh`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.struveh

:func:`struvel`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.struvel


Anger-Weber functions
...................................................

:func:`angerj`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.angerj

:func:`webere`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.webere


Lommel functions
...................................................

:func:`lommels1`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.lommels1

:func:`lommels2`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.lommels2

Airy and Scorer functions
...............................................

:func:`airyai`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.airyai(z, derivative=0, **kwargs)

:func:`airybi`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.airybi(z, derivative=0, **kwargs)

:func:`airyaizero`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.airyaizero(k, derivative=0)

:func:`airybizero`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.airybizero(k, derivative=0, complex=0)

:func:`scorergi`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.scorergi(z, **kwargs)

:func:`scorerhi`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.scorerhi(z, **kwargs)


Coulomb wave functions
...............................................

:func:`coulombf`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.coulombf(l,eta,z)

:func:`coulombg`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.coulombg(l,eta,z)

:func:`coulombc`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.coulombc(l,eta)

Confluent U and Whittaker functions
...................................

:func:`hyperu`
^^^^^^^^^^^^^^
.. autofunction:: mpmath.hyperu(a, b, z)

:func:`whitm`
^^^^^^^^^^^^^^
.. autofunction:: mpmath.whitm(k,m,z)

:func:`whitw`
^^^^^^^^^^^^^^
.. autofunction:: mpmath.whitw(k,m,z)

Parabolic cylinder functions
.................................

:func:`pcfd`
^^^^^^^^^^^^^^
.. autofunction:: mpmath.pcfd(n,z,**kwargs)

:func:`pcfu`
^^^^^^^^^^^^^^
.. autofunction:: mpmath.pcfu(a,z,**kwargs)

:func:`pcfv`
^^^^^^^^^^^^^^
.. autofunction:: mpmath.pcfv(a,z,**kwargs)

:func:`pcfw`
^^^^^^^^^^^^^^
.. autofunction:: mpmath.pcfw(a,z,**kwargs)
