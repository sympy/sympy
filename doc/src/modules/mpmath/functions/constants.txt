Mathematical constants
----------------------

Mpmath supports arbitrary-precision computation of various common (and less common) mathematical constants. These constants are implemented as lazy objects that can evaluate to any precision. Whenever the objects are used as function arguments or as operands in arithmetic operations, they automagically evaluate to the current working precision. A lazy number can be converted to a regular ``mpf`` using the unary ``+`` operator, or by calling it as a function::

    >>> from mpmath import *
    >>> mp.dps = 15
    >>> pi
    <pi: 3.14159~>
    >>> 2*pi
    mpf('6.2831853071795862')
    >>> +pi
    mpf('3.1415926535897931')
    >>> pi()
    mpf('3.1415926535897931')
    >>> mp.dps = 40
    >>> pi
    <pi: 3.14159~>
    >>> 2*pi
    mpf('6.283185307179586476925286766559005768394338')
    >>> +pi
    mpf('3.141592653589793238462643383279502884197169')
    >>> pi()
    mpf('3.141592653589793238462643383279502884197169')

Exact constants
...............

The predefined objects :data:`j` (imaginary unit), :data:`inf` (positive infinity) and :data:`nan` (not-a-number) are shortcuts to :class:`mpc` and :class:`mpf` instances with these fixed values.

Pi (``pi``)
....................................

.. autoattribute:: mpmath.mp.pi

Degree (``degree``)
....................................

.. autoattribute:: mpmath.mp.degree

Base of the natural logarithm (``e``)
.....................................

.. autoattribute:: mpmath.mp.e

Golden ratio (``phi``)
......................

.. autoattribute:: mpmath.mp.phi


Euler's constant (``euler``)
............................

.. autoattribute:: mpmath.mp.euler

Catalan's constant (``catalan``)
................................

.. autoattribute:: mpmath.mp.catalan

Apery's constant (``apery``)
............................

.. autoattribute:: mpmath.mp.apery

Khinchin's constant (``khinchin``)
..................................

.. autoattribute:: mpmath.mp.khinchin

Glaisher's constant (``glaisher``)
..................................

.. autoattribute:: mpmath.mp.glaisher

Mertens constant (``mertens``)
..................................

.. autoattribute:: mpmath.mp.mertens

Twin prime constant (``twinprime``)
...................................

.. autoattribute:: mpmath.mp.twinprime
