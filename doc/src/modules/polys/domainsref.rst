.. _polys-domainsref:

===================================
Reference docs for the Poly Domains
===================================

This page lists the reference documentation for the domains in the polys
module. For a general introduction to the polys module it is recommended to
read :ref:`polys-basics` instead. For an introductory explanation of the
what the domain system is and how it is used it is recommended to read
:ref:`polys-domainsintro`. This page lists the reference docs for the
:py:class:`~.Domain` class and its subclasses (the specific domains such as
``ZZ``) as well as the classes that represent the domain elements.


Domains
=======

.. currentmodule:: sympy.polys.domains

Here we document the various implemented ground domains (see
:ref:`polys-domainsintro` for more of an explanation). There are three types
of :py:class:`~.Domain` subclass: abstract domains, concrete domains, and
"implementation domains". Abstract domains cannot be (usefully) instantiated
at all, and just collect together functionality shared by many other domains.
Concrete domains are those meant to be instantiated and used in the polynomial
manipulation algorithms. In some cases, there are various possible ways to
implement the data type the domain provides. For example, depending on what
libraries are available on the system, the integers are implemented either
using the python built-in integers, or using gmpy. Note that various aliases
are created automatically depending on the libraries available. As such e.g.
``ZZ`` always refers to the most efficient implementation of the integer ring
available.


Abstract Domains
================

.. autoclass:: sympy.polys.domains.domain.Domain
   :members:

.. autoclass:: sympy.polys.domains.domainelement.DomainElement
   :members:

.. autoclass:: sympy.polys.domains.field.Field
   :members:

.. autoclass:: sympy.polys.domains.ring.Ring
   :members:

.. autoclass:: sympy.polys.domains.simpledomain.SimpleDomain
   :members:

.. autoclass:: sympy.polys.domains.compositedomain.CompositeDomain
   :members:


.. _GF(p):

GF(p)
=====

.. autoclass:: FiniteField
   :members:

.. autoclass:: PythonFiniteField
   :members:

.. autoclass:: GMPYFiniteField
   :members:

.. _ZZ:


ZZ
==

The :ref:`ZZ` domain represents the `integers`_ `\mathbb{Z}` as a
:py:class:`~.Domain` in the domain system (see :ref:`polys-domainsintro`).

By default a :py:class:`~.Poly` created from an expression with integer
coefficients will have the domain :ref:`ZZ`::

  >>> from sympy import Poly, Symbol
  >>> x = Symbol('x')
  >>> p = Poly(x**2 + 1)
  >>> p
  Poly(x**2 + 1, x, domain='ZZ')
  >>> p.domain
  ZZ

The corresponding `field of fractions`_ is the domain of the rationals
:ref:`QQ`. Conversely :ref:`ZZ` is the `ring of integers`_ of :ref:`QQ`::

  >>> from sympy import ZZ, QQ
  >>> ZZ.get_field()
  QQ
  >>> QQ.get_ring()
  ZZ

When using the domain directly :ref:`ZZ` can be used as a constructor to
create instances which then support the operations ``+,-,*,**,//,%`` (true
division ``/`` should not be used with :ref:`ZZ` - see the
:py:meth:`~.Domain.exquo` domain method)::

  >>> x = ZZ(5)
  >>> y = ZZ(2)
  >>> x // y  # floor division
  2
  >>> x % y   # modulo division (remainder)
  1

The :py:meth:`~.Domain.gcd` method can be used to compute the `gcd`_ of any
two elements::

  >>> ZZ.gcd(ZZ(10), ZZ(2))
  2

There are two implementations of :ref:`ZZ` in SymPy. If ``gmpy`` or ``gmpy2``
is installed then :ref:`ZZ` will be implemented by :py:class:`GMPYIntegerRing`
and the elements will be instances of the ``gmpy.mpz`` type. Otherwise if
``gmpy`` and ``gmpy2`` are not installed then :ref:`ZZ` will be implemented by
:py:class:`PythonIntegerRing` which uses Python's standard builtin ``int``
type. With larger integers ``gmpy`` can be more efficient so it is preferred
when available.

.. autoclass:: IntegerRing
   :members:
   :exclude-members: dtype, tp

.. autoclass:: PythonIntegerRing
.. autoclass:: GMPYIntegerRing
   :members:
   :exclude-members: dtype, tp

.. _QQ:


QQ
==

The :ref:`QQ` domain represents the `rationals`_ `\mathbb{Q}` as a
:py:class:`~.Domain` in the domain system (see :ref:`polys-domainsintro`).

By default a :py:class:`~.Poly` created from an expression with rational
coefficients will have the domain :ref:`QQ`::

  >>> from sympy import Poly, Symbol
  >>> x = Symbol('x')
  >>> p = Poly(x**2 + x/2)
  >>> p
  Poly(x**2 + 1/2*x, x, domain='QQ')
  >>> p.domain
  QQ

The corresponding `ring of integers`_ is the :py:class:`~.Domain` of the
integers :ref:`ZZ`. Conversely :ref:`QQ` is the `field of fractions`_ of
:ref:`ZZ`::

  >>> from sympy import ZZ, QQ
  >>> QQ.get_ring()
  ZZ
  >>> ZZ.get_field()
  QQ

When using the domain directly :ref:`QQ` can be used as a constructor to
create instances which then support the operations ``+,-,*,**,/`` (true
division ``/`` is always possible for nonzero divisors in :ref:`QQ`)::

  >>> x = QQ(5)
  >>> y = QQ(2)
  >>> x / y  # true division
  5/2

There are two implementations of :ref:`QQ` in SymPy. If ``gmpy`` or ``gmpy2``
is installed then :ref:`QQ` will be implemented by
:py:class:`GMPYRationalField` and the elements will be instances of the
``gmpy.mpq`` type. Otherwise if ``gmpy`` and ``gmpy2`` are not installed then
:ref:`QQ` will be implemented by :py:class:`PythonRationalField` which is a
pure Python class as part of sympy. The ``gmpy`` implementation is
preferred because it is significantly faster.

.. autoclass:: RationalField
   :members:
   :exclude-members: dtype, tp

.. autoclass:: PythonRationalField
.. autoclass:: GMPYRationalField
   :members:
   :exclude-members: dtype, tp

.. autoclass:: sympy.external.pythonmpq.PythonMPQ


.. _MPQ:


MPQ
===

The ``MPQ`` type is either :py:class:`~.PythonMPQ` or otherwise the ``mpq``
type from ``gmpy2``.


Gaussian domains
================

The Gaussian domains :ref:`ZZ_I` and :ref:`QQ_I` share common superclasses
:py:class:`~.GaussianElement` for the domain elements and
:py:class:`~.GaussianDomain` for the domains themselves.

.. autoclass:: sympy.polys.domains.gaussiandomains.GaussianDomain
   :members:
.. autoclass:: sympy.polys.domains.gaussiandomains.GaussianElement
   :members:


.. _ZZ_I:


ZZ_I
====

.. autoclass:: sympy.polys.domains.gaussiandomains.GaussianIntegerRing
   :members:
.. autoclass:: sympy.polys.domains.gaussiandomains.GaussianInteger
   :members:

.. _QQ_I:


QQ_I
====

.. autoclass:: sympy.polys.domains.gaussiandomains.GaussianRationalField
   :members:
.. autoclass:: sympy.polys.domains.gaussiandomains.GaussianRational
   :members:

.. _QQ(a):


QQ<a>
=====

.. autoclass:: AlgebraicField
   :members:

.. _RR:


RR
==

.. autoclass:: RealField
   :members:

.. _CC:


CC
==

.. autoclass:: ComplexField
   :members:

.. _K[x]:


K[x]
====

.. autoclass:: PolynomialRing
   :members:

.. _K(x):


K(x)
====

.. autoclass:: FractionField
   :members:

.. _EX:


EX
==

.. autoclass:: ExpressionDomain
   :members:

.. autoclass:: sympy.polys.domains.expressiondomain::ExpressionDomain.Expression
   :members:


Quotient ring
=============

.. autoclass:: sympy.polys.domains.quotientring.QuotientRing


Sparse polynomials
==================

.. currentmodule:: sympy.polys.rings

Sparse polynomials are represented as dictionaries.

.. autofunction:: ring
.. autofunction:: xring
.. autofunction:: vring
.. autofunction:: sring

.. autoclass:: PolyRing
   :members:

.. autoclass:: PolyElement
   :members:


Sparse rational functions
=========================

.. currentmodule:: sympy.polys.fields

Sparse polynomials are represented as dictionaries.

.. autofunction:: field
.. autofunction:: xfield
.. autofunction:: vfield
.. autofunction:: sfield

.. autoclass:: FracField
   :members:

.. autoclass:: FracElement
   :members:


Dense polynomials
=================

.. currentmodule:: sympy.polys.polyclasses

.. autoclass:: DMP
   :members:

.. autoclass:: DMF
   :members:

.. autoclass:: ANP
   :members:


.. _integers: https://en.wikipedia.org/wiki/Integer
.. _rationals: https://en.wikipedia.org/wiki/Rational_number
.. _gcd: https://en.wikipedia.org/wiki/Greatest_common_divisor
.. _field of fractions: https://en.wikipedia.org/wiki/Field_of_fractions
.. _ring of integers: https://en.wikipedia.org/wiki/Ring_of_integers
