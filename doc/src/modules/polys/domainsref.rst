.. _polys-domainsref:

===================================
Reference docs for the Poly Domains
===================================

This page lists the reference documentation for the domains in the polys
module. For a general introduction to the polts module it is recommended to
read :ref:`polys-basics` instead. For an introductoray explanation of the
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
****************

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

Concrete Domains
****************

.. autoclass:: FiniteField
   :members:

.. autoclass:: IntegerRing
   :members:

.. autoclass:: PolynomialRing
   :members:

.. autoclass:: RationalField
   :members:

.. autoclass:: AlgebraicField
   :members:

.. autoclass:: FractionField
   :members:

.. autoclass:: RealField
   :members:

.. autoclass:: ExpressionDomain
   :members:

Implementation Domains
**********************

.. autoclass:: PythonFiniteField
.. autoclass:: GMPYFiniteField

.. autoclass:: PythonIntegerRing
.. autoclass:: GMPYIntegerRing

.. autoclass:: PythonRationalField
.. autoclass:: GMPYRationalField

Sparse polynomials
******************

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

Dense polynomials
*****************

.. currentmodule:: sympy.polys.polyclasses

.. autoclass:: DMP
   :members:

.. autoclass:: DMF
   :members:

.. autoclass:: ANP
   :members:
