.. _polys-internals:

===============================================
Internals of the Polynomial Manipulation Module
===============================================

The implementation of the polynomials module is structured internally in
"levels". There are four levels, called L0, L1, L2 and L3. The levels three
and four contain the user-facing functionality and were described in the
previous section. This section focuses on levels zero and one.

Level zero provides core polynomial manipulation functionality with C-like,
low-level interfaces. Level one wraps this low-level functionality into object
oriented structures. These are *not* the classes seen by the user, but rather
classes used internally throughout the polys module.

There is one additional complication in the implementation. This comes from the
fact that all polynomial manipulations are relative to a *ground domain*. For
example, when factoring a polynomial like `x^{10} - 1`, one has to decide what
ring the coefficients are supposed to belong to, or less trivially, what
coefficients are allowed to appear in the factorization. This choice of
coefficients is called a ground domain. Typical choices include the integers
`\mathbb{Z}`, the rational numbers `\mathbb{Q}` or various related rings and
fields. But it is perfectly legitimate (although in this case uninteresting)
to factorize over polynomial rings such as `k[Y]`, where `k` is some fixed
field.

Thus the polynomial manipulation algorithms (both
complicated ones like factoring, and simpler ones like addition or
multiplication) have to rely on other code to manipulate the coefficients.
In the polynomial manipulation module, such code is encapsulated in so-called
"domains". A domain is basically a factory object: it takes various
representations of data, and converts them into objects with unified interface.
Every object created by a domain has to implement the arithmetic operations
`+`, `-` and `\times`. Other operations are accessed through the domain, e.g.
as in ``ZZ.quo(ZZ(4), ZZ(2))``.

Note that there is some amount of *circularity*: the polynomial ring domains
use the level one classes, the level one classes use the level zero functions,
and level zero functions use domains. It is possible, in principle, but not in
the current implementation, to work in rings like `k[X][Y]`. This would create
even more layers. For this reason, working in the isomorphic ring `k[X, Y]`
is preferred.

Basic usage of domains
======================

Several domains are predefined and ready to be used such as ``ZZ`` and ``QQ``
which represent the ring of integers `\mathbb{Z}` and the field of rationals
`\mathbb{Q}`. The :py:class:`~.Domain` object is used to construct elements
which can then be used in ordinary arithmetic operations.::

  >>> from sympy import ZZ
  >>> z1 = ZZ(2)
  >>> z1
  2
  >>> z1 + z1
  4
  >>> type(z1)
  <class 'int'>
  >>> z1 in ZZ
  True

The basic operations ``+``, ``-``, and ``*`` for addition, subtraction and
multiplication will work for the elements of any domain and will produce new
domain elements. Division with ``/`` (Python's "true division" operator) is not
possible for all domains and should not be used with domain elements unless
the domain is known to be a field. For example dividing two elements of ``ZZ``
gives a ``float`` which is not an element of ``ZZ``::

  >>> z1 / z1
  1.0
  >>> type(z1 / z1)
  <class 'float'>
  >>> ZZ.is_Field
  False

Most domains representing non-field rings allow floor and modulo division
(remainder) with Python's floor division ``//`` and modulo division ``%``
operators. For example with ``ZZ``::

  >>> z1 // z1
  1
  >>> z1 % z1
  0

The ``QQ`` domain represents the field of rational numbers and does allow
division::

  >>> from sympy import QQ
  >>> q1 = QQ(1, 2)
  >>> q1
  1/2
  >>> q2 = QQ(2, 3)
  >>> q2
  2/3
  >>> q1 / q2
  3/4
  >>> type(q1)
  <class 'sympy.polys.domains.pythonrational.PythonRational'>

In general code that is expected to work with elements of an arbitrary domain
should not use the division operators ``/``, ``//`` and ``%``. Only the operators
``+``, ``-``, ``*`` and ``**`` (with non-negative exponents) should be assumed
to work with arbitrary domain elements. All other operations should be
accessed as functions from the :py:class:`~.Domain` object::

  >>> ZZ.quo(ZZ(5), ZZ(3))  # 5 // 3
  1
  >>> ZZ.rem(ZZ(5), ZZ(3))  # 5 % 3
  2
  >>> ZZ.div(ZZ(5), ZZ(3))  # divmod(5, 3)
  (1, 2)
  >>> QQ.div(QQ(5), QQ(3))
  (5/3, 0)

The :py:meth:`~.Domain.exquo` function is used to compute an exact quotient.
This is the analogue of ``a / b`` but where the division is expected to be exact
(with no remainder) or an error will be raised::

  >>> QQ.exquo(QQ(5), QQ(3))
  5/3
  >>> ZZ.exquo(ZZ(4), ZZ(2))
  2
  >>> ZZ.exquo(ZZ(5), ZZ(3))
  Traceback (most recent call last):
  ...
  sympy.polys.polyerrors.ExactQuotientFailed: 3 does not divide 5 in ZZ

The exact methods and attributes of the domain elements are not guaranteed in
general beyond the basic arithmetic operations. It should not be presumed that
e.g. ``ZZ`` will always be of type ``int``. If ``gmpy`` or ``gmpy2`` is
installed then the ``mpz`` or ``mpq`` types are used instead for ``ZZ`` and
``QQ``::

  >>> from sympy import ZZ, QQ
  >>> ZZ(2)  # doctest: +SKIP
  mpz(2)
  >>> QQ(2, 3)  # doctest: +SKIP
  mpq(2, 3)

The ``mpz`` type is faster than Python's standard ``int`` type for operations
with large integers although for smaller integers the difference is not so
significant. The ``mpq`` type representing rational numbers is implemented in
C rather than Python and is many times faster than the pure Python
implementation of ``QQ`` that is used when gmpy is not installed.

In general the Python type used for the elements of a domain can be checked
from the ``dtype`` attribute of the domain. The :py:meth:`~.Domain.of_type`
method can be used to check if an object is an instance of ``dtype``.::

  >>> z = ZZ(2)
  >>> type(z)
  <class 'int'>
  >>> ZZ.dtype
  <class 'int'>
  >>> isinstance(z, ZZ.dtype)
  True
  >>> ZZ.of_type(z)
  True

Domain elements vs sympy expressions
====================================

Note that domain elements are not of the same type as ordinary sympy
expressions which are subclasses of :py:class:`~.Basic` such as
:py:class:`~sympy.core.numbers.Integer`. Ordinary sympy expressions are
created with the :py:func:`~sympy.core.sympify.sympify` function.::

  >>> from sympy import sympify
  >>> z1_sympy = sympify(2)  # Normal sympy object
  >>> z1_sympy
  2
  >>> type(z1_sympy)
  <class 'sympy.core.numbers.Integer'>
  >>> from sympy import Basic
  >>> isinstance(z1_sympy, Basic)
  True

It is important when working with the domains not to mix sympy expressions
with domain elements even though it will sometimes work in simple cases. Each
domain object has the methods :py:meth:`~.Domain.to_sympy` and
:py:meth:`~.Domain.from_sympy` for converting back and forther between sympy
expressions and domain elements::

  >>> z_sympy = sympify(2)
  >>> z_zz = ZZ.from_sympy(z_sympy)
  >>> z_zz
  2
  >>> type(z_sympy)
  <class 'sympy.core.numbers.Integer'>
  >>> type(z_zz)
  <class 'int'>
  >>> ZZ.to_sympy(z_zz)
  2
  >>> type(ZZ.to_sympy(z_zz))
  <class 'sympy.core.numbers.Integer'>

Any particular domain will only be able to represent some sympy expressions so
conversion will fail if the expression can not be represented in the domain::

  >>> from sympy import sqrt
  >>> e = sqrt(2)
  >>> e
  sqrt(2)
  >>> ZZ.from_sympy(e)
  Traceback (most recent call last):
  ...
  sympy.polys.polyerrors.CoercionFailed: expected an integer, got sqrt(2)

It is important not to mix domain elements with other Python types such as
``int``, ``float``, as well as standard sympy :py:class:`~.Basic` expressions.
When working in a domain care should be taken as some Python operations will
do this implicitly. for example the ``sum`` function will use the regular
``int`` value of zero so that ``sum([a, b])`` is effectively evaluated as ``(0
+ a) + b`` where ``0`` is of type ``int``.

Every domain is at least a ring if not a field and as such is guaranteed to
have two elements in particular corresponding to `1` and `0`.  The domain
object provides domain elements for these as the attributes ``one`` and
``zero``. These are useful for something like Python's ``sum`` function which
allows to provide an alternative object as the "zero"::

  >>> ZZ.one
  1
  >>> ZZ.zero
  0
  >>> sum([ZZ(1), ZZ(2)])  # don't do this (even it sometimes works)
  3
  >>> sum([ZZ(1), ZZ(2)], ZZ.zero) # provide the zero from the domain
  3

A standard pattern then for performing calculations in a domain is:

#. Start with sympy :py:class:`~.Basic` instances representing expressions.
#. Choose an appropriate domain that can represent the expressions.
#. Convert all expressions to domain elements using
   :py:meth:`~.Domain.from_sympy`.
#. Perform the calculation with the domain elements.
#. Convert back to :py:class:`~.Basic` with :py:meth:`~.Domain.to_sympy`.

Here is an implementation of the ``sum`` function that illustrates these
steps and sums some integers but performs the calculation using the domain
elements rather than standard sympy expressions::

  def sum_domain(expressions_sympy):
      """Sum sympy expressions but performing calculations in domain ZZ"""

      # Convert to domain
      expressions_dom = [ZZ.from_sympy(e) for e in expressions_sympy]

      # Perform calculations in the domain
      result_dom = ZZ.zero
      for e_dom in expressions_dom:
          result_dom += e_dom

      # Convert the result back to Basic
      result_sympy = ZZ.to_sympy(result_dom)
      return result_sympy

Converting between domains
==========================

It is often useful to combine calculations performed over different domains.
However just as it is important to avoid mixing domain elements with normal
sympy expressions and other Python types it is also important to avoid mixing
elements from different domains. The :py:meth:`~.Domain.convert` method is
used to convert elements from one domain into elements of another domain::

  >>> num_zz = ZZ(3)
  >>> ZZ.of_type(num_zz)
  True
  >>> num_qq = QQ.convert(num_zz, ZZ)
  >>> ZZ.of_type(num_qq)
  False
  >>> QQ.of_type(num_qq)
  True

When looking to combine elements from two different domains and we want to
perform mixed calculations with them we need to

#. Choose a new domain that can represent all elements of both.
#. Convert all elements to the new domain.
#. Perform the calculation in the new domain.

The key question arising from point 1. is how to choose a domain that can
represent the elements of both domains. For this there is the
:py:meth:`~.Domain.unify` method::

  >>> x1, K1 = ZZ(2), ZZ
  >>> y2, K2 = QQ(3, 2), QQ
  >>> K1
  ZZ
  >>> K2
  QQ
  >>> K3 = K1.unify(K2)
  >>> K3
  QQ
  >>> x3 = K3.convert(x1, K1)
  >>> y3 = K3.convert(y2, K2)
  >>> x3 + y3
  7/2

The :py:meth:`~.Domain.unify` method will find a domain that encompasses both
domains so in this example ``ZZ.unify(QQ)`` gives ``QQ`` because every element
of ``ZZ`` can be represented as an element of ``QQ``. This means that all
inputs (``x1`` and ``y2``) can be converted to the elements of the common
domain ``K3`` (as ``x3`` and ``y3``). Once in the common domain we can safely
use arithmetic operations like ``+``. In this example one domain is a superset
of the other and we see that ``K1.unify(K2) == K2`` so it is not actually
necessary to convert ``y2``. In general though ``K1.unify(K2)`` can give a new
domain that is not equal to either ``K1`` or ``K2``.

The :py:meth:`~.Domain.convert` method can be called without specifying the
source domain as the second ergument e.g.::

  >>> QQ.convert(ZZ(2))
  2

This works because :py:meth:`~.Domain.convert` can check the type of ``ZZ(2)``
and can try to work out what domain (``ZZ``) it is an element of. Certain
domains like ``ZZ`` and ``QQ`` are treated as special cases to make this work.
Elements of more complicated domains are instances of subclass of
:py:class:`~.DomainElement` which has a :py:meth:`~.DomainElement.parent`
method that can identify the domain that the element belongs to. For example
in the polynomial ring ``ZZ[x]`` we have::

  >>> from sympy import ZZ, Symbol
  >>> x = Symbol('x')
  >>> K = ZZ[x]
  >>> K
  ZZ[x]
  >>> p = K(x) + K.one
  >>> p
  x + 1
  >>> type(p)
  <class 'sympy.polys.rings.PolyElement'>
  >>> p.parent()
  ZZ[x]
  >>> p.parent() == K
  True

It is more efficient though to call :py:meth:`~.Domain.convert` with
the source domain specified as the second argument::

  >>> QQ.convert(ZZ(2), ZZ)
  2

Polynomial ring domains
=======================

So far we have seen the domains ``ZZ`` and ``QQ`` which are the simplest
possible examples (and the most commonly used). The complicated domain system
was not created just to represent simple examples such as this. A more
interesting example that was introduced in the last section is that of a
polynomial ring like ``ZZ[x]`` which is the domain of polynomials in the
generator ``x`` with coefficients over ``ZZ``::

  >>> from sympy import ZZ, symbols
  >>> x = symbols('x')
  >>> K = ZZ[x]
  >>> K
  ZZ[x]
  >>> x_dom = K(x)
  >>> x_dom + K.one
  x + 1

All the operations discussed before will work with elements of a polynomial
ring::

  >>> p = x_dom + K.one
  >>> p
  x + 1
  >>> p + p
  2*x + 2
  >>> p - p
  0
  >>> p * p
  x**2 + 2*x + 1
  >>> p ** 3
  x**3 + 3*x**2 + 3*x + 1
  >>> K.exquo(x_dom**2 - K.one, x_dom - K.one)
  x + 1

The internal representation of elements of ``K[x]`` is different from the way
that ordinary sympy (:py:class:`~.Basic`) expressions are represented. The
:py:class:`~.Basic` representation of any expression is as a tree e.g.::

  >>> from sympy import srepr
  >>> K = ZZ[x]
  >>> p_basic = x**2 + 2*x + 1
  >>> p_basic
  x**2 + 2*x + 1
  >>> srepr(p_basic)
  "Add(Pow(Symbol('x'), Integer(2)), Mul(Integer(2), Symbol('x')), Integer(1))"

Here the expression is a tree where the top node is an :py:class:`~.Add` and
its children nodes are :py:class:`~.Pow` etc. This tree representation makes
is possible to represent equivalent expressions in different ways e.g.::

  >>> x = symbols('x')
  >>> p_basic = x*(x + 1) + x
  >>> p_basic
  x*(x + 1) + x
  >>> p_basic.expand()
  x**2 + 2*x

By contrast the domain ``ZZ[x]`` represents only polynomials and does so by
simply storing the non-zero coefficients of the expanded polynomial (the
"sparse" polynomial representation). In particular elements of ``ZZ[x]`` are
represented as a Python ``dict``. Their type is :py:class:`~.PolyElement`
which is a subclass of ``dict``. Converting to a normal dict shows the
internal representation::

  >>> x = symbols('x')
  >>> K = ZZ[x]
  >>> x_dom = K(x)
  >>> p_dom = K(3)*x_dom**2 + K(2)*x_dom + K(7)
  >>> p_dom
  3*x**2 + 2*x + 7
  >>> dict(p_dom)
  {(0,): 7, (1,): 2, (2,): 3}

This internal form makes it impossible to represent unexpanded multiplications
so any multiplication of elements of ``ZZ[x]`` will always be
expanded::

  >>> x = symbols('x')
  >>> K = ZZ[x]
  >>> x_dom = K(x)
  >>> p_basic = x * (x + 1) + x
  >>> p_basic
  x*(x + 1) + x
  >>> p_dom = x_dom * (x_dom + K.one) + x_dom
  >>> p_dom
  x**2 + 2*x

These same considerations apply to powers::

  >>> (x + 1) ** 2
  (x + 1)**2
  >>> (x_dom + K.one) ** 2
  x**2 + 2*x + 1

We can also construct multivariate polynomial rings::

  >>> x, y = symbols('x, y')
  >>> K = ZZ[x,y]
  >>> xk = K(x)
  >>> yk = K(y)
  >>> xk**2*yk + xk + yk
  x**2*y + x + y

The :py:meth:`~.Domain.unify` method understands how to combine different
polynomial ring domains and how to unify the base domain::

  >>> ZZ[x].unify(ZZ[y])
  ZZ[x,y]
  >>> ZZ[x,y].unify(ZZ[y])
  ZZ[x,y]
  >>> ZZ[x].unify(QQ)
  QQ[x]

It is also possible to construct nested polynomial rings (although it is less
efficient). The ring ``K[x][y]`` is formally equivalent to ``K[x,y]`` although
the implementations in sympy are different::

  >>> K = ZZ[x][y]
  >>> p = K(x**2 + x*y + y**2)
  >>> p
  y**2 + x*y + x**2
  >>> dict(p)
  {(0,): x**2, (1,): x, (2,): 1}

Here the coefficients like ``x**2`` are instances of :py:class:`~.PolyElement`
as well so this is a ``dict`` where the values are also dicts. The
multivariate ring domain ``ZZ[x,y]`` has a more efficient representation as a
single flattened ``dict``::

  >>> K = ZZ[x,y]
  >>> p = K(x**2 + x*y + y**2)
  >>> p
  x**2 + x*y + y**2
  >>> dict(p)
  {(0, 2): 1, (1, 1): 1, (2, 0): 1}

The difference in efficiency between these representations grows as the number
of generators increases i.e. ``ZZ[x,y,z,t,...]`` vs ``ZZ[x][y][z][t]...``.

Algebraic number fields
=======================

A polynomial ring such as ``ZZ[x]`` is an example of a transcendental
extension but algebraic extensions are also possible. An algebraic extension
of ``QQ`` is known as an *algebraic number field* and these are implemented in
sympy. The natural syntax for these would be something like ``QQ(sqrt(2))``
however ``QQ()`` is already overloaded as the constructor for elements of
``QQ``. These domains are instead created as e.g.
``QQ.algebraic_field(sqrt(2))``. The resulting field will be an instance of
:py:class:`~.AlgebraicField`.

The printing support for these is less developed but we can use
:py:meth:`~.Domain.to_sympy` to take advantage of the corresponding
:py:class:`~.Basic` printing support::

  >>> K = QQ.algebraic_field(sqrt(2))
  >>> K
  QQ<sqrt(2)>
  >>> b = K.one + K.convert(sqrt(2))
  >>> b
  ANP([1, 1], [1, 0, -2], QQ)
  >>> K.to_sympy(b)
  1 + sqrt(2)
  >>> b ** 2
  ANP([2, 3], [1, 0, -2], QQ)
  >>> K.to_sympy(b**2)
  2*sqrt(2) + 3

The raw printed display immediately shows the internal representation of the
elements as :py:class:`~.ANP` instances. The field
`\mathbb{Q}(\sqrt{2})` consists of numbers of the form
`a+b\sqrt{2}` where `a` and `b` are rational numbers. Consequently every
number in this field can be represented as a pair ``(a, b)`` of elements of
``QQ``. The domain element stores these two in a list and also stores a list
representation of the *minimal polynomial* for the extension element
`\sqrt{2}`. There is a sympy function :py:func:`~.minpoly` that can compute
the minimal polynomial of any algebraic expression over the rationals::

  >>> from sympy import minpoly
  >>> minpoly(sqrt(2), x)
  x**2 - 2

In the dense polynomial representation as a list of coefficients this
polynomial is represented as ``[1, 0, -2]`` as seen in the :py:class:`~.ANP`
display for the elements of ``QQ<sqrt(2)>`` above.

It is also possible to create an algebraic number field with multiple
generators such as `\mathbb{Q}(\sqrt{2},\sqrt{3})`::

  >>> K = QQ.algebraic_field(sqrt(2), sqrt(3))
  >>> K
  QQ<sqrt(2) + sqrt(3)>
  >>> sqrt2 = K.convert(sqrt(2))
  >>> sqrt3 = K.convert(sqrt(3))
  >>> p = (K.one + sqrt2) * (K.one + sqrt3)
  >>> p
  ANP([1/2, 1, -3/2], [1, 0, -10, 0, 1], QQ)
  >>> K.to_sympy(p)
  1 + sqrt(2) + sqrt(3) + sqrt(6)
  >>> K.to_sympy(p**2)
  4*sqrt(6) + 6*sqrt(3) + 8*sqrt(2) + 12

Here the algebraic extension `\mathbb{Q}(\sqrt{2},\sqrt{3})` is converted to
the (isomorphic) `\mathbb{Q}(\sqrt{2}+\sqrt{3})` with a single generator
`\sqrt{2}+\sqrt{3}`. It is always possible to find a single generator like
this due to the *primitive element theorem*. There is a sympy function
:py:func:`~.primitive_element` that can compute the minimal polynomial for a
primitive element of an extension::

  >>> from sympy import primitive_element, minpoly
  >>> e = primitive_element([sqrt(2), sqrt(3)], x)
  >>> e[0]
  x**4 - 10*x**2 + 1
  >>> e[0].subs(x, sqrt(2) + sqrt(3)).expand()
  0

The minimal polynomial ``x**4 - 10*x**2 + 1`` has the dense list representation
``[1, 0, -10, 0, 1]`` as seen in the :py:class:`~.ANP` output above. What the
primitive element theorem means is that all algebraic number fields can be
represented as an extension of the rationals by a single generator with some
minimal polynomial. Calculations over the algebraic number field only need to
take advantage of the minimal polynomial and that makes it possible to compute
all arithmetic operations and also to carry out higher level operations like
factorisation of polynomials.

Rational function fields
========================

Some domains are classified as fields and others are not. It is generally
possible to convert any domain to a field that contains that domain with the
:py:meth:`~.Domain.get_field` method::

  >>> from sympy import ZZ, QQ, symbols
  >>> x, y = symbols('x, y')
  >>> ZZ.is_Field
  False
  >>> QQ.is_Field
  True
  >>> QQ[x]
  QQ[x]
  >>> QQ[x].is_Field
  False
  >>> QQ[x].get_field()
  QQ(x)
  >>> QQ[x].get_field().is_Field
  True
  >>> QQ.frac_field(x)
  QQ(x)

This introduces a new kind of domain ``QQ(x)``. It is not possible to
construct the domain ``QQ(x)`` with the same syntax so the easiest ways to
create it are using either ``QQ.frac_field(x)`` or ``QQ[x].get_field()`` where
the :py:meth:`~.Domain.frac_field` method is the more direct approach.

The rational function field ``QQ(x)`` is an instance of
:py:class:`~.RationalField`. This domain represents functions of the form
`p(x) / q(x)` for polynomials `p` and `q`. The domain elements are
represented as a pair of polynomials in ``QQ[x]``::

  >>> K = QQ.frac_field(x)
  >>> xk = K(x)
  >>> f = xk / (K.one + xk**2)
  >>> f
  x/(x**2 + 1)
  >>> f.numer
  x
  >>> f.denom
  x**2 + 1
  >>> QQ[x].of_type(f.numer)
  True
  >>> QQ[x].of_type(f.denom)
  True

Cancellation between the numerator and denominator is automatic in this
field::

  >>> p1 = xk**2 - 1
  >>> p2 = xk - 1
  >>> p1
  x**2 - 1
  >>> p2
  x - 1
  >>> p1 / p2
  x + 1

Computing this cancellation can be slow which makes rational function fields
potentially slower than polynomial rings or algebraic fields.

Domains
=======

.. currentmodule:: sympy.polys.domains

Here we document the various implemented ground domains. There are three
types: abstract domains, concrete domains, and "implementation domains".
Abstract domains cannot be (usefully) instantiated at all, and just collect
together functionality shared by many other domains. Concrete domains are
those meant to be instantiated and used in the polynomial manipulation
algorithms. In some cases, there are various possible ways to implement the
data type the domain provides. For example, depending on what libraries are
available on the system, the integers are implemented either using the python
built-in integers, or using gmpy. Note that various aliases are created
automatically depending on the libraries available. As such e.g. ``ZZ`` always
refers to the most efficient implementation of the integer ring available.

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

Level One
=========

.. currentmodule:: sympy.polys.polyclasses

.. autoclass:: DMP
   :members:

.. autoclass:: DMF
   :members:

.. autoclass:: ANP
   :members:

Level Zero
==========

Level zero contains the bulk code of the polynomial manipulation module.

Manipulation of dense, multivariate polynomials
***********************************************

These functions can be used to manipulate polynomials in `K[X_0, \ldots, X_u]`.
Functions for manipulating multivariate polynomials in the dense representation
have the prefix ``dmp_``. Functions which only apply to univariate polynomials
(i.e. `u = 0`)
have the prefix ``dup__``. The ground domain `K` has to be passed explicitly.
For many multivariate polynomial manipulation functions also the level `u`,
i.e. the number of generators minus one, has to be passed.
(Note that, in many cases, ``dup_`` versions of functions are available, which
may be slightly more efficient.)

**Basic manipulation:**

.. currentmodule:: sympy.polys.densebasic

.. autofunction:: dmp_LC
.. autofunction:: dmp_TC
.. autofunction:: dmp_ground_LC
.. autofunction:: dmp_ground_TC
.. autofunction:: dmp_true_LT
.. autofunction:: dmp_degree
.. autofunction:: dmp_degree_in
.. autofunction:: dmp_degree_list
.. autofunction:: dmp_strip
.. autofunction:: dmp_validate
.. autofunction:: dup_reverse
.. autofunction:: dmp_copy
.. autofunction:: dmp_to_tuple
.. autofunction:: dmp_normal
.. autofunction:: dmp_convert
.. autofunction:: dmp_from_sympy
.. autofunction:: dmp_nth
.. autofunction:: dmp_ground_nth
.. autofunction:: dmp_zero_p
.. autofunction:: dmp_zero
.. autofunction:: dmp_one_p
.. autofunction:: dmp_one
.. autofunction:: dmp_ground_p
.. autofunction:: dmp_ground
.. autofunction:: dmp_zeros
.. autofunction:: dmp_grounds
.. autofunction:: dmp_negative_p
.. autofunction:: dmp_positive_p
.. autofunction:: dmp_from_dict
.. autofunction:: dmp_to_dict
.. autofunction:: dmp_swap
.. autofunction:: dmp_permute
.. autofunction:: dmp_nest
.. autofunction:: dmp_raise
.. autofunction:: dmp_deflate
.. autofunction:: dmp_multi_deflate
.. autofunction:: dmp_inflate
.. autofunction:: dmp_exclude
.. autofunction:: dmp_include
.. autofunction:: dmp_inject
.. autofunction:: dmp_eject
.. autofunction:: dmp_terms_gcd
.. autofunction:: dmp_list_terms
.. autofunction:: dmp_apply_pairs
.. autofunction:: dmp_slice
.. autofunction:: dup_random

**Arithmetic operations:**

.. currentmodule:: sympy.polys.densearith

.. autofunction:: dmp_add_term
.. autofunction:: dmp_sub_term
.. autofunction:: dmp_mul_term
.. autofunction:: dmp_add_ground
.. autofunction:: dmp_sub_ground
.. autofunction:: dmp_mul_ground
.. autofunction:: dmp_quo_ground
.. autofunction:: dmp_exquo_ground
.. autofunction:: dup_lshift
.. autofunction:: dup_rshift
.. autofunction:: dmp_abs
.. autofunction:: dmp_neg
.. autofunction:: dmp_add
.. autofunction:: dmp_sub
.. autofunction:: dmp_add_mul
.. autofunction:: dmp_sub_mul
.. autofunction:: dmp_mul
.. autofunction:: dmp_sqr
.. autofunction:: dmp_pow
.. autofunction:: dmp_pdiv
.. autofunction:: dmp_prem
.. autofunction:: dmp_pquo
.. autofunction:: dmp_pexquo
.. autofunction:: dmp_rr_div
.. autofunction:: dmp_ff_div
.. autofunction:: dmp_div
.. autofunction:: dmp_rem
.. autofunction:: dmp_quo
.. autofunction:: dmp_exquo
.. autofunction:: dmp_max_norm
.. autofunction:: dmp_l1_norm
.. autofunction:: dmp_expand

**Further tools:**

.. currentmodule:: sympy.polys.densetools

.. autofunction:: dmp_integrate
.. autofunction:: dmp_integrate_in
.. autofunction:: dmp_diff
.. autofunction:: dmp_diff_in
.. autofunction:: dmp_eval
.. autofunction:: dmp_eval_in
.. autofunction:: dmp_eval_tail
.. autofunction:: dmp_diff_eval_in
.. autofunction:: dmp_trunc
.. autofunction:: dmp_ground_trunc
.. autofunction:: dup_monic
.. autofunction:: dmp_ground_monic
.. autofunction:: dup_content
.. autofunction:: dmp_ground_content
.. autofunction:: dup_primitive
.. autofunction:: dmp_ground_primitive
.. autofunction:: dup_extract
.. autofunction:: dmp_ground_extract
.. autofunction:: dup_real_imag
.. autofunction:: dup_mirror
.. autofunction:: dup_scale
.. autofunction:: dup_shift
.. autofunction:: dup_transform
.. autofunction:: dmp_compose
.. autofunction:: dup_decompose
.. autofunction:: dmp_lift
.. autofunction:: dup_sign_variations
.. autofunction:: dmp_clear_denoms
.. autofunction:: dmp_revert

Manipulation of dense, univariate polynomials with finite field coefficients
****************************************************************************
.. currentmodule:: sympy.polys.galoistools

Functions in this module carry the prefix ``gf_``, referring to the classical
name "Galois Fields" for finite fields. Note that many polynomial
factorization algorithms work by reduction to the finite field case, so having
special implementations for this case is justified both by performance, and by
the necessity of certain methods which do not even make sense over general
fields.

.. autofunction:: gf_crt
.. autofunction:: gf_crt1
.. autofunction:: gf_crt2
.. autofunction:: gf_int
.. autofunction:: gf_degree
.. autofunction:: gf_LC
.. autofunction:: gf_TC
.. autofunction:: gf_strip
.. autofunction:: gf_trunc
.. autofunction:: gf_normal
.. autofunction:: gf_from_dict
.. autofunction:: gf_to_dict
.. autofunction:: gf_from_int_poly
.. autofunction:: gf_to_int_poly
.. autofunction:: gf_neg
.. autofunction:: gf_add_ground
.. autofunction:: gf_sub_ground
.. autofunction:: gf_mul_ground
.. autofunction:: gf_quo_ground
.. autofunction:: gf_add
.. autofunction:: gf_sub
.. autofunction:: gf_mul
.. autofunction:: gf_sqr
.. autofunction:: gf_add_mul
.. autofunction:: gf_sub_mul
.. autofunction:: gf_expand
.. autofunction:: gf_div
.. autofunction:: gf_rem
.. autofunction:: gf_quo
.. autofunction:: gf_exquo
.. autofunction:: gf_lshift
.. autofunction:: gf_rshift
.. autofunction:: gf_pow
.. autofunction:: gf_pow_mod
.. autofunction:: gf_gcd
.. autofunction:: gf_lcm
.. autofunction:: gf_cofactors
.. autofunction:: gf_gcdex
.. autofunction:: gf_monic
.. autofunction:: gf_diff
.. autofunction:: gf_eval
.. autofunction:: gf_multi_eval
.. autofunction:: gf_compose
.. autofunction:: gf_compose_mod
.. autofunction:: gf_trace_map
.. autofunction:: gf_random
.. autofunction:: gf_irreducible
.. autofunction:: gf_irreducible_p
.. autofunction:: gf_sqf_p
.. autofunction:: gf_sqf_part
.. autofunction:: gf_sqf_list
.. autofunction:: gf_Qmatrix
.. autofunction:: gf_Qbasis
.. autofunction:: gf_berlekamp
.. autofunction:: gf_zassenhaus
.. autofunction:: gf_shoup
.. autofunction:: gf_factor_sqf
.. autofunction:: gf_factor
.. autofunction:: gf_value
.. autofunction:: gf_csolve

Manipulation of sparse, distributed polynomials and vectors
***********************************************************

Dense representations quickly require infeasible amounts of storage and
computation time if the number of variables increases. For this reason,
there is code to manipulate polynomials in a *sparse* representation.



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

In commutative algebra, one often studies not only polynomials, but also
*modules* over polynomial rings. The polynomial manipulation module provides
rudimentary low-level support for finitely generated free modules. This is
mainly used for Groebner basis computations (see there), so manipulation
functions are only provided to the extend needed. They carry the prefix
``sdm_``. Note that in examples, the generators of the free module are called
`f_1, f_2, \ldots`.

.. currentmodule:: sympy.polys.distributedmodules

.. autofunction:: sdm_monomial_mul
.. autofunction:: sdm_monomial_deg
.. autofunction:: sdm_monomial_divides
.. autofunction:: sdm_LC
.. autofunction:: sdm_to_dict
.. autofunction:: sdm_from_dict
.. autofunction:: sdm_add
.. autofunction:: sdm_LM
.. autofunction:: sdm_LT
.. autofunction:: sdm_mul_term
.. autofunction:: sdm_zero
.. autofunction:: sdm_deg
.. autofunction:: sdm_from_vector
.. autofunction:: sdm_to_vector

Polynomial factorization algorithms
***********************************

Many variants of Euclid's algorithm:

.. currentmodule:: sympy.polys.euclidtools

Classical remainder sequence
----------------------------

Let `K` be a field, and consider the ring `K[X]` of polynomials in a single
indeterminate `X` with coefficients in `K`. Given two elements `f` and `g`
of `K[X]` with `g\neq 0` there are unique polynomials `q` and `r` such that
`f = qg + r` and `\deg(r) < \deg(g)` or `r = 0`.
They are denoted by `\mathrm{quo}(f,g)`
(*quotient*) and `\mathrm{rem}(f,g)` (*remainder*), so we have
the *division identity*

.. math::

  f = \mathrm{quo}(f,g)g + \mathrm{rem}(f,g).

It follows that every ideal `I` of `K[X]` is a principal ideal, generated by
any element `\neq 0` of minimum degree (assuming `I` non-zero). In fact,
if `g` is such a polynomial and `f` is any element of `I`,
`\mathrm{rem}(f,g)` belongs to `I` as a linear combination of `f` and `g`,
hence must be zero; therefore `f` is a multiple of `g`.

Using this result it is possible to find a `greatest common
divisor <https://en.wikipedia.org/wiki/Greatest_common_divisor>`_
(gcd) of any polynomials `f,g,\ldots` in `K[X]`.
If `I` is the ideal formed by all linear combinations of the given polynomials
with coefficients in `K[X]`, and `d` is its generator,
then every common divisor of the polynomials also divides `d`.
On the other hand, the given polynomials are multiples of the generator `d`;
hence `d` is a gcd of the polynomials, denoted `\mathrm{gcd}(f,g,\ldots)`.

An algorithm for the gcd of two polynomials `f` and `g` in `K[X]` can
now be obtained as follows.
By the division identity, `r = \mathrm{rem}(f,g)` is in the ideal generated
by `f` and `g`, as well as `f` is in the ideal generated by `g` and `r`.
Hence the ideals generated by the pairs `(f,g)` and `(g,r)` are the same.
Set `f_0 = f`, `f_1 = g`, and define recursively
`f_i = \mathrm{rem}(f_{i-2},f_{i-1})` for `i\ge 2`.
The recursion ends after a finite number of steps with `f_{k+1}=0`,
since the degrees of the polynomials are strictly decreasing.
By the above remark, all the pairs `(f_{i-1},f_i)` generate the same ideal.
In particular, the ideal generated by `f` and `g` is generated by `f_k`
alone as `f_{k+1} = 0`. Hence `d = f_k` is a gcd of `f` and `g`.

The sequence of polynomials `f_0`, `f_1,\ldots, f_k` is called the
*Euclidean polynomial remainder sequence* determined by `(f,g)` because
of the analogy with the classical `Euclidean algorithm
<https://en.wikipedia.org/wiki/Euclidean_algorithm>`_ for the gcd of
natural numbers.

The algorithm may be extended to obtain an expression for `d` in terms of
`f` and `g` by using the full division identities
to write recursively each `f_i` as a linear combination of `f` and `g`.
This leads to an equation

.. math::

   d = uf + vg\qquad (u,v \in K[X])

analogous to `Bézout's identity
<https://en.wikipedia.org/wiki/B%C3%A9zout%27s_identity>`_
in the case of integers.

.. autofunction:: dmp_half_gcdex
.. autofunction:: dmp_gcdex
.. autofunction:: dmp_invert
.. autofunction:: dmp_euclidean_prs

Simplified remainder sequences
------------------------------

Assume, as is usual, that the coefficient field `K` is
the field of fractions of an integral domain `A`.
In this case the coefficients (numerators and denominators)
of the polynomials in the Euclidean remainder sequence
tend to grow very fast.

If `A` is a unique factorization domain, the coefficients may be
reduced by cancelling common factors of numerators and denominators.
Further reduction is possible noting that a gcd of polynomials in
`K[X]` is not unique:
it may be multiplied by any (non-zero) constant factor.

Any polynomial `f` in `K[X]` can be simplified by extracting
the denominators and common factors of the numerators of its coefficients.
This yields the representation `f = cF` where `c\in K` is
the *content* of `f` and `F` is a *primitive* polynomial, i.e.,
a polynomial in `A[X]` with coprime coefficients.

It is possible to start the algorithm by replacing the given polynomials
`f` and `g` with their primitive parts. This will only modify
`\mathrm{rem}(f,g)` by a constant factor.
Replacing it with its primitive part and continuing recursively
we obtain all the primitive parts of the polynomials in
the Euclidean remainder sequence, including the primitive
`\mathrm{gcd}(f,g)`.

This sequence is the *primitive polynomial remainder sequence*.
It is an example of *general polynomial remainder sequences* where
the computed remainders are modified by constant multipliers (or divisors)
in order to simplify the results.

.. autofunction:: dmp_primitive_prs

Subresultant sequence
---------------------

The coefficients of the primitive polynomial sequence do not grow
exceedingly, but the computation of the primitive parts requires
extra processing effort. Besides, the method only works with fraction fields
of unique factorization domains, excluding, for example, the general number
fields.

Collins [Collins67] realized that the so-called *subresultant polynomials*
of a pair of polynomials also form a generalized remainder sequence.
The coefficients of these polynomials
are expressible as determinants in the coefficients of the given
polynomials. Hence (the logarithm of) their size only grows linearly.
In addition, if the coefficients of the given polynomials
are in the subdomain `A`, so are those
of the subresultant polynomials. This means that the subresultant
sequence is comparable to the primitive remainder sequence without
relying on unique factorization in `A`.

To see how subresultants are associated with remainder sequences
recall that all polynomials `h` in the sequence are linear combinations of
the given polynomials `f` and `g`

.. math::

   h = uf+vg

with polynomials `u` and `v` in `K[X]`. Moreover, as is seen from the
extended Euclidean algorithm, the degrees of `u` and `v` are relatively
low, with limited growth from step to step.

Let `n = \deg(f)`, and `m = \deg(g)`, and assume `n\ge m`.
If `\deg(h) = j < m`, the coefficients of the powers `X^k` (`k > j`)
in the products `uf` and `vg` cancel each other. In particular, the
products must have the same degree, say, `l`.
Then `\deg(u) = l - n` and `\deg(v) = l - m` with a total of `2l -n - m + 2`
coefficients to be determined.

On the other hand, the equality `h = uf + vg` implies that `l - j`
linear combinations of the coefficients are zero, those associated with
the powers `X^i` (`j < i \leq l`), and one has a given non-zero value,
namely the leading coefficient of `h`.

To satisfy these `l - j + 1` linear equations the total number of
coefficients to be determined cannot be lower than `l - j + 1`, in general.
This leads to the inequality `l \ge n + m - j - 1`.
Taking `l = n + m - j - 1`, we obtain `\deg(u) = m - j - 1` and
`\deg(v) = n - j - 1`.

In the case `j = 0` the matrix of the resulting system of linear equations
is the `Sylvester matrix <https://en.wikipedia.org/wiki/Sylvester_matrix>`_
`S(f,g)` associated to `f` and `g`,
an `(n+m)\times (n+m)` matrix with coefficients of `f` and `g` as entries.
Its determinant is the `resultant <https://en.wikipedia.org/wiki/Resultant>`_
`\mathrm{res}(f,g)` of the pair `(f,g)`.
It is non-zero if and only if `f` and `g` are relatively prime.

For any `j` in the interval from `0` to `m` the matrix of the linear system is
an `(n+m-2j)\times (n+m-2j)` submatrix of the Sylvester matrix.
Its determinant `s_j(f,g)`
is called the `j` th *scalar subresultant* of `f` and `g`.

If `s_j(f,g)` is not zero, the associated equation `h = uf + vg` has
a unique solution where `\deg(h) = j` and the leading coefficient
of `h` has any given value; the one with leading coefficient
`s_j(f,g)` is the `j` th *subresultant polynomial* or, briefly,
*subresultant* of the pair `(f,g)`, and denoted `S_j(f,g)`.
This choice guarantees that the remainining coefficients
are also certain subdeterminants of the Sylvester matrix.
In particular, if `f` and `g` are in `A[X]`, so is `S_j(f,g)` as well.
This construction of subresultants applies to any `j` between
`0` and `m` regardless of the value of `s_j(f,g)`; if it is zero, then
`\deg(S_j(f,g)) < j`.

The properties of subresultants are as follows. Let `n_0 = \deg(f)`,
`n_1 = \deg(g)`, `n_2, \ldots, n_k` be the decreasing sequence of
degrees of polynomials in a remainder sequence.
Let `0 \le j \le n_1`; then

- `s_j(f,g)\ne 0` if and only if `j = n_i` for some `i`.

- `S_j(f,g)\ne 0` if and only if `j = n_i` or `j = n_i - 1` for some `i`.

Normally, `n_{i-1} - n_i = 1` for `1 < i \le k`. If `n_{i-1} - n_i > 1`
for some `i` (the *abnormal* case), then `S_{n_{i-1}-1}(f,g)` and
`S_{n_i}(f,g)` are constant multiples of each other.
Hence either one could be included in the polynomial remainder sequence.
The former is given by smaller determinants,
so it is expected to have smaller coefficients.

Collins defined the *subresultant remainder sequence* by setting

.. math::

   f_i = S_{n_{i-1}-1}(f,g) \qquad (2\le i \le k).

In the normal case, these are the same as the `S_{n_i}(f,g)`. He also
derived expressions for the constants `\gamma_i` in the remainder
formulas

.. math::

   \gamma_i f_i = \mathrm{rem}(f_{i-2},f_{i-1})

in terms of the leading coefficients of `f_1,\ldots,f_{i-1}`, working
in the field `K`.

Brown and Traub [BrownTraub71] later developed a recursive procedure
for computing the coefficients `\gamma_i`. Their algorithm deals with elements
of the domain `A` exclusively (assuming `f,g\in A[X]`). However, in the
abnormal case there was a problem, a division in `A`
which could only be conjectured to be exact.

This was subsequently justified by Brown [Brown78] who showed that
the result of the division is, in fact, a scalar subresultant.
More specifically, the constant appearing in the computation of `f_i` is
`s_{n_{i-2}}(f,g)` (Theorem 3).
The implication of this discovery is that the scalar subresultants
are computed as by-products of the algorithm, all but `s_{n_k}(f,g)`
which is not needed after finding `f_{k+1} = 0`.
Completing the last step we obtain all non-zero scalar subresultants,
including the last one which is the resultant if this does not vanish.

.. autofunction:: dmp_inner_subresultants
.. autofunction:: dmp_subresultants
.. autofunction:: dmp_prs_resultant
.. autofunction:: dmp_zz_modular_resultant
.. autofunction:: dmp_zz_collins_resultant
.. autofunction:: dmp_qq_collins_resultant
.. autofunction:: dmp_resultant
.. autofunction:: dmp_discriminant
.. autofunction:: dmp_rr_prs_gcd
.. autofunction:: dmp_ff_prs_gcd
.. autofunction:: dmp_zz_heu_gcd
.. autofunction:: dmp_qq_heu_gcd
.. autofunction:: dmp_inner_gcd
.. autofunction:: dmp_gcd
.. autofunction:: dmp_lcm
.. autofunction:: dmp_content
.. autofunction:: dmp_primitive
.. autofunction:: dmp_cancel

Polynomial factorization in characteristic zero:

.. currentmodule:: sympy.polys.factortools

.. autofunction:: dmp_trial_division
.. autofunction:: dmp_zz_mignotte_bound
.. autofunction:: dup_zz_hensel_step
.. autofunction:: dup_zz_hensel_lift
.. autofunction:: dup_zz_zassenhaus
.. autofunction:: dup_zz_irreducible_p
.. autofunction:: dup_cyclotomic_p
.. autofunction:: dup_zz_cyclotomic_poly
.. autofunction:: dup_zz_cyclotomic_factor
.. autofunction:: dup_zz_factor_sqf
.. autofunction:: dup_zz_factor
.. autofunction:: dmp_zz_wang_non_divisors
.. autofunction:: dmp_zz_wang_test_points
.. autofunction:: dmp_zz_wang_lead_coeffs
.. autofunction:: dmp_zz_diophantine
.. autofunction:: dmp_zz_wang_hensel_lifting
.. autofunction:: dmp_zz_wang
.. autofunction:: dmp_zz_factor
.. autofunction:: dmp_ext_factor
.. autofunction:: dup_gf_factor
.. autofunction:: dmp_factor_list
.. autofunction:: dmp_factor_list_include
.. autofunction:: dmp_irreducible_p

Groebner basis algorithms
*************************

Groebner bases can be used to answer many problems in computational
commutative algebra. Their computation in rather complicated, and very
performance-sensitive. We present here various low-level implementations of
Groebner basis computation algorithms; please see the previous section of the
manual for usage.

.. currentmodule:: sympy.polys.groebnertools

.. autofunction:: groebner
.. autofunction:: spoly
.. autofunction:: red_groebner
.. autofunction:: is_groebner
.. autofunction:: is_minimal
.. autofunction:: is_reduced

.. currentmodule:: sympy.polys.fglmtools

.. autofunction:: matrix_fglm

Groebner basis algorithms for modules are also provided:

.. currentmodule:: sympy.polys.distributedmodules

.. autofunction:: sdm_spoly
.. autofunction:: sdm_ecart
.. autofunction:: sdm_nf_mora
.. autofunction:: sdm_groebner

Options
=======

.. automodule:: sympy.polys.polyoptions

.. autoclass:: sympy.polys.polyoptions.Options
.. autofunction:: sympy.polys.polyoptions.build_options

Configuration
=============

.. automodule:: sympy.polys.polyconfig

.. autofunction:: sympy.polys.polyconfig.setup

Exceptions
==========

These are exceptions defined by the polynomials module.

TODO sort and explain

.. currentmodule:: sympy.polys.polyerrors

.. autoclass:: BasePolynomialError

.. autoclass:: ExactQuotientFailed
.. autoclass:: OperationNotSupported
.. autoclass:: HeuristicGCDFailed
.. autoclass:: HomomorphismFailed
.. autoclass:: IsomorphismFailed
.. autoclass:: ExtraneousFactors
.. autoclass:: EvaluationFailed
.. autoclass:: RefinementFailed
.. autoclass:: CoercionFailed
.. autoclass:: NotInvertible
.. autoclass:: NotReversible
.. autoclass:: NotAlgebraic
.. autoclass:: DomainError
.. autoclass:: PolynomialError
.. autoclass:: UnificationFailed
.. autoclass:: GeneratorsNeeded
.. autoclass:: ComputationFailed
.. autoclass:: GeneratorsError
.. autoclass:: UnivariatePolynomialError
.. autoclass:: MultivariatePolynomialError
.. autoclass:: PolificationFailed
.. autoclass:: OptionError
.. autoclass:: FlagError

Reference
=========

Modular GCD
***********

.. currentmodule:: sympy.polys.modulargcd

.. autofunction:: modgcd_univariate
.. autofunction:: modgcd_bivariate
.. autofunction:: modgcd_multivariate
.. autofunction:: _modgcd_multivariate_p
.. autofunction:: func_field_modgcd

Undocumented
============

Many parts of the polys module are still undocumented, and even where there is
documentation it is scarce. Please contribute!
