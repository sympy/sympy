.. _polys-domains:

======================================================
Understanding and using the Domains of the poly module
======================================================

This page introduces the idea of the "domains" that are used in SymPy's
``polys`` module. The emphasis is on introducing how to use the domains
directly and on understanding how they are used internally as part of the
:py:class:`~.Poly` class. This is a relatively advanced topic so for a more
introductory understanding of the :py:class:`~.Poly` class and the ``polys``
module it is recommended to read :ref:`polys-basics` instead.

What are the domains?
=====================

For most users the domains are only really noticable in the printed output of
a :py:class:`~.Poly`::

  >>> from sympy import Symbol, Poly
  >>> x = Symbol('x')
  >>> Poly(x**2 + x)
  Poly(x**2 + x, x, domain='ZZ')
  >>> Poly(x**2 + x/2)
  Poly(x**2 + 1/2*x, x, domain='QQ')

We see here that one :py:class:`~.Poly` has domain ``ZZ`` representing the
integers and the other has domain ``QQ`` representing the rationals. These
indicate the "domain" from which the coefficients of the polynomial are drawn.

From a high-level the domains represent formal concepts such as the set of
integers `\mathbb{Z}` or rationals `\mathbb{Q}`. The word "domain" here is a
reference to the mathematical concept of an `integral domain`_.

.. _integral domain: https://en.wikipedia.org/wiki/Integral_domain

Internally the domains correspond to different computational implementations
and representations of the expressions that the polynomials correspond to.
The :py:class:`~.Poly` object itself has an internal representation as a
``list`` of coefficients and a ``domain`` attribute representing the
implementation of those coefficients::

  >>> p = Poly(x**2 + x/2)
  >>> p
  Poly(x**2 + 1/2*x, x, domain='QQ')
  >>> p.domain
  QQ
  >>> p.rep
  DMP([1, 1/2, 0], QQ, None)
  >>> p.rep.rep
  [1, 1/2, 0]
  >>> type(p.rep.rep[0])
  <class 'sympy.polys.domains.pythonrational.PythonRational'>

Here the domain is ``QQ`` which represents the implementation of the rational
numbers in the domain system. The :py:class:`~.Poly` instance itself has a
``domain`` attribute ``QQ`` and then a list of ``PythonRational`` coefficients
where ``PythonRational`` is the class that implements the elements of the
``QQ`` domain. The list of coefficients ``[1, 1/2, 0]`` gives a standardised
low-level representation of the polynomial expression ``(1)*x**2 + (1/2)*x +
(0)``.

This page looks at the different domains that are defined in SymPy, how they
are implemented and how they can be used. It introduces how to use the domains
and domain elements directly and explains how they are used internally as part
of :py:class:`~.Poly` objects. This information is more relevant for
development in SymPy than it is for users of the ``polys`` module.

Representing expressions symbolically
=====================================

There are many different ways that a mathematical expression can be
represented symbolically. The purpose of the polynomial domains is to provide
suitable implementations for different classes of expressions. This section
considers the basic approaches to the symbolic representation of mathematical
expressions: "tree", "dense polynomial"  and "sparse polynomial".

The most general representation is as a `tree`_ and this is the representation
used for most ordinary SymPy expressions. We can see this representation using
the :py:func:`~.srepr` function::

  >>> from sympy import Symbol, srepr
  >>> x = Symbol('x')
  >>> e = 1 + 1/(2 + x**2)
  >>> e
  1 + 1/(x**2 + 2)
  >>> print(srepr(e))
  Add(Integer(1), Pow(Add(Pow(Symbol('x'), Integer(2)), Integer(2)), Integer(-1)))

Here the expression ``e`` is represented as an :py:class:`~.Add` node which
has two children ``1`` and ``1/(x**2 + 2)``. The child ``1`` is represented as
``Integer(1)`` and the other child is represented as a :py:class:`~.Pow` with
base ``x**2 + 2`` and exponent ``1``. Then ``x**2 + 2`` is represented as an
:py:class:`~.Add` with children ``x**2`` and ``2`` and so on. In this way the
expression is represented as a tree where the internal nodes are operations
like :py:class:`~.Add`, :py:class:`~.Mul`, :py:class:`~.Pow` and so on and the
leaf nodes are atomic expressions like ``Integer(1)`` or ``Symbol('x')``. See
:ref:`tutorial-manipulation` for more about this representation.

The tree representation is core to the architecture of :py:class:`~.Basic` in
SymPy. It is a highly flexible representation that can represent a very wide
range of possible expressions. It can also represent equivalent expressions in
different ways e.g.::

  >>> e = x*(x + 1)
  >>> e
  x*(x + 1)
  >>> e.expand()
  x**2 + x

These two expression although equivalent have different tree representations::

  >>> print(srepr(e))
  Mul(Symbol('x'), Add(Symbol('x'), Integer(1)))
  >>> print(srepr(e.expand()))
  Add(Pow(Symbol('x'), Integer(2)), Symbol('x'))

Being able to represent the same expression in different ways is both a
strength and a weakness. It is useful to be able to convert an
expression in to different forms for different tasks but having non-unique
representations makes it hard to tell when two expressions are equivalent
which is in fact very important for many computational algorithms. The most
important task is being able to tell when an expression is equal to zero which
is undecidable in general (`Richardon's theorem`_) but is decidable in many
important special cases.

Restricting the set of allowed expressions to special cases allows for much
more efficient symbolic representations. As we already saw :py:class:`~.Poly`
can represent a polynomial as a list of coefficients. This means that an
expression like ``x**4 + x + 1`` could be represented simply as ``[1, 0, 0, 1,
1]``. This list of coefficients representation of a polynomial expression is
known as the "dense univariate polynomial" (DUP) representation. Working
within that representation algorithms for multiplication, addition and
crucially zero-testing can be much more efficient than with the corresponding
tree representations. We can see this representation from a :py:class:`~.Poly`
instance by looking it its ``rep.rep`` attribute::

  >>> p = Poly(x**4 + x + 1)
  >>> p.rep.rep
  [1, 0, 0, 1, 1]

In the DUP representation it is not possible to represent the same expression
in different ways. There is no distinction between ``x*(x + 1)`` and ``x**2 +
x`` because both are just ``[1, 1, 0]``. This means that comparing two
expressions is easy: they are equal if and only if all of their coefficients
are equal. Zero-testing is particularly easy: the polynomial is zero if and
only if all coefficients are zero (of course we need to have easy zero-testing
for the coefficients themselves).

We can make functions that operate on the DUP representation much more
efficiently than functions that operate on the tree representation. Many
operations with standard sympy expressions are in fact computed by converting
to a polynomial representation and then performing the calculation. An example
is the :py:func:`~.factor` function::

  >>> from sympy import factor
  >>> e = 2*x**3 + 10*x**2 + 16*x + 8
  >>> e
  2*x**3 + 10*x**2 + 16*x + 8
  >>> factor(e)
  2*(x + 1)*(x + 2)**2

Internally :py:func:`~.factor` will convert the expression from the tree
representation into the DUP representation and then use the function
``dup_factor_list``::

  >>> from sympy import ZZ
  >>> from sympy.polys.factortools import dup_factor_list
  >>> p = [ZZ(2), ZZ(10), ZZ(16), ZZ(8)]
  >>> p
  [2, 10, 16, 8]
  >>> dup_factor_list(p, ZZ)
  (2, [([1, 1], 1), ([1, 2], 2)])

There are many more examples of functions with ``dup_*`` names for operating
on the DUP representation that are documented in :ref:`polys-internals`. There
are also functions with the ``dmp_*`` prefix for operating on multivariate
polynomials.

A multivariate polynomial (a polynomial in multiple variables) can be
represented as a polynomial with coefficients that are themselves polynomials.
For example ``x**2*y + x**2 + x*y + y + 1`` can be represented as polynomial
in ``x`` where the coefficients are themselves polynomials in ``y`` i.e.: ``(y
+ 1)*x**2 + (y)*x + (y+1)``. Since we can represent a polynomial with a list
of coefficients a multivariate polynomial can be represented with a list of
lists of coefficients::

  >>> from sympy import symbols
  >>> x, y = symbols('x, y')
  >>> p = Poly(x**2*y + x**2 + x*y + y + 1)
  >>> p
  Poly(x**2*y + x**2 + x*y + y + 1, x, y, domain='ZZ')
  >>> p.rep.rep
  [[1, 1], [1, 0], [1, 1]]

This list of lists of (lists of...) coefficients representation is known as
the "dense multivariate polynomial" (DMP) representation. Instead of lists we
can use a dict mapping nonzero monomial terms to their coefficients. This is
known as the "sparse polynomial" representation. We can see what this would
look like using the :py:meth:`~.Poly.as_dict` method::

  >>> Poly(7*x**20 + 8*x + 9).rep.rep
  [7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 9]
  >>> Poly(7*x**20 + 8*x + 9).as_dict()
  {(0,): 9, (1,): 8, (20,): 7}

The keys of this dict are the exponents of the powers of ``x`` and the values
are the coefficients so e.g. ``7*x**20`` becomes ``(20,): 7`` in the dict. The
key is a tuple so that in the multivariate case something like ``4*x**2*y**3``
can be represented as ``(2, 3): 4``. The sparse representation can be more
efficient as it avoids the need to store and manipulate the zero coefficients.
With a large number of generators (variables) the dense representation becomes
particularly inefficient and it is better to use the sparse representation::

  >>> from sympy import prod
  >>> gens = symbols('x:10')
  >>> gens
  (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9)
  >>> p = Poly(prod(gens))
  >>> p
  Poly(x0*x1*x2*x3*x4*x5*x6*x7*x8*x9, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, domain='ZZ')
  >>> p.rep.rep
  [[[[[[[[[[1, 0], []], [[]]], [[[]]]], [[[[]]]]], [[[[[]]]]]], [[[[[[]]]]]]], [[[[[[[]]]]]]]], [[[[[[[[]]]]]]]]], [[[[[[[[[]]]]]]]]]]
  >>> p.as_dict()
  {(1, 1, 1, 1, 1, 1, 1, 1, 1, 1): 1}

The dict representation shown in the last output maps from the monomial which
is represented as a tuple of powers (``(1, 1, 1, ...)`` i.e. ``x0**1 * x1**1,
...``) to the coefficient ``1``. Compared to the DMP representation we have a
much more flattened data structure: it is a ``dict`` with only one key and
value. Algorithms for working with sparse representations would likely be
much more efficient than dense algorithms for this particular example
polynomial.

SymPy's polynomial module has implementations of polynomial expressions based
on both the dense and sparse representations. There are also other
implementations of different special classes of expressions that can be used
as the coefficients of those polynomials. The rest of this page discusses what
those representations are and how to use them.

.. _tree: https://en.wikipedia.org/wiki/Tree_(data_structure)
.. _Richardon's theorem: https://en.wikipedia.org/wiki/Richardson%27s_theorem

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
``+``, ``-``, ``*`` and ``**`` (with nonnegative integer exponent) should be
assumed to work with arbitrary domain elements. All other operations should be
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
:py:meth:`~.Domain.from_sympy` for converting back and forth between sympy
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

We have already seen that in some cases we can use the domain object itself as
a constructor e.g. ``QQ(2)``. This will generally work provided the arguments
given are valid for the ``dtype`` of the domain. Although it is convenient to
use this in interactive sessions and in demonstrations it is generally better
to use the :py:meth:`~.Domain.from_sympy` method for constructing domain elements
from sympy expressions (or from objects that can be sympified to sympy
expressions).

It is important not to mix domain elements with other Python types such as
``int``, ``float``, as well as standard sympy :py:class:`~.Basic` expressions.
When working in a domain, care should be taken as some Python operations will
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

Gaussian integers and Gaussian rationals
========================================

The two example domains that we have seen so far are ``ZZ`` and ``QQ``
representing the integers and the rationals respectively. There are other
simple domains such as ``ZZ_I`` and ``QQ_I`` representing the
`Gaussian integers`_ and `Gaussian rationals`_. The Gaussian integers are
numbers of the form `a\sqrt{-1} + b` where `a` and `b` are integers. The
Gaussian rationals are defined similarly except that `a` and `b` can be
rationals. We can use the Gaussian domains like::

  >>> from sympy import ZZ_I, QQ_I, I
  >>> z = ZZ_I.from_sympy(1 + 2*I)
  >>> z
  (1 + 2*I)
  >>> z**2
  (-3 + 4*I)

Note the constrast with the way this calculation works in the tree
representation where :py:func:`~.expand` is needed to get the reduced form::

  >>> from sympy import expand, I
  >>> z = 1 + 2*I
  >>> z**2
  (1 + 2*I)**2
  >>> expand(z**2)
  -3 + 4*I

The ``ZZ_I`` and ``QQ_I`` domains are implemented by the classes
``GaussianIntegerRing`` and ``GaussianRationalField`` and their elements by
``GaussianInteger`` and ``GaussianRational`` respectively. The internal
representation for an element of ``ZZ_I`` or ``QQ_I`` is simply as a pair
``(a, b)`` of elements of ``ZZ`` or ``QQ`` respectively. The domain ``ZZ_I``
is a ring with similar properties to ``ZZ`` whereas ``QQ_I`` is a field
much like ``QQ``::

  >>> ZZ.is_Field
  False
  >>> QQ.is_Field
  True
  >>> ZZ_I.is_Field
  False
  >>> QQ_I.is_Field
  True

Since ``QQ_I`` is a field division by nonzero elements is always possible
whereas in ``ZZ_I`` we have the important concept of the greatest common
divisor (GCD)::

  >>> e1 = QQ_I.from_sympy(1+I)
  >>> e2 = QQ_I.from_sympy(2-I/2)
  >>> e1/e2
  (6/17 + 10/17*I)
  >>> ZZ_I.gcd(ZZ_I(5), ZZ_I.from_sympy(1+2*I))
  (1 + 2*I)

.. _Gaussian integers: https://en.wikipedia.org/wiki/Gaussian_integer
.. _Gaussian rationals: https://en.wikipedia.org/wiki/Gaussian_rational

Finite fields
=============

So far we have seen the domains ``ZZ``, ``QQ``, ``ZZ_I``, and ``QQ_I``. There
are also domains representing the `Finite fields`_ although the implementation
of these is incomplete. A finite field of *prime* order can be constructed
with ``FF`` or ``GF``. A domain for the finite field or prime order `p` can be
constructed with ``FF(p)``::

  >>> from sympy import FF
  >>> K = FF(5)
  >>> two = K(2)
  >>> two
  2 mod 5
  >>> two ** 2
  4 mod 5
  >>> two ** 3
  3 mod 5

Finite fields of order `p^n` where `n \ne 1` are not implemented. It is
possible to use e.g. ``FF(6)`` or ``FF(9)`` but the resulting domain is *not*
a field. It would be good to have a proper implementation of prime-power order
finite fields but this is not yet available in SymPy (contributions welcome!).

.. _Finite fields: https://en.wikipedia.org/wiki/Finite_field

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
it possible to represent equivalent expressions in different ways e.g.::

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

It is also possible to construct nested polynomial rings (although it is less
efficient). The ring ``K[x][y]`` is formally equivalent to ``K[x,y]`` although
their implementations in sympy are different::

  >>> K = ZZ[x][y]
  >>> p = K(x**2 + x*y + y**2)
  >>> p
  y**2 + x*y + x**2
  >>> dict(p)
  {(0,): x**2, (1,): x, (2,): 1}

Here the coefficients like ``x**2`` are instances of :py:class:`~.PolyElement`
as well so this is a ``dict`` where the values are also dicts. The full
representation is more like::

  >>> {k: dict(v) for k, v in p.items()}
  {(0,): {(2,): 1}, (1,): {(1,): 1}, (2,): {(0,): 1}}

The multivariate ring domain ``ZZ[x,y]`` has a more efficient representation
as a single flattened ``dict``::

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
  >>> b = K.one + K.from_sympy(sqrt(2))
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
  >>> sqrt2 = K.from_sympy(sqrt(2))
  >>> sqrt3 = K.from_sympy(sqrt(3))
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

Some domains are classified as fields and others are not. The principal
difference between a field and a non-field domain is that in a field it is
always possible to divide any element by any nonzero element. It is usually
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

Choosing a domain
=================

In the workflow described above the idea is to start with some sympy
expressions, choose a domain and convert all the expressions into that domain
in order to perform some calculation. The obvious question that arises is how
to choose an appropriate domain to represent some sympy expressions. For this
there is a function :py:func:`~.construct_domain` which takes a list of
expressions and will choose a domain and convert all of the expressions to
that domain::

  >>> from sympy import construct_domain, Integer
  >>> elements_sympy = [Integer(3), Integer(2)]  # elements as Basic instances
  >>> elements_sympy
  [3, 2]
  >>> K, elements_K = construct_domain(elements_sympy)
  >>> K
  ZZ
  >>> elements_K
  [3, 2]
  >>> type(elements_sympy[0])
  <class 'sympy.core.numbers.Integer'>
  >>> type(elements_K[0])
  <class 'int'>

In this example we see that the two integers ``3`` and ``2`` can be
represented in the domain ``ZZ``. The expressions have been converted to
elements of that domain which in this case means the ``int`` type rather than
instances of :py:class:`~.Basic`. It is not necessary to explicitly create
:py:class:`~.Basic` instances when the inputs can be sympified so e.g.
``construct_domain([3, 2])`` would give the same output as above.

Given more complicated inputs :py:func:`~.construct_domain` will choose more
complicated domains::

  >>> from sympy import Rational, symbols
  >>> x, y = symbols('x, y')
  >>> construct_domain([Rational(1, 2), Integer(3)])[0]
  QQ
  >>> construct_domain([2*x, 3])[0]
  ZZ[x]
  >>> construct_domain([x/2, 3])[0]
  QQ[x]
  >>> construct_domain([2/x, 3])[0]
  ZZ(x)
  >>> construct_domain([x, y])[0]
  ZZ[x,y]

If any noninteger rational numbers are found in the inputs then the ground
domain will be ``QQ`` rather than ``ZZ``. If any symbol is found in the
inputs then a :py:class:`~.PolynomialRing` will be created. A multivariate
polynomial ring such as ``QQ[x,y]`` can also be created if there are multiple
symbols in the inputs. If any symbols appear in the denominators then a
:py:class:`~.RationalField` like ``QQ(x)`` will be created instead (and the
ground domain ``ZZ`` will be promoted to the field ``QQ`` even if all
coeficients are integers).

Some of the domains above are fields and others are (non-field) rings. In some
contexts it is necessary to have a field domain so that division is possible
and for this :py:func:`~.construct_domain` has an option ``field=True`` which
will force the construction of a field domain even if the expressions can all
be represented in a non-field ring::

  >>> construct_domain([1, 2], field=True)[0]
  QQ
  >>> construct_domain([2*x, 3], field=True)[0]
  ZZ(x)
  >>> construct_domain([x/2, 3], field=True)[0]
  ZZ(x)
  >>> construct_domain([2/x, 3], field=True)[0]
  ZZ(x)
  >>> construct_domain([x, y], field=True)[0]
  ZZ(x,y)

By default :py:func:`~.construct_domain` will not construct an algebraic
extension field and will instead use the ``EX`` domain
(:py:class:`~.ExpressionDomain`). The keyword argument ``extension=True`` can
be used to construct an :py:class:`~.AlgebraicField` if the inputs are
irrational but algebraic::

  >>> from sympy import sqrt
  >>> construct_domain([sqrt(2)])[0]
  EX
  >>> construct_domain([sqrt(2)], extension=True)[0]
  QQ<sqrt(2)>
  >>> construct_domain([sqrt(2), sqrt(3)], extension=True)[0]
  QQ<sqrt(2) + sqrt(3)>

When there are algebraically independent transcendentals in the inputs a
:py:class:`~.PolynomialRing` or :py:class:`~.RationalField` will be
constructed treating those transcendentals as generators::

  >>> from sympy import sin, cos
  >>> construct_domain([sin(x), y])[0]
  ZZ[y,sin(x)]

However if there is a possibility that the inputs are not algebraically
independent then the domain will be ``EX``::

  >>> construct_domain([sin(x), cos(x)])[0]
  EX

Converting elements between different domains
=============================================

It is often useful to combine calculations performed over different domains.
However just as it is important to avoid mixing domain elements with normal
sympy expressions and other Python types it is also important to avoid mixing
elements from different domains. The :py:meth:`~.Domain.convert_from` method
is used to convert elements from one domain into elements of another domain::

  >>> num_zz = ZZ(3)
  >>> ZZ.of_type(num_zz)
  True
  >>> num_qq = QQ.convert_from(num_zz, ZZ)
  >>> ZZ.of_type(num_qq)
  False
  >>> QQ.of_type(num_qq)
  True

The :py:meth:`~.Domain.convert` method can be called without specifying the
source domain as the second ergument e.g.::

  >>> QQ.convert(ZZ(2))
  2

This works because :py:meth:`~.Domain.convert` can check the type of ``ZZ(2)``
and can try to work out what domain (``ZZ``) it is an element of. Certain
domains like ``ZZ`` and ``QQ`` are treated as special cases to make this work.
Elements of more complicated domains are instances of subclasses of
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

It is more efficient though to call :py:meth:`~.Domain.convert_from` with the
source domain specified as the second argument::

  >>> QQ.convert_from(ZZ(2), ZZ)
  2

Unifying domains
================

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
  >>> x3 = K3.convert_from(x1, K1)
  >>> y3 = K3.convert_from(y2, K2)
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

The :py:meth:`~.Domain.unify` method understands how to combine different
polynomial ring domains and how to unify the base domain::

  >>> ZZ[x].unify(ZZ[y])
  ZZ[x,y]
  >>> ZZ[x,y].unify(ZZ[y])
  ZZ[x,y]
  >>> ZZ[x].unify(QQ)
  QQ[x]

It is also possible to unify algebraic fields and rational function fields as
well::

  >>> K1 = QQ.algebraic_field(sqrt(2))[x]
  >>> K2 = QQ.algebraic_field(sqrt(3))[y]
  >>> K1
  QQ<sqrt(2)>[x]
  >>> K2
  QQ<sqrt(3)>[y]
  >>> K1.unify(K2)
  QQ<sqrt(2) + sqrt(3)>[x,y]
  >>> QQ.frac_field(x).unify(ZZ[y])
  ZZ(x,y)
