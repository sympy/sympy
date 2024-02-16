.. _assumptions-guide:

=============
 Assumptions
=============

This page outlines the core assumptions system in SymPy. It explains what the
core assumptions system is, how the assumptions system is used and what the
different assumptions predicates mean.

.. note:: This page describes the core assumptions system also often referred
   to as the "old assumptions" system. There is also a "new assumptions"
   system which is described elsewhere. Note that the system described here is
   actually the system that is widely used in SymPy. The "new assumptions"
   system is not really used anywhere in SymPy yet and the "old assumptions"
   system will not be removed. At the time of writing (SymPy 1.7) it is
   still recommended for users to use the old assumption system.

Firstly we consider what happens when taking the square root of the square of
a concrete integer such as $2$ or $-2$:

    >>> from sympy import sqrt
    >>> sqrt(2**2)
    2
    >>> sqrt((-2)**2)
    2
    >>> x = 2
    >>> sqrt(x**2)
    2
    >>> sqrt(x**2) == x
    True
    >>> y = -2
    >>> sqrt(y**2) == y
    False
    >>> sqrt(y**2) == -y
    True

What these examples demonstrate is that for a positive number $x$ we have
$\sqrt{x^2} = x$ whereas for a negative number we would instead have
$\sqrt{x^2} = -x$. That may seem obvious but the situation can be more
surprising when working with a symbol rather then an explicit number. For
example

    >>> from sympy import Symbol, simplify
    >>> x = Symbol('x')
    >>> sqrt(x**2)
    sqrt(x**2)

It might look as if that should simplify to ``x`` but it does not even if
:func:`~.simplify` is used:

    >>> simplify(sqrt(x**2))
    sqrt(x**2)

This is because SymPy will refuse to simplify this expression if the
simplification is not valid for *every* possible value of ``x``. By default
the symbol ``x`` is considered only to represent something roughly like an
arbitrary complex number and the obvious simplification here is only valid for
positive real numbers. Since ``x`` is not known to be positive or even real no
simplification of this expression is possible.

We can tell SymPy that a symbol represents a positive real number when
creating the symbol and then the simplification will happen automatically:

    >>> y = Symbol('y', positive=True)
    >>> sqrt(y**2)
    y

This is what is meant by "assumptions" in SymPy. If the symbol ``y``
is created with ``positive=True`` then SymPy will *assume* that it represents
a positive real number rather than an arbitrary complex or possibly infinite
number. That *assumption* can make it possible to simplify expressions or
might allow other manipulations to work. It is usually a good idea to be as
precise as possible about the assumptions on a symbol when creating it.


The (old) assumptions system
============================

There are two sides to the assumptions system. The first side is that we can
declare assumptions on a symbol when creating the symbol. The other side is
that we can query the assumptions on any expression using the corresponding
``is_*`` attribute. For example:

    >>> x = Symbol('x', positive=True)
    >>> x.is_positive
    True

We can query assumptions on any expression not just a symbol:

    >>> x = Symbol('x', positive=True)
    >>> expr = 1 + x**2
    >>> expr
    x**2 + 1
    >>> expr.is_positive
    True
    >>> expr.is_negative
    False

The values given in an assumptions query use three-valued "fuzzy" logic. Any
query can return ``True``, ``False``, or ``None`` where ``None`` should be
interpreted as meaning that the result is *unknown*.

    >>> x = Symbol('x')
    >>> y = Symbol('y', positive=True)
    >>> z = Symbol('z', negative=True)
    >>> print(x.is_positive)
    None
    >>> print(y.is_positive)
    True
    >>> print(z.is_positive)
    False

.. note:: We need to use ``print`` in the above examples because the special
          value ``None`` does not display by default in the Python
          interpreter.

There are several reasons why an assumptions query might give ``None``. It is
possible that the query is *unknowable* as in the case of ``x`` above. Since
``x`` does not have any assumptions declared it roughly represents an
arbitrary complex number. An arbitrary complex number *might* be a positive
real number but it also might *not* be. Without further information there is
no way to resolve the query ``x.is_positive``.

Another reason why an assumptions query might give ``None`` is that there does
in many cases the problem of determining whether an expression is e.g.
positive is *undecidable*. That means that there does not exist an algorithm
for answering the query in general. For some cases an algorithm or at least a
simple check would be possible but has not yet been implemented although it
could be added to SymPy.

The final reason that an assumptions query might give ``None`` is just that
the assumptions system does not try very hard to answer complicated queries.
The system is intended to be fast and uses simple heuristic methods to
conclude a ``True`` or ``False`` answer in common cases. For example any sum
of positive terms is positive so:

    >>> from sympy import symbols
    >>> x, y = symbols('x, y', positive=True)
    >>> expr = x + y
    >>> expr
    x + y
    >>> expr.is_positive
    True

The last example is particularly simple so the assumptions system is able to
give a definite answer. If the sum involved a mix of positive or negative
terms it would be a harder query:

    >>> x = Symbol('x', real=True)
    >>> expr = 1 + (x - 2)**2
    >>> expr
    (x - 2)**2 + 1
    >>> expr.is_positive
    True
    >>> expr2 = expr.expand()
    >>> expr2
    x**2 - 4*x + 5
    >>> print(expr2.is_positive)
    None

Ideally that last example would give ``True`` rather than ``None`` because the
expression is always positive for any real value of ``x`` (and ``x`` has been
assumed real). The assumptions system is intended to be efficient though: it
is expected many more complex queries will not be fully resolved. This is
because assumptions queries are primarily used internally by SymPy as part of
low-level calculations. Making the system more comprehensive would slow SymPy
down.

Note that in fuzzy logic giving an indeterminate result ``None`` is never a
contradiction. If it is possible to infer a definite ``True`` or ``False``
result when resolving a query then that is better than returning ``None``.
However a result of ``None`` is not a *bug*. Any code that uses the
assumptions system needs to be prepared to handle all three cases for any
query and should not presume that a definite answer will always be given.

The assumptions system is not just for symbols or for complex expressions. It
can also be used for plain SymPy integers and other objects. The assumptions
predicates are available on any instance of :class:`~.Basic` which is the superclass
for most classes of SymPy objects. A plain Python :class:`int` is not a
:class:`~.Basic` instance and can not be used to query assumptions
predicates. We can "sympify" regular Python objects to become SymPy objects
with :func:`~.sympify` or ``S`` (:class:`~.SingletonRegistry`) and then the
assumptions system can be used:

    >>> from sympy import S
    >>> x = 2
    >>> x.is_positive
    Traceback (most recent call last):
    ...
    AttributeError: 'int' object has no attribute 'is_positive'
    >>> x = S(2)
    >>> type(x)
    <class 'sympy.core.numbers.Integer'>
    >>> x.is_positive
    True


Gotcha: symbols with different assumptions
==========================================

In SymPy it is possible to declare two symbols with different names and they
will implicitly be considered equal under *structural equality*:

    >>> x1 = Symbol('x')
    >>> x2 = Symbol('x')
    >>> x1
    x
    >>> x2
    x
    >>> x1 == x2
    True

However if the symbols have different assumptions then they will be considered
to represent distinct symbols:

    >>> x1 = Symbol('x', positive=True)
    >>> x2 = Symbol('x')
    >>> x1
    x
    >>> x2
    x
    >>> x1 == x2
    False

One way to simplify an expression is to use the :func:`~.posify` function
which will replace all symbols in an expression with symbols that have the
assumption ``positive=True`` (unless that contradicts any existing assumptions
for the symbol):

    >>> from sympy import posify, exp
    >>> x = Symbol('x')
    >>> expr = exp(sqrt(x**2))
    >>> expr
    exp(sqrt(x**2))
    >>> posify(expr)
    (exp(_x), {_x: x})
    >>> expr2, rep = posify(expr)
    >>> expr2
    exp(_x)

The :func:`~.posify` function returns the expression with all symbols replaced
(which might lead to simplifications) and also a dict which maps the new
symbols to the old that can be used with :py:meth:`~.Basic.subs`. This is
useful because otherwise the new expression with the new symbols having the
``positive=True`` assumption will not compare equal to the old:

    >>> expr2
    exp(_x)
    >>> expr2 == exp(x)
    False
    >>> expr2.subs(rep)
    exp(x)
    >>> expr2.subs(rep) == exp(x)
    True


Applying assumptions to string inputs
=====================================

We have seen how to set assumptions when  :class:`~.Symbol` or
:func:`~.symbols` explicitly. A natural question to ask is in what other
situations can we assign assumptions to an object?

It is common for users to use strings as input to SymPy functions (although
the general feeling among SymPy developers is that this should be discouraged)
e.g.:

    >>> from sympy import solve
    >>> solve('x**2 - 1')
    [-1, 1]

When creating symbols explicitly it would be possible to assign assumptions
that would affect the behaviour of :func:`~.solve`:

    >>> x = Symbol('x', positive=True)
    >>> solve(x**2 - 1)
    [1]

When using string input SymPy will create the expression and create all of the
symbolc implicitly so the question arises how can the assumptions be
specified? The answer is that rather than depending on implicit string
conversion it is better to use the :func:`~.parse_expr` function explicitly
and then it is possible to provide assumptions for the symbols e.g.:

    >>> from sympy import parse_expr
    >>> parse_expr('x**2 - 1')
    x**2 - 1
    >>> eq = parse_expr('x**2 - 1', {'x':Symbol('x', positive=True)})
    >>> solve(eq)
    [1]

.. note:: The :func:`~.solve` function is unusual as a high level API in that it
          actually checks the assumptions on any input symbols (the unknowns)
          and uses that to tailor its output. The assumptions system otherwise
          affects low-level evaluation but is not necessarily handled
          explicitly by high-level APIs.


Predicates
==========

There are many different predicates that can be assumed for a symbol or can be
queried for an expression. It is possible to combine multiple predicates when
creating a symbol. Predicates are logically combined using *and* so if a
symbol is declared with ``positive=True`` and also with ``integer=True`` then
it is both positive *and* integer:

    >>> x = Symbol('x', positive=True, integer=True)
    >>> x.is_positive
    True
    >>> x.is_integer
    True

The full set of known predicates for a symbol can be accessed using the
:attr:`~.Basic.assumptions0` attribute:

    >>> x.assumptions0
    {'algebraic': True,
     'commutative': True,
     'complex': True,
     'extended_negative': False,
     'extended_nonnegative': True,
     'extended_nonpositive': False,
     'extended_nonzero': True,
     'extended_positive': True,
     'extended_real': True,
     'finite': True,
     'hermitian': True,
     'imaginary': False,
     'infinite': False,
     'integer': True,
     'irrational': False,
     'negative': False,
     'noninteger': False,
     'nonnegative': True,
     'nonpositive': False,
     'nonzero': True,
     'positive': True,
     'rational': True,
     'real': True,
     'transcendental': False,
     'zero': False}

We can see that there are many more predicates listed than the two that were
used to create ``x``. This is because the assumptions system can infer some
predicates from combinations of other predicates. For example if a symbol is
declared with ``positive=True`` then it is possible to infer that it should
have ``negative=False`` because a positive number can never be negative.
Similarly if a symbol is created with ``integer=True`` then it is possible to
infer that is should have ``rational=True`` because every integer is a
rational number.

A full table of the possible predicates and their definitions is given below.

.. _assumptions-guide-predicates:

.. list-table:: Assumptions predicates for the (old) assumptions
    :widths: 20, 45, 35
    :header-rows: 1

    * - Predicate
      - Definition
      - Implications

    * - ``commutative``
      - A commutative expression. A ``commutative`` expression commutes with
        all other expressions under multiplication. If an expression ``a`` has
        ``commutative=True`` then ``a * b == b * a`` for any other expression
        ``b`` (even if ``b`` is not ``commutative``).  Unlike all other
        assumptions predicates ``commutative`` must always be ``True`` or
        ``False`` and can never be ``None``. Also unlike all other predicates
        ``commutative`` defaults to ``True`` in e.g.  ``Symbol('x')``.
        [commutative]_
      -

    * - ``infinite``
      - An infinite expression such as ``oo``, ``-oo`` or ``zoo``.
        [infinite]_
      - | ``== !finite``

    * - ``finite``
      - A finite expression. Any expression that is not ``infinite`` is
        considered ``finite``.
        [infinite]_
      - | ``== !infinite``

    * - ``hermitian``
      - An element of the field of Hermitian operators.
        [antihermitian]_
      -

    * - ``antihermitian``
      - An element of the field of antihermitian operators.
        [antihermitian]_
      -

    * - ``complex``
      - A complex number, $z\in\mathbb{C}$. Any number of the form $x + iy$
        where $x$ and $y$ are ``real`` and $i = \sqrt{-1}$. All ``complex``
        numbers are ``finite``. Includes all ``real`` numbers.
        [complex]_
      - | ``-> commutative``
        | ``-> finite``

    * - ``algebraic``
      - An algebraic number, $z\in\overline{\mathbb{Q}}$. Any number that is a root
        of a non-zero polynomial $p(z)\in\mathbb{Q}[z]$ having rational
        coefficients. All ``algebraic`` numbers are ``complex``. An
        ``algebraic`` number may or may not be ``real``. Includes all
        ``rational`` numbers.
        [algebraic]_
      - | ``-> complex``

    * - ``transcendental``
      - A complex number that is not algebraic,
        $z\in\mathbb{C}-\overline{\mathbb{Q}}$. All ``transcendental`` numbers are
        ``complex``. A ``transcendental`` number may or may not be ``real``
        but can never be ``rational``.
        [transcendental]_
      - | ``== (complex & !algebraic)``

    * - ``extended_real``
      - An element of the extended real number line,
        $x\in\overline{\mathbb{R}}$ where
        $\overline{\mathbb{R}}=\mathbb{R}\cup\{-\infty,+\infty\}$. An
        ``extended_real`` number is either ``real`` or $\pm\infty$. The
        relational operators ``<``, ``<=``, ``>=`` and ``>`` are defined only
        for expressions that are ``extended_real``.
        [extended_real]_
      - | ``-> commutative``

    * - ``real``
      - A real number, $x\in\mathbb{R}$. All ``real`` numbers are ``finite``
        and ``complex`` (the set of reals is a subset of the set of complex
        numbers).  Includes all ``rational`` numbers. A ``real`` number is
        either ``negative``, ``zero`` or ``positive``.
        [real]_
      - | ``-> complex``
        | ``== (extended_real & finite)``
        | ``== (negative | zero | positive)``
        | ``-> hermitian``

    * - ``imaginary``
      - An imaginary number, $z\in\mathbb{I}-\{0\}$. A number of the form $z=yi$
        where $y$ is real, $y\ne 0$ and $i=\sqrt{-1}$. All ``imaginary``
        numbers are ``complex`` and not ``real``. Note in particular that
        ``zero`` is `not` considered ``imaginary`` in SymPy.
        [imaginary]_
      - | ``-> complex``
        | ``-> antihermitian``
        | ``-> !extended_real``

    * - ``rational``
      - A rational number, $q\in\mathbb{Q}$. Any number of the form
        $\frac{a}{b}$ where $a$ and $b$ are integers and $b \ne 0$. All
        ``rational`` numbers are ``real`` and ``algebraic``.  Includes all
        ``integer`` numbers.
        [rational]_
      - | ``-> real``
        | ``-> algebraic``

    * - ``irrational``
      - A real number that is not rational, $x\in\mathbb{R}-\mathbb{Q}$.
        [irrational]_
      - | ``== (real & !rational)``

    * - ``integer``
      - An integer, $a\in\mathbb{Z}$. All integers are ``rational``.  Includes
        ``zero`` and all ``prime``, ``composite``, ``even`` and ``odd`` numbers.
        [integer]_
      - | ``-> rational``

    * - ``noninteger``
      - An extended real number that is not an integer,
        $x\in\overline{\mathbb{R}}-\mathbb{Z}$.
      - | ``== (extended_real & !integer)``

    * - ``even``
      - An even number, $e\in\{2k: k\in\mathbb{Z}\}$. All ``even`` numbers are
        ``integer`` numbers. Includes ``zero``.
        [parity]_
      - | ``-> integer``
        | ``-> !odd``

    * - ``odd``
      - An odd number, $o\in\{2k + 1: k\in\mathbb{Z}\}$. All ``odd`` numbers are
        ``integer`` numbers.
        [parity]_
      - | ``-> integer``
        | ``-> !even``

    * - ``prime``
      - A prime number, $p\in\mathbb{P}$. All ``prime`` numbers are
        ``positive`` and ``integer``.
        [prime]_
      - | ``-> integer``
        | ``-> positive``

    * - ``composite``
      - A composite number, $c\in\mathbb{N}-(\mathbb{P}\cup\{1\})$. A positive
        integer that is the product of two or more primes. A ``composite``
        number is always a ``positive`` ``integer`` and is not ``prime``.
        [composite]_
      - | ``-> (integer & positive & !prime)``
        | ``!composite -> (!positive | !even | prime)``

    * - ``zero``
      - The number $0$. An expression with ``zero=True`` represents the
        number ``0`` which is an ``integer``.
        [zero]_
      - | ``-> even & finite``
        | ``== (extended_nonnegative & extended_nonpositive)``
        | ``== (nonnegative & nonpositive)``

    * - ``nonzero``
      - A nonzero real number, $x\in\mathbb{R}-\{0\}$. A ``nonzero`` number
        is always ``real`` and can not be ``zero``.
      - | ``-> real``
        | ``== (extended_nonzero & finite)``

    * - ``extended_nonzero``
      - A member of the extended reals that is not zero,
        $x\in\overline{\mathbb{R}}-\{0\}$.
      - | ``== (extended_real & !zero)``

    * - ``positive``
      - A positive real number, $x\in\mathbb{R}, x>0$. All ``positive``
        numbers are ``finite`` so ``oo`` is not ``positive``.
        [positive]_
      - | ``== (nonnegative & nonzero)``
        | ``== (extended_positive & finite)``

    * - ``nonnegative``
      - A nonnegative real number, $x\in\mathbb{R}, x\ge 0$. All ``nonnegative``
        numbers are ``finite`` so ``oo`` is not ``nonnegative``.
        [positive]_
      - | ``== (real & !negative)``
        | ``== (extended_nonnegative & finite)``

    * - ``negative``
      - A negative real number, $x\in\mathbb{R}, x<0$. All ``negative``
        numbers are ``finite`` so ``-oo`` is not ``negative``.
        [negative]_
      - | ``== (nonpositive & nonzero)``
        | ``== (extended_negative & finite)``

    * - ``nonpositive``
      - A nonpositive real number, $x\in\mathbb{R}, x\le 0$. All ``nonpositive``
        numbers are ``finite`` so ``-oo`` is not ``nonpositive``.
        [negative]_
      - | ``== (real & !positive)``
        | ``== (extended_nonpositive & finite)``

    * - ``extended_positive``
      - A positive extended real number, $x\in\overline{\mathbb{R}}, x>0$.
        An ``extended_positive`` number is either ``positive`` or ``oo``.
        [extended_real]_
      - | ``== (extended_nonnegative & extended_nonzero)``

    * - ``extended_nonnegative``
      - A nonnegative extended real number, $x\in\overline{\mathbb{R}}, x\ge 0$.
        An ``extended_nonnegative`` number is either ``nonnegative`` or ``oo``.
        [extended_real]_
      - | ``== (extended_real & !extended_negative)``

    * - ``extended_negative``
      - A negative extended real number, $x\in\overline{\mathbb{R}}, x<0$.
        An ``extended_negative`` number is either ``negative`` or ``-oo``.
        [extended_real]_
      - | ``== (extended_nonpositive & extended_nonzero)``

    * - ``extended_nonpositive``
      - A nonpositive extended real number, $x\in\overline{\mathbb{R}}, x\le 0$.
        An ``extended_nonpositive`` number is either ``nonpositive`` or ``-oo``.
        [extended_real]_
      - | ``== (extended_real & !extended_positive)``

References for the above definitions
------------------------------------

.. [commutative] https://en.wikipedia.org/wiki/Commutative_property
.. [infinite] https://en.wikipedia.org/wiki/Infinity
.. [antihermitian] https://en.wikipedia.org/wiki/Skew-Hermitian_matrix
.. [complex] https://en.wikipedia.org/wiki/Complex_number
.. [algebraic] https://en.wikipedia.org/wiki/Algebraic_number
.. [transcendental] https://en.wikipedia.org/wiki/Transcendental_number
.. [extended_real] https://en.wikipedia.org/wiki/Extended_real_number_line
.. [real] https://en.wikipedia.org/wiki/Real_number
.. [imaginary] https://en.wikipedia.org/wiki/Imaginary_number
.. [rational] https://en.wikipedia.org/wiki/Rational_number
.. [irrational] https://en.wikipedia.org/wiki/Irrational_number
.. [integer] https://en.wikipedia.org/wiki/Integer
.. [parity] https://en.wikipedia.org/wiki/Parity_%28mathematics%29
.. [prime] https://en.wikipedia.org/wiki/Prime_number
.. [composite] https://en.wikipedia.org/wiki/Composite_number
.. [zero] https://en.wikipedia.org/wiki/0
.. [positive] https://en.wikipedia.org/wiki/Positive_real_numbers
.. [negative] https://en.wikipedia.org/wiki/Negative_number


Implications
============

The assumptions system uses the inference rules to infer new predicates beyond
those immediately specified when creating a symbol:

    >>> x = Symbol('x', real=True, negative=False, zero=False)
    >>> x.is_positive
    True

Although ``x`` was not explicitly declared ``positive`` it can be inferred
from the predicates that were given explicitly. Specifically one of the
inference rules is ``real == negative | zero | positive`` so if ``real`` is
``True`` and both ``negative`` and ``zero`` are ``False`` then ``positive``
must be ``True``.

In practice the assumption inference rules mean that it is not necessary to
include redundant predicates for example a positive real number can be simply
be declared as positive:

    >>> x1 = Symbol('x1', positive=True, real=True)
    >>> x2 = Symbol('x2', positive=True)
    >>> x1.is_real
    True
    >>> x2.is_real
    True
    >>> x1.assumptions0 == x2.assumptions0
    True

Combining predicates that are inconsistent will give an error:

    >>> x = Symbol('x', commutative=False, real=True)
    Traceback (most recent call last):
    ...
    InconsistentAssumptions: {
          algebraic: False,
          commutative: False,
          complex: False,
          composite: False,
          even: False,
          extended_negative: False,
          extended_nonnegative: False,
          extended_nonpositive: False,
          extended_nonzero: False,
          extended_positive: False,
          extended_real: False,
          imaginary: False,
          integer: False,
          irrational: False,
          negative: False,
          noninteger: False,
          nonnegative: False,
          nonpositive: False,
          nonzero: False,
          odd: False,
          positive: False,
          prime: False,
          rational: False,
          real: False,
          transcendental: False,
          zero: False}, real=True


Interpretation of the predicates
================================

Although the predicates are defined in the table above it is worth taking some
time to think about how to interpret them. Firstly many of the concepts
referred to by the predicate names like "zero", "prime", "rational" etc have
a basic meaning in mathematics but can also have more general meanings. For
example when dealing with matrices a matrix of all zeros might be referred to
as "zero". The predicates in the assumptions system do not allow any
generalizations such as this. The predicate ``zero`` is strictly reserved for
the plain number $0$. Instead matrices have an
:py:meth:`~.MatrixBase.is_zero_matrix` property for this purpose (although
that property is not strictly part of the assumptions system):

    >>> from sympy import Matrix
    >>> M = Matrix([[0, 0], [0, 0]])
    >>> M.is_zero
    False
    >>> M.is_zero_matrix
    True

Similarly there are generalisations of the integers such as the Gaussian
integers which have a different notion of prime number. The ``prime``
predicate in the assumptions system does not include those and strictly refers
only to the standard prime numbers $\mathbb{P} = \{2, 3, 5, 7, 11, \cdots\}$.
Likewise ``integer`` only means the standard concept of the integers
$\mathbb{Z} = \{0, \pm 1, \pm 2, \cdots\}$, ``rational`` only means the
standard concept of the rational numbers $\mathbb{Q}$ and so on.

The predicates set up schemes of subsets such as the chain beginning with the
complex numbers which are considered as a superset of the reals which are in
turn a superset of the rationals and so on. The chain of subsets

.. math::

    \mathbb{Z} \subset \mathbb{Q} \subset \mathbb{R} \subset \mathbb{C}

corresponds to the chain of implications in the assumptions system

.. code-block:: C

    integer -> rational -> real -> complex

A "vanilla" symbol with no assumptions explicitly attached is not known to
belong to any of these sets and is not even known to be finite:

    >>> x = Symbol('x')
    >>> x.assumptions0
    {'commutative': True}
    >>> print(x.is_commutative)
    True
    >>> print(x.is_rational)
    None
    >>> print(x.is_complex)
    None
    >>> print(x.is_real)
    None
    >>> print(x.is_integer)
    None
    >>> print(x.is_finite)
    None

It is hard for SymPy to know what it can do with such a symbol that is not
even known to be finite or complex so it is generally better to give some
assumptions to the symbol explicitly. Many parts of SymPy will implicitly
treat such a symbol as complex and in some cases SymPy will permit
manipulations that would not strictly be valid given that ``x`` is not known
to be finite. In a formal sense though very little is known about a vanilla
symbol which makes manipulations involving it difficult.

Defining *something* about a symbol can make a big difference. For example if
we declare the symbol to be an integer then this implies a suite of other
predicates that will help in further manipulations:

    >>> n = Symbol('n', integer=True)
    >>> n.assumptions0
    {'algebraic': True,
     'commutative': True,
     'complex': True,
     'extended_real': True,
     'finite': True,
     'hermitian': True,
     'imaginary': False,
     'infinite': False,
     'integer': True,
     'irrational': False,
     'noninteger': False,
     'rational': True,
     'real': True,
     'transcendental': False}

These assumptions can lead to very significant simplifications e.g.
``integer=True`` gives:

    >>> from sympy import sin, pi
    >>> n1 = Symbol('n1')
    >>> n2 = Symbol('n2', integer=True)
    >>> sin(n1 * pi)
    sin(pi*n1)
    >>> sin(n2 * pi)
    0

Replacing a whole expression with $0$ is about as good as simplification can
get!

It is normally advisable to set as many assumptions as possible on any symbols
so that expressions can be simplified as much as possible. A common
misunderstanding leads to defining a symbol with a ``False`` predicate e.g.:

    >>> x = Symbol('x', negative=False)
    >>> print(x.is_negative)
    False
    >>> print(x.is_nonnegative)
    None
    >>> print(x.is_real)
    None
    >>> print(x.is_complex)
    None
    >>> print(x.is_finite)
    None

If the intention is to say that ``x`` is a real number that is not positive
then that needs to be explicitly stated. In the context that the symbol is
known to be real, the predicate ``positive=False`` becomes much more
meaningful:

    >>> x = Symbol('x', real=True, negative=False)
    >>> print(x.is_negative)
    False
    >>> print(x.is_nonnegative)
    True
    >>> print(x.is_real)
    True
    >>> print(x.is_complex)
    True
    >>> print(x.is_finite)
    True

A symbol declared as ``Symbol('x', real=True, negative=False)`` is equivalent
to a symbol declared as ``Symbol('x', nonnegative=True)``. Simply declaring a
symbol as ``Symbol('x', positive=False)`` does not allow the assumptions
system to conclude much about it because a vanilla symbol is not known to be
finite or even complex.

A related confusion arises with ``Symbol('x', complex=True)`` and
``Symbol('x', real=False)``. Often when either of these is used neither is
what is actually wanted. The first thing to understand is that all real
numbers are complex so a symbol created with ``real=True`` will also have
``complex=True`` and a symbol created with ``complex=True`` will not have
``real=False``. If the intention was to create a complex number that is not
a real number then it should be ``Symbol('x', complex=True, real=False)``. On
the other hand declaring ``real=False`` alone is not sufficient to conclude
that ``complex=True`` because knowing that it is not a real number does not
tell us whether it is finite or whether or not it is some completely different
kind of object from a complex number.

A vanilla symbol is defined by not knowing whether it is ``finite`` etc but
there is no clear definition of what it *should* actually represent. It is
tempting to think of it as an "arbitrary complex number or possibly one of the
infinities" but there is no way to query an arbitrary (non-symbol) expression
in order to determine if it meets those criteria. It is important to bear in
mind that within the SymPy codebase and potentially in downstream libraries
many other kinds of mathematical objects can be found that might also have
``commutative=True`` while being something very different from an ordinary
number (in this context even SymPy's standard infinities are considered
"ordinary").

The only predicate that is applied by default for a symbol is ``commutative``.
We can also declare a symbol to be *noncommutative* e.g.:

    >>> x, y = symbols('x, y', commutative=False)
    >>> z = Symbol('z')  # defaults to commutative=True
    >>> x*y + y*x
    x*y + y*x
    >>> x*z + z*x
    2*z*x

Note here that since ``x`` and ``y`` are both noncommutative ``x`` and ``y``
do not commute so ``x*y != y*x``. On the other hand since ``z`` is commutative
``x`` and ``z`` commute and ``x*z == z*x`` even though ``x`` is
noncommutative.

The interpretation of what a vanilla symbol represents is unclear but the
interpretation of an expression with ``commutative=False`` is entirely
obscure. Such an expression is necessarily not a complex number or an
extended real or any of the standard infinities (even ``zoo`` is commutative).
We are left with very little that we can say about what such an expression
*does* represent.


.. _assumptions-guide-other-is-properties:

Other is_* properties
=====================

There are many properties and attributes in SymPy that that have names
beginning with ``is_`` that look similar to the properties used in the
(old) assumptions system but are not in fact part of the assumptions system.
Some of these have a similar meaning and usage as those of the assumptions
system such as the :py:meth:`~.MatrixBase.is_zero_matrix` property shown
above.  Another example is the ``is_empty`` property of sets:

    >>> from sympy import FiniteSet, Intersection
    >>> S1 = FiniteSet(1, 2)
    >>> S1
    {1, 2}
    >>> print(S1.is_empty)
    False
    >>> S2 = Intersection(FiniteSet(1), FiniteSet(Symbol('x')))
    >>> S2
    Intersection({1}, {x})
    >>> print(S2.is_empty)
    None

The ``is_empty`` property gives a fuzzy-bool indicating whether or not a
:class:`~.Set` is the empty set. In the example of ``S2`` it is not possible to know
whether or not the set is empty without knowing whether or not ``x`` is equal
to ``1`` so ``S2.is_empty`` gives ``None``. The ``is_empty`` property for sets
plays a similar role to the ``is_zero`` property for numbers in the
assumptions system: ``is_empty`` is normally only ``True`` for the
:class:`~.EmptySet` object but it is still useful to be able to distinguish between
the cases where ``is_empty=False`` and ``is_empty=None``.

Although ``is_zero_matrix`` and ``is_empty`` are used for similar purposes to
the assumptions properties such as ``is_zero`` they are not part of the
(old) assumptions system. There are no associated inference rules connecting
e.g.  ``Set.is_empty`` and ``Set.is_finite_set`` because the inference rules
are part of the (old) assumptions system which only deals with the predicates
listed in the table above. It is not possible to declare a
:class:`~.MatrixSymbol` with e.g. ``zero_matrix=False`` and there is no
``SetSymbol`` class but if there was it would not have a system for
understanding predicates like ``empty=False``.

The properties :py:meth:`~.MatrixBase.is_zero_matrix` and ``is_empty`` are
similar to those of the assumptions system because they concern *semantic*
aspects of an expression. There are a large number of other properties that
focus on *structural* aspects such as ``is_Number``,
:py:meth:`~.Expr.is_number`, :py:meth:`~.Basic.is_comparable`. Since these
properties refer to structural aspects of an expression they will always give
``True`` or ``False`` rather than a fuzzy bool that also has the possibility
of being ``None``. Capitalised properties such as ``is_Number`` are usually
shorthand for ``isinstance`` checks e.g.:

    >>> from sympy import Number, Rational
    >>> x = Rational(1, 2)
    >>> isinstance(x, Number)
    True
    >>> x.is_Number
    True
    >>> y = Symbol('y', rational=True)
    >>> isinstance(y, Number)
    False
    >>> y.is_Number
    False

The :class:`~.Number` class is the superclass for :class:`~.Integer`,
:class:`~.Rational` and :class:`~.Float` so any instance of :class:`~.Number`
represents a concrete number with a known value. A symbol such as ``y`` that
is declared with ``rational=True`` might represent the same value as ``x`` but
it is not a concrete number with a known value so this is a structural rather
than a semantic distinction.  Properties like ``is_Number`` are sometimes used
in SymPy in place of e.g. ``isinstance(obj, Number)`` because they do not have
problems with circular imports and checking ``x.is_Number`` can be faster than
a call to ``isinstance``.

The :attr:`~.Expr.is_number` (lower-case) property is very different from
``is_Number``. The :attr:`~.Expr.is_number` property is ``True`` for any
expression that can be numerically evaluated to a floating point complex
number with :py:meth:`~.EvalfMixin.evalf`:

    >>> from sympy import I
    >>> expr1 = I + sqrt(2)
    >>> expr1
    sqrt(2) + I
    >>> expr1.is_number
    True
    >>> expr1.evalf()
    1.4142135623731 + 1.0*I
    >>> x = Symbol('x')
    >>> expr2 = 1 + x
    >>> expr2
    x + 1
    >>> expr2.is_number
    False
    >>> expr2.evalf()
    x + 1.0

The primary reason for checking ``expr.is_number`` is to predict whether a
call to :py:meth:`~.EvalfMixin.evalf` will fully evaluate. The
:py:meth:`~.Basic.is_comparable` property is similar to
:py:meth:`~.Expr.is_number` except that if ``is_comparable`` gives ``True``
then the expression is guaranteed to numerically evaluate to a *real*
:class:`~.Float`.  When ``a.is_comparable`` and ``b.is_comparable`` the
inequality ``a < b`` should be resolvable as something like ``a.evalf() <
b.evalf()``.

The full set of ``is_*`` properties, attributes and methods in SymPy is
large. It is important to be clear though that only those that are listed in
the table of predicates above are actually part of the assumptions system. It
is only those properties that are involved in the *mechanism* that implements
the assumptions system which is explained below.


Implementing assumptions handlers
=================================

We will now work through an example of how to implement a SymPy symbolic
function so that we can see how the old assumptions are used internally. SymPy
already has an ``exp`` function which is defined for all complex numbers but
we will define an ``expreal`` function which is restricted to real arguments.

    >>> from sympy import Function
    >>> from sympy.core.logic import fuzzy_and, fuzzy_or
    >>>
    >>> class expreal(Function):
    ...     """exponential function E**x restricted to the extended reals"""
    ...
    ...     is_extended_nonnegative = True
    ...
    ...     @classmethod
    ...     def eval(cls, x):
    ...         # Validate the argument
    ...         if x.is_extended_real is False:
    ...             raise ValueError("non-real argument to expreal")
    ...         # Evaluate for special values
    ...         if x.is_zero:
    ...             return S.One
    ...         elif x.is_infinite:
    ...             if x.is_extended_negative:
    ...                 return S.Zero
    ...             elif x.is_extended_positive:
    ...                 return S.Infinity
    ...
    ...     @property
    ...     def x(self):
    ...         return self.args[0]
    ...
    ...     def _eval_is_finite(self):
    ...         return fuzzy_or([self.x.is_real, self.x.is_extended_nonpositive])
    ...
    ...     def _eval_is_algebraic(self):
    ...         if fuzzy_and([self.x.is_rational, self.x.is_nonzero]):
    ...             return False
    ...
    ...     def _eval_is_integer(self):
    ...         if self.x.is_zero:
    ...             return True
    ...
    ...     def _eval_is_zero(self):
    ...         return fuzzy_and([self.x.is_infinite, self.x.is_extended_negative])

The ``Function.eval`` method is used to pick up on special values of the function so
that we can return a different object if it would be a simplification. When
``expreal(x)`` is called the ``expreal.__new__`` class method (defined in the
superclass ``Function``) will call ``expreal.eval(x)``. If ``expreal.eval``
returns something other than ``None`` then that will be returned instead of an
unevaluated ``expreal(x)``:

    >>> from sympy import oo
    >>> expreal(1)
    expreal(1)
    >>> expreal(0)
    1
    >>> expreal(-oo)
    0
    >>> expreal(oo)
    oo

Note that the ``expreal.eval`` method does not compare the argument using
``==``. The special values are verified using the assumptions system to query
the properties of the argument. That means that the ``expreal`` method can
also evaluate for different forms of expression that have matching properties
e.g.

    >>> x = Symbol('x', extended_negative=True, infinite=True)
    >>> x
    x
    >>> expreal(x)
    0

Of course the assumptions system can only resolve a limited number of special
values so most ``eval`` methods will also check against some special values
with ``==`` but it is preferable to check e.g. ``x.is_zero`` rather than
``x==0``.

Note also that the ``expreal.eval`` method validates that the argument is
real. We want to allow $\pm\infty$ as arguments to ``expreal`` so we check for
``extended_real`` rather than ``real``. If the argument is not extended real
then we raise an error:

    >>> expreal(I)
    Traceback (most recent call last):
    ...
    ValueError: non-real argument to expreal

Importantly we check ``x.is_extended_real is False`` rather than ``not
x.is_extended_real`` which means that we only reject the argument if it is
*definitely* not extended real: if ``x.is_extended_real`` gives ``None`` then
the argument will not be rejected. The first reason for allowing
``x.is_extended_real=None`` is so that a vanilla symbol can be used with
``expreal``. The second reason is that an assumptions query can always give
``None`` even in cases where an argument is definitely real e.g.:

    >>> x = Symbol('x')
    >>> print(x.is_extended_real)
    None
    >>> expreal(x)
    expreal(x)
    >>> expr = (1 + I)/sqrt(2) + (1 - I)/sqrt(2)
    >>> print(expr.is_extended_real)
    None
    >>> expr.expand()
    sqrt(2)
    >>> expr.expand().is_extended_real
    True
    >>> expreal(expr)
    expreal(sqrt(2)*(1 - I)/2 + sqrt(2)*(1 + I)/2)

Validating the argument in ``expreal.eval`` does mean that it will not be
validated when ``evaluate=False`` is passed but there is not really a better
place to perform the validation:

    >>> expreal(I, evaluate=False)
    expreal(I)

The ``extended_nonnegative`` class attribute and the ``_eval_is_*`` methods on
the ``expreal`` class implement queries in the assumptions system for
instances of ``expreal``:

    >>> expreal(2)
    expreal(2)
    >>> expreal(2).is_finite
    True
    >>> expreal(2).is_integer
    False
    >>> expreal(2).is_rational
    False
    >>> expreal(2).is_algebraic
    False
    >>> z = expreal(-oo, evaluate=False)
    >>> z
    expreal(-oo)
    >>> z.is_integer
    True
    >>> x = Symbol('x', real=True)
    >>> expreal(x)
    expreal(x)
    >>> expreal(x).is_nonnegative
    True

The assumptions system resolves queries like ``expreal(2).is_finite`` using
the corresponding handler ``expreal._eval_is_finite`` and *also* the
implication rules. For example it is known that ``expreal(2).is_rational`` is
``False`` because ``expreal(2)._eval_is_algebraic`` returns ``False`` and
there is an implication rule ``rational -> algebraic``. This means that an
``is_rational`` query can be resolved in this case by the
``_eval_is_algebraic`` handler. It is actually better not to implement
assumptions handlers for every possible predicate but rather to try and
identify a minimal set of handlers that can resolve as many queries as
possible with as few checks as possible.

Another point to note is that the ``_eval_is_*`` methods only make assumptions
queries on the argument ``x`` and do not make any assumptions queries on
``self``. Recursive assumptions queries on the same object will interfere with
the assumptions implications resolver potentially leading to non-deterministic
behaviour so they should not be used (there are examples of this in the SymPy
codebase but they should be removed).

Many of the ``expreal`` methods implicitly return ``None``. This is a common
pattern in the assumptions system. The ``eval`` method and the ``_eval_is_*``
methods can all return ``None`` and often will. A Python function that ends
without reaching a ``return`` statement will implicitly return ``None``. We
take advantage of this by leaving out many of the ``else`` clauses from the
``if`` statements and allowing ``None`` to be returned implicitly. When
following the control flow of these methods it is important to bear in mind
firstly that any queried property can give ``True``, ``False`` or ``None`` and
also that any function will implicitly return ``None`` if all of the
conditionals fail.


Mechanism of the assumptions system
===================================

.. note:: This section describes internal details that could change in a
          future SymPy version.

This section will explain the inner workings of the assumptions system. It is
important to understand that these inner workings are implementation details
and could change from one SymPy version to another. This explanation is
written as of SymPy 1.7. Although the (old) assumptions system has many
limitations (discussed in the next section) it is a mature system that is used
extensively in SymPy and has been well optimised for its current usage. The
assumptions system is used implicitly in most SymPy operations to control
evaluation of elementary expressions.

There are several stages in the implementation of the assumptions system
within a SymPy process that lead up to the evaluation of a single query in the
assumptions system. Briefly these are:

1. At import time the assumptions rules defined in
   ``sympy/core/assumptions.py`` are processed into a canonical form ready for
   efficiently applying the implication rules. This happens once when SymPy is
   imported before even the :class:`~.Basic` class is defined.

2. The ``Basic.__init_subclass__`` method will post-process every
   :class:`~.Basic` subclass to add the relevant properties needed for
   assumptions queries. This also adds the ``default_assumptions`` attribute
   to the class. This happens each time a :class:`~.Basic` subclass is
   defined (when its containing module is imported).

3. Every :class:`~.Basic` instance initially uses the ``default_assumptions`` class
   attribute. When an assumptions query is made on a :class:`~.Basic` instance
   in the first instance the query will be answered from the
   ``default_assumptions`` for the class.

4. If there is no cached value for the assumptions query in the
   ``default_assumptions`` for the class then the default assumptions will be
   copied to make an assumptions cache for the instance. Then the ``_ask()``
   function is called to resolve the query which will firstly call the
   relevant instance handler ``_eval_is`` method. If the handler returns
   non-None then the result will be cached and returned.

5. If the handler does not exist or gives None then the implications resolver
   is tried. This will enumerate (in a randomised order) all possible
   combinations of predicates that could potentially be used to resolve the
   query under the implication rules. In each case the handler ``_eval_is``
   method will be called to see if it gives non-None. If any combination of
   handlers and implication rules leads to a definitive result for the query
   then that result is cached in the instance cache and returned.

6. Finally if the implications resolver failed to resolve the query then the
   query is considered unresolvable. The value of ``None`` for the query is
   cached in the instance cache and returned.

The assumptions rules defined in ``sympy/core/assumptions.py`` are given in
forms like ``real ==  negative | zero | positive``. When this module is
imported these are converted into a ``FactRules`` instance called
``_assume_rules``. This preprocesses the implication rules into the form of
"A" and "B" rules that can be used for the implications resolver. This is
explained in the code in ``sympy/core/facts.py``. We can access this internal
object directly like (full output omitted):

    >>> from sympy.core.assumptions import _assume_rules
    >>> _assume_rules.defined_facts   # doctest: +SKIP
    {'algebraic',
     'antihermitian',
     'commutative',
     'complex',
     'composite',
     'even',
     ...
    >>> _assume_rules.full_implications   # doctest: +SKIP
    defaultdict(set,
                {('extended_positive', False): {('composite', False),
      ('positive', False),
      ('prime', False)},
     ('finite', False): {('algebraic', False),
      ('complex', False),
      ('composite', False),
      ...

The ``Basic.__init_subclass__`` method will inspect the attributes of each
``Basic`` class to see if any assumptions related attributes are defined. An
example of these is the ``is_extended_nonnegative = True`` attribute defined
in the ``expreal`` class. The implications of any such attributes will be
used to precompute any statically knowable assumptions. For example
``is_extended_nonnegative=True`` implies ``real=True`` etc. A ``StdFactKB``
instance is created for the class which stores those assumptions whose values
are known at this stage. The ``StdFactKB`` instance is assigned as the class
attribute ``default_assumptions``. We can see this with

    >>> from sympy import Expr
    ...
    >>> class A(Expr):
    ...     is_positive = True
    ...
    ...     def _eval_is_rational(self):
    ...         # Let's print something to see when this method is called...
    ...         print('!!! calling _eval_is_rational')
    ...         return True
    ...
    >>> A.is_positive
    True
    >>> A.is_real  # inferred from is_positive
    True

Although only ``is_positive`` was defined in the class ``A`` it also has
attributes such as ``is_real`` which are inferred from ``is_positive``. The
set of all such assumptions for class ``A`` can be seen in
``default_assumptions`` which looks like a ``dict`` but is in fact a
``StdFactKB`` instance:

    >>> type(A.default_assumptions)
    <class 'sympy.core.assumptions.StdFactKB'>
    >>> A.default_assumptions
    {'commutative': True,
     'complex': True,
     'extended_negative': False,
     'extended_nonnegative': True,
     'extended_nonpositive': False,
     'extended_nonzero': True,
     'extended_positive': True,
     'extended_real': True,
     'finite': True,
     'hermitian': True,
     'imaginary': False,
     'infinite': False,
     'negative': False,
     'nonnegative': True,
     'nonpositive': False,
     'nonzero': True,
     'positive': True,
     'real': True,
     'zero': False}

When an instance of any :class:`~.Basic` subclass is created ``Basic.__new__`` will
assign its ``_assumptions`` attribute which will initially be a reference to
``cls.default_assumptions`` shared amongst all instances of the same class.
The instance will use this to resolve any assumptions queries until that fails
to give a definitive result at which point a copy of
``cls.default_assumptions`` will be created and assigned to the instance's
``_assumptions`` attribute. The copy will be used as a cache to store any
results computed for the instance by its ``_eval_is`` handlers.

When the ``_assumptions`` attribute fails to give the relevant result it is
time to call the ``_eval_is`` handlers. At this point the ``_ask()`` function
is called. The ``_ask()`` function will initially try to resolve a query such
as ``is_rational`` by calling the corresponding method i.e.
``_eval_is_rational``. If that gives non-None then the result is stored in
``_assumptions`` and any implications of that result are computed and stored
as well. At that point the query is resolved and the value returned.

    >>> a = A()
    >>> a._assumptions is A.default_assumptions
    True
    >>> a.is_rational
    !!! calling _eval_is_rational
    True
    >>> a._assumptions is A.default_assumptions
    False
    >>> a._assumptions   # rational now shows as True
    {'algebraic': True,
     'commutative': True,
     'complex': True,
     'extended_negative': False,
     'extended_nonnegative': True,
     'extended_nonpositive': False,
     'extended_nonzero': True,
     'extended_positive': True,
     'extended_real': True,
     'finite': True,
     'hermitian': True,
     'imaginary': False,
     'infinite': False,
     'irrational': False,
     'negative': False,
     'nonnegative': True,
     'nonpositive': False,
     'nonzero': True,
     'positive': True,
     'rational': True,
     'real': True,
     'transcendental': False,
     'zero': False}

If e.g. ``_eval_is_rational`` does not exist or gives ``None`` then ``_ask()``
will try all possibilities to use the implication rules and any other handler
methods such as ``_eval_is_integer``, ``_eval_is_algebraic`` etc that might
possibly be able to give an answer to the original query. If any method leads
to a definite result being known for the original query then that is returned.
Otherwise once all possibilities for using a handler and the implication rules
to resolve the query are exhausted ``None`` will be cached and returned.

    >>> b = A()
    >>> b.is_algebraic    # called _eval_is_rational indirectly
    !!! calling _eval_is_rational
    True
    >>> c = A()
    >>> print(c.is_prime)   # called _eval_is_rational indirectly
    !!! calling _eval_is_rational
    None
    >>> c._assumptions   # prime now shows as None
    {'algebraic': True,
     'commutative': True,
     'complex': True,
     'extended_negative': False,
     'extended_nonnegative': True,
     'extended_nonpositive': False,
     'extended_nonzero': True,
     'extended_positive': True,
     'extended_real': True,
     'finite': True,
     'hermitian': True,
     'imaginary': False,
     'infinite': False,
     'irrational': False,
     'negative': False,
     'nonnegative': True,
     'nonpositive': False,
     'nonzero': True,
     'positive': True,
     'prime': None,
     'rational': True,
     'real': True,
     'transcendental': False,
     'zero': False}

.. note:: In the ``_ask()`` function the handlers are called in a randomised
          order which can mean that execution at this point is
          non-deterministic. Provided all of the different handler methods
          are consistent (i.e. there are no bugs) then the end result will
          still be deterministic. However a bug where two handlers are
          inconsistent can manifest in non-deterministic behaviour because
          this randomisation might lead to the handlers being called in
          different orders when the same program is run multiple times.

Limitations
===========

Combining predicates with or
----------------------------

In the old assumptions we can easily combine predicates with *and* when
creating a Symbol e.g.:

    >>> x = Symbol('x', integer=True, positive=True)
    >>> x.is_positive
    True
    >>> x.is_integer
    True

We can also easily query whether two conditions are jointly satisfied with

    >>> fuzzy_and([x.is_positive, x.is_integer])
    True
    >>> x.is_positive and x.is_integer
    True

However there is no way in the old assumptions to create a :class:`~.Symbol`
with assumptions predicates combined with *or*. For example if we wanted to
say that "x is positive or x is an integer" then it is not possible to create
a :class:`~.Symbol` with those assumptions.

It is also not possible to ask an assumptions query based on *or* e.g. "is
expr an expression that is positive or an integer". We can use e.g.

    >>> fuzzy_or([x.is_positive, x.is_integer])
    True

However if all that is known about ``x`` is that it is possibly positive or
otherwise a negative integer then both queries ``x.is_positive`` and
``x.is_integer`` will resolve to ``None``.  That means that the query becomes

    >>> fuzzy_or([None, None])

which then also gives ``None``.

Relations between different symbols
-----------------------------------

A fundamental limitation of the old assumptions system is that all explicit
assumptions are properties of an individual symbol. There is no way in this
system to make an assumption about the *relationship* between two symbols. One
of the most common requests is the ability to assume something like ``x < y``
but there is no way to even specify that in the old assumptions.

The new assumptions have the theoretical capability that relational
assumptions can be specified. However the algorithms to make use of that
information are not yet implemented and the exact API for specifying
relational assumptions has not been decided upon.

..
  Titles That can be added:
  1. Dynamic changing assumptions
  2. Selectively controlling evaluation
  3. Extensibility
  4. New assumptions
