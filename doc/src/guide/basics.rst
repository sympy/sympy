
.. include:: ../definitions.def

Basics of expressions in SymPy
==============================

SymPy is all about construction and manipulation of *expressions*. By the
term expression we mean mathematical expressions represented in the Python
language using SymPy's classes and objects. Expressions may consist of
symbols, numbers, functions and function applications (and many other) and
operators binding them together (addiction, subtraction, multiplication,
division, exponentiation).

Suppose we want to construct an expression for `x + 1`::

    >>> x = Symbol('x')

    >>> x + 1
    x + 1

    >>> type(_)
    <class 'sympy.core.add.Add'>

Entering ``x + 1`` gave us an instance of :class:`Add` class. This expression
consists of a symbol (``x``), a number (``1``) and addition operator, which
is represented by the topmost class (:class:`Add`). This was the simplest way
of entering an expression for `x + 1`. We could also enter::

    >>> y = Symbol('y')

    >>> x - y + 17 + y - 16 + sin(pi)
    x + 1

In this case SymPy automatically rewrote the input expression and gave its
canonical form, which is ``x + 1`` once again. This is a very important
behavior: all expressions are subject to automatic evaluation, during which
SymPy tries to find a canonical form for expressions, but it doesn't apply
"heroic" measures to achieve this goal. For example the following expression::

    >>> (x**2 - 1)/(x - 1)
     2
    x  - 1
    ──────
    x - 1

is left unsimplified. This is because automatic canonicalization would
lose important information about this expression (`x \not= 1`). We can
use :func:`cancel` remove common factors from the numerator and the
denominator::

    >>> cancel(_)
    x + 1

SymPy never applies any transformations automatically that could cause
information loss or that would result in results that are valid only
almost everywhere. Consider the following expression::

    >>> log(x*y)
    log(x⋅y)

We know that `\log(x y)` is equivalent to `\log x + \log y` and there
is :func:`expand` that is supposed be able to do this::

    >>> expand(_)
    log(x⋅y)

Unfortunately nothing interesting happened. This is because the formula
stated above is not universally valid, e.g.::

    >>> log((-2)*(-3))
    log(6)
    >>> log(-2) + log(-3)
    log(2) + log(3) + 2⋅ⅈ⋅π

It is possible to ignore such cases and expand forcibly::

    >>> expand(log(x*y), force=True)
    log(x) + log(y)

Many other expression manipulation function also support ``force`` option.
Usually a better way is to assign additional knowledge with an expression::

    >>> var('a,b', positive=True)
    (a, b)

    >>> log(a*b)
    log(a⋅b)

    >>> expand(_)
    log(a) + log(b)

In this case ``force=True`` wasn't necessary, because we gave sufficient
information to :func:`expand` so that it was able to decide that the
expansion rule is valid universally for this expression.

Arithmetic operators
--------------------

Arithmetic operators ``+``, ``-``, ``*``, ``/``, ``**`` are mapped to
combinations of three core SymPy's classes: :class:`Add`, :class:`Mul`
and :class:`Pow`, and work the following way:

* ``x + y`` uses :class:`Add` class and ``__add__`` method::

    >>> x + y
    x + y
    >>> type(_)
    <class 'sympy.core.add.Add'>

    >>> x.__add__(y)
    x + y
    >>> type(_)
    <class 'sympy.core.add.Add'>

    >>> Add(x, y)
    x + y
    >>> type(_)
    <class 'sympy.core.add.Add'>

* ``x - y`` uses :class:`Add` and :class:`Mul` classes, and ``__sub__`` method::

    >>> x - y
    x - y
    >>> type(_)
    <class 'sympy.core.add.Add'>
    >>> __.args
    (-y, x)
    >>> type(_[0])
    <class 'sympy.core.mul.Mul'>

    >>> x.__sub__(y)
    x - y
    >>> type(_)
    <class 'sympy.core.add.Add'>
    >>> __.args
    (-y, x)
    >>> type(_[0])
    <class 'sympy.core.mul.Mul'>

    >>> Add(x, -y))
    x - y
    >>> type(_)
    <class 'sympy.core.add.Add'>
    >>> __.args
    (-y, x)
    >>> type(_[0])
    <class 'sympy.core.mul.Mul'>

* ``x*y`` uses :class:`Mul` class and ``__mul__`` method::

    >>> x*y
    x*y
    >>> type(_)
    <class 'sympy.core.mul.Mul'>

    >>> x.__mul__(y)
    x*y
    >>> type(_)
    <class 'sympy.core.mul.Mul'>

    >>> Mul(x, y)
    x*y
    >>> type(_)
    <class 'sympy.core.mul.Mul'>

* ``x/y`` uses :class:`Pow` and :class:`Mul` classes and ``__div__`` method::

    >>> x/y
    x
    ─
    y
    >>> type(_)
    <class 'sympy.core.mul.Mul'>
    >>> __.args
    ⎛   1⎞
    ⎜x, ─⎟
    ⎝   y⎠
    >>> type(_[1])
    <class 'sympy.core.pow.Pow'>

    >>> x.__div__(y)
    x
    ─
    y
    >>> type(_)
    <class 'sympy.core.mul.Mul'>
    >>> __.args
    ⎛   1⎞
    ⎜x, ─⎟
    ⎝   y⎠
    >>> type(_[1])
    <class 'sympy.core.pow.Pow'>

    >>> Mul(x, 1/y)
    x
    ─
    y
    >>> type(_)
    <class 'sympy.core.mul.Mul'>
    >>> __.args
    ⎛   1⎞
    ⎜x, ─⎟
    ⎝   y⎠
    >>> type(_[1])
    <class 'sympy.core.pow.Pow'>

* ``x**y`` uses :class:`Pow` class and ``__pow__`` method::

    >>> x**y
     y
    x
    >>> type(_)
    <class 'sympy.core.pow.Pow'>

    >>> x.__pow__(y)
     y
    x
    >>> type(_)
    <class 'sympy.core.pow.Pow'>

    >>> Pow(x, y)
     y
    x
    >>> type(_)
    <class 'sympy.core.pow.Pow'>

When the first argument is not an instance SymPy's class, e.g. as in ``1 - x``,
then Python falls back to ``__r*__`` methods, which are also implemented in all
SymPy's classes::

    >>> (1).__sub__(x)
    NotImplemented

    >>> x.__rsub__(1)
    -x + 1
    >>> 1 - x
    -x + 1

Tasks
~~~~~

1. Construct an expression for `1 + x + x^2 + \ldots + x^{10}`. Can you
   construct this expression in a different way? Write a function that
   could generate an expression for `1 + x + x^2 + \ldots + x^n` for any
   integer `n >= 0`. Extend this function to allow `n < 0`.

   (:ref:`solution <solution_arith_op_1>`)

2. Write a function that can compute nested powers, e.g. `x^x`, `x^{x^x}` and
   so on. The function should take two parameters: an expression and a positive
   integer `n` that specifies the depth.

   (:ref:`solution <solution_arith_op_2>`)

Building blocks of expressions
------------------------------

Expressions can consist of instances of subclasses of :class:`Expr` class. This
includes:

* numbers::

    >>> Integer(2)
    2
    >>> Rational(1, 2)
    1/2
    >>> Float("1e-1000")
    1.00000000000000e-1000

* symbols::

    >>> Symbol('x')
    x
    >>> Dummy('y')
    y

* numer symbols::

    >>> pi
    π
    >>> E
    ℯ
    >>> Catalan
    Catalan

* functions::

    >>> Function('f')
    f
    >>> sin
    sin
    >>> cos
    cos

* function applications::

    >>> Function('f')(x)
    f(x)
    >>> sin(x)
    sin(x)
    >>> cos(x)
    cos(x)

* operators::

    >>> Add(x, y, z)
    x + y + z
    >>> Mul(x, y, z)
    x⋅y⋅z
    >>> Pow(x, y)
     y
    x
    >>> Or(x, y, z)
    x ∨ y ∨ z
    >>> And(x, y, z)
    x ∧ y ∧ z

* unevaluated operators::

    >>> Derivative(1/x, x)
    d ⎛1⎞
    ──⎜─⎟
    dx⎝x⎠
    >>> Integral(1/x, x)
    ⌠
    ⎮ 1
    ⎮ ─ dx
    ⎮ x
    ⌡
    >>> Sum(1/k, (k, 1, n))
      n
     ___
     \  `
      \   1
       )  ─
      /   k
     /__,
    k = 1

* other::

    >>> Poly(x**2 + y, x)
    Poly(x**2 + y, x, domain='ZZ[y]')
    >>> RootOf(z**5 + z + 3, 2)
          ⎛ 5           ⎞
    RootOf⎝z  + z + 3, 2⎠

This list isn't at all complete and we included only few classes that SymPy
implements that can be used as expression building blocks. Besides those,
SymPy has also very many classes that represent entities that can't be used
for constructing expressions, but can be useful as containers of expressions
or as utilities for expression building blocks.

Tasks
~~~~~

1. Expressions implement a :func:`doit` method. For most types expressions
   it doesn't do anything useful, but in the case of unevaluated operators,
   it executes an action assigned to to an unevaluated operator (it
   differentiates, integrates, etc.). Take advantage of :func:`doit` and
   write a function that generates integral tables for a few polynomials,
   rational functions and elementary functions.

   (:ref:`solution <solution_blocks_1>`)

Foreign types in SymPy
----------------------

SymPy internally expects that all objects it works with are instances of
subclasses of :class:`Basic` class. So why ``x + 1`` works without raising
an exception? The number ``1`` is not a SymPy's type, but::

    >>> type(1)
    <type 'int'>

it's a built-in type. SymPy implements :func:`sympify` function for the task
of converting foreign types to SymPy's types (yes, Python's built-in types
are also considered as foreign). All SymPy's classes, methods and functions
use :func:`sympify` and this is the reason why you can safely write ``x + 1``
instead of more verbose and less convenient ``x + Integer(1)``. Note that
not all functions return instances of SymPy's types. Usually, if a function
is supposed to return a property of an expression, it will use built-in
Python's types, e.g.::

    >>> Poly(x**2 + y).degree(y)
    1
    >>> type(_)
    <type 'int'>

Now see what :func:`sympify` can do. Let's start with built-ins::

    >>> sympify(1)
    1
    >>> type(_)
    <class 'sympy.core.numbers.One'>

    >>> sympify(117)
    117
    >>> type(_)
    <class 'sympy.core.numbers.Integer'>

    >>> sympify(0.5)
    0.500000000000000
    >>> type(_)
    <class 'sympy.core.numbers.Float'>

    >>> from fractions import Fraction

    >>> sympify(Fraction(1, 2))
    1/2
    >>> type(_)
    <class 'sympy.core.numbers.Rational'>

SymPy implements explicit sympification rules, heuristics based on ``__int__``,
``__float__`` and other attributes, and in the worst case scenario it falls
back to parsing string representation of an object. This usually works fine,
but sometimes :func:`sympify` can be wrong::

    >>> from gmpy import mpz, mpq

    >>> sympify(mpz(117))
    117.000000000000
    >>>> type(_)
    <class 'sympy.core.numbers.Float'>

    >>> sympify(mpq(1, 2))
    0.500000000000000
    >>>> type(_)
    <class 'sympy.core.numbers.Float'>

This happens because :func:`sympify` doesn't know about either ``mpz`` or
``mpq``, and it first looks for ``__float__`` attribute, which is implemented
by both those types. Getting float for exact value isn't very useful so let's
extend :func:`sympify` and add support for ``mpz``. The way to achieve this
is to add a new entry to ``converter`` dictionary. ``converter`` takes types
as keys and sympification functions as values. Before we extend this ``dict``,
we have to resolve a little problem with ``mpz``::

    >>> mpz
    <built-in function mpz>

which isn't a type but a function. We can use a little trick here and take
the type of some ``mpz`` object::

    >>> type(mpz(1))
    <type 'mpz'>

Let's now add an entry to ``converter`` for ``mpz``::

    >>> from sympy.core.sympify import converter

    >>> def mpz_to_Integer(obj):
    ...     return Integer(int(obj))
    ...
    ...

    >>> converter[type(mpz(1))] = mpz_to_Integer

We could use ``lambda`` as well. Now we can sympify ``mpz``::

    >>> sympify(mpz(117))
    117
    >>> type(_)
    <class 'sympy.core.numbers.Integer'>

Similar things should be done for ``mpq``. Let's try one more type::

    >>> import numpy

    >>> ar = numpy.array([1, 2, 3])
    >>> sympify(ar)

    >>> sympify(ar)
    Traceback (most recent call last):
    ...
    SympifyError: SympifyError: "could not parse u'[1 2 3]'"

:func:`sympify` isn't aware of ``numpy.ndarray`` and heuristics didn't work,
so it computed string representation of ``ar`` using :func:`str` and tried
to parse is, which failed because::

    >>> str(ar)
    [1 2 3]

We might be tempted to add support for ``numpy.ndarray`` to :func:`sympify`
by treating NumPy's arrays (at least a subset of) as SymPy's matrices, but
matrices aren't sympifiable::

    >>> Matrix(3, 3, lambda i, j: i + j)
    ⎡0  1  2⎤
    ⎢       ⎥
    ⎢1  2  3⎥
    ⎢       ⎥
    ⎣2  3  4⎦
    >>> sympify(_)
    Traceback (most recent call last):
    ...
    SympifyError: SympifyError: 'Matrix cannot be sympified'

We will explain this odd behavior later.

Tasks
~~~~~

1. Add support for ``mpq`` to :func:`sympify`.

   (:ref:`solution <solution_foreign_1>`)

2. SymPy implements :class:`Tuple` class, which provides functionality of
   Python's built-in ``tuple``, but is a subclass of :class:`Basic`. Take
   advantage of this and make :func:`sympify` work for row NumPy arrays,
   for which it should return instances of :class:`Tuple`. Raise
   :exc:`SympifyError` for other classes of arrays.

   (:ref:`solution <solution_foreign_2>`)

The role of symbols
-------------------

Let's now talk about the most important part of expressions: symbols. Symbols
are placeholders, abstract entities that can be filled in with whatever
content we want (unless there are explicit restrictions given). For example
in expression ``x + 1`` we have one symbol ``x``. Let's start fresh Python's
interpreter and issue::

    >>> from sympy import *
    >>> init_printing()

We want to start work with our very advanced ``x + 1`` expression, so we
may be tempted to simply write::

    >>> x + 1
    Traceback (most recent call last):
    ...
    NameError: name 'x' is not defined

For users that come from other symbolic mathematics systems, this behavior
may seem odd, because in those systems, symbols are constructed implicitly
when necessary. In general purpose programming language like Python, we
have to define all objects we want to use before we actually use them. So,
the first thing we have to always do is to construct symbols and assign
them to Python's variables::

    >>> x = Symbol('x')

    >>> x + 1
    x + 1

Now it worked. Symbols are independent of variables, so nothing prevents
you from issuing::

    >>> t = Symbol('a')

Well, besides taste. It's also perfectly valid to create symbols containing
special characters::

    >>> Symbol('+')
    +

``_`` and ``^`` characters in symbols have special meaning and are used to
denote subscripts and superscripts, respectively::

    >>> Symbol('x_1')
    x₁
    >>> Symbol('x^1')
    x¹

If you need more symbols in your expression, you have to define and assign
them all before using them. Later you can reuse existing symbols for other
purposes. To make life easier, SymPy provides several methods for constructing
symbols. The most low-level method is to use :class:`Symbol` class, as we
have been doing it before. However, if you need more symbols, then your can
use :func:`symbols`::

    >>> symbols('x,y,z')
    (x, y, z)

It takes a textual specification of symbols and returns a ``tuple`` with
constructed symbols. :func:`symbols` supports several syntaxes and can make
your life much simpler, when it comes to constructing symbols. First of all,
commas can be followed by or completely replaced by whitespace::

    >>> symbols('x, y, z')
    (x, y, z)
    >>> symbols('x y z')
    (x, y, z)

If you need indexed symbols, then use range syntax::

    >>> symbols("x:5")
    (x₀, x₁, x₂, x₃, x₄)
    >>> symbols('x5:10')
    (x₅, x₆, x₇, x₈, x₉)

You can also create consecutive symbols with lexicographic syntax::

    >>> symbols('a:d')
    (a, b, c, d)

Note that range syntax simulates :func:`range`'s behavior, so it is exclusive,
lexicographic syntax is inclusive, because it makes more sense in this case.

When we issue::

    >>> symbols('u,v')
    (u, v)

we may be tempted to use ``u`` and ``v``::

    >>> u
    Traceback (most recent call last):
    ...
    NameError: name 'u' is not defined

    >>> v
    Traceback (most recent call last):
    ...
    NameError: name 'v' is not defined

We got :exc:`NameError`, because we constructed those symbols, but we didn't
assign them to any variables. This solves the problem::

    >>> u, v = symbols('u,v')
    >>> u, v
    u, v

but is a little redundant, because we have to repeat the same information
twice. To save time and typing effort, SymPy has another function :func:`var`
for constructing symbols, which has exactly the same syntax and semantics
as :func:`symbols`, but it also injects constructed symbols into the global
namespace, making this function very useful in interactive sessions::

    >>> del u, v
    >>> var('u,v)
    (u, v)

    >>> u + v
    u + v

We don't allow to use :func:`var` in SymPy's library code. There is one
more way of constructing symbols, which is related to indexed symbols.
Sometimes we don't know in advance how many symbols will be required to
solve a certain problem. For this case, SymPy has :func:`numbered_symbols`
generator::

    >>> X = numbered_symbols('x')

    >>> X.next()
    x₀

    >>> [ X.next() for i in xrange(5) ]
    [x₁, x₂, x₃, x₄, x₅]

Tasks
~~~~~

1. Implement a function that would generate an expression for `x_1^1 +
   x_2^2 + \ldots + x_n^n`. This function would take two arguments: base
   name for indexed symbols and integer exponent `n >= 1`. What's the
   best approach among the four presented above?

   (:ref:`solution <solution_symbols_1>`)

Obtaining parts of expressions
------------------------------

We already know how to construct expressions, but how to get parts of complex
expressions? The most basic and low-level way of decomposing expressions is to
use ``args`` property::

    >>> x + y + 1
    x + y + 1
    >>> _.args
    (1, y, x)
    >>> map(type)
    [<class 'sympy.core.numbers.One'>, <class 'sympy.core.symbol.Symbol'>, <class 'sympy.core.symbol.Symbol'>]

``args`` always gives a ``tuple`` of instances of SymPy's classes. One should
notice the weird order of elements, which doesn't match printing order. This
happens for classes that in which order of arguments is insignificant. The
most notable examples of such class are :class:`Add` and :class:`Mul` (for
commutative part). In this particular case we can use :func:`as_ordered_terms`
method to get ``args`` in printing order::

    >>> (x + y + 1).as_ordered_terms()
    [x, y, 1]

When dealing which classes that have fixed order of arguments, printing
order and ``args`` order match::

    >>> Derivative(sin(x), x, x)
       2
      d
    ─────(sin(x))
    dx dx

    >>> _.args
    (sin(x), x, x)

Lets suppose that :class:`Cls` represents any SymPy's class and ``expr``
is an instance of this class (``expr = Cls()``). Then the following holds::

    Cls(*expr.args) == expr

This is very useful invariant, because we can easily decompose, modify and
rebuild expressions of various kinds in SymPy exactly the same way. This
invariant is being used in all functions that manipulation expressions.

Let's now use ``args`` to something a little more interesting than simple
decomposition of expressions. Working with expressions, one may be interested
in the depth of such expressions. By viewing expressions as n-ary trees, by
depth we understand the longest path in a tree.

Trees consist of branches and leafs. In SymPy, leafs of expressions are
instances of subclasses of :class:`Atom` class (numbers, symbols, special
constants)::

    >>> Integer(10)
    10
    >>> isinstance(_, Atom)
    True

    >>> pi
    π
    >>> isinstance(_, Atom)
    True

Atoms can be also recognized by the fact that their ``args`` are empty.
Note, however, that this is an implementation detail, and one should use
either :func:`isinstance` built-in function or ``is_Atom`` property to
recognize atoms properly. Everything else than an :class:`Atom` is a
branch.

Let's implement :func:`depth` function:

.. literalinclude:: python/depth.py

The implementation is straightforward. First we check if the input
expression is an atom. In this case we return ``1`` and terminate
recursion. Otherwise :func:`depth` recurses for every argument of
``expr`` and returns ``1`` plus maximum of depths of all branches.

Let's see :func:`depth` in action::

    >>> depth(x)
    1
    >>> depth(x + 1)
    2
    >>> depth(x + sin(x))
    3
    >>> depth(x + sin(x) + sin(cos(x)))
    4

All those examples work as expected. However, not everything is perfect
with this function. Let's look at the following phenomenon::

    >>> depth(Integer(117))
    1
    >>> depth(117)
    Traceback (most recent call last):
    ...
    AttributeError: 'int' object has no attribute 'args'

``117`` is an instance of Python's built-in type :class:`int`, but this type
is not a subclass of :class:`Atom`, so Python choses the other branch in
:func:`depth` and this must fail. Before the last example we pass only
instances of SymPy's expression to :func:`depth`. If we want :func:`depth` to
work also for non-SymPy types, we have to sympify ``expr`` with :func:`sympify`
before using it.

Tasks
~~~~~

1. Change :func:`depth` so that it sympifies its input argument. Rewrite
   :func:`depth` so that is calls :func:`sympify` only once.

   (:ref:`solution <solution_parts_1>`)

2. Add support for iterable containers to :func:`depth`. Containers should
   be treated as branches and have depth defined the same way.

   (:ref:`solution <solution_parts_2>`)

Immutability of expressions
---------------------------

Expressions in SymPy are immutable and cannot be modified by an in-place
operation. This means that a function will always return an object, and
the original expression will not be modified. Consider the following
code::

    >>> var('x,y,a,b')
    (x, y, a, b)

    >>> original = 3*x + 4*y
    >>> modified = original.subs({x: a, y: b})

    >>> original
    3*x + 4*y
    >>> modified
    3*a + 4*b

The output shows that the :func:`subs` method gave a new expression with
symbol ``x`` replaced with symbol ``a`` and symbol ``y`` replaced with
symbol ``b``. The original expression wasn't modified. This behavior
applies to all classes that are subclasses of :class:`Basic`. An exception
to immutability rule is :class:`Matrix`, which allows in-place modifications,
but it is not a subclass of :class:`Basic`::

    >>> Matrix.mro()
    [<class 'sympy.matrices.matrices.Matrix'>, <type 'object'>]

Be also aware of the fact that SymPy's symbols aren't Python's variables (they
just can be assigned to Python's variables), so if you issue::

    >>> u = Symbol('u')
    >>> v = u
    >>> v += 1
    >>> v
    u + 1

then in-place operator ``+=`` constructed an new instance of :class:`Add` and
left the original expression stored in variable ``u`` unchanged::

    >>> u
    u

For efficiency reason, any in-place operator used on elements of a matrix,
modifies the matrix in-place and doesn't waste memory for unnecessary copies.

Tasks
~~~~~

1. This is the first time we used :func:`subs`. This is a very important method
   and we will talk more about it later. However, we can also use :func:`subs`
   to generate some cool looking expressions. Start with ``x**x`` expression
   and substitute in it ``x**x`` for ``x``. What do you get? (make sure you
   use pretty printer) Can you achieve the same effect without :func:`subs`?

   (:ref:`solution <solution_immutability_1>`)

Comparing expressions with ``==``
---------------------------------

Consider the following two expressions::

    >>> f = (x + 1)**2
    >>> f
           2
    (x + y)

    >>> g = x**2 + 2*x + 1
    >>> g
     2
    x  + 2⋅x + 1

We should remember from calculus 101 that those two expressions are
equivalent, because we can use binomial theorem to expand ``f`` and
we will get ``g``. However in SymPy::

    >>> f == g
    False

This is correct result, because SymPy implements structural understanding
of ``==`` operator, not semantic. So, for SymPy ``f`` and ``g`` are very
different expressions.

What to do if we have two variables and we want to know if their contents
are equivalent, but not necessarily structurally equal? There is no simple
answer to this question in general. In the particular case of ``f`` and
``g``, it is sufficient to issue::

    >>> expand(f) == expand(g)
    True

or, based on `f = g \equiv f - g = 0` equivalence::

    >>> expand(f - g) == 0
    True

In case of more complicated expression, e.g. those involving elementary or
special functions, this approach may be insufficient. For example::

    >>> u = sin(x)**2 - 1
    >>> v = cos(x)**2

    >>> u == v
    False
    >>> expand(u - v) == 0
    False

In this case we have to use more advanced term rewriting function::

    >>> simplify(u - v) == 0
    True

The meaning of expressions
--------------------------

Expressions don't have any meaning assigned to them by default. Thus `x + 1`
is simply an expression, not a function or a univariate polynomial. Meaning
is assigned when we use expressions in a context, e.g.::

    >>> div(x**2 - y, x - y)
    ⎛        2    ⎞
    ⎝x + y, y  - y⎠

In this case, ``x**2 - y`` and ``x - y`` where treated as multivariate
polynomials in variables ``x`` and ``y`` (in this order). We could change
this understanding and ask explicitly for polynomials in variables ``y``
and ``x``. This makes :func:`div` return a different result::

    >>> div(x**2 - y, x - y, y, x)
    ⎛    2    ⎞
    ⎝1, x  - x⎠

Quite often SymPy is capable of deriving the most useful understanding of
expressions in a given context. However, there are situations when expressions
simply don't carry enough information to make SymPy perform computations without
telling it explicitly what to do::

    >>> roots(x**2 - y)
    Traceback (most recent call last):
    ...
    PolynomialError: multivariate polynomials are not supported

Here we have to tell :func:`roots` in which variable roots should be computed::

    >>> roots(x**2 - y, x)
    ⎧   ⎽⎽⎽       ⎽⎽⎽   ⎫
    ⎨-╲╱ y : 1, ╲╱ y : 1⎬
    ⎩                   ⎭

Of course the choice of ``y`` is also a valid one, assuming that this is what
you really want. This of course doesn't apply only to polynomials.

Turning strings into expressions
--------------------------------

Suppose we saved the following expression::

    >>> var('x,y')

    >>> expr = x**2 + sin(y) + S(1)/2
    >>> expr
     2            1
    x  + sin(y) + ─
                  2

by printing it with :func:`sstr` printer and storing to a file::

    >>> sstr(expr)
    x**2 + sin(y) + 1/2

    >>> with open("expression.txt", "w") as f:
    ...     f.write(_)
    ...
    ...

We used this kind of printer because we wanted the file to be fairly readable.
Now we want to restore the original expression. First we have to read the text
form from the file::

    >>> with open("expression.txt") as f:
    ...     text_form = f.read()
    ...
    ...

    >>> text_form
    x**2 + sin(y) + 1/2
    >>> type(_)
    <type 'str'>

We could try to try to use :func:`eval` on ``text_form`` but this doesn't give
expected results::

    >>> eval(text_form)
     2
    x  + sin(y) + 0.5

This happens because ``1/2`` isn't understood by Python as rational number
and is equivalent to a problem we had when entering expressions of this kind
in interactive sessions.

To overcome this problem we have to use :func:`sympify`, which implements
:mod:`tokenize`--based parser that allows us to handle this issue::

    >>> sympify(text_form)
     2            1
    x  + sin(y) + ─
                  2
    >>> _ == expr
    True

Let's now consider a more interesting problem. Suppose we define our own function::

    >>> class my_func(Function):
    ...     """Returns zero for integer values. """
    ...
    ...     @classmethod
    ...     def eval(cls, arg):
    ...         if arg.is_Number:
    ...             return 2*arg
    ...
    ...

This function gives twice the input argument if the argument is a number and
doesn't do anything for all other classes of arguments::

    >>> my_func(117)
    234
    >>> my_func(S(1)/2)
    1

    >>> my_func(x)
    my_func(x)
    >>> _.subs(x, 2.1)
    4.20000000000000

    >>> my_func(1) + 1
    3

Let's create an expression that contains :func:`my_func`::

    >>> expr = my_func(x) + 1
    >>> expr
    my_func(x) + 1

    >>> _.subs(x, 1)
    3

Now we will print it using :func:`sstr` printer and sympify the result::

    >>> sympified = sympify(sstr(expr))
    >>> sympified
    my_func(x) + 1

We can use :func:`subs` method to quickly verify the expression is correct::

    >>> sympified.subs(x, 1)
    my_func(1) + 1

This is not exactly what we expected. This happens because::

    >>> expr == sympified
    False

    >>> expr.args
    (1, my_func(x))
    >>> type(_[1]) is my_func
    True

    >>> sympified.args
    (1, my_func(x))
    >>> type(_[1]) is my_func
    False

:func:`sympify` evaluates the given string in the context of ``from sympy import *``
and is not aware of user defined names. We can explicitly pass a mapping between
names and values to it::

    >>> sympify(sstr(expr), {'my_func': my_func})
    my_func(x) + 1
    >>> _.subs(x, 1)
    3

This time we got the desired result. This shows that we have to be careful when
working with expressions encoded as strings. This happens to be even more tricky
when we put assumptions on symbols. Do you remember the example in which we
tried to expand `\log(a b)`? Lets do it once again::

    >>> var('a,b', positive=True)
    (a, b)
    >>> log(a*b).expand()
    log(a) + log(b)

This worked as previously. However, let's now print `\log(a b)`, sympify the
resulting string and expand the restored expression::

    >>> sympify(sstr(log(a*b))).expand()
    log(a⋅b)

This didn't work, because :func:`sympify` doesn't know what ``a`` and ``b``
are, so it assumed that those are symbols and it created them implicitly.
This issue is similar to what we already experienced with :func:`my_func`.

The most reliable approach to storing expression is to use :mod:`pickle`
module. In the case of `\log(a b)` it works like this::

    >>> import pickle
    >>> pickled = pickle.dumps(log(a*b))
    >>> expr = pickle.loads(pickled)
    >>> expr.expand()
    log(a) + log(b)

Unfortunately, due to :mod:`pickle`'s limitations, this doesn't work for
user defined functions like :func:`my_func`::

    >>> pickle.dumps(my_func(x))
    Traceback (most recent call last):
    ...
    PicklingError: Can't pickle my_func: it's not found as __main__.my_func

Tasks
~~~~~

1. Construct a polynomial of degree, let's say, 1000. Use both techniques
   to save and restore this expression. Compare speed of those approaches.
   Verify that the result is correct.

   (:ref:`solution <solution_sympify_1>`)
