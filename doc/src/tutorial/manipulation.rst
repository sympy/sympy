.. _tutorial-manipulation:

==================================
 Advanced Expression Manipulation
==================================

In this section, we discuss some ways that we can perform advanced
manipulation of expressions.

Understanding Expression Trees
==============================

.. sidebar :: Quick Tip

   To play with the ``srepr`` form of expressions in the SymPy Live shell,
   change the output format to ``Repr`` in the settings.

Before we can do this, we need to understand how expressions are represented
in SymPy.  A mathematical expression is represented as a tree.  Let us take
the expression `x^2 + xy`, i.e., ``x**2 + x*y``.  We can see what this
expression looks like internally by using ``srepr``

    >>> from sympy import *
    >>> x, y, z = symbols('x y z')

    >>> expr = x**2 + x*y
    >>> srepr(expr)
    "Add(Pow(Symbol('x'), Integer(2)), Mul(Symbol('x'), Symbol('y')))"

The easiest way to tear this apart is to look at a diagram of the expression
tree:

.. This comes from dotprint(x**2 + x*y, labelfunc=srepr)

.. graphviz::

    digraph{

    # Graph style
    "rankdir"="TD"

    #########
    # Nodes #
    #########

    "Symbol(x)_(0, 0)" ["color"="black", "label"="Symbol('x')", "shape"="ellipse"];
    "Integer(2)_(1, 1)" ["color"="black", "label"="Integer(2)", "shape"="ellipse"];
    "Symbol(y)_(0, 1)" ["color"="black", "label"="Symbol('y')", "shape"="ellipse"];
    "Symbol(x)_(1, 0)" ["color"="black", "label"="Symbol('x')", "shape"="ellipse"];
    "Mul(Symbol(x), Symbol(y))_(0,)" ["color"="black", "label"="Mul", "shape"="ellipse"];
    "Pow(Symbol(x), Integer(2))_(1,)" ["color"="black", "label"="Pow", "shape"="ellipse"];
    "Add(Mul(Symbol(x), Symbol(y)), Pow(Symbol(x), Integer(2)))_()" ["color"="black", "label"="Add", "shape"="ellipse"];

    #########
    # Edges #
    #########

    "Mul(Symbol(x), Symbol(y))_(0,)" -> "Symbol(x)_(0, 0)";
    "Mul(Symbol(x), Symbol(y))_(0,)" -> "Symbol(y)_(0, 1)";
    "Pow(Symbol(x), Integer(2))_(1,)" -> "Symbol(x)_(1, 0)";
    "Pow(Symbol(x), Integer(2))_(1,)" -> "Integer(2)_(1, 1)";
    "Add(Mul(Symbol(x), Symbol(y)), Pow(Symbol(x), Integer(2)))_()" -> "Mul(Symbol(x), Symbol(y))_(0,)";
    "Add(Mul(Symbol(x), Symbol(y)), Pow(Symbol(x), Integer(2)))_()" -> "Pow(Symbol(x), Integer(2))_(1,)";
    }

.. note::

   The above diagram was made using `Graphviz <http://www.graphviz.org/>`_ and
   the :py:meth:`dotprint <sympy.printing.dot.dotprint>` function.

First, let's look at the leaves of this tree.  Symbols are instances of the
class Symbol.  While we have been doing

    >>> x = symbols('x')

we could have also done

    >>> x = Symbol('x')

Either way, we get a Symbol with the name "x" [#symbols-fn]_.  For the number in the
expression, 2, we got ``Integer(2)``.  ``Integer`` is the SymPy class for
integers.  It is similar to the Python built-in type ``int``, except that
``Integer`` plays nicely with other SymPy types.

When we write ``x**2``, this creates a ``Pow`` object.  ``Pow`` is short for
"power".

    >>> srepr(x**2)
    "Pow(Symbol('x'), Integer(2))"

We could have created the same object by calling ``Pow(x, 2)``

    >>> Pow(x, 2)
    x**2

Note that in the ``srepr`` output, we see ``Integer(2)``, the SymPy version of
integers, even though technically, we input ``2``, a Python int.  In general,
whenever you combine a SymPy object with a non-SymPy object via some function
or operation, the non-SymPy object will be converted into a SymPy object.  The
function that does this is ``sympify`` [#sympify-fn]_.

    >>> type(2)
    <... 'int'>
    >>> type(sympify(2))
    <class 'sympy.core.numbers.Integer'>

We have seen that ``x**2`` is represented as ``Pow(x, 2)``.  What about
``x*y``?  As we might expect, this is the multiplication of ``x`` and ``y``.
The SymPy class for multiplication is ``Mul``.

    >>> srepr(x*y)
    "Mul(Symbol('x'), Symbol('y'))"

Thus, we could have created the same object by writing ``Mul(x, y)``.

    >>> Mul(x, y)
    x*y

Now we get to our final expression, ``x**2 + x*y``.  This is the addition of
our last two objects, ``Pow(x, 2)``, and ``Mul(x, y)``.  The SymPy class for
addition is ``Add``, so, as you might expect, to create this object, we use
``Add(Pow(x, 2), Mul(x, y))``.

    >>> Add(Pow(x, 2), Mul(x, y))
    x**2 + x*y

SymPy expression trees can have many branches, and can be quite deep or quite
broad.  Here is a more complicated example

    >>> expr = sin(x*y)/2 - x**2 + 1/y
    >>> srepr(expr)
    "Add(Mul(Integer(-1), Pow(Symbol('x'), Integer(2))), Mul(Rational(1, 2),
    sin(Mul(Symbol('x'), Symbol('y')))), Pow(Symbol('y'), Integer(-1)))"

Here is a diagram

.. dotprint(sin(x*y)/2 - x**2 + 1/y, labelfunc=srepr)

.. graphviz::

    digraph{

    # Graph style
    "rankdir"="TD"

    #########
    # Nodes #
    #########

    "Half()_(0, 0)" ["color"="black", "label"="Rational(1, 2)", "shape"="ellipse"];
    "Symbol(y)_(2, 0)" ["color"="black", "label"="Symbol('y')", "shape"="ellipse"];
    "Symbol(x)_(1, 1, 0)" ["color"="black", "label"="Symbol('x')", "shape"="ellipse"];
    "Integer(2)_(1, 1, 1)" ["color"="black", "label"="Integer(2)", "shape"="ellipse"];
    "NegativeOne()_(2, 1)" ["color"="black", "label"="Integer(-1)", "shape"="ellipse"];
    "NegativeOne()_(1, 0)" ["color"="black", "label"="Integer(-1)", "shape"="ellipse"];
    "Symbol(y)_(0, 1, 0, 1)" ["color"="black", "label"="Symbol('y')", "shape"="ellipse"];
    "Symbol(x)_(0, 1, 0, 0)" ["color"="black", "label"="Symbol('x')", "shape"="ellipse"];
    "Pow(Symbol(x), Integer(2))_(1, 1)" ["color"="black", "label"="Pow", "shape"="ellipse"];
    "Pow(Symbol(y), NegativeOne())_(2,)" ["color"="black", "label"="Pow", "shape"="ellipse"];
    "Mul(Symbol(x), Symbol(y))_(0, 1, 0)" ["color"="black", "label"="Mul", "shape"="ellipse"];
    "sin(Mul(Symbol(x), Symbol(y)))_(0, 1)" ["color"="black", "label"="sin", "shape"="ellipse"];
    "Mul(Half(), sin(Mul(Symbol(x), Symbol(y))))_(0,)" ["color"="black", "label"="Mul", "shape"="ellipse"];
    "Mul(NegativeOne(), Pow(Symbol(x), Integer(2)))_(1,)" ["color"="black", "label"="Mul", "shape"="ellipse"];
    "Add(Mul(Half(), sin(Mul(Symbol(x), Symbol(y)))), Mul(NegativeOne(), Pow(Symbol(x), Integer(2))), Pow(Symbol(y), NegativeOne()))_()" ["color"="black", "label"="Add", "shape"="ellipse"];

    #########
    # Edges #
    #########

    "Pow(Symbol(y), NegativeOne())_(2,)" -> "Symbol(y)_(2, 0)";
    "Pow(Symbol(x), Integer(2))_(1, 1)" -> "Symbol(x)_(1, 1, 0)";
    "Pow(Symbol(x), Integer(2))_(1, 1)" -> "Integer(2)_(1, 1, 1)";
    "Pow(Symbol(y), NegativeOne())_(2,)" -> "NegativeOne()_(2, 1)";
    "Mul(Symbol(x), Symbol(y))_(0, 1, 0)" -> "Symbol(x)_(0, 1, 0, 0)";
    "Mul(Symbol(x), Symbol(y))_(0, 1, 0)" -> "Symbol(y)_(0, 1, 0, 1)";
    "Mul(Half(), sin(Mul(Symbol(x), Symbol(y))))_(0,)" -> "Half()_(0, 0)";
    "Mul(NegativeOne(), Pow(Symbol(x), Integer(2)))_(1,)" -> "NegativeOne()_(1, 0)";
    "sin(Mul(Symbol(x), Symbol(y)))_(0, 1)" -> "Mul(Symbol(x), Symbol(y))_(0, 1, 0)";
    "Mul(NegativeOne(), Pow(Symbol(x), Integer(2)))_(1,)" -> "Pow(Symbol(x), Integer(2))_(1, 1)";
    "Mul(Half(), sin(Mul(Symbol(x), Symbol(y))))_(0,)" -> "sin(Mul(Symbol(x), Symbol(y)))_(0, 1)";
    "Add(Mul(Half(), sin(Mul(Symbol(x), Symbol(y)))), Mul(NegativeOne(), Pow(Symbol(x), Integer(2))), Pow(Symbol(y), NegativeOne()))_()" -> "Pow(Symbol(y), NegativeOne())_(2,)";
    "Add(Mul(Half(), sin(Mul(Symbol(x), Symbol(y)))), Mul(NegativeOne(), Pow(Symbol(x), Integer(2))), Pow(Symbol(y), NegativeOne()))_()" -> "Mul(Half(), sin(Mul(Symbol(x), Symbol(y))))_(0,)";
    "Add(Mul(Half(), sin(Mul(Symbol(x), Symbol(y)))), Mul(NegativeOne(), Pow(Symbol(x), Integer(2))), Pow(Symbol(y), NegativeOne()))_()" -> "Mul(NegativeOne(), Pow(Symbol(x), Integer(2)))_(1,)";
    }

This expression reveals some interesting things about SymPy expression
trees. Let's go through them one by one.

Let's first look at the term ``x**2``.  As we expected, we see ``Pow(x, 2)``.
One level up, we see we have ``Mul(-1, Pow(x, 2))``.  There is no subtraction
class in SymPy.  ``x - y`` is represented as ``x + -y``, or, more completely,
``x + -1*y``, i.e., ``Add(x, Mul(-1, y))``.

    >>> expr = x - y
    >>> srepr(x - y)
    "Add(Symbol('x'), Mul(Integer(-1), Symbol('y')))"

.. dotprint(x - y, labelfunc=srepr)

.. graphviz::

    digraph{

    # Graph style
    "rankdir"="TD"

    #########
    # Nodes #
    #########

    "Symbol(x)_(1,)" ["color"="black", "label"="Symbol('x')", "shape"="ellipse"];
    "Symbol(y)_(0, 1)" ["color"="black", "label"="Symbol('y')", "shape"="ellipse"];
    "NegativeOne()_(0, 0)" ["color"="black", "label"="Integer(-1)", "shape"="ellipse"];
    "Mul(NegativeOne(), Symbol(y))_(0,)" ["color"="black", "label"="Mul", "shape"="ellipse"];
    "Add(Mul(NegativeOne(), Symbol(y)), Symbol(x))_()" ["color"="black", "label"="Add", "shape"="ellipse"];

    #########
    # Edges #
    #########

    "Mul(NegativeOne(), Symbol(y))_(0,)" -> "Symbol(y)_(0, 1)";
    "Mul(NegativeOne(), Symbol(y))_(0,)" -> "NegativeOne()_(0, 0)";
    "Add(Mul(NegativeOne(), Symbol(y)), Symbol(x))_()" -> "Symbol(x)_(1,)";
    "Add(Mul(NegativeOne(), Symbol(y)), Symbol(x))_()" -> "Mul(NegativeOne(), Symbol(y))_(0,)";
    }

Next, look at ``1/y``.  We might expect to see something like ``Div(1, y)``,
but similar to subtraction, there is no class in SymPy for division.  Rather,
division is represented by a power of -1.  Hence, we have ``Pow(y, -1)``.
What if we had divided something other than 1 by ``y``, like ``x/y``?  Let's
see.

    >>> expr = x/y
    >>> srepr(expr)
    "Mul(Symbol('x'), Pow(Symbol('y'), Integer(-1)))"

.. dotprint(x/y, labelfunc=srepr)

.. graphviz::

    digraph{

    # Graph style
    "rankdir"="TD"

    #########
    # Nodes #
    #########

    "Symbol(x)_(0,)" ["color"="black", "label"="Symbol('x')", "shape"="ellipse"];
    "Symbol(y)_(1, 0)" ["color"="black", "label"="Symbol('y')", "shape"="ellipse"];
    "NegativeOne()_(1, 1)" ["color"="black", "label"="Integer(-1)", "shape"="ellipse"];
    "Pow(Symbol(y), NegativeOne())_(1,)" ["color"="black", "label"="Pow", "shape"="ellipse"];
    "Mul(Symbol(x), Pow(Symbol(y), NegativeOne()))_()" ["color"="black", "label"="Mul", "shape"="ellipse"];

    #########
    # Edges #
    #########

    "Pow(Symbol(y), NegativeOne())_(1,)" -> "Symbol(y)_(1, 0)";
    "Pow(Symbol(y), NegativeOne())_(1,)" -> "NegativeOne()_(1, 1)";
    "Mul(Symbol(x), Pow(Symbol(y), NegativeOne()))_()" -> "Symbol(x)_(0,)";
    "Mul(Symbol(x), Pow(Symbol(y), NegativeOne()))_()" -> "Pow(Symbol(y), NegativeOne())_(1,)";
    }

We see that ``x/y`` is represented as ``x*y**-1``, i.e., ``Mul(x, Pow(y,
-1))``.

Finally, let's look at the ``sin(x*y)/2`` term.  Following the pattern of the
previous example, we might expect to see ``Mul(sin(x*y), Pow(Integer(2),
-1))``.  But instead, we have ``Mul(Rational(1, 2), sin(x*y))``.  Rational
numbers are always combined into a single term in a multiplication, so that
when we divide by 2, it is represented as multiplying by 1/2.

Finally, one last note.  You may have noticed that the order we entered our
expression and the order that it came out from ``srepr`` or in the graph were
different.  You may have also noticed this phenonemon earlier in the
tutorial.  For example

     >>> 1 + x
     x + 1

This because in SymPy, the arguments of the commutative operations ``Add`` and
``Mul`` are stored in an arbitrary (but consistent!) order, which is
independent of the order inputted (if you're worried about noncommutative
multiplication, don't be.  In SymPy, you can create noncommutative Symbols
using ``Symbol('A', commutative=False)``, and the order of multiplication for
noncommutative Symbols is kept the same as the input).  Furthermore, as we
shall see in the next section, the printing order and the order in which
things are stored internally need not be the same either.

.. sidebar:: Quick Tip

   The way an expression is represented internally and the way it is printed
   are often not the same.

In general, an important thing to keep in mind when working with SymPy expression
trees is this:  the internal representation of an expression and the way it is
printed need not be the same.  The same is true for the input form.   If some
expression manipulation algorithm is not working in the way you expected it
to, chances are, the internal representation of the object is different from
what you thought it was.

Recursing through an Expression Tree
====================================

Now that you know how expression trees work in SymPy, let's look at how to dig
our way through an expression tree.  Every object in SymPy has two very
important attributes, ``func``, and ``args``.


func
----

``func`` is the head of the object. For example, ``(x*y).func`` is ``Mul``.
Usually it is the same as the class of the object (though there are exceptions
to this rule).

Two notes about ``func``.  First, the class of an object need not be the same
as the one used to create it.  For example

    >>> expr = Add(x, x)
    >>> expr.func
    <class 'sympy.core.mul.Mul'>

We created ``Add(x, x)``, so we might expect ``expr.func`` to be ``Add``, but
instead we got ``Mul``.  Why is that?  Let's take a closer look at ``expr``.

    >>> expr
    2*x

``Add(x, x)``, i.e., ``x + x``, was automatically converted into ``Mul(2,
x)``, i.e., ``2*x``, which is a ``Mul``.   SymPy classes make heavy use of the
``__new__`` class constructor, which, unlike ``__init__``, allows a different
class to be returned from the constructor.

Second, some classes are special-cased, usually for efficiency reasons
[#singleton-fn]_.

    >>> Integer(2).func
    <class 'sympy.core.numbers.Integer'>
    >>> Integer(0).func
    <class 'sympy.core.numbers.Zero'>
    >>> Integer(-1).func
    <class 'sympy.core.numbers.NegativeOne'>

For the most part, these issues will not bother us.  The special classes
``Zero``, ``One``, ``NegativeOne``, and so on are subclasses of ``Integer``,
so as long as you use ``isinstance``, it will not be an issue.

args
----

``args`` are the top-level arguments of the object.  ``(x*y).args`` would be
``(x, y)``.  Let's look at some examples

    >>> expr = 3*y**2*x
    >>> expr.func
    <class 'sympy.core.mul.Mul'>
    >>> expr.args
    (3, x, y**2)

From this, we can see that ``expr == Mul(3, y**2, x)``.  In fact, we can see
that we can completely reconstruct ``expr`` from its ``func`` and its
``args``.

    >>> expr.func(*expr.args)
    3*x*y**2
    >>> expr == expr.func(*expr.args)
    True

Note that although we entered ``3*y**2*x``, the ``args`` are ``(3, x, y**2)``.
In a ``Mul``, the Rational coefficient will come first in the ``args``, but
other than that, the order of everything else follows no special pattern.  To
be sure, though, there is an order.

    >>> expr = y**2*3*x
    >>> expr.args
    (3, x, y**2)

Mul's ``args`` are sorted, so that the same ``Mul`` will have the same
``args``.  But the sorting is based on some criteria designed to make the
sorting unique and efficient that has no mathematical significance.

The ``srepr`` form of our ``expr`` is ``Mul(3, x, Pow(y, 2))``.  What if we
want to get at the ``args`` of ``Pow(y, 2)``.  Notice that the ``y**2`` is in
the third slot of ``expr.args``, i.e., ``expr.args[2]``.

    >>> expr.args[2]
    y**2

So to get the ``args`` of this, we call ``expr.args[2].args``.

    >>> expr.args[2].args
    (y, 2)

Now what if we try to go deeper.  What are the args of ``y``.  Or ``2``.
Let's see.

    >>> y.args
    ()
    >>> Integer(2).args
    ()

They both have empty ``args``.  In SymPy, empty ``args`` signal that we have
hit a leaf of the expression tree.

So there are two possibilities for a SymPy expression. Either it has empty
``args``, in which case it is a leaf node in any expression tree, or it has
``args``, in which case, it is a branch node of any expression tree.  When it
has ``args``, it can be completely rebuilt from its ``func`` and its ``args``.
This is expressed in the key invariant.

.. topic:: Key Invariant

   Every well-formed SymPy expression must either have empty ``args`` or
   satisfy ``expr == expr.func(*expr.args)``.

(Recall that in Python if ``a`` is a tuple, then ``f(*a)`` means to call ``f``
with arguments from the elements of ``a``, e.g., ``f(*(1, 2, 3))`` is the same
as ``f(1, 2, 3)``.)

This key invariant allows us to write simple algorithms that walk expression
trees, change them, and rebuild them into new expressions.

Walking the Tree
----------------

With this knowledge, let's look at how we can recurse through an expression
tree.  The nested nature of ``args`` is a perfect fit for recursive functions.
The base case will be empty ``args``.  Let's write a simple function that goes
through an expression and prints all the ``args`` at each level.

    >>> def pre(expr):
    ...     print(expr)
    ...     for arg in expr.args:
    ...         pre(arg)

See how nice it is that ``()`` signals leaves in the expression tree.  We
don't even have to write a base case for our recursion; it is handled
automatically by the for loop.

Let's test our function.

    >>> expr = x*y + 1
    >>> pre(expr)
    x*y + 1
    1
    x*y
    x
    y

Can you guess why we called our function ``pre``?  We just wrote a pre-order
traversal function for our expression tree.   See if you can write a
post-order traversal function.

Such traversals are so common in SymPy that the generator functions
``preorder_traversal`` and ``postorder_traversal`` are provided to make such
traversals easy.  We could have also written our algorithm as

    >>> for arg in preorder_traversal(expr):
    ...     print(arg)
    x*y + 1
    1
    x*y
    x
    y

.. rubric:: Footnotes

.. [#symbols-fn] We have been using ``symbols`` instead of ``Symbol`` because it
  automatically splits apart strings into multiple ``Symbol``\ s.
  ``symbols('x y z')`` returns a tuple of three ``Symbol``\ s.  ``Symbol('x y
  z')`` returns a single ``Symbol`` called ``x y z``.
.. [#sympify-fn] Technically, it is an internal function called ``_sympify``,
  which differs from ``sympify`` in that it does not convert strings.  ``x +
  '2'`` is not allowed.
.. [#singleton-fn] Classes like ``One`` and ``Zero`` are singletonized, meaning
  that only one object is ever created, no matter how many times the class is
  called.  This is done for space efficiency, as these classes are very
  common.  For example, ``Zero`` might occur very often in a sparse matrix
  represented densely.  As we have seen, ``NegativeOne`` occurs any time we
  have ``-x`` or ``1/x``.  It is also done for speed efficiency because
  singletonized objects can be compared by ``is``.  The unique objects for
  each singletonized class can be accessed from the ``S`` object.
