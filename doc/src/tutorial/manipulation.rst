==================================
 Advanced Expression Manipulation
==================================

In this section, we discuss some ways that we can perform advanced
manipulation of expressions.

Understanding Expression Trees
==============================

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

.. This comes from dotprint(x**2 + x*y)

.. graphviz::

   digraph{

   # Graph style
   "rankdir"="TD"

   #########
   # Nodes #
   #########

   "Symbol(y)" ["color"="black", "label"="Symbol('y')", "shape"="ellipse"];
   "Symbol(x)1" ["color"="black", "label"="Symbol('x')", "shape"="ellipse"];
   "Symbol(x)2" ["color"="black", "label"="Symbol('x')", "shape"="ellipse"];
   "Integer(2)" ["color"="black", "label"="Integer(2)", "shape"="ellipse"];
   "Mul(Symbol(x)1, Symbol(y))" ["color"="black", "label"="Mul", "shape"="ellipse"];
   "Pow(Symbol(x)2, Integer(2))" ["color"="black", "label"="Pow", "shape"="ellipse"];
   "Add(Mul(Symbol(x)1, Symbol(y)), Pow(Symbol(x)2, Integer(2)))" ["color"="black", "label"="Add", "shape"="ellipse"];

   #########
   # Edges #
   #########

   "Mul(Symbol(x)1, Symbol(y))" -> "Symbol(x)1";
   "Mul(Symbol(x)1, Symbol(y))" -> "Symbol(y)";
   "Pow(Symbol(x)2, Integer(2))" -> "Symbol(x)2";
   "Pow(Symbol(x)2, Integer(2))" -> "Integer(2)";
   "Add(Mul(Symbol(x)1, Symbol(y)), Pow(Symbol(x)2, Integer(2)))" -> "Pow(Symbol(x)2, Integer(2))";
   "Add(Mul(Symbol(x)1, Symbol(y)), Pow(Symbol(x)2, Integer(2)))" -> "Mul(Symbol(x)1, Symbol(y))";
   }

.. note::

   The above diagram was made using `Graphviz <http://www.graphviz.org/>`_ and
   the ``dotprint`` function, in ``sympy.printing.dot``.

First, let's look at the leaves of this tree.  Symbols are instances of class
Symbol.  While we have been doing

    >>> x = symbols('x')

We could have also done

    >>> x = Symbol('x')

Either way, we get a Symbol with the name "x".  For the number in the
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
or operation, the non-SymPy object will be converted into a SymPy object.

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
``Add(Pow(x, 2), Mul(x, y)``.

    >>> Add(Pow(x, 2), Mul(x, y))
    x**2 + x*y

A SymPy expression tree might have many branches, and could be quite deep or
quite broad.  Here is a more complicated example

    >>> expr = sin(x*y)/2 - x**2 + 1/y
    >>> srepr(expr)
    "Add(Mul(Integer(-1), Pow(Symbol('x'), Integer(2))), Mul(Rational(1, 2),
    Function('sin')(Mul(Symbol('x'), Symbol('y')))), Pow(Symbol('y'),
    Integer(-1)))"

Here is a diagram

.. graphviz::

    digraph{

    # Graph style
    "rankdir"="TD"

    #########
    # Nodes #
    #########

    "Half()" ["color"="black", "label"="Rational(1, 2)", "shape"="ellipse"];
    "Symbol(y)" ["color"="black", "label"="Symbol('y')", "shape"="ellipse"];
    "Symbol(x)" ["color"="black", "label"="Symbol('x')", "shape"="ellipse"];
    "Integer(2)" ["color"="black", "label"="Integer(2)", "shape"="ellipse"];
    "Symbol(x1)" ["color"="black", "label"="Symbol('x')", "shape"="ellipse"];
    "Symbol(y1)" ["color"="black", "label"="Symbol('y')", "shape"="ellipse"];
    "NegativeOne()p" ["color"="black", "label"="Integer(-1)", "shape"="ellipse"];
    "NegativeOne()c" ["color"="black", "label"="Integer(-1)", "shape"="ellipse"];
    "Pow(Symbol(x), Integer(2))" ["color"="black", "label"="Pow", "shape"="ellipse"];
    "Mul(Symbol(x1), Symbol(y1))" ["color"="black", "label"="Mul", "shape"="ellipse"];
    "Pow(Symbol(y), NegativeOne()p)" ["color"="black", "label"="Pow", "shape"="ellipse"];
    "sin(Mul(Symbol(x1), Symbol(y1)))" ["color"="black", "label"="sin", "shape"="ellipse"];
    "Mul(Half(), sin(Mul(Symbol(x1), Symbol(y1))))" ["color"="black", "label"="Mul", "shape"="ellipse"];
    "Mul(NegativeOne()c, Pow(Symbol(x), Integer(2)))" ["color"="black", "label"="Mul", "shape"="ellipse"];
    "Add(Mul(Half(), sin(Mul(Symbol(x1), Symbol(y1)))), Mul(NegativeOne()c, Pow(Symbol(x), Integer(2))), Pow(Symbol(y), NegativeOne()p))" ["color"="black", "label"="Add", "shape"="ellipse"];

    #########
    # Edges #
    #########

    "Pow(Symbol(x), Integer(2))" -> "Symbol(x)";
    "Pow(Symbol(x), Integer(2))" -> "Integer(2)";
    "Mul(Symbol(x1), Symbol(y1))" -> "Symbol(x1)";
    "Mul(Symbol(x1), Symbol(y1))" -> "Symbol(y1)";
    "Pow(Symbol(y), NegativeOne()p)" -> "Symbol(y)";
    "Pow(Symbol(y), NegativeOne()p)" -> "NegativeOne()p";
    "Mul(Half(), sin(Mul(Symbol(x1), Symbol(y1))))" -> "Half()";
    "Mul(NegativeOne()c, Pow(Symbol(x), Integer(2)))" -> "NegativeOne()c";
    "sin(Mul(Symbol(x1), Symbol(y1)))" -> "Mul(Symbol(x1), Symbol(y1))";
    "Mul(NegativeOne()c, Pow(Symbol(x), Integer(2)))" -> "Pow(Symbol(x), Integer(2))";
    "Mul(Half(), sin(Mul(Symbol(x1), Symbol(y1))))" -> "sin(Mul(Symbol(x1), Symbol(y1)))";
    "Add(Mul(Half(), sin(Mul(Symbol(x1), Symbol(y1)))), Mul(NegativeOne()c, Pow(Symbol(x), Integer(2))), Pow(Symbol(y), NegativeOne()p))" -> "Pow(Symbol(y), NegativeOne()p)";
    "Add(Mul(Half(), sin(Mul(Symbol(x1), Symbol(y1)))), Mul(NegativeOne()c, Pow(Symbol(x), Integer(2))), Pow(Symbol(y), NegativeOne()p))" -> "Mul(Half(), sin(Mul(Symbol(x1), Symbol(y1))))";
    "Add(Mul(Half(), sin(Mul(Symbol(x1), Symbol(y1)))), Mul(NegativeOne()c, Pow(Symbol(x), Integer(2))), Pow(Symbol(y), NegativeOne()p))" -> "Mul(NegativeOne()c, Pow(Symbol(x), Integer(2)))";
    }

This expression reveals some interesting things about SymPy. Let's go through
them one by one.

Let's first look at the term ``-x**2``.  As we expected, we see ``Pow(x,
2)``.  One level up, we see we have ``Mul(-1, Pow(x, 2))``.  There is no
subtraction class in SymPy.  ``x - y`` is represented as ``x + -y``, or, more
completely, ``x + -1*y``, i.e., ``Add(x, Mul(-1, y))``.

Next, look at ``1/y``.  We might expect to see something like ``Div(1, y)``,
but similar to subtraction, there is no class in SymPy for division.  Rather,
division is represented by a power of -1.  Hence, we have ``Pow(y, -1)``.
What if we had divided something other than 1 by ``y``, like ``x/y``?  Let's
see.

    >>> expr = x/y
    >>> srepr(expr)
    "Mul(Symbol('x'), Pow(Symbol('y'), Integer(-1)))"

.. graphviz::

   digraph{

   # Graph style
   "rankdir"="TD"

   #########
   # Nodes #
   #########

   "Symbol(x)" ["color"="black", "label"="Symbol('x')", "shape"="ellipse"];
   "Symbol(y)" ["color"="black", "label"="Symbol('y')", "shape"="ellipse"];
   "NegativeOne()" ["color"="black", "label"="Integer(-1)", "shape"="ellipse"];
   "Pow(Symbol(y), NegativeOne())" ["color"="black", "label"="Pow", "shape"="ellipse"];
   "Mul(Symbol(x), Pow(Symbol(y), NegativeOne()))" ["color"="black", "label"="Mul", "shape"="ellipse"];

   #########
   # Edges #
   #########

   "Pow(Symbol(y), NegativeOne())" -> "Symbol(y)";
   "Pow(Symbol(y), NegativeOne())" -> "NegativeOne()";
   "Mul(Symbol(x), Pow(Symbol(y), NegativeOne()))" -> "Symbol(x)";
   "Mul(Symbol(x), Pow(Symbol(y), NegativeOne()))" -> "Pow(Symbol(y), NegativeOne())";
   }

We see that ``x/y`` is represented as ``x*(1/y)``, which ends up as ``Mul(x,
Pow(y, -1))``.

Finally, let's look at the ``sin(x*y)/2`` term.  Following the patterns of
before, we might expect to see ``Mul(sin(x*y), Pow(Integer(2), -1))``.  But
instead, we have ``Mul(Rational(1, 2), sin(x*y))``.  Rational numbers are
always combined into a single term in a multiplication, so that when we divide
by 2, it is represented as multiplying by 1/2.

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

Recursing through an Expression Tree
====================================

Now that you know how expression trees work in SymPy, let's look at how to dig
our way through an expression tree.  Every object in SymPy has two very
important attributes, ``func``, and ``args``.  ``func`` is the head of the
object.  Usually, it is the same as the class of the object, though there are
exceptions to this rule.  For example, ``(x*y).func`` is ``Mul``.  ``args`` are
the top-level arguments of the object.  ``(x*y).args`` would be ``(x, y)``.
Let's look at some examples

    >>> expr = 3*y**2*x
    >>> expr.func
    <class 'sympy.core.mul.Mul'>
    >>> expr.args
    (3, x, y**2)

From this, we can see that ``expr == Mul(2, y**2, x)``.  Note that although we
entered ``2*y**2*x``, the ``args`` are ``(2, x, y**2)``.  In a ``Mul``, the
Rational coefficient will come first in the ``args``, but other than that, the
order of everything else follows no special pattern.  To be sure, though,
there is an order.

    >>> expr = y**2*3*x
    >>> expr.args
    (3, x, y**2)

Mul's ``args`` are sorted, so that the same ``Mul`` will have the same
``args``.  But the sorting is based on several criteria that has no
mathematical significance.

Recall that our ``expr`` should be ``Mul(2, x, Pow(y, 2))``.  What if we want
to get at the ``args`` of ``Pow(y, 2)``.  Notice that the ``y**2`` is in the
third slot of ``expr.args``, i.e., ``expr.args[2]``.

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
hit the bottom of the expression tree.

With this knowledge, let's look at how we can recurse through an expression
tree.  The nested nature of ``args`` is a perfect fit for recursive
functions.  The base case will be empty ``args``.  Let's write a simple
function that goes through an expression and prints all the ``args``.

    >>> def pre(expr):
    ...     print expr
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
    ...     print arg
    x*y + 1
    1
    x*y
    x
    y
