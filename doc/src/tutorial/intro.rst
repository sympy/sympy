==============
 Introduction
==============

What is Symbolic Computation?
=============================

Symbolic computation deals with the computation and manipulation of
mathematical objects symbolically.  This means that the mathematical objects
are represented exactly, not approximately, and mathematical expressions with
unevaluated variables are left in symbolic form.

Let's take an example. Say we wanted to use the built-in Python functions to
compute square roots. We might do something like this

::

   >>> impot math
   >>> math.sqrt(9)
   3.0

9 is a perfect square, so we got the exact answer, 3. But suppose we computed
the square root of a number that isn't a perfect square

::

   >>> math.sqrt(8)
   2.8284271247461903

Here we got an approximate result. 2.8284271247461903 is not the exact square
root of 8 (indeed, the actual square root of 8 cannot be represented by a
finite decimal, since it is an irrational number).  If all we cared about was
the decimal form of the square root of 8, we would be done.

But suppose we want to go further. Recall that `\sqrt{8} = \sqrt{4\dot 2} =
2\sqrt{2}`.  We would have a hard time deducing this from the above result.
This is where symbolic computation comes in.  With a symbolic computation
system like SymPy, square roots of numbers that are not perfect squares are
left unevaluated by default

::

   >>> import sympy
   >>> sympy.sqrt(3)
   sqrt(3)

Furthermore---and this is where we start to see the real power of symbolic
computation---symbolic results can be symbolically simplified.

::

   >>> sympy.sqrt(8)
   2*sqrt(2)

A More Interesting Example
==========================

The above example starts to show how we can compute with exact numbers using
SymPy.  But it is much more powerful than that.  Symbolic computation systems
(which by the way, are also often called computer algebra systems, or just
CASs) such as SymPy are capable of computing with symbolic expressions with
variables.
