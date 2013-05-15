================
 Simplification
================

``simplify``
============

Now that we know how to create Symbols, let's jump in and do some interesting
mathematics.  One of the most useful features of a symbolic manipulation
system is the ability to simplify mathematical expressions.  SymPy has dozens
of functions to perform various kinds of simplification.  There is also one
general function called ``simplify`` that attempts to apply all of these
functions in an intelligent way to arrive at the simplest form of an
expression.  Here are some examples

    >>> simplify(sin(x)**2 + cos(x)**2)
    1
    >>> simplify((x**3 + x**2 - x - 1)/(x**2 + 2*x + 1))
    x - 1
    >>> simplify(gamma(x)/gamma(x - 1))
    (x - 2)*(x - 1)

Here, ``gamma(x)`` is `\Gamma(x)`, the gamma function.  We see that 
