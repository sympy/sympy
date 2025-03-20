Basic Usage
===========

Solve a simple equation:

>>> from sympy import solveset, S
>>> from sympy.abc import x
>>> solveset(x**2 - 1, x)
{-1, 1}

Solve an equation in the real domain:

>>> solveset(x**2 + 1, x, domain=S.Reals)
EmptySet

Solve a transcendental equation:

>>> from sympy import exp, sin
>>> solveset(exp(x) - 1, x)
{2*I*pi*n | n in Integers} ∪ {0}

Solve a trigonometric equation:

>>> solveset(sin(x), x)
{2*pi*n | n in Integers} ∪ {2*pi*n + pi | n in Integers}

Advanced Usage
=============

Solve an equation with parameters:

>>> from sympy.abc import a
>>> solveset(x**2 - a**2, x)
{-a, a}

Solve an equation with multiple variables:

>>> solveset((x + y - 2, x - y), (x, y))
{(1, 1)}

Solve an equation with complex solutions:

>>> solveset(x**2 + 1, x, domain=S.Complexes)
{-I, I}

Handle equations with no solutions:

>>> solveset(x + 1, x, domain={2})
EmptySet

Handle equations with infinitely many solutions:

>>> solveset(sin(x) - sin(x), x, domain=S.Reals)
Reals 