===
Map
===

.. module:: sympy.map

Introduction
------------

This module implements mathematical function as ``Basic`` object.

SymPy's ``Function`` is type object, which hinders itself being treated as mathematical object.
Using ``Map`` in this module instead, several useful features such as operations between
functions or derivative on pure function (not on expression) can be achieved.  
The ultimate goal of this module is to entirely replace ``core.function`` and ``functions``
module.

Function vs Expression
----------------------

Frequently, function $f: X \rightarrow Y$ is written as $f(x)$. Here, $f(x)$ is not the
function itself - it's an expression which represents $f$. It implies that the domain of
$f$ is one-dimensional line, and $x$ is an arbitrary point on it.  
The difference is clear when you are applying a fixed value - say, $3$ - on $f$. With
function $f$, you apply $3$ to $f$; in SymPy, it's ``f(3)``. But with expression, you have to
substitue $x$ with $3$ in $f(x)$; in SymPy, it's ``f(x).subs(x, 3)``.

Taking derivative is another important case where the difference between function and expression
is prominent. With expression, you differentiate $f(x)$ with respect to $x$, which can be
written as $D_x{\left( f(x) \right)}$. But with function, you don't differentiate $f$ with respect
to $x$. Rather, you are constructing a derivative function, which is noted with prime notation $f'$.
With expression, value of $f$'s derivative on $x=3$ is $D_x{\left( f(x) \right)}|_{\substact{x=3}}$.
With function, it's $f'(3)$.

Prime notation is only valid for unary function. For n-ary function, we often encounter the notation
$\partial_{x} f$ for $f: X \times Y \rightarrow Z$. However, this is possible only when there is an
agreement on the relation between a symbol and a set. More general notation for partial derivative
is to denote with the index of argument - such as $\partial_{1} f$ or $\partial_{2} f$.