# Writing Custom Functions

This guide will describe how to create custom function classes in SymPy.
Custom functions use the same mechanisms as the {ref}`functions <functions>`
that are included with SymPy such as the common {ref}`elementary functions
<elementary-functions>` like {class}`~.exp()` or {class}`~.sin`, {ref}`special
functions <special-functions>` like {class}`~.gamma()` or {class}`~.Si()`, and
{ref}`combinatorial functions <combinatorial-functions>` and {ref}`number
theory functions <ntheory-module>` like {class}`~.factorial()` or
{func}`~.primepi()`. Thus custom user functions are able to do the same things
as functions defined by SymPy.

## Easy Cases: Fully Symbolic or Fully Evaluated

Before digging into the more advanced functionality for custom functions, two
common cases are much simpler.

### Fully Symbolic

If your function `f` has no mathematical properties you
want to define on it, and should never evaluate on any arguments, you can
create an undefined function using `Function('f')`

```
>>> from sympy import symbols, Function
>>> x = symbols('x')
>>> f = Function('f')
>>> f(x)
f(x)
>>> f(0)
f(0)
```

This is useful, for instance, when solving [ODEs](ode-docs). This is also
useful if you only wish to create a symbol that depends on another symbol for
the purposes of differentiation. By default, SymPy assumes all symbols are
independent of one another:

```
>>> from sympy.abc import x, y
>>> y.diff(x)
0
```

To make a symbol depend on another symbol, you should make it a function that
explicitly depends on that symbol.

```
>>> y = Function('y')
>>> y(x).diff(x)
Derivative(y(x), x)
```

If you want your function to have additional behavior, for example, to have a
custom derivative, or to evaluate on certain arguments, you should create a
custom `Function` subclass as described below. However, undefined functions do
support one additional feature, which is that assumptions can be defined on
them, using the same syntax as used by symbols. This defines the assumptions on the output of the function, not the
inputs.

```
>>> g = Function('g', real=True)
>>> g(x)
g(x)
>>> g(x).is_real
True
```

To make a function's assumptions depend on its input in some way, you should
create a custom `Function` subclass as described below.

### Fully Evaluated

At the other end of the spectrum are functions that always evaluate to
something no matter what their inputs are. These functions are never left in
an unevaluated, symbolic form like `f(x)`.

In this case, you should just use a normal Python function using the `def`
keyword.

```
>>> def f(x):
...     if x == 0:
...         return 0
...     else:
...         return x + 1
>>> f(0)
0
>>> f(1)
2
>>> f(x)
x + 1
```

If you find yourself defining an `eval` method on a `Function` subclass where
you always return a value and never return `None`, you should consider just
using a normal Python function instead, as there is no benefit to using a
symbolic `Function` subclass in that case.

Note that in many cases, functions like these can be represented directly
using SymPy classes. For example, the above function can be represented
symbolically using `Piecewise`. The function can then be evaluated using
{meth}`subs() <sympy.core.basic.Basic.subs>`.

```
>>> from sympy import Piecewise, Eq, pprint
>>> f = Piecewise((0, Eq(x, 0)), (x + 1, True))
>>> pprint(f, use_unicode=True)
⎧  0    for x = 0
⎨
⎩x + 1  otherwise
>>> f.subs(x, 0)
0
>>> f.subs(x, 1)
2
```

This is often preferred as symbolic representations like Piecewise can
accurately represent symbolic values. For example, the above Python `def`
definition of `f` implicitly assumes that symbolic `x` is nonzero, but the
Piecewise version handles this case correctly.

Ultimately, the correct tool for the job depends on what you are doing and
what exact behavior you want.
