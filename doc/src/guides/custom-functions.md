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

Before digging into the more advanced functionality for custom functions, we
should mention two common cases that are much simpler, the case where the
function is fully symbolic, and the case where the function is fully
evaluated.

(custom-functions-fully-symbolic)=
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

This is useful, for instance, when solving [ODEs](ode-docs).

This is also useful if you only wish to create a symbol that depends on
another symbol for the purposes of differentiation. By default, SymPy assumes
all symbols are independent of one another:

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
them, using the same syntax as used by symbols. This defines the assumptions
of the output of the function, not the input (that is, it defines the
function's range, not its domain).

```
>>> g = Function('g', real=True)
>>> g(x)
g(x)
>>> g(x).is_real
True
```

To make a function's assumptions depend on its input in some way, you should
create a custom `Function` subclass as described below.

(custom-functions-fully-evaluated)=
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

Fully symbolic representations like Piecewise have the advantage that they
accurately represent symbolic values. For example, the above Python `def`
definition of `f`, `f(x)` would implicitly assume that `x` is nonzero. The
Piecewise version handles this case correctly.

Ultimately, the correct tool for the job depends on what you are doing and
what exact behavior you want.

## Creating a custom function

The basic way to create a custom function is to subclass {class}`~.Function`.
The name of the subclass will be the name of the function. Different methods
should be defined on the function depending on what functionality we want to
provide.

As a motivating example for this document, let's create a custom function
representing the [versine](https://en.wikipedia.org/wiki/Versine) function.
Versine is a trigonometric function which was used historically alongside some
of the more familiar trigonometric functions like sine and cosine. It is
rarely used today. Versine can be defined by the identity

$$\operatorname{versin}(x) = 1 - \cos(x).$$

Because versine is used so rarely in modern mathematics and because it is so
easily defined in terms of the more familiar cosine, it does not already have
a definition in SymPy, but this simple definition will make it easy for us to
define the various behaviors.

Let us start by subclassing `Function`.

```
>>> class versin(Function):
...     pass
```

At this point, `versin` has no behaviors defined on it. It is very similar to
an [undefined function](custom-functions-fully-symbolic). We note that
`versin` is a class, and `versin(x)` is an instance of this class.

```
>>> type(versin)
<class 'sympy.core.function.FunctionClass'>
>>> versin(x)
versin(x)
>>> isinstance(versin(x), versin)
True
```

### Defining Automatic Evaluation with `eval`


```{sidebar} Reminder
Remember that `eval` should be defined with the `@classmethod` decorator.
```

The first and most common thing we might want to define on our custom function
is to define automatic evaluation, that is, cases where it will return an
actual value instead of just remaining unevaluated as-is.

This is done by defining the class method `eval`. `eval` should take the
arguments of the function and return either a value or `None`. If it returns
`None`, the function will remain unevaluated in that case.

For our function `versin`, we might recall that $\cos(n\pi) = (-1)^n$ for
integer $n$, so $\operatorname{versin}(n\pi) = 1 - (-1)^n.$ We can make
`versin` automatically evaluate to this value when passed an integer multiple
of `pi`.

```
>>> from sympy import pi, Integer
>>> class versin(Function):
...    @classmethod
...    def eval(cls, x):
...        # If x is an integer multiple of pi, x/pi will cancel and be an Integer
...        n = x/pi
...        if isinstance(n, Integer):
...            return 1 - (-1)**n
>>> versin(pi)
2
>>> versin(2*pi)
0
>>> versin(x*pi)
versin(x*pi)
```

Note that in this example, we make use of the fact that if a Python function
does not explicitly return a value, it automatically returns `None`, so in the
cases where the `if` statement is not triggered, `eval` returns `None` and
`versin` remains unevaluated.

Two things are worth noting about evaluation:

- First, we might have been tempted to write

  ```py
  >>> from sympy import cos
  >>> class versin(Function):
  ...    @classmethod
  ...    def eval(cls, x):
  ...        return 1 - cos(x)
  ```

  However, this would make it so that `versin(x)` would *always* return `1 -
  cos(x)`, regardless of what `x` is. If all you want is a quick shorthand to
  `1 - cos(x)`, that is fine, it is much simpler and more explicit to just use
  a Python function, as [described above](custom-functions-fully-evaluated).
  If we defined `versin` like this, it would never actually return as
  `versin(x)`, and none of the other behavior we define below would matter,
  because the other behaviors we are going to define on the `versin` class
  only apply when the returned object is actually a `versin` instance. So for
  example, `versin(x).diff(x)` would actually just be `(1 - cos(x)).diff(x)`,
  instead of calling our derivative we defined below. Of course, for a
  function as simple as `versin`, this might not seem that bad, but remember
  that we just chose it as an example.

- Secondly, it is recommended to minimize what is evaluated automatically by
  `eval`. It is typically better to put more advanced simplifications in other
  methods like `doit` or `simplify` (see below). Remember that whatever we
  define for automatic evaluation will always evaluate. It is possible to skip
  automatic evaluation by using `evaluate=False`, but this is fragile, and it
  tends to be bug prone because other code may be written expecting the
  automatic evaluation to occur. As in the previous point, if we evaluate
  every value, there is little point to even having a symbolic function in the
  first place. For example, we might be tempted to evaluate some trig
  identities on `versin` in `eval`, but then these identities would always
  evaluate, and it wouldn't be possible to represent one half of the identity.

  One should also avoid doing anything in `eval` that is slow to compute.
  SymPy functions generally assume that it is cheap to create expressions, and
  if this is not true, it can lead to performance issues everywhere.

  You might have noticed that we used `isinstance(n, Integer)` instead of
  checking `n.is_integer` using the assumptions system. We could have done
  that instead, which would make `versin(n*pi)` evaluate even if `n =
  Symbol('n', integer=True)`. But this is a case where we might not always
  want evaluation to happen, and if `n` is a more complicated expression,
  `n.is_integer` might take more time to compute. Although SymPy itself
  doesn't always follow this rule, it is recommended to avoid performing
  automatic evaluation in `eval` based on assumptions. Instead, `eval` should
  typically only evaluate explicit numerical special values and return `None`
  for everything else.
