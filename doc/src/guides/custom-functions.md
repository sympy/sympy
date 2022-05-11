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

This guide describes how to define functions that map a subset of
$\mathbb{C}^n$ to $\mathbb{C}$. Functions that accept or return other kinds of
objects should subclass another class, such as {class}`~.Boolean`,
{class}`~.MatrixExpr`, {class}`~.Expr`, or {class}`~.Basic`. Much of what is
written here only applies to {class}`~.Function` subclasses and will not work
for general {class}`~.Basic` or {class}`~.Expr` subclasses.

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
versin(pi*x)
```

Note that in this example, we make use of the fact that if a Python function
does not explicitly return a value, it automatically returns `None`, so in the
cases where the `if` statement is not triggered, `eval` returns `None` and
`versin` remains unevaluated.

`eval` can take any number of arguments, including an arbitrary number with
`*args` and optional keyword arguments. The `.args` of the function will
always be the arguments that were passed in by the user. For example

```py
>>> class f(Function):
...     @classmethod
...     def eval(cls, x, y=1, *args):
...         return None
>>> f(1).args
(1,)
>>> f(1, 2).args
(1, 2)
>>> f(1, 2, 3).args
(1, 2, 3)
```

#### Best Practices for `eval`

- **Don't just return an expression.**

  In the above example, we might have been tempted to write

  ```py
  >>> from sympy import cos
  >>> class versin(Function):
  ...    @classmethod
  ...    def eval(cls, x):
  ...        return 1 - cos(x)
  ```

  However, this would make it so that `versin(x)` would *always* return `1 -
  cos(x)`, regardless of what `x` is. If all you want is a quick shorthand to
  `1 - cos(x)`, that is fine. But would be much simpler and more explicit to
  just use a Python function, as [described
  above](custom-functions-fully-evaluated). If we defined `versin` like this,
  it would never actually return as `versin(x)`, and none of the other
  behavior we define below would matter, because the other behaviors we are
  going to define on the `versin` class only apply when the returned object is
  actually a `versin` instance. So for example, `versin(x).diff(x)` would
  actually just be `(1 - cos(x)).diff(x)`, instead of calling our derivative
  we defined below. Of course, for a function as simple as `versin`, this
  might seem reasonable, but remember that we just chose it as an example.
  Again, you should think about what you actually want to achieve and whether
  you just want a shorthand function for an expression or a symbolic function.

- **Avoid too much automatic evaluation.**

  It is recommended to minimize what is evaluated automatically by
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

- **Allow `None` input assumptions.**

  Our example function $\operatorname{versin}$ is a function from $\mathbb{C}$
  to $\mathbb{C}$, so it can accept any input. But suppose we had a function
  that only made sense with certain inputs. As a second example, let's define
  a function `divides` as

  (custom-functions-divides-example)=

  $$\operatorname{divides}(m, n) = \begin{cases} 1 & \text{for}\: m \mid n \\
  0 & \text{for}\: m\not\mid n  \end{cases}.$$

  That is, `divides(m, n)` will be `1` if `m` divides `n` and `0` otherwise.
  `divides` clearly only makes sense if `m` and `n` are integers.

  We might want to define the `eval` method for `divides` like this.

  ```py
  >>> class divides(Function):
  ...     @classmethod
  ...     def eval(cls, m, n):
  ...         # Evaluate for explicit integer m and n. This part is fine.
  ...         if isinstance(m, Integer) and isinstance(n, Integer):
  ...             return int(n % m != 0)
  ...         # For symbolic arguments, require m and n to be integer.
  ...         # If we write the logic this way, we will run into trouble.
  ...         if not m.is_integer or not n.is_integer:
  ...             raise TypeError("m and n should be integers")
  ```

  The problem here is that by using `if not m.is_integer`, we are requiring
  `m.is_integer` to be `True`. If it is `None`, it will fail (see also the
  [guide on booleans and three-valued logic](booleans-guide)). This is
  problematic for two reasons. Firstly, it forces the user to define
  assumptions on any input variable. If a user omits them, it will fail:

  ```
  >>> n, m = symbols('n m')
  >>> divides(m, n)
  Traceback (most recent call last):
  ...
  TypeError: m and n should be integers
  ```

  Instead they have to write

  ```
  >>> n, m = symbols('n m', integer=True)
  >>> divides(m, n)
  divides(m, n)
  ```

  This may seem like an acceptable restriction, but there is a bigger problem.
  Sometimes, SymPy's assumptions system cannot deduce an assumption, even
  though it is mathematically true. In this case, it will give `None` (`None`
  is used to mean both "unknown" and "cannot compute" in SymPy's assumptions).
  For example

  ```
  >>> # n and m are still defined as integer=True as above
  >>> divides(2, (m**2 + m)/2)
  Traceback (most recent call last):
  ...
  TypeError: m and n should be integers
  ```

  Here the expression `(m**2 + m)/2` is always an integer, but SymPy's
  assumptions system is not able to deduce this

  ```
  >>> print(((m**2 + m)/2).is_integer)
  None
  ```

  Consequently, one should always test *negated* assumptions for input variables, that is,
  fail if the assumption is `False` but allow the assumption to be `None`.


  ```py
  >>> class divides(Function):
  ...     @classmethod
  ...     def eval(cls, m, n):
  ...         # Evaluate for explicit integer m and n. This part is fine.
  ...         if isinstance(m, Integer) and isinstance(n, Integer):
  ...             return int(n % m != 0)
  ...         # For symbolic arguments, require m and n to be integer.
  ...         # This is the better way to write this logic.
  ...         if m.is_integer is False or n.is_integer is False:
  ...             raise TypeError("m and n should be integers")
  ```

  This still disallows non-integer inputs as desired:

  ```
  >>> divides(1.5, 1)
  Traceback (most recent call last):
  ...
  TypeError: m and n should be integers
  ```

  But it does not fail in cases where the assumption is `None`:

  ```
  >>> divides(2, (m**2 + m)/2)
  divides(2, m**2/2 + m/2)
  >>> _.subs(m, 2)
  1
  >>> x, y = symbols('x y')
  >>> divides(y, x)
  divides(y, x)
  ```

### Assumptions

The next thing we might want to define are the assumptions on our function.
The [guide on the assumptions system](assumptions-guide) goes into the
assumptions system in great detail. Here we will show the basic ways we can
define the assumptions of a custom function.

The simplest case is if a function always has a given assumption regardless of
its input. For example, the `divides` example we [showed
above](custom-functions-divides-example) is always an integer, because its value is always
either 0 or 1. In this case, you can define <span
class="pre">is_*assumption*</span> directly on the class.

```py
>>> class divides(Function):
...     is_integer = True
>>> divides(m, n).is_integer
True
```

```{note}
From here on out in this guide, in the interest of space, we will omit the
previous method definitions in the examples unless they are needed for the
given example to work. A [complete example](custom-functions-complete-example)
with all methods is given below.
```

In general, however, the assumptions of a function depend on the assumptions
of its inputs. In this case, you should define an <span class="pre">_eval_*assumption*</span>
method.

For our $\operatorname{versin}$ function, the function is always between 0 and
2 for real valued inputs, and it is 0 whenever the input is an even multiple
of π. So we `versin(x)` should be *nonnegative* whenever `x` is *real*
(remember that unless you restricted the assumptions of your function, a
Function is assumed to be complex valued with complex valued inputs) and
*positive* whenever `x` is not an even multiple of π.

In the method, as in all methods, we can access the arguments of the function
using `.args`.

```py
>>> from sympy.core.logic import fuzzy_and, fuzzy_not
>>> class versin(Function):
...     def _eval_is_nonnegative(self):
...         x = self.args[0]
...         return x.is_real
...
...     def _eval_is_positive(self):
...         x = self.args[0]
...         return fuzzy_and([x.is_real, fuzzy_not((x/pi).is_even)])
>>> versin(1).is_nonnegative
True
```

Note the use of `fuzzy_` functions in the more complicated `_eval_is_positive`
handler. It is important when working with assumptions to always be careful
about handling [three-valued logic](booleans-guide) correctly.

Note that it is not necessary to define `is_real` because those are deduced
automatically from the other assumptions, since `nonnegative -> real`. In
general, you should avoid defining assumptions that the assumptions system can
deduce automatically given its [known facts](assumptions-guide-predicates).

```py
>>> versin(1).is_real
True
```

The assumptions system is sometimes able to deduce more than you might think.
For example, from the above, it can deduce that `versin(2*n*pi)` is zero when
`n` is an integer.

```py
>>> n = symbols('n', integer=True)
>>> versin(2*n*pi).is_zero
True
```

It's worth checking if the assumptions system can deduce something
automatically before manually coding it.

Finally, a word of warning to be careful about correctness when coding
assumptions. Be careful to use the exact
[definitions](assumptions-guide-predicates) of the various assumptions, and be
careful to handle `None` cases correctly. Incorrect or inconsistent
assumptions can lead to subtle bugs.

### Numerical Evaluation with `evalf`

### Rewriting and Simplification

### Differentiation

### Series Expansions

### Printing

(custom-functions-complete-example)=
## Complete Example
