 (custom-functions)=
# Writing Custom Functions

<!-- Note to contributors: if you update one of the examples in this guide, be
     sure to also update the complete example at the end. -->

This guide will describe how to create custom function classes in SymPy.
Custom functions use the same mechanisms as the {ref}`functions <functions>`
that are included with SymPy such as the common {ref}`elementary functions
<elementary-functions>` like {class}`~.exp()` or {class}`~.sin()`,
{ref}`special functions <special-functions>` like {class}`~.gamma()` or
{class}`~.Si()`, and {ref}`combinatorial functions <combinatorial-functions>`
and {ref}`number theory functions <ntheory-module>` like
{class}`~.factorial()` or {func}`~.primepi()`.

This guide describes how to define functions that map a subset of
$\mathbb{C}^n$ to $\mathbb{C}$. Functions that accept or return other kinds of
objects should subclass another class, such as {class}`~.Boolean`,
{class}`~.MatrixExpr`, {class}`~.Expr`, or {class}`~.Basic`. Much of what is
written here only applies to {class}`~.Function` subclasses and will not work
for general {class}`~.Basic` or {class}`~.Expr` subclasses.

## Easy Cases: Fully Symbolic or Fully Evaluated

Before digging into the more advanced functionality for custom functions, we
should mention two common cases, the case where the function is fully
symbolic, and the case where the function is fully evaluated. Both of these
cases have much simpler alternatives than the full mechanisms described in ths
guide.

(custom-functions-fully-symbolic)=
### The Fully Symbolic Case

If your function `f` has no mathematical properties you
want to define on it, and should never evaluate on any arguments, you can
create an undefined function using `Function('f')`

```py
>>> from sympy import symbols, Function
>>> x = symbols('x')
>>> f = Function('f')
```

```py
>>> f(x)
f(x)
>>> f(0)
f(0)
```

This is useful, for instance, when solving [ODEs](ode-docs).

This is also useful if you only wish to create a symbol that depends on
another symbol for the purposes of differentiation. By default, SymPy assumes
all symbols are independent of one another:

```py
>>> from sympy.abc import x, y
>>> y.diff(x)
0
```

To make a symbol that depends on another symbol, you can use a function that
explicitly depends on that symbol.

```py
>>> y = Function('y')
>>> y(x).diff(x)
Derivative(y(x), x)
```

If you want your function to have additional behavior, for example, to have a
custom derivative, or to evaluate on certain arguments, you should create a
custom `Function` subclass as [described
below](custom-functions-function-subclass). However, undefined functions do
support one additional feature, which is that assumptions can be defined on
them, using the same syntax as used by symbols. This defines the assumptions
of the output of the function, not the input (that is, it defines the
function's range, not its domain).

```py
>>> g = Function('g', real=True)
```

```py
>>> g(x)
g(x)
>>> g(x).is_real
True
```

To make a function's assumptions depend on its input in some way, you should
create a custom `Function` subclass and define assumptions handlers as
[described below](custom-functions-assumptions).

(custom-functions-fully-evaluated)=
### The Fully Evaluated Case

At the other end of the spectrum are functions that always evaluate to
something no matter what their inputs are. These functions are never left in
an unevaluated, symbolic form like `f(x)`.

In this case, you should just use a normal Python function using the `def`
keyword.

```py
>>> def f(x):
...     if x == 0:
...         return 0
...     else:
...         return x + 1
```

```py
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
symbolic `Function` subclass in that case (see the
{any}`custom-functions-eval-best-practices` section below)

Note that in many cases, functions like these can be represented directly
using SymPy classes. For example, the above function can be represented
symbolically using {class}`~.Piecewise`. The `Piecewise` expression can be
evaluated for specific values of `x` using {meth}`subs()
<sympy.core.basic.Basic.subs>`.

```py
>>> from sympy import Piecewise, Eq, pprint
>>> f = Piecewise((0, Eq(x, 0)), (x + 1, True))
```

```py
>>> pprint(f, use_unicode=True)
⎧  0    for x = 0
⎨
⎩x + 1  otherwise
>>> f.subs(x, 0)
0
>>> f.subs(x, 1)
2
```

Fully symbolic representations like `Piecewise` have the advantage that they
accurately represent symbolic values. For example, the above Python `def`
definition of `f`, `f(x)` implicitly assumes that `x` is nonzero. The
`Piecewise` version handles this case correctly and won't evaluate to the $x
\neq 0$ case unless `x` is known to not be zero.

Ultimately, the correct tool for the job depends on what you are doing and
what exact behavior you want.

(custom-functions-function-subclass)=
## Creating a custom function

The basic way to create a custom function is to subclass {class}`~.Function`.
The name of the subclass will be the name of the function. Different methods
should then e defined on this subclass, depending on what functionality you
want to provide.

As a motivating example for this document, let's create a custom function
class representing the [versine
function](https://en.wikipedia.org/wiki/Versine). Versine is a trigonometric
function which was used historically alongside some of the more familiar
trigonometric functions like sine and cosine. It is rarely used today. Versine
can be defined by the identity

$$\operatorname{versin}(x) = 1 - \cos(x).$$

SymPy does not already include versine because it is used so rarely in modern
mathematics and because it is so easily defined in terms of the more familiar
cosine.

Let us start by subclassing `Function`.

```py
>>> class versin(Function):
...     pass
```

At this point, `versin` has no behaviors defined on it. It is very similar to
the [undefined functions](custom-functions-fully-symbolic) we discussed above.
We note that `versin` is a class, and `versin(x)` is an instance of this
class.

```py
>>> type(versin)
<class 'sympy.core.function.FunctionClass'>
>>> versin(x)
versin(x)
>>> isinstance(versin(x), versin)
True
```

```{note}

All the methods described below are optional. They can be included if you want
to define the given behavior, but if they are omitted, SymPy will default to
leaving things unevaluated. For example, if you do not define
[differentiation](custom-functions-differentiation), {func}`~.diff` will just
return an unevaluated {class}`~.Derivative`.

```

(custom-functions-eval)=
### Defining Automatic Evaluation with `eval`

```{sidebar} Reminder
Remember that `eval` should be defined with the `@classmethod` decorator.
```

The first and most common thing we might want to define on our custom function
is automatic evaluation, that is, the cases where it will return an actual
value instead of just remaining unevaluated as-is.

This is done by defining the class method `eval`. `eval` should take the
arguments of the function and return either a value or `None`. If it returns
`None`, the function will remain unevaluated in that case.

For our function `versin`, we might recall that $\cos(n\pi) = (-1)^n$ for
integer $n$, so $\operatorname{versin}(n\pi) = 1 - (-1)^n.$ We can make
`versin` automatically evaluate to this value when passed an integer multiple
of `pi`:

(custom-functions-versin-eval-example)=
```py
>>> from sympy import pi, Integer
>>> class versin(Function):
...    @classmethod
...    def eval(cls, x):
...        # If x is an integer multiple of pi, x/pi will cancel and be an Integer
...        n = x/pi
...        if isinstance(n, Integer):
...            return 1 - (-1)**n
```

```py
>>> versin(pi)
2
>>> versin(2*pi)
0
>>> versin(x*pi)
versin(pi*x)
```

Here we make use of the fact that if a Python function does not explicitly
return a value, it automatically returns `None`, so in the cases where the `if
isinstance(n, Integer)` statement is not triggered, `eval` returns `None` and
`versin` remains unevaluated.

```{note}

`Function` subclasses should not redefine `__new__` or `__init__`. If you want
to implement behavior that isn't possible with `eval`, it might make more
sense to subclass {class}`~.Expr` rather than `Function`.

```

`eval` can take any number of arguments, including an arbitrary number with
`*args` and optional keyword arguments. The `.args` of the function will
always be the arguments that were passed in by the user. For example

```py
>>> class f(Function):
...     @classmethod
...     def eval(cls, x, y=1, *args):
...         return None
```

```py
>>> f(1).args
(1,)
>>> f(1, 2).args
(1, 2)
>>> f(1, 2, 3).args
(1, 2, 3)
```

Finally, note that automatic evaluation on floating-point inputs happens
automatically once `evalf` is defined. See the section on [defining numerical
evaluation with `evalf`](custom-functions-evalf) below.

(custom-functions-eval-best-practices)=
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
  `1 - cos(x)`, that is fine, but would be much simpler and more explicit to
  just [use a Python function as described
  above](custom-functions-fully-evaluated). If we defined `versin` like this,
  it would never actually be represented as `versin(x)`, and none of the other
  behavior we define below would matter, because the other behaviors we are
  going to define on the `versin` class only apply when the returned object is
  actually a `versin` instance. So for example, `versin(x).diff(x)` would
  actually just be `(1 - cos(x)).diff(x)`, instead of calling [the `fdiff`
  method we define below](custom-functions-differentiation).

  **The point of `eval` is not to define what the function *is*,
  mathematically, but rather to specify on what inputs it should automatically
  evaluate.** The mathematical definition of a function is determined through
  specifying various mathematical properties with the methods outlined below,
  like [numerical evaluation](custom-functions-evalf),
  [differentiation](custom-functions-differentiation), and so on.

  If you find yourself doing this, you should think about what you actually
  want to achieve. If you just want a shorthand function for an expression, it
  will be simpler to just [define a Python
  function](custom-functions-fully-evaluated). If you really do want a
  symbolic function, think about when you want it to evaluate to something
  else and when you want it to stay unevaluated.

- **Avoid too much automatic evaluation.**

  It is recommended to minimize what is evaluated automatically by `eval`. It
  is typically better to put more advanced simplifications in [other methods
  like `doit` or `simplify`](custom-functions-rewriting-and-simplification).
  Remember that whatever you define for automatic evaluation will *always*
  evaluate. It is technically possible to bypass automatic evaluation by using
  `evaluate=False`, but this is fragile and tends to be bug prone because
  other code may be written expecting the automatic evaluation to occur. As in
  the previous point, if we evaluate every value, there is little point to
  even having a symbolic function in the first place. For example, we might be
  tempted to evaluate some trig identities on `versin` in `eval`, but then
  these identities would always evaluate, and it wouldn't be possible to
  represent one half of the identity.

  One should also avoid doing anything in `eval` that is slow to compute.
  SymPy generally assumes that it is cheap to create expressions, and if this
  is not true, it can lead to performance issues.

  Finally, it is recommended to avoid performing automatic evaluation in
  `eval` based on assumptions. Instead, `eval` should typically only evaluate
  explicit numerical special values and return `None` for everything else. You
  might have noticed in [the example
  above](custom-functions-versin-eval-example) that we used `isinstance(n,
  Integer)` instead of checking `n.is_integer` using the assumptions system.
  We could have done that instead, which would make `versin(n*pi)` evaluate
  even if `n = Symbol('n', integer=True)`. But this is a case where we might
  not always want evaluation to happen, and if `n` is a more complicated
  expression, `n.is_integer` might take more time to compute.

- **When restricting the input domain: allow `None` input assumptions.**

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
  `m.is_integer` to be `True`. If it is `None`, it will fail (see the [guide
  on booleans and three-valued logic](booleans-guide) for details on what it
  means for an assumption to be `None`). This is problematic for two reasons.
  Firstly, it forces the user to define assumptions on any input variable. If
  the user omits them, it will fail:

  ```py
  >>> n, m = symbols('n m')
  >>> divides(m, n)
  Traceback (most recent call last):
  ...
  TypeError: m and n should be integers
  ```

  Instead they have to write

  ```py
  >>> n, m = symbols('n m', integer=True)
  >>> divides(m, n)
  divides(m, n)
  ```

  This may seem like an acceptable restriction, but there is a bigger problem.
  Sometimes, SymPy's assumptions system cannot deduce an assumption, even
  though it is mathematically true. In this case, it will give `None` (`None`
  means "undefined, "unknown", and "cannot compute" in SymPy's assumptions).
  For example

  ```py
  >>> # n and m are still defined as integer=True as above
  >>> divides(2, (m**2 + m)/2)
  Traceback (most recent call last):
  ...
  TypeError: m and n should be integers
  ```

  Here the expression `(m**2 + m)/2` is always an integer, but SymPy's
  assumptions system is not able to deduce this

  ```py
  >>> print(((m**2 + m)/2).is_integer)
  None
  ```

  SymPy's assumptions system is always improving, but there will always be
  cases like this that it cannot deduce, due to the fundamental computation
  complexity of the problem.

  Consequently, one should always test *negated* assumptions for input
  variables, that is, fail if the assumption is `False` but allow the
  assumption to be `None`.


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

  ```py
  >>> divides(1.5, 1)
  Traceback (most recent call last):
  ...
  TypeError: m and n should be integers
  ```

  But it does not fail in cases where the assumption is `None`:

  ```py
  >>> divides(2, (m**2 + m)/2)
  divides(2, m**2/2 + m/2)
  >>> _.subs(m, 2)
  1
  >>> n, m = symbols('n m') # Redefine n and m without the integer assumption
  >>> divides(m, n)
  divides(m, n)
  ```

(custom-functions-assumptions)=
### Assumptions

The next thing we might want to define are the assumptions on our function.
The [guide on the assumptions system](assumptions-guide) goes into the
assumptions system in great detail. It is recommended to read through that
guide first to understand what the different assumptions mean and how the
assumptions system works.

The simplest case is a function that always has a given assumption regardless
of its input. For example, the `divides` example we [showed
above](custom-functions-divides-example) is always an integer, because its
value is always either 0 or 1. In this case, you can define <code
class="docutils literal notranslate"><span
class="pre">is_*assumption*</span></code> directly on the class.

```py
>>> class divides(Function):
...     is_integer = True
```

```py
>>> divides(m, n).is_integer
True
```

```{note}

From here on out in this guide, in the interest of space, we will omit the
previous method definitions in the examples unless they are needed for the
given example to work. There are [complete
examples](custom-functions-complete-examples) with all methods included at the
end of this guide.

```

In general, however, the assumptions of a function depend on the assumptions
of its inputs. In this case, you should define an <code class="docutils literal notranslate"><span
class="pre">\_eval\_*assumption*</span> method.

For our $\operatorname{versin}$ function, the function is always between 0 and
2 for real valued inputs, and it is 0 whenever the input is an even multiple
of π. So `versin(x)` should be *nonnegative* whenever `x` is *real* (remember
that by default, a function's domain is all of $\mathbb{C}$, and indeed our
`versin` will make sense with non-real `x`) and *positive* whenever `x` is
*real* and not an *even* multiple of π.

In the assumptions handler methods, as in all methods, we can access the
arguments of the function using `.args`.

```py
>>> from sympy.core.logic import fuzzy_and, fuzzy_not
>>> class versin(Function):
...     def _eval_is_nonnegative(self):
...         x = self.args[0]
...         return x.is_real
...
...     def _eval_is_positive(self):
...         x = self.args[0]
...         # versin(x) is positive if x is real and not an even multiple of pi
...         return fuzzy_and([x.is_real, fuzzy_not((x/pi).is_even)])
```

```py
>>> versin(1).is_nonnegative
True
```

Note the use of `fuzzy_` functions in the more complicated `_eval_is_positive`
handler. It is important when working with assumptions to always be careful
about [handling three-valued logic correctly](booleans-guide).

It is not necessary to define `_eval_is_real` because it is deduced automatically
from the other assumptions, since `nonnegative -> real`. In general, you
should avoid defining assumptions that the assumptions system can deduce
automatically given its [known facts](assumptions-guide-predicates).

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

Finally, a word of warning: be very careful about correctness when coding
assumptions. Make sure to use the exact
[definitions](assumptions-guide-predicates) of the various assumptions, and
always check that you're handling `None` cases correctly with the fuzzy
three-valued logic functions. Incorrect or inconsistent assumptions can lead
to subtle bugs.

(custom-functions-evalf)=
### Numerical Evaluation with `evalf`

<!-- TODO: this goes over the basics, but it might be useful to have a
     separate guide dedicated to numerical evaluation. -->

Here we show how to define how to define how a function should  numerically
evaluate to a floating point {class}`~.Float` value, for instance, via `evalf`.

If your function has the same name as a function in
[mpmath](https://mpmath.org/doc/), as is the case for most functions in SymPy,
numerical evaluation will happen automatically and you do not need to do
anything.

If this is not the case, numerical evaluation can be specified by defining the
method `_eval_evalf(self, prec)`, where `prec` binary precision of the input.
The method should return the expression evaluated to the given precision, or
`None` if this is not possible.

```{note}

The `prec` argument to `_eval_evalf` is the *binary* precision, that is, the
number of bits in the floating-point representation. This differs from the
first argument to the `evalf` method, which is the *decimal* precision, or
`dps`. For example, the default binary precision of `Float` is 53,
corresponding to a decimal precision of 15. Therefore, when your `_eval_evalf`
method recursively calls evalf on another expression, it should call
`expr._eval_evalf(prec)` rather than `expr.evalf(prec)`, as the latter will
incorrectly use `prec` as the decimal precision.

```

We can define numerical evaluation for our example $\operatorname{versin}(x)$
function by recursively evaluating $1 - \cos(x)$.

```py
>>> class versin(Function):
...     def _eval_evalf(self, prec):
...         return (1 - cos(self.args[0]))._eval_evalf(prec)
```

```
>>> versin(1).evalf()
0.459697694131860
```

Once `_eval_evalf` is defined, this enables the automatic evaluation of
floating-point inputs. It is not required to implement this manually in
[`eval`](custom-functions-eval).

```py
>>> versin(1.)
0.459697694131860
```

<!-- TODO: Mention _should_evalf() here? Seems like something most people
     shouldn't mess with. -->

Note that `evalf` may be passed any expression, not just one that can be
evaluated numerically. In this case, it is expected that the numerical parts
of an expression will be evaluated. A general pattern to follow is to
recursively call `_eval_evalf` on the arguments of the function.

Whenever possible, it's best to reuse the evalf functionality defined in
existing SymPy functions. However, in some cases it will be necessary to call
out to `mpmath` directly.

(custom-functions-rewriting-and-simplification)=
### Rewriting and Simplification

#### Rewrite

#### Simplify

#### doit

(custom-functions-differentiation)=
### Differentiation

To define differentiation via {func}`~.diff`, define a method `fdiff` on the
class with an argument `argindex`. `fdiff` should return the derivative of the
function without considering the chain rule, with respect to the `argindex`-th
variable. `argindex` is indexed starting at `1`.

That is, `f(x1, ..., xi, ..., xn).fdiff(i)` should return $\frac{d}{d x_i}
f(x_1, \ldots, x_i, \ldots, x_n)$, where $x_k$ are independent of one another.
`diff` will automatically apply the chain rule using the result of `fdiff`.
User code should use `diff` and not call `fdiff` directly.

```{note}

`Function` subclasses should define differentiation using `fdiff`. Subclasses
of {class}`~.Expr` that aren't `Function` subclasess will need to define
`_eval_derivative` instead. It is not recommended to redefine
`_eval_derivative` on a `Function` subclass.

```

For our example function $\operatorname{versin}(x) = 1 - \cos(x)$, the
derivative is $\sin(x)$.

```py
>>> from sympy import sin
>>> class versin(Function):
...     def fdiff(self, argindex=1):
...         return sin(self.args[0])
```

```py
>>> versin(x).diff(x)
sin(x)
>>> versin(x**2).diff(x)
2*x*sin(x**2)
```

As an example of a function that has multiple arguments, let's define a
function for [fused
multiply-add](https://en.wikipedia.org/w/index.php?title=Fused_multiply_add):
$\operatorname{FMA}(x, y, z) = xy + z$.

We have

$$\frac{d}{dx} \operatorname{FMA}(x, y, z) = y,$$
$$\frac{d}{dy} \operatorname{FMA}(x, y, z) = x,$$
$$\frac{d}{dz} \operatorname{FMA}(x, y, z) = 1.$$

So the `fdiff` method for `FMA` would look like this

```py
>>> from sympy import Number, symbols
>>> x, y, z = symbols('x y z')
>>> class FMA(Function):
...     """
...     FMA(x, y, z) = x*y + z
...     """
...     def fdiff(self, argindex):
...         x, y, z = self.args
...         if argindex == 1:
...             return y
...         elif argindex == 2:
...             return x
...         elif argindex == 3:
...             return 1
```

```py
>>> FMA(x, y, z).diff(x)
y
>>> FMA(x, y, z).diff(y)
x
>>> FMA(x, y, z).diff(z)
1
>>> FMA(x**2, x + 1, y).diff(x)
x**2 + 2*x*(x + 1)
```

A [more complete example for `FMA`](custom-functions-example-fma) is given
below.

To leave a derivative unevaluated, raise
`sympy.core.function.ArgumentIndexError(self, argindex)`. This is the default
behavior if `fdiff` is not defined. Here is an example function $f(x, y)$ that
is linear in the first argument and has an unevaluated derivative on the
second argument.

```py
>>> from sympy.core.function import ArgumentIndexError
>>> class f(Function):
...    @classmethod
...    def eval(cls, x, y):
...        pass
...
...    def fdiff(self, argindex):
...        if argindex == 1:
...           return 1
...        raise ArgumentIndexError(self, argindex)
```

```py
>>> f(x, y).diff(x)
1
>>> f(x, y).diff(y)
Derivative(f(x, y), y)
```

### Series Expansions

### Printing

You can define how a function prints itself with the varions
[printers](module-printing) such as the {class}`string printer
<sympy.printing.str.StrPrinter>`, {func}`pretty printers
<sympy.printing.pretty.pretty.PrettyPrinter>`, and {func}`LaTeX printer
<sympy.printing.latex.LatexPrinter>`, as well as code printers for various
languages such as {class}`C <sympy.printing.c.C99CodePrinter>` and
{class}`Fortran <sympy.printing.fortran.FCodePrinter>`.


In most cases, you will not need to define any printing methods. The default
behavior is to print functions using their name. However, in some cases we may
want to define special printing for a function.

For example, for our [divides example
above](custom-functions-divides-example), we may want the LaTeX printer to
print a more mathematical expression. Let's make it print `divides(m, n)` as
`\left [ m \middle | n \right ]`, which looks like $\left [ m \middle | n
\right ]$ (here $[P]$ is the [Iverson
bracket](https://en.wikipedia.org/wiki/Iverson_bracket), which is $1$ if $P$
is true and $0$ if $P$ is false).

There are two primary ways to define printing for SymPy objects. One is to
define a printer on the printer class. Most classes that are part of the SymPy
library should use this method, by defining the printers on the respective
classes in sympy.printing. For user code, this may be preferable if you are
defining a custom printer, or if you have many custom functions that you want
to define printing for. See [](printer_example) for an example of thos to
define a printer in this way.

The other method is to define the printing on a method on the class. To do
this, first look up the `printmethod` attribute on the printer you want to
define the printing for. This is the name of the method you should define for
that printer. For the LaTeX printer, {attr}`.LatexPrinter.printmethod` is
`'_latex'`. So you should define a function `_latex(self, printer)` on the
class. For our `divides` example, it might look like

```py
>>> from sympy import latex
>>> class divides(Function):
...     def _latex(self, printer):
...         m, n = self.args
...         _m, _n = printer._print(m), printer._print(n)
...         return r'\left [ %s \middle | %s \right ]' % (_m, _n)
```

```py
>>> print(latex(divides(m, n)))
\left [ m \middle | n \right ]
```

See [](printer_method_example) for more details on how to define printer
methods and some pitfalls to avoid. Most importantly, you should always use
`printer._print()` to recursively print the arguments of the function inside
of a custom printer.

### Other Methods

#### `inverse`

#### `as_real_imag`

#### Miscellaneous `_eval_*` methods


(custom-functions-complete-examples)=
## Complete Examples

(custom-functions-example-fma)=
### FMA

```py
>>> from sympy import Number, symbols
>>> x, y, z = symbols('x y z')
>>> class FMA(Function):
...     """
...     FMA(x, y, z) = x*y + z
...     """
...     @classmethod
...     def eval(cls, x, y, z):
...         # Number is the base class of Integer, Rational, and Float
...         if all(isinstance(i, Number) for i in [x, y, z]):
...            return x*y + z
...
...     def fdiff(self, argindex):
...         x, y, z = self.args
...         if argindex == 1:
...             return y
...         elif argindex == 2:
...             return x
...         elif argindex == 3:
...             return 1
```

## Additional Tips

- SymPy includes dozens of functions. These can serve as useful examples for
  how to write a custom function, especially if the function is similar to one
  that is already implemented. Remember that everything in this guide applies
  equally well to functions that are included with SymPy and user defined
  functions.

- If you have many custom functions that share common logic, you can use a
  common base class to contain this shared logic. For an example of this, see
  the source code for the trigonometric functions in SymPy, which use
  `TrigonometricFunction`, `InverseTrigonometricFunction`, and
  `ReciprocalTrigonometricFunction` base classes with some shared logic.

- As with any code, it is a good idea to write extensive tests for your
  function. The SymPy test suite is a good resource for examples of how to
  write tests for such functions.
