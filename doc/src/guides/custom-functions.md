(custom-functions)=
# Writing Custom Functions

<!-- Note to contributors: if you update one of the examples in this guide, be
     sure to also update the complete example at the end. -->

This guide will describe how to create custom function classes in SymPy.
Custom user defined functions use the same mechanisms as the {ref}`functions
<functions>` that are included with SymPy such as the common {ref}`elementary
functions <elementary-functions>` like {class}`~.exp()` or {class}`~.sin()`,
{ref}`special functions <special-functions>` like {class}`~.gamma()` or
{class}`~.Si()`, and {ref}`combinatorial functions <combinatorial-functions>`
and {ref}`number theory functions <ntheory-module>` like
{class}`~.factorial()` or {func}`~.primepi()`. Consequently, this guide serves
both as a guide to end users who want to define their own custom functions and
to SymPy developers wishing to extend the functions included with SymPy.

This guide describes how to define complex valued functions, that is functions
that map a subset of $\mathbb{C}^n$ to $\mathbb{C}$. Functions that accept or
return other kinds of objects than complex numbers should subclass another
class, such as {class}`~.Boolean`, {class}`~.MatrixExpr`, {class}`~.Expr`, or
{class}`~.Basic`. Some of what is written here will apply to general
{class}`~.Basic` or {class}`~.Expr` subclasses, but much of it only applies to
{class}`~.Function` subclasses.

## Easy Cases: Fully Symbolic or Fully Evaluated

Before digging into the more advanced functionality for custom functions, we
should mention two common cases, the case where the function is fully
symbolic, and the case where the function is fully evaluated. Both of these
cases have much simpler alternatives than the full mechanisms described in this
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

In this case, you should use a normal Python function using the `def`
keyword:

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

If you find yourself defining an [`eval()`](custom-functions-eval) method on a
`Function` subclass where you always return a value and never return `None`,
you should consider just using a normal Python function instead, as there is
no benefit to using a symbolic `Function` subclass in that case (see the
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
accurately represent symbolic values. For example, in the above Python `def`
definition of `f`, `f(x)` implicitly assumes that `x` is nonzero. The
`Piecewise` version handles this case correctly and won't evaluate to the $x
\neq 0$ case unless `x` is known to not be zero.

Another option, if you want a function that not only evaluates, but always
evaluates to a numerical value, is to use {func}`~.lambdify`. This will
convert a SymPy expression into a function that can be evaluated using NumPy.

```
>>> from sympy import lambdify
>>> func = lambdify(x, Piecewise((0, Eq(x, 0)), (x + 1, True)))
>>> import numpy as np # doctest: +SKIP
>>> func(np.arange(5)) # doctest: +SKIP
array([0., 2., 3., 4., 5.])
```

Ultimately, the correct tool for the job depends on what you are doing and
what exact behavior you want.

(custom-functions-function-subclass)=
## Creating a Custom Function

The first step to creating a custom function is to subclass
{class}`~.Function`. The name of the subclass will be the name of the
function. Different methods should then be defined on this subclass, depending
on what functionality you want to provide.

As a motivating example for this document, let's create a custom function
class representing the [versine
function](https://en.wikipedia.org/wiki/Versine). Versine is a trigonometric
function which was used historically alongside some of the more familiar
trigonometric functions like sine and cosine. It is rarely used today. Versine
can be defined by the identity

(custom-functions-versine-definition)=

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
Note that `versin` is a class, and `versin(x)` is an instance of this class.

```py
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
### Defining Automatic Evaluation with `eval()`

```{sidebar} Reminder
Remember that `eval()` should be defined with the `@classmethod` decorator.
```

The first and most common thing we might want to define on our custom function
is automatic evaluation, that is, the cases where it will return an actual
value instead of just remaining unevaluated as-is.

This is done by defining the class method `eval()`. `eval()` should take the
arguments of the function and return either a value or `None`. If it returns
`None`, the function will remain unevaluated in that case. This also serves to
define the signature of the function (by default, without an `eval()` method, a
`Function` subclass will accept any number of arguments).

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
```

Here we make use of the fact that if a Python function does not explicitly
return a value, it automatically returns `None`. So in the cases where the `if
isinstance(n, Integer)` statement is not triggered, `eval()` returns `None`
and `versin` remains unevaluated.

```py
>>> versin(x*pi)
versin(pi*x)
```

```{note}

`Function` subclasses should not redefine `__new__` or `__init__`. If you want
to implement behavior that isn't possible with `eval()`, it might make more
sense to subclass {class}`~.Expr` rather than `Function`.

```

`eval()` can take any number of arguments, including an arbitrary number with
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
automatically once [`evalf()` is defined](custom-functions-evalf), so you do
not need to handle it explicitly in `eval()`.

(custom-functions-eval-best-practices)=
#### Best Practices for `eval()`

Certain antipatterns are common when defining `eval()` methods and should be
avoided.

- **Don't just return an expression.**

  In the above example, we might have been tempted to write

  ```py
  >>> from sympy import cos
  >>> class versin(Function):
  ...     @classmethod
  ...     def eval(cls, x):
  ...         # !! Not actually a good eval() method !!
  ...         return 1 - cos(x)
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
  actually just be `(1 - cos(x)).diff(x)`, instead of calling [the `fdiff()`
  method we define below](custom-functions-differentiation).

  ```{admonition} Key Point

  **The purpose of `eval()` is not to define what the function *is*,
  mathematically, but rather to specify on what inputs it should automatically
  evaluate.** The mathematical definition of a function is determined through
  the specification of various mathematical properties with the methods
  outlined below, like [numerical evaluation](custom-functions-evalf),
  [differentiation](custom-functions-differentiation), and so on.

  ```

  If you find yourself doing this, you should think about what you actually
  want to achieve. If you just want a shorthand function for an expression, it
  will be simpler to just [define a Python
  function](custom-functions-fully-evaluated). If you really do want a
  symbolic function, think about when you want it to evaluate to something
  else and when you want it to stay unevaluated. One option is to make your
  function unevaluated in `eval()` and define a [`doit()`
  method](custom-functions-doit) to evaluate it.

(custom-functions-automatic-evaluation)=
- **Avoid too much automatic evaluation.**

  It is recommended to minimize what is evaluated automatically by `eval()`.
  It is typically better to put more advanced simplifications in [other
  methods](custom-functions-rewriting-and-simplification), like
  [`doit()`](custom-functions-doit). Remember that whatever you define for
  automatic evaluation will *always* evaluate.[^evaluate-footnote] As in the
  previous point, if you evaluate every value, there is little point to even
  having a symbolic function in the first place. For example, we might be
  tempted to evaluate some trig identities on `versin` in `eval()`, but then
  these identities would always evaluate, and it wouldn't be possible to
  represent one half of the identity.

  [^evaluate-footnote]: While it is technically possible to bypass automatic
      evaluation by using `evaluate=False`, this is recommended against for
      two reasons. Firstly, `evaluate=False` is fragile because any function
      that rebuilds the expression from its `.args` will not keep the
      `evaluate=False` flag, causing it to evaluate. Secondly,
      `evaluate=False` tends to be bug prone, because other code may be
      written expecting the invariants from the automatic evaluation to hold.
      It is much better to not evaluate such cases at all in `eval()`, and
      move such simplifications to [`doit()`](custom-functions-doit) instead.

  One should also avoid doing anything in `eval()` that is slow to compute.
  SymPy generally assumes that it is cheap to create expressions, and if this
  is not true, it can lead to performance issues.

  Finally, it is recommended to avoid performing automatic evaluation in
  `eval()` based on assumptions. Instead, `eval()` should typically only
  evaluate explicit numerical special values and return `None` for everything
  else. You might have noticed in [the example
  above](custom-functions-versin-eval-example) that we used `isinstance(n,
  Integer)` instead of checking `n.is_integer` using the assumptions system.
  We could have done that instead, which would make `versin(n*pi)` evaluate
  even if `n = Symbol('n', integer=True)`. But this is a case where we might
  not always want evaluation to happen, and if `n` is a more complicated
  expression, `n.is_integer` might be more expensive to compute.

  Let's consider an example. Using the identity $\cos(x + y) = \cos(x)\cos(y) - \sin(x)\sin(y)$, we can derive the identity

  $$\operatorname{versin}(x + y) =
  \operatorname{versin}(x)\operatorname{versin}(y) -
  \operatorname{versin}(x) - \operatorname{versin}(y) - \sin(x)\sin(y) + 1.$$

  Suppose we decided to automatically expand this in `eval()`:

  ```
  >>> from sympy import Add, sin
  >>> class versin(Function):
  ...     @classmethod
  ...     def eval(cls, x):
  ...         # !! Not actually a good eval() method !!
  ...         if isinstance(x, Add):
  ...             a, b = x.as_two_terms()
  ...             return (versin(a)*versin(b) - versin(a) - versin(b)
  ...                     - sin(a)*sin(b) + 1)
  ```

  This method recursively splits `Add` terms into two parts and applies the
  above identity.

  ```
  >>> x, y, z = symbols('x y z')
  >>> versin(x + y)
  -sin(x)*sin(y) + versin(x)*versin(y) - versin(x) - versin(y) + 1
  ```

  But now it's impossible to represent `versin(x + y)` without it expanding.
  This will affect other methods too. For example, suppose we define
  [differentiation (see below)](custom-functions-differentiation):

  ```
  >>> class versin(Function):
  ...     @classmethod
  ...     def eval(cls, x):
  ...         # !! Not actually a good eval() method !!
  ...         if isinstance(x, Add):
  ...             a, b = x.as_two_terms()
  ...             return (versin(a)*versin(b) - versin(a) - versin(b)
  ...                     - sin(a)*sin(b) + 1)
  ...
  ...     def fdiff(self, argindex=1):
  ...         return sin(self.args[0])
  ```

  We would expect `versin(x + y).diff(x)` to return `sin(x + y)`, and indeed,
  if we hadn't expanded this identity in `eval()`, [it
  would](custom-functions-differentiation-examples). But with this version,
  `versin(x + y)` gets automatically expanded before `diff()` gets called,
  instead we get a more complicated expression:

  ```
  >>> versin(x + y).diff(x)
  sin(x)*versin(y) - sin(x) - sin(y)*cos(x)
  ```

  And things are even worse than that. Let's try an `Add` with three terms:

  ```
  >>> versin(x + y + z)
  (-sin(y)*sin(z) + versin(y)*versin(z) - versin(y) - versin(z) +
  1)*versin(x) - sin(x)*sin(y + z) + sin(y)*sin(z) - versin(x) -
  versin(y)*versin(z) + versin(y) + versin(z)
  ```

  We can see that things are getting out of control quite quickly. In fact,
  `versin(Add(*symbols('x:100')))` (`versin()` on an `Add` with 100 terms)
  takes over a second to evaluate, and that's just to *create* the expression,
  without even doing anything with it yet.

  Identities like this are better left out of `eval` and implemented in other
  methods instead (in the case of this identity,
  [`expand_trig()`](custom-functions-expand)).

- **When restricting the input domain: allow `None` input assumptions.**

  Our example function $\operatorname{versin}(x)$ is a function from
  $\mathbb{C}$ to $\mathbb{C}$, so it can accept any input. But suppose we had
  a function that only made sense with certain inputs. As a second example,
  let's define a function `divides` as

  (custom-functions-divides-definition)=

  $$\operatorname{divides}(m, n) = \begin{cases} 1 & \text{for}\: m \mid n \\
  0 & \text{for}\: m\not\mid n  \end{cases}.$$

  That is, `divides(m, n)` will be `1` if `m` divides `n` and `0` otherwise.
  `divides` clearly only makes sense if `m` and `n` are integers.

  We might be tempted to define the `eval()` method for `divides` like this:

  ```py
  >>> class divides(Function):
  ...     @classmethod
  ...     def eval(cls, m, n):
  ...         # !! Not actually a good eval() method !!
  ...
  ...         # Evaluate for explicit integer m and n. This part is fine.
  ...         if isinstance(m, Integer) and isinstance(n, Integer):
  ...             return int(n % m == 0)
  ...
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
  >>> print(n.is_integer)
  None
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
  means both "undefined" and "cannot compute" in SymPy's assumptions). For
  example

  ```py
  >>> # n and m are still defined as integer=True as above
  >>> divides(2, (m**2 + m)/2)
  Traceback (most recent call last):
  ...
  TypeError: m and n should be integers
  ```

  Here the expression `(m**2 + m)/2` is always an integer, but SymPy's
  assumptions system is not able to deduce this:

  ```py
  >>> print(((m**2 + m)/2).is_integer)
  None
  ```

  SymPy's assumptions system is always improving, but there will always be
  cases like this that it cannot deduce, due to the fundamental computational
  complexity of the problem, and the fact that the general problem is
  [often](https://en.wikipedia.org/wiki/Hilbert%27s_tenth_problem)
  [undecidable](https://en.wikipedia.org/wiki/Richardson%27s_theorem).

  Consequently, one should always test *negated* assumptions for input
  variables, that is, fail if the assumption is `False` but allow the
  assumption to be `None`.

  (custom-functions-divides-eval)=
  ```py
  >>> class divides(Function):
  ...     @classmethod
  ...     def eval(cls, m, n):
  ...         # Evaluate for explicit integer m and n. This part is fine.
  ...         if isinstance(m, Integer) and isinstance(n, Integer):
  ...             return int(n % m == 0)
  ...
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
  0
  >>> n, m = symbols('n m') # Redefine n and m without the integer assumption
  >>> divides(m, n)
  divides(m, n)
  ```

  ```{note}

  This rule of allowing `None` assumptions only applies to instances
  where an exception would be raised, such as type checking an input domain.
  In cases where simplifications or other operations are done, one should
  treat a `None` assumption as meaning "can be either `True` or `False`" and
  not perform an operation that might not be mathematically valid.

  ```

(custom-functions-assumptions)=
### Assumptions

The next thing you might want to define are the assumptions on our function.
The assumptions system allows defining what mathematical properties your
function has given its inputs, for example, "$f(x)$ is *positive* when $x$ is
*real*."


The [guide on the assumptions system](assumptions-guide) goes into the
assumptions system in great detail. It is recommended to read through that
guide first to understand what the different assumptions mean and how the
assumptions system works.

The simplest case is a function that always has a given assumption regardless
of its input. In this case, you can define <code class="docutils
literal notranslate"><span class="pre">is_*assumption*</span></code> directly
on the class.

For example, our [example `divides`
function](custom-functions-divides-definition) is always an integer, because
its value is always either 0 or 1:


```{sidebar} Note

From here on out in this guide, in the interest of space, we will omit the
previous method definitions in the examples unless they are needed for the
given example to work. There are [complete
examples](custom-functions-complete-examples) at the end of this guide with
all the methods.

```

```py
>>> class divides(Function):
...     is_integer = True
...     is_negative = False
```

```py
>>> divides(m, n).is_integer
True
>>> divides(m, n).is_nonnegative
True
```

In general, however, the assumptions of a function depend on the assumptions
of its inputs. In this case, you should define an <code class="docutils literal notranslate"><span
class="pre">\_eval\_*assumption*</span></code> method.

For our [$\operatorname{versin}(x)$
example](custom-functions-versine-definition), the function is always in $[0,
2]$ when $x$ is real, and it is 0 exactly when $x$ is an even multiple of
$\pi$. So `versin(x)` should be *nonnegative* whenever `x` is *real* and
*positive* whenever `x` is *real* and not an *even* multiple of π. Remember
that by default, a function's domain is all of $\mathbb{C}$, and indeed
`versin(x)` makes perfect sense with non-real `x`.

To see if `x` is an even multiple of `pi`, we can use {meth}`~.as_independent`
to match `x` structurally as `coeff*pi`. Pulling apart subexpressions
structurally like this in assumptions handlers is preferable to using
something like `(x/pi).is_even`, because that will create a new expression
`x/pi`. The creation of a new expression is much slower. Furthermore, whenever
an expression is created, the constructors that are called when creating the
expression will often themselves cause assumptions to be queried. If you are
not careful, this can lead to infinite recursion. So a good general rule for
assumptions handlers is, **never create a new expression in an assumptions
handler**. Always pull apart the args of the function using structural methods
like `as_independent`.


Note that $\operatorname{versin}(x)$ can be
nonnegative for nonreal $x$, for example:

```py
>>> from sympy import I
>>> 1 - cos(pi + I*pi)
1 + cosh(pi)
>>> (1 - cos(pi + I*pi)).evalf()
12.5919532755215
```

So for the `_eval_is_nonnegative` handler, we want to return `True` if
`x.is_real` is `True` but `None` if `x.is_real` is either `False` or `None`.
It is left as an exercise to the reader to handle the cases for nonreal `x`
that make `versin(x)` nonnegative, using similar logic from the
`_eval_is_positive` handler.

In the assumptions handler methods, as in all methods, we can access the
arguments of the function using `self.args`.

```py
>>> from sympy.core.logic import fuzzy_and, fuzzy_not
>>> class versin(Function):
...     def _eval_is_nonnegative(self):
...         # versin(x) is nonnegative if x is real
...         x = self.args[0]
...         if x.is_real is True:
...             return True
...
...     def _eval_is_positive(self):
...         # versin(x) is positive if x is real and not an even multiple of pi
...         x = self.args[0]
...
...         # x.as_independent(pi, as_Add=False) will split x as a Mul of the
...         # form coeff*pi
...         coeff, pi_ = x.as_independent(pi, as_Add=False)
...         # If pi_ = pi, x = coeff*pi. Otherwise x is not (structurally) of
...         # the form coeff*pi.
...         if pi_ == pi:
...             return fuzzy_and([x.is_real, fuzzy_not(coeff.is_even)])
...         elif x.is_real is False:
...             return False
...         # else: return None. We do not know for sure whether x is an even
...         # multiple of pi
```

```py
>>> versin(1).is_nonnegative
True
>>> versin(2*pi).is_positive
False
>>> versin(3*pi).is_positive
True
```

Note the use of `fuzzy_` functions in the more complicated
`_eval_is_positive()` handler, and the careful handling of the `if`/`elif`. It
is important when working with assumptions to always be careful about
[handling three-valued logic correctly](booleans-guide). This ensures that the
method returns the correct answer when `x.is_real` or `coeff.is_even` are
`None`.

```{warning}

Never define <code class="docutils literal notranslate"><span
class="pre">is_*assumption*</span></code> as a `@property` method. Doing so
will break the automatic deduction of other assumptions. <code class="docutils
literal notranslate"><span class="pre">is_*assumption*</span></code> should
only ever be defined as a class variable equal to `True` or `False`. If the
assumption depends on the `.args` of the function somehow, define the <code
class="docutils literal notranslate"><span
class="pre">\_eval\_*assumption*</span></code> method.

```

In this example, it is not necessary to define `_eval_is_real()` because it is
deduced automatically from the other assumptions, since `nonnegative -> real`.
In general, you should avoid defining assumptions that the assumptions system
can deduce automatically given its [known
facts](assumptions-guide-predicates).

```py
>>> versin(1).is_real
True
```

The assumptions system is often able to deduce more than you might think.
For example, from the above, it can deduce that `versin(2*n*pi)` is zero when
`n` is an integer.

```py
>>> n = symbols('n', integer=True)
>>> versin(2*n*pi).is_zero
True
```

It's always worth checking if the assumptions system can deduce something
automatically before manually coding it.

Finally, a word of warning: be very careful about correctness when coding
assumptions. Make sure to use the exact
[definitions](assumptions-guide-predicates) of the various assumptions, and
always check that you're handling `None` cases correctly with the fuzzy
three-valued logic functions. Incorrect or inconsistent assumptions can lead
to subtle bugs. It's recommended to use unit tests to check all the various
cases whenever your function has a nontrivial assumption handler. All
functions defined in SymPy itself are required to be extensively tested.

(custom-functions-evalf)=
### Numerical Evaluation with `evalf()`

<!-- TODO: this goes over the basics, but it might be useful to have a
     separate guide dedicated to numerical evaluation. -->

Here we show how to define how a function should numerically evaluate to a
floating point {class}`~.Float` value, for instance, via `evalf()`.
Implementing numerical evaluation enables several behaviors in SymPy. For
example, once `evalf()` is defined, you can plot your function, and things
like inequalities can evaluate to explicit values.

If your function has the same name as a function in
[mpmath](https://mpmath.org/doc/current/), which is the case for most
functions included with SymPy, numerical evaluation will happen automatically
and you do not need to do anything.

If this is not the case, numerical evaluation can be specified by defining the
method `_eval_evalf(self, prec)`, where `prec` is the binary precision of the
input. The method should return the expression evaluated to the given
precision, or `None` if this is not possible.

```{note}

The `prec` argument to `_eval_evalf()` is the *binary* precision, that is, the
number of bits in the floating-point representation. This differs from the
first argument to the `evalf()` method, which is the *decimal* precision, or
`dps`. For example, the default binary precision of `Float` is 53,
corresponding to a decimal precision of 15. Therefore, if your `_eval_evalf()`
method recursively calls evalf on another expression, it should call
`expr._eval_evalf(prec)` rather than `expr.evalf(prec)`, as the latter will
incorrectly use `prec` as the decimal precision.

```

We can define numerical evaluation for [our example $\operatorname{versin}(x)$
function](custom-functions-versine-definition) by recursively evaluating
$2\sin^2\left(\frac{x}{2}\right)$, which is a more numerically stable way of writing $1 -
\cos(x)$.

```py
>>> from sympy import sin
>>> class versin(Function):
...     def _eval_evalf(self, prec):
...         return (2*sin(self.args[0]/2)**2)._eval_evalf(prec)
```

```
>>> versin(1).evalf()
0.459697694131860
```

Once `_eval_evalf()` is defined, this enables the automatic evaluation of
floating-point inputs. It is not required to implement this manually in
[`eval()`](custom-functions-eval).

```py
>>> versin(1.)
0.459697694131860
```

<!-- TODO: Mention _should_evalf() here? Seems like something most people
     shouldn't mess with. -->

Note that `evalf()` may be passed any expression, not just one that can be
evaluated numerically. In this case, it is expected that the numerical parts
of an expression will be evaluated. A general pattern to follow is to
recursively call `_eval_evalf(prec)` on the arguments of the function.

Whenever possible, it's best to reuse the evalf functionality defined in
existing SymPy functions. However, in some cases it will be necessary to use
mpmath directly.

(custom-functions-rewriting-and-simplification)=
### Rewriting and Simplification

Various simplification functions and methods allow specifying their behavior
on custom subclasses. Not every function in SymPy has such hooks. See the
documentation of each individual function for details.

(custom-functions-rewrite)=
#### `rewrite()`

The {meth}`~.rewrite` method allows rewriting an expression in terms of a
specific function or rule. For example,

```
>>> sin(x).rewrite(cos)
cos(x - pi/2)
```

To implement rewriting, define a method `_eval_rewrite(self, rule, args,
**hints)`, where

- `rule` is the *rule* passed to the `rewrite()` method. Typically `rule` will
  be the class of the object to be rewritten to, although for more complex
  rewrites, it can be anything. Each object that defines `_eval_rewrite()`
  defines what rule(s) it supports. Many SymPy functions rewrite to common
  classes, like `expr.rewrite(Add)`, to perform simplifications or other
  computations.

- `args` are the arguments of the function to be used for rewriting. This
  should be used instead of `self.args` because any recursive expressions in
  the args will be rewritten in `args` (assuming the caller used
  `rewrite(deep=True)`, which is the default).

- `**hints` are additional keyword arguments which may be used to specify the
  behavior of the rewrite. Unknown hints should be ignored as they may be
  passed to other `_eval_rewrite()` methods. If you recursively call rewrite,
  you should pass the `**hints` through.


The method should return a rewritten expression, using `args` as the
arguments to the function, or `None` if the expression should be unchanged.

For our [`versin` example](custom-functions-versine-definition), an obvious
rewrite we can implement is rewriting `versin(x)` as `1 - cos(x)`:

```py
>>> class versin(Function):
...     def _eval_rewrite(self, rule, args, **hints):
...         if rule == cos:
...             return 1 - cos(*args)
>>> versin(x).rewrite(cos)
1 - cos(x)
```

Once we've defined this, {func}`~.simplify` is now able to simplify some
expressions containing `versin`:

```
>>> from sympy import simplify
>>> simplify(versin(x) + cos(x))
1
```

<!--

TODO: simplify() can be customized with _eval_simplify(), but it isn't very
powerful right now (see https://github.com/sympy/sympy/issues/19281). Most
users are better off just defining rewrite(), doit(), and/or expand().

#### Simplify

-->

(custom-functions-doit)=
#### `doit()`

The {meth}`doit() <sympy.core.basic.Basic.doit>` method is used to evaluate
"unevaluated" functions. To define `doit()` implement `doit(self, deep=True,
**hints)`. If `deep=True`, `doit()` should recursively call `doit()` on the
arguments. `**hints` will be any other keyword arguments passed to the user,
which should be passed to any recursive calls to `doit()`. You can use `hints`
to allow the user to specify specific behavior for `doit()`.

The typical usage of `doit()` in custom `Function` subclasses is to perform more
advanced evaluation which is not performed in [`eval()`](custom-functions-eval).

For example, for our [`divides` example](custom-functions-divides-definition),
there are several instances that could be simplified using some identities.
For example, we defined `eval()` to evaluate on explicit integers, but we might
also want to evaluate examples like `divides(k, k*n)` where the divisibility
is symbolically true. One of the [best practices for
`eval()`](custom-functions-eval-best-practices) is to avoid too much automatic
evaluation. Automatically evaluating in this case might be considered too
much, as it would make use of the assumptions system, which could be
expensive. Furthermore, we might want to be able to represent `divides(k,
k*n)` without it always evaluating.

The solution is to implement these more advanced evaluations in `doit()`. That
way, we can explicitly perform them by calling `expr.doit()`, but they won't
happen by default. An example `doit()` for `divides` that performs this
simplification (along with the [above definition of
`eval()`](custom-functions-divides-eval)) might look like this:

```{note}
If  `doit()` returns a Python `int` literal, convert it to an `Integer` so
that the returned object is a SymPy type.
```

```
>>> from sympy import Integer
>>> class divides(Function):
...     # Define evaluation on basic inputs, as well as type checking that the
...     # inputs are not nonintegral.
...     @classmethod
...     def eval(cls, m, n):
...         # Evaluate for explicit integer m and n.
...         if isinstance(m, Integer) and isinstance(n, Integer):
...             return int(n % m == 0)
...
...         # For symbolic arguments, require m and n to be integer.
...         if m.is_integer is False or n.is_integer is False:
...             raise TypeError("m and n should be integers")
...
...     # Define doit() as further evaluation on symbolic arguments using
...     # assumptions.
...     def doit(self, deep=False, **hints):
...         m, n = self.args
...         # Recursively call doit() on the args whenever deep=True.
...         # Be sure to pass deep=True and **hints through here.
...         if deep:
...            m, n = m.doit(deep=deep, **hints), n.doit(deep=deep, **hints)
...
...         # divides(m, n) is 1 iff n/m is an integer. Note that m and n are
...         # already assumed to be integers because of the logic in eval().
...         isint = (n/m).is_integer
...         if isint is True:
...             return Integer(1)
...         elif isint is False:
...             return Integer(0)
...         else:
...             return divides(m, n)
```

(Note that this uses the
[convention](https://en.wikipedia.org/wiki/Divisor#Definition) that $k \mid 0$
for all $k$ so that we do not need to check if `m` or `n` are nonzero. If we
used a different convention we would need to check if `m.is_zero` and
`n.is_zero` before performing the simplification.)

```
>>> n, m, k = symbols('n m k', integer=True)
>>> divides(k, k*n)
divides(k, k*n)
>>> divides(k, k*n).doit()
1
```

Another common way to implement `doit()` is for it to always return another
expression. This effectively treats the function as an "unevaluated" form of
another expression.

(custom-functions-fma-definition)=

For example, let's define a function for [fused
multiply-add](https://en.wikipedia.org/w/index.php?title=Fused_multiply_add):
$\operatorname{FMA}(x, y, z) = xy + z$. It may be useful to express this
function as a distinct function, e.g., for the purposes of code generation,
but it would also be useful in some contexts to "evaluate" `FMA(x, y, z)` to
`x*y + z` so that it can properly simplify with other expressions.

```
>>> from sympy import Number
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
...     def doit(self, deep=True, **hints):
...         x, y, z = self.args
...         # Recursively call doit() on the args whenever deep=True.
...         # Be sure to pass deep=True and **hints through here.
...         if deep:
...             x = x.doit(deep=deep, **hints)
...             y = y.doit(deep=deep, **hints)
...             z = z.doit(deep=deep, **hints)
...         return x*y + z
```

```
>>> x, y, z = symbols('x y z')
>>> FMA(x, y, z)
FMA(x, y, z)
>>> FMA(x, y, z).doit()
x*y + z
```

Most custom functions will not want to define `doit()` in this way. However,
this can provide a happy medium between having a function that always
evaluates and a function that never evaluates, producing a function that
doesn't evaluate by default but can be evaluated on demand (see the
[discussion above](custom-functions-eval-best-practices)).

(custom-functions-expand)=
#### `expand()`

The {func}`~.expand()` function "expands" an expression in various ways. It is
actually a wrapper around several sub-expansion hints. Each function
corresponds to a hint to the `expand()` function/method. A specific expand
*hint* can be defined in a custom function by defining <code class="docutils
literal notranslate"><span class="pre">\_eval_expand_<em>hint</em>(self,
**hints)</span></code>. See the documentation of {func}`~.expand` for details
on which hints are defined and the documentation for each specific <code
class="docutils literal notranslate"><span
class="pre">expand_*hint*()</span></code> function (e.g.,
{func}`~.expand_trig`) for details on what each hint is designed to do.

The `**hints` keyword arguments are additional hints that may be passed to the
expand function to specify additional behavior (these are separate from the
predefined *hints* described in the previous paragraph). Unknown hints should
be ignored as they may apply to other functions' custom `expand()` methods. A
common hint to define is `force`, where `force=True` would force an expansion
that might not be mathematically valid for all the given input assumptions.
For example, `expand_log(log(x*y), force=True)` produces `log(x) + log(y)`
even though this identity is not true for all complex `x` and `y` (typically
`force=False` is the default).

Note that `expand()` automatically takes care of recursively expanding
expressions using its own `deep` flag, so `_eval_expand_*` methods should not
recursively call expand on the arguments of the function.

For our [`versin` example](custom-functions-versine-definition), we can define
rudimentary `trig` expansion by defining an `_eval_expand_trig` method,
which recursively calls `expand_trig()` on `1 - cos(x)`:

```
>>> from sympy import expand_trig
>>> y = symbols('y')
>>> class versin(Function):
...    def _eval_expand_trig(self, **hints):
...        x = self.args[0]
...        return expand_trig(1 - cos(x))
>>> versin(x + y).expand(trig=True)
sin(x)*sin(y) - cos(x)*cos(y) + 1
```

A more sophisticated implementation might attempt to rewrite the result of
`expand_trig(1 - cos(x))` back into `versin` functions. This is left as an
exercise for the reader.

(custom-functions-differentiation)=
### Differentiation

To define differentiation via {func}`~.diff`, define a method `fdiff(self,
argindex)`. `fdiff()` should return the derivative of the function, without
considering the chain rule, with respect to the `argindex`-th variable.
`argindex` is indexed starting at `1`.

That is, `f(x1, ..., xi, ..., xn).fdiff(i)` should return $\frac{d}{d x_i}
f(x_1, \ldots, x_i, \ldots, x_n)$, where $x_k$ are independent of one another.
`diff()` will automatically apply the chain rule using the result of
`fdiff()`. User code should use `diff()` and not call `fdiff()` directly.

```{note}

`Function` subclasses should define differentiation using `fdiff()`. Subclasses
of {class}`~.Expr` that aren't `Function` subclasses will need to define
`_eval_derivative()` instead. It is not recommended to redefine
`_eval_derivative()` on a `Function` subclass.

```

For our [$\operatorname{versin}$ example
function](custom-functions-versine-definition), the derivative is $\sin(x)$.

```py
>>> class versin(Function):
...     def fdiff(self, argindex=1):
...         # argindex indexes the args, starting at 1
...         return sin(self.args[0])
```

(custom-functions-differentiation-examples)=
```py
>>> versin(x).diff(x)
sin(x)
>>> versin(x**2).diff(x)
2*x*sin(x**2)
>>> versin(x + y).diff(x)
sin(x + y)
```

As an example of a function that has multiple arguments, consider the [fused
multiply-add (FMA) example](custom-functions-fma-definition) defined above
($\operatorname{FMA}(x, y, z) = xy + z$).

We have

$$\frac{d}{dx} \operatorname{FMA}(x, y, z) = y,$$
$$\frac{d}{dy} \operatorname{FMA}(x, y, z) = x,$$
$$\frac{d}{dz} \operatorname{FMA}(x, y, z) = 1.$$

So the `fdiff()` method for `FMA` would look like this:

```py
>>> from sympy import Number, symbols
>>> x, y, z = symbols('x y z')
>>> class FMA(Function):
...     """
...     FMA(x, y, z) = x*y + z
...     """
...     def fdiff(self, argindex):
...         # argindex indexes the args, starting at 1
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

To leave a derivative unevaluated, raise
`sympy.core.function.ArgumentIndexError(self, argindex)`. This is the default
behavior if `fdiff()` is not defined. Here is an example function $f(x, y)$ that
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

<!-- TODO: series expansion will need its own guide. See
     https://github.com/sympy/sympy/issues/23625. -->

<!-- TODO: integration will need its own guide. See
     https://github.com/sympy/sympy/issues/23624. -->

### Printing

You can define how a function prints itself with the varions
[printers](module-printing) such as the {class}`string printer
<sympy.printing.str.StrPrinter>`, {func}`pretty printers
<sympy.printing.pretty.pretty.PrettyPrinter>`, and {func}`LaTeX printer
<sympy.printing.latex.LatexPrinter>`, as well as code printers for various
languages like {class}`C <sympy.printing.c.C99CodePrinter>` and
{class}`Fortran <sympy.printing.fortran.FCodePrinter>`.

In most cases, you will not need to define any printing methods. The default
behavior is to print functions using their name. However, in some cases we may
want to define special printing for a function.

For example, for our [divides example
above](custom-functions-divides-definition), we may want the LaTeX printer to
print a more mathematical expression. Let's make the LaTeX printer represent
`divides(m, n)` as `\left [ m \middle | n \right ]`, which looks like $\left [
m \middle | n \right ]$ (here $[P]$ is the [Iverson
bracket](https://en.wikipedia.org/wiki/Iverson_bracket), which is $1$ if $P$
is true and $0$ if $P$ is false).

There are two primary ways to define printing for SymPy objects. One is to
define a printer on the printer class. Most classes that are part of the SymPy
library should use this method, by defining the printers on the respective
classes in `sympy.printing`. For user code, this may be preferable if you are
defining a custom printer, or if you have many custom functions that you want
to define printing for. See [](printer_example) for an example of how to
define a printer in this way.

The other way is to define the printing as a method on the function class. To
do this, first look up the `printmethod` attribute on the printer you want to
define the printing for. This is the name of the method you should define for
that printer. For the LaTeX printer, {attr}`.LatexPrinter.printmethod` is
`'_latex'`. The print method always takes one argument, `printer`.
`printer._print` should be used to recursively print any other expressions,
including the arguments of the function.

So to define our `divides` LaTeX printer, we will define the function
`_latex(self, printer)` on the class, like this:

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

Several other methods can be defined on custom functions to specify various behaviors.

#### `inverse()`

The `inverse(self, argindex=1)` method can be defined to specify the inverse
of the function. This is used by {func}`~.solve` and {func}`~.solveset`. The
`argindex` argument is the argument of the function, starting at 1 (similar to
the same argument name for the [`fdiff()`
method](custom-functions-differentiation)).

`inverse()` should return a function (not an expression) for the inverse. If the
inverse is a larger expression than a single function, it can return a
`lambda` function.

`inverse()` should only be defined for functions that are one-to-one. In other
words, `f(x).inverse()` is the [left
inverse](https://en.wikipedia.org/wiki/Inverse_function#Left_and_right_inverses)
of `f(x)`. Defining `inverse()` on a function that is not one-to-one may
result in `solve()` not giving all possible solutions to an expression
containing the function.

Our [example versine function](custom-functions-versine-definition) is not
one-to-one (because cosine is not), but its inverse $\operatorname{arcversin}$
is. We may define it as follows (using the same naming convention as other
inverse trig functions in SymPy):

```py
>>> class aversin(Function):
...     def inverse(self, argindex=1):
...         return versin
```

This makes `solve()` work on `aversin(x)`:

```
>>> from sympy import solve
>>> solve(aversin(x) - y, x)
[versin(y)]
```

#### `as_real_imag()`

The method {meth}`as_real_imag() <sympy.core.expr.Expr.as_real_imag>` method
defines how to split a function into its real and imaginary parts. It is used
by various SymPy functions that operate on the real and imaginary parts of an
expression separately.

`as_real_imag(self, deep=True, **hints)` should return a 2-tuple containing
the real part and imaginary part of the function. That is
`expr.as_real_imag()` returns `(re(expr), im(expr))`, where
`expr == re(expr) + im(expr)*I`, and `re(expr)` and `im(expr)` are real.

If `deep=True`, it should recursively call `as_real_imag(deep=True, **hints)`
on its arguments. As with [`doit()`](custom-functions-doit) and [the
`_eval_expand_*()` methods](custom-functions-expand), `**hints` may be any
hints to allow the user to specify the behavior of the method. Unknown hints
should be ignored and passed through on any recursive calls in case they are
meant for other `as_real_imag()` methods.

For our [`versin` example](custom-functions-versine-definition), we can
recursively use the `as_real_imag()` that is already defined for `1 - cos(x)`.

```py
>>> class versin(Function):
...     def as_real_imag(self, deep=True, **hints):
...         return (1 - cos(self.args[0])).as_real_imag(deep=deep, **hints)
>>> versin(x).as_real_imag()
(-cos(re(x))*cosh(im(x)) + 1, sin(re(x))*sinh(im(x)))
```

Defining `as_real_imag()` also automatically makes {func}`~.expand_complex`
work.

```py
>>> versin(x).expand(complex=True)
I*sin(re(x))*sinh(im(x)) - cos(re(x))*cosh(im(x)) + 1
```

#### Miscellaneous `_eval_*` methods

There are many other functions in SymPy whose behavior can be defined on
custom functions via a custom `_eval_*` method, analogous to the ones
described above. See the documentation of the specific function for details on
how to define each method.

(custom-functions-complete-examples)=
## Complete Examples

Here are complete examples for the example functions defined in this guide.
See the above sections for details on each method.

(custom-functions-versine-full-example)=
### Versine

The versine (versed sine) function is defined as

$$\operatorname{versin}(x) = 1 - \cos(x).$$

Versine is an example of a simple function defined for all complex numbers.
The mathematical definition is simple, which makes it straightforward to
define all the above methods on it (in most cases we can just reuse the
existing SymPy logic defined on `1 - cos(x)`).

#### Definition

```
>>> from sympy import Function, cos, expand_trig, Integer, pi, sin
>>> from sympy.core.logic import fuzzy_and, fuzzy_not
>>> class versin(Function):
...     r"""
...     The versine function.
...
...     $\operatorname{versin}(x) = 1 - \cos(x) = 2\sin(x/2)^2.$
...
...     Geometrically, given a standard right triangle with angle x in the
...     unit circle, the versine of x is the positive horizontal distance from
...     the right angle of the triangle to the rightmost point on the unit
...     circle. It was historically used as a more numerically accurate way to
...     compute 1 - cos(x), but it is rarely used today.
...
...     References
...     ==========
...
...     .. [1] https://en.wikipedia.org/wiki/Versine
...     .. [2] https://blogs.scientificamerican.com/roots-of-unity/10-secret-trig-functions-your-math-teachers-never-taught-you/
...     """
...     # Define evaluation on basic inputs.
...     @classmethod
...     def eval(cls, x):
...         # If x is an explicit integer multiple of pi, x/pi will cancel and
...         # be an Integer.
...         n = x/pi
...         if isinstance(n, Integer):
...             return 1 - (-1)**n
...
...     # Define numerical evaluation with evalf().
...     def _eval_evalf(self, prec):
...         return (2*sin(self.args[0]/2)**2)._eval_evalf(prec)
...
...     # Define basic assumptions.
...     def _eval_is_nonnegative(self):
...         # versin(x) is nonnegative if x is real
...         x = self.args[0]
...         if x.is_real is True:
...             return True
...
...     def _eval_is_positive(self):
...         # versin(x) is positive if x is real and not an even multiple of pi
...         x = self.args[0]
...
...         # x.as_independent(pi, as_Add=False) will split x as a Mul of the
...         # form n*pi
...         coeff, pi_ = x.as_independent(pi, as_Add=False)
...         # If pi_ = pi, x = coeff*pi. Otherwise pi_ = 1 and x is not
...         # (structurally) of the form n*pi.
...         if pi_ == pi:
...             return fuzzy_and([x.is_real, fuzzy_not(coeff.is_even)])
...         elif x.is_real is False:
...             return False
...         # else: return None. We do not know for sure whether x is an even
...         # multiple of pi
...
...     # Define the behavior for various simplification and rewriting
...     # functions.
...     def _eval_rewrite(self, rule, args, **hints):
...         if rule == cos:
...             return 1 - cos(*args)
...         elif rule == sin:
...             return 2*sin(x/2)**2
...
...     def _eval_expand_trig(self, **hints):
...         x = self.args[0]
...         return expand_trig(1 - cos(x))
...
...     def as_real_imag(self, deep=True, **hints):
...         # reuse _eval_rewrite(cos) defined above
...         return self.rewrite(cos).as_real_imag(deep=deep, **hints)
...
...     # Define differentiation.
...     def fdiff(self, argindex=1):
...         return sin(self.args[0])
```

#### Examples

**Evaluation:**

```
>>> x, y = symbols('x y')
>>> versin(x)
versin(x)
>>> versin(2*pi)
0
>>> versin(1.0)
0.459697694131860
```

**Assumptions:**

```
>>> n = symbols('n', integer=True)
>>> versin(n).is_real
True
>>> versin((2*n + 1)*pi).is_positive
True
>>> versin(2*n*pi).is_zero
True
>>> print(versin(n*pi).is_positive)
None
>>> r = symbols('r', real=True)
>>> print(versin(r).is_positive)
None
>>> nr = symbols('nr', real=False)
>>> print(versin(nr).is_nonnegative)
None
```

**Simplification:**

```
>>> a, b = symbols('a b', real=True)
>>> from sympy import I
>>> versin(x).rewrite(cos)
1 - cos(x)
>>> versin(x).rewrite(sin)
2*sin(x/2)**2
>>> versin(2*x).expand(trig=True)
2 - 2*cos(x)**2
>>> versin(a + b*I).expand(complex=True)
I*sin(a)*sinh(b) - cos(a)*cosh(b) + 1
```

**Differentiation:**

```
>>> versin(x).diff(x)
sin(x)
```

**Solving:**

<!-- Note: doit() is commented out in the example above because it makes the
     below return 1 - cos(x) instead of versin(x). Most people shouldn't
     implement doit() like that anyway, though. -->

(a more general version of `aversin` would have all the above methods defined
as well)

```
>>> class aversin(Function):
...     def inverse(self, argindex=1):
...         return versin
>>> from sympy import solve
>>> solve(aversin(x**2) - y, x)
[-sqrt(versin(y)), sqrt(versin(y))]
```

(custom-functions-divides-full-example)=
### divides

divides is a function defined by

$$\operatorname{divides}(m, n) = \begin{cases} 1 & \text{for}\: m \mid n \\
  0 & \text{for}\: m\not\mid n  \end{cases},$$

that is, `divides(m, n)` is 1 if `m` divides `n` and `0` if `m` does not
divide `m`. It is only defined for integer `m` and `n`. For the sake of
simplicity, we use the convention that $m \mid 0$ for all integer $m$.

`divides` is an example of a function that is only defined for certain input
values (integers). `divides` also gives an example of defining a custom
printer (`_latex()`).

#### Definition

```
>>> from sympy import Function, Integer
>>> from sympy.core.logic import fuzzy_not
>>> class divides(Function):
...     r"""
...     $$\operatorname{divides}(m, n) = \begin{cases} 1 & \text{for}\: m \mid n \\ 0 & \text{for}\: m\not\mid n  \end{cases}.$$
...
...     That is, ``divides(m, n)`` is ``1`` if ``m`` divides ``n`` and ``0``
...     if ``m`` does not divide ``n`. It is undefined if ``m`` or ``n`` are
...     not integers. For simplicity, the convention is used that
...     ``divides(m, 0) = 1`` for all integers ``m``.
...
...     References
...     ==========
...
...     .. [1] https://en.wikipedia.org/wiki/Divisor#Definition
...     """
...     # Define evaluation on basic inputs, as well as type checking that the
...     # inputs are not nonintegral.
...     @classmethod
...     def eval(cls, m, n):
...         # Evaluate for explicit integer m and n.
...         if isinstance(m, Integer) and isinstance(n, Integer):
...             return int(n % m == 0)
...
...         # For symbolic arguments, require m and n to be integer.
...         if m.is_integer is False or n.is_integer is False:
...             raise TypeError("m and n should be integers")
...
...     # Define basic assumptions.
...
...     # divides is always either 0 or 1.
...     is_integer = True
...     is_negative = False
...
...     # Whether divides(m, n) is 0 or 1 depends on m and n. Note that this
...     # method only makes sense because we don't automatically evaluate on
...     # such cases, but instead simplify these cases in doit() below.
...     def _eval_is_zero(self):
...         m, n = self.args
...         if m.is_integer and n.is_integer:
...              return fuzzy_not((n/m).is_integer)
...
...     # Define doit() as further evaluation on symbolic arguments using
...     # assumptions.
...     def doit(self, deep=False, **hints):
...         m, n = self.args
...         # Recursively call doit() on the args whenever deep=True.
...         # Be sure to pass deep=True and **hints through here.
...         if deep:
...            m, n = m.doit(deep=deep, **hints), n.doit(deep=deep, **hints)
...
...         # divides(m, n) is 1 iff n/m is an integer. Note that m and n are
...         # already assumed to be integers because of the logic in eval().
...         isint = (n/m).is_integer
...         if isint is True:
...             return Integer(1)
...         elif isint is False:
...             return Integer(0)
...         else:
...             return divides(m, n)
...
...     # Define LaTeX printing for use with the latex() function and the
...     # Jupyter notebook.
...     def _latex(self, printer):
...         m, n = self.args
...         _m, _n = printer._print(m), printer._print(n)
...         return r'\left [ %s \middle | %s \right ]' % (_m, _n)
...
```

#### Examples

**Evaluation**

```
>>> from sympy import symbols
>>> n, m, k = symbols('n m k', integer=True)
>>> divides(3, 10)
0
>>> divides(3, 12)
1
>>> divides(m, n).is_integer
True
>>> divides(k, 2*k)
divides(k, 2*k)
>>> divides(k, 2*k).is_zero
False
>>> divides(k, 2*k).doit()
1
```

**Printing:**

```py
>>> str(divides(m, n)) # This is using the default str printer
'divides(m, n)'
>>> print(latex(divides(m, n)))
\left [ m \middle | n \right ]
```


(custom-functions-fma-full-example)=
### Fused Multiply-Add (FMA)

[Fused Multiply-Add
(FMA)](https://en.wikipedia.org/wiki/Multiply%E2%80%93accumulate_operation#Fused_multiply%E2%80%93add)
is a multiplication followed by an addition:

$$\operatorname{FMA}(x, y, z) = xy + z.$$

It is often implemented in hardware as a single floating-point operation that
has better rounding and performance than the equivalent combination of
multiplication and addition operations.

FMA is an example of a custom function that is defined as an unevaluated
"shorthand" to another function. This is because the
[`doit()`](custom-functions-doit) method is defined to return `x*y + z`,
meaning the `FMA` function can easily be evaluated to the expression is
represents, but the [`eval()`](custom-functions-eval) method does *not* return
anything (except when `x`, `y`, and `z` are all explicit numeric values),
meaning that it stays unevaluated by default.

Contrast this with the
[versine](custom-functions-versine-full-example) example, which treats
`versin` as a first-class function in its own regard. Even though `versin(x)`
can be expressed in terms of other functions (`1 - cos(x)`) it does not
evaluate on general symbolic inputs in `versin.eval()`, and `versin.doit()` is
not defined at all.

`FMA` also represents an example of a continuous function defined on multiple
variables, which demonstrates how `argindex` works in the
[`fdiff`](custom-functions-differentiation) example.

Finally, `FMA` shows an example of defining some code printers for `C` and
`C++` (using the method names from {attr}`.C99CodePrinter.printmethod` and
{attr}`.CXX11CodePrinter.printmethod`), since that is a typical use-case for
this function.

The mathematical definition of FMA is very simple and it would be easy to
define every method on it, but only a handful are shown here. The
[versine](custom-functions-versine-full-example) and
[divides](custom-functions-divides-full-example) examples show how to define
the other important methods discussed in this guide.

Note that if you want to actually use fused-multiply add for code generation,
there is already a version in SymPy `sympy.codegen.cfunctions.fma()` which is
supported by the existing code printers. The version here is only designed to
serve as an example.

#### Definition

```py
>>> from sympy import Number, symbols, Add, Mul
>>> x, y, z = symbols('x y z')
>>> class FMA(Function):
...     """
...     FMA(x, y, z) = x*y + z
...
...     FMA is often defined as a single operation in hardware for better
...     rounding and performance.
...
...     FMA can be evaluated by using the doit() method.
...
...     References
...     ==========
...
...     .. [1] https://en.wikipedia.org/wiki/Multiply%E2%80%93accumulate_operation#Fused_multiply%E2%80%93add
...     """
...     # Define automatic evaluation on explicit numbers
...     @classmethod
...     def eval(cls, x, y, z):
...         # Number is the base class of Integer, Rational, and Float
...         if all(isinstance(i, Number) for i in [x, y, z]):
...            return x*y + z
...
...     # Define numerical evaluation with evalf().
...     def _eval_evalf(self, prec):
...         return self.doit(deep=False)._eval_evalf(prec)
...
...     # Define full evaluation to Add and Mul in doit(). This effectively
...     # treats FMA(x, y, z) as just a shorthand for x*y + z that is useful
...     # to have as a separate expression in some contexts and which can be
...     # evaluated to its expanded form in other contexts.
...     def doit(self, deep=True, **hints):
...         x, y, z = self.args
...         # Recursively call doit() on the args whenever deep=True.
...         # Be sure to pass deep=True and **hints through here.
...         if deep:
...             x = x.doit(deep=deep, **hints)
...             y = y.doit(deep=deep, **hints)
...             z = z.doit(deep=deep, **hints)
...         return x*y + z
...
...     # Define FMA.rewrite(Add) and FMA.rewrite(Mul).
...     def _eval_rewrite(self, rule, args, **hints):
...         x, y, z = self.args
...         if rule in [Add, Mul]:
...             return self.doit()
...
...     # Define differentiation.
...     def fdiff(self, argindex):
...         # argindex indexes the args, starting at 1
...         x, y, z = self.args
...         if argindex == 1:
...             return y
...         elif argindex == 2:
...             return x
...         elif argindex == 3:
...             return 1
...
...     # Define code printers for ccode() and cxxcode()
...     def _ccode(self, printer):
...         x, y, z = self.args
...         _x, _y, _z = printer._print(x), printer._print(y), printer._print(z)
...         return "fma(%s, %s, %s)" % (_x, _y, _z)
...
...     def _cxxcode(self, printer):
...         x, y, z = self.args
...         _x, _y, _z = printer._print(x), printer._print(y), printer._print(z)
...         return "std::fma(%s, %s, %s)" % (_x, _y, _z)
```

#### Examples

**Evaluation:**

```
>>> x, y, z = symbols('x y z')
>>> FMA(2, 3, 4)
10
>>> FMA(x, y, z)
FMA(x, y, z)
>>> FMA(x, y, z).doit()
x*y + z
>>> FMA(x, y, z).rewrite(Add)
x*y + z
>>> FMA(2, pi, 1).evalf()
7.28318530717959
```

**Differentiation**

```
>>> FMA(x, x, y).diff(x)
2*x
>>> FMA(x, y, x).diff(x)
y + 1
```

**Code Printers**

```
>>> from sympy import ccode, cxxcode
>>> ccode(FMA(x, y, z))
'fma(x, y, z)'
>>> cxxcode(FMA(x, y, z))
'std::fma(x, y, z)'
```

## Additional Tips

- SymPy includes dozens of functions. These can serve as useful examples for
  how to write a custom function, especially if the function is similar to one
  that is already implemented. Remember that everything in this guide applies
  equally well to functions that are included with SymPy and user-defined
  functions. Indeed, this guide is designed to serve as both a developer guide
  for contributors to SymPy and a guide for end-users of SymPy.

- If you have many custom functions that share common logic, you can use a
  common base class to contain this shared logic. For an example of this, see
  the [source code for the trigonometric functions in
  SymPy](https://github.com/sympy/sympy/blob/master/sympy/functions/elementary/trigonometric.py),
  which use `TrigonometricFunction`, `InverseTrigonometricFunction`, and
  `ReciprocalTrigonometricFunction` base classes with some shared logic.

- As with any code, it is a good idea to write extensive tests for your
  function. The [SymPy test
  suite](https://github.com/sympy/sympy/tree/master/sympy/functions/elementary/tests)
  is a good resource for examples of how to write tests for such functions.
  All code included in SymPy itself is required to be tested. Functions
  included in SymPy should also always contain a docstring with references, a
  mathematical definition, and doctest examples.
