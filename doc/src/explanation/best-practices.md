# Best Practices

This page outlines some of the best practices for users of SymPy. The best
practices here will help avoid some common bugs and pitfalls that can occur
when using SymPy.

(best-practices-interactive-vs-programmatic)=
## An Important Aside: Interactive vs. Programmatic Usage

There are two primary modes of usage for SymPy, and Python in general,
interactive and programmatic.

Interactive usage refers to live usage in an
interactive Python session such as a Python or [IPython](https://ipython.org/)
terminal session or [Jupyter Notebook](https://jupyter.org/). Interactive
sessions focus on live experimentation and fast iteration. It typically
involves reaching a solution in small iterative steps, viewing outputs often,
and going back and redefining parts of a calculation. Interactive sessions are
ephemeral in nature. The session history is either not saved at all, or even
if it is saved, it is not written or edited with any intention of rerunning
the code at an unspecified point in the future.

To contrast programmatic usage is any usage where the code saved with the
intention of being run again in the future. Typical interactive usage of SymPy
involves code that is saved to a `.py` file, or to a Jupyter notebook file
that is intended to be rerun. In programmatic usage, iterative steps are
consolidated to make code more readable, and abandoned ideas are deleted
entirely.

The distinction between interactive and programmatic usage matters because
many of the best practices outlined here may be ignored in interactive usage.
This is true of many general coding best practices, not just in SymPy. For
example, it is a best practice in Python to avoid wildcard imports, like `from
sympy import *`. These wildcard imports save typing, but they make it harder
for someone reading the code to understand where things come from, and they
can lead to bugs if wildcard imports from multiple libraries are used. But for
interactive usage, a wildcard import is acceptable, because it does save on
typing, and the ambiguity about what is actually imported can be resolved by
examining the namespace live in the interactive session.

SymPy has similar "shortcuts", which are perfectly acceptable in interactive
usage, but should be avoided for programmatic usage. Such instances will be
notated where appropriate. Other best practices listed here, which are not
notated as such, should be followed in all cases, as they can prevent pitfalls
and bugs even in interactive usage.

## Basic Usage

(best-practices-defining-symbols)=
### Defining Symbols

- **The best way to define symbols is using the {func}`~.symbols` function.**
  The `symbols()` function supports creating one or multiple symbols at a
  time:

  ```py
  >>> from sympy import symbols
  >>> x = symbols('x')
  >>> a, b, c = symbols('a b c')
  ```

  Additionally, it supports adding assumptions to symbols

  ```py
  >>> i, j, k = symbols('i j k', integer=True)
  ```

  and defining {class}`~.Function` objects:

  ```py
  >>> from sympy import Function
  >>> f, g, h = symbols('f g h', cls=Function)
  ```

  It also supports shorthands for defining many numbered symbols at once:

  ```py
  >>> symbols('x:10')
  (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9)
  ```

  There also exist some alternatives to {func}`~.symbols`. Two alternatives
  are using the {class}`~.Symbol` constructor directly and {mod}`sympy.abc`.
  These are both fine, but they are less general than {func}`~.symbols`, so
  just using `symbols()` is generally preferred. For example,
  {class}`~.Symbol` only supports creating one symbol at a time, and
  {mod}`sympy.abc` allows importing common single letter symbol names and
  doesn't allow including assumptions.

  The {func}`~.var` function should be avoided, except when working
  interactively. It works like the {func}`~.symbols` function, except it
  automatically injects symbol names into the calling namespace. This is
  designed only for typing convenience interactively, and shouldn't be used in
  programmatic environments (see
  [](best-practices-interactive-vs-programmatic) above).

- **Add assumptions to symbols when they are known.**
  [Assumptions](assumptions-guide) can be added by passing the relevant
  keywords to {func}`~.symbols`. The most common assumptions are `real=True`,
  `positive=True` (or `nonnegative=True`), and `integer=True`.

  Assumptions are never required, but it is always recommended to include them
  if they are known because it will allow certain operations to simplify. If
  no assumptions are provided, symbols are assumed to be general complex
  numbers, and simplifications will not be made unless they are true for all
  complex numbers.

  For example:


  ```py
  >>> from sympy import integrate, exp, oo
  >>> a = symbols('a') # no assumptions
  >>> integrate(exp(-a*x), (x, 0, oo))
  Piecewise((1/a, Abs(arg(a)) < pi/2), (Integral(exp(-a*x), (x, 0, oo)), True))
  ```

  ```py
  >>> a = symbols('a', positive=True)
  >>> integrate(exp(-a*x), (x, 0, oo))
  1/a
  ```

  Here, $\int_0^\infty e^{-ax}\,dx$ gives a piecewise result when `a`
  is defined with no assumptions, because the integral only converges when `a`
  is positive. Setting `a` to be positive removes this piecewise.

See also [](best-practices-avoid-string-inputs) and
[](best-practices-dont-hardcode-symbol-names) for related best practices
around defining symbols.

(best-practices-avoid-string-inputs)=
### Avoid String Inputs


Don't use strings as input to functions. Rather, create the objects
symbolically using Symbols and the appropriate SymPy functions, and manipulate
them.

**Don't**

```py
>>> from sympy import simplify
>>> simplify("(x**2 + x)/x")
x + 1
```

**Do**

```py
>>> from sympy import symbols
>>> x = symbols('x')
>>> simplify((x**2 + x)/x)
x + 1
```

It's always best to create expressions explicitly using Python operators, but
sometimes you really do start with a string input, like if you accept an
expression from the user. If you do have a string that you are starting with,
you should parse it explicitly with
[`parse_expr()`](sympy.parsing.sympy_parser.parse_expr). It is best to parse
all strings early and only use symbolic manipulation from there on.

```
>>> from sympy import parse_expr
>>> string_input = "(x**2 + x)/x"
>>> expr = parse_expr(string_input)
>>> simplify(expr)
x + 1
```

**Reason**

Support for string input is in many ways accidental. It only happens because
functions call `sympify()` on their input to ensure that it is a SymPy object,
and `sympify()` translates strings. Support for this may go away in a future
version of SymPy.

There are many disadvantages to using strings:

- They are not explicit. They make code much harder to read.

- `sympify()` automatically turns all undefined names into Symbols of
  Functions, so if you have a typo, the string will still parse correctly, but
  the output will not be what you expect. For example

  ```py
  >>> from sympy import expand_trig
  >>> expand_trig("sine(x + y)")
  sine(x + y)
  ```

  Compare this to the explicit error you get when you don't use strings.

  ```
  >>> from sympy import sin, symbols
  >>> x, y = symbols('x y')
  >>> expand_trig(sine(x + y)) # The typo is caught by a NameError
  Traceback (most recent call last):
  ...
  NameError: name 'sine' is not defined
  >>> expand_trig(sin(x + y))
  sin(x)*cos(y) + sin(y)*cos(x)
  ```

  In the first example, `sine`, a typo for `sin` is parsed into
  `Function("sine")`, and it appears that `expand_trig` cannot handle it. In the
  second case, we immediately get an error from the undefined name `sine`, and
  fixing our typo, we see that `expand_trig` can indeed do what we want.

- The biggest gotcha when using string inputs comes from using assumptions. In
  SymPy, if two symbols have the same name but different assumptions, they are
  considered unequal:

  ```
  >>> from sympy import Symbol
  >>> z1 = Symbol('z')
  >>> z2 = Symbol('z', positive=True)
  >>> z1 == z2
  False
  >>> z1 + z2
  z + z
  ```

  It is generally recommended to avoid doing this, as it can lead to confusing
  expressions like the one above.

  However, strings always parse symbols into symbols without assumptions. So
  if you have a symbol with an assumption and later try to use the string
  version of it, you will end up with confusing results.

  ```py
  >>> from sympy import diff
  >>> z = Symbol('z', positive=True)
  >>> diff('z**2', z)
  0
  ```

  The answer here is apparently wrong, but what happened is that the `z` in
  `"z**2"` parsed to `Symbol('z')` with no assumptions, which SymPy considers
  to be a different symbol from `Symbol('z', positive=True)`. So as far as
  `diff` is concerned, the expression is constant and the result is 0.

  This sort of thing is particularly bad because it generally doesn't lead to
  any errors. It will just silently give the "wrong" answer because SymPy will
  be treating symbols that you thought were the same as different. The
  solution is to not mix string inputs with non-string inputs.

  If you are parsing strings, and you want some of the symbols in it to have
  certain assumptions, you should create those symbols and pass them to the
  dictionary to [`parse_expr()`](sympy.parsing.sympy_parser.parse_expr). For example:

  **Don't**

  ```py
  >>> a, b, c = symbols('a b c', real=True)
  >>> # a, b, and c in expr are different symbols without assumptions
  >>> expr = parse_expr('a**2 + b - c')
  >>> expr.subs({a: 1, b: 1, c: 1})
  a**2 + b - c
  ```

  **Do**

  ```
  >>> # a, b, and c are the same as the a, b, c with real=True defined above
  >>> expr = parse_expr('a**2 + b - c', {'a': a, 'b': b, 'c': c})
  >>> expr.subs({a: 1, b: 1, c: 1})
  1
  ```

- Many operations on SymPy objects work as methods, not functions. These
  methods won't work on strings, since they are not yet SymPy objects. For
  example:

  ```
  >>> "x + 1".subs("x", "y")
  Traceback (most recent call last):
  ...
  AttributeError: 'str' object has no attribute 'subs'
  ```

  vs.

  ```
  >>> x, y = symbols('x y')
  >>> (x + 1).subs(x, y)
  y + 1
  ```

- Symbol names can contain any character, including things that aren't
  valid Python. But if you use strings as input, it is impossible to use such
  symbols. For example

  ```
  >>> from sympy import solve
  >>> solve('x_{2} - 1') # doctest: +SKIP
  ValueError: Error from parse_expr with transformed code: "Symbol ('x_' ){Integer (2 )}-Integer (1 )"
  ...
  SyntaxError: invalid syntax (<string>, line 1)
  ```

  This doesn't work because `x_{2}` is not valid Python.

  Symbol names do not have to be valid Python, if they are created from the
  `Symbol` or `symbols` functions. It's often convenient to use non-Python
  notation in Symbol names so that they will print nicely in LaTeX.

  ```py
  >>> x2 = Symbol('x_{2}')
  >>> solve(x2 - 1, x2)
  [1]
  ```

  This is the best case scenario, where you get an error. It is also possible
  you might get something unexpected:

  ```py
  >>> solve('x^1_2 - 1')
  [-1, 1, -I, I, -1/2 - sqrt(3)*I/2, -1/2 + sqrt(3)*I/2, 1/2 - sqrt(3)*I/2, 1/2 + sqrt(3)*I/2, -sqrt(3)/2 - I/2, -sqrt(3)/2 + I/2, sqrt(3)/2 - I/2, sqrt(3)/2 + I/2]
  ```

  What happened here is that instead of parsing `x^1_2` as $x^1_2$, it is
  parsed as `x**12`, because `^` is converted to `**` and [`_` is ignored in
  numeric literals in Python](https://peps.python.org/pep-0515/).

  If we instead create a Symbol, the actual contents of the symbol name are
  ignored. It is always represented as a single symbol.

  ```py
  >>> x12 = Symbol('x^1_2')
  >>> solve(x12 - 1, x12)
  [1]
  ```


- If you use strings, syntax errors won't be caught until the line is run. If
  you build up the expressions, syntax errors will be caught immediately when
  Python compiles the script before any of it runs.

- When you are using strings, it is tempting to build up expressions and do
  symbolic manipulation via manipulation of strings. This is almost always a
  bad idea. It is much better to use SymPy functions and operations to do
  symbolic manipulation. You will avoid mistakes, and doing things with
  strings will soon get too complex to handle. That's because there is no
  notion of a symbolic expression in a string. To Python, `"(x + y)/z"` is no
  different from `"/x+)(y z "`, which is the same string with the characters
  in another order. But SymPy expressions know what type of expression they
  are and how to traverse the expression tree (e.g., using `expr.args`), and
  you can use operators like `+-*/`. And as to the previous point, if there is
  an error, it will be difficult to track down because it won't be noticed
  until you actually parse the string into an expression.

  For example

  **Don't**

  ```py
  >>> expression_str = '+'.join([f'{i}*x_{i}' for i in range(10)])
  >>> expr = parse_expr(expression_str)
  >>> expr
  x_1 + 2*x_2 + 3*x_3 + 4*x_4 + 5*x_5 + 6*x_6 + 7*x_7 + 8*x_8 + 9*x_9
  ```

  **Do**

  ```py
  >>> from sympy import Add, Symbol
  >>> expr = Add(*[i*Symbol(f'x_{i}') for i in range(10)])
  >>> expr
  x_1 + 2*x_2 + 3*x_3 + 4*x_4 + 5*x_5 + 6*x_6 + 7*x_7 + 8*x_8 + 9*x_9
  ```

- In code editors that do syntax highlighting, strings will be highlighted all
  one color, whereas Python expressions will be highlighted according to their
  actual content.

- As mentioned, string inputs are not officially supported or tested, and so
  may go away at any time.

### Exact Rational Numbers vs. Floats

If a number is known to be exactly equal to some quantity, avoid defining it
as a floating-point number.

For example,

**Don't**

```py
>>> expression = x**2 + 0.5*x + 1
```

**Do**

```py
>>> expression = x**2 + x/2 + 1
```

However, this isn't to say that you should never use floating-point numbers in
SymPy, only that if a more exact value is known it should be preferred. SymPy
does support [arbitrary precision floating-point
numbers](sympy.core.numbers.Float), but some operations may not perform as
well with them.

One should also take care to avoid writing `integer/integer` where both
integers are explicit integers. This is because Python will evaluate this to a floating-point value before SymPy is able to parse it.

**Don't**

```py
>>> x + 2/7 # The exact value of 2/7 is lost
x + 0.2857142857142857
```

In this case, use {class}`~.Rational` to create a rational number, or use
`S()` shorthand if you want to save on typing.

**Do**

```py
>>> from sympy import Rational, S
>>> x + Rational(2, 7)
x + 2/7
>>> x + S(2)/7 # Equivalently
x + 2/7
```

This also applies to non-rational values which can be represented exactly. For
example, one should avoid using `math.pi` and prefer `sympy.pi`, since the
former is a numerical approximation to $\pi$ and the latter is exactly $\pi$
(see also [](best-practices-separate-sympy-and-non-sympy) below; in general,
one should avoid importing `math` when using SymPy).

**Don't**

```py
>>> import math
>>> import sympy
>>> math.pi
3.141592653589793
>>> sympy.sin(math.pi)
1.22464679914735e-16
```

**Do**

```py
>>> sympy.pi
pi
>>> sympy.sin(sympy.pi)
0
```

Here `sympy.sin(math.pi)` is not exactly 0, because `math.pi` is not exactly $\pi$.

**Reason**

Exact values, if they are known, should be preferred over floats for the
following reasons:

- An exact symbolic value can often be symbolically simplified or manipulated.
  A float represents an approximation to an exact real number, and therefore
  cannot be simplified exactly. For example, in the above example,
  `sin(math.pi)` does not simplify to `0` because `math.pi` is not exactly
  $\pi$. It is just a floating-point number that approximates $\pi$ to 15
  digits (effectively, a close rational approximation to $\pi$, but not
  exactly $\pi$).

- Some algorithms will not be able to compute a result if there are
  floating-point values, but can if the values are rational numbers. This is
  because rational numbers have properties that make it easier for these
  algorithms to work with them. For instance, with floats, one can have a
  situation where a number should be 0, but due to approximation errors, does
  not equal exactly 0.

  A particularly notable example of this is with floating-point exponents. For
  example,

  ```py
  >>> from sympy import factor
  >>> factor(x**2.0 - 1)
  x**2.0 - 1
  ```

- SymPy floats can have the same loss of significance cancellation issues that
  can occur from using finite precision floating-point approximations:

  ```py
  >>> from sympy import expand
  >>> expand((x + 1.0)*(x - 1e-16)) # the coefficient of x should be slightly less than 1
  x**2 + 1.0*x - 1.0e-16
  >>> expand((x + 1)*(x - Rational('1e-16'))) # Using rational numbers gives the coefficient of x exactly
  x**2 + 9999999999999999*x/10000000000000000 - 1/10000000000000000
  ```

  It is possible to avoid these issues in SymPy in many cases by making
  careful use of `evalf` with its ability to evaluate in arbitrary precision.
  This typically involves either computing an expression with symbolic values
  and substituting them later with `expr.evalf(subs=...)`, or by starting with
  `Float` values with a precision higher than the default of 15 digits:

  ```
  >>> from sympy import Float
  >>> expand((x + 1.0)*(x - Float('1e-16', 20)))
  x**2 + 0.9999999999999999*x - 1.0e-16
  ```

A `Float` number can be converted to its exact rational equivalent by passing
it to `Rational`. Alternatively, you can use `nsimplify` to find the nicest
rational approximation. This may often reproduce the number that was intended
if the number is supposed to be rational (although again, it's best to just
start with rational numbers in the first place, if you can):

```py
>>> from sympy import nsimplify
>>> Rational(0.7)
3152519739159347/4503599627370496
>>> nsimplify(0.7)
7/10
```

### Avoid `simplify()`

{func}`~.simplify` is designed as a general purpose heuristic. It tries
various simplification algorithms on the input expression and returns the
result that seems the "simplest" based on some metric.

`simplify()` is perfectly fine for interactive use, where you just want SymPy
to do whatever it can to an expression. However, in programmatic usage, it's
better to avoid `simplify()` and use more targeted simplification functions
like {func}`~.cancel` or {func}`~.collect` instead.

There are a few reasons why this is generally preferred:

- Due to its heuristical nature, `simplify()` can potentially be slow, since
  it tries a lot of different approaches to try to find the best
  simplification.

- There are no guarantees about what form an expression will have after being
  passed through `simplify()`. It may actually end up "less simple" by
  whatever metric you were hoping for. To contrast, targeted simplification
  functions are very specific about what behaviors they have and what they
  guarantee about the output. For example,

  - {func}`~.factor` will always factor a polynomial into irreducible factors.

  - {func}`~.cancel` will always convert a rational function into the form
    $p/q$ where $p$ and $q$ are expanded polynomials with no common factors.

  - The {func}`~.trigsimp` function is similarly heuristical in nature
    (although targeted to only trigonometric functions), but the routines in
    the {mod}`sympy.simplify.fu` submodule allow applying specific
    trigonometric identities.

- A targeted simplification will not do something unexpected if the expression
  contains an unexpected form, or an unexpected subexpression. This is
  especially the case if simplification functions are applied with `deep=False`
  to only apply the simplification to the top-level expression.

The [simplify section of the tutorial](tutorial-simplify) and the
[simplify module reference](../modules/simplify/simplify.rst) list the various
targeted simplification functions.

In some cases, you may know exactly what simplification operations you wish to
apply to an expression, but there may not be an exact set of simplification
functions that do them. When this happens, you can create your own targeted
simplification using {meth}`~sympy.core.basic.Basic.replace`, or in general, manually using
[advanced expression manipulation](tutorial-manipulation).

<!-- TODO: Add an example -->


(best-practices-dont-hardcode-symbol-names)=
### Don't Hardcode Symbol Names in Functions

Instead of hard-coding {class}`~.Symbol` names inside of a function
definition, make the symbols a parameter to the function.

**Don't**

```py
def theta_operator(expr):
    t = symbols('t')
    return t*expr.diff(t)
```

**Do**

```py
def theta_operator(expr, t):
    return t*expr.diff(t)
```

A hard-coded symbol name has the disadvantage of requiring all expressions to
use that exact symbol name. In the above example, it is not possible to
compute $\theta = xD_x$ because it is hard-coded to $tD_t$. What's worse,
trying to do so silently leads to a wrong result instead of an error, since
`x` is treated as a constant expression:

```py
>>> def theta_operator(expr):
...     t = symbols('t')
...     return t*expr.diff(t)
>>> theta_operator(x**2)
0
```

This is particularly problematic if the function accepts arbitrary user input,
as the user may be using a different variable name that makes more sense in
their mathematical context. And if the user already used the symbol `t` but as
a constant, they would need to swap things around with `subs` before being
able to use the function.

The other reason this antipattern is problematic is due to the gotcha that
symbols with assumptions are considered unequal to symbols without
assumptions. If someone defined their expression using

```py
>>> t = symbols('t', positive=True)
```

for example, to make further simplifications possible (see
[](best-practices-defining-symbols) above), the function hard-coding
`Symbol('t')` without assumptions would not work:

```py
>>> theta_operator(t**2)
0
```

By making the symbol an argument, like `theta_operator(expr, t)`, these
problems all go away.

(best-practices-separate-sympy-and-non-sympy)=
### Separate SymPy and non-SymPy Numerical Code

TODO

## Advanced Usage

### Be Careful Comparing and Sorting Symbolic Objects

Be careful with programmatic code that compares numerical quantities, whether
using an inequality (`<`, `<=`, `>`, `>=`) or indirectly with something like
`sort`. The issue is that if an inequality is unknown, it will produce a
symbolic inequality, like

```
>>> x > 0
x > 0
```

A symbolic inequality will raise an exception if `bool()` is called on it, due
to the ambiguity:

```py
>>> bool(x > 0)
Traceback (most recent call last):
...
TypeError: cannot determine truth value of Relational
```

A check like

```py
if x > 0:
    ...
```

May work just fine if you only ever test it for numerical `x`. But if `x` can
ever be symbolic, the above code is wrong. It will fail with `TypeError:
cannot determine truth value of Relational`. If you ever see this exception,
it means this error has been made somewhere (sometimes the error is in SymPy
itself; if this appears to be the case, please [open an
issue](https://github.com/sympy/sympy/issues)).

The exact same issue occurs when using `sorted`, since this internally uses `>`.

```py
>>> sorted([x, 0])
Traceback (most recent call last):
...
TypeError: cannot determine truth value of Relational
```

There are a few options for fixing this issue, and the correct one to choose
depends on what you are doing:

- **Disallow symbolic inputs.** If your function cannot possibly work on
  symbolic inputs, you can explicitly disallow them. The primary benefit here
  is to give a more readable error message to users. The
  {attr}`~sympy.core.expr.Expr.is_number` can be used to check expression can
  be evaluated to a specific number with `evalf()`. If you want to only accept
  integers, you can check `isinstance(x, Integer)` (after calling `sympify`)
  (beware that `is_integer` uses the assumptions system and may be True even
  for symbolic objects, like `Symbol('x', integer=True)`).

- **Use the assumptions system.** If you do support symbolic inputs, you
  should use the assumptions system to check for things like `x > 0`, e.g.,
  using `x.is_positive`. When doing this, you should always [be aware of the
  nuances](booleans-guide) of the {term}`three-valued fuzzy logic
  <three-valued logic>` used in the assumptions system. That is, always be
  aware that an assumption could be `None`, meaning its value is unknown and
  could be either true or false. For example,

  ```py
  if x.is_positive:
      ...
  ```

  will only run the block if `x.is_positive` is `True`, but you may want to do
  something when `x.is_positive` is `None`.

- **Return a Piecewise result.** If the result of a function depends on an
  inequality or other boolean condition, you can use {class}`~.Piecewise` to
  return a result that applies to both possibilities symbolically. This is
  generally preferred when possible, as it offers the most flexibility. This
  is because the result is represented symbolically, meaning, for instance,
  one can later substitute specific values for the symbols and it will
  evaluate to the specific case, even if it is combined with other
  expressions.

  For example, instead of

  ```py
  if x > 0:
      expr = 1
  else:
      expr = 0
  ```

  this can be represented symbolically as

  ```py
  >>> from sympy import Piecewise
  >>> expr = Piecewise((1, x > 0), (0, True))
  >>> expr
  Piecewise((1, x > 0), (0, True))
  >>> expr.subs(x, 1)
  1
  >>> expr.subs(x, -1)
  0
  ```

- **Use {func}`~.ordered` to sort expressions into a canonical order.** If you
  are trying to use `sorted` because you want a canonical ordering, but you
  don't particularly care what that ordering is, you can use `ordered`.

  ```py
  >>> from sympy import ordered
  >>> list(ordered([x, 0]))
  [0, x]
  ```

  Alternatively, try to write the function in a way so that the result does not depend on
  the order that arguments are processed in.

## Custom SymPy Objects

### Args Invariants

Custom SymPy objects should always satisfy the following invariants:

1. `all(isinstance(arg, Basic) for arg in args)`
2. `expr.func(*expr.args) == expr`

The first says that all elements of {term}`args` should be instances of
{term}`Basic`. The second says that an expression should be rebuildable from
its `args` (note that {term}`func` is usually the same as `type(expr)`).

These two invariants are assumed throughout SymPy, and are essential for any
function that manipulates expressions.

For example, consider this simple function, which is a simplified version of
{meth}`~sympy.core.basic.Basic.xreplace`:

```py
>>> def replace(expr, x, y):
...     """Replace x with y in expr"""
...     newargs = []
...     for arg in expr.args:
...         if arg == x:
...             newargs.append(y)
...         else:
...             newargs.append(replace(arg, x, y))
...     return expr.func(*newargs)
>>> replace(x + sin(x - 1), x, y)
y + sin(y - 1)
```

The function works by recursively traversing the `args` of `expr`, and
rebuilding it except any instances of `x` are replaced by `y`.

It's easy to see how this function would break if the args invariants did not
hold:

1. If an expression has args that are not `Basic`, they will fail with
   `AttributeError` on a recursive call, because the non-`Basic` args will not
   have the `args` or `func` attributes.

   Usually this means calling `_sympify()` on the inputs to the class so that
   they are basic instances. If you want to store a string on a class, you
   should either use a `Symbol` or `sympy.core.symbols.Str`.

2. If an expression does not rebuild from its `args`, the line `return
   exr.func(*newargs)`, which should be a no-op if none of the args are
   changed by the replacement, will fail.

   In some cases a class may accept args in multiple equivalent forms. It is
   important that whatever form is stored in `args` is one of the ways that
   can be used to reconstruct the class.

Note that most user-defined custom functions should be defined by subclassing
`Function` (see the [guide to writing custom functions](custom-functions)).
The `Function` class automatically takes care of both of the args invariants,
so if you are using it, you do not need to worry about this.

### Avoid Too Much Automatic Evaluation

When defining a custom function, avoid doing too much automatic evaluation
(i.e., evaluation in the `eval` or `__new__` methods).

Generally, automatic evaluation should only be done in instances where it is
fast, and it is something that no one ever want to not happen. Automatic
evaluation is difficult to undo. A good rule of thumb is to evaluate on
explicit numeric values (`isinstance(x, Number)`), and leave everything else
symbolically unevaluated.
Further simplification using more advanced identities should be done in
specific simplification functions or `doit` (see the [custom functions
guide](custom-functions) for a list of common simplification routines that can
be defined on SymPy objects).

The [custom functions guide](custom-functions-automatic-evaluation) goes over
this in depth (but note that this guideline applies equally to all SymPy
objects, not just functions). But in a nutshell, the reason for this is that the only way to prevent automatic evaluation is to use
`evaluate=False`, which is fragile. Additionally, code will invariably be
written assuming the invariants that are true due to automatic evaluations,
meaning that expressions created with `evaluate=False` can lead to wrong
results from this code.

Evaluation that can potentially be expensive (for
instance, applying a symbolic identity) is itself bad because it can make
creating an expression without even doing anything with it allow. This can
also apply to checking for symbolic assumptions (like `x.is_integer`), so it
is also better to avoid this.

**Don't**

```py
class f(Function):
    @classmethod
    def eval(cls, x):
        if x.is_integer: # Bad (checking general assumptions)
            return 0
        if isinstance(x, Add): # Bad (applying symbolic identities)
            return Add(*[f(i) for i in x.args])

```

**Do**

```
class f(Function):
    @classmethod
    def eval(cls, x):
        if isinstance(x, Integer): # Good (only evaluating on explicit integers)
            return 0
    # Good (applying simplification on assumptions in doit())
    def doit(self, deep=True, **hints):
        x = self.args[0]
        if deep:
           x = x.doit(deep=deep, **hints)
        if x.is_integer:
           return S(0)
        return self
    # Good (applying symbolic identities inside of simplification functions)
    def _eval_expand_func(self, **hints):
        x = self.args[0]
        if isinstance(x, Add):
            return Add(*[f(i) for i in x.args])
        return self
```

Note that not all the classes in SymPy currently follow this guideline very
well, but it is something that we are improving.

### Don't Denest Collections

TODO

(best-practices-extra-attributes)=
### Avoid Storing Extra Attributes on an Object

A common reason that you might want to create a custom SymPy object is that
you want to store extra attributes on the object. However, doing this in a
naive way, i.e., by simply storing the data as a Python attribute on the
object, is almost always a bad idea.

SymPy does not expect objects to have extra data stored in them beyond what is
in their {term}`args`. For instance, this breaks `==` checking, which only
compares an objects `args`. See the [](best-practices-eq) section below for
why it is a bad idea to override `__eq__`. This section and that one are
closely related.

Typically, there is a better way to do what you are trying to do, depending on
the specific details of your situation:

- **Store the extra data in the object's `args`.** This is the best approach if
  the extra data you want to store is part of the *mathematical* description
  of your object.

  As long as the data is representable using other SymPy objects, it can be
  stored in `args`. Note that an object's `args` should be usable to recreate
  the object (e.g., something like `YourObject(*instance.args)` should
  recreate `instance`).

  Additionally, it should be mentioned that it is not a good idea to subclass
  `Symbol` if you plan to store anything extra in `args`. `Symbol` is designed
  around having no `args`. You are better off subclassing `Function` (see
  [](custom-functions)) or `Expr` directly. If you simply want to have two
  symbols that are distinct from one another, the best approach is often just
  to give them different names. If you are concerned about how they are
  printed, you can replace them with a more canonical name when it comes time
  to print things, or use a [custom printer](module-printing).

- **Store the data about the object separately.** This is the best approach if
  the extra data is not directly related to an objects mathematical
  properties.

  Remember that SymPy objects are hashable, so they can easily be used as
  dictionary keys. So maintaining a separate dictionary of `{object:
  extra_data}` pairs is straightforward.

  Note that some SymPy APIs already allow redefining how they operate on
  objects separately from the objects themselves. A big example of this is the
  {term}`printers <printing>`, which allow defining [custom
  printers](module-printing) that change how any SymPy object is printed
  without modifying those object themselves. Functions like {func}`~.lambdify`
  and {func}`~.init_printing` allow passing in a custom printer.

- **Represent the attribute using different subclasses.** This is often a good
  idea if there are only a few possible values for the attribute (e.g., a
  boolean flag). Code duplication can be avoided by using a common superclass.

- **If the data you want to store is a Python function**, it's best to just
  use as a method on the class. In many cases, the method may already fit into
  one of the [existing set of overridable SymPy methods](custom-functions). If
  you want to define how a function evaluates itself numerically, you can use
  {func}`~.implemented_function`.

- **Represent the information using by modifying the object's `func`.** This
  solution is much more complicated than the others, and should only be used
  when it is necessary. In some extreme cases, it is not possible to represent
  every mathematical aspect of an object using `args` alone. This can happen,
  for example, because of the limitation that `args` should only contain
  `Basic` instances. It is still possible to create custom SymPy objects in
  these situations by using a custom {term}`func` that is different from
  `type(expr)` (in this case, you would override `__eq__` on the `func`
  [rather than on the class](best-practices-eq)).

  However, this sort of situation is rare.

(best-practices-eq)=
### Don't Overwrite `__eq__`

When building a custom SymPy object, it is sometimes tempting to overwrite
`__eq__` to define custom logic for the `==` operator. This is almost always a
bad idea. Custom SymPy classes should leave `__eq__` undefined and use the
default implementation in the `Basic` superclass.

In SymPy, `==` compares objects using {term}`structural equality`. That is, `a
== b` means that `a` and `b` are exactly the same object. They have the same
type and the same {term}`args`. `==` does not perform any sort of
*mathematical* equality checking. For example,

```py
>>> x*(x - 1) == x**2 - x
False
```

`==` also always returns a boolean `True` or `False`. Symbolic equations can
be represented with {class}`Eq <sympy.core.relational.Equality>`.

There are several reasons for this

- Mathematical equality checking can be very expensive to compute, and in
  general, it is [computationally impossible to
  determine](https://en.wikipedia.org/wiki/Richardson%27s_theorem).

- Python itself automatically uses `==` in various places and assumes that it
  returns a boolean and is inexpensive to compute. For example, `a in b` uses
  `==`, where `b` is a builtin Python container like `list`, `dict`, or
  `set`.[^dict-footnote]

[^dict-footnote]: Python dicts and sets use `hash`, but fallback to using `==`
    when there is a hash collision.

- SymPy internally uses `==` all over the place, both explicitly and
  implicitly via things like `in` or dictionary keys. This
  usage all implicitly assumes that `==` operates structurally.


In affect, *structural equality* means that if `a == b` is `True`, then `a`
and `b` are for all intents and purposes the same object. This is because all
SymPy objects are {term}`immutable`. Any SymPy function may freely `a` with
`b` in any subexpression.

The default `__eq__` method on {term}`Basic` checks if the two objects have
the same type and the same `args`. There are also many parts of SymPy that
implicitly assume that if two objects are equal, then they have the same
`args`. Therefore, it is not a good idea to try to override `__eq__` as a way
to avoid storing some identifying information about an object in its `args`.
The `args` of an object should contain everything that is needed to recreate
it (see {term}`args`). Note that it is possible for an objects constructor to
accept multiple forms of arguments, so long as it accepts the form stored in
`args` (e.g., it is perfectly fine for some args to have default values).

Here are some examples of reasons you might be tempted to override `__eq__`
and the preferred alternatives:

- To make `==` apply some smarter equality check than purely structural
  equality. As noted above, this is a bad idea because too many things
  implicitly assume `==` works structurally only. Instead, use a function or
  method to implement the smarter equality (for example, the `equals` method).

  Another option is to define a {term}`canonicalization <canonicalize>` method
  that puts objects into canonical form (e.g., via `doit`), so that, for
  instance, `x.doit() == y.doit()` is true whenever `x` and `y` are
  mathematically equal. This is not always possible because not every type of
  object has a commutable canonical form, but it is a convenient approach when
  one does exist.

- To make `==` check for some additional attributes beyond those stored in the
  `args` of an expression. [](best-practices-extra-attributes) for more
  details on why it's a bad idea to directly store extra attributes on a SymPy
  object, and what the best alternatives are.

- To make `==` compare equal to some non-SymPy object. It is preferable to
  extend `sympify` to be able to convert this object into the SymPy object.
  The default `__eq__` implementation will automatically call `sympify` on the
  other argument if it isn't a `Basic` instance (e.g., `Integer(1) == int(1)`
  gives `True`). It is possible to extend `sympify` both for objects you
  control by defining a `_sympy_` method and for objects you do not control by
  extending the `converter` dictionary. See the {func}`~.sympify`
  documentation for more details.

### Avoiding Infinite Recursion from Assumptions Handlers

When writing assumptions handlers on custom functions like `_eval_is_positive`
(see the [custom
functions guide](custom-functions-assumptions) for details on how to do this),
there are two important things to keep in mind:

**Firstly, avoid creating new expressions inside of an assumption handler. You
should always pull apart the arguments of a function directly instead.** The
reason is that creating a new expression could itself result in an assumptions
query. This can easily lead to infinite recursion. And even when it doesn't,
creating a new expression which itself could lead to many recursive
assumptions queries is bad for performance compared to querying the desired
property more directly.

This generally means using methods like {meth}`~.as_independent` and checking
the `args` of expressions directly (see the [custom functions
guide](custom-functions-assumptions) for an example).

<!-- TODO: Give an example here. Can we show an infinite recursion? -->

Secondly, **do not recursively evaluate assumptions on `self` in assumptions
handlers**. Assumptions handlers should only check for assumptions on
`self.args`. The global assumptions system will automatically handle
implications between different assumptions.

For example, you may be tempted to write something like

```py
# BAD

class f(Function):
    def _eval_is_integer(self):
        # Quick return if self is not real (do not do this).
        if self.is_real is False:
            return False
        return self.args[0].is_integer
```

However, the `if self.is_real is False` check is completely unnecessary. The
assumptions system already knows that `integer` implies `real`, and it will
not bother checking `is_integer` if it already knows that `is_real` is False.

If you define the function this way, it will lead to an infinite recursion:

```py
>>> class f(Function):
...     def _eval_is_integer(self):
...         if self.is_real is False:
...             return False
...         return self.args[0].is_integer
>>> f(x).is_real
Traceback (most recent call last):
...
RecursionError: maximum recursion depth exceeded while calling a Python object
```

Instead, define the handler based on the arguments of the function only:

```
# GOOD

class f(Function):
    def _eval_is_integer(self):
        return self.args[0].is_integer
```
