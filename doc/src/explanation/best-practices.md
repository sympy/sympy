# Best Practices

This page outlines some of the best practices for users of SymPy. The best
practices here will help avoid some common bugs and pitfalls that can occur
when using SymPy.

This page primarily focuses on best practices that apply generally to all
parts of SymPy. Best practices that are specific to certain SymPy submodules
or functions are outlined in the documentation for those specific functions.

## Basic Usage

(best-practices-defining-symbols)=
### Defining Symbols

- **Define symbols with {func}`~.symbols` or {class}`~.Symbol()`.** The
  `symbols()` function is the most convenient way to create symbols. It
  supports creating one or more symbols at once:

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

  The `Symbol()` constructor may also be used directly. Unlike `symbols()`,
  `Symbol()` always creates one symbol. It is the best option if you want to
  make a symbol with unusual characters in its name or if you are creating
  symbols programmatically.

  ```py
  >>> from sympy import Symbol
  >>> x_y = Symbol('x y') # This creates a single symbol named 'x y'
  ```

  The {func}`~.var` function should be avoided, except when working
  interactively. It works like the {func}`~.symbols` function, except it
  automatically injects symbol names into the calling namespace. This function
  is designed solely for interactive typing convenience and is not recommended
  for programmatic use.

  Do not use `sympify()` or `S()` to create symbols. This may appear to work:

  ```py
  >>> from sympy import S
  >>> x = S("x") # DO NOT DO THIS
  ```

  However, `S()`/`sympify()` are not designed to create symbols. They are
  designed to parse entire expressions. This method fails if the input string
  is not valid Python. It also fails if the string parses to a larger
  expression:

  ```py
  >>> # These both fail
  >>> x = S("0x") # doctest: +SKIP
  Traceback (most recent call last):
  ...
  SyntaxError: invalid syntax (<string>, line 1)
  >>> x = S("x+") # doctest: +SKIP
  Traceback (most recent call last):
  ...
  SyntaxError: invalid syntax (<string>, line 1)
  ```

  Any Python string can be used as a valid Symbol name.

  Furthermore, all the same issues described in the
  [](best-practices-avoid-string-inputs) section below apply here.

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

  When you do use assumptions, the best practice is to always use the same
  assumptions for each symbol name. SymPy allows the same symbol name to be
  defined with different assumptions, but these symbols will be considered
  unequal to each other:

  ```py
  >>> z1 = symbols('z')
  >>> z2 = symbols('z', positive=True)
  >>> z1 == z2
  False
  >>> z1 + z2
  z + z
  ```

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
>>> from sympy import expand
>>> expand("(x**2 + x)/x")
x + 1
```

**Do**

```py
>>> from sympy import symbols
>>> x = symbols('x')
>>> expand((x**2 + x)/x)
x + 1
```

It's always best to create expressions explicitly using Python operators, but
sometimes you really do start with a string input, like if you accept an
expression from the user. If you do have a string that you are starting with,
you should parse it explicitly with
[`parse_expr()`](sympy.parsing.sympy_parser.parse_expr). It is best to parse
all strings early and only use symbolic manipulation from there on.

```py
>>> from sympy import parse_expr
>>> string_input = "(x**2 + x)/x"
>>> expr = parse_expr(string_input)
>>> expand(expr)
x + 1
```

**Reason**

There are many disadvantages to using strings as input to SymPy functions:

- It is unpythonic and makes code harder to read. See [the Zen of
  Python](https://peps.python.org/pep-0020/) "explicit is better than implicit".

- Support for string inputs in general SymPy functions is mostly accidental.
  It happens because these functions call {func}`~.sympify` on their inputs in
  order to convert things like Python `int`s into SymPy `Integer`s. However,
  `sympify()` also parses strings into SymPy expressions, unless the
  `strict=True` flag is used. Automatic parsing of strings for general SymPy
  functions (other than `sympify()` or {func}`~.parse_expr()`) [may go away in
  a future version of SymPy](https://github.com/sympy/sympy/issues/11003).

- Typos in symbol or function names can go unnoticed. This is because all
  undefined names in the string will be automatically parsed into Symbols or
  Functions. If the input has a typo, the string will still parse correctly,
  but the output will not be what was expected. For example

  ```py
  >>> from sympy import expand_trig
  >>> expand_trig("sine(x + y)")
  sine(x + y)
  ```

  Compare this to the explicit error you get when not using strings:

  ```py
  >>> from sympy import sin, symbols
  >>> x, y = symbols('x y')
  >>> expand_trig(sine(x + y)) # The typo is caught by a NameError
  Traceback (most recent call last):
  ...
  NameError: name 'sine' is not defined
  >>> expand_trig(sin(x + y))
  sin(x)*cos(y) + sin(y)*cos(x)
  ```

  In the first example, `sine`, a typo for `sin`, is parsed into
  `Function("sine")`, and it appears that `expand_trig` cannot handle it. In the
  second case, we immediately get an error from the undefined name `sine`, and
  fixing our typo, we see that `expand_trig` can indeed do what we want.

- The biggest gotcha when using string inputs comes from using assumptions. In
  SymPy, if two symbols have the same name but different assumptions, they are
  considered unequal:

  ```py
  >>> z1 = symbols('z')
  >>> z2 = symbols('z', positive=True)
  >>> z1 == z2
  False
  >>> z1 + z2
  z + z
  ```

  It is generally recommended to avoid doing this, as it can lead to confusing
  expressions like the one above (see [](best-practices-defining-symbols)
  above).

  However, string inputs will always create symbols without assumptions. So
  if you have a symbol with an assumption and later try to use the string
  version of it, you will end up with confusing results.

  ```py
  >>> from sympy import diff
  >>> z = symbols('z', positive=True)
  >>> diff('z**2', z)
  0
  ```

  The answer here is apparently wrong, but what happened is that the `z` in
  `"z**2"` parsed to `Symbol('z')` with no assumptions, which SymPy considers
  to be a different symbol from `z = Symbol('z', positive=True)`, which is
  used as the second argument to `diff()`. So as far as `diff` is concerned,
  the expression is constant and the result is 0.

  This sort of thing is particularly bad because it generally doesn't lead to
  any errors. It will just silently give the "wrong" answer because SymPy will
  be treating symbols that you thought were the same as different. The
  situation is avoided by not using string inputs.

  If you are parsing strings, and you want some of the symbols in it to have
  certain assumptions, you should create those symbols and pass them to the
  dictionary to [`parse_expr()`](sympy.parsing.sympy_parser.parse_expr). For example:

  **Don't**

  ```py
  >>> a, b, c = symbols('a b c', real=True)
  >>> # a, b, and c in expr are different symbols without assumptions
  >>> expr = parse_expr('a**2 + b - c')
  >>> expr.subs({a: 1, b: 1, c: 1}) # The substitution (apparently) doesn't work
  a**2 + b - c
  ```

  **Do**

  ```py
  >>> # a, b, and c are the same as the a, b, c with real=True defined above
  >>> expr = parse_expr('a**2 + b - c', {'a': a, 'b': b, 'c': c})
  >>> expr.subs({a: 1, b: 1, c: 1})
  1
  ```

- Many SymPy operations are defined as methods, not functions, that is, they
  are called like `sympy_obj.method_name()`. These methods won't work on
  strings, since they are not yet SymPy objects. For example:

  ```py
  >>> "x + 1".subs("x", "y")
  Traceback (most recent call last):
  ...
  AttributeError: 'str' object has no attribute 'subs'
  ```

  Contrasted with:

  ```py
  >>> x, y = symbols('x y')
  >>> (x + 1).subs(x, y)
  y + 1
  ```

- Symbol names can contain any character, including things that aren't
  valid Python. But if you use strings as input, it is impossible to use such
  symbols. For example

  ```py
  >>> from sympy import solve
  >>> solve('x_{2} - 1') # doctest: +SKIP
  ValueError: Error from parse_expr with transformed code: "Symbol ('x_' ){Integer (2 )}-Integer (1 )"
  ...
  SyntaxError: invalid syntax (<string>, line 1)
  ```

  This doesn't work because `x_{2}` is not valid Python. But it is perfectly
  possible to use this as a Symbol name:

  ```py
  >>> x2 = symbols('x_{2}')
  >>> solve(x2 - 1, x2)
  [1]
  ```

  Actually, the above is the best case scenario, where you get an error. It is
  also possible you might get something unexpected:

  ```py
  >>> solve('x^1_2 - 1')
  [-1, 1, -I, I, -1/2 - sqrt(3)*I/2, -1/2 + sqrt(3)*I/2, 1/2 - sqrt(3)*I/2, 1/2 + sqrt(3)*I/2, -sqrt(3)/2 - I/2, -sqrt(3)/2 + I/2, sqrt(3)/2 - I/2, sqrt(3)/2 + I/2]
  ```

  What happened here is that instead of parsing `x^1_2` as $x^1_2$, it is
  parsed as `x**12` (`^` is converted to `**` and [`_` is ignored in
  numeric literals in Python](https://peps.python.org/pep-0515/)).

  If we instead create a Symbol, the actual contents of the symbol name are
  ignored. It is always represented as a single symbol.

  ```py
  >>> x12 = symbols('x^1_2')
  >>> solve(x12 - 1, x12)
  [1]
  ```

- If you use strings, syntax errors won't be caught until the line is run. If
  you build up the expressions, syntax errors will be caught immediately by
  before any of it runs.

- Syntax highlighting in code editors doesn't typically recognize and
  color-code the content of strings, whereas it can recognize Python
  expressions.

### Avoid Manipulating Expressions as Strings

If you find yourself doing a lot of string or regular expression manipulations
on symbolic expressions, this is generally a sign that you are using SymPy
incorrectly. It's better to build up expressions directly with operators like
`+`, `-`, `*`, and `/` and SymPy's various functions and methods. String-based
manipulations can introduce errors, grow complex quickly, and lose the
benefits of symbolic expression structures.

The reason for this is that there is no notion of a symbolic expression in a
string. To Python, `"(x + y)/z"` is no different from `"/x+)(y z "`, which is
the same string with the characters in another order. To contrast, a SymPy
expression actually knows about what type of mathematical object it
represents. SymPy has many methods and functions for building and manipulating
expressions, and they all operate on SymPy objects, not strings.

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

See also the [previous section on avoiding string inputs to
functions](best-practices-avoid-string-inputs).

(best-practices-exact-rational-numbers-vs-floats)=
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
>>> from sympy import Rational
>>> expression = x**2 + Rational(1, 2)*x + 1
>>> expression = x**2 + x/2 + 1 # Equivalently
```

However, this isn't to say that you should never use floating-point numbers in
SymPy, only that if a more exact value is known it should be preferred. SymPy
does support [arbitrary precision floating-point
numbers](sympy.core.numbers.Float), but some operations may not perform as
well with them.

This also applies to non-rational numbers which can be represented exactly. For
example, one should avoid using `math.pi` and prefer `sympy.pi`, since the
former is a numerical approximation to $\pi$ and the latter is exactly $\pi$
(see also [](best-practices-separate-symbolic-and-numeric-code) below; in
general, one should avoid importing `math` when using SymPy).

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


One should also take care to avoid writing `integer/integer` where both
integers are explicit integers. This is because Python will evaluate this to a
floating-point value before SymPy is able to parse it.

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

**Reason**

Exact values, if they are known, should be preferred over floats for the
following reasons:

- An exact symbolic value can often be symbolically simplified or manipulated.
  A float represents an approximation to an exact real number, and therefore
  cannot be simplified exactly. For example, in the above example,
  `sin(math.pi)` does not produce `0` because `math.pi` is not exactly $\pi$.
  It is just a floating-point number that approximates $\pi$ to 15 digits
  (effectively, a close rational approximation to $\pi$, but not exactly
  $\pi$).

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
  >>> factor(x**2 - 1)
  (x - 1)*(x + 1)
  ```

- SymPy Floats have the same loss of significance cancellation issues that can
  occur from using finite precision floating-point approximations:

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
rational approximation. This can sometimes reproduce the number that was
intended if the number is supposed to be rational (although again, it's best
to just start with rational numbers in the first place, if you can):

```py
>>> from sympy import nsimplify
>>> Rational(0.7)
3152519739159347/4503599627370496
>>> nsimplify(0.7)
7/10
```

### Avoid `simplify()`

{func}`~.simplify` (not to be confused with {func}`~.sympify`) is designed as
a general purpose heuristic. It tries various simplification algorithms on the
input expression and returns the result that seems the "simplest" based on
some metric.

`simplify()` is perfectly fine for interactive use, where you just want SymPy
to do whatever it can to an expression. However, in programmatic usage, it's
better to avoid `simplify()` and use more [targeted simplification
functions](simplify-docs) instead (e.g., {func}`~.cancel`, {func}`~.expand`,
or {func}`~.collect`).

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

  The documentation for each function describes exactly what behavior it will
  have on the input expression.

- A targeted simplification will not do something unexpected if the expression
  contains an unexpected form, or an unexpected subexpression. This is
  especially the case if simplification functions are applied with `deep=False`
  to only apply the simplification to the top-level expression.

Some other simplification functions are heuristical in nature, and care should
be taken with them as well. For example, the {func}`~.trigsimp` function is a
heuristic targeted to trigonometric functions, but the routines in the
{mod}`sympy.simplify.fu` submodule allow applying specific trigonometric
identities.

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
### Don't Hardcode Symbol Names in Python Functions

Instead of hard-coding {class}`~.Symbol` names inside of a function
definition, make the symbols a parameter to the function.

For example, consider a function `theta_operator` that computes the [theta
operator](https://en.wikipedia.org/wiki/Theta_operator) $\theta =
zD_z$:

**Don't**

```py
def theta_operator(expr):
    z = symbols('z')
    return z*expr.diff(z)
```

**Do**

```py
def theta_operator(expr, z):
    return z*expr.diff(z)
```

A hard-coded symbol name has the disadvantage of requiring all expressions to
use that exact symbol name. In the above example, it is not possible to
compute $\theta = xD_x$ because it is hard-coded to $zD_z$. What's worse,
trying to do so silently leads to a wrong result instead of an error, since
`x` is treated as a constant expression:

```py
>>> def theta_operator(expr):
...     z = symbols('z')
...     return z*expr.diff(z)
>>> theta_operator(x**2) # The expected answer is 2*x**2
0
```

This is particularly problematic if the function accepts arbitrary user input,
as the user may be using a different variable name that makes more sense in
their mathematical context. And if the user already used the symbol `z` but as
a constant, they would need to swap things around with `subs` before being
able to use the function.

The other reason this antipattern is problematic is due to the gotcha that
symbols with assumptions are considered unequal to symbols without
assumptions. If someone defined their expression using

```py
>>> z = symbols('z', positive=True)
```

for example, to make further simplifications possible (see
[](best-practices-defining-symbols) above), the function hard-coding
`Symbol('z')` without assumptions would not work:

```py
>>> theta_operator(z**2)
0
```

By making the symbol an argument to the function, like `theta_operator(expr,
z)`, these problems all go away.

(best-practices-separate-symbolic-and-numeric-code)=
### Separate Symbolic and Numeric Code

SymPy sets itself apart from most of the rest of the libraries in the Python
ecosystem in that it operates symbolically, whereas other libraries, like
NumPy, operate numerically. These two paradigms are different enough that it's
always best to keep them as separate as possible.

Importantly, SymPy is not designed to work with NumPy arrays, and conversely,
NumPy will not work directly with SymPy objects.

```py
>>> import numpy as np
>>> import sympy
>>> a = np.array([0., 1., 2.])
>>> sympy.sin(a)
Traceback (most recent call last):
...
AttributeError: 'ImmutableDenseNDimArray' object has no attribute 'as_coefficient'
```

```py
>>> x = Symbol('x')
>>> np.sin(x) # NumPy functions do not know how to handle SymPy expressions
Traceback (most recent call last):
...
TypeError: loop of ufunc does not support argument 0 of type Symbol which has no callable sin method
```

If you want to use both SymPy and NumPy, you should explicitly convert your
SymPy expressions into NumPy functions using {func}`~.lambdify`. The typical
workflow in SymPy is to model your problem symbolically using SymPy, then
convert the result into a numerical function with `lambdify()` that can be
evaluated on NumPy arrays. For advanced use-cases, `lambdify()`/NumPy may not
be enough and you may instead need to use SymPy's more general [code
generation](codegen_prose) routines to generate code for other fast numerical
languages such as Fortran or C.

```python
>>> # First symbolically construct the expression you are interested in with SymPy
>>> from sympy import diff, sin, exp, lambdify, symbols
>>> x = symbols('x')
>>> expr = diff(sin(x)*exp(x**2), x)

>>> # Then convert it to a numeric function with lambdify()
>>> f = lambdify(x, expr)

>>> # Now use this function with NumPy
>>> import numpy as np
>>> a = np.linspace(0, 10)
>>> f(a) # doctest: +SKIP
[ 1.00000000e+00  1.10713341e+00  1.46699555e+00 ... -3.15033720e+44]
```

These are some antipatterns that should be generally avoided

- **Do not use `import math`.** It is virtually never necessary to use the
  [standard library `math`
  module](https://docs.python.org/3/library/math.html) alongside SymPy (or
  NumPy). Every function that is in `math` is already in SymPy. SymPy can
  compute values numerically using {term}`evalf`, which provides more
  precision and accuracy than `math`. Or better, SymPy will by default compute
  things symbolically. Functions and constants in `math` are floats, which are
  inexact. SymPy always works better with exact quantities when possible. For
  example,

  ```py
  >>> import math
  >>> math.pi # a float
  3.141592653589793
  >>> import sympy
  >>> sympy.sin(math.pi)
  1.22464679914735e-16
  ```

  The result of `sympy.sin(math.pi)` is not `0` as you might expect, because
  `math.pi` is only an approximation of $\pi$, equal to 16 digits. On the
  other hand, `sympy.pi` is *exactly* equal to $\pi$ because it is represented
  symbolically, so it is able to give the exact answer:

  ```
  >>> sympy.sin(sympy.pi)
  0
  ```

  So in general, one should [prefer symbolic
  representations](best-practices-exact-rational-numbers-vs-floats). But even
  if you actually do want a float, you are better off using SymPy's `evalf()`
  rather than `math`. This avoids the pitfall that `math` functions can only
  operate on `float` objects, not symbolic expressions

  ```py
  >>> x = Symbol('x')
  >>> math.sin(x)
  Traceback (most recent call last):
  ...
  TypeError: Cannot convert expression to float
  ```

  And furthermore, SymPy's `evalf()` is more accurate than `math`, because it
  uses arbitrary precision arithmetic, and allows you to specify any number
  of digits.

  ```py
  >>> sympy.sin(1).evalf(30)
  0.841470984807896506652502321630
  >>> math.sin(1)
  0.8414709848078965
  ```

  Even when using NumPy, `math` should be avoided. NumPy functions are faster
  than their `math` equivalents, support a larger range of numerical dtypes,
  and can operate on arrays of values, whereas `math` functions can only
  operate on a single scalar at a time.

- **Don't pass SymPy expressions to a NumPy function.** You should not pass a
  SymPy expression to a NumPy function. This includes anything in the `numpy`
  or `scipy` namespaces, as well as most functions from other Python libraries
  such as `matplotlib`. These functions are only designed to work with NumPy
  arrays with numeric values.

- **Don't pass SymPy expressions to a lambdified function.** Similar to the
  previous point, you should not pass SymPy expressions to a function created
  with `lambdify`. In effect, the functions returned by `lambdify` *are* NumPy
  functions, so the situation here is exactly the same. It is possible that in
  some cases a function created from `lambdify()` will work with a SymPy
  expression, but this is just an accident of the way it works. See [the "how
  it works" section of the `lambdify()` documentation](lambdify-how-it-works)
  for more details on why this happens.

- **Avoid storing SymPy expressions in a NumPy array.** While it is
  technically possible to store SymPy expressions inside of a NumPy array,
  doing so usually represents a mistake. A sign that this is happening is if
  the `dtype` of the NumPy array is `object` (instead of a numeric dtype like
  `float64` or `int64`).

  Just as one should avoid using NumPy when doing symbolic calculations with
  SymPy, one should stop using SymPy once the calculation have moved over to
  the numeric side of things with NumPy.

  A NumPy array that contains SymPy expressions effectively has the same
  problem as trying to call NumPy functions directly on a SymPy expression.
  They do not know how to operate on SymPy objects, so they will fail. This
  applies even if the SymPy objects are all SymPy
  [`Float`s](sympy.core.numbers.Float).

  ```
  >>> import numpy as np
  >>> import sympy
  >>> a = np.asarray([sympy.Float(1.0), sympy.Float(0.0)]) # Do not do this
  >>> print(repr(a)) # Note that the dtype is 'object'
  array([1.00000000000000, 0.0], dtype=object)
  >>> np.sin(a)
  Traceback (most recent call last):
  ...
  TypeError: loop of ufunc does not support argument 0 of type Float which has no callable sin method
  ```

  If you are doing this, you should probably either be using native NumPy
  floats, or, if you really do want to store an array of SymPy expressions,
  you should use SymPy's [`Matrix`](sympy.matrices.dense.Matrix) or
  `NDimArray` classes.

## Advanced Usage

### Be Careful Comparing and Sorting Symbolic Objects

Be careful with programmatic code that compares numerical quantities, either
directly using an inequality (`<`, `<=`, `>`, `>=`) or indirectly with
something like `sorted`. The issue is that if an inequality is unknown, the
result will be symbolic, like

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

-  **Disallow symbolic inputs.** If your function cannot possibly work on
  symbolic inputs, you can explicitly disallow them. The primary benefit here
  is to give a more readable error message to users than `TypeError:  cannot
  determine truth value of Relational`. The
  {attr}`~sympy.core.expr.Expr.is_number` attribute can be used to check if an
  expression can be evaluated to a specific number with `evalf()`. If you want
  to only accept integers, you can check `isinstance(x, Integer)` (after
  calling `sympify()` to convert Python ints). Beware that `is_integer` uses
  the assumptions system and may be True even for symbolic objects, like
  `Symbol('x', integer=True)`.

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
  return a result that represents both possibilities symbolically. This is
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
  >>> from sympy import Piecewise, pprint
  >>> expr = Piecewise((1, x > 0), (0, True))
  >>> pprint(expr, use_unicode=True)
  ⎧1  for x > 0
  ⎨
  ⎩0  otherwise
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

  Alternatively, try to write the function in a way so that the correctness of
  the result does not depend on the order that arguments are processed in.

## Custom SymPy Objects

SymPy is designed to be extended with custom classes, typically by subclassing
{term}`Basic`, {term}`Expr`, or {term}`Function <Function (class)>`. All the
symbolic classes in SymPy itself are written this way, and the points here
apply equally to them as to user-defined classes.

For an in-depth guide on how to write a `Function` subclass, see the [guide on
writing custom functions](custom-functions).

(best-practices-args-invariants)=
### Args Invariants

Custom SymPy objects should always satisfy the following invariants:

1. `all(isinstance(arg, Basic) for arg in args)`
2. `expr.func(*expr.args) == expr`

The first says that all elements of {term}`args` should be instances of
{term}`Basic`. The second says that an expression should be rebuildable from
its `args` (note that {term}`func` is usually the same as `type(expr)`, though
it may not always be).

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

1. If an expression had args that were not `Basic`, they would fail with
   `AttributeError` on a recursive call, because the non-`Basic` args would not
   have the `.args` or `.func` attributes.

2. If an expression did not rebuild from its `args`, the line `return
   exr.func(*newargs)` would fail, even in the trivial case where none of the
   args are changed by the replacement, which should effectively be a no-op.

Making all `args` instances of `Basic` usually just means calling `_sympify()`
on the inputs to the class so that they are basic instances. If you want to
store a string on a class, you should either use a `Symbol` or
`sympy.core.symbols.Str`.

In some cases a class may accept args in multiple equivalent forms. It is
important that whatever form is stored in `args` is one of the ways that can
be used to reconstruct the class. It is okay to normalize `args` as long as
that normalized form is accepted as input. For example, `Integral` always
stores the variable argument as a tuple to make things easier to process
internally, but this form is also accepted by the class constructor:

```py
>>> from sympy import Integral
>>> expr = Integral(sin(x), x)
>>> expr.args # args are normalized
(sin(x), (x,))
>>> Integral(sin(x), (x,)) # Also accepted
Integral(sin(x), x)
```

Note that most user-defined custom functions should be defined by subclassing
`Function` (see the [guide to writing custom functions](custom-functions)).
The `Function` class automatically takes care of both of the args invariants,
so if you are using it, you do not need to worry about this.

(best-practices-avoid-automatic-evaluation)=
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
objects, not just functions). But in a nutshell, the reason for this is that
the only way to prevent automatic evaluation is to use `evaluate=False`, which
is fragile. Additionally, code will invariably be written assuming the
invariants that are true due to automatic evaluations, meaning that
expressions created with `evaluate=False` can lead to wrong results from this
code. This also means that removing automatic evaluation later can be
difficult.

Evaluation that can potentially be expensive (for instance, applying a
symbolic identity) is itself bad because it can make creating an expression
without even doing anything with it allow. This also applies to checking for
symbolic assumptions (like `x.is_integer`), so this should also be avoided in
class constructors.

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

Functions and classes that accept an arbitrary number of arguments should
either accept the arguments directly, like `f(*args)`, or as a single
argument, like `f(args)`. They should not try to support both at once.

The reason is that this makes it impossible to represented nested collections.
For example, take the {class}`~.FiniteSet` class. It is constructed like
`FiniteSet(x, y, z)` (i.e., using `*args`).

```py
>>> from sympy import FiniteSet
>>> FiniteSet(1, 2, 3)
{1, 2, 3}
```

It might be tempting to also support `FiniteSet([1, 2, 3])`, to match the
built-in `set`. However, doing so would make it impossible to represent a
nested `FiniteSet` containing a single `FiniteSet`, like $\{\{1, 2, 3\}\}$:

```py
>>> FiniteSet(FiniteSet(1, 2, 3)) # We don't want this to be the same as {1, 2, 3}
FiniteSet({1, 2, 3})
```

As to whether `args` or `*args` should be used, if it is only possible for
there to be a finite number of arguments, `*args` is generally better, as this
makes things easier to deal with using the object's {term}`args`, since
`obj.args` will be the direct arguments of the class. However, if it is
possible that you might want to support a symbolic infinite collection in
addition to finite ones, like {class}`~.Integers` or {class}`~.Range`, then it
is better to use `args` as this will be impossible to do with `*args`.

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
  stored in `args`. Note that an object's `args` should be usable to [recreate
  the object](best-practices-args-invariants) (e.g., something like
  `YourObject(*instance.args)` should recreate `instance`).

  Additionally, it should be mentioned that it is not a good idea to subclass
  `Symbol` if you plan to store anything extra in `args`. `Symbol` is designed
  around having no `args`. You are better off subclassing `Function` (see
  [](custom-functions)) or `Expr` directly. If you simply want to have two
  symbols that are distinct from one another, the best approach is often just
  to give them different names. If you are concerned about how they are
  printed, you can replace them with a more canonical name when it comes time
  to print things, or use a [custom printer](module-printing).

- **Store the data about the object separately.** This is the best approach if
  the extra data is not directly related to an object's mathematical
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
  for example, because of the limitation that [`args` should only contain
  `Basic` instances](best-practices-args-invariants). It is still possible to
  create custom SymPy objects in these situations by using a custom
  {term}`func` that is different from `type(expr)` (in this case, you would
  override `__eq__` on the `func` [rather than on the
  class](best-practices-eq)).

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
  returns a boolean and is inexpensive to compute. For example, if `b` is a
  builtin Python container like `list`, `dict`, or `set`, then `a in b` uses
  `==`.[^dict-footnote]

[^dict-footnote]: Python dicts and sets use `hash`, but fallback to using `==`
    when there is a hash collision.

- SymPy internally uses `==` all over the place, both explicitly and
  implicitly via things like `in` or dictionary keys. This
  usage all implicitly assumes that `==` operates structurally.


In affect, *structural equality* means that if `a == b` is `True`, then `a`
and `b` are for all intents and purposes the same object. This is because all
SymPy objects are {term}`immutable`. When `a == `, any SymPy function may
freely replace `a` with `b` in any subexpression.

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
  method to implement the smarter equality checking (for example, the `equals`
  method).

  Another option is to define a {term}`canonicalization <canonicalize>` method
  that puts objects into canonical form (e.g., via `doit`), so that, for
  instance, `x.doit() == y.doit()` is true whenever `x` and `y` are
  mathematically equal. This is not always possible because not every type of
  object has a computable canonical form, but it is a convenient approach when
  one does exist.

- To make `==` check for some additional attributes beyond those stored in the
  `args` of an expression. See the [](best-practices-extra-attributes) section
  above for more details on why it's a bad idea to directly store extra
  attributes on a SymPy object, and what the best alternatives are.

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
creating a new expression that itself could lead to many recursive
assumptions queries is bad for performance compared to querying the desired
property more directly.

This generally means using methods like {meth}`~.as_independent` or
{meth}`~.as_coeff_mul` and checking
the `args` of expressions directly (see the [custom functions
guide](custom-functions-assumptions) for an example).

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
