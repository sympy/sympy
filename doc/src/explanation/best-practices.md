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

### Avoid `simplify()`

(best-practices-dont-hardcode-symbol-names)=
### Don't Hardcode Symbol Names in Functions

### Separate Symbolic and Non-symbolic Code

## Custom SymPy Objects

### Args Invariants

### Avoid Too Much Automatic Evaluation

### Don't Denest Collections

### Don't Store Attributes in `.args`

### Don't Overwrite `__eq__`
