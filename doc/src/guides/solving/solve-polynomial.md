# Solve (Find the Roots of) a Polynomial Algebraically 

Use SymPy to solve (find the roots of) a polynomial algebraically. For example,
solving $ax^2 + bx + c = 0$ for $x$ yields $x = \frac{-b\pm\sqrt{b^2 -
4ac}}{2a}$.

Alternatives to consider:
- If you need a numeric (rather than algebraic) solution, you can use either
    - NumPy's {external:func}`~numpy.roots`
    - SciPy's {external:func}`~scipy.optimize.root`

Here is an example of solving a polynomial algebraically:

```py
>>> from sympy import solve
>>> from sympy.abc import x, a, b, c
>>> solve(a*x**2 + b*x + c, x, dict=True)
[{x: (-b - sqrt(-4*a*c + b**2))/(2*a)}, {x: (-b + sqrt(-4*a*c + b**2))/(2*a)}]
```

This example reproduces the [quadratic
formula](https://en.wikipedia.org/wiki/Quadratic_formula).

## Guidance

### *Guidance 1*

*Guidance 1 content*

### *Guidance 2*

*Guidance 2 content*


## Solve (Find the Roots of) a Polynomial

You can solve a polynomial algebraically in several ways. The one to use depends
on whether you
- want an algebraic or numeric answer
- want the multiplicity of each root (how many times each root is a solution).
  In the `expression` below representing $(x+2)^2(x-3)$, the root -2 has a
  multiplicity of two because $x+2$ is squared, whereas 3 has a multiplicity of
  one because $x-3$ has no exponent. Similarly, for the `symbolic` expression,
  the root $-a$ has a multiplicity of two and the root $b$ has a multiplicity of
  one.

```py
>>> from sympy import solve, roots, real_roots, factor, nroots, RootOf, expand
>>> from sympy import Poly
>>> from sympy.abc import x, a, b
>>> expression = (x+2)**2 * (x-3)
>>> symbolic = (x+a)**2 * (x-b)
```

### Algebraic Solution Without Root Multiplicities

You can use SymPy's standard {func}`~.solve` function, though it will not return
the multiplicity of roots:

```py
>>> solve(expression, x, dict=True)
[{x: -2}, {x: 3}]
>>> solve(symbolic, x, dict=True)
[{x: -a}, {x: b}]
```

Refer to [](solve-equation-algebraically.md) for more about using
{func}`~.solve`.

### Algebraic Solution With Root Multiplicities

#### {func}`~.roots`

{func}`~.roots` is the most rigorous function because it can give explicit
expressions for the roots of polynomials that have symbolic coefficients (that
is, if there are symbols in the coefficients) if {func}`~.factor` does not
reveal them. However, it may fail for some polynomials. Here are examples of
{func}`~.roots`:

```py
>>> roots(expression, x)
{-2: 2, 3: 1}
>>> roots(symbolic, x)
{-a: 2, b: 1}
```

It returns results as a dictionary, where the key is the root (for example, -2)
and the value is the multiplicity of that root (for example, 2).

{func}`~.roots` function uses a combination of techniques (factorization,
decomposition, radical formulae) to find expressions in radicals if possible for
the roots. When it can find some radical expressions for the roots, it returns
them along with their multiplicity. This function will fail for most high-degree
polynomials (five or greater) because they do not have radical solutions, and
there is no guarantee that they have closed-form solutions at all, as explained
by the [Abel-Ruffini
theorem](https://en.wikipedia.org/wiki/Abel%E2%80%93Ruffini_theorem).

#### Factor the Equation

A different approach is to factor a polynomial using {func}`~.factor`, which
does not give the roots directly but can give you simpler expressions:

```py
>>> expression_expanded = expand(expression)
>>> expression_expanded
x**3 + x**2 - 8*x - 12
>>> factor(expression_expanded)
(x - 3)*(x + 2)**2
>>> symbolic_expanded = expand(symbolic)
>>> symbolic_expanded
-a**2*b + a**2*x - 2*a*b*x + 2*a*x**2 - b*x**2 + x**3
>>> factor(symbolic_expanded)
(a + x)**2*(-b + x)
```

{func}`~.factor` can also factorize a polynomial in a given [polynomial
ring](polys-ring) which can reveal roots lie in the coefficient ring. For
example, if the polynomial has rational coefficients, then {func}`~.factor` will
reveal any rational roots. If the coefficients are polynomials involving, for
example, symbol $a$ with rational coefficients then any roots that are
polynomial functions of $a$ with rational coefficients will be revealed. In this
example, {func}`~.factor` reveals that $x = a^2$ and $x = -a^3 - a$ are roots:

```py
>>> from sympy import expand, factor
>>> from sympy.abc import x, a
>>> p = expand((x - a**2)*(x + a + a**3))
>>> p
-a**5 + a**3*x - a**3 - a**2*x + a*x + x**2
>>> factor(p)
(-a**2 + x)*(a**3 + a + x)
```

### Exact Numeric Solution With Root Multiplicities

#### `real_roots`

If the roots to your polynomial are real, using {func}`~.real_roots` ensures
that only real (not complex or imaginary) roots will be returned.

```py
>>> real_roots(expression)
[-2, -2, 3]
```

{func}`~.real_roots` calls {func}`~sympy.polys.rootoftools.RootOf`, so you can
get the same results by iterating over the number of roots of your equation:

```py
>>> [RootOf(expression, n) for n in range(0,3)]
[-2, -2, 3]
```

### Approximate Numeric Solution With Root Multiplicities

#### `nroots`

{func}`~.nroots` gives an approximate numerical approximation to the roots of a
polynomial. This example demonstrates that it can include numerical noise, for
example a (negligible) imaginary component in what should be a real root:

```py
>>> nroots(expression)
[3.0, -2.0 - 4.18482169793536e-14*I, -2.0 + 4.55872552179222e-14*I]
```

{func}`~.nroots`  is analogous to NumPy's {external:func}`~numpy.roots`
function. Usually the difference between these two is that {func}`~.nroots` is
more accurate but slower.

A major advantage of {func}`~.nroots` is that it can compute numerical
approximations of the roots of any polynomial whose coefficients can be
numerically evaluated with {func}`~sympy.core.evalf` (that is, they do not have
free symbols). Contrarily, symbolic solutions may not be possible for
higher-order (fifth or greater) polynomials as explained by the [Abel-Ruffini
theorem](https://en.wikipedia.org/wiki/Abel%E2%80%93Ruffini_theorem). Even if
closed-form solutions are available, they may have so many terms that they are
not useful in practice. You may therefore want to use {func}`~.nroots` to find
approximate numeric solutions even if closed-form symbolic solutions are
available.

{func}`~.nroots` can fail sometimes for polynomials that are numerically ill
conditioned, for example [Wilkinson's
polynomial](https://en.wikipedia.org/wiki/Wilkinson%27s_polynomial).

## Complex Roots

For complex roots, similar functions can be used, for example {func}`~.solve`:

```py
>>> from sympy import solve, roots, nroots, real_roots, expand, RootOf, CRootOf, Symbol
>>> from sympy import Poly
>>> from sympy.abc import x
>>> expression_complex = (x**2+4)**2 * (x-3)
>>> solve(expression_complex, x, dict=True)
[{x: 3}, {x: -2*I}, {x: 2*I}]
```

If the constants are symbolic, you may need to specify their domain for SymPy to
recognize that the solutions are not real. For example, specifying that $a$ is
positive leads to imaginary roots:

```py
>>> a = Symbol("a", positive=True)
>>> symbolic_complex = (x**2+a)**2 * (x-3)
>>> solve(symbolic_complex, x, dict=True)
[{x: 3}, {x: -I*sqrt(a)}, {x: I*sqrt(a)}]
```

{func}`~.roots` will also find imaginary or complex roots:

```py
>>> roots(expression_complex, x)
{3: 1, -2*I: 2, 2*I: 2}
```

{func}`~sympy.polys.rootoftools.RootOf` will also return complex roots:

```py
>>> [RootOf(expression_complex, n) for n in range(0,3)]
[3, -2*I, -2*I]
```

{func}`~.real_roots` will only return the real roots and give no indication that
there are complex roots, so use it with caution if your equation could have
complex roots:

```py
>>> real_roots(expression_complex)
[3]
```

If you make the expression into a polynomial class {class}`~.Poly`, you can use
its {meth}`~sympy.polys.polytools.Poly.all_roots` method to find the roots:

```py
>>> expression_complex_poly = Poly(expression_complex)
>>> expression_complex_poly.all_roots()
[3, -2*I, -2*I, 2*I, 2*I]
```

## Use the Solution Result

### *Usage Method 1*

*Usage method 1 content*

### *Usage Method 2*

*Usage method 2 content*

## *Tradeoffs (speed vs. accuracy, etc.) for function*

### *Tradeoff 1*

*Tradeoff 1 content*

### *Tradeoff 2*

*Tradeoff 2 content*

## Not All Equations Can Be Solved

### Equations With No Solution

*Equations with no solution content*

### Equations With No Closed-Form Solution

As mentioned above, higher-order polynomials (fifth or greater) are unlikely to
have closed-form solutions, so you may have to use a numerical method such as
[`nroots` as described above](#nroots).

Please post the problem on the [mailing
list](https://groups.google.com/g/sympy), or open an issue on [SymPy's GitHub
page](https://github.com/sympy/sympy/issues). Until the issue is resolved, you
can *workaround*.
