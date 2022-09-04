# Solve (Find the Roots of) a Polynomial Algebraically 

Use SymPy to solve (find the roots of) a polynomial algebraically. For example,
solving $ax^2 + bx + c = 0$ for $x$ yields $x = \frac{-b\pm\sqrt{b^2 -
4ac}}{2a}$.

Alternatives to consider:
- *alternative 1*
- *alternative 2*

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


## Solve (Find the Roots of) a Polynomial Algebraically

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

#### {func}`~.factor`

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

### Exact Numeric Solution Without Root Multiplicities

### Exact Numeric Solution With Root Multiplicities

#### {func}`~.real_roots`

If the roots to your polynomial are real, using {func}`~.real_roots` ensures
that only real (not complex) roots will be returned.

```py
>>> real_roots(expression)
[-2, -2, 3]
```

### Approximate Numeric Solution With Root Multiplicities

#### {func}`~.nroots`

{func}`~.nroots` gives an approximate numerical approximation to the roots of a
polynomial. This example demonstrates that it can include numerical noise, for
example a (negligible) imaginary component in what should be a real root:

```py
>>> nroots(expression)
[3.0, -2.0 - 4.18482169793536e-14*I, -2.0 + 4.55872552179222e-14*I]
```

nroots(expression)
# [3.00000000000000,
#  -2.0 - 4.18482169793536e-14*I,
#  -2.0 + 4.55872552179222e-14*I]
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

### Equations With No Analytical Solution

*Equations with no analytical solution content*

Please post the problem on the [mailing
list](https://groups.google.com/g/sympy), or open an issue on [SymPy's GitHub
page](https://github.com/sympy/sympy/issues). Until the issue is resolved, you
can *workaround*.
