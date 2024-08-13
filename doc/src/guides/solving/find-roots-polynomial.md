# Find the Roots of a Polynomial Algebraically or Numerically

Use SymPy to find the roots of a univariate polynomial algebraically. For
example, finding the roots of $ax^2 + bx + c$ for $x$ yields $x =
\frac{-b\pm\sqrt{b^2 - 4ac}}{2a}$.

## Alternatives to Consider

- If you need a numeric (rather than algebraic) solution, you can use either
    - NumPy's {external:func}`~numpy.roots`
    - SciPy's {external:func}`~scipy.optimize.root`
- If you need to solve systems of polynomial equations algebraically, use
  {func}`~.solve`

## Example of Finding the Roots of a Polynomial Algebraically

Here is an example of finding the roots of a polynomial algebraically:

```py
>>> from sympy import roots
>>> from sympy.abc import x, a, b, c
>>> roots(a*x**2 + b*x + c, x)
{-b/(2*a) - sqrt(-4*a*c + b**2)/(2*a): 1,
 -b/(2*a) + sqrt(-4*a*c + b**2)/(2*a): 1}
```

This example reproduces the [quadratic
formula](https://en.wikipedia.org/wiki/Quadratic_formula).

## Functions to Find the Roots of a Polynomial

There are several functions that you can use to find the roots of a polynomial:
- {func}`~.solve` is a general solving function which can find roots, though is
  less efficient than {meth}`~sympy.polys.polytools.Poly.all_roots` and is the
  only function in this list that does not convey the multiplicity of roots;
  {func}`~.solve` also works on [non-polynomial
  equations](solve-equation-algebraically.md) and [systems of non-polynomial
  equations](solve-system-of-equations-algebraically.md)
- {func}`~.roots` computes the symbolic roots of a univariate polynomial; will
fail for most high-degree polynomials (five or greater)
- {func}`~.nroots` computes numerical approximations of the roots of any
polynomial whose coefficients can be numerically evaluated, whether the
coefficients are rational or irrational
- {func}`~sympy.polys.rootoftools.RootOf` can represent all the roots exactly of
a polynomial of arbitrarily large degree, as long as the coefficients are
rational numbers. {func}`~sympy.polys.rootoftools.RootOf` can avoid both
ill-conditioning and returning spurious complex parts because it uses a more
exact, but much slower, numerical algorithm based on isolating intervals. The
  following two functions use {func}`~sympy.polys.rootoftools.RootOf` so they
  have the same properties:
    - {func}`~.real_roots` can find all the real roots exactly of a polynomial
of arbitrarily large degree; because it finds only the real roots, it can be
  more efficient than functions that find all roots.
    - {meth}`~sympy.polys.polytools.Poly.all_roots` can find all the roots
exactly of a polynomial of arbitrarily large degree
- {func}`~.factor` factors a polynomial into irreducibles and can reveal that
  roots lie in the coefficient ring

Each will be used on this page.

## Guidance

Refer to
[](solving-guidance.md#include-the-variable-to-be-solved-for-in-the-function-call)
and [](solving-guidance.md#use-exact-values).

## Find the Roots of a Polynomial

You can find the roots of a polynomial algebraically in several ways. The one to
use depends on whether you
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

{func}`~.solve` will first try using {func}`~.roots`; if that doesn't work, it
will try using {meth}`~sympy.polys.polytools.Poly.all_roots`. For cubics
(third-degree polynomials) and quartics (fourth-degree polynomials), that means
that {func}`~.solve` will use radical formulae from roots rather than
{func}`~sympy.polys.rootoftools.RootOf` even if RootOf is possible. The cubic
and quartic formulae often give very complex expressions that are not useful in
practice. As a result, you may want to set the {func}`~.solve` parameter
`cubics` or `quartics` to `False` to return
{func}`~sympy.polys.rootoftools.RootOf` results:

```py
>>> from sympy import solve
>>> from sympy.abc import x
>>> # By default, solve() uses the radical formula, yielding very complex terms
>>> solve(x**4 - x + 1, x)
[-sqrt(2/(3*(1/16 + sqrt(687)*I/144)**(1/3)) + 2*(1/16 + sqrt(687)*I/144)**(1/3))/2 - sqrt(-2*(1/16 + sqrt(687)*I/144)**(1/3) - 2/sqrt(2/(3*(1/16 + sqrt(687)*I/144)**(1/3)) + 2*(1/16 + sqrt(687)*I/144)**(1/3)) - 2/(3*(1/16 + sqrt(687)*I/144)**(1/3)))/2,
 sqrt(2/(3*(1/16 + sqrt(687)*I/144)**(1/3)) + 2*(1/16 + sqrt(687)*I/144)**(1/3))/2 - sqrt(-2*(1/16 + sqrt(687)*I/144)**(1/3) + 2/sqrt(2/(3*(1/16 + sqrt(687)*I/144)**(1/3)) + 2*(1/16 + sqrt(687)*I/144)**(1/3)) - 2/(3*(1/16 + sqrt(687)*I/144)**(1/3)))/2,
 sqrt(-2*(1/16 + sqrt(687)*I/144)**(1/3) - 2/sqrt(2/(3*(1/16 + sqrt(687)*I/144)**(1/3)) + 2*(1/16 + sqrt(687)*I/144)**(1/3)) - 2/(3*(1/16 + sqrt(687)*I/144)**(1/3)))/2 - sqrt(2/(3*(1/16 + sqrt(687)*I/144)**(1/3)) + 2*(1/16 + sqrt(687)*I/144)**(1/3))/2,
 sqrt(-2*(1/16 + sqrt(687)*I/144)**(1/3) + 2/sqrt(2/(3*(1/16 + sqrt(687)*I/144)**(1/3)) + 2*(1/16 + sqrt(687)*I/144)**(1/3)) - 2/(3*(1/16 + sqrt(687)*I/144)**(1/3)))/2 + sqrt(2/(3*(1/16 + sqrt(687)*I/144)**(1/3)) + 2*(1/16 + sqrt(687)*I/144)**(1/3))/2]
>>> # If you set quartics=False, solve() uses RootOf()
>>> solve(x**4 - x + 1, x, quartics=False)
[CRootOf(x**4 - x + 1, 0),
 CRootOf(x**4 - x + 1, 1),
 CRootOf(x**4 - x + 1, 2),
 CRootOf(x**4 - x + 1, 3)]
```

Writing the first root from {func}`~.solve` in standard mathematical notation
emphasizes how complex it is:

$$- \frac{\sqrt{\frac{2}{3 \sqrt[3]{\frac{1}{16} + \frac{\sqrt{687} i}{144}}} +
2 \sqrt[3]{\frac{1}{16} + \frac{\sqrt{687} i}{144}}}}{2} - \frac{\sqrt{- 2
\sqrt[3]{\frac{1}{16} + \frac{\sqrt{687} i}{144}} - \frac{2}{\sqrt{\frac{2}{3
\sqrt[3]{\frac{1}{16} + \frac{\sqrt{687} i}{144}}} + 2 \sqrt[3]{\frac{1}{16} +
\frac{\sqrt{687} i}{144}}}} - \frac{2}{3 \sqrt[3]{\frac{1}{16} +
\frac{\sqrt{687} i}{144}}}}}{2}$$

Further, there is no general radical formula for quintics (fifth degree) or
higher polynomials, so their {func}`~sympy.polys.rootoftools.RootOf`
representations may be the best option.

Refer to [](solve-equation-algebraically.md) for more about using
{func}`~.solve`.

### Algebraic Solution With Root Multiplicities

#### `roots`

{func}`~.roots` can give explicit expressions for the roots of polynomials that
have symbolic coefficients (that is, if there are symbols in the coefficients)
if {func}`~.factor` does not reveal them. However, it may fail for some
polynomials. Here are examples of {func}`~.roots`:

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
>>> from sympy import real_roots
>>> from sympy.abc import x
>>> cubed = x**3 - 1
>>> # roots() returns real and complex roots
>>> roots(cubed)
{1: 1, -1/2 - sqrt(3)*I/2: 1, -1/2 + sqrt(3)*I/2: 1}
>>> # real_roots() returns only real roots
>>> real_roots(cubed)
[1]
```

{func}`~.real_roots` calls {func}`~sympy.polys.rootoftools.RootOf`, so for
equations whose roots are all real, you can get the same results by iterating
over the number of roots of your equation:

```py
>>> [RootOf(expression, n) for n in range(3)]
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

If you want numeric approximations of the real roots, but you want to know
exactly which roots are real, then the best method is {func}`~.real_roots` with
{func}`~sympy.core.evalf`:

```py
>>> [r.n(2) for r in real_roots(expression)]
[-2.0, -2.0, 3.0]
>>> [r.is_real for r in real_roots(expression)]
[True, True, True]
```

{func}`~.nroots` is analogous to NumPy's {external:func}`~numpy.roots` function.
Usually the difference between these two is that {func}`~.nroots` is more
accurate but slower.

A major advantage of {func}`~.nroots` is that it can compute numerical
approximations of the roots of any polynomial whose coefficients can be
numerically evaluated with {func}`~sympy.core.evalf` (that is, they do not have
free symbols). Contrarily, symbolic solutions may not be possible for
higher-order (fifth or greater) polynomials as explained by the [Abel-Ruffini
theorem](https://en.wikipedia.org/wiki/Abel%E2%80%93Ruffini_theorem). Even if
closed-form solutions are available, they may have so many terms that they are
not useful in practice. You may therefore want to use {func}`~.nroots` to find
approximate numeric solutions even if closed-form symbolic solutions are
available. For example, the closed-form roots of a fourth-order (quartic)
polynomial may be rather complicated:

```py
>>> rq0, rq1, rq2, rq3 = roots(x**4 + 3*x**2 + 2*x + 1)
>>> rq0
sqrt(-4 - 2*(-1/8 + sqrt(237)*I/36)**(1/3) + 4/sqrt(-2 + 7/(6*(-1/8 + sqrt(237)*I/36)**(1/3)) + 2*(-1/8 + sqrt(237)*I/36)**(1/3)) - 7/(6*(-1/8 + sqrt(237)*I/36)**(1/3)))/2 - sqrt(-2 + 7/(6*(-1/8 + sqrt(237)*I/36)**(1/3)) + 2*(-1/8 + sqrt(237)*I/36)**(1/3))/2
```

so you may prefer an approximate numerical solution:

```py
>>> rq0.n()
-0.349745826211722 - 0.438990337475312*I
```

{func}`~.nroots` can fail sometimes for polynomials that are numerically ill
conditioned, for example [Wilkinson's
polynomial](https://en.wikipedia.org/wiki/Wilkinson%27s_polynomial). Using
{func}`~sympy.polys.rootoftools.RootOf` and {func}`~sympy.core.evalf` as
described in [](#numerically-evaluate-crootof-roots) can avoid both
ill-conditioning and returning spurious complex parts because it uses a more
exact, but much slower, numerical algorithm based on isolating intervals.

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

{func}`~.real_roots` will return only the real roots.

```py
>>> real_roots(expression_complex)
[3]
```

An advantage of {func}`~.real_roots` is that it can be more efficient than
generating all the roots: {func}`~sympy.polys.rootoftools.RootOf` can be slow
for complex roots.

If you make the expression into a polynomial class {class}`~.Poly`, you can use
its {meth}`~sympy.polys.polytools.Poly.all_roots` method to find the roots:

```py
>>> expression_complex_poly = Poly(expression_complex)
>>> expression_complex_poly.all_roots()
[3, -2*I, -2*I, 2*I, 2*I]
```

## Use the Solution Result

The way to extract solutions from the result depends on the form of the result.

### List (`all_roots`, `real_roots`, `nroots`)

You can use standard Python list traversal techniques such as looping. Here, we
substitute each root into the expression to verify that the result is $0$:

```py
>>> expression = (x+2)**2 * (x-3)
>>> my_real_roots = real_roots(expression)
>>> my_real_roots
[-2, -2, 3]
>>> for root in my_real_roots:
...         print(f"expression({root}) = {expression.subs(x, root)}")
expression(-2) = 0
expression(-2) = 0
expression(3) = 0
```

### List of dictionaries (`solve`)

Refer to [](solve-equation-algebraically.md#use-the-solution-result).

### Dictionary (`roots`)

You can use standard Python list traversal techniques such as looping through
the keys and values in a dictionary. Here we print the value and multiplicity of
each root:

```py
>>> my_roots = roots(expression)
>>> my_roots
{-2: 2, 3: 1}
>>> for root, multiplicity in my_roots.items():
...     print(f"Root {root} has multiplicity of {multiplicity}")
Root 3 has multiplicity of 1
Root -2 has multiplicity of 2
```

### Expression (`factor`)

You can manipulate an algebraic expression using various SymPy techniques, for
example substituting in a symbolic or numeric value for $x$:

```py
>>> from sympy.abc import y
>>> factored = factor(expression_expanded)
>>> factored
(x - 3)*(x + 2)**2
>>> factored.subs(x, 2*y)
(2*y - 3)*(2*y + 2)**2
>>> factored.subs(x, 7)
324
```

## Tradeoffs

### Mathematical Exactness, Completeness of List of Roots, and Speed

Consider the high-order polynomial $x^5 - x + 1 = 0$. {func}`~.nroots` returns
numerical approximations to all five roots:

```py
>>> from sympy import roots, solve, real_roots, nroots
>>> from sympy.abc import x
>>> fifth_order = x**5 - x + 1
>>> nroots(fifth_order)
[-1.16730397826142,
 -0.181232444469875 - 1.08395410131771*I,
 -0.181232444469875 + 1.08395410131771*I,
 0.764884433600585 - 0.352471546031726*I,
 0.764884433600585 + 0.352471546031726*I]
```

{func}`~.roots` can sometimes return only a subset of the roots or nothing if it
can't express any roots in radicals. In this case, it returns no roots (an empty
set):

```py
>>> roots(fifth_order, x)
{}
```

But if you set the flag `strict=True`, {func}`~.roots` will inform you that all
roots cannot be returned:

```py
>>> roots(x**5 - x + 1, x, strict=True)
Traceback (most recent call last):
...
sympy.polys.polyerrors.UnsolvableFactorError: Strict mode: some factors cannot be solved in radicals, so a complete
list of solutions cannot be returned. Call roots with strict=False to
get solutions expressible in radicals (if there are any).
```

#### Get All Roots, Perhaps Implicitly
{func}`~.solve` will return all five roots as `CRootOf`
({func}`~sympy.polys.rootoftools.ComplexRootOf`) class members:

```py
>>> fifth_order_solved = solve(fifth_order, x, dict=True)
>>> fifth_order_solved
[{x: CRootOf(x**5 - x + 1, 0)},
{x: CRootOf(x**5 - x + 1, 1)},
{x: CRootOf(x**5 - x + 1, 2)},
{x: CRootOf(x**5 - x + 1, 3)},
{x: CRootOf(x**5 - x + 1, 4)}]
```

where the second argument in each `CRootOf` is the index of the root.

#### Numerically Evaluate `CRootOf` Roots
You can then numerically evaluate those `CRootOf` roots using `n` from
{func}`~sympy.core.evalf`:

```py
>>> for root in fifth_order_solved:
...     print(root[x].n(10))
-1.167303978
-0.1812324445 - 1.083954101*I
-0.1812324445 + 1.083954101*I
0.7648844336 - 0.352471546*I
0.7648844336 + 0.352471546*I
```

If you are only interested in the sole real root, it is faster to use
{func}`~.real_roots` because it will not attempt to find the complex roots:

```py
>>> real_root = real_roots(fifth_order, x)
>>> real_root
[CRootOf(x**5 - x + 1, 0)]
>>> real_root[0].n(10)
-1.167303978
```

### Representing Roots

{func}`~sympy.polys.rootoftools.RootOf`, {func}`~.real_roots`, and
{meth}`~sympy.polys.polytools.Poly.all_roots` can find all the roots exactly of
a polynomial of arbitrarily large degree despite the [Abel-Ruffini
theorem](https://en.wikipedia.org/wiki/Abel%E2%80%93Ruffini_theorem). Those
functions allow the roots to be categorized precisely and manipulated
symbolically.

```py
>>> from sympy import init_printing
>>> init_printing()
>>> real_roots(fifth_order)
        / 5           \
[CRootOf\x  - x + 1, 0/]
>>> r = r0, r1, r2, r3, r4 = Poly(fifth_order, x).all_roots(); r
        / 5           \         / 5           \         / 5           \         / 5           \         / 5           \
[CRootOf\x  - x + 1, 0/, CRootOf\x  - x + 1, 1/, CRootOf\x  - x + 1, 2/, CRootOf\x  - x + 1, 3/, CRootOf\x  - x + 1, 4/]
>>> r0
       / 5           \
CRootOf\x  - x + 1, 0/
```

Now that the roots have been found exactly, their properties can be determined
free of numerical noise. For example, we can tell whether roots are real or not.
If we request the {meth}`~sympy.core.expr.Expr.conjugate` (same real part and
imaginary part with opposite sign) of a root, for example `r1`, and that is
exactly equal to another root `r2`, that root `r2` will be returned:

```py
>>> r0.n()
-1.16730397826142
>>> r0.is_real
True
>>> r1.n()
-0.181232444469875 - 1.08395410131771*I
>>> r2.n()
-0.181232444469875 + 1.08395410131771*I
>>> r1
        / 5           \
CRootOf\x  - x + 1, 1/
>>> r1.conjugate()
        / 5           \
CRootOf\x  - x + 1, 2/
>>> r1.is_real
False
```

{func}`~.solve` will also give the complex roots where possible but it is less
efficient than using {meth}`~sympy.polys.polytools.Poly.all_roots` directly.

{func}`~sympy.polys.rootoftools.RootOf` exactly represents the root in a way
that can be manipulated symbolically, and computed to arbitrary precision. The
{func}`~sympy.polys.rootoftools.RootOf` representation makes it possible to
precisely:

- Compute all roots of a polynomial with exact rational coefficients.
- Decide exactly the multiplicity of every root.
- Determine exactly whether roots are real or not.
- Order the real and complex roots precisely.
- Know which roots are complex conjugate pairs of each other.
- Determine precisely which roots are rational vs irrational.
- Represent every possible algebraic number exactly.

The other numerical methods such NumPy's {external:func}`~numpy.roots`,
{func}`~.nroots`, and {func}`~.nsolve` cannot do any of these things robustly,
if at all. Similarly, when numerically evaluated using
{func}`~sympy.core.evalf`, the radical expressions returned by {func}`~.solve`
or {func}`~.roots` cannot do these things robustly.

## Not All Equations Can Be Solved

### Equations With No Closed-Form Solution

As mentioned above, higher-order polynomials (fifth or greater) are unlikely to
have closed-form solutions, so you may have to represent them using, for
example, [`RootOf` as described above](#representing-roots), or use a numerical
method such as [`nroots` as described above](#nroots).

## Report a Bug

If you encounter a bug with these commands, please post the problem on the
[SymPy mailing list](https://groups.google.com/g/sympy). Until the issue is
resolved, you can use another of the
[](#functions-to-find-the-roots-of-a-polynomial) or try one of the
[](#alternatives-to-consider).
