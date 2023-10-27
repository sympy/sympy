# Solve One or a System of Equations Numerically

Use SymPy to numerically solve a system of one or more equations. For example,
numerically solving $\cos(x) = x $ returns $ x \approx 0.739085133215161$.

Solving numerically is useful if:
- You only need a numeric solution, not a symbolic one
- A closed-form solution is not available or is overly complicated; refer to
  [](solving-guidance.md#when-you-might-prefer-a-numeric-solution)

{func}`~.solve` and {func}`~.solveset` will not try to find a numeric solution,
only a mathematically-exact symbolic solution. So if you want a numeric
solution, use {func}`~.nsolve`.

SymPy is designed for symbolic mathematics. If you do not need to do symbolic
operations, then for numerical operations you can use another free and
open-source package such as NumPy or SciPy which will be faster, work with
arrays, and have more algorithms implemented. The main reasons to use SymPy (or
its dependency [mpmath](https://mpmath.org/)) for numerical calculations are:
- to do a simple numerical calculation within the context of a symbolic
  calculation using SymPy
- if you need the arbitrary precision capabilities to get more digits of
  precision than you would get from float64.

## Alternatives to Consider

- SciPy's {external:func}`scipy.optimize.fsolve` can solve a system of
  (non-linear) equations
- NumPy's {external:func}`numpy.linalg.solve` can solve a system of linear
  scalar equations
- mpmath's {external:func}`~mpmath.findroot`, which {func}`~.nsolve` calls and
  can pass parameters to


## Example of Numerically Solving an Equation

Here is an example of numerically solving one equation:

```py
>>> from sympy import cos, nsolve, Symbol
>>> x = Symbol('x')
>>> nsolve(cos(x) - x, x, 1)
0.739085133215161
```

## Guidance

Overdetermined systems of equations are supported.

### Find Complex Roots of a Real Function

To solve for complex roots of real functions, specify a nonreal (either purely
imaginary, or complex) initial point:

```py
>>> from sympy import nsolve
>>> from sympy.abc import x
>>> nsolve(x**2 + 2, 1) # Real initial point returns no root
Traceback (most recent call last):
    ...
ValueError: Could not find root within given tolerance. (4.18466446988997098217 > 2.16840434497100886801e-19)
Try another starting point or tweak arguments.
>>> from sympy import I
>>> nsolve(x**2 + 2, I) # Imaginary initial point returns a complex root
1.4142135623731*I
>>> nsolve(x**2 + 2, 1 + I) # Complex initial point returns a complex root
1.4142135623731*I
```

### Ensure the Root Found is in a Given Interval

It is not guaranteed that {func}`~.nsolve` will find the root closest to the
initial point. Here, even though the root `-1` is closer to the initial point of
`-0.1`, {func}`~.nsolve` finds the root `1`:

```py
>>> from sympy import nsolve
>>> from sympy.abc import x
>>> nsolve(x**2 - 1, -0.1)
1.00000000000000
```

You can ensure the root found is in a given interval, if such a root exists,
using `solver='bisect'` by specifying the interval in a tuple. Here, specifying
the interval `(-10, 0)` ensures that the root `-1` is found:

```py
>>> from sympy import nsolve
>>> from sympy.abc import x
>>> nsolve(x**2 - 1, (-10, 0), solver='bisect')
-1.00000000000000
```

### Solve a System of Equations Numerically

To solve a system of multidimensional functions, supply a tuple of
- functions `(f1, f2)`
- variables to solve for `(x1, x2)`
- starting values `(-1, 1)`

```py
>>> from sympy import Symbol, nsolve
>>> x1 = Symbol('x1')
>>> x2 = Symbol('x2')
>>> f1 = 3 * x1**2 - 2 * x2**2 - 1
>>> f2 = x1**2 - 2 * x1 + x2**2 + 2 * x2 - 8
>>> print(nsolve((f1, f2), (x1, x2), (-1, 1)))
Matrix([[-1.19287309935246], [1.27844411169911]])
```

### Increase Precision of the Solution

You can increase the precision of the solution using `prec`:

```py
>>> from sympy import Symbol, nsolve
>>> x1 = Symbol('x1')
>>> x2 = Symbol('x2')
>>> f1 = 3 * x1**2 - 2 * x2**2 - 1
>>> f2 = x1**2 - 2 * x1 + x2**2 + 2 * x2 - 8
>>> print(nsolve((f1, f2), (x1, x2), (-1, 1), prec=25))
Matrix([[-1.192873099352460791205211], [1.278444111699106966687122]])
```

### Create a Function That Can Be Solved With SciPy

As noted above, SymPy focuses on symbolic computation and is not optimized for
numerical calculations. If you need to make many calls to a numerical solver, it
can be much faster to use a solver optimized for numerical calculations such as
SciPy's {external:func}`~scipy.optimize.root_scalar`. A recommended workflow is:
1. use SymPy to generate (by symbolically simplifying or solving an equation)
  the mathematical expression
2. convert it to a lambda function using {func}`~.lambdify`
3. use a numerical library such as SciPy to generate numerical solutions

```py
>>> from sympy import simplify, cos, sin, lambdify
>>> from sympy.abc import x, y
>>> from scipy.optimize import root_scalar
>>> expr = cos(x * (x + x**2)/(x*sin(y)**2 + x*cos(y)**2 + x))
>>> simplify(expr) # 1. symbolically simplify expression
cos(x*(x + 1)/2)
>>> lam_f = lambdify(x, cos(x*(x + 1)/2)) # 2. lambdify
>>> sol = root_scalar(lam_f, bracket=[0, 2]) # 3. numerically solve using SciPy
>>> sol.root
1.3416277185114782
```

## Use the Solution Result

### Substitute the Result Into an Expression

The best practice is to use {func}`~sympy.core.evalf` to substitute numerical
values into expressions. The following code demonstrates that the numerical
value is not an exact root because substituting it back into the expression
produces a result slightly different from zero:

```py
>>> from sympy import cos, nsolve, Symbol
>>> x = Symbol('x')
>>> f = cos(x) - x
>>> x_value = nsolve(f, x, 1); x_value
0.739085133215161
>>> f.evalf(subs={x: x_value})
-5.12757857962640e-17
```

Using [`subs`](sympy.core.basic.Basic.subs) can give an incorrect result due to
precision errors, here effectively rounding `-5.12757857962640e-17` to zero:

```py
>>> f.subs(x, x_value)
0
```

When substituting in values, you can also leave some symbols as variables:

```py
>>> from sympy import cos, nsolve, Symbol
>>> x = Symbol('x')
>>> f = cos(x) - x
>>> x_value = nsolve(f, x, 1); x_value
0.739085133215161
>>> y = Symbol('y')
>>> z = Symbol('z')
>>> g = x * y**2
>>> values = {x: x_value, y: 1}
>>> (x + y - z).evalf(subs=values)
1.73908513321516 - z
```

## Not all Equations Can be Solved

{func}`~.nsolve` is a numerical solving function, so it can often provide a
solution for equations which cannot be solved algebraically.

### Equations With no Solution

Some equations have no solution, in which case SymPy may return an error. For
example, the equation $e^x = 0$ (`exp(x)` in SymPy) has no solution:

```py
>>> from sympy import nsolve, exp
>>> from sympy.abc import x
>>> nsolve(exp(x), x, 1, prec=20)
Traceback (most recent call last):
...
ValueError: Could not find root within given tolerance. (5.4877893607115270300540019e-18 > 1.6543612251060553497428174e-24)
Try another starting point or tweak arguments.
```

## Report a Bug

If you find a bug with {func}`~.nsolve`, please post the problem on the [SymPy mailing
list](https://groups.google.com/g/sympy). Until the issue is resolved, you can
use a different method listed in [](#alternatives-to-consider).
