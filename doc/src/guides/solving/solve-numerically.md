# Solve one or a system of equations numerically

Use SymPy to numerically solve a system of one or more equations. For example,
numerically solving $\cos(x) = x $ returns $ x \approx 0.739085133215161$.

Solving numerically is useful if:
- You only need a numeric solution, not a symbolic one
- A closed-form solution is not available or is overly complicated

{func}`~.solve` and {func}`~.solveset` will not try to find a numeric solution,
only a mathematically exact symbolic solution. So if you want a numeric
solution, use {func}`~.nsolve`.

Alternatives to consider:
- [NumPy](https://numpy.org/doc/stable/reference/generated/numpy.linalg.solve.html?highlight=solve#numpy.linalg.solve)
and
[SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve.html#scipy.linalg.solve)
can each solve a system of linear scalar equations

Here is a simple example of numerically solving one equation:

```py
>>> from sympy import cos, nsolve, Symbol
>>> x = Symbol('x')
>>> nsolve(cos(x) - x, x, 1)
0.739085133215161
```

{func}`~.nsolve` calls, and can pass parameters to,
[mpmath.findroot](https://mpmath.org/doc/current/calculus/optimization.html#root-finding-findroot).

Overdetermined systems of equations are supported.

## Guidance

### Find complex roots of real function

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

### Ensure the root found is in a given interval

It is not guaranteed that nsolve will find the root closest to the initial
point. Here, even though the root `-1` is closer to the initial point of `-0.1`,
nsolve finds the root `1`:

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

## Solve a system of equations numerically

### Solve multidimensional functions

To solve multidimensional functions, supply a tuple of
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

### Use SciPy on a lambda function *what is the advantage of this method?*

You can *description*

```py
>>> from sympy.abc import x
>>> from sympy.utilities.lambdify import implemented_function
>>> from sympy import lambdify
>>> from scipy import optimize
>>> f = implemented_function('f', lambda x: x**3 - 1)
>>> lam_f = lambdify(x, f(x))
>>> sol = optimize.root_scalar(lam_f, bracket=[0, 2], method='brentq')
>>> sol.root
1.0
```

## Use the solution result

### Substitute the result into an expression

You can use [`subs`](sympy.core.basic.Basic.subs) to substitute numerical values
into expressions to numerically evaluate them fully:

```py
>>> from sympy import cos, nsolve, Symbol, diff
>>> x = Symbol('x')
>>> f = cos(x) - x
>>> x_value = nsolve(f, x, 1); x_value
0.739085133215161
>>> f.subs(x, x_value) # Verify that root is correct
0
>>> derivative = diff(f, x)
>>> derivative.subs(x, x_value) # Calculate the derivative where the function value is zero
-1.67361202918321
```

or leave some symbols as variables:

```py
>>> from sympy import cos, nsolve, Symbol, diff
>>> x = Symbol('x')
>>> f = cos(x) - x
>>> x_value = nsolve(f, x, 1); x_value
0.739085133215161
>>> y = Symbol('y')
>>> g = x * y**2
>>> g.subs(x, x_value)
0.739085133215161*y**2
>>> derivative = diff(g, y)
>>> derivative.subs(x, x_value)
1.47817026643032*y
```

## Tradeoffs

### Do not use `verify` for functions which are very steep near the root

For functions which are very steep near the root, the verification of the
solution may fail. In this case you should use the flag `verify=False` and
independently verify the solution.

```py
>>> from sympy import cos, cosh, nsolve, Symbol
>>> x = Symbol('x')
>>> f = cos(x)*cosh(x) - 1
>>> nsolve(f, 3.14*100)
Traceback (most recent call last):
...
ValueError: Could not find root within given tolerance. (1.39267e+230 > 2.1684e-19)
>>> ans = nsolve(f, 3.14*100, verify=False); ans
312.588469032184
>>> f.subs(x, ans).n(2)
2.1e+121
>>> (f/f.diff(x)).subs(x, ans).n(2)
7.4e-15
```

## Not all equations can be solved

{func}`~.nsolve` is a numerical solving function, so it is often the solution to
equations which cannot be solved algebraically.

### Equations with no solution

Some equations have no solution, in which case SymPy may return an empty set.
For example, the equation $x - 7 - x - 2 = 0$ reduces to $-9 = 0$, which has no
solution because no value of $x$ will make it true:

```py
>>> from sympy import nsolve
>>> from sympy.abc import x
>>> nsolve(x - 7 - x - 2, x, 1)
Traceback (most recent call last):
...
ValueError:
expected a one-dimensional and numerical function
```

SymPy reports that the function to be solved is not one-dimensional, because
SymPy simplifies $x - x$ to $0$, leaving a zero-dimensional function.

## Report a problem

If you find a problem with {func}`~.nsolve`, please post it on the [mailing
list](https://groups.google.com/g/sympy), or open an issue on [SymPy's GitHub
page](https://github.com/sympy/sympy/issues). Until the issue is resolved, you
can use a
[NumPy](https://numpy.org/doc/stable/reference/generated/numpy.linalg.solve.html?highlight=solve#numpy.linalg.solve)
or
[SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve.html#scipy.linalg.solve)
solver.