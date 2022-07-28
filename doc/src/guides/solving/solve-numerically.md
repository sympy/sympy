# Solve a nonlinear equation system numerically

Use SymPy to numerically solve a system of one or more equations. For example, numerically solving $\cos(x) = x $ returns $ x \approx 0.739085133215161$.

Alternatives to consider:
- [NumPy](https://numpy.org/doc/stable/reference/generated/numpy.linalg.solve.html?highlight=solve#numpy.linalg.solve)
and
[SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve.html#scipy.linalg.solve)
can each solve a system of linear scalar equations

Here is a simple example of numerically solving an equation:

```py
>>> from sympy import solve, cos, Eq, nsolve, Symbol
>>> x = Symbol('x')
>>> nsolve(cos(x) - x, x, 1)
0.739085133215161
```

{func}`~.nsolve` is based on [mpmath.findroot](https://mpmath.org/doc/current/calculus/optimization.html#root-finding-findroot), and can pass parameters to it.

## Guidance

### Finding complex roots of real function

To solve for complex roots of real functions, a nonreal (either purely imaginary, or complex) initial point must be specified:

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

### *Guidance 2*

*Guidance 2 content*


## *Title*

You can *title* in several ways. 

### *Method 1*

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

### *Method 2*

*Method 2 content*

## Use the solution result

### *Usage method 1*

*Usage method 1 content*

### *Usage method 2*

*Usage method 2 content*

## *Tradeoffs (speed vs. accuracy, etc.) for function*

### *Tradeoff 1*

*Speed-up option 1 content*

### *Speed-up option 2*

*Speed-up option 2 content*

## Not all equations can be solved

### Equations with no solution

*Equations with no solution content*

### Equations with no analytical solution

*Equations with no analytical solution content*

### Equations which have an analytical solution, and SymPy cannot solve

*Equations which have an analytical solution, and SymPy cannot solve content*

Please post the problem on the 
[mailing list](https://groups.google.com/g/sympy), or open an issue on 
[SymPy's GitHub page](https://github.com/sympy/sympy/issues). Until the issue 
is resolved, you can *workaround*.
