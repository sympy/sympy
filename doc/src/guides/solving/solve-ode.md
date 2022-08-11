# Solve an Ordinary Differential Equation Algebraically

Use SymPy to solve an ordinary differential equation algebraically. For example,
solving $y''(x) + 9y(x)=0 $ yields $ y(x)=C_{1} \sin(3x)+ C_{2} \cos(3x)$.

Alternatives to consider:
- To solve a system of ordinary differential equations, use
  {func}`~.dsolve_system`
- *alternative 2*

Here is an example of solving an ordinary differential equation algebraically:

```py
>>> from sympy import Function, dsolve, Derivative
>>> from sympy.abc import x
>>> f = Function('f')
>>> result = dsolve(Derivative(f(x), x, x) + 9*f(x), f(x)); result
Eq(f(x), C1*sin(3*x) + C2*cos(3*x))
```

You can then use SymPy to verify that the solution is correct:

```py
>>> from sympy import checkodesol
>>> solution = result.rhs; solution
C1*sin(3*x) + C2*cos(3*x)
>>> checkodesol(Derivative(f(x), x, x) + 9*f(x), solution)
(True, 0)
```

The output of {func}`~.checkodesol` is a tuple where the first item, a boolean,
tells whether the substitution results in `0`.

## Guidance

### *Guidance 1*

*Guidance 1 content*

### *Guidance 2*

*Guidance 2 content*


## *Title*

You can *title* in several ways. 

### *Method 1*

*Method 1 content*

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

Please post the problem on the [mailing
list](https://groups.google.com/g/sympy), or open an issue on [SymPy's GitHub
page](https://github.com/sympy/sympy/issues). Until the issue is resolved, you
can *workaround*.
