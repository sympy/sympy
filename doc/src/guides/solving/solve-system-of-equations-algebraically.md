# Solve a system of equations algebraically

Use SymPy to solve a system of equations algebraically. For example, $x^2 + y = 2, x - y = 4$ yields $\{(x = -3, y = -7), (x = 2, y = 2)\}$.

Alternatives to consider:
- *alternative 1*
- *alternative 2*

Here is a simple example of *title*:

```py
>>> from sympy import solve
>>> from sympy.abc import x, y
>>> solve([x**2 + y - 2, x - y - 4], x, y, dict=True)
[{x: -3, y: -7}, {x: 2, y: -2}]
```

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

Please post the problem on the 
[mailing list](https://groups.google.com/g/sympy), or open an issue on 
[SymPy's GitHub page](https://github.com/sympy/sympy/issues). Until the issue 
is resolved, you can *workaround*.
