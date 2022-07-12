# Solve a system of inequalities algebraically

Use SymPy to solve a system of inequalities algebraically. For example, solving $x^2 < \pi$, $x > 0$ yields $0 < x < \sqrt{pi}$.

Alternatives to consider:
- *alternative 1*
- *alternative 2*

Here is a simple example of solving a system of inequalities algebraically:

```py
>>> from sympy import symbols, solve, pi
>>> from sympy.core.relational import Relational
>>> x = symbols('x')
>>> solve([x >= 0, x**2 <= pi])
(0 <= x) & (x <= sqrt(pi))
```

## Guidance

### *Guidance 1*

*Guidance 1 content*

### *Guidance 2*

*Guidance 2 content*


## Examples

### Using exact values

If you want to preserve the exact mathematical value of symbols, for example fractions, define them so that SymPy can interpret them symbolically, for example:

```py
>>> from sympy import symbols, solve, pi, Integer
>>> x = symbols('x')
>>> one_third = Integer(1)/Integer(3)
>>> solve([x >= one_third, x**2 <= pi])
```

If you use Python's standard division operator between the numerator and denominator, for example `1/3`, Python will perform the division and pass a numerical value to SymPy, leading to an inexact, numerical solution for all parts:

```py
>>> from sympy import symbols, solve, pi, Integer
>>> x = symbols('x')
>>> one_third = 1/3
>>> solve([x >= one_third, x**2 <= pi])
(0.333333333333333 <= x) & (x <= 1.77245385090552)
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

SymPy has implemented algorithms to solve inequalities involving only one symbol ("variable"), so it cannot solve a set of inequalities involving more than one symbol:

```py
>>> from sympy import solve, symbols
>>> x, y = symbols('x y')
>>> from sympy.abc import x, y
>>> solve([x**2 <= y])

```

Please post the problem on the 
[mailing list](https://groups.google.com/g/sympy), or open an issue on 
[SymPy's GitHub page](https://github.com/sympy/sympy/issues). Until the issue 
is resolved, you can *workaround*.
