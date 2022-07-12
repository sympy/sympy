# Solve a system of inequalities algebraically

Use SymPy to solve a system of inequalities algebraically. For example, solving $x^2 < \pi$ and $x > 0$ yields $0 < x < \sqrt{pi}$.

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

### Systems of inequalities with no solution

If the system of inequalities has incompatible conditions, for example $x < 0$ and $x > \pi$, SymPy will return `False`:

```py
>>> from sympy import symbols, solve, pi
>>> x = symbols('x')
>>> solve([x < 0, x > pi])
False
```

### Equations with no analytical solution

SymPy may reflect that your equation has no solutions that can be expressed algebraically (symbolically) by returning an error such as `NotImplementedError`:

```py
>>> from sympy import symbols, solve, cos
>>> x = symbols('x')
>>> solve([cos(x) - x > 0, x > 0])
Traceback (most recent call last):
    ...
NotImplementedError: 
The inequality, -x + cos(x) > 0, cannot be solved using
solve_univariate_inequality.
```

so you may have to 
{func}`solve your equation numerically <sympy.solvers.solvers.nsolve>` 
instead.

### Equations which have an analytical solution, and SymPy cannot solve

SymPy has implemented algorithms to solve inequalities involving only one symbol ("variable"), so it cannot solve a set of inequalities involving more than one symbol:

```py
>>> from sympy import solve, symbols
>>> x, y = symbols('x y')
>>> from sympy.abc import x, y
>>> solve([x**2 < y, x > 0])
Traceback (most recent call last):
    ...
NotImplementedError: 
inequality has more than one symbol of interest.
```

You can use [Wolfram
Alpha](https://www.wolframalpha.com/input?i2d=true&i=solve%5C%2840%29Power%5Bx%2C2%5D+%3C+y+and+++x%3E0%5C%2844%29x%5C%2841%29)
to solve this problem.

If you encounter a problem with SymPy, please post the problem on the [mailing
list](https://groups.google.com/g/sympy), or open an issue on [SymPy's GitHub
page](https://github.com/sympy/sympy/issues). Until the issue is resolved, you
may be able to use [Wolfram
Alpha](https://www.wolframalpha.com/input?i2d=true&i=solve%5C%2840%29Power%5Bx%2C2%5D+%3C+y+and+++x%3E0%5C%2844%29x%5C%2841%29)
to solve the problem.
