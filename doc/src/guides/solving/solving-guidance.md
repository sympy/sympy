# Solving Guidance

Here are guidelines that apply to many types of solving.

## Guidance

### Use exact values

If you want to preserve the exact mathematical values of symbols such as
[fractions](tutorial-gotchas-final-notes) and [square
roots](symbolic-computation), define them so that SymPy can interpret them
symbolically, for example define one-third as a {class}`rational object <sympy.core.numbers.Rational>` using `Rational(1, 3)`:


```py
>>> from sympy import symbols, solve, pi, Rational
>>> x = symbols('x')
>>> one_third = Rational(1, 3)
>>> solve((x >= one_third, x**2 <= pi), x)
(1/3 <= x) & (x <= sqrt(pi))
```

If you use the division operator between standard Python numbers, for example
`1/3`, [Python will perform the division numerically](python-vs-sympy-numbers) and pass that inexact value to SymPy,
leading to an inexact, numerical expression for all relations:

```py
>>> from sympy import symbols, solve, pi, Integer
>>> x = symbols('x')
>>> one_third = 1/3
>>> solve([x >= one_third, x**2 <= pi], x)
(0.333333333333333 <= x) & (x <= 1.77245385090552)
```

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
