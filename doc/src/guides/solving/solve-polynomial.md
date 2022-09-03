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


## *Title*

You can *title* in several ways. 

### *Method 1*

*Method 1 content*

### *Method 2*

*Method 2 content*

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
