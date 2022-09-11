# Solve a Diophantine Equation Algebraically

Use SymPy to solve a Diophantine equation (find integer solutions to a
polynomial equation) algebraically, returning a parameterized general solution
if possible. For example, solving $a^2 + b^2 = c^2$ yields $(a=2pq, b=p^2-q^2,
c=p^2-q^2)$.

## Alternatives to Consider

- *solve() but not very useful, e.g. solve(a**2 + b**2 - c**2, [a, b, c])
  returns  [(-sqrt(-b**2 + c**2), b, c), (sqrt(-b**2 + c**2), b, c)]*
- *numerical method--numpy or scipy?*

Here is an example of solving a Diophantine equation, specifically the
[Pythagorean theorem](https://en.wikipedia.org/wiki/Pythagorean_theorem) $a^2 +
b^2 = c^2$, using {func}`~.diophantine`:

```py
>>> from sympy.solvers.diophantine import diophantine
>>> from sympy import symbols
>>> a, b, c = symbols("a, b, c", integer=True)
>>> diophantine(a**2 + b**2 - c**2)
>>> {(2*p*q, p**2 - q**2, p**2 + q**2)}
```

## Guidance

### Limitations

Currently, following five types of Diophantine equations can be solved using
{meth}`~sympy.solvers.diophantine.diophantine.diophantine` and other helper
functions of the Diophantine module.

- Linear Diophantine equations: $a_1x_1 + a_2x_2 + \ldots + a_nx_n = b$
- General binary quadratic equation: $ax^2 + bxy + cy^2 + dx + ey + f = 0$
- Homogeneous ternary quadratic equation: $ax^2 + by^2 + cz^2 + dxy + eyz + fzx
  = 0$
- Extended Pythagorean equation: $a_{1}x_{1}^2 + a_{2}x_{2}^2 + \ldots +
  a_{n}x_{n}^2 = a_{n+1}x_{n+1}^2$
- General sum of squares: $x_{1}^2 + x_{2}^2 + \ldots + x_{n}^2 = k$

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

### Equations With No Closed-Form Solution

*Equations with no analytical solution content*

## Report a Problem

If you find a problem with *function*, please post the problem on the [mailing
list](https://groups.google.com/g/sympy), or open an issue on [SymPy's GitHub
page](https://github.com/sympy/sympy/issues). Until the issue is resolved, you
can use a different method listed in [](#alternatives-to-consider).
