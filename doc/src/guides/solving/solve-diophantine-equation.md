# Solve a Diophantine Equation Algebraically

Use SymPy to solve a Diophantine equation (find integer solutions to a
polynomial equation) algebraically, returning a parameterized general solution
if possible. For example, solving the [Pythagorean
theorem](https://en.wikipedia.org/wiki/Pythagorean_theorem) $a^2 + b^2 = c^2$
yields $(a=2pq, b=p^2-q^2, c=p^2-q^2)$.

## Alternatives to Consider

There are few alternatives for finding a parameterized general solution a
Diophantine equation. {func}`~.solve` simply solves for one variable in terms of
the others. For example, attempting to solve $a^2 + b^2 = c^2$ for $a$, $b$, and
$c$ can only reveal that $a = \pm \sqrt{c^2-b^2}$:

```py
>>> from sympy import solve
>>> from sympy import symbols
>>> a, b, c = symbols("a, b, c", integer=True)
>>> solve(a**2 + b**2 - c**2, [a, b, c], dict=True)
[{a: -sqrt(-b**2 + c**2)}, {a: sqrt(-b**2 + c**2)}]
```

## Example of Solving a Diophantine Equation

Here is an example of solving a Diophantine equation, specifically $a^2 + b^2 =
c^2$, using {func}`~.diophantine`:

```py
>>> from sympy.solvers.diophantine import diophantine
>>> from sympy import symbols
>>> a, b, c = symbols("a, b, c", integer=True)
>>> diophantine(a**2 + b**2 - c**2, syms=(a, b, c))
{(2*p*q, p**2 - q**2, p**2 + q**2)}
```

Refer to the [Diophantine API reference](../../modules/solvers/diophantine.rst)
for more examples of solving various types of Diophantine equations.

## Guidance

### *expr = 0*

### Specify the Order of Symbols in the Result

We recommend you specify the order of symbols in the result to avoid confusion.
Use the `syms` parameter and pass it a tuple or list of symbols to ensure the
result will be in that order.

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

## Use the Solution Result

### Extract Expressions From the Result

{func}`~.diophantine` returns results as a set of tuples, where each element in
a tuple is an expression for a variable in your equation. For example, in

```py
>>> from sympy.solvers.diophantine import diophantine
>>> from sympy import symbols
>>> a, b, c, p, q = symbols("a, b, c, p, q", integer=True)
>>> d = diophantine(a**2 + b**2 - c**2, syms=(a, b, c))
>>> d
{(2*p*q, p**2 - q**2, p**2 + q**2)}
```

the result is a set containing one tuple where the expressions correspond to (a,
b, c). That is, the tuple represents `a = 2*p*q, b = p**2 - q**2, c =
p**2-q**2`.

Because you cannot extract an element (here, a tuple) from a set by subscripting
the set, you can convert the set to a list, and then subscript the list:

```py
>>> solution_list = list(d)
>>> solution_list
[(2*p*q, p**2 - q**2, p**2 + q**2)]
>>> solution_list[0] # Extract a tuple corresponding to a solution
(2*p*q, p**2 - q**2, p**2 + q**2)
>>> solution_list[0][0] # Extract an expression for one variable, here a
2*p*q
```

You can also create a dictionary of symbol-expression pairs to extract an
expression by its symbol:

```py
>>> from sympy.solvers.diophantine import diophantine
>>> from sympy import symbols
>>> a, b, c = symbols("a, b, c", integer=True)
>>> my_syms = (a, b, c)
>>> solution, = diophantine(a**2 + b**2 - c**2, syms=(a, b, c))
>>> solution
(2*p*q, p**2 - q**2, p**2 + q**2)
>>> solution_dict = dict(zip(my_syms, solution))
>>> solution_dict
{a: 2*p*q, b: p**2 - q**2, c: p**2 + q**2}
>>> solution_dict[a]
2*p*q
```

### Work With Parameters

You can manipulate parameters such as `p` and `q`, which are generated
automatically by {func}`~.diophantine`, by creating them as symbols. For
example, to find a particular set of values that satisfies the Diophantine
equation, you can substitute in values for the parameters by
1. creating the parameters as symbols
2. substituting in their values using {meth}`~sympy.core.basic.Basic.subs`.

```py
>>> p, q = symbols("p, q", integer=True)
>>> [var.subs({p:4, q:3}) for var in solution_list[0]]
[24, 7, 25]
```

### Programmatically Extract Parameter Symbols

If you want to programmatically obtain the set of auto-generated parameters for
one solution, you can use the following code:

```py
>>> solution, = diophantine(a**2 + b**2 - c**2, syms=(a, b, c))
>>> solution
(2*p*q, p**2 - q**2, p**2 + q**2)
>>> set().union(*(s.free_symbols for s in solution))
{p, q}
```

## *Tradeoffs (speed vs. accuracy, etc.) for function*

### *Are there any tradeoffs to mention--permute parameter?*

*Tradeoff 1 content*

## Not All Equations Can Be Solved

### Equations With No Solution

Some Diophantine equations have no solution, in which case {func}`~.diophantine`
will return an empty set, `set()`. For example, in this equation, the
coefficients are both even ($2$ and $4$), so the sum of the terms ($2x + 4y$)
can only be even. However, the constant $3$ is odd, so there is no solution.

```py
>>> from sympy.solvers.diophantine import diophantine
>>> from sympy import symbols
>>> x, y = symbols("x, y", integer=True)
>>> diophantine(2*x + 4*y - 3, syms=(x, y))
set()
```

## Report a Problem

If you find a problem with {func}`~.diophantine`, please post the problem on the
[mailing list](https://groups.google.com/g/sympy), or open an issue on [SymPy's
GitHub page](https://github.com/sympy/sympy/issues). Until the issue is
resolved, you can use a different method listed in
[](#alternatives-to-consider).
