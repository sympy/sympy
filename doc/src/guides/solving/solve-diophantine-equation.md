(solving-guide-diophantine)=
# Solve a Diophantine Equation Algebraically

Use SymPy to solve a [Diophantine
equation](https://en.wikipedia.org/wiki/Diophantine_equation) (find integer
solutions to a polynomial equation) algebraically, returning a parameterized
general solution if possible. For example, solving the [Pythagorean
equation](https://en.wikipedia.org/wiki/Pythagorean_theorem) $a^2 + b^2 = c^2$
yields $(a=2pq, b=p^2-q^2, c=p^2+q^2)$. Here, $p$ and $q$ are new parameters
introduced in the solution. $p$ and $q$ can take on any integer value to
parameterize the full set of solutions. More formally, $p,q \in \mathbb{Z}$
parameterize the infinite set of [Pythagorean
triples](https://en.wikipedia.org/wiki/Pythagorean_triple).

## Alternatives to Consider

There are few alternatives for finding a parameterized general solution a
Diophantine equation.
- Numerical alternatives:
    - [Sage's EllipticCurve
  command](https://doc.sagemath.org/html/en/constructions/elliptic_curves.html)
  may be able to find a set of relative numerical values for each variable
    - You can test explicit integer values, for example using a nested for loop
  of ranges of values. This is inefficient, but fine if you are only interested
  in solutions that are relatively small.
- {func}`~.solve` treats the variables as real or complex numbers, and simply
  solves for one variable in terms of the others, which produces a different
type of solution. For example, attempting to solve $a^2 + b^2 = c^2$ for $a$,
$b$, and $c$ can only reveal that $a = \pm \sqrt{c^2-b^2}$.

## Example of Solving a Diophantine Equation

Here is an example of solving a Diophantine equation, specifically $a^2 + b^2 =
c^2$, using {func}`~.diophantine`:

```py
>>> from sympy.solvers.diophantine import diophantine
>>> from sympy import symbols, Eq
>>> a, b, c = symbols("a, b, c", integer=True)
>>> my_syms = (a, b, c)
>>> pythag_eq = Eq(a**2 + b**2, c**2)
>>> # Solve Diophantine equation
>>> d = diophantine(pythag_eq, syms=my_syms)
>>> d
{(2*p*q, p**2 - q**2, p**2 + q**2)}
```

Refer to the [Diophantine API reference](../../modules/solvers/diophantine.rst)
for more examples of solving various types of Diophantine equations.

## Guidance

### Diophantine Equation Can be Expressed as Expression That Equals Zero

If you already have an expression that equals zero, you can solve that
expression. For example, expressing the Pythagorean equation as $a^2 + b^2 -
c^2$ is also valid:

```py
>>> from sympy.solvers.diophantine import diophantine
>>> from sympy import symbols
>>> a, b, c = symbols("a, b, c", integer=True)
>>> my_syms = (a, b, c)
>>> pythag = a**2 + b**2 - c**2
>>> diophantine(pythag, syms=my_syms)
{(2*p*q, p**2 - q**2, p**2 + q**2)}
```

### Specify the Order of Symbols in the Result

We recommend you specify the order of symbols in the result to avoid confusion.
Use the `syms` parameter and pass it a tuple or list of symbols to ensure the
result will be in that order, for example `syms=my_syms`, as in the examples on
this page.

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
a tuple is an expression for a variable in your equation. For example, for the
Pythogorean equation, the result is a set containing one tuple where the
expressions correspond to (a, b, c). That is, the tuple represents `a = 2*p*q, b
= p**2 - q**2, c = p**2-q**2`. Because you cannot extract an element (here, a
tuple) from a set by subscripting the set, you can create a dictionary of
symbol-expression pairs to extract an expression by its symbol:

```py
>>> from sympy.solvers.diophantine import diophantine
>>> from sympy import symbols
>>> a, b, c = symbols("a, b, c", integer=True)
>>> my_syms = (a, b, c)
>>> pythag = a**2 + b**2 - c**2
>>> solution, = diophantine(pythag, syms=my_syms)
>>> solution
(2*p*q, p**2 - q**2, p**2 + q**2)
>>> # Convert set to list
>>> solution_dict = dict(zip(my_syms, solution))
>>> solution_dict
{a: 2*p*q, b: p**2 - q**2, c: p**2 + q**2}
>>> # Extract an expression for one variable using its symbol, here a
>>> solution_dict[a]
2*p*q
```

Less elegantly, you can convert the set to a list, and then subscript the list.
It is a common mistake to forget the order of parameters, so this method is more
prone to errors:

```py
>>> from sympy.solvers.diophantine import diophantine
>>> from sympy import symbols
>>> a, b, c, p, q = symbols("a, b, c, p, q", integer=True)
>>> my_syms = (a, b, c)
>>> pythag = a**2 + b**2 - c**2
>>> d = diophantine(pythag, syms=my_syms)
>>> d
{(2*p*q, p**2 - q**2, p**2 + q**2)}
>>> # Convert set to list
>>> solution_list = list(d)
>>> solution_list
[(2*p*q, p**2 - q**2, p**2 + q**2)]
>>> # Extract a tuple corresponding to a solution
>>> solution_first = solution_list[0]
>>> solution_first
(2*p*q, p**2 - q**2, p**2 + q**2)
>>> # Extract an expression for one variable using its order, here a is element number zero
>>> solution_first[0]
2*p*q
```

### Work With Parameters

You can manipulate parameters such as `p` and `q`, which are generated
automatically by {func}`~.diophantine`, by creating them as symbols. For
example, to find a particular set of values that satisfies the Diophantine
equation, you can substitute in values for the parameters by
1. creating the parameters as symbols
2. substituting in their values using {meth}`~sympy.core.basic.Basic.subs`.

Here, we express the set of values as a dictionary to associate each variable
($a, b, c$) with its example value:

```py
>>> from sympy.solvers.diophantine import diophantine
>>> from sympy import symbols
>>> my_syms = (a, b, c)
>>> pythag = a**2 + b**2 - c**2
>>> d = diophantine(pythag, syms=my_syms)
>>> solution_list = list(d)
>>> solution_list
[(2*p*q, p**2 - q**2, p**2 + q**2)]
>>> p, q = symbols("p, q", integer=True)
>>> # Substitute in values as the dictionary is created
>>> solution_p4q3 = dict(zip(my_syms, [var.subs({p:4, q:3}) for var in solution_list[0]]))
>>> solution_p4q3
{a: 24, b: 7, c: 25}
```

Note that you need to include the `integer=True` assumption for the generated
parameters (`p` and `q`) to substitute numerical values for them. Conversely,
you do not need to include the `integer=True` assumption for the symbols in the
original equation (`a`, `b`, and `c`), although it is a good practice.

To iterate the set of solutions, you can iterate over value of the parameters
(`p` and `q`) in a nested loop:

```py
>>> from sympy.solvers.diophantine import diophantine
>>> from sympy import symbols
>>> a, b, c, p, q = symbols("a, b, c, p, q", integer=True)
>>> my_syms = (a, b, c)
>>> pythag = a**2 + b**2 - c**2
>>> d = diophantine(pythag, syms=my_syms)
>>> solution_list = list(d)
>>> # Iterate over the value of parameters p and q
>>> for p_val in range(-1,2):
...     for q_val in range(-1,2):
...         # Substitute in the values of p and q
...         pythag_vals = dict(zip(my_syms, [var.subs({p:p_val, q:q_val}) for var in solution_list[0]]))
...         # Print out the values of the generated parameters, and the Pythagorean triple a, b, c
...         print(f"p: {p_val}, q: {q_val} -> {pythag_vals}")
p: -1, q: -1 -> {a: 2, b: 0, c: 2}
p: -1, q: 0 -> {a: 0, b: 1, c: 1}
p: -1, q: 1 -> {a: -2, b: 0, c: 2}
p: 0, q: -1 -> {a: 0, b: -1, c: 1}
p: 0, q: 0 -> {a: 0, b: 0, c: 0}
p: 0, q: 1 -> {a: 0, b: -1, c: 1}
p: 1, q: -1 -> {a: -2, b: 0, c: 2}
p: 1, q: 0 -> {a: 0, b: 1, c: 1}
p: 1, q: 1 -> {a: 2, b: 0, c: 2}
```

### Verify a Solution

You can verify a solution is correct by substituting its integer values back
into the original equation (expression which equals zero) and checking that the
result is zero, either by using the dictionary approach from
[](#work-with-parameters), or by manually substituting in values determined by
any procedure:

```py
>>> from sympy.solvers.diophantine import diophantine
>>> from sympy import symbols
>>> a, b, c, p, q = symbols("a, b, c, p, q", integer=True)
>>> my_syms = (a, b, c)
>>> pythag = a**2 + b**2 - c**2
>>> d = diophantine(pythag, syms=my_syms)
>>> solution_list = list(d)
>>> solution_p4q3 = dict(zip(my_syms, [var.subs({p:4, q:3}) for var in solution_list[0]]))
>>> # Substitute values in using a dictionary
>>> pythag.subs({a: solution_p4q3[a], b: solution_p4q3[b], c: solution_p4q3[c]})
0
>>> # Manually substitute in values
>>> pythag.subs({a: 24, b: 7, c: 25})
0
```

### Programmatically Extract Parameter Symbols

If you want to programmatically obtain the set of auto-generated parameters for
one solution, you can use the following code:

```py
>>> from sympy.solvers.diophantine import diophantine
>>> from sympy import symbols
>>> a, b, c, p, q = symbols("a, b, c, p, q", integer=True)
>>> my_syms = (a, b, c)
>>> pythag = a**2 + b**2 - c**2
>>> # Solve Diophantine equation
>>> solution, = diophantine(pythag, syms=my_syms)
>>> solution
(2*p*q, p**2 - q**2, p**2 + q**2)
>>> # Extract parameter symbols
>>> set().union(*(s.free_symbols for s in solution))
{p, q}
```

## Not All Equations Can Be Solved

### Equations With No Solution

Some Diophantine equations have no solution, in which case {func}`~.diophantine`
will return an empty set, `set()`. For example, in the expression $2x + 4y - 3$
(which we will try to set to zero), the coefficients are both even ($2$ and
$4$), so the sum of the terms $(2x + 4y)$ can only be even. However, the
constant $3$ is odd, so there is no solution.

```py
>>> from sympy.solvers.diophantine import diophantine
>>> from sympy import symbols
>>> x, y = symbols("x, y", integer=True)
>>> diophantine(2*x + 4*y - 3, syms=(x, y))
set()
```

## Report a Bug

If you find a bug with {func}`~.diophantine`, please post the problem on the
[SymPy mailing list](https://groups.google.com/g/sympy). Until the issue is
resolved, you can use a different method listed in
[](#alternatives-to-consider).
