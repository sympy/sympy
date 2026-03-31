# Solve a Matrix Equation Algebraically

Use SymPy to solve a matrix (linear) equation. For example, solving $
\left[\begin{array}{cc} c & d\\1 & -e\end{array}\right] \left[\begin{array}{cc}
x\\y\end{array}\right] = \left[\begin{array}{cc} 2\\0\end{array}\right] $ yields
$ \left[\begin{array}{cc} x\\y\end{array}\right] = \left[\begin{array}{cc}
\frac{2e}{ce+d}\\\frac{2}{ce+d}\end{array}\right]$.

## Alternatives to Consider
- If your matrix and constant vector contain only numbers, not symbols, for
  example $\left[\begin{array}{cc} 1 & 2\\3 & 4\end{array}\right]
\left[\begin{array}{cc} x\\y\end{array}\right] = \left[\begin{array}{cc}
  2\\0\end{array}\right]$, you can use one of these other free and open-source
  packages instead of SymPy:
    - NumPy's {external:func}`numpy.linalg.solve`
    - SciPy's {external:func}`scipy.linalg.solve`
    - mpmath's
      [lu_solve()](https://mpmath.org/doc/current/matrices.html#linear-equations)
- Solving a matrix equation is equivalent to solving a system of linear
  equations, so if you prefer you can
  [](solve-system-of-equations-algebraically.md)
- If you formulated your problem as a system of linear equations, and want to
  convert it to matrix form, you can use {func}`~.linear_eq_to_matrix` and then
  follow the procedures in this guide.

## Solve a Matrix Equation

Here is an example of solving a matrix equation with SymPy's
{meth}`sympy.matrices.matrixbase.MatrixBase.solve`. We use the standard matrix
equation formulation $Ax=b$ where
- $A$ is the matrix representing the coefficients in the linear equations
- $x$ is the column vector of unknowns to be solved for
- $b$ is the column vector of constants, where each row is the value of an
  equation

```py
>>> from sympy import init_printing
>>> init_printing(use_unicode=True)
```

```py
>>> from sympy import symbols
>>> from sympy.matrices import Matrix
>>> c, d, e = symbols("c, d, e")
>>> A = Matrix([[c,d], [1, -e]])
>>> A
⎡c  d ⎤
⎢     ⎥
⎣1  -e⎦
>>> b = Matrix([2, 0])
>>> b
⎡2⎤
⎢ ⎥
⎣0⎦
>>> A.solve(b)
⎡  2⋅e  ⎤
⎢───────⎥
⎢c⋅e + d⎥
⎢       ⎥
⎢   2   ⎥
⎢───────⎥
⎣c⋅e + d⎦
```

## Guidance

### Matrix Usually Must Be Square

The matrix $A$ usually must be square to represent a system of linear equations
with the same number of unknowns as equations. If not, SymPy will give the error
``ShapeError: `self` and `rhs` must have the same number of rows.``

The exception to the requirement that a matrix be square comes from SymPy's use
of the {any}`Moore-Penrose pseudoinverse
<sympy.matrices.matrixbase.MatrixBase.pinv>`.

### Methods for Solving Matrix Equations

SymPy's matrix solving method, {meth}`sympy.matrices.matrixbase.MatrixBase.solve`,
can use several different methods, which are listed at that API reference link.
Depending on the nature of the matrix, a given method may be more efficient. By
default, [Gauss-Jordan
elimination](https://en.wikipedia.org/wiki/Gaussian_elimination) will be used.

Specifying a method in solve is equivalent to using a specialized solving
function. For example, using `solve` with `method='LU'` calls
{meth}`~sympy.matrices.matrixbase.MatrixBase.LUsolve`.

### Solving Several Matrix Equations With the Same Matrix

If you need to repeatedly solve matrix equations with the same matrix $A$ but
different constant vectors $b$, it is more efficient to use one of the following
methods.

You can use [LU decomposition](https://en.wikipedia.org/wiki/LU_decomposition)
via {meth}`~sympy.matrices.matrixbase.MatrixBase.LUsolve`:

```py
>>> from sympy import symbols, Matrix, eye, simplify
>>> c, d, e = symbols("c, d, e")
>>> A = Matrix([[c,d], [1, -e]])
>>> A
⎡c  d ⎤
⎢     ⎥
⎣1  -e⎦
>>> b = Matrix([2, 0])
>>> b
    ⎡2⎤
    ⎢ ⎥
    ⎣0⎦
>>> solution = A.LUsolve(b)
>>> solution
    ⎡  2⋅e  ⎤
    ⎢───────⎥
    ⎢c⋅e + d⎥
    ⎢       ⎥
    ⎢   2   ⎥
    ⎢───────⎥
    ⎣c⋅e + d⎦
>>> # Demonstrate that solution is correct
>>> simplify(A * solution)
    ⎡2⎤
    ⎢ ⎥
    ⎣0⎦
>>> b2 = Matrix([4, 0])
>>> b2
    ⎡4⎤
    ⎢ ⎥
    ⎣0⎦
>>> solution2 = A.LUsolve(b2)
>>> solution2
    ⎡  4⋅e  ⎤
    ⎢───────⎥
    ⎢c⋅e + d⎥
    ⎢       ⎥
    ⎢   4   ⎥
    ⎢───────⎥
    ⎣c⋅e + d⎦
>>> # Demonstrate that solution2 is correct
>>> simplify(A * solution2)
    ⎡4⎤
    ⎢ ⎥
    ⎣0⎦
```

Another approach is to compute the inverse matrix, but this is almost always
slower, and significantly slower for larger matrices. If efficient computation
is not a priority, you can use
{meth}`~sympy.matrices.matrixbase.MatrixBase.inv`:

```py
>>> from sympy import symbols, Matrix, simplify
>>> c, d, e = symbols("c, d, e")
>>> A = Matrix([[c,d], [1, -e]])
>>> b = Matrix([2, 0])
>>> b
    ⎡2⎤
    ⎢ ⎥
    ⎣0⎦
>>> b2 = Matrix([4, 0])
>>> b2
    ⎡4⎤
    ⎢ ⎥
    ⎣0⎦
>>> inv = A.inv()
>>> inv
    ⎡   e        d   ⎤
    ⎢───────  ───────⎥
    ⎢c⋅e + d  c⋅e + d⎥
    ⎢                ⎥
    ⎢   1       -c   ⎥
    ⎢───────  ───────⎥
    ⎣c⋅e + d  c⋅e + d⎦
>>> # Solves Ax = b for x
>>> solution = inv * b
>>> solution
    ⎡  2⋅e  ⎤
    ⎢───────⎥
    ⎢c⋅e + d⎥
    ⎢       ⎥
    ⎢   2   ⎥
    ⎢───────⎥
    ⎣c⋅e + d⎦
>>> # Demonstrate that solution is correct
>>> simplify(A * solution)
    ⎡2⎤
    ⎢ ⎥
    ⎣0⎦
>>> # Solves Ax = b2 for x
>>> solution2 = inv * b2
>>> solution2
    ⎡  4⋅e  ⎤
    ⎢───────⎥
    ⎢c⋅e + d⎥
    ⎢       ⎥
    ⎢   4   ⎥
    ⎢───────⎥
    ⎣c⋅e + d⎦
>>> # Demonstrate that solution2 is correct
>>> simplify(A * solution2)
    ⎡4⎤
    ⎢ ⎥
    ⎣0⎦
```

Determining the inverse of a large symbolic matrix may not be computationally
tractable.

### Work With Symbolic Matrices

The computational complexity of manipulating symbolic matrices can increase
rapidly with matrix size. For example, the number of terms in the determinant of
a symbolic matrix increases with the factorial of the matrix dimension. As a
result, the maximum dimensionality of matrices that can be solved is more
limited than for numerical matrices. For example, the determinant of this 4x4
symbolic matrix has 24 terms with four elements in each term:

```py
>>> from sympy import MatrixSymbol
>>> A = MatrixSymbol('A', 4, 4).as_explicit()
>>> A
⎡A₀₀  A₀₁  A₀₂  A₀₃⎤
⎢                  ⎥
⎢A₁₀  A₁₁  A₁₂  A₁₃⎥
⎢                  ⎥
⎢A₂₀  A₂₁  A₂₂  A₂₃⎥
⎢                  ⎥
⎣A₃₀  A₃₁  A₃₂  A₃₃⎦
>>> A.det()
A₀₀⋅A₁₁⋅A₂₂⋅A₃₃ - A₀₀⋅A₁₁⋅A₂₃⋅A₃₂ - A₀₀⋅A₁₂⋅A₂₁⋅A₃₃ + A₀₀⋅A₁₂⋅A₂₃⋅A₃₁ +
A₀₀⋅A₁₃⋅A₂₁⋅A₃₂ - A₀₀⋅A₁₃⋅A₂₂⋅A₃₁ - A₀₁⋅A₁₀⋅A₂₂⋅A₃₃ + A₀₁⋅A₁₀⋅A₂₃⋅A₃₂ +
A₀₁⋅A₁₂⋅A₂₀⋅A₃₃ - A₀₁⋅A₁₂⋅A₂₃⋅A₃₀ - A₀₁⋅A₁₃⋅A₂₀⋅A₃₂ + A₀₁⋅A₁₃⋅A₂₂⋅A₃₀ +
A₀₂⋅A₁₀⋅A₂₁⋅A₃₃ - A₀₂⋅A₁₀⋅A₂₃⋅A₃₁ - A₀₂⋅A₁₁⋅A₂₀⋅A₃₃ + A₀₂⋅A₁₁⋅A₂₃⋅A₃₀ +
A₀₂⋅A₁₃⋅A₂₀⋅A₃₁ - A₀₂⋅A₁₃⋅A₂₁⋅A₃₀ - A₀₃⋅A₁₀⋅A₂₁⋅A₃₂ + A₀₃⋅A₁₀⋅A₂₂⋅A₃₁ +
A₀₃⋅A₁₁⋅A₂₀⋅A₃₂ - A₀₃⋅A₁₁⋅A₂₂⋅A₃₀ - A₀₃⋅A₁₂⋅A₂₀⋅A₃₁ + A₀₃⋅A₁₂⋅A₂₁⋅A₃₀
```

and solving a matrix equation of it takes about a minute, whereas the analogous
3x3 matrix takes less than one second. The more unrelated, symbolic entries in a
matrix, the more likely it is to be slow to manipulate. This example, finding a
general solution to a matrix where all elements are independent symbols, is the
extreme case and thus the slowest for a matrix of its size.

### Speed up Solving Matrix Equations
Here are some suggestions:
- If matrix elements are zero, ensure that they are recognized as zero. You can
  do this by either making them zero or by applying
  [assumptions](assumptions_module).
- Selecting a solve method suited to the properties of the matrix, for example
  hermitian, symmetric, or triangular. Refer to
  [](#methods-for-solving-matrix-equations).
- Use the {class}`~.DomainMatrix` class, which can be faster to operate on
because it limits the domain of matrix elements.

## Use the Solution Result

### Use the Solution as a Vector

You can use the solution result as a vector. For example, to prove that the
solution $x$ is correct, you can multiply it the matrix $A$ and verify that it
produces the constants vector $b$:

```py
>>> from sympy import symbols, simplify
>>> from sympy.matrices import Matrix
>>> c, d, e = symbols("c, d, e")
>>> A = Matrix([[c,d], [1, -e]])
>>> b = Matrix([2, 0])
>>> solution = A.solve(b)
>>> solution
    ⎡  2⋅e  ⎤
    ⎢───────⎥
    ⎢c⋅e + d⎥
    ⎢       ⎥
    ⎢   2   ⎥
    ⎢───────⎥
    ⎣c⋅e + d⎦
>>> # Not immediately obvious whether this result is a zeroes vector
>>> (A * solution) - b
    ⎡ 2⋅c⋅e      2⋅d      ⎤
    ⎢─────── + ─────── - 2⎥
    ⎢c⋅e + d   c⋅e + d    ⎥
    ⎢                     ⎥
    ⎣          0          ⎦
>>> # simplify reveals that this result is a zeroes vector
>>> simplify((A * solution) - b)
    ⎡0⎤
    ⎢ ⎥
    ⎣0⎦
```

Note that we had to use {func}`~sympy.simplify.simplify.simplify` to make SymPy
simplify the expression in a matrix element to make it immediately obvious that
the solution is correct.

### Extract Elements From the Solution

Because you can iterate through the elements in a column vector, you can extract
its elements using standard Python techniques. For example, you can create a
list of the elements using list comprehension

```py
>>> [element for element in solution]
    ⎡  2⋅e       2   ⎤
    ⎢───────, ───────⎥
    ⎣c⋅e + d  c⋅e + d⎦
```

or you can extract individual elements by subscripting

```py
>>> solution[0]
      2⋅e
    ───────
    c⋅e + d
```

## Equations With No Solution

If the determinant of a matrix is zero, matrix equations with it have no
solution:

```py
>>> from sympy import symbols
>>> from sympy.matrices import Matrix
>>> c, d, e = symbols("c, d, e")
>>> A = Matrix([[c*e**2, d*e], [c*e, d]])
>>> A
    ⎡   2     ⎤
    ⎢c⋅e   d⋅e⎥
    ⎢         ⎥
    ⎣c⋅e    d ⎦
>>> b = Matrix([2, 0])
>>> A.LUsolve(b)
Traceback (most recent call last):
    ...
NonInvertibleMatrixError: Matrix det == 0; not invertible.
```

## Report a Bug

If you find a bug with matrix-solving functions, please post the problem on the
[SymPy mailing list](https://groups.google.com/g/sympy). Until the issue is resolved,
you can use a different method listed in [](#alternatives-to-consider).
