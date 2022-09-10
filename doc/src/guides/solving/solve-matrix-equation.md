# Solve a Matrix Equation

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
  equations, so if you prefer you can use *link to solve a system of equations
  once published*
- If you formulated your problem as a system of linear equations, and want to
  convert it to matrix form, you can use {func}`~.linear_eq_to_matrix` and then
  follow the procedures in this guide.

## Solve a Matrix Equation

Here is an example of solving a matrix equation with SymPy's
{meth}`sympy.matrices.matrices.MatrixBase.solve`. We use the standard matrix
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

### Matrix Must Be Square

The matrix $A$ must be square to represent a system of linear equations with the
same number of unknowns as equations. If not, SymPy will give the error
``ShapeError: `self` and `rhs` must have the same number of rows.``

### Solving Several Matrix Equations With the Same Matrix

If you need to repeatedly solve matrix equations with the same matrix $A$ but
different constant vectors $b$, it is more efficient to use one of the following
methods: use [LU decomposition](https://en.wikipedia.org/wiki/LU_decomposition)
via {meth}`~sympy.matrices.matrices.MatrixBase.LUdecomposition`

```py
>>> from sympy import symbols, Matrix, eye
>>> c, d, e = symbols("c, d, e")
>>> A = Matrix([[c,d], [1, -e]])
>>> A
⎡c  d ⎤
⎢     ⎥
⎣1  -e⎦
>>> L, U, perm = A.LUdecomposition()
>>> L
    ⎡1  0⎤
    ⎢    ⎥
    ⎣c  1⎦
>>> U
    ⎡1    -e   ⎤
    ⎢          ⎥
    ⎣0  c⋅e + d⎦
>>> perm
[[0, 1]]
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
>>> P = eye(A.rows).permuteFwd(perm)
>>> P
    ⎡0  1⎤
    ⎢    ⎥
    ⎣1  0⎦
>>> y = L.solve(P*b) # Step-by-step approach, step 1
>>> y
    ⎡0⎤
    ⎢ ⎥
    ⎣2⎦
>>> U.solve(y) # Step-by-step approach, step 2
    ⎡  2⋅e  ⎤
    ⎢───────⎥
    ⎢c⋅e + d⎥
    ⎢       ⎥
    ⎢   2   ⎥
    ⎢───────⎥
    ⎣c⋅e + d⎦
>>> U.solve(L.solve(P*b)) # One-line approach
    ⎡  2⋅e  ⎤
    ⎢───────⎥
    ⎢c⋅e + d⎥
    ⎢       ⎥
    ⎢   2   ⎥
    ⎢───────⎥
    ⎣c⋅e + d⎦
>>> U.solve(L.solve(P*b2)) # Repeating one-line approach for b2
    ⎡  4⋅e  ⎤
    ⎢───────⎥
    ⎢c⋅e + d⎥
    ⎢       ⎥
    ⎢   4   ⎥
    ⎢───────⎥
    ⎣c⋅e + d⎦
```

or compute the inverse matrix using
{meth}`~sympy.matrices.matrices.MatrixBase.inv`:

```py
>>> from sympy import symbols, Matrix
>>> c, d, e = symbols("c, d, e")
>>> A = Matrix([[c,d], [1, -e]])
>>> b = Matrix([2, 0])
>>> b2 = Matrix([4, 0])
>>> inv = A.inv()
>>> inv
    ⎡   e        d   ⎤
    ⎢───────  ───────⎥
    ⎢c⋅e + d  c⋅e + d⎥
    ⎢                ⎥
    ⎢   1       -c   ⎥
    ⎢───────  ───────⎥
    ⎣c⋅e + d  c⋅e + d⎦
>>> inv * b # Solves to Ax = b for x
    ⎡  2⋅e  ⎤
    ⎢───────⎥
    ⎢c⋅e + d⎥
    ⎢       ⎥
    ⎢   2   ⎥
    ⎢───────⎥
    ⎣c⋅e + d⎦
>>> inv * b2 # Solves to Ax = b2 for x
    ⎡  4⋅e  ⎤
    ⎢───────⎥
    ⎢c⋅e + d⎥
    ⎢       ⎥
    ⎢   4   ⎥
    ⎢───────⎥
    ⎣c⋅e + d⎦
```

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
    A₀₀⋅A₁₁⋅A₂₂⋅A₃₃ - A₀₀⋅A₁₁⋅A₂₃⋅A₃₂ - A₀₀⋅A₁₂⋅A₂₁⋅A₃₃ + A₀₀⋅A₁₂⋅A₂₃⋅A₃₁ + A₀₀⋅A₁
    ₃⋅A₂₁⋅A₃₂ - A₀₀⋅A₁₃⋅A₂₂⋅A₃₁ - A₀₁⋅A₁₀⋅A₂₂⋅A₃₃ + A₀₁⋅A₁₀⋅A₂₃⋅A₃₂ + A₀₁⋅A₁₂⋅A₂₀⋅
    A₃₃ - A₀₁⋅A₁₂⋅A₂₃⋅A₃₀ - A₀₁⋅A₁₃⋅A₂₀⋅A₃₂ + A₀₁⋅A₁₃⋅A₂₂⋅A₃₀ + A₀₂⋅A₁₀⋅A₂₁⋅A₃₃ - 
    A₀₂⋅A₁₀⋅A₂₃⋅A₃₁ - A₀₂⋅A₁₁⋅A₂₀⋅A₃₃ + A₀₂⋅A₁₁⋅A₂₃⋅A₃₀ + A₀₂⋅A₁₃⋅A₂₀⋅A₃₁ - A₀₂⋅A₁
    ₃⋅A₂₁⋅A₃₀ - A₀₃⋅A₁₀⋅A₂₁⋅A₃₂ + A₀₃⋅A₁₀⋅A₂₂⋅A₃₁ + A₀₃⋅A₁₁⋅A₂₀⋅A₃₂ - A₀₃⋅A₁₁⋅A₂₂⋅
    A₃₀ - A₀₃⋅A₁₂⋅A₂₀⋅A₃₁ + A₀₃⋅A₁₂⋅A₂₁⋅A₃₀
```

and solving a matrix equation of it takes about a minute, whereas the analogous
3x3 matrix takes less than one second. The more unrelated, symbolic entries in a
matrix, the more likely it is to be slow to manipulate.

It is often much faster to solve as a system of linear
[](solvers-algebraic-equations) rather than as a symbolic matrix equation for
large matrix dimensions, so you may want to convert your matrix equation into a
system of linear symbolic equations.

One way to speed up matrix operations is by using the {class}`~.DomainMatrix`
class, which can be faster to operate on because it limits the domain of matrix
elements.

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

>>> (A * solution) - b # Not immediately obvious whether this result is a zeroes vector
    ⎡ 2⋅c⋅e      2⋅d      ⎤
    ⎢─────── + ─────── - 2⎥
    ⎢c⋅e + d   c⋅e + d    ⎥
    ⎢                     ⎥
    ⎣          0          ⎦
>>> simplify((A * solution) - b) # simplify reveals that this result is a zeroes vector
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

or you can extract individual elements by subscripting ("slicing")

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

## Report a Problem

If you know your matrix equation has a solution, and SymPy cannot find it,
please post the problem on the [mailing
list](https://groups.google.com/g/sympy), or open an issue on [SymPy's GitHub
page](https://github.com/sympy/sympy/issues). Until the issue is resolved, you
can try one of the [](#alternatives-to-consider).
