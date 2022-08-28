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
>>> from sympy import symbols
>>> from sympy.matrices import Matrix
>>> c, d, e = symbols("c, d, e")
>>> A = Matrix([[c,d], [1, -e]])
>>> A
Matrix([
[c,  d],
[1, -e]])
>>> b = Matrix([2, 0])
>>> b
Matrix([
[2],
[0]])
>>> A.solve(b)
Matrix([
[2*e/(c*e + d)],
[  2/(c*e + d)]])
```

## Guidance

### Matrix Must Be Square

The matrix $A$ must be square to represent a system of linear equations with the
same number of unknowns as equations. If not, SymPy will give an error:

```py
ShapeError: `self` and `rhs` must have the same number of rows.
```

### *Guidance 2*

*Guidance 2 content*


## *Title*

You can *title* in several ways. 

### *Method 1*

*Method 1 content*

### *Method 2*

*Method 2 content*

## Working With Symbolic Matrices

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

## Equations With No Solution

If the determinant of a matrix is zero, matrix equations with it have no
solution:

```py
>>> from sympy import symbols
>>> from sympy.matrices import Matrix
>>> c, d, e = symbols("c, d, e")
>>> A = Matrix([[c*e**2, d*e], [c*e, d]])
>>> A
Matrix([
[c*e**2, d*e],
[   c*e,   d]])
>>> b = Matrix([2, 0])
>>> A.LUsolve(b)
Traceback (most recent call last):
    ...
NonInvertibleMatrixError: Matrix det == 0; not invertible.
```

## If You Encounter a Problem With SymPy Matrix Equation Solving

If you know your matrix equation has a solution, and SymPy cannot find it,
please post the problem on the [mailing
list](https://groups.google.com/g/sympy), or open an issue on [SymPy's GitHub
page](https://github.com/sympy/sympy/issues). Until the issue is resolved, you
can try one of the [](#alternatives-to-consider).
