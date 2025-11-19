Matrix derivatives
==================

Matrix and array expressions in SymPy
-------------------------------------

SymPy supports matrices and N-dimensional arrays.  Matrices and N-dimensional
arrays have their basic component-explicit modules, i.e. ``sympy.matrices`` and
``sympy.tensor.array`` respectively.  By component-explicit it is meant that
the object contains the values of all components.

Besides, SymPy has support for matrix expressions and array expressions, these
are modules that define matrices and arrays as symbols, with only their
identifying name (such as the string `X` or `M`) and their shape tuple (in case
of matrices always a tuple of two elements).  The types of the shape tuples are
either integers, in case the shape is determined, or any kind mathematical
expression, whose integer value is not known yet.  Matrix symbols can undergo
the usual arithmetic operations, but unlike component-explicit matrices,
where operations happen elementwise to produce a resulting matrix, matrix
symbols will be stored in expression trees freezing their operations.  Matrix
symbols is the convention commonly used in mathematics textbooks.  In this case
expressions such as ``M*N.T*P`` (or `\mathbf{M} \mathbf{N}' \mathbf{P}`)
represent a matrix multiplication.

A further approach is to explicitly show indices after the symbol, in SymPy
both matrix symbols and array symbols support overloading for the ``[ ]``
operator to create matrix elements and array elements, respectively.  So you
can define ``M[i, j]`` to represent the value at the $(i+1)$-th row and
$(j+1)$-th column (SymPy index offset starts at zero, not at unit as commonly
assumed in mathematics).  Given that $i$ and $j$ are undefined variables, this
matrix element convention can be used to represent the whole span of
matrix/array values, and henceforth the whole matrix.  Still notice that SymPy
supports any kind of expression as index, no restrictions to numbers.

Matrix derivatives
------------------

A matrix derived by another matrix is a 4-dimensional array, combining the
components of the derivative of a matrix expression elements by the elements of
the matrix symbol.

The most intuitive understanding of matrix derivatives can be shown using
index-notation, if you have a matrix expression $A_{ij}$ derived by matrix
$X_{mn}$, the resulting derivative will be

$$\frac{\partial}{\partial X_{mn}} \Big ( A_{ij} \Big ) $$.

In SymPy, we use the convention of putting the derivation dimensions in front
of the derivand dimensions, convention in common with Wolfram Mathematica but
unlike most libraries such as PyTorch and NumPy. Therefore, when indices are
not explicitly stated, the derivative of $A_{ij}$ by matrix $X_{mn}$ will be in
the order ${}_{\{mnij\}}$.

In order to keep track of the index reading order (which may determine
transpositions of axes), we will use the convention of writing the index order
mapping to the full expression, so:

$${}_{\{mnij\}} \Longrightarrow \frac{\partial}{\partial X_{mn}} \Big ( A_{ij} \Big ) .$$

Notice that in this convention some expressions have multiple representations,
for example we can even write transpositions without using the operator,

$${}_{\{ij\}} \Longrightarrow \Big(M'\Big)_{ij} = M_{ji},$$

as transposition is indeed just the swapping of index reading order of two
indices.

As an example, the derivative of matrix $X$ by itself is given by

$${}_{\{mnij\}} \Longrightarrow \frac{\partial}{\partial X_{mn}} \Big ( X_{ij} \Big ) = \delta_{mi} \delta_{nj} = \Big( \mathbf{I} \otimes \mathbf{I} \Big)_{minj}, $$

where $\delta$ is the Kronecker delta symbol. Indeed, matrix $\mathbf{X}$ is
made of element variables that are to be considered different from one another,
so $X_{mn}$ is the same variable as $X_{ij}$ if and only if $m=i$ and $n=j$.
The previous expression also shows how the product of two Kronecker deltas on
four different free indices may be viewed as the tensor product of two identity
matrices, with the free index order properly permuted.

This matrix derivative returns a 4-dimensional array expression if computed in SymPy

>>> X.diff(X)
PermuteDims(ArrayTensorProduct(I, I), (3)(1 2))

The permutation `(3)(1 2)` maps ${}_{\{mnij\}}$ into ${}_{\{minj\}}$. Remember
that SymPy uses zero offset and the trivial cycle `(3)` is just a trick to
display the maximum size of the permutation.

Derivative is a matrix
----------------------

In SymPy, if you derive a matrix expression by a matrix symbol you will
generally get an array expression, as this has dimension 4.  In some cases, the
presence of trivial dimensions (i.e. axes of unit size) allows 4-dimensional
arrays to be represented as an equivalent matrix expression.

For example, if `x` is a matrix of shape `(k, 1)`, that is a vector-shaped
matrix, the derivative of $x$ by $x$, that is

$${}_{\{mnij\}} \Longrightarrow \frac{\partial}{\partial x_{mn}}x_{ij} = \delta_{mi} \delta_{nj} = \big( \mathbf{I}_{[k]} \otimes \mathbf{I}_{[1]} \big) \delta_{mi} = \big(\mathbf{I}_{[k]}\big)_{mi}$$

is an array of shape `(k, 1, k, 1)` that is equivalent to the identity matrix
$\mathbf{I}_{[k]}$ of size $k$. Indeed, in the index representation, we notice
that $\delta_{nj}$ is the scalar unit if trivial dimensions are dropped.
Similarly $\frac{\partial}{\partial \mathbf{x}}\mathbf{x}'$ has shape
`(k, 1, 1, k)`, as is still equivalent to the identity matrix $I_k$ of shape (k, k).

Matrix expressions and array expressions
----------------------------------------

The matrix expression module in SymPy reflects the common convention used in
mathematics to represent matrices when indices are not explicit. Therefore,
``M*N`` or `\mathbf{M}\cdot\mathbf{N}` is the matrix multiplication, or in
index-explicit form $\sum_j M_{ij} N_{jk}$ with final free indices $ik$.
Common operators such as determinant, trace, hadamard product and hadamard
power are provided.  Applying a function to a matrix is meant in a mathematical
way, not elementwise.  Matrix expressions has the ``.applyfunc`` method for
element-wise functions acting on their components.

The array expressions module on the other hand is meant to represent general
operations on arrays.

+---------------------------+-------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| operation name            | representation                            | index-explicit equivalent                                                                                 |
+===========================+===========================================+===========================================================================================================+
| tensor product            | $A \otimes B$                             | $A_{ij} B_{kl}$                                                                                           |
+---------------------------+-------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| contraction               | $A$ on axes $a$, $b$                      | $A_{i_{1} i_{2} \ldots i_{a} \ldots i_{b} \ldots } \Rightarrow \sum_j A_{i_{1} \ldots j \ldots j \ldots}$ |
+---------------------------+-------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| diagonal                  | $A$ on axes $a$, $b$                      | $A_{i_{1} i_{2} \ldots i_{a} \ldots i_{b} \ldots } \Rightarrow A_{i_{1} \ldots j \ldots j \ldots}$        |
+---------------------------+-------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| permutation of dimensions | permutation $\sigma$ on axes of $A$       | $A_{i_{1} i_{2} \ldots} \Rightarrow A_{i_{\sigma(1)} i_{\sigma(2) \ldots}}$                               |
+---------------------------+-------------------------------------------+-----------------------------------------------------------------------------------------------------------+

These operators on array expressions are handled in SymPy by expression tree
nodes ``ArrayTensorProduct``, ``ArrayContraction``, ``ArrayDiagonal``,
``PermuteDims``.

Matrix expressions can be represented by array expressions using these
operator, but the converse is not always true.

+---------------------------+-------------------------------------------+------------------------------------------------+----------------------------------------------------------------------------------+
| matrix operation          | matrix expression form                    | index form                                     | array expression form                                                            |
+===========================+===========================================+================================================+==================================================================================+
| matrix multiplication     | $\mathbf{M} \mathbf{N}$                   | ${}_{\{ij\}} \Rightarrow \sum_k M_{ik} N_{kj}$ | contraction: $M \otimes N$ on 2nd and 3rd axes                                   |
+---------------------------+-------------------------------------------+------------------------------------------------+----------------------------------------------------------------------------------+
| trace                     | $\mbox{tr}(\mathbf{M})$                   | ${}_{\{\}} \Rightarrow \sum_i M_{ii}$          | contraction: $M$ on 1nd and 2rd axes                                             |
+---------------------------+-------------------------------------------+------------------------------------------------+----------------------------------------------------------------------------------+
| diagonal                  | $\mbox{diag}(\mathbf{M})$                 | ${}_{\{i\}} \Rightarrow  M_{ii}$               | diagonalize: $M$ on 1nd and 2rd axes                                             |
+---------------------------+-------------------------------------------+------------------------------------------------+----------------------------------------------------------------------------------+
| transposition             | $\mathbf{M}'$                             | ${}_{\{ij\}} \Rightarrow M_{ji}$               | permutation: $M$ on 1nd and 2rd axes                                             |
+---------------------------+-------------------------------------------+------------------------------------------------+----------------------------------------------------------------------------------+
| Hadamard product          | $\mathbf{M} \circ \mathbf{N}$             | ${}_{\{ij\}} \Rightarrow M_{ij} N_{ij}$        | diagonalize: $M \otimes N$ on 1st-3rd axes and 2nd-4th axes                      |
+---------------------------+-------------------------------------------+------------------------------------------------+----------------------------------------------------------------------------------+

How matrix derivatives work
---------------------------

The matrix derivation algorithm takes a matrix expression, converts it to an
array expression equivalent using array operators.  At this point, the
derivation algorithm is straightforward by using a variation of the chain rule
that also keeps track of free indices and index order at each step.

After computing the derivative, an attempt at finding a matrix expression
equivalent to the derivative is made, dropping possible trivial dimensions.  If
no such matrix expression has been found, the array expression will be
returned.

Symbolic derivation algorithm finds a tensor expression that is equivalent to the derivative,
it does not perform the derivative in a way other platforms do with the chain rule.

Array expression to matrix expression conversion
------------------------------------------------

The core complexity of the algorithm of matrix derivation lies in the conversion back to matrix expression of the derivative.
The array derivative returns the derivative as an array expression with axes properly contracted, diagonalized and permuted.
In order to rewrite the output of a matrix expression, you need to identify the equivalent matrix operations on top of these expressions.
The most common of these operations in the matrix multiplication, but other operations such as traces, diag-expansion, hadamard products are also common.

The first limitation is that matrix derivatives.

References
----------

* The Matrix Cookbook, by Kaare Brandt Petersen and Michael Syskind Pedersen, https://www.math.uwaterloo.ca/~hwolkowi/matrixcookbook.pdf
