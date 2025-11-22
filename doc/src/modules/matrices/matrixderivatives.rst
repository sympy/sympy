Matrix derivatives
==================

Matrix and array expressions in SymPy
-------------------------------------

SymPy provides matrices and general N-dimensional arrays through the
modules ``sympy.matrices`` and ``sympy.tensor.array``, respectively. Each
of these two modules is component-explicit, meaning that every element is stored and directly
accessible, so the object explicitly contains the full set of component values
rather than representing them implicitly.

SymPy also supports symbolic matrix expressions and array expressions through
dedicated modules. These represent matrices and arrays abstractly as symbols,
defined solely by an identifying name (such as $X$ or $M$) and a shape tuple. For
matrices, the shape is always a 2-element tuple, representing the number of rows and columns. The dimensions within this
tuple can be either fixed integers (for known sizes) or arbitrary symbolic
expressions (when dimensions are undetermined).

Unlike component-explicit matrices, which perform elementwise calculations
immediately, symbolic matrix expressions retain operations unevaluated in
expression trees. Arithmetic operations (addition, multiplication, transpose,
etc.) are formally applied but frozen as symbolic representations rather than
computed results. This aligns with conventional mathematical notation in
textbooks, where expressions like ``M*N.T*P`` (also written as
$\mathbf{M}\mathbf{N}^{T} \mathbf{P}$) represent abstract matrix multiplication rather
than explicit numerical results.

Array expressions follow analogous principles, with the distinction that their
shape tuples support arbitrary sizes ($N$ dimensions) rather than only two.

In SymPy both matrix
symbols and array symbols support Python's ``[ ]`` operator to reference individual
elements. For example:

* ``M[i, j]`` denotes the element at row $i$ (0-indexed) and column $j$ of matrix symbol ``M``,
* ``A[k, l, m]`` accesses an element of a 3D array symbol ``A``.

SymPy uses zero-based indexing (starting at 0) for all positions, contrasting
with the one-based indexing (starting at 1) common in mathematical literature.
The indices themselves can be arbitrary symbolic expressions, not limited to
integers—allowing representations like ``M[2*i, j+1]``.  When combined with
symbolic indices (e.g., undefined variables $i$ and $j$), these element
references collectively span the entire matrix/array, facilitating
component-wise operations while preserving the abstract nature of the symbolic
expression.

Using SymPy code:

>>> from sympy import symbols, MatrixSymbol
>>> from sympy.tensor.array.expressions import ArraySymbol
>>> i, j, k, l, m = symbols("i j k l m")
>>> d = symbols("d")
>>> M = MatrixSymbol("M", d, d)
>>> M.shape
(d, d)
>>> M[i, j]
M[i, j]
>>> A = ArraySymbol("A", (d, d, d))
>>> A.shape
(d, d, d)
>>> A[k, l, m]
A[k, l, m]

Operations on matrix and array expressions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The matrix expression module in SymPy reflects the common convention used in
mathematics to represent matrices when indices are not explicit. Therefore,
``M*N`` or `\mathbf{M}\cdot\mathbf{N}` is the matrix multiplication, or in
index-explicit form $\sum_j M_{ij} N_{jk}$ with final free indices ${}_{\{ik\}}$.
Common operators such as determinant, trace, Hadamard product and Hadamard
power are provided.  Applying a function to a matrix is meant in a mathematical
way, not elementwise. Therefore ``exp(M)`` or `\exp(\mathbf{M})` represents the matrix
exponential of matrix $\mathbf{M}$. For the elementwise function application, matrix expressions have the ``.applyfunc`` method for
element-wise functions acting on their components, so applying the exponential
to every element of matrix $\mathbf{M}$ is achieved by calling ``M.applyfunc(exp)``,
commonly represented as $\exp_{\circ}(\mathbf{M})$ to distinguish it from the matrix function.

The array expressions module on the other hand is meant to represent general
operations on arrays.

+---------------------------+-------------------------+-------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| operation name            | SymPy object            | representation                            | index-explicit equivalent                                                                                 |
+===========================+=========================+===========================================+===========================================================================================================+
| tensor product            | ``ArrayTensorProduct``  | $A \otimes B$                             | $A_{ij} B_{kl}$                                                                                           |
+---------------------------+-------------------------+-------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| contraction               | ``ArrayContraction``    | $A$ on axes $a$, $b$                      | $A_{i_{1} i_{2} \ldots i_{a} \ldots i_{b} \ldots } \Rightarrow \sum_j A_{i_{1} \ldots j \ldots j \ldots}$ |
+---------------------------+-------------------------+-------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| diagonal                  | ``ArrayDiagonal``       | $A$ on axes $a$, $b$                      | $A_{i_{1} i_{2} \ldots i_{a} \ldots i_{b} \ldots } \Rightarrow A_{i_{1} \ldots j \ldots j \ldots}$        |
+---------------------------+-------------------------+-------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| permutation of dimensions | ``PermuteDims``         | permutation $\sigma$ on axes of $A$       | $A_{i_{1} i_{2} \ldots} \Rightarrow A_{i_{\sigma(1)} i_{\sigma(2) \ldots}}$                               |
+---------------------------+-------------------------+-------------------------------------------+-----------------------------------------------------------------------------------------------------------+

These operators on array expressions are handled in SymPy by expression tree
nodes ``ArrayTensorProduct``, ``ArrayContraction``, ``ArrayDiagonal``,
``PermuteDims``.

There is a canonical form for these three operators: where all these four occur, they should be nested
in these sequence (from outermost to innermost): permutation, diagonalization, contraction, tensor product.
Creating these objects itself does not resort the order, but it can be achieved by calling ``.doit()``.

For example, nested contractions can be flattened:

>>> from sympy import symbols
>>> from sympy.tensor.array.expressions import ArraySymbol, ArrayContraction, PermuteDims
>>> k = symbols("k")
>>> A = ArraySymbol("A", (k, k, k, k, k, k, k, k))
>>> expr = ArrayContraction(ArrayContraction(A, (1, 6), (2, 5)), (0, 3))
>>> expr
ArrayContraction(ArrayContraction(A, (1, 6), (2, 5)), (0, 3))
>>> expr.doit()
ArrayContraction(A, (0, 7), (1, 6), (2, 5))

>>> expr = ArrayContraction(PermuteDims(A, [6, 3, 0, 7, 1, 5, 4, 2]), (2, 4), (1, 7))
>>> expr
ArrayContraction(PermuteDims(A, (0 6 4 1 3 7 2)), (1, 7), (2, 4))
>>> expr.doit()
PermuteDims(ArrayContraction(A, (0, 1), (2, 3)), (0 2 1 3))

Matrix expressions can be represented by array expressions using these
operator, but the converse is not always true.

+---------------------------+-------------------------------------------+------------------------------------------------------------+----------------------------------------------------------------------------------+
| matrix operation          | matrix expression form                    | index form                                                 | array expression form                                                            |
+===========================+===========================================+============================================================+==================================================================================+
| matrix multiplication     | $\mathbf{M} \mathbf{N}$                   | ${}_{\{ij\}} \Rightarrow \sum_k M_{ik} N_{kj}$             | contraction: $M \otimes N$ on 2nd and 3rd axes                                   |
+---------------------------+-------------------------------------------+------------------------------------------------------------+----------------------------------------------------------------------------------+
| trace                     | $\mbox{tr}(\mathbf{M})$                   | ${}_{\{\}} \Rightarrow \sum_i M_{ii}$                      | contraction: $M$ on 1st and 2nd axes                                             |
+---------------------------+-------------------------------------------+------------------------------------------------------------+----------------------------------------------------------------------------------+
| diagonal                  | $\mbox{diag}(\mathbf{M})$                 | ${}_{\{i\}} \Rightarrow  M_{ii}$                           | diagonalize: $M$ on 1st and 2nd axes                                             |
+---------------------------+-------------------------------------------+------------------------------------------------------------+----------------------------------------------------------------------------------+
| transposition             | $\mathbf{M}'$                             | ${}_{\{ij\}} \Rightarrow M_{ji}$                           | permutation: $M$ on 1st and 2nd axes                                             |
+---------------------------+-------------------------------------------+------------------------------------------------------------+----------------------------------------------------------------------------------+
| Hadamard product          | $\mathbf{M} \circ \mathbf{N}$             | ${}_{\{ij\}} \Rightarrow M_{ij} N_{ij}$                    | diagonalize: $M \otimes N$ on 1st-3rd axes and 2nd-4th axes                      |
+---------------------------+-------------------------------------------+------------------------------------------------------------+----------------------------------------------------------------------------------+
| Kronecker product         | $\mathbf{M} \boxtimes \mathbf{N}$         | ${}_{\{m=id_1+k,n=jd_2+l\}} \Longrightarrow A_{ij} B_{kl}$ | permute $M \otimes N$ on 2nd and 3rd axes, then reshape                          |
+---------------------------+-------------------------------------------+------------------------------------------------------------+----------------------------------------------------------------------------------+

Kronecker product versus tensor product
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In SymPy, the Kronecker product $\boxtimes$ and tensor product $\otimes$ are two ways of seeing the same underlying operation.
Given matrices $\mathbf{A}$ and $\mathbf{B}$, their elements can be combined as

$$\mathbf{A} \otimes \mathbf{B} \Longrightarrow A_{ij} B_{kl} $$

When representing this expression as an array, there are different approaches depending on the index-order and whether to reshape the 4-dimensional
array as a 2-dimensional array by joining dimensions.

The array tensor product as implemented by ``ArrayTensorProduct`` and ``tensorproduct`` simply nests the arrays,
maintaining the index-order ${}_{\{ijkl\}}$, and produces:

$$\left[
\begin{matrix}\left[\begin{matrix}{A}_{0,0} {B}_{0,0} & {A}_{0,0} {B}_{0,1}\\{A}_{0,0} {B}_{1,0} & {A}_{0,0} {B}_{1,1}\end{matrix}\right] & \left[\begin{matrix}{A}_{0,1} {B}_{0,0} & {A}_{0,1} {B}_{0,1}\\{A}_{0,1} {B}_{1,0} & {A}_{0,1} {B}_{1,1}\end{matrix}\right]\\\left[\begin{matrix}{A}_{1,0} {B}_{0,0} & {A}_{1,0} {B}_{0,1}\\{A}_{1,0} {B}_{1,0} & {A}_{1,0} {B}_{1,1}\end{matrix}\right] & \left[\begin{matrix}{A}_{1,1} {B}_{0,0} & {A}_{1,1} {B}_{0,1}\\{A}_{1,1} {B}_{1,0} & {A}_{1,1} {B}_{1,1}\end{matrix}\right]\end{matrix}
\right] $$

The Kronecker product can be defined in terms of indices ${}_{\{mn\}}$ where $m$ spans over rows of both $\mathbf{A}$
and $\mathbf{B}$, while $n$ spans over their columns:
$$A \boxtimes B = \Big[ {}_{\{mn\}} = {}_{\{m=id_1+k,n=jd_2+l\}} \Longrightarrow A_{ij} B_{kl} \Big]$$
where $[d_1 \times d_2]$ is the shape of $\mathbf{A}$.

The Kronecker product, implemented by ``KroneckerProduct`` and ``kronecker_product``, combines the rows and columns of $\mathbf{A}$ and $\mathbf{B}$

$$\left[
\begin{matrix}{A}_{0,0} {B}_{0,0} & {A}_{0,0} {B}_{0,1} & {A}_{0,1} {B}_{0,0} & {A}_{0,1} {B}_{0,1}\\{A}_{0,0} {B}_{1,0} & {A}_{0,0} {B}_{1,1} & {A}_{0,1} {B}_{1,0} & {A}_{0,1} {B}_{1,1}\\{A}_{1,0} {B}_{0,0} & {A}_{1,0} {B}_{0,1} & {A}_{1,1} {B}_{0,0} & {A}_{1,1} {B}_{0,1}\\{A}_{1,0} {B}_{1,0} & {A}_{1,0} {B}_{1,1} & {A}_{1,1} {B}_{1,0} & {A}_{1,1} {B}_{1,1}\end{matrix}
\right]$$

In these visual representations of components, the Kronecker product and tensor product look very similar, apart from nesting.
However, using the ``Reshape`` operator to convert shape $[2 \times 2 \times 2 \times 2]$ into an array of shape $[4 \times 4]$
will preserve the inner-outer index order, therefore the inner matrices will be converted into rows display a different order
from the Kronecker product.
To convert a tensor product to a Kronecker product by reshaping
requires permuting the $j$ and $k$ indices first:

$$\left[
\begin{matrix}\left[\begin{matrix}{A}_{0,0} {B}_{0,0} & {A}_{0,0} {B}_{0,1}\\{A}_{0,1} {B}_{0,0} & {A}_{0,1} {B}_{0,1}\end{matrix}\right] & \left[\begin{matrix}{A}_{0,0} {B}_{1,0} & {A}_{0,0} {B}_{1,1}\\{A}_{0,1} {B}_{1,0} & {A}_{0,1} {B}_{1,1}\end{matrix}\right]\\\left[\begin{matrix}{A}_{1,0} {B}_{0,0} & {A}_{1,0} {B}_{0,1}\\{A}_{1,1} {B}_{0,0} & {A}_{1,1} {B}_{0,1}\end{matrix}\right] & \left[\begin{matrix}{A}_{1,0} {B}_{1,0} & {A}_{1,0} {B}_{1,1}\\{A}_{1,1} {B}_{1,0} & {A}_{1,1} {B}_{1,1}\end{matrix}\right]\end{matrix}
\right]
$$

This array can be reshaped as a $[4 \times 4]$ Kronecker product of $\mathbf{A}$ and $\mathbf{B}$, the element order of the inner matrices
is indeed equivalent to the rows of the Kronecker product.

>>> from sympy import MatrixSymbol, kronecker_product
>>> A = MatrixSymbol("A", 2, 2).as_explicit()
>>> B = MatrixSymbol("B", 2, 2).as_explicit()
>>> kronecker_product(A, B)
Matrix([
[A[0, 0]*B[0, 0], A[0, 0]*B[0, 1], A[0, 1]*B[0, 0], A[0, 1]*B[0, 1]],
[A[0, 0]*B[1, 0], A[0, 0]*B[1, 1], A[0, 1]*B[1, 0], A[0, 1]*B[1, 1]],
[A[1, 0]*B[0, 0], A[1, 0]*B[0, 1], A[1, 1]*B[0, 0], A[1, 1]*B[0, 1]],
[A[1, 0]*B[1, 0], A[1, 0]*B[1, 1], A[1, 1]*B[1, 0], A[1, 1]*B[1, 1]]])
>>> tensorproduct(A, B)
[[[[A[0, 0]*B[0, 0], A[0, 0]*B[0, 1]], [A[0, 0]*B[1, 0], A[0, 0]*B[1, 1]]], [[A[0, 1]*B[0, 0], A[0, 1]*B[0, 1]], [A[0, 1]*B[1, 0], A[0, 1]*B[1, 1]]]], [[[A[1, 0]*B[0, 0], A[1, 0]*B[0, 1]], [A[1, 0]*B[1, 0], A[1, 0]*B[1, 1]]], [[A[1, 1]*B[0, 0], A[1, 1]*B[0, 1]], [A[1, 1]*B[1, 0], A[1, 1]*B[1, 1]]]]]

Matrix derivatives
------------------

You can differentiate a matrix with respect to a scalar by differentiating each
element individually. Similarly, differentiating a scalar with respect to a
matrix yields a matrix of the same shape, where each entry is the derivative of
the scalar with respect to the corresponding matrix element.

While the trace and determinant of a matrix are scalars, their derivatives with
respect to a matrix variable can yield nontrivial expressions involving the
matrix itself.

Differentiating a matrix $\mathbf{Y}$ with respect to another matrix
$\mathbf{X}$ results in a 4-dimensional array (sometimes also called a tensor,
though tensors refer to differential geometry in SymPy), where each entry is
the partial derivative of an element of $\mathbf{Y}$ with respect to an element
of $\mathbf{X}$.

Derivatives extend naturally through index notation.  Differentiating a matrix
expression $A_{ij}$ by matrix $X_{mn}$ produces a four-dimensional array
representing:

$$D_{mnij} = \frac{\partial}{\partial X_{mn}} \Big (A_{ij} \Big ).$$

SymPy adopts a denominator-first index ordering for derivatives, positioning
differentiation variable indices ${}_{\{mn\}}$ before derivand indices ${}_{\{ij\}}$.  This
${}_{\{mnij\}}$ convention aligns with Wolfram Mathematica but differs from
PyTorch/NumPy, which prefer the ${}_{\{ijmn\}}$ index-order convention.
This structure applies universally: scalar-by-matrix
derivatives produce matrices, matrix-by-scalar derivatives produce matrices,
and matrix-by-matrix derivatives yield 4-dimensional arrays reflecting the
relationship between all component pairs.

To explicitly track index reading order, which controls axis transpositions in
multi-dimensional representations, we will use the convention of mapping the index
sequence to the full expression. For example:

$${}_{\{mnij\}} \Longrightarrow\frac{\partial}{\partial X_{mn}} \Big ( A_{ij} \Big )$$

this shows that
differentiation indices ($mn$) precede derivand indices ($ij$).  Crucially,
this notation inherently captures index permutations without explicit
operators. For example:

$${}_{\{ij\}} \Longrightarrow \Big(M^T\Big)_{ij} = M_{ji}$$

demonstrates how transposition is encoded solely through index-order
reversal, the ${}_{\{ij\}}$ mapping of $M_{ji}$ directly yields the
component expression of the transposition.

Reconstructing a matrix expression
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In SymPy, differentiating a matrix expression with respect to a matrix symbol
normally yields a 4-dimensional array.  When one or more of its axes happen to be
singleton (dimension 1), the tensor collapses to a lower-rank object; SymPy
automatically detects this and returns the simplest equivalent matrix
expression instead of keeping the full 4-dimensional form.

For example, let $\mathbf{x}$ be a matrix-vector of shape $[k \times 1]$, that is a
vector of size $[k]$ disguised as a single column matrix by adding a singleton dimension to its shape.  The derivative of $\mathbf{x}$ with
respect to $\mathbf{x}$ is then the 4-D identity array of shape $[k \times 1 \times k \times 1]$, but because the second
axis of both the numerator and the denominator has length 1, the result can be
collapsed to the $[k \times k]$ identity matrix $\mathbf{I}_{[k]}$ instead of
keeping the full tensor form. Explicitly:

$$\left[{}_{\{mnij\}} \Longrightarrow \frac{\partial}{\partial x_{mn}}x_{ij} = \delta_{mi} \delta_{nj} = \big( \mathbf{I}_{[k]} \otimes \mathbf{I}_{[1]} \big)_{minj}\right] \rightarrow \Big[ {}_{\{mi\}}\Longrightarrow \big(\mathbf{I}_{[k]}\big)_{mi}\Big]$$

Here, the second and fourth axes are singletons, appearing in the derivative as
scalar identity matrix $\mathbf{I}_{[1]}$, equivalent to scalar unit.  In index
notation they only survive in a Kronecker $\delta_{nj}$.

Likewise, $\frac{\partial}{\partial \mathbf{x}}\mathbf{x}'$ has shape
$[k \times 1 \times 1 \times k]$, and collapses again to the same $[k \times k]$ identity matrix $\mathbf{I}_{[k]}$.

Derive matrix by itself
~~~~~~~~~~~~~~~~~~~~~~~

As an example, the derivative of $[k \times l]$ matrix $\mathbf{X}$ by itself is given by identity
relationships between indices:

$${}_{\{mnij\}} \Longrightarrow \frac{\partial}{\partial X_{mn}} \Big ( X_{ij} \Big ) = \delta_{mi} \delta_{nj} = \Big( \mathbf{I}_{[k]} \otimes \mathbf{I}_{[l]} \Big)_{minj}, $$

where $\delta$ is the Kronecker delta, which satisfies $\delta_{ab} = 1$ if
$a = b$ and $\delta_{ab} = 0$ otherwise.  Indeed, matrix $\mathbf{X}$ is made of
element variables that are to be considered different from one another, so
$X_{mn}$ is the same variable as $X_{ij}$ if and only if $m=i$ and $n=j$.  The
previous expression also shows how the product of two Kronecker deltas on four
different free indices may be viewed as the tensor product of two identity
matrices, with the free index order properly permuted.

This matrix derivative returns a 4-dimensional array expression if computed in
SymPy

>>> from sympy import MatrixSymbol
>>> X = MatrixSymbol("X", 3, 3)
>>> X.diff(X)
PermuteDims(ArrayTensorProduct(I, I), (3)(1 2))

The permutation `(3)(1 2)` maps ${}_{\{mnij\}}$ into ${}_{\{minj\}}$. Remember
that SymPy uses zero offset and the trivial cycle `(3)` is just a trick to
display the maximum size of the permutation. The object ``PermuteDims`` in SymPy
takes the convention of its permutation mapping the new index order into the old index order
(the old index order of the wrapped expression). This convention is compatible with most of
the scientific python libraries, but is the opposite of the one used by Wolfram Mathematica.

Notice that the same result, ``X.diff(X)`` can be expressed as

``Reshape(KroneckerProduct(I, I), (k, l, k, l))``,

indeed some authors represent $\frac{\partial}{\partial \mathbf{X}}\mathbf{X} = \mathbf{I}_{[k]} \bar\boxtimes \mathbf{I}_{[l]}$,
where $\bar\boxtimes$ represents the 4-dimensionally reshaped Kronecker product.

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

Chain rule
~~~~~~~~~~

Chain-rule for matrix derivatives is applied in a similar manner to the
standard derivative chain rule, but in case the dimension of the intermediate
gradient increases, the axes have to properly rearranged to the right order.

Given $\mathbf{Y}$ (or $\mathbf{Z}$) as a generic matrix expression, the
derivative of a function $\partial f(\mathbf{Y})$ can be expanded with the
chain rule of $f$ to the intermediate expression in terms of $\partial
\mathbf{Y}$.  The usual non-commuting principle holds for matrices, so the
relative order of $\partial \mathbf{Y}$ and $\mathbf{Y}$ cannot be changed.

If you are deriving by a scalar, $\partial \mathbf{Y}$ will not increase the
number of dimensions, so no permutations of axes will be required.

In case you are deriving by a matrix, $\partial \mathbf{Y}$ will be a
4-dimensional array, the indices coming from $\mathbf{Y}$ will be connected to
the chain rule expression, while the deriving indices need to be brought in
front of all others (remember, we use the convention that the indices of the
deriving variable precede the indices of the expression to be derived).

+---------------------------+---------------------------------+-------------------------------------------------------------------------------------+
| operation                 | expression                      | chain rule                                                                          |
+===========================+=================================+=====================================================================================+
| matrix addition           | $\mathbf{Y} + \mathbf{Z}$       | $\partial \mathbf{Y} + \partial\mathbf{Z}$                                          |
+---------------------------+---------------------------------+-------------------------------------------------------------------------------------+
| matrix multiplication     | $\mathbf{Y}\mathbf{Z}$          | $(\partial \mathbf{Y})\mathbf{Z} + \mathbf{Y} (\partial\mathbf{Z})$                 |
+---------------------------+---------------------------------+-------------------------------------------------------------------------------------+
| Hadamard product          | $\mathbf{Y} \circ \mathbf{Z}$   | $(\partial \mathbf{Y})\circ\mathbf{Z} + \mathbf{Y}\circ(\partial\mathbf{Z})$        |
+---------------------------+---------------------------------+-------------------------------------------------------------------------------------+
| tensor product            | $\mathbf{Y}\otimes\mathbf{Z}$   | $(\partial \mathbf{Y})\otimes\mathbf{Z} + \mathbf{Y}\otimes(\partial\mathbf{Z})$    |
+---------------------------+---------------------------------+-------------------------------------------------------------------------------------+
| inverse                   | $\mathbf{Y}^{-1}$               | $-\mathbf{Y}^{-1} (\partial \mathbf{Y}) \mathbf{Y}^{-1}$                            |
+---------------------------+---------------------------------+-------------------------------------------------------------------------------------+
| trace                     | $\mbox{tr}(\mathbf{Y})$         | $\mbox{tr}(\partial\mathbf{Y})$                                                     |
+---------------------------+---------------------------------+-------------------------------------------------------------------------------------+
| determinant               | $\mbox{det}(\mathbf{Y})$        | $\mbox{det}(\mathbf{Y}) \mbox{tr}(\mathbf{Y}^{-1}\partial\mathbf{Y})$               |
+---------------------------+---------------------------------+-------------------------------------------------------------------------------------+
| transposition             | $\mathbf{Y}'$                   | $(\partial\mathbf{Y})'$                                                             |
+---------------------------+---------------------------------+-------------------------------------------------------------------------------------+

When differentiating with respect to a matrix $\mathbf{X}$, we treat the
gradient $\partial \mathbf{Y}$ as a 4-dimensional array whose first two axes
are the indices of $\mathbf{X}$ (row, column).  After the chain rule step for
the derivative is computed, these leading axes have to be permuted so that the
final tensor carries the indices in the order expected by the downstream matrix
or array expression.

Chain rule in array expressions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the operands are N-dimensional arrays, matrix multiplication has to be
rewritten as a tensor product followed by a contraction along the appropriate
axes.  The chain rule for the tensor product is the familiar one, but
contractions and diagonal extractions need extra care: their chain rules have
to shuffle the axes they act on so that the summation or diagonalization
indices still line up correctly after the derivative.

SymPy computes a matrix derivative by first lifting the whole expression to the
array level, differentiating there, and then folding the result back into the
simplest matrix form that respects the original axis structure.

Example: Derivative of the inverse matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For example, when **deriving the inverse** of matrix expression $\mathbf{Y}$ is

$$\frac{\partial \mathbf{Y}_{ij}^{-1}}{\partial \mathbf{X}_{mn}} = \left( {}_{\{mnij\}} \Rightarrow - \sum_{kl} \mathbf{Y}^{-1}_{ik} \frac{\partial \mathbf{Y}_{kl}}{\partial \mathbf{X}_{mn}} \mathbf{Y}^{-1}_{lj} \right)$$

using our array expression syntax, this can be represented as

|  PermuteDims(
|     ArrayTensorContraction(
|        ArrayTensorProduct(
|           -Inverse(Y), array_derive(Y, X), Inverse(Y)
|        ),
|        (1, 4),
|        (5, 6)
|     ),
|     index_order_old="imnj",
|     index_order_new="mnij"
|  )

here, ``array_derive(Y, X)`` returns a 4-dimensional array expression, we then proceed to create an 8-dimensional
array expression,
$$-\mathbf{Y}^{-1} \otimes \frac{\partial\mathbf{Y}}{\partial\mathbf{X}} \otimes \mathbf{Y}^{-1}$$
which is then contracted on the 2-nd and 5-th axes, i.e. (1, 4), and on the 6-th and 7-th axes, i.e. (5, 6).
$$-\frac{1}{\partial\mathbf{X}} \Big( \mathbf{Y}^{-1} (\partial\mathbf{Y}) \mathbf{Y}^{-1} \Big)$$
The contraction reproduce the structure of the matrix multiplication of the chain rule for the inverse.
As a last step, we need to take the 2-nd and 3-rd axes of the contracted expression in front of the other ones,
as they are the axes referring to the deriving variable $\mathbf{X}$.

Array expression to matrix expression conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The core complexity of the algorithm of matrix derivation lies in the conversion back to matrix expression of the derivative.
The array derivative returns the derivative as an array expression with axes properly contracted, diagonalized and permuted.
In order to rewrite the output of a matrix expression, you need to identify the equivalent matrix operations on top of these expressions.
The most common of these operations in the matrix multiplication, but other operations such as traces, diag-expansion, Hadamard products are also common.

Matrix derivation may produce array expressions that have no correspondence to closed-form matrix expressions,
therefore it will not be always possible to re-express the result as a matrix expression.

Sometimes a clear sequence of paired contractions corresponding to a matrix multiplication line will be identified,
in these cases the identification of a matrix expression is straightforward.

Otherwise, some tricks are used:

* tensor product of two matrix vectors $\mathbf{a} \otimes \mathbf{b}$, of shape $[k \times 1]$ both, can be turned into a matrix multiplication over their trivial dimension: $\mathbf{a} \cdot \mathbf{b}'$.
* identity matrices in a tensor product may be removed. Indeed they increase the dimensions of the array expression by filling the space with null values off the diagonals identified by their axes, they can thus be generally neglected.
* the triple contraction of two matrices and a matrix-vector may be reinterpreted in terms of matrix multiplication: $\sum \mathbf{A}_{ij} \mathbf{b}_{j0} \mathbf{C}_{jk} \Longrightarrow \mathbf{A}\, \mbox{diag}(\mathbf{b}) \, \mathbf{C}$.
* repeated indices without summation can be identified as Hadamard products $\mathbf{A}_{ij} \mathbf{B}_{ij} \Longrightarrow \mathbf{A} \circ \mathbf{B}$. This is often the case with ``ArrayDiagonal`` operators.
* generally, open two-paired contraction lines are matrix multiplications:

$$\sum_{j} \mathbf{A}_{ij} \mathbf{B}_{ij} \Longrightarrow \mathbf{A} \mathbf{B}' $$

* while closed two-paired contraction lines are traces:

$$\sum_{ij} \mathbf{A}_{ij} \mathbf{B}_{ij} \Longrightarrow \mbox{tr}\Big(\mathbf{A} \mathbf{B}'\Big) $$

* diagonalization of matrix and matrix-vector can be turned into matrix multiplications

$$A_{ij} b_{j0} \Longrightarrow \mathbf{A} \, \mbox{diag}(\mathbf{b})$$

* diagonalization of two matrix-vectors may be turned into a multiplication by an all-one matrix:

$$a_{i0} b_{j0} \Longrightarrow \mathbf{a}' \cdot \mathbf{1} \cdot \mathbf{b}$$

References
----------

* The Matrix Cookbook, by Kaare Brandt Petersen and Michael Syskind Pedersen, https://www.math.uwaterloo.ca/~hwolkowi/matrixcookbook.pdf
