Matrix derivatives
==================

Matrix derivatives are a powerful tool utilized in theoretical statistics and
machine learning.  SymPy provides support for computing closed-formula
derivatives of matrices and, more broadly, symbolic expressions involving
matrices.

Matrix derivatives constitute a complex concept that encompasses both the
principles of matrices and higher-dimensional arrays in general.  The
conventional matrix syntax employed in mathematics, although powerful,
occasionally proves incomplete.

In this explanation, we will introduce the framework of matrix expressions and
array expressions in SymPy, before exploring a detailed description of the
matrix derivation algorithm.

Matrix and array expressions in SymPy
-------------------------------------

SymPy provides matrices and general N-dimensional arrays through the modules
``sympy.matrices`` and ``sympy.tensor.array``, respectively. The basic objects
of these two modules are component-explicit, meaning that every element is
stored and directly accessible, so the object explicitly contains the full set
of component values rather than representing them implicitly.

SymPy also supports symbolic matrix expressions and array expressions through
the respective dedicated submodules ``sympy.matrices.expressions`` and
``sympy.tensor.array.expressions``.  These represent matrices and arrays
abstractly as symbols, defined solely by an identifying name (such as $X$ or
$M$) and a shape tuple.  For matrices, the shape is always a 2-element tuple,
representing the number of rows and columns. The dimensions within this tuple
can be either fixed integers (for known sizes) or arbitrary symbolic
expressions (when dimensions are undetermined).

Unlike component-explicit matrices, which perform elementwise calculations
immediately, symbolic matrix expressions retain operations unevaluated in
expression trees. Arithmetic operations (addition, multiplication, transpose,
etc.) are formally applied but frozen as symbolic representations rather than
computed results. This aligns with conventional mathematical notation in
textbooks, where expressions like ``M*N.T*P`` (also written as
$\mathbf{M}\mathbf{N}^{T} \mathbf{P}$) represent abstract matrix multiplication
rather than explicit numerical results.

Array expressions are analogous, the only difference is that they may have any
number of dimensions, whereas matrices are defined to have exactly two.

In SymPy both matrix symbols and array symbols support Python's ``[ ]``
operator to reference individual elements. For example:

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

The matrix expression module ``sympy.matrices.expressions`` in SymPy reflects
the common convention used in mathematics to represent matrices when indices
are not explicit. Therefore, ``M*N`` or `\mathbf{M}\cdot\mathbf{N}` is the
matrix multiplication, or in index-explicit form $\sum_j M_{ij} N_{jk}$ with
final free indices ${}_{\{ik\}}$.  Common operators such as determinant, trace,
Hadamard (i.e. elementwise) product and Hadamard (i.e. elementwise) power are
provided.

When a function is applied to a matrix, it is understood in the mathematical
sense, so ``exp(M)`` or `\exp(\mathbf{M})` denotes the matrix exponential.
Element-wise operations on matrices are provided by the ``.applyfunc`` method.
For example, to exponentiate every entry of a matrix $\mathbf{M}$ you write
``M.applyfunc(exp)``.  This is denoted as $\exp_{\circ}(\mathbf{M})$, with the
subscript ${}_{\circ}$ flagging an entry-wise application, distinguishing it
from the matrix exponential.

The array expressions module, by contrast, is designed to express arbitrary
operations on arrays of any dimension.

+---------------------------+--------------------------------+-------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| operation name            | SymPy object                   | representation                            | index-explicit equivalent                                                                                 |
+===========================+================================+===========================================+===========================================================================================================+
| tensor product            | ``ArrayTensorProduct``         | $A \boxtimes B$                           | $A_{ij} B_{kl}$                                                                                           |
+---------------------------+--------------------------------+-------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| contraction               | ``ArrayContraction``           | $A$ on axes $a$, $b$                      | $A_{i_{1} i_{2} \ldots i_{a} \ldots i_{b} \ldots } \Rightarrow \sum_j A_{i_{1} \ldots j \ldots j \ldots}$ |
+---------------------------+--------------------------------+-------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| diagonal                  | ``ArrayDiagonal``              | $A$ on axes $a$, $b$                      | $A_{i_{1} i_{2} \ldots i_{a} \ldots i_{b} \ldots } \Rightarrow A_{i_{1} \ldots j \ldots j \ldots}$        |
+---------------------------+--------------------------------+-------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| permutation of dimensions | ``PermuteDims``                | permutation $\sigma$ on axes of $A$       | $A_{i_{1} i_{2} \ldots} \Rightarrow A_{i_{\sigma(1)} i_{\sigma(2) \ldots}}$                               |
+---------------------------+--------------------------------+-------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| addition                  | ``ArrayAdd``                   | elementwise $A + B$                       | $A_{ij} + B_{ij}$                                                                                         |
+---------------------------+--------------------------------+-------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| summation                 | ``ArraySum``                   | summation $\sum_{x} A^{x}$                | $\sum_{x} A^{x}_{ij}$                                                                                     |
+---------------------------+--------------------------------+-------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| elementwise application   | ``ArrayElementwiseApplyFunc``  | $f$ applied to elements of $f_{\circ}(A)$ | $f(A^{x}_{ij})$                                                                                           |
+---------------------------+--------------------------------+-------------------------------------------+-----------------------------------------------------------------------------------------------------------+

In SymPy, array expressions are built from four primitive nodes, most
importantly ``ArrayTensorProduct``, ``ArrayContraction``, ``ArrayDiagonal``,
``PermuteDims``, that compose and manipulate arrays of arbitrary dimension.

These four operators have a canonical form where, when all four are present,
they must be nested in the following sequence (from outermost to innermost):
permutation, diagonalization, contraction, then tensor product.

Furthermore, to reduce ambiguity among mathematically identical expressions,
the contraction and diagonalization index-tuples in the canonical form are
sorted in ascending order according to the minimum index in each tuple.

The canonical order is not automatically applied upon object creation, but can
be attained by explicitly calling ``.doit()``.

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

Matrix expressions can be represented using array expressions and their
associated operators.  However, the converse is not always true, as array
expressions encompass a more general set of operations.

+---------------------------+-------------------------------------------+------------------------------------------------------------+----------------------------------------------------------------------------------+
| matrix operation          | matrix expression form                    | index form                                                 | array expression form                                                            |
+===========================+===========================================+============================================================+==================================================================================+
| matrix multiplication     | $\mathbf{M} \mathbf{N}$                   | ${}_{\{ij\}} \Rightarrow \sum_k M_{ik} N_{kj}$             | contraction: $M \boxtimes N$ on 2nd and 3rd axes                                 |
+---------------------------+-------------------------------------------+------------------------------------------------------------+----------------------------------------------------------------------------------+
| trace                     | $\mbox{tr}(\mathbf{M})$                   | ${}_{\{\}} \Rightarrow \sum_i M_{ii}$                      | contraction: $M$ on 1st and 2nd axes                                             |
+---------------------------+-------------------------------------------+------------------------------------------------------------+----------------------------------------------------------------------------------+
| diagonal                  | $\mbox{diag}(\mathbf{M})$                 | ${}_{\{i\}} \Rightarrow  M_{ii}$                           | diagonalize: $M$ on 1st and 2nd axes                                             |
+---------------------------+-------------------------------------------+------------------------------------------------------------+----------------------------------------------------------------------------------+
| transposition             | $\mathbf{M}^{T}$                          | ${}_{\{ij\}} \Rightarrow M_{ji}$                           | permutation: $M$ on 1st and 2nd axes                                             |
+---------------------------+-------------------------------------------+------------------------------------------------------------+----------------------------------------------------------------------------------+
| Hadamard product          | $\mathbf{M} \circ \mathbf{N}$             | ${}_{\{ij\}} \Rightarrow M_{ij} N_{ij}$                    | diagonalize: $M \boxtimes N$ on 1st-3rd axes and 2nd-4th axes                    |
+---------------------------+-------------------------------------------+------------------------------------------------------------+----------------------------------------------------------------------------------+
| Kronecker product         | $\mathbf{M} \otimes \mathbf{N}$         | ${}_{\{m=id_1+k,n=jd_2+l\}} \Longrightarrow A_{ij} B_{kl}$ | permute $M \boxtimes N$ on 2nd and 3rd axes, then reshape                          |
+---------------------------+-------------------------------------------+------------------------------------------------------------+----------------------------------------------------------------------------------+

Kronecker product versus tensor product
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In SymPy, the Kronecker product $\otimes$ and tensor product $\boxtimes$
represent very similar underlying operations.  Given matrices $\mathbf{A}$ and
$\mathbf{B}$, the elements of the resulting tensor product are combined as the
product of the individual elements:

$$\mathbf{A} \boxtimes \mathbf{B} \Longrightarrow A_{ij} B_{kl} $$

When representing this expression as an array, different approaches arise
depending on the index-order and the reshaping of the 4-dimensional array into
a 2-dimensional array by joining pairs of dimensions.

The array tensor product, as implemented by both ``ArrayTensorProduct`` and the
``tensorproduct`` function, performs a direct nesting of the input arrays. This
operation preserves the index order, typically resulting in a new array with
the combined index sequence:

$$\left[
\begin{matrix}\left[\begin{matrix}{A}_{0,0} {B}_{0,0} & {A}_{0,0} {B}_{0,1}\\{A}_{0,0} {B}_{1,0} & {A}_{0,0} {B}_{1,1}\end{matrix}\right] & \left[\begin{matrix}{A}_{0,1} {B}_{0,0} & {A}_{0,1} {B}_{0,1}\\{A}_{0,1} {B}_{1,0} & {A}_{0,1} {B}_{1,1}\end{matrix}\right]\\\left[\begin{matrix}{A}_{1,0} {B}_{0,0} & {A}_{1,0} {B}_{0,1}\\{A}_{1,0} {B}_{1,0} & {A}_{1,0} {B}_{1,1}\end{matrix}\right] & \left[\begin{matrix}{A}_{1,1} {B}_{0,0} & {A}_{1,1} {B}_{0,1}\\{A}_{1,1} {B}_{1,0} & {A}_{1,1} {B}_{1,1}\end{matrix}\right]\end{matrix}
\right] $$

The Kronecker product can be defined in terms of indices ${}_{\{mn\}}$ where $m$ spans over rows of both $\mathbf{A}$
and $\mathbf{B}$, while $n$ spans over their columns:
$$A \otimes B = \Big[ {}_{\{mn\}} = {}_{\{m=id_1+k,n=jd_2+l\}} \Longrightarrow A_{ij} B_{kl} \Big]$$
where $[d_1 \times d_2]$ is the shape of $\mathbf{A}$.

The Kronecker product, implemented by ``KroneckerProduct`` and ``kronecker_product``, combines the rows and columns of $\mathbf{A}$ and $\mathbf{B}$

$$\left[
\begin{matrix}{A}_{0,0} {B}_{0,0} & {A}_{0,0} {B}_{0,1} & {A}_{0,1} {B}_{0,0} & {A}_{0,1} {B}_{0,1}\\{A}_{0,0} {B}_{1,0} & {A}_{0,0} {B}_{1,1} & {A}_{0,1} {B}_{1,0} & {A}_{0,1} {B}_{1,1}\\{A}_{1,0} {B}_{0,0} & {A}_{1,0} {B}_{0,1} & {A}_{1,1} {B}_{0,0} & {A}_{1,1} {B}_{0,1}\\{A}_{1,0} {B}_{1,0} & {A}_{1,0} {B}_{1,1} & {A}_{1,1} {B}_{1,0} & {A}_{1,1} {B}_{1,1}\end{matrix}
\right]$$

In these visual representations of components, the Kronecker product and tensor
product appear very similar, differing only in their nesting structure.
However, when using a ``Reshape`` operation to convert an array of shape
$[2 \times 2 \times 2 \times 2]$ into one of shape $[4 \times 4]$, the inner-outer
index ordering is preserved.  SymPy arrays use a row-major ordering. That is,
an $N$-dimensional array is internally stored as a one-dimensional array, and a
multi-index is mapped to a single internal index by traversing the rows first.
As a result, the inner $[2\times 2]$ blocks are rearranged into rows in an
order that does not match the layout produced by the Kronecker product.
Converting a tensor product into a Kronecker product by reshaping requires
first permuting the $j$ and $k$ indices:

$$\left[
\begin{matrix}\left[\begin{matrix}{A}_{0,0} {B}_{0,0} & {A}_{0,0} {B}_{0,1}\\{A}_{0,1} {B}_{0,0} & {A}_{0,1} {B}_{0,1}\end{matrix}\right] & \left[\begin{matrix}{A}_{0,0} {B}_{1,0} & {A}_{0,0} {B}_{1,1}\\{A}_{0,1} {B}_{1,0} & {A}_{0,1} {B}_{1,1}\end{matrix}\right]\\\left[\begin{matrix}{A}_{1,0} {B}_{0,0} & {A}_{1,0} {B}_{0,1}\\{A}_{1,1} {B}_{0,0} & {A}_{1,1} {B}_{0,1}\end{matrix}\right] & \left[\begin{matrix}{A}_{1,0} {B}_{1,0} & {A}_{1,0} {B}_{1,1}\\{A}_{1,1} {B}_{1,0} & {A}_{1,1} {B}_{1,1}\end{matrix}\right]\end{matrix}
\right]
$$

This array can be reshaped into a $[4 \times 4]$ Kronecker product of
$\mathbf{A}$ and $\mathbf{B}$.  After reshaping, the ordering of the elements
in the inner blocks matches exactly the row-wise ordering of the Kronecker
product.

>>> from sympy import MatrixSymbol, kronecker_product, tensorproduct
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
element independently. Likewise, differentiating a scalar with respect to a
matrix produces a matrix of the same shape, where each entry is the derivative
of the scalar with respect to the corresponding matrix element.

While the trace and determinant of a matrix are scalars, their derivatives with
respect to a matrix variable can yield nontrivial expressions involving the
matrix itself.

Differentiating a matrix $\mathbf{Y}$ with respect to another matrix
$\mathbf{X}$ results in a 4-dimensional array, where each entry is the
drivative of an element of $\mathbf{Y}$ with respect to an element of
$\mathbf{X}$. This is sometimes informally called a tensor, although in SymPy
the term tensor refers to objects from differential geometry and is handled by
different modules.  Derivatives extend naturally through index notation.  In
index-explicit form, differentiating a matrix expression $Y_{ij}$ with respect
to a matrix $X_{mn}$ yields a four-dimensional array whose entries are given
by:

$$D_{mnij} = \frac{\partial}{\partial X_{mn}} \Big (Y_{ij} \Big ).$$

SymPy adopts a denominator-first index ordering for derivatives, placing the
indices of the differentiation variable ${}_{\{mn\}}$ before those of the
derivand ${}_{\{ij\}}$. This roughly corresponds to the mixed layout convention,
see [LayoutConventions]_.  This ${}_{\{mnij\}}$ convention aligns with Wolfram
Mathematica but differs from PyTorch/NumPy, which implicitly follow the
${}_{\{ijmn\}}$ index-order convention.

This structure also extends to expressions involving mixtures of scalars and
matrices: scalar-by-matrix derivatives produce matrices, matrix-by-scalar
derivatives produce matrices, and matrix-by-matrix derivatives yield
4-dimensional arrays reflecting the relationship between all component pairs.

To explicitly track index reading order, which controls axis transpositions in
multi-dimensional representations, we will use the convention of mapping the
index sequence to the full expression. For example:

$${}_{\{mnij\}} \Longrightarrow\frac{\partial}{\partial X_{mn}} \Big ( A_{ij} \Big )$$

this shows that differentiation indices ($mn$) precede the derivand indices
($ij$).  Importantly, this notation inherently captures index permutations,
eliminating the need for operators such as transposition. For example:

$${}_{\{ij\}} \Longrightarrow \Big(M^{T}\Big)_{ij} = M_{ji}$$

demonstrates how transposition is encoded solely through index-order
reversal, the ${}_{\{ij\}}$ mapping of $M_{ji}$ directly yields the
component expression of the transposition.

Exceptions to the derivative index-order
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Be aware that certain functions within the matrix module, such as
``.jacobian(...)``, employ their own index ordering conventions, which are
typically derivand-first.  Specifically, in the case of the Jacobian function,
the derivand axis corresponds to the columns, while the deriving variable axis
corresponds to the rows of the resulting Jacobian matrix.

Indeed, let's consider the following example:
$$X = \left[\begin{matrix}f{\left(x,y \right)}\\g{\left(x,y \right)}\end{matrix}\right]$$
and
$$Y = \left[\begin{matrix}x\\y\end{matrix}\right]$$

Using the matrix derivative convention, ``X.diff(Y).reshape(2, 2)`` returns:

$$\left[\begin{matrix}\frac{\partial}{\partial x} f{\left(x,y \right)} & \frac{\partial}{\partial x} g{\left(x,y \right)}\\\frac{\partial}{\partial y} f{\left(x,y \right)} & \frac{\partial}{\partial y} g{\left(x,y \right)}\end{matrix}\right]$$

In contrast, calling ``X.jacobian(Y)`` returns the previous matrix transposed:

$$\left[\begin{matrix}\frac{\partial}{\partial x} f{\left(x,y \right)} & \frac{\partial}{\partial y} f{\left(x,y \right)}\\\frac{\partial}{\partial x} g{\left(x,y \right)} & \frac{\partial}{\partial y} g{\left(x,y \right)}\end{matrix}\right]$$

Reconstructing a matrix expression
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In SymPy, differentiating a matrix expression with respect to a matrix symbol
normally yields a 4-dimensional array.  When one or more of its axes happen to
be singleton (dimension 1), the tensor collapses to a lower dimensional object.
SymPy automatically detects this and returns the simplest equivalent matrix
expression instead of keeping the full 4-dimensional form.

For example, let $\mathbf{x}$ be a matrix-vector of shape $[k \times 1]$, that
is a vector of size $[k]$ disguised as a single column matrix by adding a
singleton dimension to its shape.  The derivative of $\mathbf{x}$ with respect
to $\mathbf{x}$ is then the 4-dimensional identity array of shape $[k \times 1
\times k \times 1]$, but because the second axis of both the numerator and the
denominator has length 1, the result can be collapsed to the $[k \times k]$
identity matrix $\mathbf{I}_{[k]}$ instead of keeping the full tensor form.
Explicitly:

$$\left[{}_{\{mnij\}} \Longrightarrow \frac{\partial}{\partial x_{mn}}x_{ij} = \delta_{mi} \delta_{nj} = \big( \mathbf{I}_{[k]} \boxtimes \mathbf{I}_{[1]} \big)_{minj}\right] \rightarrow \Big[ {}_{\{mi\}}\Longrightarrow \big(\mathbf{I}_{[k]}\big)_{mi}\Big]$$

Here, the second and fourth axes are singletons, appearing in the derivative as
scalar identity matrix $\mathbf{I}_{[1]}$, equivalent to scalar unit.  In index
notation they only survive in a Kronecker $\delta_{nj}$.

Likewise, $\frac{\partial}{\partial \mathbf{x}}\mathbf{x}^{T}$ has shape
$[k \times 1 \times 1 \times k]$, and collapses again to the same
$[k \times k]$ identity matrix $\mathbf{I}_{[k]}$.

Derive matrix by itself
~~~~~~~~~~~~~~~~~~~~~~~

As an example, the derivative of $[k \times l]$ matrix $\mathbf{X}$ by itself
is given by identity relationships between indices:

$${}_{\{mnij\}} \Longrightarrow \frac{\partial}{\partial X_{mn}} \Big ( X_{ij} \Big ) = \delta_{mi} \delta_{nj} = \Big( \mathbf{I}_{[k]} \boxtimes \mathbf{I}_{[l]} \Big)_{minj}, $$

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

The permutation `(3)(1 2)` maps ${}_{\{mnij\}}$ into ${}_{\{minj\}}$. Note that
SymPy uses zero-offset indexing and the trivial cycle `(3)` is included simply
to indicate the maximum size of the permutation. The ``PermuteDims`` object in
SymPy interprets a permutation as mapping the new index order to the old index
order (i.e., the index order of the wrapped expression). This convention is
aligns with most scientific Python libraries but is the opposite of the
convention used by Wolfram Mathematica.

Notice that the result of ``X.diff(X)`` can also be expressed as

``Reshape(KroneckerProduct(I, I), (k, l, k, l))``,

indeed some authors represent

$$\frac{\partial}{\partial \mathbf{X}}\mathbf{X} = \mathbf{I}_{[k]} \bar\otimes \mathbf{I}_{[l]},$$

where $\bar\otimes$ represents the 4-dimensionally reshaped Kronecker
product.

Mixing symbols and elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~

SymPy distinguishes between matrix symbols ``M`` (or more generally matrix
expressions), an algebraic object whose indices are left implicit, and the
explicit scalar entries ``M[i, j]`` that carry visible indices.  While most of
SymPy focuses on derivatives of matrix expressions, it can also handle
expressions that freely mix matrix-level and index-level notation.  For
example:

$$\frac{\partial \mathbf{M}}{\partial M_{ij}} = \mathbf{E}_{ij}$$

where $\mathbf{E}_{ij}$ is the matrix unit, a matrix with only one non-zero
element with value 1 at position ${}_{\{ij\}}$.

Likewise,

$$\frac{\partial M_{ij}}{\partial \mathbf{M}} = \mathbf{E}_{ij}$$

Deriving an elementwise function of matrices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A common example in neural networks involves deriving expressions for matrices
that undergo elementwise functions. For instance, consider a matrix $\mathbf{W}$
of shape $[k \times l]$ and a matrix-vector $\mathbf{x}$ of shape $[l \times 1]$

$$\frac{\partial}{\partial \mathbf{x}}{\operatorname{atanh}}_{\circ}\left({\mathbf{W} \mathbf{x}}\right) = \mathbf{W}^{T} \operatorname{diag}\left({\left( d \mapsto \frac{1}{1 - d^{2}} \right)}_{\circ}\left({\mathbf{W} \mathbf{x}}\right)\right)$$

The symbol ${}_{\circ}$ is used to mark the function as acting elementwise.
Here, we have used the **diag** operator:

$$\mbox{diag}\left(\left[\begin{matrix}{b}_{0}\\{b}_{1}\\{b}_{2}\end{matrix}\right] \right) = \left[\begin{matrix}{b}_{0} & 0 & 0\\0 & {b}_{1} & 0\\0 & 0 & {b}_{2}\end{matrix}\right]$$

How matrix derivatives work
---------------------------

The matrix derivation algorithm takes a matrix expression and converts it into
an equivalent array expression using array operators.  From there,
differentiation becomes straightforward, using a variation of the chain rule
that also keeps track of free indices and their ordering at each step.

After computing the derivative, the algorithm attempts to reconstruct an
equivalent matrix expression, removing any singleton or unnecessary diagonal
dimensions.  If no suitable matrix expression has been found, the derivative is
returned as an array expression.

The symbolic derivation algorithm produces an array expression equivalent to
the derivative, effectively altering the operation-flow graph.  Unlike other
platforms optimized for numeric computations, which compute derivatives without
modifying the overall operation graph, this approach finds a symbolic form of
the derivative.

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

+--------------------------------+---------------------------------+--------------------------------------------------------------------------------------+
| operation                      | expression                      | chain rule                                                                           |
+================================+=================================+======================================================================================+
| matrix addition                | $\mathbf{Y} + \mathbf{Z}$       | $\partial \mathbf{Y} + \partial\mathbf{Z}$                                           |
+--------------------------------+---------------------------------+--------------------------------------------------------------------------------------+
| matrix multiplication          | $\mathbf{Y}\mathbf{Z}$          | $(\partial \mathbf{Y})\mathbf{Z} + \mathbf{Y} (\partial\mathbf{Z})$                  |
+--------------------------------+---------------------------------+--------------------------------------------------------------------------------------+
| Hadamard (elementwise) product | $\mathbf{Y} \circ \mathbf{Z}$   | $(\partial \mathbf{Y})\circ\mathbf{Z} + \mathbf{Y}\circ(\partial\mathbf{Z})$         |
+--------------------------------+---------------------------------+--------------------------------------------------------------------------------------+
| tensor product                 | $\mathbf{Y}\boxtimes\mathbf{Z}$ | $(\partial \mathbf{Y})\boxtimes\mathbf{Z} + \mathbf{Y}\boxtimes(\partial\mathbf{Z})$ |
+--------------------------------+---------------------------------+--------------------------------------------------------------------------------------+
| inverse                        | $\mathbf{Y}^{-1}$               | $-\mathbf{Y}^{-1} (\partial \mathbf{Y}) \mathbf{Y}^{-1}$                             |
+--------------------------------+---------------------------------+--------------------------------------------------------------------------------------+
| trace                          | $\mbox{tr}(\mathbf{Y})$         | $\mbox{tr}(\partial\mathbf{Y})$                                                      |
+--------------------------------+---------------------------------+--------------------------------------------------------------------------------------+
| determinant                    | $\mbox{det}(\mathbf{Y})$        | $\mbox{det}(\mathbf{Y}) \mbox{tr}(\mathbf{Y}^{-1}\partial\mathbf{Y})$                |
+--------------------------------+---------------------------------+--------------------------------------------------------------------------------------+
| transposition                  | $\mathbf{Y}^{T}$                | $(\partial\mathbf{Y})^{T}$                                                           |
+--------------------------------+---------------------------------+--------------------------------------------------------------------------------------+

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
$$-\mathbf{Y}^{-1} \boxtimes \frac{\partial\mathbf{Y}}{\partial\mathbf{X}} \boxtimes \mathbf{Y}^{-1}$$
which is then contracted on the 2-nd and 5-th axes, i.e. (1, 4), and on the 6-th and 7-th axes, i.e. (5, 6).
$$-\frac{1}{\partial\mathbf{X}} \Big( \mathbf{Y}^{-1} (\partial\mathbf{Y}) \mathbf{Y}^{-1} \Big)$$
The contraction reproduce the structure of the matrix multiplication of the chain rule for the inverse.
As a last step, we need to take the 2-nd and 3-rd axes of the contracted expression in front of the other ones,
as they are the axes referring to the deriving variable $\mathbf{X}$.

A **concrete example** involving the inverse matrix taken from [MatrixCookbook]_ (number 124):

$$\frac{\partial}{\partial \mathbf{X}} \operatorname{tr}\left(\mathbf{A} \mathbf{X}^{-1} \mathbf{B} \right) = - \left(\mathbf{X}^{T}\right)^{-1} \mathbf{A}^{T} \mathbf{B}^{T} \left(\mathbf{X}^{T}\right)^{-1} $$

>>> from sympy import MatrixSymbol, Trace, Inverse
>>> from sympy.abc import k
>>> A = MatrixSymbol("A", k, k)
>>> B = MatrixSymbol("B", k, k)
>>> X = MatrixSymbol("X", k, k)
>>> expr = Trace(A*X**(-1)*B)
>>> expr.diff(X)
-X.T**(-1)*A.T*B.T*X.T**(-1)


Example derivatives of the determinant
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some examples of the derivative of the determinant:

[MatrixCookbook]_ example 49:

$$\frac{\partial}{\partial \mathbf{X}} \mbox{det}(\mathbf{X}) = \mbox{det}(\mathbf{X}) \left(\mathbf{X}^{T}\right)^{-1}$$

>>> from sympy import Determinant
>>> expr = Determinant(X)
>>> expr.diff(X)
Determinant(X)*X.T**(-1)

[MatrixCookbook]_ example 51:

$$\frac{\partial}{\partial \mathbf{X}} \mbox{det}(\mathbf{A}\mathbf{X}\mathbf{B}) = \mbox{det}(\mathbf{A}\mathbf{X}\mathbf{B}) \left(\mathbf{X}^{T}\right)^{-1}$$

>>> expr = Determinant(A*X*B)
>>> expr.diff(X)
(Determinant(A)*Determinant(B)*Determinant(X))*X.T**(-1)

[MatrixCookbook]_ example 55:

$$\frac{\partial}{\partial \mathbf{X}} \log\,\mbox{det}(\mathbf{X}) = 2 \left(\mathbf{X}^{T}\right)^{-1}$$

>>> from sympy import log
>>> expr = log(Determinant(X.T*X))
>>> expr.diff(X)
2*X.T**(-1)


Example derivatives of the trace
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[MatrixCookbook]_ example 107:

$$\frac{\partial}{\partial \mathbf{X}} \mbox{tr}(\mathbf{X}^2\mathbf{B}) = \left( \mathbf{X} \mathbf{B} + \mathbf{B}\mathbf{X}\right)^{T}$$

>>> expr = Trace(X**2*B)
>>> expr.diff(X)
B.T*X.T + X.T*B.T

[MatrixCookbook]_ example 108:

$$\frac{\partial}{\partial \mathbf{X}} \mbox{tr}(\mathbf{X}^{T}\mathbf{B}\mathbf{X}) = \mathbf{B} \mathbf{X} + \mathbf{B}^{T} \mathbf{X}$$

>>> expr = Trace(X.T*B*X)
>>> expr.diff(X)
B*X + B.T*X

Array expression to matrix expression conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The hardest part of the algorithm of matrix derivation is converting the derivative back to a matrix expression.
The algorithm returns an array expression whose axes have been contracted, diagonalized and permuted.
To recover the matrix form, we must recognize the equivalent matrix operations those operations represent.

Matrix derivation may produce array expressions that lack a closed-form matrix counterpart,
so the result cannot always be rewritten in matrix form.

When the contractions line up in clear pairs that mirror a sequential pattern,
spotting the equivalent matrix multiplication or trace forms is straightforward.
For example:

$$\sum_{b,c,e,f} M_{cb} N_{cd} P_{ef} Q_{ba} R_{fe}
\Longrightarrow \mathbf{Q}^{T} \mathbf{M}^{T} \mathbf{N} \, \mbox{tr}\big( \mathbf{P} \mathbf{R} \big)
$$

* generally, open two-paired contraction lines are matrix multiplications:
  $$\sum_{j} \mathbf{A}_{ij} \mathbf{B}_{ij} \Longrightarrow \mathbf{A} \mathbf{B}^{T} $$
* while closed two-paired contraction lines are traces:
  $$\sum_{ij} \mathbf{A}_{ij} \mathbf{B}_{ij} \Longrightarrow \mbox{tr}\Big(\mathbf{A} \mathbf{B}^{T}\Big) $$

Matrix multiplication is the most frequent pattern, but traces, diag-expansions, and Hadamard products may also appear.

* repeated indices without summation can be identified as Hadamard products
  $$\mathbf{A}_{ij} \mathbf{B}_{ij} \Longrightarrow \mathbf{A} \circ \mathbf{B}.$$
  This is often the case with ``ArrayDiagonal`` operators.
* double single-index contraction on both axes of the same matrix or matrix expression, suppose that matrix $\mathbf{M}$ has shape $[m \times n]$, then:
  $$\sum_{i,j} M_{ij} \Longrightarrow \mathbf{1}_{[1\times m]} \, \mathbf{M} \, \mathbf{1}_{[n\times 1]}$$

Conversions from array to matrix with dimensional reduction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simplification step tries to drop any singleton or diagonal dimensions;
if this reduces the array to two dimensions, an equivalent matrix expression is sought.
A handful of tricks is used to collapse the array dimensions:

* tensor product of two matrix-vectors $\mathbf{a}$ and $\mathbf{a}$, of shape $[k \times 1]$ each, $\mathbf{a} \boxtimes \mathbf{b}$, can be turned into a matrix multiplication over their trivial dimension: $\mathbf{a} \cdot \mathbf{b}^{T}$. Indeed
  $$\mathbf{a} \boxtimes \mathbf{b} = a_{i0} b_{j0} \Longrightarrow \sum_{k=0}^0 a_{ik} b_{jk} = \mathbf{a} \mathbf{b}^{T} $$
  This has shrunk the array from $[k \times 1 \times k \times 1]$ to $[k \times k]$ by squeezing out the singleton dimensions.
* similar to the previous point, but more complex:
  $$\mathbf{a} \boxtimes \mathbf{b} \boxtimes \mathbf{x}^{T}\mathbf{x} \Longrightarrow \mathbf{a} \mathbf{x}^{T} \mathbf{x} \mathbf{b}^{T} $$
* if one term of the tensor product has shape $[1 \times 1]$, the result may be expressed as a Kronecker product $\otimes$ between matrices:
  $$\mathbf{x}^{T} \mathbf{x} \boxtimes \mathbf{A} \Longrightarrow \left( \mathbf{x}^{T} \mathbf{x} \right) \otimes \mathbf{A} $$
  This has shrunk the array from $[1 \times 1 \times k \times k]$ to $[k \times k]$.
* the triple contraction of two matrices and a matrix-vector may be reinterpreted in terms of matrix multiplication:
  $$\sum_{j} \mathbf{A}_{ij} \mathbf{b}_{j0} \mathbf{C}_{jk} \Longrightarrow \mathbf{A}\, \mbox{diag}(\mathbf{b}) \, \mathbf{C}.$$
  This has shrunk $[k\times 1 \times k]$ to $[k \times k]$.
* diagonalization of matrix and matrix-vector can be turned into matrix multiplications
  $$A_{ij} b_{j0} = \sum_{k} A_{ik} \, \mbox{diag}(\mathbf{b})_{kj} \Longrightarrow \mathbf{A} \, \mbox{diag}(\mathbf{b})$$
  This has shrunk $[k \times k \times 1]$ to $[k \times k]$.
* multi-index contractions that span up to two matrices and any number of vectors reduce to ordinary matrix–vector products once the vectors are promoted to diagonal matrices via **diag( )**:
  $$\sum_{j} A_{ij} b_{j0} c_{j0} d_{j0} E_{jk} \Longrightarrow \mathbf{A} \, \mbox{diag}(\mathbf{b}) \, \mbox{diag}(\mathbf{c}) \, \mbox{diag}(\mathbf{d}) \, \mathbf{E} $$
  Notice that an arbitrary number of matrix-vectors can be added, each contributes to the final expression with a singleton dimension,
  which can be dropped.
* multi-index contractions involving identity matrices may create redundant diagonal axes that can be simplified away:
  $$\sum_{l} A_{il} I_{lk} I_{jl} \Longrightarrow \mathbf{A}$$
  by summing up, this is equivalent to:
  $${}_{\{ijk\}} \Rightarrow A_{ik} I_{kj}$$
  or in component-explicit form:
  $$\left[\begin{matrix}\left[\begin{matrix}{A}_{00} & 0 & 0\\0 & {A}_{01} & 0\\0 & 0 & {A}_{02}\end{matrix}\right] & \left[\begin{matrix}{A}_{10} & 0 & 0\\0 & {A}_{11} & 0\\0 & 0 & {A}_{12}\end{matrix}\right] & \left[\begin{matrix}{A}_{20} & 0 & 0\\0 & {A}_{21} & 0\\0 & 0 & {A}_{22}\end{matrix}\right]\end{matrix}\right]\Longrightarrow \left[\begin{matrix}{A}_{00} & {A}_{01} & {A}_{02}\\{A}_{10} & {A}_{11} & {A}_{12}\\{A}_{20} & {A}_{21} & {A}_{22}\end{matrix}\right]$$
  This collapse has shrunk a **non-singleton dimension**, from $[k \times k \times k]$ to $[k \times k]$.
  Unlike squeezing singleton dimensions, this simplification cannot be done by a simple array reshaping.
  Similarly,
  $$\sum_{l} A_{lj} I_{il} I_{lk} \Longrightarrow \mathbf{A}^{T}$$
  summing up reduces it to
  $${}_{\{ijk\}} \Rightarrow A_{ij} I_{ik} $$
  which is equivalent, in explicit form, to:
  $$\left[\begin{matrix}\left[\begin{matrix}{A}_{00} & 0 & 0\\{A}_{01} & 0 & 0\\{A}_{02} & 0 & 0\end{matrix}\right] & \left[\begin{matrix}0 & {A}_{10} & 0\\0 & {A}_{11} & 0\\0 & {A}_{12} & 0\end{matrix}\right] & \left[\begin{matrix}0 & 0 & {A}_{20}\\0 & 0 & {A}_{21}\\0 & 0 & {A}_{22}\end{matrix}\right]\end{matrix}\right] \Longrightarrow \left[\begin{matrix}{A}_{00} & {A}_{10} & {A}_{20}\\{A}_{01} & {A}_{11} & {A}_{21}\\{A}_{02} & {A}_{12} & {A}_{22}\end{matrix}\right] $$
  this last expression has three axes, but the second ${}_{\{j\}}$ and third ${}_{\{k\}}$ are diagonal:
  with the remaining indices held fixed, they form a diagonal matrix,
  so the expression is nonzero only when $j = k$.
  One of those axes can therefore be dropped without loss of information.

Applications
------------

Linear regression
~~~~~~~~~~~~~~~~~

In linear regression, we observe a vector of responses $\mathbf{y}$ of length
$n$ and, for each of these $n$ observations, $p$ covariates (also called
predictors) whose values are assumed to explain or predict the responses.  The
covariates are collected in a design matrix $\mathbf{X}$ of size $n \times p$.

An intercept term is often included in the model, this can be represented by
adding a column of ones to the matrix $\mathbf{X}$.

The goal is to find a vector of regression coefficients $\boldsymbol{\beta}$ of
length $p$ such that the linear predictor $\mathbf{X}\boldsymbol{\beta}$
approximates the observed responses $\mathbf{y}$ as closely as possible.

In ordinary least squares (OLS) regression, the coefficients
$\boldsymbol{\beta}$ are chosen to minimize the objective function, in this
case the squared Euclidean $\ell_2$ norm of the residual vector:

$$\min_{\boldsymbol{\beta}} \Big\Vert \mathbf{y} - \mathbf{X} \boldsymbol{\beta} \Big\Vert^2$$

>>> from sympy.abc import n, p
>>> X = MatrixSymbol("X", n, p)
>>> y = MatrixSymbol("y", n, 1)
>>> beta = MatrixSymbol("beta", p, 1)
>>> obj = (y - X*beta).T * (y - X*beta)

The variable ``obj`` is the objective function to be minimized, here
represented as a dot-product,

$$\left(\mathbf{y} - \mathbf{X}\boldsymbol{\beta}\right)^{T} \left(\mathbf{y} - \mathbf{X}\boldsymbol{\beta}\right).$$

A standard way to find the minimizer is to differentiate this objective
function with respect to $\boldsymbol{\beta}$ and set the resulting gradient
equal to zero:

>>> obj.diff(beta)
-2*X.T*(-X*beta + y)

that is, computing the derivative of the objective function with respect to
$\boldsymbol{\beta}$ yields

$$- 2 \mathbf{X}^{T} \left(\mathbf{y} - \mathbf{X}\boldsymbol{\beta}\right).$$

Setting this derivative equal to zero gives the normal equations, which can
then be solved to obtain the ordinary least squares estimator.  We can
rearrange:

>>> ((X.T*X).inv() * obj.diff(beta)).expand()
-2*(X.T*X)**(-1)*X.T*y + 2*beta

that is, $- 2 \left(\mathbf{X}^{T} \mathbf{X}\right)^{-1} \mathbf{X}^{T} \mathbf{y} + 2 \boldsymbol{\beta}$.

This leads immediately to the well-known closed-form solution for the ordinary
least squares estimator:

$$\boldsymbol{\beta} = \left(\mathbf{X}^{T} \mathbf{X}\right)^{-1} \mathbf{X}^{T} \mathbf{y}$$

Ridge regression
~~~~~~~~~~~~~~~~

Ridge regression is a variation of linear regression in which the objective
function includes an additional term proportional to the squared norm of
$\boldsymbol{\beta}$:

$$\min_{\boldsymbol{\beta}} \Big( \big \Vert \mathbf{X} \boldsymbol{\beta} - \mathbf{y} \big \Vert^2 + \lambda \big \Vert \boldsymbol{\beta} \big \Vert^2 \Big)$$

This approach causes the components of $\boldsymbol{\beta}$ to remain small and
provides better behavior in the presence of multicollinearity, that is, when
the columns of $\mathbf{X}$ are correlated, which is a well-known weakness of
ordinary linear regression.

Using SymPy, the objective function can be written as:

>>> from sympy import symbols
>>> lamda = symbols("lambda")
>>> residual = y - X*beta
>>> obj = residual.T * residual + lamda * beta.T * beta
>>> obj
lambda*beta.T*beta + (-beta.T*X.T + y.T)*(-X*beta + y)

which corresponds to $\lambda \boldsymbol{\beta}^{T} \boldsymbol{\beta} + \left( \mathbf{y}^{T} - \boldsymbol{\beta}^{T} \mathbf{X}^{T} \right) \left( \mathbf{y} -
\mathbf{X} \boldsymbol{\beta} \right)$.

Taking the derivative of the objective function with respect to
$\boldsymbol{\beta}$ yields:

>>> obj.diff(beta)
(2*lambda)*beta - 2*X.T*(-X*beta + y)

which is $2 \lambda \boldsymbol{\beta} - 2 \mathbf{X}^{T} \left( \mathbf{y} - \mathbf{X} \boldsymbol{\beta} \right)$.
Regrouping the terms by $\boldsymbol{\beta}$ and dropping the constant factor
of 2, we obtain:

$$\left(\lambda \mathbf{I} + \mathbf{X}^{T} \mathbf{X}\right  ) \boldsymbol{\beta} - \mathbf{X}^{T} \mathbf{y}$$

Solving for the stationary point by setting this derivative equal to zero gives
the normal equations for ridge regression:

$$\boldsymbol{\beta} = \left(\lambda \mathbf{I} + \mathbf{X}^{T} \mathbf{X}\right)^{-1} \mathbf{X}^{T} \mathbf{y}$$

Principal component analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Given a matrix $\mathbf{X}$ whose columns are zero-mean, we seek a column
vector $\mathbf{w}$ such that the projection $\mathbf{X}\mathbf{w}$ has maximum
variance.

Since the columns of $\mathbf{X}$ have zero mean, the vector
$\mathbf{X}\mathbf{w}$ will also have zero mean, as it is simply a weighted
linear combination of the columns of $\mathbf{X}$. Consequently, the variance
of $\mathbf{X}\mathbf{w}$ is proportional to its squared norm. This observation
leads to the following optimization problem:

$$\underset{\mathbf{w}}{\operatorname{arg\,max}} \, \Big\Vert \mathbf{X} \mathbf{w} \Big\Vert^2 \\ \mbox{with } \big \Vert \mathbf{w} \big \Vert = 1$$

Expressing the norms in terms of matrix products, the problem can be written as

$$\underset{\mathbf{w}}{\operatorname{arg\,max}} \left( \mathbf{w}^{T} \mathbf{X}^{T} \mathbf{X} \mathbf{w} \right) \\ \mbox{with } \mathbf{w}^{T}\mathbf{w} = 1$$

To solve this constrained optimization problem, we introduce a Lagrange
multiplier and define the Lagrangian

$$\mathcal{L} = \mathbf{w}^{T} \mathbf{X}^{T} \mathbf{X} \mathbf{w} + \lambda ( \mathbf{w}^{T}\mathbf{w} - 1 )$$

In SymPy, this setup can be expressed as follows:

>>> from sympy import Identity
>>> lamda, n, p = symbols("lambda n p")
>>> X = MatrixSymbol("X", n, p)
>>> w = MatrixSymbol("w", p, 1)
>>> target = X*w
>>> lagr_mult = target.T*target + lamda*(w.T*w - Identity(1))
>>> lagr_mult
lambda*(-I + w.T*w) + w.T*X.T*X*w

We now proceed by computing the partial derivatives of the Lagrangian with
respect to $\mathbf{w}$ and $\lambda$:

>>> d_w = lagr_mult.diff(w)
>>> d_w
(2*lambda)*w + 2*X.T*X*w
>>> d_lambda = lagr_mult.diff(lamda)
>>> d_lambda
-I + w.T*w

and setting them equal to zero in order to obtain the stationary conditions for
this optimization problem.  That is, we are required to solve the following two
equations:

$$\frac{\partial \mathcal{L}}{\partial \mathbf{w}} = 2 \lambda \mathbf{w} + 2 \mathbf{X}^{T} \mathbf{X} \mathbf{w} = 0$$

$$\frac{\partial \mathcal{L}}{\partial \lambda} = \mathbf{w}^{T} \mathbf{w} - \mathbb{I} = 0$$

The first equation corresponds to the eigenvalue equation of the matrix
$\mathbf{X}^{T} \mathbf{X}$. Therefore, we conclude that $\mathbf{w}$ must be an
eigenvector of $\mathbf{X}^{T} \mathbf{X}$.

References
----------

.. [MatrixCookbook] The Matrix Cookbook, by Kaare Brandt Petersen and Michael Syskind Pedersen, https://www.math.uwaterloo.ca/~hwolkowi/matrixcookbook.pdf
.. [LayoutConventions] https://en.wikipedia.org/wiki/Matrix_calculus#Layout_conventions
