Matrix References
=================

Matrix User Classes
-------------------

.. currentmodule:: sympy.matrices

.. autosummary::
    :toctree: _generated/
    :template: matrix_template.rst

    ~dense.MutableDenseMatrix
    ~immutable.ImmutableDenseMatrix
    ~sparse.MutableSparseMatrix
    ~immutable.ImmutableSparseMatrix

Matrix Aliases
--------------

.. currentmodule:: sympy

.. autosummary::
    :toctree: _generated/

    Matrix
    MutableMatrix
    ImmutableMatrix
    SparseMatrix

Matrix Base Classes
-------------------

.. currentmodule:: sympy.matrices

.. autosummary::
    :toctree: _generated/
    :template: matrix_template.rst

    ~common.MatrixShaping
    ~common.MatrixSpecial
    ~common.MatrixProperties
    ~common.MatrixOperations
    ~common.MatrixArithmetic

    ~matrices.MatrixDeterminant
    ~matrices.MatrixReductions
    ~matrices.MatrixSubspaces
    ~matrices.MatrixEigen
    ~matrices.MatrixCalculus
    ~matrices.MatrixBase

    ~dense.DenseMatrix
    ~sparse.SparseMatrix

Matrix Exceptions
-----------------

.. autosummary::
    :toctree: _generated/

    ~common.MatrixError
    ~common.ShapeError
    ~common.NonSquareMatrixError
    ~common.NonInvertibleMatrixError
    ~common.NonPositiveDefiniteMatrixError

Matrix Helper Functions
-----------------------

.. autosummary::
    :toctree: _generated/

    ~dense.zeros
    ~dense.ones
    ~dense.eye
    ~dense.diag
    ~dense.jordan_cell
    ~dense.randMatrix

    ~dense.rot_axis1
    ~dense.rot_axis2
    ~dense.rot_axis3

    ~dense.matrix_multiply_elementwise

    ~dense.hessian
    ~dense.GramSchmidt
    ~dense.wronskian
    ~dense.casoratian

Numpy Utility Functions
-----------------------

.. autosummary::
    :toctree: _generated/

    ~dense.list2numpy
    ~dense.matrix2numpy
    ~dense.symarray
