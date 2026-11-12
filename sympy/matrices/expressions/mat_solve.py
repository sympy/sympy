from sympy.core.mul import Mul
from sympy.core.relational import Equality
from sympy.matrices.matrixbase import MatrixBase
from sympy.matrices.expressions.matadd import MatAdd
from sympy.matrices.expressions.matmul import MatMul
from sympy.matrices.expressions.matexpr import MatrixExpr
from sympy.matrices.expressions.special import Identity, ZeroMatrix


def _is_matrix_factor(factor):
    """Return whether factor is a matrix expression or concrete matrix."""
    return isinstance(factor, (MatrixExpr, MatrixBase))


def _split_factors(term):
    """Separate scalar factors from ordered matrix factors."""
    args = term.args if isinstance(term, MatMul) else (term,)

    scalar_factors = []
    matrix_factors = []

    for factor in args:
        if _is_matrix_factor(factor):
            matrix_factors.append(factor)
        elif factor.is_commutative is True:
            scalar_factors.append(factor)
        else:
            raise NotImplementedError(
                "Unsupported noncommutative scalar factor."
            )

    return Mul(*scalar_factors), matrix_factors


def _classify_term(term, symbol):
    """Classify a variable-dependent term as LEFT, RIGHT, or BARE."""
    if term == symbol:
        return "BARE"

    if not isinstance(term, MatMul):
        raise NotImplementedError(
            "Unsupported expression involving the unknown matrix."
        )

    _, matrix_factors = _split_factors(term)

    positions = [
        i for i, factor in enumerate(matrix_factors)
        if factor == symbol
    ]

    nested = any(
        factor != symbol and factor.has(symbol)
        for factor in matrix_factors
    )

    if nested or len(positions) != 1:
        raise NotImplementedError(
            "The unknown must occur exactly once as a direct matrix factor."
        )

    position = positions[0]

    if len(matrix_factors) == 1:
        return "BARE"
    if position == len(matrix_factors) - 1:
        return "LEFT"
    if position == 0:
        return "RIGHT"

    raise NotImplementedError(
        "Two-sided matrix multiplication is not supported."
    )


def _extract_coefficient(term, symbol, orientation):
    """Extract the coefficient multiplying symbol on the chosen side."""
    scalar, matrix_factors = _split_factors(term)
    position = matrix_factors.index(symbol)

    if orientation == "LEFT":
        factors = matrix_factors[:position]
        identity_size = symbol.rows
    elif orientation == "RIGHT":
        factors = matrix_factors[position + 1:]
        identity_size = symbol.cols
    else:
        raise ValueError("Unknown matrix-equation orientation.")

    if not factors:
        coefficient = Identity(identity_size)
    elif len(factors) == 1:
        coefficient = factors[0]
    else:
        coefficient = MatMul(*factors)

    return scalar * coefficient


def _validate_coefficient(coefficient, symbol, orientation):
    """Check that the coefficient has a compatible square shape."""
    if coefficient.rows != coefficient.cols:
        raise NotImplementedError(
            "The coefficient matrix must be square."
        )

    if orientation == "LEFT":
        if coefficient.cols != symbol.rows:
            raise NotImplementedError(
                "Incompatible dimensions for left-sided matrix solving."
            )
    elif orientation == "RIGHT":
        if coefficient.rows != symbol.cols:
            raise NotImplementedError(
                "Incompatible dimensions for right-sided matrix solving."
            )


def matrix_solve(equation, symbol):
    """Solve a supported one-sided linear matrix equation.

    The unknown must occur exactly once in each variable-dependent term,
    as a direct matrix factor at the same side of its coefficient.

    Parameters
    ==========
    equation : Equality or MatrixExpr
        The matrix equation to solve.
    symbol : MatrixExpr
        The matrix expression to solve for.

    Returns
    =======
    dict
        A dictionary mapping ``symbol`` to its formal solution.

    Notes
    =====
    The coefficient matrix is assumed to be invertible. This function
    constructs a formal inverse; it does not establish invertibility.

    Two-sided equations, nonlinear equations, and rectangular
    coefficient matrices are not supported.
    """
    if not isinstance(symbol, MatrixExpr):
        raise TypeError("The unknown must be a MatrixExpr.")

    if isinstance(equation, Equality):
        expression = equation.lhs - equation.rhs
    else:
        expression = equation

    if not isinstance(expression, MatrixExpr):
        raise TypeError("The equation must be a matrix expression.")

    terms = expression.args if isinstance(expression, MatAdd) else (expression,)

    left_terms = []
    right_terms = []
    bare_terms = []
    constant_terms = []

    for term in terms:
        if not term.has(symbol):
            constant_terms.append(term)
            continue

        orientation = _classify_term(term, symbol)

        if orientation == "LEFT":
            left_terms.append(term)
        elif orientation == "RIGHT":
            right_terms.append(term)
        else:
            bare_terms.append(term)

    if not left_terms and not right_terms and not bare_terms:
        raise ValueError("The equation does not contain the unknown.")

    if left_terms and right_terms:
        raise NotImplementedError(
            "Mixed left- and right-sided matrix equations are not supported."
        )

    if left_terms:
        orientation = "LEFT"
    elif right_terms:
        orientation = "RIGHT"
    else:
        # A standalone multiple of the unknown can be treated
        # as either Identity * symbol or symbol * Identity.
        orientation = "LEFT"

    variable_terms = left_terms + right_terms + bare_terms

    coefficients = [
        _extract_coefficient(term, symbol, orientation)
        for term in variable_terms
    ]

    coefficient = (
        coefficients[0]
        if len(coefficients) == 1
        else MatAdd(*coefficients)
    )

    constant = (
        MatAdd(*constant_terms)
        if constant_terms
        else ZeroMatrix(*symbol.shape)
    )

    _validate_coefficient(coefficient, symbol, orientation)

    if orientation == "LEFT":
        return {symbol: coefficient.I * (-constant)}

    return {symbol: (-constant) * coefficient.I}
