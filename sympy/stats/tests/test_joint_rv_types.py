from sympy import Matrix
from sympy.stats import MultivariateNormal, density

def test_singleton_matrix_pdf_equals_scalar():
    M = Matrix([[1]])
    X = MultivariateNormal('X', M, M)
    pdf = density(X)

    assert pdf(Matrix([[1]])) == pdf(1)
