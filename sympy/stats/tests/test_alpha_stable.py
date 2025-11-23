from sympy import Symbol, sqrt, pi, exp, Rational, oo, simplify
from sympy.testing.pytest import raises
from sympy.stats import AlphaStable, density, characteristic_function
def test_alpha_stable_cauchy():
    """Test Cauchy special case (alpha=1, beta=0)"""
    x = Symbol('x', real=True)
    X = AlphaStable('X', 1, 0, 1, 0)
    assert density(X)(x) == 1 / (pi * (x ** 2 + 1))


def test_alpha_stable_normal():
    """Test Normal special case (alpha=2)"""
    x = Symbol('x', real=True)
    Y = AlphaStable('Y', 2, 0, 1, 0)

    # Check PDF
    pdf = density(Y)(x)
    expected = sqrt(2) * exp(-x ** 2 / 2) / (2 * sqrt(pi))
    assert pdf == expected


def test_alpha_stable_levy():
    """Test Levy special case (alpha=0.5, beta=1)"""
    x = Symbol('x', real=True, positive=True)
    L = AlphaStable('L', Rational(1, 2), 1, 1, 0)

    # PDF should be defined for x > 0
    pdf = density(L)(x)
    assert pdf != 0


def test_alpha_stable_characteristic_function():
    """Test characteristic function"""
    t = Symbol('t', real=True)
    X = AlphaStable('X', 1.5, 0.5, 2, 1)

    # Should not raise an error
    cf = characteristic_function(X)(t)
    assert cf is not None


def test_alpha_stable_general_no_pdf():
    """Test that general case raises NotImplementedError for PDF"""
    x = Symbol('x', real=True)
    Z = AlphaStable('Z', 1.5, 0.5, 1, 0)

    raises(NotImplementedError, lambda: density(Z)(x))


def test_alpha_stable_parameters():
    """Test parameter handling"""
    X = AlphaStable('X', 1.8, 0.3, 2, 5)

    assert X.pspace.distribution.alpha == 1.8
    assert X.pspace.distribution.beta == 0.3
    assert X.pspace.distribution.scale == 2
    assert X.pspace.distribution.location == 5


def test_alpha_stable_support():
    """Test that support is (-oo, oo)"""
    X = AlphaStable('X', 1.5, 0, 1, 0)
    assert X.pspace.domain.set.left == -oo
    assert X.pspace.domain.set.right == oo


def test_alpha_stable_scaled_cauchy():
    """Test Cauchy with different scale and location"""
    x = Symbol('x', real=True)
    X = AlphaStable('X', 1, 0, 2, 3)

    # Cauchy with scale=2, location=3
    pdf = density(X)(x)
    expected = 1 / (2 * pi * (1 + ((x - 3) / 2) ** 2))
    assert simplify(pdf - expected) == 0


def test_alpha_stable_scaled_normal():
    """Test Normal with different scale and location"""
    x = Symbol('x', real=True)
    Y = AlphaStable('Y', 2, 0, 2, 1)

    # Normal with scale=2, location=1
    # Variance = 2*scale^2 = 8
    pdf = density(Y)(x)
    z = (x - 1) / 2
    expected = exp(-z ** 2 / 2) / (2 * sqrt(2 * pi))
    assert simplify(pdf - expected) == 0


def test_alpha_stable_symmetric():
    """Test symmetric case (beta=0)"""
    X = AlphaStable('X', 1.5, 0, 1, 0)

    # Should be able to create it
    assert X.pspace.distribution.beta == 0


def test_alpha_stable_skewed():
    """Test skewed case"""
    X = AlphaStable('X', 1.5, 0.8, 1, 0)

    assert X.pspace.distribution.beta == 0.8


def test_alpha_stable_parameter_validation():
    """Test that invalid parameters raise errors"""
    # alpha too small
    raises(ValueError, lambda: AlphaStable('X', 0, 0, 1, 0))

    # alpha too large
    raises(ValueError, lambda: AlphaStable('X', 2.5, 0, 1, 0))

    # beta too small
    raises(ValueError, lambda: AlphaStable('X', 1.5, -1.5, 1, 0))

    # beta too large
    raises(ValueError, lambda: AlphaStable('X', 1.5, 1.5, 1, 0))

    # negative scale
    raises(ValueError, lambda: AlphaStable('X', 1.5, 0, -1, 0))


def test_alpha_stable_symbolic_parameters():
    """Test with symbolic parameters"""
    from sympy import symbols
    a, b, s, loc = symbols('alpha beta scale location', positive=True, real=True)

    # Should be able to create with symbolic parameters
    X = AlphaStable('X', a, b, s, loc)
    assert X.pspace.distribution.alpha == a
    assert X.pspace.distribution.beta == b
    assert X.pspace.distribution.scale == s
    assert X.pspace.distribution.location == loc
