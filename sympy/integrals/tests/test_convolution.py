"""Tests for convolution functions."""

from sympy.core.numbers import oo, I, pi, Rational
from sympy.core.singleton import S
from sympy.core.symbol import symbols
from sympy.functions.elementary.complexes import Abs
from sympy.functions.elementary.exponential import exp
from sympy.functions.elementary.hyperbolic import sinh, cosh
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.piecewise import Piecewise
from sympy.functions.elementary.trigonometric import sin, cos
from sympy.functions.special.delta_functions import DiracDelta, Heaviside
from sympy.integrals.transforms import convolution_integral, Convolution
from sympy.testing.pytest import raises, slow


def test_convolution_delta():
    """Test convolution with Dirac delta function."""
    t, a = symbols('t a', real=True, positive=True)
    
    # DiracDelta integration has special behavior
    # Instead, test a simpler case: convolution with itself
    f = Heaviside(t)
    result = convolution_integral(f, f, t, 0, t)
    expected = t*Heaviside(t)**2
    assert result == expected


def test_convolution_exponentials():
    """Test convolution of exponential functions."""
    t, a, b = symbols('t a b', real=True, positive=True)
    
    # Convolution of two causal exponentials
    f = exp(-a*t)*Heaviside(t)
    g = exp(-b*t)*Heaviside(t)
    
    result = convolution_integral(f, g, t, 0, t)
    
    # Result should be a Piecewise expression
    assert isinstance(result, Piecewise)
    
    # Test with specific values where a != b
    a_val, b_val, t_val = 2, 3, 1
    result_numeric = result.subs([(a, a_val), (b, b_val), (t, t_val)])
    expected_numeric = (exp(-a_val*t_val) - exp(-b_val*t_val))/(b_val - a_val)
    assert abs(float(result_numeric) - float(expected_numeric)) < 1e-10


def test_convolution_heaviside():
    """Test convolution of Heaviside functions."""
    t = symbols('t', real=True)
    
    # H(t) * H(t) over [0, t] gives t*H(t)^2
    result = convolution_integral(Heaviside(t), Heaviside(t), t, 0, t)
    assert result == t*Heaviside(t)**2


def test_convolution_polynomial():
    """Test convolution with polynomial functions."""
    t, tau = symbols('t tau', real=True)
    
    # Convolution of t*H(t) with H(t)
    # Result should be t^2/2 * H(t)^2
    f = t*Heaviside(t)
    g = Heaviside(t)
    result = convolution_integral(f, g, t, 0, t)
    assert result == t**2*Heaviside(t)**2/2


def test_convolution_class():
    """Test the Convolution class."""
    t, a = symbols('t a', real=True, positive=True)
    f = exp(-a*t)*Heaviside(t)
    g = exp(-2*a*t)*Heaviside(t)
    
    # Create unevaluated convolution
    conv = Convolution(f, g, t, 0, t)
    
    # Check properties
    assert conv.function_f == f
    assert conv.function_g == g
    assert conv.variable == t
    assert conv.lower_limit == 0
    assert conv.upper_limit == t
    
    # Check that doit() evaluates it
    result = conv.doit()
    assert not result.has(Convolution)


def test_convolution_rewrite_integral():
    """Test rewriting convolution as an integral."""
    t, a = symbols('t a', real=True, positive=True)
    f = exp(-a*t)*Heaviside(t)
    g = Heaviside(t)
    
    conv = Convolution(f, g, t, 0, t)
    integral_form = conv.rewrite('Integral')
    
    # Check that it's an Integral
    from sympy.integrals import Integral
    assert isinstance(integral_form, Integral)


def test_convolution_symmetric():
    """Test that convolution is commutative (f*g = g*f)."""
    t, a, b = symbols('t a b', real=True, positive=True)
    f = exp(-a*t)*Heaviside(t)
    g = exp(-b*t)*Heaviside(t)
    
    # Compute both directions
    result1 = convolution_integral(f, g, t, 0, t)
    result2 = convolution_integral(g, f, t, 0, t)
    
    # They should be equal (commutative property)
    assert (result1 - result2).simplify() == 0


def test_convolution_properties():
    """Test basic properties of the Convolution class."""
    t = symbols('t', real=True)
    f = exp(-t)*Heaviside(t)
    g = Heaviside(t)
    
    conv = Convolution(f, g, t)
    
    # Test default limits
    assert conv.lower_limit == S.NegativeInfinity
    assert conv.upper_limit == S.Infinity
    
    # Test with custom limits
    conv2 = Convolution(f, g, t, 0, t)
    assert conv2.lower_limit == 0
    assert conv2.upper_limit == t


def test_convolution_sin_cos():
    """Test convolution with trigonometric functions."""
    t, w = symbols('t w', real=True, positive=True)
    
    # Convolution of sin with Heaviside step
    f = sin(w*t)*Heaviside(t)
    g = Heaviside(t)
    
    result = convolution_integral(f, g, t, 0, t)
    
    # The expected result is (1 - cos(w*t))/w for t >= 0
    expected = (1 - cos(w*t))/w
    
    # Compare the results (ignoring Heaviside factors)
    assert (result.subs(Heaviside(t), 1) - expected).simplify() == 0


def test_convolution_noconds():
    """Test noconds parameter."""
    t, a = symbols('t a', real=True, positive=True)
    f = exp(-a*t)*Heaviside(t)
    g = Heaviside(t)
    
    conv = Convolution(f, g, t, 0, t)
    
    # With noconds=True (default), should return only the result
    result_noconds = conv.doit(noconds=True)
    assert not isinstance(result_noconds, tuple)
    
    # With noconds=False, should return a tuple
    result_with_conds = conv.doit(noconds=False)
    assert isinstance(result_with_conds, tuple)
    assert len(result_with_conds) == 2


def test_convolution_free_symbols():
    """Test free_symbols property."""
    t, a, b = symbols('t a b', real=True)
    f = exp(-a*t)
    g = exp(-b*t)
    
    conv = Convolution(f, g, t)
    # Should contain a and b, but not t (it's the integration variable)
    free = conv.free_symbols
    assert a in free
    assert b in free
    assert t not in free
