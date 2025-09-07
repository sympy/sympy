from sympy.core.mul import Mul
from sympy.core.singleton import S
from sympy.core.symbol import symbols
from sympy.core.add import Add
from sympy.core.power import Pow
from sympy.testing.pytest import raises

from sympy.physics.quantum.slidingtransform import SlidingTransform

a, b, c, d, e, f = symbols('a b c d e f', commutative=False)
x, y, z = symbols('x y z', commutative=True)


def doubler(x):
    return (x, x)

def drone(x):
    return (a,)

def mapping(lhs, rhs):
    if lhs == a and rhs == b:
        return (c,)
    elif lhs == b and rhs == c:
        return (d,)
    elif lhs == a and rhs == c:
        return None
    elif lhs == c and rhs == d:
        return (e, f)


def test_unary():
    st = SlidingTransform(unary=doubler)
    assert st(a*b*c) == Mul._from_args((a, a, b, b, c, c))

    st = SlidingTransform(unary=drone)
    assert st(a*b*c) == Mul._from_args((a, a, a))

    st = SlidingTransform(unary=lambda x: None)
    assert st(a*b*c) == Mul._from_args((a, b, c))


def test_binary():
    st = SlidingTransform(binary=mapping)

    assert st(a*b) == Mul._from_args((c,))
    assert st(b*c) == Mul._from_args((d,))
    assert st(a*c) == Mul._from_args((a,c))
    assert st(c*d) == Mul._from_args((e,f))

    assert st(a*b*c) == Mul._from_args((c, c))
    assert st(b*c*d) == Mul._from_args((d, d))
    assert st(c*d*e) == Mul._from_args((e, f, e))

    assert st(a*b*d) == Mul._from_args((e, f))


def test_reverse():
    st = SlidingTransform(binary=mapping, reverse=True)

    assert st(a*b) == Mul._from_args((c,))
    assert st(b*c) == Mul._from_args((d,))
    assert st(a*c) == Mul._from_args((a,c))
    assert st(c*d) == Mul._from_args((e,f))

    assert st(a*b*c) == Mul._from_args((a, d))
    assert st(d*b*c) == Mul._from_args((d, d))
    assert st(a*c*b*a*b) == Mul._from_args((a, e, f))


# Initialization and Property Tests

def test_init_basic():
    """Test basic initialization with default parameters."""
    st = SlidingTransform()
    assert st.unary is None
    assert st.binary is None
    assert st.reverse is False
    assert st.from_args is True


def test_init_all_parameters():
    """Test initialization with all parameters provided."""
    def dummy_unary(x): return None
    def dummy_binary(lhs, rhs): return None
    
    st = SlidingTransform(unary=dummy_unary, binary=dummy_binary, reverse=True, from_args=False)
    assert st.unary is dummy_unary
    assert st.binary is dummy_binary
    assert st.reverse is True
    assert st.from_args is False


def test_property_access():
    """Test that all property getters return the correct values after initialization."""
    def test_unary(x): return (x, x)
    def test_binary(lhs, rhs): return (lhs,)
    
    st = SlidingTransform(unary=test_unary, binary=test_binary, reverse=True, from_args=False)
    assert st.unary is test_unary
    assert st.binary is test_binary
    assert st.reverse is True
    assert st.from_args is False


# Input Validation Tests

def test_non_mul_raises_typeerror():
    """Test that calling the transform with non-Mul expressions raises TypeError."""
    st = SlidingTransform()
    
    with raises(TypeError):
        st(a)  # Symbol
    
    with raises(TypeError):
        st(a + b)  # Add
    
    with raises(TypeError):
        st(a**2)  # Pow
    
    with raises(TypeError):
        st(S(5))  # Integer


# Unary Transform Tests

def test_unary_identity_transform():
    """Test unary function that returns None for all inputs (no transformation)."""
    def identity_unary(x):
        return None
    
    st = SlidingTransform(unary=identity_unary)
    assert st(a*b*c) == Mul._from_args((a, b, c))
    assert st(a*x*b) == Mul._from_args((x, a, b))


def test_unary_expansion_transform():
    """Test unary function that expands single factors into multiple factors."""
    def tripler(x):
        return (x, x, x)
    
    st = SlidingTransform(unary=tripler)
    assert st(a*b) == Mul._from_args((a, a, a, b, b, b))
    
    # Test with expressions that remain as Mul objects
    # Create a Mul with two different non-commutative factors, then test just one
    result = st(a*c)  # This will triple both a and c
    expected = Mul._from_args((a, a, a, c, c, c))
    assert result == expected


def test_unary_replacement_transform():
    """Test unary function that replaces factors with different factors."""
    def replacer(x):
        if x == a:
            return (b, c)
        elif x == b:
            return (d,)
        return None
    
    st = SlidingTransform(unary=replacer)
    assert st(a*b*e) == Mul._from_args((b, c, d, e))
    assert st(a*e*b) == Mul._from_args((b, c, e, d))


def test_unary_zero_return():
    """Test unary function that returns (S.Zero,) for some input, causing entire expression to become zero."""
    def zero_maker(x):
        if x == a:
            return (S.Zero,)
        return None
    
    st = SlidingTransform(unary=zero_maker)
    assert st(a*b*c) == S.Zero
    assert st(b*c*d) == Mul._from_args((b, c, d))


def test_unary_with_commutative_factors():
    """Test unary transforms with mix of commutative and non-commutative factors."""
    def doubler(x):
        return (x, x)
    
    st = SlidingTransform(unary=doubler, from_args=True)
    assert st(a*x*b*y) == Mul._from_args((x, y, a, a, b, b), is_commutative=False)
    st = SlidingTransform(unary=doubler, from_args=False)
    assert st(a*x*b*y) == x*y*a**2*b**2


def test_unary_mixed_output():
    """Test unary function that returns mix of commutative and non-commutative factors."""
    def mixed_expander(x):
        if x == a:
            return (x, y, b)
        return None
    
    st = SlidingTransform(unary=mixed_expander)
    result = st(a*c)
    assert result == Mul._from_args((y, a, b, c), is_commutative=False)


# Binary Transform Tests

def test_binary_simple_pairs():
    """Test binary function on simple two-factor expressions."""
    def simple_binary(lhs, rhs):
        if lhs == a and rhs == b:
            return (c,)
        return None
    
    st = SlidingTransform(binary=simple_binary)
    assert st(a*b) == Mul._from_args((c,))
    assert st(b*a) == Mul._from_args((b, a))
    assert st(a*c) == Mul._from_args((a, c))


def test_binary_no_transformation():
    """Test binary function that returns None (no transformation) for all pairs."""
    def no_transform(lhs, rhs):
        return None
    
    st = SlidingTransform(binary=no_transform)
    assert st(a*b*c*d) == Mul._from_args((a, b, c, d))


def test_binary_cascading_effect():
    """Test that when binary transform returns multiple factors, the last factor feeds back."""
    def expanding_binary(lhs, rhs):
        if lhs == a and rhs == b:
            return (c, d, e)  # Multiple factors
        elif lhs == d and rhs == e:
            return (f,)
        return None
    
    st = SlidingTransform(binary=expanding_binary)
    # Based on the sliding transform logic:
    # a*b -> (c, d, e), output gets c, d, then e gets fed back
    # But there's no more input, so d and e remain as-is
    # The cascading only happens if there are more factors to process
    result = st(a*b)
    # The actual result should be c*d*e (no further cascading happens)
    assert result == c*d*e


def test_binary_reverse_direction():
    """Test binary processing in reverse direction (right-to-left)."""
    def directional_binary(lhs, rhs):
        if lhs == a and rhs == b:
            return (c,)
        elif lhs == b and rhs == c:
            return (d,)
        return None
    
    st_forward = SlidingTransform(binary=directional_binary, reverse=False)
    st_reverse = SlidingTransform(binary=directional_binary, reverse=True)
    
    # Forward: a*b -> c, then c*c (no match) -> (c, c)
    # Reverse: b*c -> d, then a*d (no match) -> (a, d)
    assert st_forward(a*b*c) == Mul._from_args((c, c))
    assert st_reverse(a*b*c) == Mul._from_args((a, d))


def test_binary_with_commutative_separation():
    """Test binary transform correctly skips commutative factors."""
    def nc_only_binary(lhs, rhs):
        if lhs == a and rhs == b:
            return (c,)
        return None
    
    st = SlidingTransform(binary=nc_only_binary)
    # Should skip commutative x, y and only process non-commutative pairs
    result = st(a*x*b*y*d)
    assert result == Mul._from_args((x, y, c, d), is_commutative=False)


def test_binary_zero_in_transform_output():
    """Test binary function that produces zero in commutative parts."""
    def zero_producing_binary(lhs, rhs):
        if lhs == a and rhs == b:
            return (S.Zero, c)  # Zero in commutative part
        return None
    
    st = SlidingTransform(binary=zero_producing_binary)
    # When S.Zero appears in commutative parts, it gets multiplied with other factors
    # The result will be S.Zero * c * d, which is still zero but not simplified
    assert st(a*b*d) == S.Zero


# Combined Unary-Binary Tests

def test_unary_then_binary():
    """Test full pipeline with both unary and binary transforms."""
    def unary_expander(x):
        if x == a:
            return (b, c)
        return None
    
    def binary_combiner(lhs, rhs):
        if lhs == b and rhs == c:
            return (d,)
        return None
    
    st = SlidingTransform(unary=unary_expander, binary=binary_combiner)
    # Test with a*e where e is untransformed by unary
    # a -> (b, c) via unary, then b*c -> d via binary, e remains
    assert st(a*e) == d*e


def test_unary_expands_for_binary():
    """Test scenario where unary transform creates more factors that binary processes."""
    def doubler_unary(x):
        return (x, x)
    
    def pair_reducer(lhs, rhs):
        if lhs == rhs:  # Same factor twice
            return (lhs,)  # Reduce to single
        return None
    
    st = SlidingTransform(unary=doubler_unary, binary=pair_reducer)
    # a*b -> (a, a, b, b) via unary, then a*a -> a, b*b -> b via binary
    assert st(a*b) == a*b


# Edge Cases and Complex Scenarios

def test_single_factor_expression():
    """Test expressions with single effective factors."""
    def unary_doubler(x):
        return (x, x)
    
    def binary_combiner(lhs, rhs):
        return (lhs,)
    
    st_unary = SlidingTransform(unary=unary_doubler)
    st_binary = SlidingTransform(binary=binary_combiner)
    st_both = SlidingTransform(unary=unary_doubler, binary=binary_combiner)
    
    # Test with expressions that create proper Mul objects
    assert st_unary(a*c) == Mul._from_args((a, a, c, c))  # Both factors doubled
    assert st_binary(a*c) == a  # Binary combines pairs to first element
    assert st_both(a*c) == a  # Unary doubles, binary reduces


def test_all_commutative_factors():
    """Test expressions containing only commutative factors."""
    def unary_transform(x):
        return (x, x)
    
    def binary_transform(lhs, rhs):
        return (lhs,)
    
    st = SlidingTransform(unary=unary_transform, binary=binary_transform)
    # Commutitive factors should not be transformed at all
    assert st(x*y*z) == x*y*z


def test_alternating_commutative_noncommutative():
    """Test expressions with alternating commutative and non-commutative factors."""
    def simple_binary(lhs, rhs):
        if lhs == a and rhs == b:
            return (c,)
        return None
    
    st = SlidingTransform(binary=simple_binary)
    # Should skip commutative factors and only process a*b pair
    result = st(x*a*y*b*z)
    assert result == Mul._from_args((x, y, z, c), is_commutative=False)


def test_long_chain_transformations():
    """Test expressions with many factors that undergo multiple rounds of binary transformations."""
    def chain_binary(lhs, rhs):
        # Create a chain: a*b -> c, c*d -> e, e*f -> a (cycle)
        if lhs == a and rhs == b:
            return (c,)
        elif lhs == c and rhs == d:
            return (e,)
        elif lhs == e and rhs == f:
            return (a,)
        return None
    
    st = SlidingTransform(binary=chain_binary)
    # a*b*d*f -> c*d*f -> e*f -> a
    assert st(a*b*d*f) == Mul._from_args((a,))


def test_complex_cascading():
    """Test scenario where binary transform produces many factors, creating complex cascading."""
    def complex_binary(lhs, rhs):
        if lhs == a and rhs == b:
            return (c, d, e, f)  # Produces many factors
        elif lhs == e and rhs == f:
            return (a,)
        return None
    
    st = SlidingTransform(binary=complex_binary)
    # a*b -> (c, d, e, f) where f is the last item and goes back to input
    # But with just a*b, there's no more input to process f with
    # So the result is just c*d*e*f without further cascading
    result = st(a*b)
    assert result == c*d*e*f


# from_args Parameter Tests

def test_from_args():
    """Test with from_args=True to verify Mul._from_args is used."""
    st = SlidingTransform(unary=doubler, from_args=True)
    # Should use _from_args, avoiding post-processors
    assert st(x*a*b) == Mul._from_args((x, a, a, b, b), is_commutative=False)
    st = SlidingTransform(unary=doubler, from_args=False)
    # Should use _from_args, avoiding post-processors
    assert st(x*a*b) == x*a**2*b**2


# Commutative/Non-commutative Handling Tests

def test_commutative_collection():
    """Test that commutative factors from both unary and binary transforms are properly collected."""
    def unary_mixed(x):
        if x == a:
            return (x, y)  # y is commutative
        return None
    
    def binary_mixed(lhs, rhs):
        if lhs == b and rhs == c:
            return (z, d)  # z is commutative
        return None
    
    st = SlidingTransform(unary=unary_mixed, binary=binary_mixed)
    result = st(a*b*c)
    # Should collect y and z as commutative factors
    assert result == Mul._from_args((y, z, a, d), is_commutative=False)


def test_split_cnc_functionality():
    """Test the internal _split_cnc function behavior through transform operations."""
    def mixed_output_unary(x):
        return (x, y, z, b)  # Mix of commutative (y, z) and non-commutative (x, b)
    
    st = SlidingTransform(unary=mixed_output_unary)
    result = st(a*c)  # Use proper Mul with two non-commutative factors
    # a -> (a, y, z, b) and c -> (c, y, z, b)
    # Commutative parts: y, z, y, z -> combined
    # Non-commutative parts: a, b, c, b
    # Final result: y*y*z*z * a*b*c*b
    expected = y**2*z**2*a*b*c*b
    assert result == expected


def test_mixed_result_reconstruction():
    """Test reconstruction of final Mul with both commutative and non-commutative parts."""
    def add_commutative_unary(x):
        return (x, x, y)  # Adds commutative factor
    
    st = SlidingTransform(unary=add_commutative_unary)
    result = st(a*b)
    # Should have both commutative and non-commutative parts properly combined
    expected_commutative = y * y  # y from both a and b
    assert result == Mul._from_args((expected_commutative, a, a, b, b), is_commutative=False)


# Error Handling and Robustness Tests

def test_transform_function_exceptions():
    """Test behavior when unary or binary functions raise exceptions."""
    def failing_unary(x):
        if x == a:
            raise ValueError("Test exception")
        return None
    
    def failing_binary(lhs, rhs):
        if lhs == a and rhs == b:
            raise RuntimeError("Test exception")
        return None
    
    st_unary = SlidingTransform(unary=failing_unary)
    st_binary = SlidingTransform(binary=failing_binary)
    
    with raises(ValueError):
        st_unary(a*b)
    
    with raises(RuntimeError):
        st_binary(a*b)


def test_malformed_transform_output():
    """Test behavior when transform functions return invalid output."""
    def bad_unary_string(x):
        return "not a tuple"  # Invalid return type
    
    def bad_binary_int(lhs, rhs):
        return 42  # Invalid return type
    
    st_unary = SlidingTransform(unary=bad_unary_string)
    st_binary = SlidingTransform(binary=bad_binary_int)
    
    # These should raise appropriate errors during processing
    with raises((TypeError, AttributeError)):
        st_unary(a*b)
    
    with raises((TypeError, AttributeError)):
        st_binary(a*b)
