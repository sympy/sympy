from __future__ import annotations

from sympy import Symbol, assuming, Q

def test_unified_assumptions_basic():
    x = Symbol('x')
    
    # 1. Verification of default status (None)
    assert x.is_positive is None

    # 2. Inside context manager, old property query should dynamically bridge to SAT solver
    with assuming(Q.positive(x)):
        assert x.is_positive is True

    # 3. Outside context manager, status should reset to None
    assert x.is_positive is None


def test_unified_assumptions_compound():
    x = Symbol('x')
    y = Symbol('y')
    expr = x - y

    assert expr.is_positive is None

    with assuming(Q.positive(expr)):
        assert expr.is_positive is True

    assert expr.is_positive is None


def test_unified_assumptions_cache_coherency():
    x = Symbol('x')
    
    # Pre-cache x.is_positive as None
    assert x.is_positive is None
    
    # Verify that the cache is successfully bypassed and returns True inside assuming context
    with assuming(Q.positive(x)):
        assert x.is_positive is True
        
    # Verify that exiting the context reverts back to None and does not leak the True value
    assert x.is_positive is None
