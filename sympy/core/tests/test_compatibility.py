from __future__ import annotations

import pytest

from sympy.testing.pytest import warns_deprecated_sympy


@pytest.mark.thread_unsafe(
    reason="expects a warning side effect from a one-time module import"
)
def test_compatibility_submodule():
    # Test the sympy.core.compatibility deprecation warning
    with warns_deprecated_sympy():
        import sympy.core.compatibility # noqa:F401
