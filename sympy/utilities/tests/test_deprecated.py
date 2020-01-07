from sympy.testing.pytest import warns_deprecated_sympy, XFAIL

# See https://github.com/sympy/sympy/pull/18095

def test_deprecated_utilities():
    with warns_deprecated_sympy():
        import sympy.utilities.pytest
    with warns_deprecated_sympy():
        import sympy.utilities.runtests
    with warns_deprecated_sympy():
        import sympy.utilities.randtest
    with warns_deprecated_sympy():
        import sympy.utilities.tmpfiles
    with warns_deprecated_sympy():
        import sympy.utilities.quality_unicode

# This fails because benchmarking isn't importable...
@XFAIL
def test_deprecated_benchmarking():
    with warns_deprecated_sympy():
        import sympy.utilities.benchmarking
