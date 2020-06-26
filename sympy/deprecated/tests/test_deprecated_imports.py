import sympy
from sympy.testing.pytest import warns_deprecated_sympy

def test_deprecated_imports():
    # https://github.com/sympy/sympy/pull/18245
    # Before 1.6 these names were importable with e.g.
    # from sympy import *
    # from sympy import add
    # Now sympy/__init__.py uses __all__ so these names are no longer
    # accidentally imported.  However many of the names now give a warning and
    # this test checks that they are importable but a warning is given
    from sympy import add

    with warns_deprecated_sympy():
        add.Add

    modnames = type(add)._DEPRECATED_IMPORTS

    assert len(modnames) == 80

    for modname in modnames:
        name = modname.split('.')[-1]
        mod = getattr(sympy, name)
        attr = dir(mod.mod)[0]
        with warns_deprecated_sympy():
            getattr(mod, attr)
