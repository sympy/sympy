from sympy.utilities.source import get_mod_func, get_class, source
from sympy.utilities.pytest import raises
from sympy import point
from sympy.utilities.exceptions import SymPyDeprecationWarning

def test_source():
    with raises(SymPyDeprecationWarning):
        source(point)

def test_get_mod_func():
    assert get_mod_func(
        'sympy.core.basic.Basic') == ('sympy.core.basic', 'Basic')


def test_get_class():
    _basic = get_class('sympy.core.basic.Basic')
    assert _basic.__name__ == 'Basic'
