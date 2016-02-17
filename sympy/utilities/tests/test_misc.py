from sympy.core.compatibility import unichr
from sympy.utilities.misc import translate, replace

def test_translate():
    abc = 'abc'
    translate(abc, None, 'a') == 'bc'
    translate(abc, None, '') == 'abc'
    translate(abc, {'a': 'x'}, 'c') == 'xb'
    assert translate(abc, {'a': 'bc'}, 'c') == 'bcb'
    assert translate(abc, {'ab': 'x'}, 'c') == 'x'
    assert translate(abc, {'ab': ''}, 'c') == ''
    assert translate(abc, {'bc': 'x'}, 'c') == 'ab'
    assert translate(abc, {'abc': 'x', 'a': 'y'}) == 'x'
    u = unichr(4096)
    assert translate(abc, 'a', 'x', u) == 'xbc'
    assert (u in translate(abc, 'a', u, u)) is True


def test_replace():
    assert replace('abc', ('a', 'b')) == 'bbc'
    assert replace('abc', {'a': 'Aa'}) == 'Aabc'
    assert replace('abc', ('a', 'b'), ('c', 'C')) == 'bbC'
