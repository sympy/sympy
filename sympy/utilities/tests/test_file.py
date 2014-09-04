from sympy.utilities.file import SympyDict
from sympy import symbols, srepr


def test_sympy_dict():
    a, b = symbols('a, b')
    d = SympyDict({'a': a, 'b': b})
    d.save('.test_file_sympy.swp')
    del d
    with open('.test_file_sympy.swp', 'r') as f:
        s = f.read()
    assert(s =="""{'a': "Symbol('a')",
 'b': "Symbol('b')"}""")
    d2 = SympyDict.load('test_file_sympy')
    assert srepr(d2['a']) == "Symbol('a')"
    assert srepr(d2['b']) == "Symbol('b')"
