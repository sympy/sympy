# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)
from sympy.core import Mod, symbols
from sympy.logic import And, Or
from sympy.functions import Piecewise
from sympy.printing.pycode import PythonCodePrinter

x, y, z = symbols('x y z')


def test_PythonCodePrinter():
    prntr = PythonCodePrinter()
    assert prntr.doprint(x**y) == 'x**y'
    assert prntr.doprint(Mod(x, 2)) == 'x % 2'
    assert prntr.doprint(And(x, y)) == 'x and y'
    assert prntr.doprint(Or(x, y)) == 'x or y'
    assert prntr.doprint(Piecewise((x, x > 1), (y, True))) == (
        'if x > 1:\n'
        '    return x\n'
        'else:\n'
        '    return y'
    )
    pw = Piecewise((x, x > 1), (y, x > 0))
    assert prntr.doprint(pw) == (
        'if x > 1:\n'
        '    return x\n'
        'elif x > 0:\n'
        '    return y\n'
        'else:\n'
        '    raise NotImplementedError("Unhandled condition in: %s")' % pw
    )
