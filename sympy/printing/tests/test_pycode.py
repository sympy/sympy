# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)
from sympy.core import Mod, symbols
from sympy.core.numbers import pi
from sympy.logic import And, Or
from sympy.functions import Piecewise, acos
from sympy.matrices import SparseMatrix
from sympy.printing.pycode import PythonCodePrinter, SciPyPrinter


x, y, z = symbols('x y z')


def test_PythonCodePrinter():
    prntr = PythonCodePrinter()
    assert not prntr.modules
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
    assert not prntr.modules
    assert prntr.doprint(pi) == 'math.pi'
    assert prntr.modules == {'math'}
    assert prntr.doprint(acos(x)) == 'math.acos(x)'


def test_SciPyPrinter():
    p = SciPyPrinter()
    expr = acos(x)
    assert 'numpy' not in p.modules
    assert p.doprint(expr) == 'numpy.arccos(x)'
    assert 'numpy' in p.modules
    assert not any(m.startswith('scipy') for m in p.modules)
    smat = SparseMatrix(2, 5, {(0, 1): 3})
    assert p.doprint(smat) == 'scipy.sparse.coo_matrix([3], ([0], [1]), shape=(2, 5))'
    assert 'scipy.sparse' in p.modules
