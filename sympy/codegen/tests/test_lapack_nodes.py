from __future__ import annotations

from sympy import MatrixSymbol, ccode, fcode
from sympy.codegen.lapack_nodes import dgesv_function

def test_dgesv_function__ccode():
    A = MatrixSymbol('A', 3, 3)
    b = MatrixSymbol('b', 3, 1)
    func = dgesv_function(A, b)
    code = ccode(func)
    assert 'LAPACKE_dgesv' in code
    assert 'int n = 3' in code
    assert 'int ipiv[3]' in code
    assert 'return info' in code

def test_dgesv_function__fcode():
    A = MatrixSymbol('A', 3, 3)
    b = MatrixSymbol('b', 3, 1)
    func = dgesv_function(A, b)
    code = fcode(func, standard=95, source_format='free')
    assert 'call dgesv' in code
    assert 'dimension(3)' in code
