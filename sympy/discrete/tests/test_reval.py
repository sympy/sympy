from __future__ import print_function, division

from sympy import sqrt
from sympy.core import S, Symbol, I
from sympy.core.compatibility import range
from sympy.utilities.pytest import raises
from sympy.discrete import reval_lhcc

def test_reval_lhcc():
    assert reval_lhcc([1, 1], [1, 1], n=20) == 10946
