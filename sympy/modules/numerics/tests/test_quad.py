import sys
sys.path.append(".")
import py
from sympy import *
from sympy.modules.numerics import *

def test_tanhsinh():
    t = TanhSinh()
    assert t(lambda x: x**3, -3, 2).ae(-16.25)
    assert t(lambda x: 2/(1+x**2), -1, 1).ae(pi_float())
    assert t(lambda x: 2/(1+x**2), 0, 'oo').ae(pi_float())
    assert t(lambda x: exp(-x), 0, 'oo').ae(1)
    assert t(lambda x: 2*exp(-x**2), 0, 'oo').ae(sqrt(pi_float()))

