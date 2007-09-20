import sys
sys.path.append(".")
import py
from sympy import *
from sympy.numerics import *
from sympy.numerics.functions import *
from sympy.numerics.quad import *

def test_nintegrate():
    Float.store()
    Float.setdps(20)
    pi_ = pi_float()
    assert nintegrate(lambda x: sin(x), 0, pi_).ae(2)
    assert nintegrate(lambda x: abs(sin(x)), 0, 10*pi_).ae(20)
    assert nintegrate(lambda x: sin(x), 0, 10*pi_).ae(0)
    assert nintegrate(lambda x: 4/(1+x**2), 0, 1).ae(pi_)
    assert nintegrate(lambda x: 4*sqrt(1-x**2), 0, 1).ae(pi_)
    Float.revert()

def test_nintegrate_infinite():
    Float.store()
    Float.setdps(15)
    pi_ = pi_float()
    assert nintegrate(lambda x: 4/(1+x**2), 1, oo).ae(pi_)
    A = nintegrate(lambda x: 2 * exp(-x**2), 0, oo)
    B = nintegrate(lambda x: 2 * exp(-x**2), -oo, 0)
    C = nintegrate(lambda x: 2 * exp(-x**2), -oo, oo)
    D = nintegrate(lambda x: 2 * exp(-x**2), 1, oo)
    E = nintegrate(lambda x: 2 * exp(-x**2), -1, oo)
    F = nintegrate(lambda x: 2 * exp(-x**2), -oo, -1)
    G = nintegrate(lambda x: 2 * exp(-x**2), -oo, 1)
    assert A.ae(pi_ ** 0.5)
    assert A.ae(B)
    assert C.ae(2*B)
    assert D.ae(0.27880558528066197650)
    assert E.ae(3.2661021165303700781)
    assert F.ae(D)
    assert G.ae(E)
    Float.revert()

def test_tanhsinh():
    Float.store()
    Float.setdps(15)
    assert nintegrate(lambda x: x**3, -3, 2, method=1).ae(-16.25)
    assert nintegrate(lambda x: 2/(1+x**2), -1, 1, method=1).ae(pi_float())
    assert nintegrate(lambda x: 2/(1+x**2), 0, oo, method=1).ae(pi_float())
    assert nintegrate(lambda x: exp(-x), 0, oo, method=1).ae(1)
    assert nintegrate(lambda x: 2*exp(-x**2), 0, oo, method=1).ae(sqrt(pi_float()))
    Float.revert()
