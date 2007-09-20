import sys
sys.path.append(".")
import py
from sympy import *
from sympy.numerics import *
from sympy.numerics.functions import *
from sympy.numerics.functions2 import *
import math

def test_incomplete_gamma():
    Float.store()
    Float.setprec(53)
    assert upper_gamma(-2.5,-0.5).ae(-0.9453087204829418812-5.3164237738936178621j)
    assert erf(0) == 0
    assert erf(1).ae('0.84270079294971486934')
    assert erf(3+4j).ae(-120.186991395079444098 - 27.750337293623902498j)
    assert erf(-4-3j).ae(-0.99991066178539168236 + 0.00004972026054496604j)
    assert erf(pi_float()).ae(0.99999112385363235839)
    Float.revert()

def test_gamma():
    Float.store()
    Float.setprec(53)
    assert gamma(0.25).ae('3.6256099082219083119')
    assert gamma(0.0001).ae('9999.4228832316241908')
    assert gamma(300).ae('1.0201917073881354535e612')
    assert gamma(-0.5).ae('-3.5449077018110320546')
    assert gamma(-7.43).ae('0.00026524416464197007186')
    assert gamma(Rational(1,2)) == gamma(0.5)
    assert gamma(Rational(-7,3)).ae(gamma(Float(-7)/3))
    assert gamma(1+1j).ae(0.49801566811835604271 - 0.15494982830181068512j)
    assert gamma(-1+0.01j).ae(-0.422733904013474115 + 99.985883082635367436j)
    assert gamma(20+30j).ae(-1453876687.5534810 + 1163777777.8031573j)
    # Should always give exact factorials when they can
    # be represented as Floats under the current working precision
    fact = 1
    for i in range(1, 18):
        assert gamma(i) == fact
        fact *= i
    Float.setprec(1000)
    for i in range(18, 105):
        assert gamma(i) == fact
        fact *= i
    Float.setdps(100)
    assert gamma(0.5).ae(sqrt(pi_float()))
    Float.revert()

def test_zeta():
    Float.store()
    Float.setprec(53)
    assert zeta(2).ae(pi_float()**2 / 6)
    assert zeta(2.0).ae(pi_float()**2 / 6)
    assert zeta(ComplexFloat(2)).ae(pi_float()**2 / 6)
    assert zeta(100).ae(1)
    assert zeta(0).ae(-0.5)
    assert zeta(0.5).ae('-1.46035450880958681')
    assert zeta(-1).ae(-Float(1)/12)
    assert zeta(-2).ae(0)
    assert zeta(-3).ae(Float(1)/120)
    assert zeta(-4).ae(0)
    # Zeros in the critical strip
    assert zeta(ComplexFloat(0.5, 14.1347251417346937904)).ae(0)
    assert zeta(ComplexFloat(0.5, 21.0220396387715549926)).ae(0)
    assert zeta(ComplexFloat(0.5, 25.0108575801456887632)).ae(0)
    Float.setdps(50)
    im = '236.5242296658162058024755079556629786895294952121891237'
    assert zeta(ComplexFloat(0.5, im)).ae(0, 1e-46)
    Float.revert()
