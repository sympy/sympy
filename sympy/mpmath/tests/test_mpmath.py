from sympy.mpmath.libmp import *
from sympy.mpmath import *

def test_newstyle_classes():
    for cls in [mp, fp, iv, mpf, mpc]:
        for s in cls.__class__.__mro__:
            assert isinstance(s, type)
