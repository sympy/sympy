from time import time
from sympy import Symbol, S
from sympy.tensor.tensor import (tensor_indices)
from sympy.tensor.group_factors import SuNGroupFactors, match_CijCij

"""
References

[1] P. Cvitanovic "Group Theory" version 9.0.1
"""

def test_SuN():
    N = Symbol('N')
    SN = SuNGroupFactors(N)
    C = SN.C
    g = SN.ASuN.metric
    i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14 = \
        tensor_indices('i0:15', SN.ASuN)
    theta = SN.theta

    t = C(i0,i1,i2)*C(-i0,-i2,i3)
    t = SN.rule_C2theta_all(t)
    assert t == -2*N*g(i1, i3)

    t = C(i0,i1,i2)*C(-i1,i3,i4)*C(i5,-i0,-i2)*C(-i5,-i4,-i3)
    t = SN.rule_C2theta_all(t)
    t = t.expand_coeff()
    assert t == 4*N**4 - 4*N**2

    # see eq.(1.1) and page 11 in Ref.[1] see also page 72
    # We apply rule_C_square to two parts of the tensor;
    # applying rule_C_square once to the whole tensor it is 2x slower.
    # Without applying rule_C_square it is roughly 50x slower.
    t1 = C(i0,i1,i2)*C(-i1,i5,i3)*C(-i2,i4,i7)*C(-i3,-i4,i6)
    t2 = C(-i6,i10,i11)*C(-i5,-i11,i8)*C(-i7,-i10,i9)*C(-i8,-i9,i12)
    t1 = SN.rule_C_square(t1)
    t2 = SN.rule_C_square(t2)
    t = t1*t2
    t = SN.rule_C2theta_all(t)
    t = t.expand_coeff()
    assert t == (2*N**4 + 24*N**2)*g(i0, i12)

    # without using rule_C_triangle it is 2.5x slower
    t = C(i0,i1,i2)*C(-i2,i3,i5)*C(-i1,-i3,i4)*C(-i4,-i5,i6)
    t = SN.rule_C_triangle(t)
    t = SN.rule_C2theta_all(t)
    t = t.expand_coeff()
    assert t == 2*N**2*g(i0, i6)

    # first graph in page 72 in Ref. [1]
    t = theta(i0,i1,i2,i3)*theta(-i3,-i2,-i1,-i0)
    t = SN.rule_C2theta_all(t)
    t = t.expand_coeff()
    assert (t._coeff - (N**4 - 3*N**2 + 3)*(N**2 - 1)/N**2).expand() == 0


    # Identically vanishing tensors
    # contracted according to the Kuratowski graph  eq.(6.59) in Ref. [1]
    t = C(i0,i1,i2)*C(-i1,i3,i4)*C(-i3,i7,i5)*C(-i2,-i5,i6)*C(-i4,-i6,i8)
    t1 = t.canon_bp()
    assert t1 == 0

    #  contracted according to the Peterson graph eq.(6.60) in Ref. [1]
    t = C(i0,i1,i2)*C(-i1,i3,i4)*C(-i2,i5,i6)*C(-i3,i7,i8)*C(-i6,-i7,i9)*\
        C(-i8,i10,i13)*C(-i5,-i10,i11)*C(-i4,-i11,i12)*C(-i9,-i12,i14)
    t1 = t.canon_bp()
    assert t1 == 0
