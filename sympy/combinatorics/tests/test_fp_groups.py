# -*- coding: utf-8 -*-
from sympy.combinatorics.fp_groups import FpGroup, CosetTable, coset_enumeration_r
from sympy.combinatorics.free_group import free_group

"""
References
==========

[1] Holt, D., Eick, B., O'Brien, E.
"Handbook of Computational Group Theory"

[2] John J. Cannon; Lucien A. Dimino; George Havas; Jane M. Watson
Mathematics of Computation, Vol. 27, No. 123. (Jul., 1973), pp. 463-490.
"Implementation and Analysis of the Todd-Coxeter Algorithm"

"""


def test_scan_2():
    # Example 5.1 [1]
    F, x, y = free_group("x, y")
    f = FpGroup(F, [x**3, y**3, x**-1*y**-1*x*y])
    c = CosetTable(f, [x])

    c.scan_and_fill(0, x)
    assert c.table == [[0, 0, None, None]]
    assert c.p == [0]
    assert c.n == 1
    assert c.omega == [0]

    c.scan_and_fill(0, x**3)
    assert c.table == [[0, 0, None, None]]
    assert c.p == [0]
    assert c.n == 1
    assert c.omega == [0]

    c.scan_and_fill(0, y**3)
    assert c.table == [[0, 0, 1, 2], [None, None, 2, 0], [None, None, 0, 1]]
    assert c.p == [0, 1, 2]
    assert c.n == 3
    assert c.omega == [0, 1, 2]

    c.scan_and_fill(0, x**-1*y**-1*x*y)
    assert c.table == [[0, 0, 1, 2], [None, None, 2, 0], [2, 2, 0, 1]]
    assert c.p == [0, 1, 2]
    assert c.n == 3
    assert c.omega == [0, 1, 2]

    c.scan_and_fill(1, x**3)
    assert c.table == [[0, 0, 1, 2], [3, 4, 2, 0], [2, 2, 0, 1], \
            [4, 1, None, None], [1, 3, None, None]]
    assert c.p == [0, 1, 2, 3, 4]
    assert c.n == 5
    assert c.omega == [0, 1, 2, 3, 4]

    c.scan_and_fill(1, y**3)
    assert c.table == [[0, 0, 1, 2], [3, 4, 2, 0], [2, 2, 0, 1], \
            [4, 1, None, None], [1, 3, None, None]]
    assert c.p == [0, 1, 2, 3, 4]
    assert c.n == 5
    assert c.omega == [0, 1, 2, 3, 4]

    c.scan_and_fill(1, x**-1*y**-1*x*y)
    assert c.table == [[0, 0, 1, 2], [1, 1, 2, 0], [2, 2, 0, 1], \
            [None, 1, None, None], [1, 3, None, None]]
    assert c.p == [0, 1, 2, 1, 1]
    assert c.n == 3
    assert c.omega == [0, 1, 2, 1, 1]

    # Example 5.2 [1]
    f = FpGroup(F, [x**2, y**3, (x*y)**3])
    c = CosetTable(f, [x*y])

    c.scan_and_fill(0, x*y)
    assert c.table == [[1, None, None, 1], [None, 0, 0, None]]
    assert c.p == [0, 1]
    assert c.n == 2
    assert c.omega == [0, 1]

    c.scan_and_fill(0, x**2)
    assert c.table == [[1, 1, None, 1], [0, 0, 0, None]]
    assert c.p == [0, 1]
    assert c.n == 2
    assert c.omega == [0, 1]

    c.scan_and_fill(0, y**3)
    assert c.table == [[1, 1, 2, 1], [0, 0, 0, 2], [None, None, 1, 0]]
    assert c.p == [0, 1, 2]
    assert c.n == 3
    assert c.omega == [0, 1, 2]

    c.scan_and_fill(0, (x*y)**3)
    assert c.table == [[1, 1, 2, 1], [0, 0, 0, 2], [None, None, 1, 0]]
    assert c.p == [0, 1, 2]
    assert c.n == 3
    assert c.omega == [0, 1, 2]

    c.scan_and_fill(1, x**2)
    assert c.table == [[1, 1, 2, 1], [0, 0, 0, 2], [None, None, 1, 0]]
    assert c.p == [0, 1, 2]
    assert c.n == 3
    assert c.omega == [0, 1, 2]

    c.scan_and_fill(1, y**3)
    assert c.table == [[1, 1, 2, 1], [0, 0, 0, 2], [None, None, 1, 0]]
    assert c.p == [0, 1, 2]
    assert c.n == 3
    assert c.omega == [0, 1, 2]

    c.scan_and_fill(1, (x*y)**3)
    assert c.table == [[1, 1, 2, 1], [0, 0, 0, 2], [3, 4, 1, 0], [None, 2, 4, None], [2, None, None, 3]]
    assert c.p == [0, 1, 2, 3, 4]
    assert c.n == 5
    assert c.omega == [0, 1, 2, 3, 4]

    c.scan_and_fill(2, x**2)
    assert c.table == [[1, 1, 2, 1], [0, 0, 0, 2], [3, 3, 1, 0], [2, 2, 3, 3], [2, None, None, 3]]
    assert c.p == [0, 1, 2, 3, 3]
    assert c.n == 4
    assert c.omega == [0, 1, 2, 3, 3]


def test_HLT_method():
    # Example 5.1 [1]
    F, x, y = free_group("x, y")
    f = FpGroup(F, [x**3, y**3, x**-1*y**-1*x*y])
    C = coset_enumeration_r(f, [x])
    c_1 = [[0, 0, 1, 2], [1, 1, 2, 0], [2, 2, 0, 1]]
    for i in range(len(C.p)):
        if C.p[i] == i:
            assert C.table[i] == c_1[i]

    # Group denoted by E₁
    F, r, s, t = free_group("r, s, t")
    E1 = FpGroup(F, [t**-1*r*t*r**-2, r**-1*s*r*s**-2, s**-1*t*s*t**-2])
    C = coset_enumeration_r(E1, [r])
    c_E1 = [0, 0, 0, 0, 0, 0]
    for i in range(len(C.p)):
        if C.p[i] == i:
            C.table[i] == c_E1[i]

    # Group denoted by Cox
    F, a, b = free_group("a, b")
    Cox = FpGroup(F, [a**6, b**6, (a*b)**2, (a**2*b**2)**2, (a**3*b**3)**5])
    C = coset_enumeration_r(Cox, [a])
    index_cox = 0
    for i in range(len(C.p)):
        if C.p[i] == i:
            index_cox += 1
    assert index_cox == 500

    # Group denoted by B₁,₂
    F, a, b = free_group("a, b")
    B_2_4 = FpGroup(F, [a**4, b**4, (a*b)**4, (a**-1*b)**4, (a**2*b)**4, (a*b**2)**4, (a**2*b**2)**4, (a**-1*b*a*b)**4, (a*b**-1*a*b)**4])
    C = coset_enumeration_r(B_2_4, [a])
    index_b_2_4 = 0
    for i in range(len(C.p)):
        if C.p[i] == i:
            index_b_2_4 += 1
    assert index_b_2_4 == 1024
