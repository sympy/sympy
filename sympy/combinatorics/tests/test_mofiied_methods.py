from sympy.combinatorics.fp_groups import FpGroup
from sympy.combinatorics.coset_table import (CosetTable, coset_enumeration_r)
from sympy.combinatorics.free_groups import free_group

def test_modified_methods():
    '''
    Tests for modified coset table methods.
    Example 5.7 from [1] Holt, D., Eick, B., O'Brien
    "Handbook of Computational Group Theory".

    '''
    F, x, y = free_group("x, y")
    f = FpGroup(F, [x**3, y**5, (x*y)**2])
    H = [x*y, x**-1*y**-1*x*y*x]
    C = CosetTable(f, H)
    C.modified_define(0, x)
    identity = C._grp.identity
    a_0 = C._grp.generators[0]
    a_1 = C._grp.generators[1]

    assert C.P == [[identity, None, None, None],
                    [None, identity, None, None]]
    assert C.table == [[1, None, None, None],
                        [None, 0, None, None]]

    C.modified_define(1, x)
    assert C.table == [[1, None, None, None],
                        [2, 0, None, None],
                        [None, 1, None, None]]
    assert C.P == [[identity, None, None, None],
                    [identity, identity, None, None],
                    [None, identity, None, None]]

    C.modified_scan(0, x**3, C._grp.identity, fill=False)
    assert C.P == [[identity, identity, None, None],
                     [identity, identity, None, None],
                     [identity, identity, None, None]]
    assert C.table == [[1, 2, None, None],
                        [2, 0, None, None],
                        [0, 1, None, None]]

    C.modified_scan(0, x*y, C._grp.generators[0], fill=False)
    assert C.P == [[identity, identity, None, a_0**-1],
                    [identity, identity, a_0, None],
                    [identity, identity, None, None]]
    assert C.table == [[1, 2, None, 1],
                        [2, 0, 0, None],
                        [0, 1, None, None]]

    C.modified_define(2, y**-1)
    assert C.table == [[1, 2, None, 1],
                        [2, 0, 0, None],
                        [0, 1, None, 3],
                        [None, None, 2, None]]
    assert C.P == [[identity, identity, None, a_0**-1],
                    [identity, identity, a_0, None],
                    [identity, identity, None, identity],
                    [None, None, identity, None]]

    C.modified_scan(0, x**-1*y**-1*x*y*x, C._grp.generators[1])
    assert C.table == [[1, 2, None, 1],
                        [2, 0, 0, None],
                        [0, 1, None, 3],
                        [3, 3, 2, None]]
    assert C.P == [[identity, identity, None, a_0**-1],
                    [identity, identity, a_0, None],
                    [identity, identity, None, identity],
                    [a_1, a_1**-1, identity, None]]

    C.modified_scan(2, (x*y)**2, C._grp.identity)
    assert C.table == [[1, 2, 3, 1],
                        [2, 0, 0, None],
                        [0, 1, None, 3],
                        [3, 3, 2, 0]]
    assert C.P == [[identity, identity, a_1**-1, a_0**-1],
                    [identity, identity, a_0, None],
                    [identity, identity, None, identity],
                    [a_1, a_1**-1, identity, a_1]]

    C.modified_define(2, y)
    assert C.table == [[1, 2, 3, 1],
                        [2, 0, 0, None],
                        [0, 1, 4, 3],
                        [3, 3, 2, 0],
                        [None, None, None, 2]]
    assert C.P == [[identity, identity, a_1**-1, a_0**-1],
                    [identity, identity, a_0, None],
                    [identity, identity, identity, identity],
                    [a_1, a_1**-1, identity, a_1],
                    [None, None, None, identity]]

    C.modified_scan(0, y**5, C._grp.identity)
    assert C.table == [[1, 2, 3, 1], [2, 0, 0, 4], [0, 1, 4, 3], [3, 3, 2, 0], [None, None, 1, 2]]
    assert C.P == [[identity, identity, a_1**-1, a_0**-1],
                    [identity, identity, a_0, a_0*a_1**-1],
                    [identity, identity, identity, identity],
                    [a_1, a_1**-1, identity, a_1],
                    [None, None, a_1*a_0**-1, identity]]

    C.modified_scan(1, (x*y)**2, C._grp.identity)
    assert C.table == [[1, 2, 3, 1],
                        [2, 0, 0, 4],
                        [0, 1, 4, 3],
                        [3, 3, 2, 0],
                        [4, 4, 1, 2]]
    assert C.P == [[identity, identity, a_1**-1, a_0**-1],
                    [identity, identity, a_0, a_0*a_1**-1],
                    [identity, identity, identity, identity],
                    [a_1, a_1**-1, identity, a_1],
                    [a_0*a_1**-1, a_1*a_0**-1, a_1*a_0**-1, identity]]
