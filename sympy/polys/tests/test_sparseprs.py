from __future__ import annotations

from sympy.polys.domains import ZZ
from sympy.polys.sparseprs import smp_prs_resultant, smp_subresultants


def test_smp_subresultants_zero_cases():
    f = {(2,): ZZ.one, (0,): ZZ(-1)}
    one = {(0,): ZZ.one}

    assert smp_subresultants({}, {}, 0, 1, ZZ) == [{}, {}]
    assert smp_subresultants(f, {}, 0, 1, ZZ) == [f, one]
    assert smp_subresultants({}, f, 0, 1, ZZ) == [f, one]


def test_smp_prs_resultant_zero_and_common_factor():
    f = {(2,): ZZ.one, (0,): ZZ(-1)}
    g = {(1,): ZZ.one, (0,): ZZ(-1)}

    assert smp_prs_resultant({}, g, 0, 1, ZZ) == ({}, [])
    assert smp_prs_resultant(f, {}, 0, 1, ZZ) == ({}, [])

    result, prs = smp_prs_resultant(f, g, 0, 1, ZZ)

    assert result == {}
    assert prs == [f, g]


def test_smp_prs_resultant_nonzero():
    f = {(2,): ZZ.one, (0,): ZZ(-1)}
    g = {(1,): ZZ.one, (0,): ZZ(-2)}

    result, prs = smp_prs_resultant(f, g, 0, 1, ZZ)

    assert result == {(0,): ZZ(3)}
    assert prs[0] == f
    assert prs[1] == g
