from sympy.combinatorics.named_groups import SymmetricGroup
from sympy.combinatorics.testutil import _verify_bsgs

def test_verify_bsgs():
    S = SymmetricGroup(5)
    S.schreier_sims()
    base = S.base
    strong_gens = S.strong_gens
    gens = S.generators
    assert _verify_bsgs(S, base, strong_gens) == True
    assert _verify_bsgs(S, base[:-1], strong_gens) == False
    assert _verify_bsgs(S, base, S.generators) == False

