from sympy.core.symbol import symbols
from sympy.core.function import Function
from sympy.assumptions.ask import Q
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.logic.inference import satisfiable
from sympy.logic.algorithms.dpll2 import dpll_satisfiable
from sympy.logic.algorithms.euf_theory_solver import EUFTheorySolver, TRUE
from sympy.logic.boolalg import And, Or
from sympy.logic.algorithms.euf_theory import EUFCongruenceClosure
from sympy.core.random import choice, randint, sample
from sympy.testing.pytest import raises
from itertools import product

f, g = symbols('f g', cls=Function)
a, b, c, d, x, y = symbols('a b c d x y')


def _encode(prop):
    enc = EncodedCNF()
    enc.from_cnf(CNF.from_prop(prop))
    return enc


def _solver(prop):
    return EUFTheorySolver.from_encoded_cnf(_encode(prop))


def test_from_encoded_cnf_collects_equality_atoms():
    enc = _encode(Q.eq(a, b) & Q.ne(c, d))
    euf, conflicts = EUFTheorySolver.from_encoded_cnf(enc)
    assert conflicts == []
    assert euf.atom_id_to_equality == {
        enc.encoding[Q.eq(a, b)]: (a, b, True),
        enc.encoding[Q.ne(c, d)]: (c, d, False)}


def test_from_encoded_cnf_decides_reflexive_atoms():
    enc = _encode(Q.eq(a, a) & Q.ne(b, b))
    euf, conflicts = EUFTheorySolver.from_encoded_cnf(enc)
    # x = x holds and x != x does not, so neither reaches the closure
    assert euf.atom_id_to_equality == {}
    assert sorted(conflicts) == sorted([[enc.encoding[Q.eq(a, a)]],
                                        [-enc.encoding[Q.ne(b, b)]]])


def test_predicate_atoms_are_reified():
    # Q.positive(x) becomes the equation Q.positive(x) = TRUE.
    enc = _encode(Q.positive(x) & Q.eq(a, b))
    euf, conflicts = EUFTheorySolver.from_encoded_cnf(enc)
    assert conflicts == []
    assert euf.atom_id_to_equality[enc.encoding[Q.positive(x)]] == (
        Q.positive(x), TRUE, True)
    assert euf.assert_lit(enc.encoding[Q.positive(x)]) is None
    assert euf.check()[0] is True


def test_predicate_congruence():
    # a = b forces Q.positive(a) and Q.positive(b) to agree.
    enc = _encode(Q.eq(a, b) & Q.positive(a) & Q.positive(b))
    euf, _ = EUFTheorySolver.from_encoded_cnf(enc)
    ab, pa, pb = (enc.encoding[p] for p in
                  (Q.eq(a, b), Q.positive(a), Q.positive(b)))
    euf.assert_lit(pa)
    euf.assert_lit(-pb)
    res = euf.assert_lit(ab)
    assert res[0] is False
    assert set(res[1]) == {-ab, -pa, pb}


def test_different_predicates_do_not_collide():
    # Q.positive and Q.negative are distinct heads, and two predicates that
    # are merely both false say nothing about each other.
    enc = _encode(Q.eq(a, b) & Q.positive(a) & Q.negative(b))
    euf, _ = EUFTheorySolver.from_encoded_cnf(enc)
    for prop in (Q.eq(a, b), Q.positive(a)):
        assert euf.assert_lit(enc.encoding[prop]) is None
    assert euf.assert_lit(-enc.encoding[Q.negative(b)]) is None
    assert euf.check()[0] is True


def test_predicate_congruence_end_to_end():
    assert satisfiable(Q.eq(a, b) & Q.positive(a) & ~Q.positive(b),
                       use_euf_theory=True) is False
    assert satisfiable(Q.eq(a, b) & Q.positive(a) & ~Q.negative(b),
                       use_euf_theory=True) is not False
    # nested inside a function application
    assert satisfiable(Q.eq(a, b) & Q.positive(f(a)) & ~Q.positive(f(b)),
                       use_euf_theory=True) is False


def test_assert_lit_reports_transitivity_conflict():
    enc = _encode(Q.eq(a, b) & Q.eq(b, c) & Q.ne(a, c))
    euf, _ = EUFTheorySolver.from_encoded_cnf(enc)
    ab, bc, ac = (enc.encoding[p] for p in (Q.eq(a, b), Q.eq(b, c), Q.ne(a, c)))
    assert euf.assert_lit(ab) is None
    assert euf.assert_lit(ac) is None
    res = euf.assert_lit(bc)
    assert res[0] is False
    assert set(res[1]) == {-ab, -bc, -ac}
    assert euf.check() == res


def test_conflict_uses_only_the_relevant_equalities():
    enc = _encode(Q.eq(a, b) & Q.eq(x, y) & Q.ne(a, b))
    euf, _ = EUFTheorySolver.from_encoded_cnf(enc)
    ab, xy, ne = (enc.encoding[p] for p in (Q.eq(a, b), Q.eq(x, y), Q.ne(a, b)))
    euf.assert_lit(xy)
    euf.assert_lit(ab)
    res = euf.assert_lit(ne)
    # the x = y branch is disjoint and must not show up in the clause
    assert set(res[1]) == {-ab, -ne}


def test_conflict_from_congruence():
    # a = b alone forces f(a) = f(b); the disequality is what breaks.
    enc = _encode(Q.eq(a, b) & Q.ne(f(a), f(b)))
    euf, _ = EUFTheorySolver.from_encoded_cnf(enc)
    ab, ne = enc.encoding[Q.eq(a, b)], enc.encoding[Q.ne(f(a), f(b))]
    assert euf.assert_lit(ne) is None
    res = euf.assert_lit(ab)
    assert set(res[1]) == {-ab, -ne}


def test_negated_equality_atom_is_a_disequality():
    # ~Q.eq(a, b) reaches the solver as a negative literal, not a Q.ne atom.
    enc = _encode(Q.eq(a, b) & (Q.eq(b, c) | Q.eq(c, d)))
    euf, _ = EUFTheorySolver.from_encoded_cnf(enc)
    ab, bc = enc.encoding[Q.eq(a, b)], enc.encoding[Q.eq(b, c)]
    euf.assert_lit(ab)
    euf.assert_lit(-bc)
    assert euf.check()[0] is True
    # asserting a = c on top makes b = c derivable, contradicting ~(b = c)
    enc2 = _encode(Q.eq(a, b) & Q.eq(a, c) & Q.eq(b, c))
    euf2, _ = EUFTheorySolver.from_encoded_cnf(enc2)
    ab2, ac2, bc2 = (enc2.encoding[p] for p in
                     (Q.eq(a, b), Q.eq(a, c), Q.eq(b, c)))
    euf2.assert_lit(ab2)
    euf2.assert_lit(ac2)
    res = euf2.assert_lit(-bc2)
    assert set(res[1]) == {-ab2, -ac2, bc2}


def test_check_and_reset():
    enc = _encode(Q.eq(a, b) & Q.ne(a, b))
    euf, _ = EUFTheorySolver.from_encoded_cnf(enc)
    ab, ne = enc.encoding[Q.eq(a, b)], enc.encoding[Q.ne(a, b)]
    assert euf.check()[0] is True
    euf.assert_lit(ab)
    is_sat, classes = euf.check()
    assert is_sat is True
    assert classes[a] == classes[b]
    assert euf.assert_lit(ne)[0] is False
    euf.reset()
    assert euf.check() == (True, {})
    assert euf.assert_lit(ne) is None
    assert euf.check()[0] is True


def test_satisfiable_end_to_end():
    assert satisfiable(Q.eq(a, b) & Q.eq(b, c) & Q.ne(a, c),
                       use_euf_theory=True) is False
    assert satisfiable(Q.eq(a, b) & Q.ne(f(a), f(b)),
                       use_euf_theory=True) is False
    assert satisfiable(Q.eq(a, b) & Q.ne(f(a), g(b)),
                       use_euf_theory=True) is not False
    # the solver has to pick the disjunct that keeps the closure consistent
    model = satisfiable((Q.eq(a, b) | Q.eq(a, c)) & Q.ne(a, b) & Q.eq(c, d),
                        use_euf_theory=True)
    assert model[Q.eq(a, c)] is True


def test_satisfiable_without_the_theory_misses_the_conflict():
    # Without EUF the equalities are just opaque booleans.
    prop = Q.eq(a, b) & Q.eq(b, c) & Q.ne(a, c)
    assert satisfiable(prop) is not False
    assert satisfiable(prop, use_euf_theory=True) is False


def test_dpll_satisfiable_flag():
    enc = _encode(Q.eq(a, b) & Q.ne(a, b))
    assert dpll_satisfiable(enc, use_euf_theory=True) is False


def test_theories_cannot_be_combined():
    raises(ValueError, lambda: satisfiable(Q.eq(a, b), use_lra_theory=True,
                                           use_euf_theory=True))
    raises(ValueError, lambda: satisfiable(Q.eq(a, b), algorithm="dpll",
                                           use_euf_theory=True))


def _euf_consistent(assignment):
    """Brute-force reference: is this assignment of equality atoms possible?"""
    cc = EUFCongruenceClosure([])
    for (lhs, rhs), value in assignment.items():
        if value:
            cc.merge(lhs, rhs)
    return all(not cc.are_congruent(lhs, rhs)
               for (lhs, rhs), value in assignment.items() if not value)


def test_random_formulas_match_brute_force():
    # Enumerate every assignment of the equality atoms, keep the ones the
    # closure allows, and see whether any of them satisfies the formula.  A
    # conflict clause with the wrong polarity shows up here as a wrong UNSAT.
    terms = [a, b, c, f(a), f(b), g(a, b)]
    for _ in range(25):
        pairs = list(dict.fromkeys(tuple(sample(terms, 2))
                                   for _ in range(randint(2, 4))))
        clauses = []
        for _ in range(randint(1, 4)):
            clauses.append(Or(*[choice([Q.eq, Q.ne])(lhs, rhs)
                                for lhs, rhs in sample(pairs, randint(1, len(pairs)))]))
        prop = And(*clauses)

        expected = False
        for bits in product([True, False], repeat=len(pairs)):
            assignment = dict(zip(pairs, bits))
            if not _euf_consistent(assignment):
                continue
            subs = {}
            for (lhs, rhs), value in assignment.items():
                subs[Q.eq(lhs, rhs)] = value
                subs[Q.ne(lhs, rhs)] = not value
            if prop.subs(subs) == True:
                expected = True
                break
        assert (satisfiable(prop, use_euf_theory=True) is not False) is expected
