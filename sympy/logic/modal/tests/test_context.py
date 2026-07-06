import pytest
import warnings
from sympy.core.symbol import Symbol
from sympy.logic.boolalg import Implies, Or, Not
from sympy.logic.modal.frames import KripkeFrame, Axiom
from sympy.logic.modal.context import ProofContext, Strategy
from sympy.logic.modal.errors import ProofFailure
from sympy.logic.modal.kernel import ProofTerm

def test_context_assume_and_discharge():
    ctx = ProofContext(KripkeFrame.K())
    p = Symbol('p')

    hyp = ctx.assume(p)
    assert hyp in ctx.hypotheses
    assert hyp.source == 'hypothesis'

    # Fake proof of something from p
    q = Symbol('q')
    pt_q = ProofTerm(q, hypotheses=[p])

    impl = ctx.discharge(hyp, pt_q)
    assert hyp not in ctx.hypotheses
    assert impl.formula == Implies(p, q)
    assert p not in impl.hypotheses

def test_context_save_restore():
    ctx = ProofContext(KripkeFrame.K())
    p = Symbol('p')

    ctx.save()
    ctx.assume(p)
    assert len(ctx.hypotheses) == 1

    ctx.restore()
    assert len(ctx.hypotheses) == 0

def test_classical_warnings():
    ctx = ProofContext(KripkeFrame.K())
    p = Symbol('p')

    # Law of excluded middle
    lem = Or(p, Not(p))
    with pytest.warns(UserWarning, match="Law of Excluded Middle"):
        ctx.assume(lem)

    # Double negation elimination
    nn = Not(Not(p, evaluate=False), evaluate=False)
    dne = Implies(nn, p, evaluate=False)
    with pytest.warns(UserWarning, match="Double Negation Elimination"):
        ctx.assume(dne)

def test_context_prove_lob():
    # Löb's theorem should be provable in GL
    ctx = ProofContext(frame=KripkeFrame.GL())
    p = Symbol('p')

    lob = Implies(ctx.Box(Implies(ctx.Box(p), p)), ctx.Box(p))

    result = ctx.prove(lob)
    assert isinstance(result, ProofTerm)
    assert result.formula == lob
    assert ctx.kernel.check_term(result, lob)

    # Same theorem should fail in S4
    ctx_s4 = ProofContext(frame=KripkeFrame.S4())
    failure = ctx_s4.prove(lob)
    assert isinstance(failure, ProofFailure)
    assert Axiom.Lob in failure.missing_axioms
def test_classical_logic_opt_in():
    from sympy.logic.boolalg import Or, Not
    from sympy.logic.modal.context import ProofContext
    from sympy.logic.modal.frames import KripkeFrame
    from sympy.core.symbol import Symbol

    p = Symbol('p')
    lem = Or(p, Not(p))

    # Normally it should warn
    ctx = ProofContext(frame=KripkeFrame.K())
    import warnings
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        ctx.assume(lem)
        assert len(w) > 0
        assert "Law of Excluded Middle" in str(w[-1].message)

    # With allow_classical=True, it shouldn't warn
    ctx_classical = ProofContext(frame=KripkeFrame.K(), allow_classical=True)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        ctx_classical.assume(lem)
        # Assuming no other warnings happen
        assert not any("Law of Excluded Middle" in str(warn.message) for warn in w)
def test_tactic_rewrite():
    from sympy.logic.modal.context import ProofContext
    from sympy.logic.modal.frames import KripkeFrame
    from sympy.core.symbol import Symbol
    from sympy.logic.boolalg import Implies, Equivalent
    from sympy.logic.modal.kernel import ProofTerm
    import pytest

    ctx = ProofContext(frame=KripkeFrame.K())
    p = Symbol('p')
    q = Symbol('q')

    pt = ctx.assume(Implies(p, q))

    # Try rewriting without a valid equivalence proof
    invalid_eq_proof = ProofTerm(Implies(p, q)) # Not an equivalence
    with pytest.raises(ValueError):
        ctx.tactic_rewrite(pt, p, q, invalid_eq_proof)

    # Rewrite with a valid equivalence proof
    valid_eq_proof = ctx.assume(Equivalent(p, q))
    pt_rewritten = ctx.tactic_rewrite(pt, p, q, valid_eq_proof)
    assert pt_rewritten.formula == Implies(q, q)

def test_smt_integration():
    from sympy.logic.modal.context import ProofContext
    from sympy.logic.modal.frames import KripkeFrame
    from sympy.core.symbol import Symbol
    from sympy.logic.boolalg import Implies, And, Or, Not

    # Needs classical logic for SMT integration to trigger on propositional formulas
    ctx = ProofContext(frame=KripkeFrame.K(), allow_classical=True)
    p = Symbol('p')
    q = Symbol('q')

    # A complex tautology that SMT can solve instantly but backward/forward might struggle with
    # (p -> q) v (q -> p)
    tautology = Or(Implies(p, q), Implies(q, p))
    proof = ctx.prove(tautology)

    from sympy.logic.modal.kernel import ProofTerm
    assert isinstance(proof, ProofTerm)
    assert proof.derivation[0] == "SMT_Solver"
    assert proof.formula == tautology

def test_classical_injection():
    from sympy.logic.modal.context import ProofContext
    from sympy.logic.modal.frames import KripkeFrame
    from sympy.core.symbol import Symbol
    from sympy.logic.boolalg import Implies, And, Or, Not

    ctx = ProofContext(frame=KripkeFrame.K(), allow_classical=True)
    p = Symbol('p')

    # Try to prove LEM
    lem = Or(p, Not(p))
    proof = ctx.prove(lem)

    from sympy.logic.modal.kernel import ProofTerm
    assert isinstance(proof, ProofTerm)
    # The SMT solver might pick this up first depending on the implementation details,
    # but the point is we can prove classical axioms directly if allow_classical=True.
