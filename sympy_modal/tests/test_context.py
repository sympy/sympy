import pytest
import warnings
from sympy.core.symbol import Symbol
from sympy.logic.boolalg import Implies, Or, Not
from sympy_modal.frames import KripkeFrame, Axiom
from sympy_modal.context import ProofContext, Strategy
from sympy_modal.errors import ProofFailure
from sympy_modal.kernel import ProofTerm

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
