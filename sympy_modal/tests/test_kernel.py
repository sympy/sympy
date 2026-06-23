import pytest
from sympy.core.symbol import Symbol
from sympy.logic.boolalg import Implies
from sympy_modal.frames import KripkeFrame
from sympy_modal.kernel import TrustedKernel, ProofTerm, ModusPonens
from sympy_modal.operators import Box
from sympy_modal.errors import NotAnAxiomError, InvalidInferenceError, NecessitationError

def test_kernel_check_axiom():
    kernel = TrustedKernel(frame=KripkeFrame.GL())
    p = Symbol('p')
    # Tautology P -> P
    thm = kernel.check_axiom(Implies(p, p))
    assert thm.formula == Implies(p, p)

    with pytest.raises(NotAnAxiomError):
        kernel.check_axiom(p)

def test_kernel_verify_rule_modus_ponens():
    kernel = TrustedKernel(frame=KripkeFrame.GL())
    p = Symbol('p')
    q = Symbol('q')

    # Fake proofs for testing
    ab = ProofTerm(Implies(p, q))
    a = ProofTerm(p)

    b = kernel.verify_rule(ModusPonens, [ab, a])
    assert b.formula == q

    # order independent
    b2 = kernel.verify_rule(ModusPonens, [a, ab])
    assert b2.formula == q

    with pytest.raises(InvalidInferenceError):
        kernel.verify_rule(ModusPonens, [a, a])

    with pytest.raises(InvalidInferenceError):
        kernel.verify_rule("SomeOtherRule", [a, ab])

def test_kernel_necessitate():
    kernel = TrustedKernel(frame=KripkeFrame.GL())
    p = Symbol('p')

    thm = kernel.check_axiom(Implies(p, p))
    box = kernel.necessitate(thm)
    assert isinstance(box.formula, Box)
    assert box.formula.args[0] == thm.formula

    hyp = ProofTerm(p, source='hypothesis', hypotheses=[p])
    with pytest.raises(NecessitationError):
        kernel.necessitate(hyp)

def test_kernel_check_term():
    kernel = TrustedKernel(frame=KripkeFrame.GL())
    p = Symbol('p')
    thm = kernel.check_axiom(Implies(p, p))

    assert kernel.check_term(thm, Implies(p, p))
    assert not kernel.check_term(thm, p)
