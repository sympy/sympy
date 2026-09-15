from __future__ import annotations
import pytest
from sympy.core.symbol import Symbol
from sympy.logic.boolalg import Implies
from sympy.logic.modal.frames import KripkeFrame
from sympy.logic.modal.kernel import TrustedKernel, ProofTerm, ModusPonens
from sympy.logic.modal.operators import Box
from sympy.logic.modal.errors import NotAnAxiomError, InvalidInferenceError, NecessitationError

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
def test_export_lean_coq():
    from sympy.core.symbol import Symbol
    from sympy.logic.boolalg import Implies
    from sympy.logic.modal.kernel import ProofTerm

    p = Symbol('p')
    q = Symbol('q')
    formula = Implies(p, q)
    pt = ProofTerm(formula)

    lean_export = pt.export_lean()
    assert "theorem my_thm" in lean_export
    assert "->" in lean_export or "Implies" not in lean_export
    assert "sorry" in lean_export

    coq_export = pt.export_coq()
    assert "Theorem my_thm" in coq_export
    assert "->" in coq_export or "Implies" not in coq_export
    assert "admit" in coq_export
def test_robust_export():
    from sympy.core.symbol import Symbol
    from sympy.logic.boolalg import Implies, Not
    from sympy.logic.modal.kernel import ProofTerm

    # Create a symbol with "Not" in the name to test collision
    NotASymbol = Symbol('NotASymbol')
    formula = Implies(Not(NotASymbol), NotASymbol)
    pt = ProofTerm(formula)

    lean_export = pt.export_lean()
    assert "¬NotASymbol" in lean_export

    coq_export = pt.export_coq()
    assert "~NotASymbol" in coq_export
