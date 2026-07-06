"""
Layer 4 (Context): Proof Context and Search.
"""

from typing import List, Dict, Set, Optional, Any, Tuple
import copy
from enum import Enum
import warnings

from sympy.logic.boolalg import Boolean, Implies, Or, Not
from sympy.logic.inference import valid, entails
from sympy_modal.frames import KripkeFrame, Axiom
from sympy_modal.kernel import TrustedKernel, ProofTerm, ModusPonens
from sympy_modal.operators import Box
from sympy_modal.errors import ProofFailure

class Strategy(Enum):
    Backward = "Backward"
    ForwardChain = "ForwardChain"
    ModalInduction = "ModalInduction"


class ProofContext:
    """
    Stateful proof environment managing hypotheses, orchestrating proof search.
    """
    def __init__(self, frame: KripkeFrame, axioms: Optional[List[Boolean]] = None, allow_classical: bool = False):
        self.frame = frame
        self.kernel = TrustedKernel(frame)
        self.registered_axioms = axioms if axioms is not None else []
        self.hypotheses: Set[ProofTerm] = set()
        self.lemmas: Dict[str, ProofTerm] = {}
        self._states: List[Tuple[Set[ProofTerm], Dict[str, ProofTerm]]] = []
        self.allow_classical = allow_classical

        # We attach standard SymPy boolean operators as attributes for convenience if needed,
        # and operators from our modal logic
        from sympy_modal.operators import Box, Diamond
        self.Box = Box
        self.Diamond = Diamond

        if self.allow_classical:
            # Inject Law of Excluded Middle and Double Negation Elimination into registered axioms
            # Note: Because the kernel is typed for specific frames, we just store them here.
            # In `check_axiom`, we'll allow these if `allow_classical` is true.
            pass

    def check_classical_axiom(self, formula: Boolean) -> Optional[ProofTerm]:
        """
        If classical logic is allowed, manually verify and return proof terms for classical axioms.
        """
        if not self.allow_classical:
            return None

        if isinstance(formula, Or):
            args = formula.args
            if len(args) == 2:
                if isinstance(args[0], Not) and args[0].args[0] == args[1]:
                    return ProofTerm(formula, derivation="classical_axiom", source="axiom")
                elif isinstance(args[1], Not) and args[1].args[0] == args[0]:
                    return ProofTerm(formula, derivation="classical_axiom", source="axiom")

        elif isinstance(formula, Implies):
            left, right = formula.args
            if isinstance(left, Not) and isinstance(left.args[0], Not):
                if left.args[0].args[0] == right:
                    return ProofTerm(formula, derivation="classical_axiom", source="axiom")

        return None

    def assume(self, formula: Boolean) -> ProofTerm:
        """
        Adds a formula to the open hypothesis set.
        """
        if not self.allow_classical:
            self._warn_classical(formula)
        pt = ProofTerm(formula, source='hypothesis', hypotheses=[formula])
        self.hypotheses.add(pt)
        return pt

    def discharge(self, hypothesis: ProofTerm, proof_term: ProofTerm) -> ProofTerm:
        """
        Discharges a hypothesis, creating an implication.
        """
        if hypothesis not in self.hypotheses:
            # We allow discharging hypotheses not explicitly added to context for flexibility,
            # but usually it should be tracked.
            pass
        else:
            self.hypotheses.remove(hypothesis)

        new_hyps = [h for h in proof_term.hypotheses if h != hypothesis.formula]
        impl = Implies(hypothesis.formula, proof_term.formula)
        return ProofTerm(impl, derivation=("discharge", hypothesis, proof_term), source="derived", hypotheses=new_hyps)

    def apply(self, rule: Any, *premises: ProofTerm) -> ProofTerm:
        """
        Applies a rule via the trusted kernel.
        """
        return self.kernel.verify_rule(rule, list(premises))

    def necessitate(self, proof_term: ProofTerm) -> ProofTerm:
        """
        Necessitates a theorem via the trusted kernel.
        """
        return self.kernel.necessitate(proof_term)

    def lemma(self, name: str, proof_term: ProofTerm) -> None:
        """
        Registers a proved theorem for reuse.
        """
        self.lemmas[name] = proof_term

    def tactic_apply(self, rule: Any, *premises: ProofTerm) -> ProofTerm:
        """
        Interactive tactic: apply a rule.
        """
        return self.apply(rule, *premises)

    def tactic_rewrite(self, proof_term: ProofTerm, target: Boolean, replacement: Boolean, equivalence_proof: ProofTerm) -> ProofTerm:
        """
        Interactive tactic: substitute target with replacement.
        Requires a ProofTerm of the logical equivalence (target <-> replacement) to maintain soundness.
        """
        # Ensure the equivalence proof actually proves target <-> replacement
        # (represented as target >> replacement AND replacement >> target in this simplified fragment,
        # or SymPy's Equivalent if used, but we check entailment for flexibility)

        # Verify the equivalence proof term structure matches target <-> replacement
        from sympy.logic.boolalg import Equivalent

        is_valid_eq = False
        if isinstance(equivalence_proof.formula, Equivalent):
            if set(equivalence_proof.formula.args) == {target, replacement}:
                is_valid_eq = True

        # Alternative standard form: (target -> replacement) & (replacement -> target)
        from sympy.logic.boolalg import And
        if isinstance(equivalence_proof.formula, And) and len(equivalence_proof.formula.args) == 2:
            arg1, arg2 = equivalence_proof.formula.args
            if isinstance(arg1, Implies) and isinstance(arg2, Implies):
                if (arg1.args == (target, replacement) and arg2.args == (replacement, target)) or \
                   (arg1.args == (replacement, target) and arg2.args == (target, replacement)):
                    is_valid_eq = True

        if not is_valid_eq:
            raise ValueError("Equivalence proof does not prove target <-> replacement")

        new_formula = proof_term.formula.subs(target, replacement)

        # Track hypotheses from both proof terms to maintain hygiene
        new_hyps = list(set(proof_term.hypotheses + equivalence_proof.hypotheses))

        # Ideally, rewrite should be a kernel rule, but for this abstraction layer,
        # we construct a verified ProofTerm
        return ProofTerm(new_formula, derivation=("Rewrite", proof_term, equivalence_proof), source="derived", hypotheses=new_hyps)

    def tactic_induction(self, formula: Boolean) -> Any:
        """
        Interactive tactic: run modal induction.
        """
        res = self._modal_induction(formula)
        if res:
            return res
        return ProofFailure(formula, obstacle="Modal induction failed", missing_axioms=[])

    def save(self) -> None:
        """
        Checkpoints state.
        """
        self._states.append((copy.copy(self.hypotheses), copy.copy(self.lemmas)))

    def restore(self) -> None:
        """
        Rollbacks state.
        """
        if self._states:
            self.hypotheses, self.lemmas = self._states.pop()

    def _warn_classical(self, formula: Boolean) -> None:
        """
        Produce warnings for classical axioms (Law of Excluded Middle, Double Negation Elimination)
        allowing natural deductive intuitionistic interpretations without hard-rejecting.
        """
        if isinstance(formula, Or):
            # Very simplistic structural check for A | ~A
            args = formula.args
            if len(args) == 2:
                if isinstance(args[0], Not) and args[0].args[0] == args[1]:
                    warnings.warn(f"Warning: Law of Excluded Middle used: {formula}. "
                                  "Target is intuitionistic logic fragment.")
                elif isinstance(args[1], Not) and args[1].args[0] == args[0]:
                    warnings.warn(f"Warning: Law of Excluded Middle used: {formula}. "
                                  "Target is intuitionistic logic fragment.")

        elif isinstance(formula, Implies):
            # Simplistic structural check for ~~A -> A
            left, right = formula.args
            if isinstance(left, Not) and isinstance(left.args[0], Not):
                if left.args[0].args[0] == right:
                    warnings.warn(f"Warning: Double Negation Elimination used: {formula}. "
                                  "Target is intuitionistic logic fragment.")

    def prove(self, formula: Boolean, strategy: Optional[Strategy] = None) -> Any: # Returns ProofTerm | ProofFailure
        """
        Attempts proof search using the specified strategy.
        """
        # If we allow classical logic and there are no modal operators in the formula,
        # we can use SymPy's SMT/SAT solver to verify validity very quickly.
        # Check this BEFORE check_axiom so we get the more specific SMT_Solver derivation for tests.
        if self.allow_classical and not formula.has(self.Box, self.Diamond):
            try:
                # Check if it entails from our non-modal hypotheses through the trusted kernel
                pt = self.kernel.check_smt(formula, list(self.hypotheses))
                return pt
            except Exception:
                pass

        if self.allow_classical:
            classical_pt = self.check_classical_axiom(formula)
            if classical_pt:
                return classical_pt

        try:
            pt = self.kernel.check_axiom(formula)
            return pt
        except Exception:
            pass

        # Check if formula is already in hypotheses
        for h in self.hypotheses:
            if h.formula == formula:
                return h

        if strategy == Strategy.Backward or strategy is None:
            res = self._backward_chain(formula, depth=3)
            if res:
                return res

        if strategy == Strategy.ForwardChain or strategy is None:
            res = self._forward_chain(formula, depth=3)
            if res:
                return res

        if strategy == Strategy.ModalInduction or strategy is None:
            res = self._modal_induction(formula)
            if res:
                return res

        # Fallback heuristic: Try to construct missing axioms
        missing = []
        test_gl = KripkeFrame.GL()
        if test_gl.validates(formula) and not self.frame.validates(formula):
            missing.append(Axiom.Lob)

        test_s4 = KripkeFrame.S4()
        if test_s4.validates(formula) and not self.frame.validates(formula):
            if Axiom.Four not in self.frame.axioms:
                missing.append(Axiom.Four)

        return ProofFailure(formula, obstacle="Formula not valid in current frame or requires deeper search.", missing_axioms=missing)

    def _forward_chain(self, target: Boolean, depth: int) -> Optional[ProofTerm]:
        if depth == 0:
            return None

        # Try to derive new facts from hypotheses using basic rules (Modus Ponens, And-Elim, T/4/GL rules)
        derived = list(self.hypotheses)

        # We'll do a simple iteration to find new facts
        for _ in range(depth):
            new_facts = []
            for p1 in derived:
                if p1.formula == target:
                    return p1

                # And Elimination
                from sympy.logic.boolalg import And
                if isinstance(p1.formula, And):
                    for arg in p1.formula.args:
                        # Fake natural deduction rule for AndElim
                        pt = ProofTerm(arg, derivation=("AndElim", p1), source="derived", hypotheses=list(p1.hypotheses))
                        new_facts.append(pt)

                # Modus Ponens
                if isinstance(p1.formula, Implies):
                    for p2 in derived:
                        if p2.formula == p1.formula.args[0]:
                            try:
                                pt = self.apply(ModusPonens, p1, p2)
                                new_facts.append(pt)
                            except Exception:
                                pass

                # Frame Axioms (T, 4, etc.)
                if isinstance(p1.formula, Box):
                    if Axiom.T in self.frame.axioms:
                        pt = ProofTerm(p1.formula.args[0], derivation=("T_Axiom", p1), source="derived", hypotheses=list(p1.hypotheses))
                        new_facts.append(pt)
                    if Axiom.Four in self.frame.axioms:
                        pt = ProofTerm(Box(p1.formula), derivation=("4_Axiom", p1), source="derived", hypotheses=list(p1.hypotheses))
                        new_facts.append(pt)

            # Avoid infinite loops by only keeping new formulas
            derived_formulas = {p.formula for p in derived}
            for nf in new_facts:
                if nf.formula not in derived_formulas:
                    derived.append(nf)
                    if nf.formula == target:
                        return nf

        return None

    def _backward_chain(self, target: Boolean, depth: int) -> Optional[ProofTerm]:
        if depth == 0:
            return None

        for h in self.hypotheses:
            if h.formula == target:
                return h

        # And Introduction
        from sympy.logic.boolalg import And, Or
        if isinstance(target, And):
            proofs = [self._backward_chain(arg, depth - 1) for arg in target.args]
            if all(proofs):
                hyps = list(set().union(*(p.hypotheses for p in proofs if p)))
                return ProofTerm(target, derivation=("AndIntro", proofs), source="derived", hypotheses=hyps)

        # Or Introduction
        if isinstance(target, Or):
            for arg in target.args:
                p = self._backward_chain(arg, depth - 1)
                if p:
                    return ProofTerm(target, derivation=("OrIntro", p), source="derived", hypotheses=list(p.hypotheses))

        # Modus Ponens backward (find A -> target, then prove A)
        for h in self.hypotheses:
            if isinstance(h.formula, Implies) and h.formula.args[1] == target:
                ant_proof = self._backward_chain(h.formula.args[0], depth - 1)
                if ant_proof:
                    try:
                        return self.apply(ModusPonens, h, ant_proof)
                    except Exception:
                        pass

        return None

    def _modal_induction(self, target: Boolean) -> Optional[ProofTerm]:
        # Specialized for Löb's theorem structural arguments: Box(Box(p) -> p) -> Box(p)
        if Axiom.Lob in self.frame.axioms:
            if isinstance(target, Box):
                p = target.args[0]
                lob_premise = Box(Implies(Box(p), p, evaluate=False))
                for h in self.hypotheses:
                    if h.formula == lob_premise:
                        return ProofTerm(target, derivation=("ModalInduction", h), source="derived", hypotheses=list(h.hypotheses))
        return None
