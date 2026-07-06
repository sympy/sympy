"""
Layer 2 (Frames): Kripke Frame Semantics.
"""

from typing import Dict, Set, List, Any
from enum import Enum

from sympy.logic.boolalg import Boolean, Implies, And, Or, Not
from sympy.logic.modal.operators import Box, Diamond
from sympy.logic.modal.errors import FrameViolationError

class Axiom(Enum):
    K = "K"         # □(A → B) → (□A → □B)
    T = "T"         # □A → A (reflexivity)
    B = "B"         # A → □◇A (symmetry)
    Four = "4"      # □A → □□A (transitivity)
    Five = "5"      # ◇A → □◇A (euclidean)
    Lob = "Löb"     # □(□A → A) → □A (converse well-foundedness)
    D = "D"         # □A → ◇A (seriality)

class KripkeFrame:
    """
    A Kripke frame governing inference.
    """
    def __init__(self, worlds: Set[str], accessibility: Dict[str, Set[str]], axioms: List[Axiom]):
        self.worlds = worlds
        self.accessibility = accessibility
        self.axioms = axioms

    @classmethod
    def K(cls) -> 'KripkeFrame':
        """Minimal normal modal logic."""
        return cls({"w"}, {"w": set()}, [Axiom.K])

    @classmethod
    def T(cls) -> 'KripkeFrame':
        """Reflexive (alethic)"""
        return cls({"w"}, {"w": {"w"}}, [Axiom.K, Axiom.T])

    @classmethod
    def S4(cls) -> 'KripkeFrame':
        """Reflexive and transitive"""
        return cls({"w"}, {"w": {"w"}}, [Axiom.K, Axiom.T, Axiom.Four])

    @classmethod
    def S5(cls) -> 'KripkeFrame':
        """Reflexive, transitive, symmetric (equivalence relation)"""
        return cls({"w"}, {"w": {"w"}}, [Axiom.K, Axiom.T, Axiom.Four, Axiom.Five, Axiom.B])

    @classmethod
    def GL(cls) -> 'KripkeFrame':
        """Gödel-Löb provability logic: transitive and converse well-founded."""
        # Represents an irreflexive transitive relation
        # (Converse well-foundedness is harder to represent simply with a single world,
        # but axioms suffice for inference checking).
        return cls({"w1", "w2"}, {"w1": {"w2"}, "w2": set()}, [Axiom.K, Axiom.Four, Axiom.Lob])

    @classmethod
    def D(cls) -> 'KripkeFrame':
        """Deontic: serial accessibility"""
        return cls({"w"}, {"w": {"w"}}, [Axiom.K, Axiom.D])

    @classmethod
    def K45(cls) -> 'KripkeFrame':
        """Epistemic: transitive and euclidean"""
        return cls({"w"}, {"w": {"w"}}, [Axiom.K, Axiom.Four, Axiom.Five])

    def _evaluate(self, formula: Boolean, world: str, valuation: Dict[Any, Set[str]]) -> bool:
        """
        Evaluate a formula at a given world in a given valuation.
        valuation maps atomic propositions to the set of worlds where they are true.
        """
        if isinstance(formula, Box):
            from sympy.logic.boolalg import Boolean
            inner = formula.args[0]
            if not isinstance(inner, Boolean):
                return False
            # True if inner is true in ALL accessible worlds
            for w in self.accessibility.get(world, set()):
                if not self._evaluate(inner, w, valuation):
                    return False
            return True

        elif isinstance(formula, Diamond):
            from sympy.logic.boolalg import Boolean
            inner = formula.args[0]
            if not isinstance(inner, Boolean):
                return False
            # True if inner is true in SOME accessible world
            for w in self.accessibility.get(world, set()):
                if self._evaluate(inner, w, valuation):
                    return True
            return False

        elif isinstance(formula, Implies):
            left, right = formula.args
            return (not self._evaluate(left, world, valuation)) or self._evaluate(right, world, valuation)

        elif isinstance(formula, And):
            return all(self._evaluate(arg, world, valuation) for arg in formula.args)

        elif isinstance(formula, Or):
            return any(self._evaluate(arg, world, valuation) for arg in formula.args)

        elif isinstance(formula, Not):
            return not self._evaluate(formula.args[0], world, valuation)

        else:
            # Atomic proposition or unbound variable
            return world in valuation.get(formula, set())

    def validates(self, formula: Boolean) -> bool:
        """
        Checks whether a formula is valid in all worlds of this frame
        under all possible valuations (for a simple finite model check).

        Since rigorous universal validation over all possible valuations for any frame
        is PSPACE-hard, we do a basic structural check based on the frame's axioms
        for the standard axioms, and fallback to evaluating over a small set of test valuations
        for the specific frame instances provided above.
        """
        # First, structural check for common axioms:
        if self._is_axiom_instance(formula):
            return True

        # Fallback to model checking on the explicitly defined worlds:
        # Generate all valuations for the atoms in the formula
        atoms = formula.atoms()

        # Limit to basic propositional/modal atoms
        # Exclude types/universes etc
        atoms = {a for a in atoms if isinstance(a, Boolean) and not isinstance(a, (Box, Diamond, Implies, And, Or, Not))}

        # If formula is already SymPy true/false boolean literal:
        from sympy.logic.boolalg import BooleanTrue, BooleanFalse
        if isinstance(formula, BooleanTrue):
            return True
        if isinstance(formula, BooleanFalse):
            return False

        import itertools
        # Each atom can be true in any subset of worlds
        power_set_worlds = []
        for i in range(len(self.worlds) + 1):
            for subset in itertools.combinations(self.worlds, i):
                power_set_worlds.append(set(subset))

        # Generator of all possible valuations
        valuations = []
        for val_tuple in itertools.product(power_set_worlds, repeat=len(atoms)):
            val_dict = dict(zip(atoms, val_tuple))
            valuations.append(val_dict)

        for val in valuations:
            for w in self.worlds:
                if not self._evaluate(formula, w, val):
                    return False
        return True

    def _is_axiom_instance(self, formula: Boolean) -> bool:
        """Heuristically check if formula matches an allowed axiom."""
        # This is a very simplified pattern matcher for testing.
        # A full kernel proof checker handles arbitrary derivations.
        if isinstance(formula, Implies):
            left, right = formula.args
            if isinstance(left, Box) and Axiom.T in self.axioms:
                if left.args[0] == right:
                    return True # Box(A) -> A

            if isinstance(left, Box) and isinstance(right, Box) and Axiom.Four in self.axioms:
                if isinstance(right.args[0], Box) and right.args[0].args[0] == left.args[0]:
                    return True # Box(A) -> Box(Box(A))

            if isinstance(left, Box) and isinstance(right, Box) and Axiom.Lob in self.axioms:
                # Box(Box(A) -> A) -> Box(A)
                if isinstance(left.args[0], Implies) and isinstance(left.args[0].args[0], Box):
                    if left.args[0].args[0].args[0] == right.args[0] and left.args[0].args[1] == right.args[0]:
                        return True

        from sympy.logic.modal.operators import ForAllPredicates
        if isinstance(formula, ForAllPredicates):
            # Verify the inner formula
            return self._is_axiom_instance(formula.formula)

        return False

    def is_valid_inference(self, premises: List[Boolean], conclusion: Boolean) -> bool:
        """
        Checks whether the inference is sound in this frame.
        """
        # A simple check: Construct (P1 AND P2 ...) -> C and check validity
        if not premises:
            return self.validates(conclusion)

        if len(premises) == 1:
            return self.validates(Implies(premises[0], conclusion))

        return self.validates(Implies(And(*premises), conclusion))
