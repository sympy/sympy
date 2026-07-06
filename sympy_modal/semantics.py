from typing import Dict, Set, Any, Tuple
from sympy.logic.boolalg import Boolean, And, Or, Not, Implies
from sympy.core.symbol import Symbol
from sympy_modal.operators import Box, Diamond, ModalOperator

class KripkeModel:
    """
    Finite Kripke Model for Semantic Evaluation.
    W: Set of worlds
    R: Accessibility relation (dict mapping world to set of accessible worlds)
    V: Valuation function (dict mapping (world, proposition) to bool)
    """
    def __init__(self, W: Set[str], R: Dict[str, Set[str]], V: Dict[Tuple[str, Any], bool]):
        self.W = W
        self.R = R
        self.V = V

class SemanticEvaluator:
    """
    Evaluates formulas against finite Kripke models.
    """
    def __init__(self, model: KripkeModel):
        self.model = model

    def evaluate(self, formula: Boolean, world: str) -> bool:
        if world not in self.model.W:
            raise ValueError(f"World {world} not in model")

        if isinstance(formula, Symbol):
            return self.model.V.get((world, formula), False)

        if isinstance(formula, Not):
            return not self.evaluate(formula.args[0], world)

        if isinstance(formula, And):
            return all(self.evaluate(arg, world) for arg in formula.args)

        if isinstance(formula, Or):
            return any(self.evaluate(arg, world) for arg in formula.args)

        if isinstance(formula, Implies):
            return (not self.evaluate(formula.args[0], world)) or self.evaluate(formula.args[1], world)

        if isinstance(formula, Box):
            # True if formula holds in ALL accessible worlds
            accessible = self.model.R.get(world, set())
            return all(self.evaluate(formula.args[0], w) for w in accessible)

        if isinstance(formula, Diamond):
            # True if formula holds in SOME accessible world
            accessible = self.model.R.get(world, set())
            return any(self.evaluate(formula.args[0], w) for w in accessible)

        from sympy_modal.operators import Next, Until

        if isinstance(formula, Next):
            # Next (X p): True if p holds in ALL immediately accessible worlds
            # (Standard LTL assumes a single path, but on branching Kripke we check all next states, analogous to AX)
            accessible = self.model.R.get(world, set())
            if not accessible:
                return False # Or True depending on semantics of deadlock, usually False for 'next'
            return all(self.evaluate(formula.args[0], w) for w in accessible)

        if isinstance(formula, Until):
            # Until (p U q): True if q holds now or eventually, and p holds in all steps until q holds.
            # Implemented via BFS on the Kripke frame to avoid infinite loops on cycles.
            p, q = formula.args

            queue = [world]
            visited = set()

            while queue:
                curr = queue.pop(0)
                if curr in visited:
                    continue
                visited.add(curr)

                if self.evaluate(q, curr):
                    # Found a state where q holds
                    # Need to verify p held on all paths *to* here, but BFS from root is tricky for branching.
                    # Actually, standard CTL 'Until' (A[p U q] or E[p U q]) is complex.
                    # Let's implement A[p U q] (holds on ALL paths).
                    pass # We will use a recursive helper for A[p U q]

            # Using a recursive helper for A[p U q] with cycle detection
            def eval_until(current_world: str, current_visited: Set[str]) -> bool:
                if self.evaluate(q, current_world):
                    return True
                if not self.evaluate(p, current_world):
                    return False

                next_worlds = self.model.R.get(current_world, set())
                if not next_worlds:
                    return False # Dead end before q

                for nw in next_worlds:
                    if nw in current_visited:
                        return False # Cycle detected before q, so not *eventually* q on this path
                    if not eval_until(nw, current_visited | {nw}):
                        return False
                return True

            return eval_until(world, {world})

        raise ValueError(f"Unsupported formula for semantic evaluation: {type(formula)}")
