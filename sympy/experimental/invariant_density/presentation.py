from sympy import Eq, Function
from collections import deque


def _apply_axiom(expr, ax):
    """
    Apply a single axiom using Wild-based pattern matching.
    Returns a set of rewritten expressions.
    """
    results = set()

    # lhs -> rhs
    m = expr.match(ax.lhs)
    if m is not None:
        results.add(ax.rhs.subs(m))

    # rhs -> lhs (symmetry)
    m = expr.match(ax.rhs)
    if m is not None:
        results.add(ax.lhs.subs(m))

    return results


class Presentation:
    """
    Finite equational presentation P = (G, R)
    with bounded derivability.
    """

    def __init__(self, generators, relations, associative_ops=None):
        self.generators = list(generators)
        self.relations = list(relations)
        self.associative_ops = set(associative_ops or [])

    # --------------------------------------------------
    # Associativity handling (canonical right-association)
    # --------------------------------------------------

    def _right_associate(self, expr, op):
        if not expr.args or expr.func != op:
            return expr

        a, b = expr.args
        if a.func == op:
            x, y = a.args
            return self._right_associate(op(x, op(y, b)), op)

        return op(a, b)

    # --------------------------------------------------
    # Bounded normalization (NF_R^m)
    # --------------------------------------------------

    def normalize(self, expr, axioms=None, max_steps=15):
        if axioms is None:
            axioms = self.relations

        seen = set()
        current = expr

        for _ in range(max_steps):
            if current in seen:
                break
            seen.add(current)

            new = current

            # apply associativity canonically
            for op in self.associative_ops:
                new = self._right_associate(new, op)

            # attempt one-step rewrites
            for ax in axioms:
                rewrites = _apply_axiom(new, ax)
                if rewrites:
                    new = next(iter(rewrites))
                    break

            if new == current:
                break

            current = new

        return current

    # --------------------------------------------------
    # Bounded derivability R ⊢_{n,m} t = u
    # --------------------------------------------------

    def derive(self, eq, axioms=None, max_depth=5, max_steps=15):
        if axioms is None:
            axioms = self.relations

        lhs = self.normalize(eq.lhs, axioms, max_steps)
        rhs = self.normalize(eq.rhs, axioms, max_steps)

        if lhs == rhs:
            return True

        frontier = deque([lhs])
        seen = {lhs}

        for _ in range(max_depth):
            next_frontier = deque()

            while frontier:
                t = frontier.popleft()

                for ax in axioms:
                    for u0 in _apply_axiom(t, ax):
                        u = self.normalize(u0, axioms, max_steps)

                        if u == rhs:
                            return True

                        if u not in seen:
                            seen.add(u)
                            next_frontier.append(u)

            if not next_frontier:
                break

            frontier = next_frontier

        return False

    # --------------------------------------------------
    # Redundancy checking
    # --------------------------------------------------

    def is_redundant(self, index, max_depth=5, max_steps=15):
        target = self.relations[index]
        others = self.relations[:index] + self.relations[index + 1 :]

        return self.derive(
            target,
            axioms=others,
            max_depth=max_depth,
            max_steps=max_steps,
        )

    # --------------------------------------------------
    # Greedy reduction
    # --------------------------------------------------

    def reduce(self, max_depth=5, max_steps=15):
        rels = list(self.relations)

        for i in range(len(rels) - 1, -1, -1):
            test = Presentation(
                self.generators, rels, self.associative_ops
            )
            if test.is_redundant(
                i,
                max_depth=max_depth,
                max_steps=max_steps,
            ):
                del rels[i]

        return Presentation(self.generators, rels, self.associative_ops)
