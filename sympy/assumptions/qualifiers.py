from sympy import satisfiable, Expr, And, sympify, Basic
from sympy.logic.boolalg import Boolean


class ContradictionError(Exception):
    def __init__(self, contradiction, *args, **kwargs):
        """Raised when contradictory assumptions are encountered"""
        super().__init__(*args, **kwargs)
        self.contradiction = contradiction

def qualified(value,axiom=True):
    return Qualified(value,axiom)

def unqualified(value):
    return Qualified(value).unqualified

Basic.qualified = property(qualified)
Basic.unqualified = property(unqualified)

class Qualified(Expr, Boolean):
    def __new__(cls, term, axiom=True):
        new_term = sympify(term)
        new_axiom = sympify(axiom)

        def replacer(*args):
            return args[0]

        (new_term, replacements) = new_term.replace(Qualified, replacer, map=True)
        new_axiom = And(new_axiom, *[k.axiom for k in replacements.keys()])

        if not satisfiable(new_axiom):
            raise ContradictionError(new_axiom)
        else:
            obj = Expr.__new__(cls, new_term, new_axiom)
        return obj

    @property
    def unqualified(self):
        return self.args[0]

    @property
    def axiom(self):
        return self.args[1]

    def __eq__(self, other):
        if not isinstance(other, Qualified):
            return False
        return self.axiom == other.axiom and self.unqualified == other.unqualified

    def __hash__(self):
        return hash(super())

    def _hashable_content(self):
        return self.axiom, self.unqualified

