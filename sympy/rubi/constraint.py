from matchpy import Constraint, substitute
from sympy import sympify

class cons(Constraint):
    def __init__(self, expr, vars):
        self.expr = expr
        self.vars = frozenset(v.name for v in vars)

    def __call__(self, substitution):
        if isinstance(self.expr, bool): # handle rules without constraints
            return self.expr
        sub = substitute(self.expr, substitution)
        try:
            return sympify(str(sub))
        except:
            #print(('Unable to sympify: {}').format(sub))
            return False

    @property
    def variables(self):
        return self.vars

    def with_renamed_vars(self, renaming):
        if isinstance(self.expr, bool):
            copy = cons(self.expr, [])
        else:
            copy = cons(self.expr.with_renamed_vars(renaming), [])
        copy.vars = frozenset(renaming.get(v, v) for v in self.vars)
        return copy

    def __eq__(self, other):
        return isinstance(other, cons) and other.vars == self.vars and other.expr == self.expr

    def __hash__(self):
        return hash((self.vars, self.expr))
