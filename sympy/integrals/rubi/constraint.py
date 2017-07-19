from sympy.external import import_module
matchpy = import_module("matchpy")

if matchpy:
    Constraint, substitute = matchpy.Constraint, matchpy.substitute
else:
    raise ImportError('MatchPy could not be imported')

from sympy.logic.boolalg import BooleanTrue
from sympy.integrals.rubi.matchpy2sympy import matchpy2sympy

class cons(Constraint):
    def __init__(self, expr, vars):
        self.expr = expr
        self.vars = frozenset(v.name for v in vars)

    def __call__(self, substitution):

        if isinstance(self.expr, bool): #handle rules without constraints
            return self.expr

        sub = substitute(self.expr, substitution)
        try:
            res = matchpy2sympy(sub)
        except:
            return False

        if isinstance(res, BooleanTrue) or res == True:
            return True
        else:
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

    def __str__(self):
        return str(self.expr)

class FreeQ(Constraint):
    def __init__(self, expr, var):
        self.expr = expr if isinstance(expr, str) else expr.name
        self.var = var if isinstance(var, str) else var.name

    def __call__(self, substitution):
        return substitution[self.var] not in substitution[self.expr]

    @property
    def variables(self):
        return frozenset([self.expr, self.var])

    def with_renamed_vars(self, renaming):
        return FreeQ(
            renaming.get(self.expr, self.expr),
            renaming.get(self.var, self.var)
        )

    def __eq__(self, other):
        return isinstance(other, FreeQ) and other.expr == self.expr and other.var == self.var

    def __hash__(self):
        return hash((self.expr, self.var))

    def __str__(self):
        return str('FreeQ({}, {})'.format(self.expr, self.var))
