from __future__ import annotations

from sympy import Expr, Symbol, cancel, S

class KarrElement:
    """Represents an element in a Karr difference field."""
    def __init__(self, expr: Expr, field: KarrField):
        self.expr = cancel(expr)
        self.field = field

    def shift(self) -> KarrElement:
        """Applies the difference field automorphism (shift operator) sigma."""
        subs_dict = {}
        for t, info in self.field.generators.items():
            g_type, param = info
            if g_type == 'base':
                subs_dict[t] = t + 1
            elif g_type == 'sigma':
                subs_dict[t] = t + param
            elif g_type == 'pi':
                subs_dict[t] = param * t
        return KarrElement(self.expr.subs(subs_dict), self.field)

    def __add__(self, other: KarrElement | Expr) -> KarrElement:
        other_expr = other.expr if isinstance(other, KarrElement) else other
        return KarrElement(self.expr + other_expr, self.field)

    def __radd__(self, other: KarrElement | Expr) -> KarrElement:
        return self.__add__(other)

    def __sub__(self, other: KarrElement | Expr) -> KarrElement:
        other_expr = other.expr if isinstance(other, KarrElement) else other
        return KarrElement(self.expr - other_expr, self.field)

    def __rsub__(self, other: KarrElement | Expr) -> KarrElement:
        other_expr = other.expr if isinstance(other, KarrElement) else other
        return KarrElement(other_expr - self.expr, self.field)

    def __mul__(self, other: KarrElement | Expr) -> KarrElement:
        other_expr = other.expr if isinstance(other, KarrElement) else other
        return KarrElement(self.expr * other_expr, self.field)

    def __rmul__(self, other: KarrElement | Expr) -> KarrElement:
        return self.__mul__(other)

    def __truediv__(self, other: KarrElement | Expr) -> KarrElement:
        other_expr = other.expr if isinstance(other, KarrElement) else other
        return KarrElement(self.expr / other_expr, self.field)

    def __rtruediv__(self, other: KarrElement | Expr) -> KarrElement:
        other_expr = other.expr if isinstance(other, KarrElement) else other
        return KarrElement(other_expr / self.expr, self.field)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, KarrElement):
            return cancel(self.expr - other.expr) == S.Zero
        return cancel(self.expr - other) == S.Zero

    def __repr__(self) -> str:
        return f"KarrElement({self.expr})"


class KarrField:
    """Represents a difference field tower Q(t_1, ..., t_n)."""
    def __init__(self):
        self.generators = {}  # Symbol -> (type, param)
        self.ordered_generators = []  # List of Symbols (order of tower)

    def add_base(self, t: Symbol):
        self.generators[t] = ('base', S.One)
        self.ordered_generators.append(t)

    def add_sigma(self, t: Symbol, beta: Expr):
        self.generators[t] = ('sigma', cancel(beta))
        self.ordered_generators.append(t)

    def add_pi(self, t: Symbol, alpha: Expr):
        self.generators[t] = ('pi', cancel(alpha))
        self.ordered_generators.append(t)

    def element(self, expr: Expr) -> KarrElement:
        return KarrElement(expr, self)
