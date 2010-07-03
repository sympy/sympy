"""The anti-commutator: {A,B} = A*B + B*A."""

from sympy import S, Expr, Mul, Integer
from sympy.printing.pretty.stringpict import prettyForm

from sympy.physics.quantum.qexpr import split_commutative_parts
from sympy.physics.quantum.operator import Operator
from sympy.physics.quantum.dagger import Dagger

__all__ = [
    'AntiCommutator'
]

#-----------------------------------------------------------------------------
# Anti-commutator
#-----------------------------------------------------------------------------


class AntiCommutator(Expr):
    """The standard anticommutator, in an unevaluated state.

    The commutator is defined [1] as: {A, B} = A*B + B*A, but in this class
    the anticommutator is initially unevaluated. To expand the anticommutator
    out, use the ``doit`` method.

    The arguments of the anticommutator are put into canonical order using
    ``__cmp__``, so that {B,A} becomes {A,B}.

    Parameters
    ==========
    A : Expr
        The first argument of the anticommutator {A,B}.
    B : Expr
        The second argument of the anticommutator {A,B}.

    Examples
    ========

        >>> from sympy import symbols
        >>> from sympy.physics.quantum import AntiCommutator
        >>> from sympy.physics.quantum import Operator, Dagger
        >>> x, y = symbols('x,y')
        >>> A = Operator('A')
        >>> B = Operator('B')

    Create an anticommutator and use ``doit`` to multiply them out.

        >>> ac = AntiCommutator(A,B); ac
        {A,B}
        >>> ac.doit()
        A*B + B*A

    The commutator orders it arguments in canonical order::

        >>> ac = AntiCommutator(B,A); ac
        {A,B}

    Scalar constants are factored out::

        >>> AntiCommutator(3*x*A,x*y*B)
        3*y*x**2*{A,B}

    Dagger is alto handled::

        >>> Dagger(AntiCommutator(A,B))
        {Dagger(A),Dagger(B)}

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Commutator
    """

    def __new__(cls, A, B, **old_assumptions):
        r = cls.eval(A, B)
        if r is not None:
            return r
        obj = Expr.__new__(cls, *(A, B), **{'commutative': False})
        return obj

    @classmethod
    def eval(cls, a, b):
        """The Commutator [A,B] is on canonical form if A < B.
        """
        if not (a and b): return S.Zero
        if a == b: return Integer(2)*a**2
        if a.is_commutative or b.is_commutative:
            return Integer(2)*a*b

        # [xA,yB]  ->  xy*[A,B]
        # from sympy.physics.qmul import QMul
        c_part = []
        nc_part = []
        nc_part2 = []
        if isinstance(a, Mul):
            c_part, nc_part = split_commutative_parts(a)
        if isinstance(b, Mul):
            c_part2, nc_part2 = split_commutative_parts(b)
            c_part.extend(c_part2)
        if c_part:
            a = nc_part or [a]
            b = nc_part2 or [b]
            return Mul(Mul(*c_part), cls(Mul(*a), Mul(*b)))

        # Canonical ordering of arguments
        if a.compare(b) == 1:
            return cls(b,a)

    def _eval_expand_anticommutator(self, **hints):
        # No changes, so return self
        return self

    def doit(self, **hints):
        A = self.args[0]
        B = self.args[1]
        if isinstance(A, Operator) and isinstance(B, Operator):
            try:
                comm = A._eval_anticommutator(B, **hints)
            except NotImplementedError:
                try:
                    comm = B._eval_anticommutator(A, **hints)
                except NotImplementedError:
                    comm = None
            if comm is not None:
                return comm.doit(**hints)
        return (A*B + B*A).doit(**hints)

    def _eval_dagger(self):
        return AntiCommutator(Dagger(self.args[0]), Dagger(self.args[1]))

    def _sympyrepr(self, printer, *args):
        return "%s(%s,%s)" % (self.__class__.__name__, self.args[0],\
        self.args[1])

    def _sympystr(self, printer, *args):
        return "{%s,%s}" % (self.args[0], self.args[1])

    def _pretty(self, printer, *args):
        pform = printer._print(self.args[0], *args)
        pform = prettyForm(*pform.right((prettyForm(u','))))
        pform = prettyForm(*pform.right((printer._print(self.args[1], *args))))
        pform = prettyForm(*pform.parens(left='{', right='}'))
        return pform

    def _latex(self, printer, *args):
        return "\\left{}%s,%s\\right}" % tuple([
            printer._print(arg, *args) for arg in self.args])

