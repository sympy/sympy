"""The commutator: [A,B] = A*B - B*A."""

from sympy import S, Expr, Mul, Add
from sympy.printing.pretty.stringpict import prettyForm

from sympy.physics.quantum.qexpr import split_commutative_parts
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.operator import Operator


__all__ = [
    'Commutator'
]

#-----------------------------------------------------------------------------
# Commutator
#-----------------------------------------------------------------------------



class Commutator(Expr):
    """The standard commutator, in an unevaluated state.

    The commutator is defined [1] as: [A, B] = A*B - B*A, but in this class
    the commutator is initially unevaluated. To expand the commutator out,
    use the ``doit`` method.

    The arguments of the commutator are put into canonical order using
    ``__cmp__``, so that [B,A] becomes -[A,B].

    Parameters
    ==========
    A : Expr
        The first argument of the commutator [A,B].
    B : Expr
        The second argument of the commutator [A,B].

    Examples
    ========

        >>> from sympy import symbols
        >>> from sympy.physics.quantum import Commutator, Dagger
        >>> x, y = symbols('x,y')
        >>> A, B, C = symbols('A,B,C', commutative=False)

    Create some commutators and use ``doit`` to multiply them out.

        >>> comm = Commutator(A,B); comm
        [A,B]
        >>> comm.doit()
        A*B - B*A

    The commutator orders it arguments in canonical order::

        >>> comm = Commutator(B,A); comm
        -[A,B]

    Scalar constants are factored out::

        >>> Commutator(3*x*A,x*y*B)
        3*y*x**2*[A,B]

    Using ``expand(commutator=True)``, the standard commutator expansion rules
    can be applied::

        >>> Commutator(A+B,C).expand(commutator=True)
        [A,C] + [B,C]
        >>> Commutator(A,B+C).expand(commutator=True)
        [A,B] + [A,C]
        >>> Commutator(A*B,C).expand(commutator=True)
        A*[B,C] + [A,C]*B
        >>> Commutator(A,B*C).expand(commutator=True)
        B*[A,C] + [A,B]*C

    Commutator works with Dagger::

        >>> Dagger(Commutator(A,B))
        -[Dagger(A),Dagger(B)]

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
        if a == b: return S.Zero
        if a.is_commutative or b.is_commutative:
            return S.Zero

        # [xA,yB]  ->  xy*[A,B]
        # from sympy.physics.qmul import QMul
        c_part = c_part2 = []
        nc_part = nc_part2 = []
        if isinstance(a, Mul):
            c_part, nc_part = split_commutative_parts(a)
        if isinstance(b, Mul):
            c_part2, nc_part2 = split_commutative_parts(b)
            c_part.extend(c_part2)
        if c_part:
            a = nc_part or [a]
            b = nc_part2 or [b]
            return Mul(*c_part)*cls(Mul(*a),Mul(*b))

        # Canonical ordering of arguments
        if a.compare(b) == 1:
            return S.NegativeOne*cls(b,a)

    def _eval_expand_commutator(self, **hints):
        A = self.args[0].expand(**hints)
        B = self.args[1].expand(**hints)

        result = None

        if isinstance(A, Add):
            # [A+B,C]  ->  [A,C] + [B,C]
            result = Add(
                *[Commutator(term,B).expand(**hints)\
                  for term in A.args]
            )
        elif isinstance(B, Add):
            # [A,B+C]  ->  [A,B] + [A,C]
            result = Add(
                *[Commutator(A,term).expand(**hints)\
                  for term in B.args]
            )
        elif isinstance(A, Mul):
            # [A*B,C] -> A*[B,C] + [A,C]*B
            a = A.args[0]
            b = Mul(*A.args[1:])
            c = B
            comm1 = Commutator(b,c).expand(**hints)
            comm2 = Commutator(a,c).expand(**hints)
            first = Mul(a, comm1)
            second = Mul(comm2, b)
            result = Add(first, second)
        elif isinstance(B, Mul):
            # [A,B*C] -> [A,B]*C + B*[A,C]
            a = A
            b = B.args[0]
            c = Mul(*B.args[1:])
            comm1 = Commutator(a,b).expand(**hints)
            comm2 = Commutator(a,c).expand(**hints)
            first = Mul(comm1, c)
            second = Mul(b, comm2)
            result = Add(first, second)

        if result is None:
            # No changes, so return self
            return self
        else:
            return result

    def doit(self, **hints):
        A = self.args[0]
        B = self.args[1]
        if isinstance(A, Operator) and isinstance(B, Operator):
            try:
                comm = A._eval_commutator(B, **hints)
            except NotImplementedError:
                try:
                    comm = -1*B._eval_commutator(A, **hints)
                except NotImplementedError:
                    comm = None
            if comm is not None:
                return comm.doit(**hints)
        return (A*B - B*A).doit(**hints)

    def _eval_dagger(self):
        return Commutator(Dagger(self.args[1]), Dagger(self.args[0]))

    def _sympyrepr(self, printer, *args):
        return "%s(%s,%s)" % (self.__class__.__name__, self.args[0],\
        self.args[1])

    def _sympystr(self, printer, *args):
        return "[%s,%s]" % (self.args[0], self.args[1])

    def _pretty(self, printer, *args):
        pform = printer._print(self.args[0], *args)
        pform = prettyForm(*pform.right((prettyForm(u','))))
        pform = prettyForm(*pform.right((printer._print(self.args[1], *args))))
        pform = prettyForm(*pform.parens(left='[', right=']'))
        return pform

    def _latex(self, printer, *args):
        return "\\left[%s,%s\\right]" % tuple([
            printer._print(arg, *args) for arg in self.args])

