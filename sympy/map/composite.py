from sympy.core import Expr
from sympy.core.operations import AssocOp
from .map import Map, IdentityMap

__all__ = ["CompositeMap",]

class CompositeMap(Map, AssocOp):
    """
    A class for general composite mappings.

    Explanation
    ===========

    Function composition is an operation that takes two functions $f$ and $g$
    and produces a function $h$ such that $h(x)=g(f(x))$ [1]. The composition
    of functions is always associative [1].

    Examples
    ========

    >>> from sympy.map import Map, CompositeMap
    >>> from sympy.abc import x
    >>> class F(Map):
    ...     def eval(self, x):
    ...        return 2*x
    >>> f = F()

    >>> CompositeMap(f, f)(x).doit()
    4*x
    >>> CompositeMap(f, f.inv())(x).doit()
    (x,)

    Notes
    =====

    Only unary maps can be arguments of CompositeMap. If you need multivariate
    map, make maps that take tuple withoug argument unpacking. See
    test_composite.py for example.

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Function_composition

    """
    def __new__(cls, *args, evaluate=False, **options):
        domain = args[0].domain
        codomain = args[-1].codomain

        if not evaluate:
            return Expr.__new__(cls, *args)
        else:
            _, args, _ = cls.flatten(list(args))

            if len(args) == 0:
                # all mappings are cancelled out by inverse composition
                return IdentityMap(domain=domain)
            elif len(args) == 1:
                return args[0]
            else:
                return Expr.__new__(cls, *args)

    @classmethod
    def flatten(cls, seq):
        # denest the composite mappings
        _, seq, _ = super().flatten(seq)

        # cancel out composition of inverse mappings
        oldseq = []
        while len(oldseq) != len(seq):
            oldseq = seq
            seq = []
            m1 = None
            for m2 in oldseq:
                if m1 is None:
                    m1 = m2
                else:
                    are_inverse = m1.doit() == m2.inv().doit()
                    if are_inverse:
                        m1 = None
                    else:
                        seq.append(m1)
                        m1 = m2
            if m1 is not None:
                seq.append(m1)

        return [], seq, None

    def eval(self, arg):
        result = arg
        for a in reversed(self.args):
            result = a(result, evaluate=True)
        return result

    def _eval_inverse(self):
        maps = reversed([a.inverse() for a in self.args])
        return self.func(*maps)
