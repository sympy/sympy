from sympy.core import Expr
from sympy.core.operations import AssocOp
from sympy.core.sympify import _sympify
from .map import Map, IdentityMap, AppliedMap

__all__ = ["CompositeMap", "CompositionalMapPow"]

class CompositeMap(Map, AssocOp):
    """
    A class for general composite mappings.

    .. note::
       CompositeMap is unary. Pass tuple for n-dimensional argument.

    Explanation
    ===========

    Function composition is an operation that takes two functions $f$ and $g$
    and produces a function $h$ such that $h(x)=g(f(x))$ when codomain of $f$
    is subset of domain of $g$ [1]. The composition of functions is always
    associative [1].

    Examples
    ========

    >>> from sympy.map import Map, CompositeMap
    >>> from sympy.abc import x
    >>> class F(Map):
    ...     def eval(self, x):
    ...        return 2*x
    >>> f = F()

    >>> CompositeMap(f, f)
    F@F
    >>> CompositeMap(f, f)(x).doit()
    4*x
    >>> CompositeMap(f, f.inv()).doit()
    IdentityMap

   @ operator returns evaluated composition

    >>> f@f
    F@F
    >>> f@f.inv()
    IdentityMap

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Function_composition

    """
    def __new__(cls, *args, evaluate=False, **options):

        cls.check_domain(args)

        domain = args[0].domain

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
    def check_domain(cls, seq):
        for i in range(len(seq)-1):
            g, f = seq[i], seq[i+1]
            if not f.codomain.is_subset(g.domain):
                raise TypeError(
            "%s's codomain %s is not subset of %s's domain %s" % (f, f.codomain, g, g.domain))

    @classmethod
    def flatten(cls, seq):
        #1.  denest the composite mappings
        _, seq, _ = super().flatten(seq)

        #2.  find if composition is defined
        # composition with inverse or identity are handled here
        oldseq = []
        while len(oldseq) != len(seq):
            oldseq = seq
            seq = []
            m1 = None
            for m2 in oldseq:
                if m1 is None:
                    m1 = m2
                else:
                    comp_val = m1._eval_composite(m2)
                    if comp_val is None:
                        seq.append(m1)
                        m1 = m2
                    else:
                        m1 = comp_val
            if m1 is not None:
                seq.append(m1)

        return [], seq, None

    @property
    def domain(self):
        return self.args[-1].domain

    @property
    def codomain(self):
        return self.args[0].codomain

    def eval(self, arg):

        #1. denest mappings

        self_args = self.flatten([*self.args])[1]

        #2. f(g(h(x))) -> CompositeMap(f,g,h)(x)

        mappings = [] # store unresolved mappings
        innermost_arg = arg
        result = arg
        for m in reversed(self_args):
            m_eval = m(result, evaluate=True)
            if isinstance(m_eval, AppliedMap):
                mappings.insert(0, m)
            else:
                mappings = []
                innermost_arg = m_eval
            result = m_eval

        if mappings:
            return self.func(*mappings)(innermost_arg)
        else:
            return result

    def _eval_inverse(self):
        maps = reversed([a.inverse() for a in self.args])
        return self.func(*maps)

class CompositionalMapPow(Map):
    """
    A class for functional power.

    Explanation
    ===========

    The n-th functional power is a function composed with itself n times. This
    class has no direct relation with productional power, and should not be
    confused with n-fold product function which is defined in function ring [1].

    References
    ==========
    .. [1] https://en.wikipedia.org/wiki/Function_composition

    """
    def __new__(cls, b, e, evaluate=False):

        if not b.codomain.is_subset(b.domain):
            raise TypeError(
        "%s's codomain %s is not subset of %s's domain %s" % (b, b.codomain, b, b.domain))

        e = _sympify(e)

        if evaluate:
            result = b._eval_compositionalpow(e)
            if result is not None:
                return result
        return super().__new__(cls, b, e)

    @property
    def base(self):
        return self.args[0]

    @property
    def exp(self):
        return self.args[1]
