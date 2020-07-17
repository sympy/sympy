from sympy.assumptions import ask, Q
from sympy.core import Expr
from sympy.core.operations import AssocOp
from sympy.core.sympify import _sympify
from .map import Map, IdentityMap, AppliedMap

__all__ = ["CompositeMap", "IteratedMap"]

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
    >>> f@(f.inv())
    IdentityMap

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Function_composition

    """
    def __new__(cls, *args, evaluate=False, **options):

        cls.check_domain(args)

        domain = args[0].domain

        if evaluate:
            _, args, _ = cls.flatten(list(args))

            if len(args) == 0:
                # all mappings are cancelled out by inverse composition
                return IdentityMap(domain=domain)
            elif len(args) == 1:
                return args[0]

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

        oldseq = []
        while oldseq != seq:
            #1. Denest the composite mappings.
            # This must be in while loop since the user may define
            # composite of maps to return another composite of maps.
            _, oldseq, _ = super().flatten(seq)
            seq = []

            m1 = None
            for m2 in oldseq:
                if m1 is None:
                    m1 = m2
                    continue
                #2. Find if composition is defined
                comp_val = m1._eval_composite(m2)
                if comp_val is not None:
                    m1 = comp_val
                    continue
                #3. Check if m1 and m2 are inverse.
                # Unevaluated substructures are not evaluated here
                # to check if m1 and m2 are inverse for performance reason.
                # Calling doit(deep=True) will deal with such case.
                if m1 == m2.inv(evaluate=True):
                    m1 = IdentityMap(m1.codomain)
                    continue
                #4. Deal with identity maps
                # IdentityMaps cannot be blindly removed because
                # the domain and codomain of m1 and m2 can be different.
                if isinstance(m1, IdentityMap) and m1.domain == m2.codomain:
                    m1 = m2
                    continue
                if isinstance(m2, IdentityMap) and m1.domain == m2.codomain:
                    continue
                #5. No evaluation can be done with m1 and m2
                seq.append(m1)
                m1 = m2
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

class IteratedMap(Map):
    """
    A class for n-th iterated function.

    Explanation
    ===========

    The n-th iterated function is a function composed with itself n times. This
    class has no direct relation with productional power [1]. should not be
    confused with n-fold product function which is defined in function ring [2].

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Iterated_function
    .. [2] https://en.wikipedia.org/wiki/Function_composition

    """
    def __new__(cls, b, n, evaluate=False):

        if not b.codomain.is_subset(b.domain):
            raise TypeError(
        "%s's codomain %s is not subset of %s's domain %s" % (b, b.codomain, b, b.domain))

        n = _sympify(n)

        if evaluate:

            if n == 0:
                return IdentityMap(b.domain)
            if n == 1:
                return b
            if ask(Q.negative(n)):
                return cls(b.inv(evaluate=True), -n, evaluate=True)

            b_base, b_iternum = b.as_base_iternum()
            if b_iternum != 1:
                return cls(b.base, b.iternum*n, evaluate=True)

            result = b._eval_iteration(n)
            if result is not None:
                return result
        return super().__new__(cls, b, n)

    @property
    def base(self):
        return self.args[0]

    @property
    def iternum(self):
        return self.args[1]

    @property
    def domain(self):
        return self.base.domain

    @property
    def codomain(self):
        return self.base.codomain

    def eval(self, arg):
        new_self = self.doit(deep=False)

        if not isinstance(new_self, self.func):
            return new_self(arg, evaluate=True)

        b, n = new_self.args
        if not n.free_symbols and ask(Q.integer(n)):
            # n is non-abstract integer
            result = arg
            for i in range(n):
                result = b(result, evaluate=True)
                if isinstance(result, AppliedMap):
                    # cannot be evaluated: abort evaluation
                    return new_self(arg, evaluate=False)
            return result

    def as_base_iternum(self):
        return self.args
