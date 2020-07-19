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

    >>> from sympy.map import Map
    >>> from sympy.abc import x
    >>> class F(Map):
    ...     def eval(self, x):
    ...         return 2*x
    >>> f = F()

    >>> f.composite(f)
    F@F
    >>> f.composite(f)(x, evaluate=True)
    4*x
    >>> f.composite(f.inv(), evaluate=True)
    IdentityMap

   @ operator returns evaluated composition

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
                #3. Convert repeated composition to IteratedMap
                m1_base, m1_iternum = m1.as_base_iternum()
                m2_base, m2_iternum = m2.as_base_iternum()
                if m1_base == m2_base:
                    m1 = IteratedMap(m1_base, m1_iternum + m2_iternum, evaluate=True)
                    continue
                #4. Deal with inverse maps
                if m1.inv(evaluate=True) == m2 or m1 == m2.inv(evaluate=True):
                    m1 = IdentityMap(m1.domain)
                    continue
                #5. Deal with identity maps
                # IdentityMaps cannot be blindly removed because
                # the domain and codomain of m1 and m2 can be different.
                if isinstance(m1, IdentityMap) and m1.domain == m2.codomain:
                    m1 = m2
                    continue
                if isinstance(m2, IdentityMap) and m1.domain == m2.codomain:
                    continue
                #6. No evaluation can be done with m1 and m2
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

    def _map_content(self):
        return self.func, tuple(a._map_content() for a in self.args)

class IteratedMap(Map):
    """
    A class for n-th iterated function.

    Explanation
    ===========

    The n-th iterated function is a function composed with itself n times. This
    class has no direct relation with productional power [1]. should not be
    confused with n-fold product function which is defined in function ring [2].

    Parameters
    ==========

    f : Map
        Iterated map

    n : number of iteration
        Non-integers are allowed

    Examples
    ========

    >>> from sympy import S, Symbol
    >>> from sympy.map import Map
    >>> from sympy.abc import x
    >>> class F(Map):
    ...     def eval(self, x):
    ...         return 2*x
    >>> f = F()

    >>> f.composite(f, evaluate=True) == f.iterate(2)
    True
    >>> f.iterate(2)(x, evaluate=True)
    4*x

    Cannot evaluate if n is not Integer

    >>> f.iterate(S.One/2)(x, evaluate=True)
    IteratedMap(F, 1/2)(x)
    >>> f.iterate(Symbol('n', integer=True))(x, evaluate=True)
    IteratedMap(F, n)(x)

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

            result = b._eval_iterate(n)
            if result is not None:
                return result

            if n == 0:
                return IdentityMap(b.domain)
            if n == 1:
                return b
            if n == -1:
                return b.inv(evaluate=True)
            if ask(Q.negative(n)):
                return cls(b.inv(evaluate=True), -n, evaluate=True)

            b_base, b_iternum = b.as_base_iternum()
            if b_iternum != 1 and b_iternum != -1:
                return cls(b_base, b_iternum*n, evaluate=True)

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
        b, n = self.args
        if not n.free_symbols and ask(Q.integer(n)):
            # n is non-abstract integer
            result = arg
            for i in range(n):
                result = b(result, evaluate=True)
                if isinstance(result, AppliedMap):
                    # cannot be evaluated: abort evaluation
                    return self(arg, evaluate=False)
            return result

    def as_base_iternum(self):
        return self.args

    def _map_content(self):
        return self.func, self.base._map_content(), self.iternum
