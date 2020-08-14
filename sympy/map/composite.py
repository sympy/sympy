from sympy.assumptions import ask, Q
from sympy.core import Expr, S
from sympy.core.operations import AssocOp
from sympy.core.sympify import _sympify
from .map import Map, IdentityMap, AppliedMap

__all__ = [
    "CompositeMap", "IteratedMap",
]

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
    F@F : UniversalSet -> UniversalSet
    >>> f.composite(f)(x, evaluate=True)
    4*x
    >>> f.composite(f.inv(), evaluate=True)
    id : UniversalSet -> UniversalSet

    @ operator returns evaluated composition

    >>> f@(f.inv())
    id : UniversalSet -> UniversalSet

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Function_composition

    """
    def __new__(cls, *args, evaluate=False, **options):

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
    def flatten(cls, seq):

        newseq = []

        while seq:
            o2 = seq.pop(0)

            # 1. Denest the composite mappings
            if isinstance(o2, cls):
                args = [*o2.args]
                args.extend(seq)
                seq = args
                continue

            # newseq = [..., o1], seq = [o2, ...]
            o1 = newseq.pop() if newseq else None

            if o1 is None:
                newseq.append(o2)
                continue

            # 2. Find if composition is defined
            comp_val = o1._eval_composite(o2)
            if comp_val is not None:
                seq.insert(0, comp_val)
                continue

            # 3. Converted repeated composition to IteratedMap
            # InverseMaps are dealt here
            o1_base, o1_iternum = o1.as_base_iternum()
            o2_base, o2_iternum = o2.as_base_iternum()
            if o1_base == o2_base:
                result = IteratedMap(o1_base, o1_iternum + o2_iternum, evaluate=True)
                seq.insert(0, result)
                continue

            # 4. Deal with evaluated inverse maps
            if o1.inv(evaluate=True) == o2 or o1 == o2.inv(evaluate=True):
                seq.insert(0, IdentityMap(o1.domain))
                continue

            #5. Deal with identity maps
            # IdentityMaps cannot be blindly removed because
            # the domain and codomain of m1 and m2 can be different.
            if isinstance(o1, IdentityMap) and o1.domain == o2.codomain:
                seq.insert(0, o2)
                continue
            if isinstance(o2, IdentityMap) and o1.domain == o2.codomain:
                seq.insert(0, o1)
                continue

            #6. No evaluation can be done with m1 and m2
            newseq.append(o1)
            newseq.append(o2)

        return [], newseq, None

    @property
    def domain(self):
        return self.args[-1].domain

    @property
    def codomain(self):
        return self.args[0].codomain

    def eval(self, arg):

        #1. denest mappings

        _, self_args, _ = self.flatten([*self.args])

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

        n = _sympify(n)

        if evaluate:

            result = b._eval_iterate(n)
            if result is not None:
                return result

            if n is S.Zero:
                return IdentityMap(b.domain)
            if n is S.One:
                return b
            if n is S.NegativeOne:
                return b.inv(evaluate=True)
            if ask(Q.negative(n)):
                return cls(b.inv(evaluate=True), -n, evaluate=True)

            if isinstance(b, IdentityMap):
                return b

            b_base, b_iternum = b.as_base_iternum()
            if b_iternum is not S.One and b_iternum is not S.NegativeOne:
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
