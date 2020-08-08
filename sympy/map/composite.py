from sympy.assumptions import ask, Q
from sympy.core import Expr, S, Tuple
from sympy.core.operations import AssocOp
from sympy.core.sympify import _sympify
from .map import function_set, Map, IdentityMap, AppliedMap, InverseMap
from .operator import (
    BinaryOperator, ExponentOperator, ExponentElement
)

__all__ = [
    "CompositionOperator", "composite_op", "CompositeMap",
    "IterationOperator", "IteratedMap",
]

class CompositionOperator(BinaryOperator):
    """
    A class for function composition operator.

    Explanation
    ===========

    Function composition is an operation that takes two functions $f$ and $g$
    and produces a function $h$ such that $h(x)=g(f(x))$ when codomain of $f$
    is subset of domain of $g$ [1]. The composition of functions is always
    associative [1].

    Parameters
    ==========

    domain, codomain : Set, optional
        Domain and codomain of the operator.
        Default is set of any map.

    Examples
    ========

    >>> from sympy import Map, FunctionSet, CompositionOperator, S
    >>> fs = FunctionSet(S.UniversalSet, S.UniversalSet)
    >>> op = CompositionOperator(fs, fs)
    >>> class F(Map):
    ...     name = 'f'
    ...     def eval(self, x):
    ...         return 2*x
    >>> f = F()

    >>> op(f, f)
    f@f : UniversalSet -> UniversalSet

    Using @ operator returns evaluated composition

    >>> f@(f.inv())
    id : UniversalSet -> UniversalSet
    >>> op(f, f.inv(), evaluate=True) == f@(f.inv())
    True

    """
    latex_name = '\\circ'
    pretty_name = '∘'
    str_name = '@'

    is_associative = True
    is_commutative = False

    def __new__(cls, domain=None, codomain=None, **kwargs):
        if domain is None:
            domain = function_set**2
        if codomain is None:
            codomain = function_set        
        return super().__new__(cls, domain, codomain)

    @property
    def domain(self):
        return self.args[0]

    @property
    def codomain(self):
        return self.args[1]

    def apply(self, *args, **kwargs):
        evaluate = kwargs.get("evaluate", False)

        self.check_domain(args)

        domain = args[0].domain
        if evaluate:
            args = self.process_args(args)
            if not args:
                return IdentityMap(domain=domain)
            elif len(args) == 1:
                return args[0]

            result = self.eval(*args)
            if result is not None:
                return result

        args = Tuple(*args)
        return super(AppliedMap, CompositeMap).__new__(CompositeMap, self, args)

    def __call__(self, *args, evaluate=False, **kwargs):
        return CompositeMap(self, args, evaluate=evaluate)

    @staticmethod
    def check_domain(seq):
        for i in range(len(seq)-1):
            g, f = seq[i], seq[i+1]
            if not f.codomain.is_subset(g.domain):
                raise TypeError(
            "%s's codomain %s is not subset of %s's domain %s" % (f, f.codomain, g, g.domain))

    def process_args(self, seq):
        seq = self.flatten(seq)
        seq = self.remove_identity(seq)
        seq = self.cancel(seq)
        return seq

    def remove_identity(self, seq):
        # Remove IdentityMap from seq
        result = []
        o1 = None
        for o2 in seq:
            if o1 is None:
                o1 = o2
                continue
            # IdentityMaps are not removed unless their domain or codomain
            # exactly match with adjacent operator.
            if isinstance(o1, IdentityMap) and o1.domain == o2.codomain:
                o1 = o2
                continue
            if isinstance(o2, IdentityMap) and o1.domain == o2.codomain:
                continue
            result.append(o1)
            o1 = o2
        if o1 is not None:
            result.append(o1)
        return result

    def _binary_cancel(self, a, b):
        #1. Deal with inverse map & iterated map
        a_base, a_exp = a.as_base_exp(self)
        b_base, b_exp = b.as_base_exp(self)
        exp_op = self.exponent_operator()
        if a_base == b_base:
            exp = a_exp + b_exp
            if exp == 0:
                return []
            else:
                return [exp_op(a_base, exp, evaluate=True)]

        #2. Deal with evaluated inverse
        if a.inv(evaluate=True) == b or b.inv(evaluate=True) == a:
            return []

        #3. Find if composition is defined
        comp_val = a._eval_composite(b)
        if comp_val is not None:
            return [comp_val]

        return [a, b]

    def exponent_operator(self):
        return IterationOperator(self)

# General composition operator
composite_op = CompositionOperator(function_set, function_set)

class CompositeMap(AppliedMap, Map):
    """
    A class for the unevaluated result of function composition.

    .. note::
        ``CompositeMap`` instance is unary map. Pass tuple
        for n-dimensional argument.

    Examples
    ========

    >>> from sympy import Map
    >>> class F(Map):
    ...     def eval(self, x):
    ...         return 2*x
    >>> class G(Map):
    ...     def eval(self, x):
    ...         return x+3
    >>> f, g = F(), G()

    >>> (f@g)(1, evaluate=True)
    8

    """

    @property
    def domain(self):
        return self.arguments[-1].domain

    @property
    def codomain(self):
        return self.arguments[0].codomain

    def eval(self, arg):
        # recursively apply the argument.
        # if unevaluatable, avoid nested AppliedMaps.
        # i.e. (f@g)(x) returns just (f@g)(x), not f(g(x))
        mappings = [] # store unresolved mappings
        innermost_arg = arg
        result = arg
        for m in reversed(self.arguments):
            m_eval = m(result, evaluate=True)
            if isinstance(m_eval, AppliedMap):
                mappings.insert(0, m)
            else:
                mappings = []
                innermost_arg = m_eval
            result = m_eval

        if mappings:
            return self._new_rawargs(*mappings)(innermost_arg, evaluate=False)
        else:
            return result

    def _eval_inverse(self):
        maps = reversed([a.inverse() for a in self.arguments])
        return self._new_rawargs(*maps)

class IterationOperator(ExponentOperator):
    """
    A class for function iteration operator.

    Explanation
    ===========

    The n-th iterated function is a function composed with itself n times.

    .. note::
        This class has no direct relation with productional power [1].
        It should not be confused with n-fold product function which is
        defined in function ring [2].

    Parameters
    ==========

    function_domain : Set
        Set of the function that can be subject to iteration

    codomain : Set

    Examples
    ========

    >>> from sympy import Map, composite_op
    >>> from sympy.abc import x
    >>> class F(Map):
    ...     name = 'f'
    ...     def eval(self, x):
    ...         return 2*x
    >>> f = F()
    >>> iterate_op = composite_op.exponent_operator()

    >>> iterate_op(f, 2)
    f**2 : UniversalSet -> UniversalSet

    If n is negative, evaluating returns -n iteration of inverse map.

    >>> iterate_op(f, -2, evaluate=True)
    InverseMap(f)**2 : UniversalSet -> UniversalSet

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Iterated_function
    .. [2] https://en.wikipedia.org/wiki/Function_composition

    """

    def __call__(self, x, n, evaluate=False, **kwargs):
        return IteratedMap(self, (x, n), evaluate=evaluate)

    def apply(self, x, n, **kwargs):
        if not x.codomain.is_subset(x.domain):
            raise TypeError(
        "%s's codomain %s is not subset of %s's domain %s" % (x, x.codomain, x, x.domain))

        if kwargs.get("evaluate", False):
            n = _sympify(n)
            base, exp = x.as_base_exp(self.base_op)
            if exp != 1:
                # collect x**2**3 to x**6
                x, n = base, exp*n
                return self(x, n, **kwargs)
            return self.eval(x, n)   

    def eval(self, x, n):
        result = x._eval_iterate(n)
        if result is not None:
            return result

        if n == 0:
            return IdentityMap(x.domain)
        if n == 1:
            return x
        if n == -1:
            return x.inv(evaluate=True)
        if ask(Q.negative(n)):
            inv_x = x.inv(evaluate=True)
            return self(inv_x, -n, evaluate=True)

        if isinstance(x, IdentityMap):
            return x

class IteratedMap(ExponentElement, Map):
    """
    A class for the unevaluated result of function iteration.

    .. note::
        ``IteratedMap`` instance is unary map. Pass tuple
        for n-dimensional argument.

    Examples
    ========

    >>> from sympy import Map, composite_op, Symbol, S
    >>> from sympy.abc import x
    >>> class F(Map):
    ...     name = 'f'
    ...     def eval(self, x):
    ...         return 2*x
    >>> f = F()
    >>> iterate_op = composite_op.exponent_operator()

    >>> (iterate_op(f,2))(x, evaluate=True)
    4*x
    >>> (iterate_op(f,3))(x, evaluate=True)
    8*x

    Cannot evaluate if n is not integer

    >>> f.iterate(S.One/2)(x, evaluate=True)
    (f**(1/2))(x)
    >>> f.iterate(Symbol('n', integer=True))(x, evaluate=True)
    (f**n)(x)

    """

    def eval(self, arg):
        b, n = self.arguments
        if not n.free_symbols and ask(Q.integer(n)):
            # n is non-abstract integer
            result = arg
            for i in range(n):
                result = b(result, evaluate=True)
                if isinstance(result, AppliedMap):
                    # cannot be evaluated: abort evaluation
                    return self(arg, evaluate=False)
            return result
