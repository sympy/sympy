from sympy import Basic, S, Expr, Symbol, Tuple, And
from sympy.core.sets import FiniteSet

class Domain(Basic):
    """
    Represents a set of variables and the values which they can take
    """

    is_ProductDomain = False
    is_finite = False
    is_continuous = False
    def __new__(cls, symbols, *args):
        symbols = FiniteSet(*symbols)
        return Basic.__new__(cls, symbols, *args)

    @property
    def symbols(self):
        return self.args[0]
    @property
    def set(self):
        return self.args[1]

    def __contains__(self, other):
        raise NotImplementedError()
    def integrate(self, expr):
        raise NotImplementedError()

class SingleDomain(Domain):
    """
    A single variable and its domain
    """
    def __new__(cls, symbol, set):
        assert symbol.is_Symbol
        symbols = FiniteSet(symbol)
        return Domain.__new__(cls, symbols, set)

    @property
    def symbol(self):
        return tuple(self.symbols)[0]

    def __contains__(self, other):
        if len(other)!=1:
            return False
        sym, val = tuple(other)[0]
        return self.symbol == sym and val in self.set

class ConditionalDomain(Domain):
    """
    A Domain with an attached condition
    """
    def __new__(cls, fulldomain, condition):
        condition = condition.subs(dict((rs,rs.symbol)
            for rs in random_symbols(condition)))
        return Domain.__new__(cls, fulldomain.symbols, fulldomain, condition)

    @property
    def fulldomain(self):
        return self.args[1]
    @property
    def condition(self):
        return self.args[2]

    def as_boolean(self):
        return And(self.fulldomain.as_boolean(), self.condition)

class PSpace(Basic):
    """
    A Probability Space

    Probability Spaces encode processes that equal different values
    probabalistically. These underly Random Symbols which occur in SymPy
    expressions and contain the mechanics to evaluate statistical statements.

    """

    is_finite = None
    is_continuous = None
    @property
    def domain(self):
        return self.args[0]
    @property
    def density(self):
        return self.args[1]
    @property
    def values(self):
        return frozenset(RandomSymbol(self, sym) for sym in self.domain.symbols)
    @property
    def symbols(self):
        return self.domain.symbols
    def where(self, condition):
        raise NotImplementedError()
    def compute_density(self, expr):
        raise NotImplementedError()
    def sample(self):
        raise NotImplementedError()
    def P(self, condition):
        raise NotImplementedError()
    def integrate(self, expr):
        raise NotImplementedError()

    _count = 0
    _name = 'space'
    @classmethod
    def create_symbol(cls):
        cls._count += 1
        return Symbol('%s%d'%(cls._name, cls._count),
                real=True, finite=True, bounded=True)

class RandomSymbol(Symbol):
    """
    Random Symbols represent ProbabilitySpaces in SymPy Expressions
    In principle they can take on any value that their symbol can take on
    within the associated PSpace with probability determined by the PSpace
    Density.
    """

    is_bounded=True
    is_finite=True
    def __new__(cls, *args):
        return Basic.__new__(cls, *args)
    @property
    def pspace(self):
        return self.args[0]
    @property
    def symbol(self):
        return self.args[1]
    @property
    def name(self):
        return self.symbol.name
    @property
    def is_commutative(self):
        return self.symbol.is_commutative
    def _hashable_content(self):
        return self.args



class ProductPSpace(PSpace):
    """
    A probability space resulting from the merger of two independent probability
    spaces.

    Often created using the function, pspace
    """

    def __new__(cls, *spaces):
        from sympy.statistics.frv import ProductFinitePSpace
        from sympy.statistics.crv import ProductContinuousPSpace
        rs_space_dict = {}
        for space in spaces:
            for value in space.values:
                rs_space_dict[value] = space
        symbols = FiniteSet(val.symbol for val in rs_space_dict.keys())

        if all(space.is_finite for space in spaces):
            cls = ProductFinitePSpace
        if all(space.is_continuous for space in spaces):
            cls = ProductContinuousPSpace

        obj = Basic.__new__(cls, symbols, FiniteSet(*spaces))
        obj.rs_space_dict = rs_space_dict

        return obj

    @property
    def spaces(self):
        return self.args[1]

    @property
    def values(self):
        return sumsets(space.values for space in self.spaces)

    def integrate(self, expr, rvs=None, **kwargs):
        rvs = rvs or self.values
        rvs = frozenset(rvs)
        for space in self.spaces:
            expr = space.integrate(expr, rvs & space.values, **kwargs)
        return expr

    @property
    def domain(self):
        return ProductDomain(*[space.domain for space in self.spaces])
    @property
    def density(self):
        raise NotImplementedError("Density not available for ProductSpaces")

class ProductDomain(Domain):
    """
    A domain resulting from the merger of two independent domains
    """
    is_ProductDomain = True
    def __new__(cls, *domains):

        symbolslist = sumsets([domain.symbols for domain in domains])
        symbols = frozenset(symbolslist)
        if len(symbols) != len(symbolslist):
            raise ValueError("Overlapping Domains")

        # Flatten any product of products
        domains2 = []
        for domain in domains:
            if not domain.is_ProductDomain:
                domains2.append(domain)
            else:
                domains2.extend(domain.domains)
        domains2 = FiniteSet(domains2)

        sym_domain_dict = {}
        for domain in domains2:
            for symbol in domain.symbols:
                sym_domain_dict[symbol] = domain
        if all(domain.is_finite for domain in domains2):
            from sympy.statistics.frv import ProductFiniteDomain
            cls = ProductFiniteDomain
        if all(domain.is_continuous for domain in domains2):
            from sympy.statistics.crv import ProductContinuousDomain
            cls = ProductContinuousDomain
        obj = Domain.__new__(cls, symbols, domains2)
        obj.sym_domain_dict = sym_domain_dict
        return obj

    @property
    def domains(self):
        return self.args[1]

    @property
    def set(self):
        raise NotImplementedError("Product Sets not implemented")
        # return ProductSet(domain.set for domain in self.domains)

    def __contains__(self, other):
        # Split event into each subdomain
        for domain in self.domains:
            # Collect the parts of this event which associate to this domain
            elem = frozenset([item for item in other
                if item[0] in domain.symbols])
            # Test this sub-event
            if elem not in domain:
                return False
        # All subevents passed
        return True

    def as_boolean(self):
        return And(*[domain.as_boolean() for domain in self.domains])

def is_random(x):
    return isinstance(x, RandomSymbol)
def random_symbols(expr):
    """
    Returns all RandomSymbols within a SymPy Expression
    """
    try:
        return [s for s in expr.free_symbols if is_random(s)]
    except:
        return []


def pspace(expr):
    """
    Returns the underlying Probability Space of a random expression
    >>> from sympy.statistics import pspace, Normal
    >>> from sympy.statistics.rv import ProductPSpace
    >>> X, Y = Normal(0, 1), Normal(0, 1)
    >>> pspace(2*X + 1) == X.pspace
    True

    For internal use.
    """

    rvs = random_symbols(expr)
    if not rvs:
        return None
    # If only one space present
    if all(rv.pspace == rvs[0].pspace for rv in rvs):
        return rvs[0].pspace
    # Otherwise make a product space
    return ProductPSpace(*[rv.pspace for rv in rvs])

def sumsets(sets):
    """
    Union of sets
    """
    return reduce(frozenset.union, sets, frozenset())

def rs_swap(a,b):
    """
    Build a dictionary to swap RandomSymbols based on their underlying symbol
    i.e.
    if    X = ('x', pspace1)
    and   Y = ('x', pspace2)
    then X and Y match and the key, value pair
    {X:Y} will appear in the result
    """
    d = {}
    for rsa in a:
        d[rsa] = [rsb for rsb in b if rsa.symbol==rsb.symbol][0]
    return d

def Given(expr, given=None, **kwargs):
    """
    From a random expression and a condition on that expression creates a new
    probability space from the condition and returns the same expression on that
    conditional probability space.

    >>> from sympy.statistics import Given, Density, Die
    >>> X = Die(6)
    >>> Y = Given(X, X>3)
    >>> Density(Y)
    {4: 1/3, 5: 1/3, 6: 1/3}

    """

    if given is None:
        return expr

    # Get full probability space of both the expression and the condition
    fullspace = pspace(Tuple(expr, given))
    # Build new space given the condition
    space = fullspace.conditional_space(given, **kwargs)
    # Dictionary to swap out RandomSymbols in expr with new RandomSymbols
    # That point to the new conditional space
    swapdict = rs_swap(fullspace.values, space.values)
    # Swap
    expr = expr.subs(swapdict)
    return expr

def E(expr, given=None, **kwargs):
    """
    Returns the expected value of a random expression (optionally given a
    condition)

    >>> from sympy.statistics import E, Die
    >>> X = Die(6)
    >>> E(X)
    7/2
    >>> E(2*X + 1)
    8

    >>> E(X, X>3) # Expectation of X given that it is above 3
    5
    """

    if not random_symbols(expr): # expr isn't random?
        return expr
    # Create new expr and recompute E
    if given is not None: # If there is a condition
        return E(Given(expr, given, **kwargs), **kwargs)
    # Otherwise case is simple, pass work off to the ProbabilitySpace
    return pspace(expr).integrate(expr, **kwargs)


def P(condition, given=None, **kwargs):
    """
    Probability that a condition is true, optionally given a second condition

    >>> from sympy.statistics import P, Die
    >>> from sympy import Eq
    >>> X, Y = Die(6), Die(6)
    >>> P(X>3)
    1/2
    >>> P(Eq(X, 5), X>2) # Probability that X == 5 given that X > 2
    1/4
    >>> P(X>Y)
    5/12

    """

    if given is not None: # If there is a condition
        # Recompute on new conditional expr
        return P(Given(condition, given, **kwargs), **kwargs)

    # Otherwise pass work off to the ProbabilitySpace
    return pspace(condition).P(condition, **kwargs)

def Density(expr, given=None, **kwargs):
    """
    Probability Density of a random expression, optionally given a second
    condition

    This density will take on different forms for different types of probability
    spaces.
    Discrete RV's produce Dicts
    Continuous RV's produce a Tuple with expression representing the PDF and
    a symbol designating the active variable

    >>> from sympy.statistics import Density, Die, Normal
    >>> from sympy import Symbol

    >>> D = Die(6)
    >>> X = Normal(0, 1, symbol=Symbol('x'))

    >>> Density(D)
    {1: 1/6, 2: 1/6, 3: 1/6, 4: 1/6, 5: 1/6, 6: 1/6}
    >>> Density(2*D)
    {2: 1/6, 4: 1/6, 6: 1/6, 8: 1/6, 10: 1/6, 12: 1/6}
    >>> Density(X)
    (x, 2**(1/2)*exp(-x**2/2)/(2*pi**(1/2)))

    """
    if given is not None: # If there is a condition
        # Recompute on new conditional expr
        return Density(Given(expr, given, **kwargs), **kwargs)

    # Otherwise pass work off to the ProbabilitySpace
    return pspace(expr).compute_density(expr, **kwargs)

def Where(condition, given=None, **kwargs):
    """
    Returns the domain where a condition is True

    >>> from sympy.statistics import Where, Die, Normal
    >>> from sympy import symbols, And

    >>> x, a, b = symbols('x a b')
    >>> D1, D2 = Die(6, symbol=a), Die(6, symbol=b)
    >>> X = Normal(0, 1, symbol=x)

    >>> Where(X**2<1)
    Domain: And(-1 < x, x < 1)

    >>> Where(And(D1<=D2 , D2<3))
    Domain: Or(And(a == 1, b == 1), And(a == 1, b == 2), And(a == 2, b == 2))


    """
    if given is not None: # If there is a condition
        # Recompute on new conditional expr
        return Where(Given(condition, given, **kwargs), **kwargs)

    # Otherwise pass work off to the ProbabilitySpace
    return pspace(condition).where(condition, **kwargs)
