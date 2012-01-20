"""
Main Random Variables Module

Defines abstract random variable type.
Contains interfaces for probability space object (PSpace) as well as standard
operators, P, E, Sample, Density, Where

See Also
========
sympy.stats.crv
sympy.stats.frv
sympy.stats.rv_interface
"""

from sympy import Basic, S, Expr, Symbol, Tuple, And, Add, Eq
from sympy.core.sets import FiniteSet, ProductSet

class RandomDomain(Basic):
    """
    Represents a set of variables and the values which they can take

    See Also
    ========
    sympy.stats.crv.ContinuousDomain
    sympy.stats.frv.FiniteDomain
    """

    is_ProductDomain = False
    is_Finite = False
    is_Continuous = False

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

class SingleDomain(RandomDomain):
    """
    A single variable and its domain

    See Also
    ========
    sympy.stats.crv.SingleContinuousDomain
    sympy.stats.frv.SingleFiniteDomain
    """
    def __new__(cls, symbol, set):
        assert symbol.is_Symbol
        symbols = FiniteSet(symbol)
        return RandomDomain.__new__(cls, symbols, set)

    @property
    def symbol(self):
        return tuple(self.symbols)[0]

    def __contains__(self, other):
        if len(other)!=1:
            return False
        sym, val = tuple(other)[0]
        return self.symbol == sym and val in self.set

class ConditionalDomain(RandomDomain):
    """
    A RandomDomain with an attached condition

    See Also
    ========
    sympy.stats.crv.ConditionalContinuousDomain
    sympy.stats.frv.ConditionalFiniteDomain
    """
    def __new__(cls, fulldomain, condition):
        condition = condition.subs(dict((rs,rs.symbol)
            for rs in random_symbols(condition)))
        return RandomDomain.__new__(
                cls, fulldomain.symbols, fulldomain, condition)

    @property
    def fulldomain(self):
        return self.args[1]

    @property
    def condition(self):
        return self.args[2]

    @property
    def set(self):
        raise NotImplementedError("Set of Conditional Domain not Implemented")

    def as_boolean(self):
        return And(self.fulldomain.as_boolean(), self.condition)

class PSpace(Basic):
    """
    A Probability Space

    Probability Spaces encode processes that equal different values
    probabalistically. These underly Random Symbols which occur in SymPy
    expressions and contain the mechanics to evaluate statistical statements.

    See Also
    ========
    sympy.stats.crv.ContinuousPSpace
    sympy.stats.frv.FinitePSpace
    """

    is_Finite = None
    is_Continuous = None

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

class SinglePSpace(PSpace):
    """
    Represents the probabilities of a set of random events that can be
    attributed to a single variable/symbol.
    """

    @property
    def value(self):
        return tuple(self.values)[0]

class RandomSymbol(Symbol):
    """
    Random Symbols represent ProbabilitySpaces in SymPy Expressions
    In principle they can take on any value that their symbol can take on
    within the associated PSpace with probability determined by the PSpace
    Density.

    Random Symbols contain pspace and symbol properties.
    The pspace property points to the represented Probability Space
    The symbol is a standard SymPy Symbol that is used in that probability space
    for example in defining a density.

    You can form normal SymPy expressions using RandomSymbols and operate on
    those expressions with the Functions

    E - Expectation of a random expression
    P - Probability of a condition
    Density - Probability Density of an expression
    Given - A new random expression (with new random symbols) given a condition

    An object of the RandomSymbol type should almost never be created by the
    user. They tend to be created instead by the PSpace class's value method.
    Traditionally a user doesn't even do this but instead calls one of the
    convenience functions Normal, Exponential, Coin, Die, FiniteRV, etc....
    """

    is_bounded=True
    is_finite=True

    def __new__(cls, *args):
        obj = Basic.__new__(cls)
        obj.pspace = args[0]
        obj.symbol = args[1]
        return obj

    @property
    def name(self):
        return self.symbol.name

    @property
    def is_commutative(self):
        return self.symbol.is_commutative

    def _hashable_content(self):
        return self.pspace, self.symbol

class ProductPSpace(PSpace):
    """
    A probability space resulting from the merger of two independent probability
    spaces.

    Often created using the function, pspace
    """

    def __new__(cls, *spaces):
        rs_space_dict = {}
        for space in spaces:
            for value in space.values:
                rs_space_dict[value] = space

        symbols = FiniteSet(val.symbol for val in rs_space_dict.keys())

        # Overlapping symbols
        if len(symbols) < sum(len(space.symbols) for space in spaces):
            raise ValueError("Overlapping Random Variables")

        if all(space.is_Finite for space in spaces):
            from sympy.stats.frv import ProductFinitePSpace
            cls = ProductFinitePSpace
        if all(space.is_Continuous for space in spaces):
            from sympy.stats.crv import ProductContinuousPSpace
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

    def sample(self):
        return dict([(k,v) for space in self.spaces
            for k,v in space.sample().items()])

class ProductDomain(RandomDomain):
    """
    A domain resulting from the merger of two independent domains

    See Also
    ========
    sympy.stats.crv.ProductContinuousDomain
    sympy.stats.frv.ProductFiniteDomain
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

        if all(domain.is_Finite for domain in domains2):
            from sympy.stats.frv import ProductFiniteDomain
            cls = ProductFiniteDomain
        if all(domain.is_Continuous for domain in domains2):
            from sympy.stats.crv import ProductContinuousDomain
            cls = ProductContinuousDomain

        obj = RandomDomain.__new__(cls, symbols, domains2)
        obj.sym_domain_dict = sym_domain_dict
        return obj

    @property
    def domains(self):
        return self.args[1]

    @property
    def set(self):
        return ProductSet(domain.set for domain in self.domains)

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
    Returns all RandomSymbols within a SymPy Expression.
    """
    try:
        return [s for s in expr.free_symbols if is_random(s)]
    except:
        return []

def pspace(expr):
    """
    Returns the underlying Probability Space of a random expression.

    For internal use.

    Examples
    ========

    >>> from sympy.stats import pspace, Normal
    >>> from sympy.stats.rv import ProductPSpace
    >>> X, Y = Normal(0, 1), Normal(0, 1)
    >>> pspace(2*X + 1) == X.pspace
    True
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
    Build a dictionary to swap RandomSymbols based on their underlying symbol.

    i.e.
    if    X = ('x', pspace1)
    and   Y = ('x', pspace2)
    then X and Y match and the key, value pair
    {X:Y} will appear in the result

    Inputs: collections a and b of random variables which share common symbols
    Output: dict mapping RVs in a to RVs in b
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

    Examples
    ========

    >>> from sympy.stats import Given, Density, Die
    >>> X = Die(6)
    >>> Y = Given(X, X>3)
    >>> Density(Y)
    {4: 1/3, 5: 1/3, 6: 1/3}
    """

    if not random_symbols(given) or pspace_independent(expr, given):
        return expr

    # Get full probability space of both the expression and the condition
    fullspace = pspace(Tuple(expr, given))
    # Build new space given the condition
    space = fullspace.conditional_space(given, **kwargs)
    # Dictionary to swap out RandomSymbols in expr with new RandomSymbols
    # That point to the new conditional space
    swapdict = rs_swap(fullspace.values, space.values)
    # Swap random variables in the expression
    expr = expr.subs(swapdict)
    return expr

def E(expr, given=None, numsamples=None, **kwargs):
    """
    Returns the expected value of a random expression

    (optionally given a condition)

    Examples
    ========

    >>> from sympy.stats import E, Die
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
    if numsamples: # Computing by monte carlo sampling?
        return sampling_E(expr, given, numsamples=numsamples, **kwargs)

    # Create new expr and recompute E
    if given is not None: # If there is a condition
        return E(Given(expr, given, **kwargs), **kwargs)

    # A few known statements for efficiency
    if expr.is_Add:
        return Add(*[E(arg) for arg in expr.args]) # E is Linear

    # Otherwise case is simple, pass work off to the ProbabilitySpace
    return pspace(expr).integrate(expr, **kwargs)


def P(condition, given=None, numsamples=None,  **kwargs):
    """
    Probability that a condition is true, optionally given a second condition

    Examples
    ========

    >>> from sympy.stats import P, Die
    >>> from sympy import Eq
    >>> X, Y = Die(6), Die(6)
    >>> P(X>3)
    1/2
    >>> P(Eq(X, 5), X>2) # Probability that X == 5 given that X > 2
    1/4
    >>> P(X>Y)
    5/12
    """

    if numsamples:
        return sampling_P(condition, given, numsamples=numsamples, **kwargs)
    if given is not None: # If there is a condition
        # Recompute on new conditional expr
        return P(Given(condition, given, **kwargs), **kwargs)

    # Otherwise pass work off to the ProbabilitySpace
    return pspace(condition).P(condition, **kwargs)

def Density(expr, given=None, **kwargs):
    """
    Probability Density of a random expression

    Optionally given a second condition

    This density will take on different forms for different types of probability
    spaces.
    Discrete RV's produce Dicts
    Continuous RV's produce a Tuple with expression representing the PDF and
    a symbol designating the active variable

    Examples
    ========

    >>> from sympy.stats import Density, Die, Normal
    >>> from sympy import Symbol

    >>> D = Die(6)
    >>> X = Normal(0, 1, symbol=Symbol('x'))

    >>> Density(D)
    {1: 1/6, 2: 1/6, 3: 1/6, 4: 1/6, 5: 1/6, 6: 1/6}
    >>> Density(2*D)
    {2: 1/6, 4: 1/6, 6: 1/6, 8: 1/6, 10: 1/6, 12: 1/6}
    >>> Density(X)
    (x, sqrt(2)*exp(-x**2/2)/(2*sqrt(pi)))
    """
    if given is not None: # If there is a condition
        # Recompute on new conditional expr
        return Density(Given(expr, given, **kwargs), **kwargs)

    # Otherwise pass work off to the ProbabilitySpace
    return pspace(expr).compute_density(expr, **kwargs)

def CDF(expr, given=None, **kwargs):
    """
    Cumulative Distribution Function of a random expression.

    optionally given a second condition

    This density will take on different forms for different types of probability
    spaces.
    Discrete RV's produce list of tuples.
    Continuous RV's produce a Tuple with expression representing the PDF and
    a symbol designating the active variable.

    Examples
    ========

    >>> from sympy.stats import Density, Die, Normal, CDF
    >>> from sympy import Symbol

    >>> D = Die(6)
    >>> X = Normal(0, 1)

    >>> Density(D)
    {1: 1/6, 2: 1/6, 3: 1/6, 4: 1/6, 5: 1/6, 6: 1/6}
    >>> CDF(D)
    {1: 1/6, 2: 1/3, 3: 1/2, 4: 2/3, 5: 5/6, 6: 1}
    >>> CDF(3*D, D>2)
    {9: 1/4, 12: 1/2, 15: 3/4, 18: 1}

    >>> CDF(X)
    (_z, erf(sqrt(2)*_z/2)/2 + 1/2)
    """
    if given is not None: # If there is a condition
        # Recompute on new conditional expr
        return CDF(Given(expr, given, **kwargs), **kwargs)

    # Otherwise pass work off to the ProbabilitySpace
    return pspace(expr).compute_cdf(expr, **kwargs)

def Where(condition, given=None, **kwargs):
    """
    Returns the domain where a condition is True.

    Examples
    ========

    >>> from sympy.stats import Where, Die, Normal
    >>> from sympy import symbols, And

    >>> x, a, b = symbols('x a b')
    >>> D1, D2 = Die(6, symbol=a), Die(6, symbol=b)
    >>> X = Normal(0, 1, symbol=x)

    >>> Where(X**2<1)
    Domain: And(-1 < x, x < 1)

    >>> Where(X**2<1).set
    (-1, 1)

    >>> Where(And(D1<=D2 , D2<3))
    Domain: Or(And(a == 1, b == 1), And(a == 1, b == 2), And(a == 2, b == 2))
    """
    if given is not None: # If there is a condition
        # Recompute on new conditional expr
        return Where(Given(condition, given, **kwargs), **kwargs)

    # Otherwise pass work off to the ProbabilitySpace
    return pspace(condition).where(condition, **kwargs)

def Sample(expr, given=None, **kwargs):
    """
    A realization of the random expression

    Examples
    ========

    >>> from sympy.stats import Die, Sample
    >>> X, Y, Z = Die(6), Die(6), Die(6)

    >>> die_roll = Sample(X+Y+Z) # A random realization of three dice
    """
    return sample_iter(expr, given, numsamples=1).next()

def sample_iter(expr, given=None, numsamples=S.Infinity, **kwargs):
    """
    Returns an iterator of realizations from the expression given a condition

    expr: Random expression to be realized
    given: A conditional expression (optional)
    numsamples: Length of the iterator (defaults to infinity)

    See Also
    ========
    Sample
    sampling_P
    sampling_E
    """

    if given:
        ps = pspace(Tuple(expr, given))
    else:
        ps = pspace(expr)

    count = 0
    while count < numsamples:
        d = ps.sample() # a dictionary that maps RVs to values

        if given: # Check that these values satisfy the condition
            gd = given.subs(d)
            if not isinstance(gd, bool):
                raise ValueError("Conditions must not contain free symbols")
            if gd == False: # If the values don't satisfy then try again
                continue

        ed = expr.subs(d)

        yield ed
        count += 1

def sampling_P(condition, given=None, numsamples=1, **kwargs):
    """
    Sampling version of P

    See Also
    ========
    P
    sampling_E
    """

    count_true = 0
    count_false = 0

    samples = sample_iter(condition, given, numsamples=numsamples, **kwargs)

    for x in samples:
        if not isinstance(x, bool):
            raise ValueError("Conditions must not contain free symbols")

        if x==True:
            count_true += 1
        else:
            count_false += 1

    return S(count_true) / numsamples

def sampling_E(condition, given=None, numsamples=1, **kwargs):
    """
    Sampling version of E

    See Also
    ========
    P
    sampling_P
    """

    samples = sample_iter(condition, given, numsamples=numsamples, **kwargs)

    return Add(*samples) / numsamples

def dependent(a, b):
    """
    Dependence of two random expressions

    Two expressions are independent if knowledge of one does not change
    computations on the other.

    Examples
    ========

    >>> from sympy.stats import Normal, dependent, Given
    >>> from sympy import Tuple, Eq

    >>> X, Y = Normal(0, 1), Normal(0, 1)
    >>> dependent(X, Y)
    False
    >>> dependent(2*X + Y, -Y)
    True
    >>> X, Y = Given(Tuple(X, Y), Eq(X+Y,3))
    >>> dependent(X, Y)
    True

    See Also
    ========
    independent
    """
    if pspace_independent(a,b):
        return False

    z = Symbol('z', real=True)
    # Dependent if density is unchanged when one is given information about
    # the other
    return (Density(a, Eq(b, z)) != Density(a) or
            Density(b, Eq(a, z)) != Density(b))

def independent(a, b):
    """
    Independence of two random expressions

    Two expressions are independent if knowledge of one does not change
    computations on the other.

    Examples
    ========

    >>> from sympy.stats import Normal, independent, Given
    >>> from sympy import Tuple, Eq

    >>> X, Y = Normal(0, 1), Normal(0, 1)
    >>> independent(X, Y)
    True
    >>> independent(2*X + Y, -Y)
    False
    >>> X, Y = Given(Tuple(X, Y), Eq(X+Y,3))
    >>> independent(X, Y)
    False

    See Also
    ========
    dependent
    """
    return not dependent(a, b)

def pspace_independent(a,b):
    """
    Tests for independence between a and b by checking if their PSpaces have
    overlapping symbols. This is a sufficient but not necessary condition for
    independence and is intended to be used internally.
    Note:
    pspace_independent(a,b) implies independent(a,b)
    independent(a,b) does not imply pspace_independent(a,b)
    """
    a_symbols = pspace(b).symbols
    b_symbols = pspace(a).symbols
    if len(a_symbols.intersect(b_symbols)) == 0:
        return True
    return None

def rv_subs(expr, symbols=None):
    """
    Given a random expression replace all random variables with their symbols.

    If symbols keyword is given restrict the swap to only the symbols listed.
    """
    if symbols is None:
        symbols = random_symbols(expr)
    swapdict = dict([(rv, rv.symbol) for rv in symbols])
    return expr.subs(swapdict)

