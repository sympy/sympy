from sympy import Basic, S, Expr, Symbol, Tuple
from sympy.core.sets import FiniteSet

class Domain(Basic):
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
    def __new__(cls, fulldomain, condition):
        condition = condition.subs({rs:rs.symbol
            for rs in random_symbols(condition)})
        return Domain.__new__(cls, fulldomain.symbols, fulldomain, condition)

    @property
    def fulldomain(self):
        return self.args[1]
    @property
    def condition(self):
        return self.args[2]

class PSpace(Basic):
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

    def __new__(cls, *spaces):
        from sympy.statistics.frv import ProductFinitePSpace
        from sympy.statistics.crv import ProductContinuousPSpace
        rs_space_dict = {}
        for space in spaces:
            for value in space.values:
                rs_space_dict[value] = space
        symbols = frozenset((val.symbol for val in rs_space_dict.keys()))

        if all(space.is_finite for space in spaces):
            cls = ProductFinitePSpace
        if all(space.is_continuous for space in spaces):
            cls = ProductContinuousPSpace

        obj = Basic.__new__(cls, symbols, spaces)
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


def is_random(x):
    return isinstance(x, RandomSymbol)
def is_random_expr(expr):
    return any(is_random(sym) for sym in expr.free_symbols)
def random_symbols(expr):
    try:
        return [s for s in expr.free_symbols if is_random(s)]
    except:
        return []


def pspace(expr):
    rvs = random_symbols(expr)
    if not rvs:
        return None
    # If only one space present
    if all(rv.pspace == rvs[0].pspace for rv in rvs):
        return rvs[0].pspace
    # Otherwise make a product space
    return ProductPSpace(*[rv.pspace for rv in rvs])

def sumsets(sets):
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
    if not random_symbols(expr): # expr isn't random?
        return expr
    # Create new expr and recompute E
    if given is not None: # If there is a condition
        return E(Given(expr, given, **kwargs), **kwargs)
    # Otherwise case is simple, pass work off to the ProbabilitySpace
    return pspace(expr).integrate(expr, **kwargs)


def P(condition, given=None, **kwargs):
    if given is not None: # If there is a condition
        # Recompute on new conditional expr
        return P(Given(condition, given, **kwargs), **kwargs)

    # Otherwise pass work off to the ProbabilitySpace
    return pspace(condition).P(condition, **kwargs)

def Density(expr, given=None, **kwargs):
    if given is not None: # If there is a condition
        # Recompute on new conditional expr
        return Density(Given(expr, given, **kwargs), **kwargs)

    # Otherwise pass work off to the ProbabilitySpace
    return pspace(expr).compute_density(expr, **kwargs)

def Where(condition, given=None, **kwargs):
    if given is not None: # If there is a condition
        # Recompute on new conditional expr
        return Where(Given(condition, given, **kwargs), **kwargs)

    # Otherwise pass work off to the ProbabilitySpace
    return pspace(condition).where(condition, **kwargs)
