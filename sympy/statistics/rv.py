from sympy import Basic, S, Expr, Symbol
from sympy.core.sets import FiniteSet

class Domain(Basic):
    is_ProductDomain = False
    def __new__(cls, symbols, *args):
        symbols = FiniteSet(*symbols)
        return Basic.__new__(cls, symbols, *args)
    @property
    def symbols(self):
        return self.args[0]
    def __contains__(self, other):
        raise NotImplementedError()
    def integrate(self, expr):
        raise NotImplementedError()

class PSpace(Basic):
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
    def computeDensity(self, expr):
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
        return Symbol('%s%d'%(cls._name, cls._count), real=True)

class RandomSymbol(Symbol):
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

class ProductPSpace(PSpace):

    def __new__(cls, *spaces):
        from sympy.statistics.frv import ProductFinitePSpace
        rs_space_dict = {}
        for space in spaces:
            for value in space.values:
                rs_space_dict[value] = space
        symbols = frozenset((val.symbol for val in rs_space_dict.keys()))

        if all(space.is_finite for space in spaces):
            cls = ProductFinitePSpace
        #if all(space.is_continuous for space in spaces):
        #    cls = ContinuousProductPSpace

        obj = Basic.__new__(cls, symbols, spaces)
        obj.rs_space_dict = rs_space_dict

        return obj

    @property
    def spaces(self):
        return self.args[1]

    @property
    def values(self):
        return sumsets(space.values for space in self.spaces)

    def integrate(self, expr, rvs=None):
        rvs = rvs or self.values
        rvs = frozenset(rvs)
        for space in self.spaces:
            expr = space.integrate(expr, rvs & space.values)
        return expr

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

        obj = Domain.__new__(cls, symbols, domains2)
        obj.sym_domain_dict = sym_domain_dict
        return obj

    @property
    def domains(self):
        return self.args[1]

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
    return ProductPSpace(*[rv.pspace for rv in rvs])

def sumsets(sets):
    return reduce(frozenset.union, sets, frozenset())

def E(expr, given=None):
    if given == None:
        space = pspace(expr)
    else:
        space = pspace(expr+given)
        space = space.conditional_space(given)
    return space.integrate(expr)

def P(condition, given=None):
    if given == None:
        space = pspace(condition)
    else:
        space = pspace(condition+given)
        space = space.conditional_space(given)
    return space.P(condition)

def Density(expr, given=None):
    if given == None:
        space = pspace(expr)
    else:
        space = pspace(expr+given)
        space = conditional_space(given)
    return space.computeDensity(expr)

def Where(condition, given=None):
    if given == None:
        space = pspace(condition)
    else:
        space = pspace(condition+given)
        space = conditional_space(given)
    return space.where(condition)
