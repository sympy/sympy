from sympy import Basic, Mul, Symbol, Indexed, integrate, Sum, Dummy, Lambda
from sympy.stats.rv import PSpace, pspace, is_random, RandomSymbol, NamedArgsMixin, random_symbols
from sympy.stats.crv import ContinuousDistribution, SingleContinuousPSpace
from sympy.stats.drv import DiscreteDistribution, SingleDiscretePSpace
from sympy.stats.frv import SingleFiniteDistribution, SingleFinitePSpace
from sympy.stats.crv_types import ContinuousDistributionHandmade
from sympy.stats.drv_types import DiscreteDistributionHandmade
from sympy.stats.frv_types import FiniteDistributionHandmade

class CompoundPSpace(PSpace):

    def __new__(cls, s, comp_distribution):
        if isinstance(s, str):
            s = Symbol(s)
        if not isinstance(s, Symbol):
            raise TypeError("s should have been string or Symbol")
        parent_dist = comp_distribution.args[0]
        return _get_newpspace(s, parent_dist, _get_pdf(parent_dist))

    @property
    def value(self):
        return RandomSymbol(self.symbol, self)

    @property
    def symbol(self):
        return self.args[0]

    @property
    def distribution(self):
        return self.args[1]

    @property
    def pdf(self):
        return self.distribution.pdf(self.symbol)

class CompoundDistribution(Basic, NamedArgsMixin):
    is_Finite = None
    is_Continuous = None
    is_Discrete = None

    def __new__(cls, dist):
        if isinstance(dist, ContinuousDistribution):
            cls.is_Continuous = True
        elif isinstance(dist, DiscreteDistribution):
            cls.is_Discrete = True
        elif isinstance(dist, SingleFiniteDistribution):
            cls.is_Finite = True
        else:
            message = "Compound Distribution for %s is not implemeted yet" % str(dist)
            raise NotImplementedError(message)
        return Basic.__new__(cls, dist)

    @property
    def set(self):
        return self.args[0].set

    def pdf(self, x):
        res = _get_pdf(self.args[0])
        print(res(x))
        return res

def _get_pdf(dist):
    randoms = []
    for arg in dist.args:
        randoms.extend(random_symbols(arg))
    if len(randoms) > 1:
        raise NotImplementedError("Compound Distributions for more than one random"
            " argument is not implemeted yet.")
    rand_sym = randoms[0]
    x = Dummy('x')
    if isinstance(dist, SingleFiniteDistribution):
        _pdf = dist.pmf(x)
    else:
        _pdf = dist.pdf(x)
    rand_dens = rand_sym.pspace.distribution.pdf(rand_sym)
    rand_sym_dom = rand_sym.pspace.domain.set
    if rand_sym.pspace.is_Discrete or rand_sym.pspace.is_Finite:
        _pdf = Sum(_pdf*rand_dens, (rand_sym, rand_sym_dom._inf, rand_sym_dom._sup)).doit()
    else:
        _pdf = integrate(_pdf*rand_dens, (rand_sym, rand_sym_dom._inf, rand_sym_dom._sup))
    return Lambda(x, _pdf)

def _get_newpspace(sym, dist, pdf):
    pdf = Lambda(sym, pdf(sym))
    _set = dist.set
    if isinstance(dist, ContinuousDistribution):
        return SingleContinuousPSpace(sym, ContinuousDistributionHandmade(pdf, _set))
    elif isinstance(dist, DiscreteDistribution):
        return SingleDiscretePSpace(sym, DiscreteDistributionHandmade(pdf, _set))
    elif isinstance(dist, SingleFiniteDistribution):
        dens = dict((k, pdf(k)) for k in _set)
        return SingleFinitePSpace(sym, FiniteDistributionHandmade(dens))
    raise NotImplementedError("CompoundDistribution for %s is not implemeted" %str(dist))
