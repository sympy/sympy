from sympy import Basic, integrate, Sum, Dummy, Lambda
from sympy.stats.rv import NamedArgsMixin, random_symbols, _symbol_converter
from sympy.stats.crv import ContinuousDistribution, SingleContinuousPSpace
from sympy.stats.drv import DiscreteDistribution, SingleDiscretePSpace
from sympy.stats.frv import SingleFiniteDistribution, SingleFinitePSpace
from sympy.stats.crv_types import ContinuousDistributionHandmade
from sympy.stats.drv_types import DiscreteDistributionHandmade
from sympy.stats.frv_types import FiniteDistributionHandmade


def compound_pspace(s, comp_distribution):
    """
    A temporary Probability Space function for the Compound Distribution. After
    Marginalization, this returns the corresponding Probability Space of the
    parent distribution.
    """
    s = _symbol_converter(s)
    if not isinstance(comp_distribution, CompoundDistribution):
        raise ValueError("%s should be an isinstance of "
                        "CompoundDistribution"%(comp_distribution))
    x = Dummy('x')
    parent_dist = comp_distribution.args[0]
    new_pspace = _get_newpspace(s, parent_dist,
                        Lambda(x, comp_distribution.pdf(x)))
    if new_pspace is not None:
        return new_pspace
    message = ("Compound Distribution for %s is not implemeted yet" % str(parent_dist))
    raise NotImplementedError(message)


def _get_newpspace(sym, dist, pdf):
    """
    This function returns the new pspace of the distribution using handmade
    Distributions and their corresponding pspace.
    """
    pdf = Lambda(sym, pdf(sym))
    _set = dist.set
    if isinstance(dist, ContinuousDistribution):
        return SingleContinuousPSpace(sym, ContinuousDistributionHandmade(pdf, _set))
    elif isinstance(dist, DiscreteDistribution):
        return SingleDiscretePSpace(sym, DiscreteDistributionHandmade(pdf, _set))
    elif isinstance(dist, SingleFiniteDistribution):
        dens = dict((k, pdf(k)) for k in _set)
        return SingleFinitePSpace(sym, FiniteDistributionHandmade(dens))


class CompoundDistribution(Basic, NamedArgsMixin):
    """
    Class for Compound Distributions.

    Parameters
    ==========

    dist : Distribution
        Distribution must contain a random parameter

    Examples
    ========

    >>> from sympy.stats.compound_rv import CompoundDistribution
    >>> from sympy.stats.crv_types import NormalDistribution
    >>> from sympy.stats import Normal
    >>> from sympy.abc import x
    >>> X = Normal('X', 2, 4)
    >>> N = NormalDistribution(X, 4)
    >>> C = CompoundDistribution(N)
    >>> C.set
    Interval(-oo, oo)
    >>> C.pdf(x).simplify()
    exp(-x**2/64 + x/16 - 1/16)/(8*sqrt(pi))

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Compound_probability_distribution

    """

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
        if not cls._compound_check(dist):
            raise ValueError("There are no random arguments for considering it as "
                            "Compound Distribution.")
        return Basic.__new__(cls, dist)

    @property
    def set(self):
        return self.args[0].set

    def pdf(self, x):
        dist = self.args[0]
        randoms = []
        for arg in dist.args:
            randoms.extend(random_symbols(arg))
        if len(randoms) > 1:
            raise NotImplementedError("Compound Distributions for more than"
                " one random argument is not implemeted yet.")
        rand_sym = randoms[0]
        if isinstance(dist, SingleFiniteDistribution):
            y = Dummy('y', integer=True, negative=False)
            _pdf = dist.pmf(y)
        else:
            y = Dummy('y')
            _pdf = dist.pdf(y)
        if isinstance(rand_sym.pspace.distribution, SingleFiniteDistribution):
            rand_dens = rand_sym.pspace.distribution.pmf(rand_sym)
        else:
            rand_dens = rand_sym.pspace.distribution.pdf(rand_sym)
        rand_sym_dom = rand_sym.pspace.domain.set
        if rand_sym.pspace.is_Discrete or rand_sym.pspace.is_Finite:
            _pdf = Sum(_pdf*rand_dens, (rand_sym, rand_sym_dom._inf,
                    rand_sym_dom._sup)).doit()
        else:
            _pdf = integrate(_pdf*rand_dens, (rand_sym, rand_sym_dom._inf,
                    rand_sym_dom._sup))
        return Lambda(y, _pdf)(x)

    @classmethod
    def _compound_check(self, dist):
        """
        Checks if the given distribution contains random parameters.
        """
        randoms = []
        for arg in dist.args:
            randoms.extend(random_symbols(arg))
        if len(randoms) == 0:
            return False
        return True
