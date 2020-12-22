from sympy import S, exp, log, symbols, Add, sqrt
from sympy.stats import Normal, LogNormal

# This import is not optional, it will register the main SymPy nodes with the MatchPy tree walker:
from sympy.utilities.matchpy_connector import matchpy, WildDot, WildPlus, WildStar

if matchpy:
    from matchpy import ManyToOneReplacer, ReplacementRule, Pattern, Operation

from sympy.stats.rv import RandomSymbol
from sympy.stats.crv import SingleContinuousPSpace, SingleContinuousDistribution


def _check_reparse(n):
    import re
    m = re.match(r"RandomSymbol<(.+)>", str(n))
    if m is None:
        return n
    return S(m.group(1))


def _name_combine(n1, n2, operator):
    n1 = _check_reparse(n1)
    n2 = _check_reparse(n2)
    return f"RandomSymbol<{operator(n1, n2)}>"


class _GetMatchpyReplacerForStats:

    replacer = None

    def __new__(cls):
        obj = object.__new__(cls)
        if cls.replacer is None:
            cls.register_operations()
            cls.initialize_replacer()
        return obj

    @classmethod
    def register_operations(cls):
        # MatchPy needs every node of the expression tree to be registered:
        Operation.register(RandomSymbol)
        Operation.register(SingleContinuousPSpace)
        Operation.register(SingleContinuousDistribution)

    @classmethod
    def initialize_replacer(cls):

        # I believe `WildDot` is limited to matching symbols only, this can be fixed:
        name_ = WildDot("name")
        mu_ = WildDot("mu")
        sigma_ = WildDot("sigma")

        wplus = WildPlus("wplus")
        wstar = WildStar("wstar")

        # Define the replacer:
        cls.replacer = ManyToOneReplacer()

        # Recognize `exp(Normal)` as the LogNormal distribution:
        rule1 = Pattern(exp(Normal(name_, mu_, sigma_)))
        repl1 = lambda name, mu, sigma: LogNormal(f"exp({name.name})", mu, sigma)
        cls.replacer.add(ReplacementRule(rule1, repl1))

        # Recognize `exp(Normal)` as the LogNormal distribution:
        rule2 = Pattern(log(LogNormal(name_, mu_, sigma_)))
        repl2 = lambda name, mu, sigma: Normal(f"log({name.name})", mu, sigma)
        cls.replacer.add(ReplacementRule(rule2, repl2))

        # Addition of two normal distributions:
        n1, n2, mu1, mu2, sigma1, sigma2 = symbols("n1 n2 mu1 mu2 sigma1 sigma2", cls=WildDot)
        rule3 = Pattern(wstar + Normal(n1, mu1, sigma1) + Normal(n2, mu2, sigma2))
        repl3 = lambda wstar, n1, n2, mu1, mu2, sigma1, sigma2: Add.fromiter(wstar) + Normal(_name_combine(n1, n2, Add), mu1+mu2, sqrt(sigma1**2 + sigma2**2))
        cls.replacer.add(ReplacementRule(rule3, repl3))


def probability_simplify(expr):
    replacer = _GetMatchpyReplacerForStats().replacer
    return replacer.replace(expr)
