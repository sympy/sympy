"""Definitions of common exceptions for `polys` module. """

class OperationNotSupported(Exception):

    def __init__(self, poly, func):
        self.poly = poly
        self.func = func

    def __str__(self): # pragma: no cover
        return "`%s` operation not supported by %s representation" % (self.func, self.poly.rep.__class__.__name__)

class ExactQuotientFailed(Exception):
    pass

class HeuristicGCDFailed(Exception):
    pass

class HomomorphismFailed(Exception):
    pass

class IsomorphismFailed(Exception):
    pass

class ExtraneousFactors(Exception):
    pass

class EvaluationFailed(Exception):
    pass

class RefinementFailed(Exception):
    pass

class CoercionFailed(Exception):
    pass

class NotInvertible(Exception):
    pass

class NotReversible(Exception):
    pass

class NotAlgebraic(Exception):
    pass

class DomainError(Exception):
    pass

class PolynomialError(Exception):
    pass

class UnificationFailed(Exception):
    pass

class GeneratorsNeeded(Exception):
    pass

class GeneratorsError(Exception):
    pass

class UnivariatePolynomialError(PolynomialError):
    pass

class MultivariatePolynomialError(PolynomialError):
    pass

class OptionError(Exception):
    pass

class PolificationFailed(PolynomialError):

    def __init__(self, origs, exprs):
        self.origs = tuple(origs)
        self.exprs = tuple(exprs)

    def __str__(self): # pragma: no cover
        return "can't construct polynomials from %s" % ', '.join(map(str, self.origs))
