"""Definitions of common exceptions for `polys` module. """

class OperationNotSupported(Exception):

    def __init__(self, poly, func):
        self.poly = poly
        self.func = func

    def __str__(self): # pragma: no cover
        return "`%s` operation not supported by %s representation" % (self.func, self.poly.rep.__class__.__name__)

class ExactQuotientFailed(Exception):

    def __init__(self, f, g, dom=None):
        self.f, self.g, self.dom = f, g, dom

    def __str__(self): # pragma: no cover
        from sympy.printing.str import sstr

        if self.dom is None:
            return "%s does not divide %s" % (sstr(self.g), sstr(self.f))
        else:
            return "%s does not divide %s in %s" % (sstr(self.g), sstr(self.f), sstr(self.dom))

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

class ComputationFailed(Exception):

    def __init__(self, func, nargs, exc):
        self.func = func
        self.nargs = nargs
        self.exc = exc

    def __str__(self):
        return "%s(%s) failed without generators" % (self.func, ', '.join(map(str, self.exc.exprs[:self.nargs])))

class GeneratorsError(Exception):
    pass

class UnivariatePolynomialError(PolynomialError):
    pass

class MultivariatePolynomialError(PolynomialError):
    pass

class OptionError(Exception):
    pass

class FlagError(OptionError):
    pass

class PolificationFailed(PolynomialError):

    def __init__(self, origs, exprs, seq=False):
        if not seq:
            self.orig = origs
            self.expr = exprs
            self.origs = [origs]
            self.exprs = [exprs]
        else:
            self.origs = origs
            self.exprs = exprs

        self.seq = seq

    def __str__(self): # pragma: no cover
        if not self.seq:
            return "can't construct a polynomial from %s" % str(self.orig)
        else:
            return "can't construct polynomials from %s" % ', '.join(map(str, self.origs))

