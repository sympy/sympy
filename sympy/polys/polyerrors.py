"""Definitions of common exceptions for `polys` module. """

class BasePolynomialError(Exception):
    """Base class for polynomial related exceptions. """

    def new(self, *args):
        raise NotImplementedError("abstract base class")

class ExactQuotientFailed(BasePolynomialError):

    def __init__(self, f, g, dom=None):
        self.f, self.g, self.dom = f, g, dom

    def __str__(self): # pragma: no cover
        from sympy.printing.str import sstr

        if self.dom is None:
            return "%s does not divide %s" % (sstr(self.g), sstr(self.f))
        else:
            return "%s does not divide %s in %s" % (sstr(self.g), sstr(self.f), sstr(self.dom))

    def new(self, f, g):
        return self.__class__(f, g, self.dom)

class OperationNotSupported(BasePolynomialError):

    def __init__(self, poly, func):
        self.poly = poly
        self.func = func

    def __str__(self): # pragma: no cover
        return "`%s` operation not supported by %s representation" % (self.func, self.poly.rep.__class__.__name__)

class HeuristicGCDFailed(BasePolynomialError):
    pass

class HomomorphismFailed(BasePolynomialError):
    pass

class IsomorphismFailed(BasePolynomialError):
    pass

class ExtraneousFactors(BasePolynomialError):
    pass

class EvaluationFailed(BasePolynomialError):
    pass

class RefinementFailed(BasePolynomialError):
    pass

class CoercionFailed(BasePolynomialError):
    pass

class NotInvertible(BasePolynomialError):
    pass

class NotReversible(BasePolynomialError):
    pass

class NotAlgebraic(BasePolynomialError):
    pass

class DomainError(BasePolynomialError):
    pass

class PolynomialError(BasePolynomialError):
    pass

class UnificationFailed(BasePolynomialError):
    pass

class GeneratorsNeeded(BasePolynomialError):
    pass

class ComputationFailed(BasePolynomialError):

    def __init__(self, func, nargs, exc):
        self.func = func
        self.nargs = nargs
        self.exc = exc

    def __str__(self):
        return "%s(%s) failed without generators" % (self.func, ', '.join(map(str, self.exc.exprs[:self.nargs])))

class GeneratorsError(BasePolynomialError):
    pass

class UnivariatePolynomialError(PolynomialError):
    pass

class MultivariatePolynomialError(PolynomialError):
    pass

class PolificationFailed(PolynomialError):

    def __init__(self, opt, origs, exprs, seq=False):
        if not seq:
            self.orig = origs
            self.expr = exprs
            self.origs = [origs]
            self.exprs = [exprs]
        else:
            self.origs = origs
            self.exprs = exprs

        self.opt = opt
        self.seq = seq

    def __str__(self): # pragma: no cover
        if not self.seq:
            return "can't construct a polynomial from %s" % str(self.orig)
        else:
            return "can't construct polynomials from %s" % ', '.join(map(str, self.origs))

class OptionError(BasePolynomialError):
    pass

class FlagError(OptionError):
    pass
