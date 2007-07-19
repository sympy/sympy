
from basic import Basic
from function import DefinedFunction, Lambda, Function, Apply

class OrthogonalPolynomial(Function):

    def __new__(cls, index, **assumptions):
        index = Basic.sympify(index)
        if index.is_negative:
            raise ValueError("%s index must be nonnegative integer (got %r)" % (cls.__name__, index))
        obj = Basic.__new__(cls, index, **assumptions)
        obj.name = obj.__class__.__name__
        return obj

    def _hashable_content(self):
        return (self.name, self.index)

    @property
    def index(self):
        return self._args[0]

    def torepr(self):
        return '%s(%r)' % (self.__class__.__name__, self.index)

    precedence = Basic.Apply_precedence

    def tostr(self, level=0):
        r = '%s(%s)' % (self.name, self.index)
        p = self.precedence
        if p<=level: r = '(%s)' % r
        return r

    def weight(self):
        raise NotImplementedError('%s.weight()' % (self))



class Chebyshev(OrthogonalPolynomial):
    """
    Chebyshev polynomials of the first kind.
    """

    nofargs = 1
    _cache = {}
    kind = 1

    def as_lambda(self):
        i = self.index
        x = Basic.Symbol('x',dummy=True)
        cls = self.__class__
        if not isinstance(i, Basic.Integer):
            return Lambda((2*x*cls(i-1)(x)-cls(i-2)(x)),x)
        try:
            return cls._cache[i]
        except KeyError:
            pass
        if i==0:
            r = Lambda(Basic.One(),x)
        elif i==1:
            r = Lambda(self.kind * x,x)
        elif i>1:
            r = Lambda((2*x*cls(i-1).as_lambda()(x) - cls(i-2).as_lambda()(x)).expand(), x)
        else:
            raise ValueError("%s index must be nonnegative integer (got %r)" % (self, i))
        cls._cache[i] = r
        return r

    def _eval_expand(self):
        return self.as_lambda()

    def _eval_apply(self, x):
        if isinstance(x, Basic.Number):
            return self.as_lambda()(x)

        if isinstance(x, ApplyChebyshev): # nesting property
            return self.__class__(self.index * x.func.index)(x.args[0])

    def weight(self):
        x = Basic.Symbol('x',dummy=True)
        return Lambda(Basic.Sqrt()(1-x**2),x)

class Chebyshev2(Chebyshev):
    """
    Chebyshev polynomials of the second kind.
    """

    _cache = {}
    kind = 2

    def _eval_apply(self, x):
        if isinstance(x, Basic.Number):
            return self.as_lambda()(x)

class ApplyChebyshev(Apply):

    def expand(self):
        return self.func.as_lambda()(self.args[0]).expand()

class ApplyChebyshev2(ApplyChebyshev):

    pass

Basic.singleton['Chebyshev'] = lambda: Chebyshev
Basic.singleton['Chebyshev2'] = lambda: Chebyshev2
