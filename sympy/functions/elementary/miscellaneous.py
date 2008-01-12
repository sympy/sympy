
from sympy.core.basic import Basic, S, cache_it, cache_it_immutable
from sympy.core.function import Lambda, Function

###############################################################################
############################# SQUARE ROOT FUNCTION ############################
###############################################################################

class sqrt(Function):

    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            s = Basic.Symbol('x',dummy=True)
            return Lambda(S.Half * s**(-S.Half),s)
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        s = Basic.Symbol('x', dummy=True)
        return Lambda(s**2, s)

    def tostr(self, p=None):
        return "sqrt(%s)" % self[0]

    @classmethod
    def canonize(cls, arg):
        arg = Basic.sympify(arg)
        return arg**S.Half

    def _eval_apply_subs(self, x, old, new):
        base, exp = old.as_base_exp()
        if base==x:
            return new ** (exp/2)

    def _eval_is_zero(self):
        return isinstance(self[0], Basic.Zero)

    def as_base_exp(self):
        return self[0], S.Half

###############################################################################
############################# MINIMUM and MAXIMUM #############################
###############################################################################

class max_(Function):

    nargs = 2

    def canonize(cls, x, y):
        if isinstance(x, Basic.Number) and isinstance(y, Basic.Number):
            return max(x, y)
        if x.is_positive:
            if y.is_negative:
                return x
            if y.is_positive:
                if x.is_unbounded:
                    if y.is_unbounded:
                        return
                    return x
        elif x.is_negative:
            if y.is_negative:
                if y.is_unbounded:
                    if x.is_unbounded:
                        return
                    return x

class min_(Function):

    nargs = 2

    def canonize(cls, x, y):
        if isinstance(x, Basic.Number) and isinstance(y, Basic.Number):
            return min(x, y)
