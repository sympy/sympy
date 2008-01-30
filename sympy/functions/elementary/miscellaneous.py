
from sympy.core.basic import Basic, S, C
from sympy.core.function import Lambda, Function

###############################################################################
############################# SQUARE ROOT FUNCTION ############################
###############################################################################

def sqrt(arg):
    arg = Basic.sympify(arg)
    return arg**S.Half

###############################################################################
############################# MINIMUM and MAXIMUM #############################
###############################################################################

class max_(Function):

    nargs = 2

    def canonize(cls, x, y):
        if isinstance(x, C.Number) and isinstance(y, C.Number):
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
        if isinstance(x, C.Number) and isinstance(y, C.Number):
            return min(x, y)
