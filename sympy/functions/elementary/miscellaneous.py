from sympy.core.basic import S, sympify
from sympy.core.function import Function

###############################################################################
############################# SQUARE ROOT FUNCTION ############################
###############################################################################

def sqrt(arg):
    arg = sympify(arg)
    return arg**S.Half

###############################################################################
############################# MINIMUM and MAXIMUM #############################
###############################################################################

class max_(Function):

    nargs = 2

    @classmethod
    def eval(cls, x, y):
        if x.is_Number and y.is_Number:
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

    @classmethod
    def eval(cls, x, y):
        if x.is_Number and y.is_Number:
            return min(x, y)
