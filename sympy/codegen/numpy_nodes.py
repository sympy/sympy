from sympy.core.function import Add, ArgumentIndexError, Function
from sympy.core.singleton import S
from sympy.functions.elementary.exponential import exp, log
from sympy.utilities import default_sort_key


def _logaddexp(x1, x2, *, evaluate=True):
    return log(Add(exp(x1, evaluate=evaluate), exp(x2, evaluate=evaluate), evaluate=evaluate))


class logaddexp(Function):
    """ Logarithm of the sum of exponentiations of the inputs.

    Helper class for use with e.g. numpy.logaddexp
    See: https://docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.logaddexp.html
    """
    nargs = 2

    def __new__(cls, *args):
        return Function.__new__(cls, *sorted(args, key=default_sort_key))

    def fdiff(self, argindex=1):
        """
        Returns the first derivative of this function.
        """
        if argindex == 1:
            wrt, other = self.args
        elif argindex == 2:
            other, wrt = self.args
        else:
            raise ArgumentIndexError(self, argindex)
        return S.One/(S.One + exp(other-wrt))

    def _eval_rewrite_as_log(self, x1, x2, **kwargs):
        return _logaddexp(x1, x2)

    def _eval_simplify(self, *args, **kwargs):
        a, b = map(lambda x: x.simplify(**kwargs), self.args)
        candidate = _logaddexp(a, b)
        if candidate != _logaddexp(a, b, evaluate=False):
            return candidate
        else:
            return logaddexp(a, b)
