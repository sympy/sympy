from sympy.core.function import ArgumentIndexError, Function
from sympy.functions.elementary.exponential import exp, log
from sympy.utilities import default_sort_key


def _logaddexp(x1, x2):
    return log(exp(x1) + exp(x2))


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
        return 1/(1 + exp(other-wrt))

    def _eval_rewrite_as_log(self, x1, x2, **kwargs):
        return _logaddexp(x1, x2)

    def _eval_simplify(self, *args, **kwargs):
        a, b = map(lambda x: x.simplify(**kwargs), self.args)
        candidate = _logaddexp(a, b)
        if candidate.is_Function and candidate.func == log and len(candidate.args[0].args) == 2:
            for arg in candidate.args[0].args:
                if arg.is_Function and arg.func == exp:
                    continue
                else:
                    break
            else:
                return logaddexp(a, b)  # no simplification made

        return candidate  # simplified
