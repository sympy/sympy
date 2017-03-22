"""
Functions with corresponding implementations in Fortran.

The functions defined in this module allows the user to express functions such as ``dsign``
as a SymPy function for symbolic manipulation.

"""
from sympy.core.function import Function
from sympy.core.numbers import Float

class FFunction(Function):
    _required_standard = 77

    def _fcode(self, printer):
        name = self.__class__.__name__
        if printer._settings['standard'] < self._required_standard:
            raise NotImplementedError("%s requires Fortran %d or newer" %
                                      (name, self._required_standard))
        return '{0}({1})'.format(name, ', '.join(map(printer._print, self.args)))

class F95Function(FFunction):
    _required_standard = 95


class isign(FFunction):
    """ Fortran sign intrinsic with for integer arguments. """
    nargs = 2


class dsign(FFunction):
    """ Fortran sign intrinsic with for double precision arguments. """
    nargs = 2


class cmplx(FFunction):
    """ Fortran complex conversion function. """
    nargs = 2  # may be extended to (2, 3) at a later point


class kind(FFunction):
    """ Fortran kind function. """
    nargs = 1


class merge(F95Function):
    """ Fortran merge function """
    nargs = 3


class _literal(Float):
    _token = None
    _decimals = None

    def _fcode(self, printer):
        mantissa, sgnd_ex = ('%.{0}e'.format(self._decimals) % self).split('e')
        mantissa = mantissa.strip('0').rstrip('.')
        ex_sgn, ex_num = sgnd_ex[0], sgnd_ex[1:].lstrip('0')
        ex_sgn = '' if ex_sgn == '+' else ex_sgn
        return (mantissa or '0') + self._token + ex_sgn + (ex_num or '0')


class literal_sp(_literal):
    """ Fortran single precision real literal """
    _token = 'e'
    _decimals = 9


class literal_dp(_literal):
    """ Fortran double precision real literal """
    _token = 'd'
    _decimals = 17
