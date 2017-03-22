"""
Functions with corresponding implementations in Fortran.

The functions defined in this module allows the user to express functions such as ``dsign``
as a SymPy function for symbolic manipulation.

"""
from sympy.core.function import Function


class isign(Function):
    """ Fortran sign intrinsic with for integer arguments. """
    nargs = 2

class dsign(Function):
    """ Fortran sign intrinsic with for double precision arguments. """
    nargs = 2

class cmplx(Function):
    """ Fortran complex conversion function. """
    nargs = 2  # may be extended to (2, 3) at a later point

class kind(Function):
    """ Fortran kind function. """
    nargs = 1
