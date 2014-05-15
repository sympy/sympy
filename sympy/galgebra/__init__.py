"""
A module for geometric algebra

"""
from sympy.galgebra.ga import MV, Com, DD, Format, Nga, ReciprocalFrame, ScalarFunction, \
    cross, dual, ga_print_off, ga_print_on, inv, proj, refl, rot, rotor

from sympy.galgebra.debug import oprint, ostr

from sympy.galgebra.precedence import define_precedence, GAeval

from sympy.galgebra.printing import enhance_print, Get_Program, Print_Function, latex, xdvi

from sympy.galgebra.manifold import Manifold
