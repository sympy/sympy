"""
A module for geometric algebra

"""
from ga import MV, Com, DD, Format, GAeval, Nga, ReciprocalFrame, ScalarFunction, \
               cross, define_precedence, dual, ga_print_off, ga_print_on, inv, proj, \
               refl, rot, rotor, ONE, ZERO

from ga_debug import oprint,ostr

from ga_print import enhance_print , Get_Program, Print_Function, latex, xdvi

from manifold import Manifold
