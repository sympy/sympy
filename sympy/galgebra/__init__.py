"""
A module for geometric algebra

"""
from ga import MV, Com, DD, Format, Nga, ReciprocalFrame, ScalarFunction, \
    cross, dual, ga_print_off, ga_print_on, inv, proj, refl, rot, rotor

from debug import oprint, ostr

from precedence import define_precedence, GAeval

from printing import enhance_print, Get_Program, Print_Function, latex, xdvi

from manifold import Manifold
