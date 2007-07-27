"""
This module contains functions and classes for arbitrary-precision
numerical computations.

>>> from sympy.modules.numerics import *
>>> Float.setdps(50)
>>> print pi_float()
3.1415926535897932384626433832795028841971693993751
>>> print exp(1)
2.7182818284590452353602874713526624977572470937000
>>> print log(2)
0.69314718055994530941723212145817656807550013436025
>>> print sqrt(2)
1.4142135623730950488016887242096980785696718753769

"""

from float_ import *
from constants import *
from functions import *
from functions2 import *
from quad import TanhSinh
