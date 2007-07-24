try:
    from ctypes import *
except:
    raise ImportError("ctypes is required for plotting.\n'easy_install ctypes' or visit http://sourceforge.net/projects/ctypes/")

import sys
import os
try:
    libdir = os.path.abspath(os.path.dirname(__file__))
    sys.path.insert(0, libdir)
except:
    pass

from cartesian import CartesianFunction
from polar import PolarFunction
from parametric import ParametricFunction
from spherical import SphericalFunction

from plot import Plot
