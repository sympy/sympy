__all__ = []

# The following pattern is used below for importing sub-modules:
#
# 1. "from foo import *".  This imports all the names from foo.__all__ into
#    this module. But, this does not put those names into the __all__ of
#    this module. This enables "from sympy.physics.mechanics import kinematics" to
#    work.
# 2. "import foo; __all__.extend(foo.__all__)". This adds all the names in
#    foo.__all__ to the __all__ of this module. The names in __all__
#    determine which names are imported when
#    "from sympy.physics.mechanics import *" is done.

from . import kane
from .kane import *
__all__.extend(kane.__all__)

from . import rigidbody
from .rigidbody import *
__all__.extend(rigidbody.__all__)

from . import functions
from .functions import *
__all__.extend(functions.__all__)

from . import particle
from .particle import *
__all__.extend(particle.__all__)

from . import lagrange
from .lagrange import *
__all__.extend(lagrange.__all__)


#Import essential elements from physics.vector module
from sympy.physics.vector.frame import ReferenceFrame, CoordinateSym
from sympy.physics.vector.dyadic import Dyadic
from sympy.physics.vector.vector import Vector
from sympy.physics.vector.printing import (
    VectorStrPrinter as MechanicsStrPrinter,
    VectorLatexPrinter as MechanicsLatexPrinter,
    VectorPrettyPrinter as MechanicsPrettyPrinter)
from sympy.physics.vector.point import Point
from sympy.physics.vector.functions import (cross, dot, express,
                                            time_derivative, outer,
                                            kinematic_equations,
                                            get_motion_params,
                                            partial_velocity,
                                            dynamicsymbols)
from sympy.physics.vector.functions import (
     time_derivative_printing as mechanics_printing,
     vprint as mprint, vsprint as msprint,
     vpprint as mpprint, vlatex as mlatex)

#essentialnames contains all names to be imported from vector package
essentialnames = ['ReferenceFrame', 'CoordinateSym',
                  'Dyadic', 'Vector', 'MechanicsStrPrinter',
                  'MechanicsLatexPrinter',
                  'MechanicsPrettyPrinter', 'dynamicsymbols',
                  'Point', 'cross', 'dot', 'express',
                  'time_derivative', 'outer', 'kinematic_equations',
                  'get_motion_params', 'partial_velocity',
                  'mechanics_printing', 'mprint', 'msprint',
                  'mpprint', 'mlatex']

__all__.extend(essentialnames)
