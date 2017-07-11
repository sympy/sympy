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

from sympy.physics import vector
from sympy.physics.vector import *

from . import body, functions, kane, lagrange, linearize, particle, \
    rigidbody, system
from .body import *
from .functions import *
from .kane import *
from .lagrange import *
from .linearize import *
from .particle import *
from .rigidbody import *
from .system import *

__all__.extend(kane.__all__)

__all__.extend(rigidbody.__all__)

__all__.extend(functions.__all__)

__all__.extend(particle.__all__)

__all__.extend(lagrange.__all__)

__all__.extend(vector.__all__)

__all__.extend(linearize.__all__)

__all__.extend(body.__all__)

__all__.extend(system.__all__)
