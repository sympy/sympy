"""controlsys is a package for control systems theory. It currently only
supports linear, time invariant systems via the module 'lti'. You can
create linear control systems in state space or transfer funciton model
representation, transform theese representatin into oneanother,
interconnect systems, evaluate the systems etc.
"""
__all__ = ['lti']

# The following pattern is used below for importing sub-modules:
#
# 1. "from foo import *".  This imports all the names from foo.__all__ into
#    this module. But, this does not put those names into the __all__ of
#    this module. This enables "from sympy.physics.mechanics import kinematics" to
#    work.
# 2. "import foo; __all__.extend(foo.__all__)". This adds all the names in
#    foo.__all__ to the __all__ of this module. The names in __all__
#    determine which names are imported when
#    "from sympy.controlsys import *" is done.


from sympy.controlsys import lti
__all__.extend([lti.__all__])
