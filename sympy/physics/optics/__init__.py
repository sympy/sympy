__all__ = []

# The following pattern is used below for importing sub-modules:
#
# 1. "from foo import *".  This imports all the names from foo.__all__ into
#    this module. But, this does not put those names into the __all__ of
#    this module. This enables "from sympy.physics.optics import TWave" to
#    work.
# 2. "import foo; __all__.extend(foo.__all__)". This adds all the names in
#    foo.__all__ to the __all__ of this module. The names in __all__
#    determine which names are imported when
#    "from sympy.physics.optics import *" is done.

from . import gaussopt, medium, utils, waves
from .gaussopt import BeamParameter, CurvedMirror, CurvedRefraction, \
    FlatMirror, FlatRefraction, FreeSpace, GeometricRay, RayTransferMatrix, \
    ThinLens, conjugate_gauss_beams, gaussian_conj, geometric_conj_ab, \
    geometric_conj_af, geometric_conj_bf, rayleigh2waist, waist2rayleigh
from .medium import Medium
from .utils import deviation, hyperfocal_distance, lens_formula, \
    lens_makers_formula, mirror_formula, refraction_angle
from .waves import TWave

__all__.extend(waves.__all__)


__all__.extend(gaussopt.__all__)


__all__.extend(medium.__all__)


__all__.extend(utils.__all__)
