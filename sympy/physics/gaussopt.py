from sympy.physics.optics.gaussopt import BeamParameter, CurvedMirror, \
    CurvedRefraction, FlatMirror, FlatRefraction, FreeSpace, GeometricRay, \
    RayTransferMatrix, ThinLens, conjugate_gauss_beams, gaussian_conj, \
    geometric_conj_ab, geometric_conj_af, geometric_conj_bf, rayleigh2waist, \
    waist2rayleigh
from sympy.utilities.exceptions import SymPyDeprecationWarning

SymPyDeprecationWarning(feature="Module sympy.physics.gaussopt",
        useinstead="sympy.physics.optics.gaussopt",
        deprecated_since_version="0.7.6", issue=7659).warn()
