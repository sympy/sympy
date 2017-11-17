"""A functions module, includes all the standard functions.

Combinatorial - factorial, fibonacci, harmonic, bernoulli...
Elementary - hyperbolic, trigonometric, exponential, floor and ceiling, sqrt...
Special - gamma, zeta,spherical harmonics...
"""

__all__ = []

from .combinatorial.factorials import (
    factorial, factorial2, rf, ff,
    binomial, RisingFactorial, FallingFactorial, subfactorial
)
__all__ += [
    "factorial", "factorial2", "rf", "ff",
    "binomial", "RisingFactorial", "FallingFactorial", "subfactorial"
]

from .combinatorial.numbers import (
    fibonacci, lucas, harmonic,
    bernoulli, bell, euler, catalan, genocchi
)
__all__ += [
    "fibonacci", "lucas", "harmonic",
    "bernoulli", "bell", "euler", "catalan", "genocchi"
]

from .elementary.miscellaneous import (
    sqrt, root, Min, Max,
    Id, real_root, cbrt
)
__all__ += [
    "sqrt", "root", "Min", "Max",
    "Id", "real_root", "cbrt"
]

from .elementary.complexes import (
    re, im, sign, Abs, conjugate, arg,
    polar_lift, periodic_argument, unbranched_argument,
    principal_branch, transpose, adjoint, polarify, unpolarify
)
__all__ += [
    "re", "im", "sign", "Abs", "conjugate", "arg",
    "polar_lift", "periodic_argument", "unbranched_argument",
    "principal_branch", "transpose", "adjoint", "polarify", "unpolarify"
]

from .elementary.trigonometric import (
    sin, cos, tan, sec, csc, cot, sinc, asin,
    acos, atan, asec, acsc, acot, atan2
)
__all__ += [
    "sin", "cos", "tan", "sec", "csc", "cot", "sinc", "asin",
    "acos", "atan", "asec", "acsc", "acot", "atan2"
]

from .elementary.exponential import exp_polar, exp, log, LambertW
__all__ += ["exp_polar", "exp", "log", "LambertW"]

from .elementary.hyperbolic import (
    sinh, cosh, tanh, coth, sech, csch, asinh,
    acosh, atanh, acoth, asech, acsch
)
__all__ += [
    "sinh", "cosh", "tanh", "coth", "sech", "csch", "asinh",
    "acosh", "atanh", "acoth", "asech", "acsch"
]

from .elementary.integers import floor, ceiling, frac
__all__ += ["floor", "ceiling", "frac"]

from .elementary.piecewise import Piecewise, piecewise_fold
__all__ += ["Piecewise", "piecewise_fold"]

from .special.error_functions import (
    erf, erfc, erfi, erf2, erfinv, erfcinv, erf2inv,
    Ei, expint, E1, li, Li, Si, Ci, Shi, Chi,
    fresnels, fresnelc
)
__all__ += [
    "erf", "erfc", "erfi", "erf2", "erfinv", "erfcinv", "erf2inv",
    "Ei", "expint", "E1", "li", "Li", "Si", "Ci", "Shi", "Chi",
    "fresnels", "fresnelc"
]

from .special.gamma_functions import (
    gamma, lowergamma, uppergamma, polygamma,
    loggamma, digamma, trigamma
)
__all__ += [
    "gamma", "lowergamma", "uppergamma", "polygamma",
    "loggamma", "digamma", "trigamma"
]

from .special.zeta_functions import (
    dirichlet_eta, zeta, lerchphi,
    polylog, stieltjes
)
__all__ += [
    "dirichlet_eta", "zeta", "lerchphi",
    "polylog", "stieltjes"
]

from .special.tensor_functions import Eijk, LeviCivita, KroneckerDelta
__all__ += ["Eijk", "LeviCivita", "KroneckerDelta"]

from .special.singularity_functions import SingularityFunction
__all__ += ["SingularityFunction"]

from .special.delta_functions import DiracDelta, Heaviside
__all__ += ["DiracDelta", "Heaviside"]

from .special.bsplines import bspline_basis, bspline_basis_set
__all__ += ["bspline_basis", "bspline_basis_set"]

from sympy.functions.special.bessel import (
    besselj, bessely, besseli, besselk, hankel1, hankel2,
    jn, yn, jn_zeros, hn1, hn2,
    airyai, airybi, airyaiprime, airybiprime
)
__all__ += [
    "besselj", "bessely", "besseli", "besselk", "hankel1", "hankel2",
    "jn", "yn", "jn_zeros", "hn1", "hn2",
    "airyai", "airybi", "airyaiprime", "airybiprime"
]

from .special.hyper import hyper, meijerg
__all__ += ["hyper", "meijerg"]

from .special.polynomials import (
    legendre, assoc_legendre, hermite, chebyshevt, chebyshevu,
    chebyshevu_root, chebyshevt_root, laguerre, assoc_laguerre,
    gegenbauer, jacobi, jacobi_normalized
)
__all__ += [
    "legendre", "assoc_legendre", "hermite", "chebyshevt", "chebyshevu",
    "chebyshevu_root", "chebyshevt_root", "laguerre", "assoc_laguerre",
    "gegenbauer", "jacobi", "jacobi_normalized"
]

from .special.spherical_harmonics import Ynm, Ynm_c, Znm
__all__ += ["Ynm", "Ynm_c", "Znm"]

from .special.elliptic_integrals import (
    elliptic_k, elliptic_f,
    elliptic_e, elliptic_pi
)
__all__ += [
    "elliptic_k", "elliptic_f",
    "elliptic_e", "elliptic_pi"
]

from .special.beta_functions import beta
__all__ += ["beta"]

from .special.mathieu_functions import (
    mathieus, mathieuc,
    mathieusprime, mathieucprime
)
__all__ += [
    "mathieus", "mathieuc",
    "mathieusprime", "mathieucprime"
]

ln = log
__all__ += ["ln"]
