""" This module contains the elliptic integrals. """

from sympy.core import S, C, pi, I
from sympy.core.function import Function, ArgumentIndexError
from sympy.functions.elementary.miscellaneous import sqrt, root
from sympy.functions.special.hyper import hyper, meijerg

###############################################################################
################### INCOMPLETE ELLIPTIC INTEGRALS #############################
###############################################################################

#
# Incomplete elliptic integral of the first kind: F(z, m)
#

class ellipticfinc(Function):
    r"""
    The incomplete elliptic integral F.

    This function is defined as:

    :math:`\mathrm{F}(z, m) = \int_0^z \frac{1}{\sqrt{1 - m \sin(t)^2}} \, \mathrm{d}t`

    Examples
    ========

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Elliptic_integrals
    .. [2] http://functions.wolfram.com/EllipticIntegrals/EllipticF
    """

    nargs = 2

#
# Incomplete elliptic integral of the second kind: E(z, m)
#

class ellipticeinc(Function):
    r"""
    The incomplete elliptic integral E.

    This function is defined as:

    :math:`\mathrm{E}(z, m) = \int_0^z \sqrt{1 - m \sin(t)^2} \, \mathrm{d}t`

    Examples
    ========

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Elliptic_integrals
    .. [2] http://functions.wolfram.com/EllipticIntegrals/EllipticE2
    """

    nargs = 2

#
# Incomplete elliptic integral of the third kind: Pi(n, z, m)
#

class ellipticpiinc(Function):
    r"""
    The incomplete elliptic integral Pi.

    This function is defined as:

    :math:`\Pi(n, z, m) = \int_0^z \frac{1}
           {\left(1 - n \sin(t)^2\right) \sqrt{1 - m \sin(t)^2} \, \mathrm{d}t`

    Examples
    ========

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Elliptic_integrals
    .. [2] http://functions.wolfram.com/EllipticIntegrals/EllipticPi3
    """

    nargs = 3

###############################################################################
##################### COMPLETE ELLIPTIC INTEGRALS #############################
###############################################################################

#
# Complete elliptic integral of the first kind: K(z)
#

class elliptick(Function):
    r"""
    The complete elliptic integral K.

    This function is defined as:

    :math:`\mathrm{K}(z) = \mathrm{F}\left(\frac{\pi}{2}, z\right)`

    Examples
    ========

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Elliptic_integrals
    .. [2] http://functions.wolfram.com/EllipticIntegrals/EllipticK
    """

    nargs = 1

#
# Complete elliptic integral of the second kind: E(z)
#

class elliptice(Function):
    r"""
    The complete elliptic integral E.

    This function is defined as:

    :math:`\mathrm{E}(z) = \mathrm{E}\left(\frac{\pi}{2}, z\right)`

    Examples
    ========

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Elliptic_integrals
    .. [2] http://functions.wolfram.com/EllipticIntegrals/EllipticE
    """

    nargs = 1

#
# Complete elliptic integral of the third kind: Pi(n, m)
#

class ellipticpi(Function):
    r"""
    The complete elliptic integral Pi.

    This function is defined as:

    :math:`\Pi(n, m) = \Pi\left(n, \frac{\pi}{2}, m\right)`

    Examples
    ========

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Elliptic_integrals
    .. [2] http://functions.wolfram.com/EllipticIntegrals/EllipticPi
    """

    nargs = 2
