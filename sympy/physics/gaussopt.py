# -*- encoding: utf-8 -*-
"""
Gaussian optics.

The module implements:
    Ray transfer matrices for geometrical and gaussian optics
     See RayTransferMatrix, GeometricRay and BeamParameter
    Conjugation relations for geometrical and gaussian optics
     See geometric_conj*, gauss_conj and conjugate_gauss_beams

The conventions for the distances are as follows:
    focal distance - positive for convergent lenses
    object distance - positive for real objects
    image distance - positive for real images

Module uses dual formalism to represent ray transfer matrix and geometric
ray. One of formalisms uses unimodular matrices (det(M) = 1) in which
geometric ray is represented by two paramaeters: distance to optical axis
and optical-direction cosine (which equals n*u, n - refractive index,
u - angle to optical axis). The second formalism uses distance and angle
to optical axis. The first formalism is set by default.

See Also
========

[1] Gerrard, Anthony; Burch, James M. (1994).
    Introduction to matrix methods in optics. Courier Dover.
"""

from sympy import (atan2, Expr, I, im, Matrix, oo, pi, re, sqrt, sympify,
    together)
from sympy.utilities.misc import filldedent


def isunimodular(function):
    """Review kwargs argument in function and deside
    wheather unimodular or non-unimodular formalism needed
    """
    def _isunimodular(*args, **kwargs):
        if len(kwargs) > 1:
            raise ValueError('''To many named aguments only one needed''')
        elif len(kwargs) == 1:
            if 'unimodular' in kwargs:
                if kwargs['unimodular'] == True:
                    pass
                elif kwargs['unimodular'] == False:
                    pass
                else:
                    raise ValueError('''Unimodular can be only True of False''')
            else:
                raise ValueError('''Only unimodular can be named agument''')
        elif len(kwargs) == 0:
            kwargs['unimodular'] = True
        res = function(*args, **kwargs)
        return res
    return _isunimodular

###
# A, B, C, D matrices
###

class RayTransferMatrix(Matrix):
    """
    Base class for a Ray Transfer Matrix.

    It should be used if there isn't already a more specific subclass mentioned
    in See Also.

    Parameters
    ==========

    parameters A, B, C and D or 2x2 matrix (Matrix(2, 2, [A, B, C, D]))

    Examples
    =======

    >>> from sympy.physics.gaussopt import RayTransferMatrix, ThinLens
    >>> from sympy import Symbol, Matrix

    >>> mat = RayTransferMatrix(1, 2, 3, 4)
    >>> mat
    [1, 2]
    [3, 4]

    >>> RayTransferMatrix(Matrix([[1, 2], [3, 4]]), unimodular=False)
    [1, 2]
    [3, 4]

    >>> mat.A
    1

    >>> f = Symbol('f')
    >>> lens = ThinLens(f)
    >>> lens
    [   1, 0]
    [-1/f, 1]

    >>> lens.C
    -1/f

    See Also
    ========

    GeometricRay, BeamParameter,
    FreeSpace, FlatRefraction, CurvedRefraction,
    FlatMirror, CurvedMirror, ThinLens

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Ray_transfer_matrix_analysis
    """

    _unimodular = True

    def __new__(cls, *args, **kwargs):
        if len(args) == 4:
            temp = ((args[0], args[1]), (args[2], args[3]))
        elif len(args) == 1 \
             and isinstance(args[0], Matrix) \
             and args[0].shape == (2, 2):
            temp = args[0]
        else:
            raise ValueError(filldedent('''
                Expecting 2x2 Matrix or the 4 elements of
                the Matrix but got %s''' % str(args)))
        return super(RayTransferMatrix, cls).__new__(cls, temp)

    @isunimodular
    def __init__(self, *args, **kwargs):
        self._unimodular = kwargs['unimodular']

    def __mul__(self, other):
        if isinstance(other, RayTransferMatrix):
            if self.unimodular != other.unimodular:
                raise ValueError('''Only unimodular or non-unimodular formalism
                can be applied''')
            return RayTransferMatrix(Matrix.__mul__(self, other),
                                     unimodular=self.unimodular)
        elif isinstance(other, GeometricRay):
            if self.unimodular != other.unimodular:
                raise ValueError('''Only unimodular or non-unimodular formalism
                can be applied''')
            return GeometricRay(Matrix.__mul__(self, other),
                                unimodular=self.unimodular)
        elif isinstance(other, BeamParameter):
            temp = self*Matrix(((other.q,), (1,)))
            q = (temp[0]/temp[1]).expand(complex=True)
            return BeamParameter(other.wavelen, \
                                 together(re(q)), \
                                 z_r = together(im(q)))
        else:
            return Matrix.__mul__(self, other)

    @property
    def unimodular(self):
        return self._unimodular

    @property
    def A(self):
        """
        The A parameter of the Matrix.

        Examples
        ========

        >>> from sympy.physics.gaussopt import RayTransferMatrix
        >>> mat = RayTransferMatrix(1, 2, 3, 4)
        >>> mat.A
        1
        """
        return self[0, 0]

    @property
    def B(self):
        """
        The B parameter of the Matrix.

        Examples
        ========

        >>> from sympy.physics.gaussopt import RayTransferMatrix
        >>> mat = RayTransferMatrix(1, 2, 3, 4)
        >>> mat.B
        2
        """
        return self[0, 1]

    @property
    def C(self):
        """
        The C parameter of the Matrix.

        Examples
        ========

        >>> from sympy.physics.gaussopt import RayTransferMatrix
        >>> mat = RayTransferMatrix(1, 2, 3, 4)
        >>> mat.C
        3
        """
        return self[1, 0]

    @property
    def D(self):
        """
        The D parameter of the Matrix.

        Examples
        ========

        >>> from sympy.physics.gaussopt import RayTransferMatrix
        >>> mat = RayTransferMatrix(1, 2, 3, 4)
        >>> mat.D
        4
        """
        return self[1, 1]

class FreeSpace(RayTransferMatrix):
    """
    Ray Transfer Matrix for free space.

    Parameters
    ==========

    distance
    n index or refraction (oly for unimodular case)
    and
    unimodular (default = True)

    See Also
    ========

    RayTransferMatrix

    Examples
    ========

    >>> from sympy.physics.gaussopt import FreeSpace
    >>> from sympy import symbols
    >>> d,n = symbols('d n')
    >>> FreeSpace(d)
    [1, d]
    [0, 1]

    >>> FreeSpace(d, n)
    [1, d/n]
    [0,   1]

    >>> FreeSpace(d, n).unimodular
    True

    >>> FreeSpace(d, unimodular=False).unimodular
    False
    """

    @isunimodular
    def __new__(cls, d, n=1, **kwargs):
        d = sympify((d))
        if kwargs['unimodular'] == False:
            return RayTransferMatrix.__new__(cls, 1, d, 0, 1,
                                             unimodular=kwargs['unimodular'])
        n = sympify((n))
        return RayTransferMatrix.__new__(cls, 1, d/n, 0, 1,
                                         unimodular=kwargs['unimodular'])

class FlatRefraction(RayTransferMatrix):
    """
    Ray Transfer Matrix for refraction.

    Parameters
    ==========

    n1: refractive index of one medium
    n2: refractive index of other medium
    and
    unimodular (default = True)

    See Also
    ========

    RayTransferMatrix

    Examples
    ========

    >>> from sympy.physics.gaussopt import FlatRefraction
    >>> from sympy import symbols
    >>> n1, n2 = symbols('n1 n2')
    >>> FlatRefraction(n1, n2, unimodular=False)
    [1,     0]
    [0, n1/n2]

    >>> FlatRefraction(n1, n2, unimodular=False).unimodular
    False
    """
    @isunimodular
    def __new__(cls, n1, n2, **kwargs):
        if kwargs['unimodular'] == False:
            n1, n2 = sympify((n1, n2))
            return RayTransferMatrix.__new__(cls, 1, 0, 0, n1/n2,
                                             unimodular=kwargs['unimodular'])
        return RayTransferMatrix.__new__(cls, 1, 0, 0, 1,
                                         unimodular=kwargs['unimodular'])

class CurvedRefraction(RayTransferMatrix):
    """
    Ray Transfer Matrix for refraction on curved interface.

    Parameters
    ==========

    R: radius of curvature (positive for concave),
    n1: refractive index of one medium
    n2: refractive index of other medium
    or
    P: optical power ((n2 - n1) / R = 1/f, where f - focal length)
    and
    unimodular (default = True)

    See Also
    ========

    RayTransferMatrix

    Examples
    ========

    >>> from sympy.physics.gaussopt import CurvedRefraction
    >>> from sympy import symbols
    >>> R, n1, n2 = symbols('R n1 n2')
    >>> CurvedRefraction(R, n1, n2)
    [          1, 0]
    [(n1 - n2)/R, 1]

    >>> CurvedRefraction(R, n1, n2, unimodular=False)
    [               1,     0]
    [(n1 - n2)/(R*n2), n1/n2]
    """

    @isunimodular
    def __new__(cls, R, n1, n2, **kwargs):
        if kwargs['unimodular'] == False:
            return RayTransferMatrix.__new__(cls, 1, 0, (n1-n2)/R/n2, n1/n2,
                                             unimodular=kwargs['unimodular'])
        return RayTransferMatrix.__new__(cls, 1, 0, (n1-n2)/R, 1,
                                         unimodular=kwargs['unimodular'])

class FlatMirror(RayTransferMatrix):
    """
    Ray Transfer Matrix for reflection.

    Parameters
    ==========
    unimodular (default = True)

    See Also
    ========

    RayTransferMatrix

    Examples
    ========

    >>> from sympy.physics.gaussopt import FlatMirror
    >>> FlatMirror()
    [1, 0]
    [0, 1]
    """

    @isunimodular
    def __new__(cls, **kwargs):
        return RayTransferMatrix.__new__(cls, 1, 0, 0, 1,
                                         unimodular=kwargs['unimodular'])

class CurvedMirror(RayTransferMatrix):
    """
    Ray Transfer Matrix for reflection from curved surface.

    Parameters
    ==========

    radius of curvature (positive for concave)
    and
    unimodular (default = True)

    See Also
    ========

    RayTransferMatrix

    Examples
    ========

    >>> from sympy.physics.gaussopt import CurvedMirror
    >>> from sympy import symbols
    >>> R = symbols('R')
    >>> CurvedMirror(R)
    [  1, 0]
    [2/R, 1]
    """

    @isunimodular
    def __new__(cls, R, **kwargs):
        R = sympify(R)
        if kwargs['unimodular'] == False:
            return RayTransferMatrix.__new__(cls, 1, 0, -2/R, 1,
                                             unimodular=kwargs['unimodular'])
        return RayTransferMatrix.__new__(cls, 1, 0, 2/R, 1,
                                         unimodular=kwargs['unimodular'])

class ThinLens(RayTransferMatrix):
    """
    Ray Transfer Matrix for a thin lens.

    Parameters
    ==========

    the focal distance

    See Also
    ========

    RayTransferMatrix

    Examples
    ========

    >>> from sympy.physics.gaussopt import ThinLens
    >>> from sympy import symbols
    >>> f = symbols('f')
    >>> ThinLens(f)
    [   1, 0]
    [-1/f, 1]
    """

    @isunimodular
    def __new__(cls, f, **kwargs):
        f = sympify(f)
        return RayTransferMatrix.__new__(cls, 1, 0, -1/f, 1,
                                         unimodular=kwargs['unimodular'])

class ThinPrism(Matrix):
    """Ray Transfer matrix for thin prism

    Adds delta = -n1 * angle * (n2 - n1) to optical direction cosine
    [1, 0] + [                      0]
    [0, 1]   [-n1 * angle * (n2 - n1)]

    Parameters
    ==========

    angle: top angle of the prism
    n1: refractive index of one medium (default=1.0)
    n2: refractive index of prism medium

    See Also
    ========

    Matrix

    [1] Matrix Methods for Optical Layouts, Gerhaed Kloos,
        Bellingham, Washington USA, ISBN 978-0-8194-6780-5

    Examples
    ========
    >>> from sympy import symbols
    >>> from sympy.physics.gaussopt import ThinPrism
    >>> alpha, n2, n1 = symbols('alpha n2 n1')
    >>> ThinPrism(alpha, n2, n1)
    [                   0]
    [-alpha*n1*(-n1 + n2)]

    >>> from sympy.physics.gaussopt import GeometricRay
    >>> h, angle, n = symbols('h angle n')
    >>> ThinPrism(alpha, n2)*GeometricRay(h, angle*n)
    [                        h]
    [alpha*(-n2 + 1) + angle*n]
    >>> ThinPrism(alpha, n2, n1)*GeometricRay(h, angle*n)
    [                             h]
    [-alpha*n1*(-n1 + n2) + angle*n]
    """

    _unimodular = True

    @isunimodular
    def __new__(cls, angle, n2, n1=1, **kwargs):
        angle, n2, n1 = sympify((angle, n2, n1))
        if kwargs['unimodular'] == False:
            return Matrix.__new__(cls, [[0,], [(n1 - n2) * angle,]])
        return Matrix.__new__(cls, [[0,], [-n1 * (n2 - n1) * angle,]])

    @isunimodular
    def __init__(self, *args, **kwargs):
        self._unimodular = kwargs['unimodular']

    def __mul__(self, other):
        if isinstance(other, Matrix) and other.shape == (2, 1):
            if isinstance(other, GeometricRay):
                return GeometricRay(Matrix.__add__(self,other),
                                    unimodular=self.unimodular)
            else:
                return Matrix.__add__(self,other)
        else:
            raise ValueError("Only 2x1 Matrix can be multiplied with ThinPrism")

    @property
    def unimodular(self):
        return self._unimodular

# TODO Add "Duct" (radially variad index and gain)
# TODO Add class OpticalSystem which will contain optical elements

###
# Representation for geometric ray
###

class GeometricRay(Matrix):
    """
    Representation for a geometric ray in the Ray Transfer Matrix formalism.

    Parameters
    ==========

    n : refraction index (default = 1.0), needed only for unimodular formalism
    height: distance to optical axis
    angle: optical direction angle; equals n*u in case of unimodular formalism,
    where u - angle to optical axis)
    or
    2x1 matrix (Matrix(2, 1, [height, angle]))

    Examples
    =======

    >>> from sympy.physics.gaussopt import GeometricRay, FreeSpace
    >>> from sympy import symbols, Matrix
    >>> d, h, angle, n = symbols('d h angle n')

    >>> GeometricRay(h, angle*n)
    [      h]
    [angle*n]

    >>> FreeSpace(d)*GeometricRay(h, angle*n)
    [angle*d*n + h]
    [      angle*n]

    >>> GeometricRay( Matrix( ((h,), (angle*n,)) ))
    [      h]
    [angle*n]

    See Also
    ========

    RayTransferMatrix

    """

    _unimodular = True

    @isunimodular
    def __new__(cls, *args, **kwargs):
        if len(args) == 1 and isinstance(args[0], Matrix) \
                          and args[0].shape == (2, 1):
            temp = args[0]
        elif len(args) == 2:
            temp = ((args[0],), (args[1],))
        else:
            raise ValueError(filldedent('''
                Expecting 2x1 Matrix or the 2 elements of
                the Matrix but got %s''' % str(args)))
        return Matrix.__new__(cls, temp)

    @isunimodular
    def __init__(self, *args, **kwargs):
        self._unimodular = kwargs['unimodular']

    @property
    def height(self):
        """
        The distance from the optical axis.

        Examples
        ========

        >>> from sympy.physics.gaussopt import GeometricRay
        >>> from sympy import symbols
        >>> h, angle = symbols('h, angle')
        >>> gRay = GeometricRay(h, angle, unimodular=False)
        >>> gRay.height
        h
        """
        return self[0]

    @property
    def angle(self):
        """
        The angle with the optical axis.

        Examples
        ========

        >>> from sympy.physics.gaussopt import GeometricRay
        >>> from sympy import symbols
        >>> h, angle = symbols('h, angle')
        >>> gRay = GeometricRay(h, angle, unimodular=False)
        >>> gRay.angle
        angle
        """
        return self[1]

    @property
    def unimodular(self):
        return self._unimodular


###
# Representation for gauss beam
###

class BeamParameter(Expr):
    """
    Representation for a gaussian ray in the Ray Transfer Matrix formalism.

    Parameters
    ==========

    wavelength, distance to waist, and w (waist) or z_r (rayleigh range)

    Examples
    ========

    >>> from sympy.physics.gaussopt import BeamParameter
    >>> p = BeamParameter(530e-9, 1, w=1e-3)
    >>> p.q
    1 + 1.88679245283019*I*pi

    >>> p.q.n()
    1.0 + 5.92753330865999*I
    >>> p.w_0.n()
    0.00100000000000000
    >>> p.z_r.n()
    5.92753330865999

    >>> from sympy.physics.gaussopt import FreeSpace
    >>> fs = FreeSpace(10)
    >>> p1 = fs*p
    >>> p.w.n()
    0.00101413072159615
    >>> p1.w.n()
    0.00210803120913829

    >>> from sympy.physics.gaussopt import CurvedMirror
    >>> from sympy import symbols
    >>> R, t, n = 1.5, 2.0, 1.76
    >>> Resonator = FreeSpace(t,n)*CurvedMirror(R)*\
    FreeSpace(t,n)*CurvedMirror(-R)
    >>> outp = (Resonator**10)*BeamParameter(530e-9, 1, w=1e-3)
    >>> outp.w.n()
    0.00260674329962508

    See Also
    ========

    RayTransferMatrix

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Complex_beam_parameter
    """
    #TODO A class Complex may be implemented. The BeamParameter may
    # subclass it. See:
    # https://groups.google.com/d/topic/sympy/7XkU07NRBEs/discussion

    __slots__ = ['z', 'z_r', 'wavelen']

    def __new__(cls, wavelen, z, **kwargs):
        wavelen, z = sympify((wavelen, z))
        inst = Expr.__new__(cls, wavelen, z, **kwargs)
        inst.wavelen = wavelen
        inst.z = z
        if len(kwargs) !=1:
            raise ValueError('Constructor expects exactly one named argument.')
        elif 'z_r' in kwargs:
            inst.z_r = sympify(kwargs['z_r'])
        elif 'w' in kwargs:
            inst.z_r = waist2rayleigh(sympify(kwargs['w']), wavelen)
        else:
            raise ValueError('The constructor needs named argument w or z_r')
        return inst

    @property
    def q(self):
        """
        The complex parameter representing the beam.

        Examples
        ========

        >>> from sympy.physics.gaussopt import BeamParameter
        >>> p = BeamParameter(530e-9, 1, w=1e-3)
        >>> p.q
        1 + 1.88679245283019*I*pi
        """
        return self.z + I*self.z_r

    @property
    def radius(self):
        """
        The radius of curvature of the phase front.

        Examples
        ========

        >>> from sympy.physics.gaussopt import BeamParameter
        >>> p = BeamParameter(530e-9, 1, w=1e-3)
        >>> p.radius
        0.2809/pi**2 + 1
        """
        return self.z*(1+(self.z/self.z_r)**2)

    @property
    def w(self):
        """
        The beam radius at 1/e^2 intensity.

        See Also
        ========

        w_0: minimal radius of beam

        Examples
        ========

        >>> from sympy.physics.gaussopt import BeamParameter
        >>> p = BeamParameter(530e-9, 1, w=1e-3)
        >>> p.w
        0.001*sqrt(0.2809/pi**2 + 1)
        """
        return self.w_0*sqrt(1+(self.z/self.z_r)**2)

    @property
    def w_0(self):
        """
        The beam waist (minimal radius).


        See Also
        ========

        w: beam radius at 1/e^2 intensity

        Examples
        ========

        >>> from sympy.physics.gaussopt import BeamParameter
        >>> p = BeamParameter(530e-9, 1, w=1e-3)
        >>> p.w_0
        0.00100000000000000
        """
        return sqrt(self.z_r/pi*self.wavelen)

    @property
    def divergence(self):
        """
        Half of the total angular spread.

        Examples
        ========

        >>> from sympy.physics.gaussopt import BeamParameter
        >>> p = BeamParameter(530e-9, 1, w=1e-3)
        >>> p.divergence
        0.00053/pi
        """
        return self.wavelen/pi/self.w_0

    @property
    def gouy(self):
        """
        The Gouy phase.

        Examples
        ========

        >>> from sympy.physics.gaussopt import BeamParameter
        >>> p = BeamParameter(530e-9, 1, w=1e-3)
        >>> p.gouy
        atan(0.53/pi)
        """
        return atan2(self.z, self.z_r)

    @property
    def waist_approximation_limit(self):
        """
        The minimal waist for which the gauss beam approximation is valid.

        The gauss beam is a solution to the paraxial equation. For curvatures
        that are too great it is not a valid approximation.

        Examples
        ========

        >>> from sympy.physics.gaussopt import BeamParameter
        >>> p = BeamParameter(530e-9, 1, w=1e-3)
        >>> p.waist_approximation_limit
        1.06e-6/pi
        """
        return 2*self.wavelen/pi


###
# Utilities
###

def waist2rayleigh(w, wavelen):
    """
    Calculate the rayleigh range from the waist of a gaussian beam.

    See Also
    ========

    rayleigh2waist, BeamParameter

    Examples
    ========

    >>> from sympy.physics.gaussopt import waist2rayleigh
    >>> from sympy import symbols
    >>> w, wavelen = symbols('w wavelen')
    >>> waist2rayleigh(w, wavelen)
    pi*w**2/wavelen
    """
    w, wavelen = sympify((w, wavelen))
    return w**2*pi/wavelen

def rayleigh2waist(z_r, wavelen):
    """Calculate the waist from the rayleigh range of a gaussian beam.

    See Also
    ========

    waist2rayleigh, BeamParameter

    Examples
    ========

    >>> from sympy.physics.gaussopt import rayleigh2waist
    >>> from sympy import symbols
    >>> z_r, wavelen = symbols('z_r wavelen')
    >>> rayleigh2waist(z_r, wavelen)
    sqrt(wavelen*z_r)/sqrt(pi)
    """
    z_r, wavelen = sympify((z_r, wavelen))
    return sqrt(z_r/pi*wavelen)


def geometric_conj_ab(a, b):
    """
    Conjugation relation for geometrical beams under paraxial conditions.

    Takes the distances to the optical element and returns the needed
    focal distance.

    See Also
    ========

    geometric_conj_af, geometric_conj_bf

    Examples
    ========

    >>> from sympy.physics.gaussopt import geometric_conj_ab
    >>> from sympy import symbols
    >>> a, b = symbols('a b')
    >>> geometric_conj_ab(a, b)
    a*b/(a + b)
    """
    a, b = sympify((a, b))
    if abs(a) == oo or abs(b) == oo:
        return a if abs(b) == oo else b
    else:
        return a*b/(a+b)

def geometric_conj_af(a, f):
    """
    Conjugation relation for geometrical beams under paraxial conditions.

    Takes the object distance (for geometric_conj_af) or the image distance
    (for geometric_conj_bf) to the optical element and the focal distance.
    Then it returns the other distance needed for conjugation.

    See Also
    ========

    geometric_conj_ab

    Examples
    ========

    >>> from sympy.physics.gaussopt import geometric_conj_af, geometric_conj_bf
    >>> from sympy import symbols
    >>> a, b, f = symbols('a b f')
    >>> geometric_conj_af(a, f)
    a*f/(a - f)
    >>> geometric_conj_bf(b, f)
    b*f/(b - f)
    """
    a, f = sympify((a, f))
    return -geometric_conj_ab(a, -f)

geometric_conj_bf = geometric_conj_af

def gaussian_conj(s_in, z_r_in, f):
    """
    Conjugation relation for gaussian beams.

    Parameters
    ==========

    s_in: distance to optical element from the waist
    z_r_in: the rayleigh range of the incident beam
    f: the focal length of the optical element

    Returns
    =======

    A tuple containing (s_out, z_r_out, m)
     - s_out - distance between the new waist and the optical element
     - z_r_out - rayleigh range of the emergent beam
     - m - the ration between the new and the old waists

    Examples
    ========

    >>> from sympy.physics.gaussopt import gaussian_conj
    >>> from sympy import symbols
    >>> s_in, z_r_in, f = symbols('s_in z_r_in f')

    >>> gaussian_conj(s_in, z_r_in, f)[0]
    1/(-1/(s_in + z_r_in**2/(-f + s_in)) + 1/f)

    >>> gaussian_conj(s_in, z_r_in, f)[1]
    z_r_in/(1 - s_in**2/f**2 + z_r_in**2/f**2)

    >>> gaussian_conj(s_in, z_r_in, f)[2]
    1/sqrt(1 - s_in**2/f**2 + z_r_in**2/f**2)
    """
    s_in, z_r_in, f = sympify((s_in, z_r_in, f))
    s_out = 1 / ( -1/(s_in + z_r_in**2/(s_in-f)) + 1/f )
    m = 1/sqrt((1-(s_in/f)**2) + (z_r_in/f)**2)
    z_r_out = z_r_in / ((1-(s_in/f)**2) + (z_r_in/f)**2)
    return (s_out, z_r_out, m)

def conjugate_gauss_beams(wavelen, waist_in, waist_out, **kwargs):
    """
    Find the optical setup conjugating the object/image waists.

    Parameters
    ==========

    wavelen: the wavelength of the beam
    waist_in and waist_out: the waists to be conjugated
    f: the focal distance of the element used in the conjugation

    Returns
    =======

    A tuple containing (s_in, s_out, f)
     - s_in - distance before the optical element
     - s_out - distance after the optical element
     - f -  focal distance of the optical element

    Examples
    ========

    >>> from sympy.physics.gaussopt import conjugate_gauss_beams
    >>> from sympy import symbols, factor
    >>> l, w_i, w_o, f = symbols('l w_i w_o f')

    >>> conjugate_gauss_beams(l, w_i, w_o, f=f)[0]
    f*(-sqrt(w_i**2/w_o**2 - pi**2*w_i**4/(f**2*l**2)) + 1)

    >>> factor(conjugate_gauss_beams(l, w_i, w_o, f=f)[1])
    f*w_o**2*(w_i**2/w_o**2 - sqrt(w_i**2/w_o**2 - pi**2*w_i**4/(f**2*l**2)))/w_i**2

    >>> conjugate_gauss_beams(l, w_i, w_o, f=f)[2]
    f
    """
    #TODO add the other possible arguments
    wavelen, waist_in, waist_out = sympify((wavelen, waist_in, waist_out))
    m = waist_out / waist_in
    z = waist2rayleigh(waist_in, wavelen)
    if len(kwargs) != 1:
        raise ValueError("The function expects only one named argument")
    elif 'dist' in kwargs:
        raise NotImplementedError(filldedent('''
            Currently only focal length is supported as a parameter'''))
    elif 'f' in kwargs:
        f = sympify(kwargs['f'])
        s_in = f * (1 - sqrt(1/m**2 - z**2/f**2))
        s_out = gaussian_conj(s_in, z, f)[0]
    elif 's_in' in kwargs:
        raise NotImplementedError(filldedent('''
            Currently only focal length is supported as a parameter'''))
    else:
        raise ValueError(filldedent('''
            The functions expects the focal length as a named argument'''))
    return (s_in, s_out, f)

#TODO
#def plot_beam():
#    """Plot the beam radius as it propagates in space."""
#    pass

#TODO
#def plot_beam_conjugation():
#    """
#    Plot the intersection of two beams.
#
#    Represents the conjugation relation.
#
#    See Also
#    ========
#
#    conjugate_gauss_beams
#    """
#    pass
