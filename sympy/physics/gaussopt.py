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
"""

from sympy import (atan2, cos, Expr, I, im, Matrix, oo, pi, re, sqrt, sympify,
    together, symbols)
from sympy.utilities.misc import filldedent

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
    [1,  2]
    [3,  4]

    >>> RayTransferMatrix(Matrix([[1, 2], [3, 4]]))
    [1,  2]
    [3,  4]

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

    def __new__(cls, *args):

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
        return Matrix.__new__(cls, temp)

    def __mul__(self, other):
        if isinstance(other, RayTransferMatrix):
            return RayTransferMatrix(Matrix.__mul__(self, other))
        elif isinstance(other, GeometricRay):
            return GeometricRay(Matrix.__mul__(self, other))
        elif isinstance(other, BeamParameter):
            temp = self*Matrix(((other.q,), (1,)))
            q = (temp[0]/temp[1]).expand(complex=True)
            return BeamParameter(other.wavelen, \
                                 together(re(q)), \
                                 z_r = together(im(q)))
        else:
            return Matrix.__mul__(self, other)

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

    See Also
    ========

    RayTransferMatrix

    Examples
    ========

    >>> from sympy.physics.gaussopt import FreeSpace
    >>> from sympy import symbols
    >>> d = symbols('d')
    >>> FreeSpace(d)
    [1, d]
    [0, 1]
    """
    def __new__(cls, d):
        return RayTransferMatrix.__new__(cls, 1, d, 0, 1)

class FlatRefraction(RayTransferMatrix):
    """
    Ray Transfer Matrix for refraction.

    Parameters
    ==========

    n1: refractive index of one medium
    n2: refractive index of other medium

    See Also
    ========

    RayTransferMatrix

    Examples
    ========

    >>> from sympy.physics.gaussopt import FlatRefraction
    >>> from sympy import symbols
    >>> n1, n2 = symbols('n1 n2')
    >>> FlatRefraction(n1, n2)
    [1,     0]
    [0, n1/n2]
    """
    def __new__(cls, n1, n2):
        n1, n2 = sympify((n1, n2))
        return RayTransferMatrix.__new__(cls, 1, 0, 0, n1/n2)

class CurvedRefraction(RayTransferMatrix):
    """
    Ray Transfer Matrix for refraction on curved interface.

    Parameters
    ==========

    R: radius of curvature (positive for concave),
    n1: refractive index of one medium
    n2: refractive index of other medium

    See Also
    ========

    RayTransferMatrix

    Examples
    ========

    >>> from sympy.physics.gaussopt import CurvedRefraction
    >>> from sympy import symbols
    >>> R, n1, n2 = symbols('R n1 n2')
    >>> CurvedRefraction(R, n1, n2)
    [               1,     0]
    [(n1 - n2)/(R*n2), n1/n2]
    """
    def __new__(cls, R, n1, n2):
        R, n1 , n2 = sympify((R, n1, n2))
        return RayTransferMatrix.__new__(cls, 1, 0, (n1-n2)/R/n2, n1/n2)

class FlatMirror(RayTransferMatrix):
    """
    Ray Transfer Matrix for reflection.

    See Also: RayTransferMatrix

    Examples
    ========

    >>> from sympy.physics.gaussopt import FlatMirror
    >>> FlatMirror()
    [1, 0]
    [0, 1]
    """
    def __new__(cls):
        return RayTransferMatrix.__new__(cls, 1, 0, 0, 1)

class CurvedMirror(RayTransferMatrix):
    """
    Ray Transfer Matrix for reflection from curved surface.

    Parameters
    ==========

    radius of curvature (positive for concave)

    See Also
    ========

    RayTransferMatrix

    Examples
    ========

    >>> from sympy.physics.gaussopt import CurvedMirror
    >>> from sympy import symbols
    >>> R = symbols('R')
    >>> CurvedMirror(R)
    [   1, 0]
    [-2/R, 1]
    """
    def __new__(cls, R):
        R = sympify(R)
        return RayTransferMatrix.__new__(cls, 1, 0, -2/R, 1)

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
    def __new__(cls, f):
        f = sympify(f)
        return RayTransferMatrix.__new__(cls, 1, 0, -1/f, 1)


class SingleRightAnglePrism(RayTransferMatrix) :
    """
    Ray Transfer Matrix for a single right andle prism.
    
    Parameters
    ==========

    n: refraction index
    d: prism pass length
    k: beam expansion factor

    See Also
    ========

    RayTransferMatrix
    
    BeamExansionFactor

    Examples
    ========

    >>> from sympy.physics.gaussopt import SingleRightAndlePrism
    >>> from sympy import symbols
    >>> n,d,k = symbols('n d k')
    >>> SingleRightAnglePrism(n,d,k)
    [   k, d/(nk)]
    [   0, 1/k   ]
    """    
    def __new__(cls, n, d, k) :
        n, d, k = sympify((n, d, k))
        return RayTransferMatrix.__new__(cls, k, d/n/k, 0, 1/k)
###
# Representation for geometric ray
###

class GeometricRay(Matrix):
    """
    Representation for a geometric ray in the Ray Transfer Matrix formalism.

    Parameters
    ==========

    height and angle or 2x1 matrix (Matrix(2, 1, [height, angle]))

    Examples
    =======

    >>> from sympy.physics.gaussopt import GeometricRay, FreeSpace
    >>> from sympy import symbols, Matrix
    >>> d, h, angle = symbols('d, h, angle')

    >>> GeometricRay(h, angle)
    [    h]
    [angle]

    >>> FreeSpace(d)*GeometricRay(h, angle)
    [angle*d + h]
    [      angle]

    >>> GeometricRay( Matrix( ((h,), (angle,)) ) )
    [    h]
    [angle]

    See Also
    ========

    RayTransferMatrix

    """

    def __new__(cls, *args):
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

    @property
    def height(self):
        """
        The distance from the optical axis.

        Examples
        ========

        >>> from sympy.physics.gaussopt import GeometricRay
        >>> from sympy import symbols
        >>> h, angle = symbols('h, angle')
        >>> gRay = GeometricRay(h, angle)
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
        >>> gRay = GeometricRay(h, angle)
        >>> gRay.angle
        angle
        """
        return self[1]

class BeamExpansionFactor(Expr) :
    """ Definition of beam expansion factor k = cos(phi) * cos(psi)
    
    Parameters
    ==========
    phi: the angle of incidence
    psi: the angle of refraction
    
    Examples
    ========
    >>> from sympy.physics.gaussopt import BeamExpansionFactor, SingleRightAnglePrism
    >>> from sympy import symbols
    >>> phi, psi = symbols('phi psi')
    >>> bExp = BeamExapansionFactor(phi, psi)
    >>> n,d = symbols('n d')
    >>> SinleRightAnglePrism(n,d,bExp)
    [   bExp, d/(nk)]
    [   0, 1/bExp   ]
    """
    def __new__(cls, phi, psi) :
        phi, psi = sympify((phi, psi))
        return cos(phi)*cos(psi)
    
    @property
    def incidence(self) :
        """
        The angle of icidence

        Examples
        ========

        >>> from sympy.physics.gaussopt import BeamExpansionFactor
        >>> from sympy import symbols
        >>> phi, psi = symbols('phi psi')
        >>> bExp = BeamExpansionFactor(psi, phi)
        >>> gRay.incidence
        phi
        """
        return self.phi
    
    @property
    def refraction(self) :
        """
        The angle of refraction 

        Examples
        ========

        >>> from sympy.physics.gaussopt import BeamExpansionFactor
        >>> from sympy import symbols
        >>> phi, psi = symbols('phi psi')
        >>> bExp = BeamExpansionFactor(psi, phi)
        >>> gRay.refraction
        psi        
        """
        return self.psi
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

def plot_beam(beam, **kwargs):
    """Plot the beam radius as it propagates in space.
    Uses pyglet and ctype libraries

    Parameters
    ==========
    
    beam : BeamParameter for gaussian beam
    z_range : plot range for the beam propagation coordinate

    See Also
    ========    
    
    BeamParameter
    plot_beam2
    
    Examples
    ========
    
    >>> from sympy.physics.gaussopt import BeamParameter, plot_beam
    >>> b = BeamParameter(530e-9, 1, w=1e-5)
    >>> plot_beam(b,z_range=2*b.z_r)    
    """
    
    if len(kwargs) != 1:
        raise ValueError("The function expects only one named argument")    
    elif 'z_range' in kwargs :
        z_range = sympify(kwargs['z_range'])
    else :
        raise ValueError(filldedent('''
            The functions expects the z_range as a named argument'''))

    x = symbols('x')
    # TODO beam.w_0 * 
    # pyglet needs to have better zoom adjustment for geting a normal view
    weist_d = sqrt(1+(x/beam.z_r)**2)
    
    p = Plot(visible=False)
    
    p[1] = weist_d, 'color=black', [x, -z_range, z_range, int(beam.w_0**-1)]
    p[2] = -weist_d, 'color=black', [x, -z_range, z_range, int(beam.w_0**-1)]
    
    p.adjust_all_bounds()
    p.show()
    
def plot_beam2(beam, x_n, **kwargs) :
    """Plot the beam radius as it propagates in space.
    Uses external matplotlib library.

    Parameters
    ==========
    
    beam : BeamParameter for gaussian beam
    x_n : number of samples to plot function
    z_range : plot range for the beam propagation coordinate

    See Also
    ========    
    
    BeamParameter
    plot_beam
    
    Examples
    ========
        
    >>> from sympy.physics.gaussopt import BeamParameter, plot_beam2
    >>> b = BeamParameter(530e-9,1,w=1e-5)
    >>> plot_beam2(b,100,z_range=2*b.z_r)    
    """
    
    if len(kwargs) != 1:
        raise ValueError("The function expects only one named argument")    
    elif 'z_range' in kwargs :
        z_range = sympify(kwargs['z_range'])
    else :
        raise ValueError(filldedent('''
            The functions expects the z_range as a named argument'''))
    
    x = symbols('x')
    
    weist_d = beam.w_0*sqrt(1+(x/beam.z_r)**2)
    
    mplot2d([weist_d, -weist_d], (x, -z_range, z_range, x_n))

def mplot2d(f, var, show=True):
    """
    Plot a 2d function using matplotlib/Tk.
    """

    import warnings
    warnings.filterwarnings("ignore", "Could not match \S")
    
    p = import_module('pylab')
    if not p:
        sys.exit("Matplotlib is required to use mplot2d.")

    if not is_sequence(f):
        f = [f,]
    
    for f_i in f:
        x, y = sample2d(f_i, var)
        p.plot(x, y,'black')

    p.draw()
    
    p.ylabel("Transverse beam cordinate")
    p.xlabel("Propagation coordinate")
    p.grid(True)
    
    if show:
        p.show()

def sample2d(f, x_args):
    """
    Samples a 2d function f over specified intervals and returns two
    arrays (X, Y) suitable for plotting with matlab (matplotlib)
    syntax. See examples\mplot2d.py.

    f is a function of one variable, such as x**2.
    x_args is an interval given in the form (var, min, max, n)
    """
    try:
        f = sympify(f)
    except SympifyError:
        raise ValueError("f could not be interpretted as a SymPy function")
    try:
        x, x_min, x_max, x_n = x_args
    except AttributeError:
        raise ValueError("x_args must be a tuple of the form (var, min, max, n)")

    x_l = float(x_max - x_min)
    x_d = x_l/float(x_n)
    X = arange(float(x_min), float(x_max)+x_d, x_d)

    Y = empty(len(X))
    for i in range(len(X)):
        try:
            Y[i] = float(f.subs(x, X[i]))
        except TypeError:
            Y[i] = None
    return X, Y

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