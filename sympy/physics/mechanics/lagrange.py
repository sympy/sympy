from __future__ import print_function, division

__all__ = ['LagrangesMethod']

from sympy import diff, zeros, Matrix, eye, sympify, Symbol
from sympy.physics.vector import (dynamicsymbols, ReferenceFrame, Point)
from sympy.physics.mechanics.functions import _mat_inv_mul
from sympy.utilities.exceptions import SymPyDeprecationWarning
import warnings

warnings.simplefilter("always", SymPyDeprecationWarning)

class LagrangesMethod(object):
    """Lagrange's method object.

    This object generates the equations of motion in a two step procedure. The
    first step involves the initialization of LagrangesMethod by supplying the
    Lagrangian and a list of the generalized coordinates, at the bare minimum.
    If there are any constraint equations, they can be supplied as keyword
    arguments. The Lagrangian multipliers are automatically generated and are
    equal in number to the constraint equations.Similarly any non-conservative
    forces can be supplied in a list (as described below and also shown in the
    example) along with a ReferenceFrame. This is also discussed further in the
    __init__ method.

    Attributes
    ==========

    mass_matrix : Matrix
        The system's mass matrix

    forcing : Matrix
        The system's forcing vector

    mass_matrix_full : Matrix
        The "mass matrix" for the qdot's, qdoubledot's, and the
        lagrange multipliers (lam)

    forcing_full : Matrix
        The forcing vector for the qdot's, qdoubledot's and
        lagrange multipliers (lam)

    Examples
    ========

    This is a simple example for a one degree of freedom translational
    spring-mass-damper.

    In this example, we first need to do the kinematics.$
    This involves creating generalized coordinates and its derivative.
    Then we create a point and set its velocity in a frame::

        >>> from sympy.physics.mechanics import LagrangesMethod, Lagrangian
        >>> from sympy.physics.mechanics import ReferenceFrame, Particle, Point
        >>> from sympy.physics.mechanics import dynamicsymbols, kinetic_energy
        >>> from sympy import symbols
        >>> q = dynamicsymbols('q')
        >>> qd = dynamicsymbols('q', 1)
        >>> m, k, b = symbols('m k b')
        >>> N = ReferenceFrame('N')
        >>> P = Point('P')
        >>> P.set_vel(N, qd * N.x)

    We need to then prepare the information as required by LagrangesMethod to
    generate equations of motion.
    First we create the Particle, which has a point attached to it.
    Following this the lagrangian is created from the kinetic and potential
    energies.
    Then, a list of nonconservative forces/torques must be constructed, where
    each entry in is a (Point, Vector) or (ReferenceFrame, Vector) tuple, where
    the Vectors represent the nonconservative force or torque.

        >>> Pa = Particle('Pa', P, m)
        >>> Pa.set_potential_energy(k * q**2 / 2.0)
        >>> L = Lagrangian(N, Pa)
        >>> fl = [(P, -b * qd * N.x)]

     Finally we can generate the equations of motion.
     First we create the LagrangesMethod object.To do this one must supply an
     the Lagrangian, the list of generalized coordinates. Also supplied are the
     constraint equations, the forcelist and the inertial frame, if relevant.
     Next we generate Lagrange's equations of motion, such that:
     Lagrange's equations of motion = 0.
     We have the equations of motion at this point.

        >>> l = LagrangesMethod(L, [q], forcelist = fl, frame = N)
        >>> print(l.form_lagranges_equations())
        Matrix([[b*Derivative(q(t), t) + 1.0*k*q(t) + m*Derivative(q(t), t, t)]])

    We can also solve for the states using the 'rhs' method.

        >>> print(l.rhs())
        Matrix([[Derivative(q(t), t)], [(-b*Derivative(q(t), t) - 1.0*k*q(t))/m]])

    Please refer to the docstrings on each method for more details.

    """

    def __init__(self, Lagrangian, q_list, coneqs=None, forcelist=None, frame=None):
        """Supply the following for the initialization of LagrangesMethod

        Lagrangian : Sympifyable

        q_list : list
            A list of the generalized coordinates

        coneqs : list
            A list of the holonomic and non-holonomic constraint equations.
            VERY IMPORTANT NOTE- The holonomic constraints must be
            differentiated with respect to time and then included in coneqs.

        forcelist : list
            Takes a list of (Point, Vector) or (ReferenceFrame, Vector) tuples
            which represent the force at a point or torque on a frame. This
            feature is primarily to account for the nonconservative forces
            amd/or moments.

        frame : ReferenceFrame
            Supply the inertial frame. This is used to determine the
            generalized forces due to non-sonservative forces.

        """

        self._L = sympify(Lagrangian)
        self.eom = None  # initializing the eom Matrix
        self._m_cd = Matrix([])  # Mass Matrix of differentiated coneqs
        self._m_d = Matrix([])  # Mass Matrix of dynamic equations
        self._f_cd = Matrix([])  # Forcing part of the diff coneqs
        self._f_d = Matrix([])  # Forcing part of the dynamic equations
        self.lam_coeffs = Matrix([])  # Initializing the coeffecients of lams

        self.forcelist = forcelist
        self.inertial = frame

        self.lam_vec = Matrix([])

        self._term1 = Matrix([])
        self._term2 = Matrix([])
        self._term3 = Matrix([])
        self._term4 = Matrix([])

        # Creating the qs, qdots and qdoubledots

        q_list = list(q_list)
        if not isinstance(q_list, list):
            raise TypeError('Generalized coords. must be supplied in a list')
        self._q = q_list
        self._qdots = [diff(i, dynamicsymbols._t) for i in self._q]
        self._qdoubledots = [diff(i, dynamicsymbols._t) for i in self._qdots]

        self.coneqs = coneqs

    def form_lagranges_equations(self):
        """Method to form Lagrange's equations of motion.

        Returns a vector of equations of motion using Lagrange's equations of
        the second kind.

        """

        q = self._q
        qd = self._qdots
        qdd = self._qdoubledots
        n = len(q)

        #Putting the Lagrangian in a Matrix
        L = Matrix([self._L])

        #Determining the first term in Lagrange's EOM
        self._term1 = L.jacobian(qd)
        self._term1 = ((self._term1).diff(dynamicsymbols._t)).transpose()

        #Determining the second term in Lagrange's EOM
        self._term2 = (L.jacobian(q)).transpose()

        #term1 and term2 will be there no matter what so leave them as they are

        if self.coneqs is not None:
            coneqs = self.coneqs
            #If there are coneqs supplied then the following will be created
            coneqs = list(coneqs)
            if not isinstance(coneqs, list):
                raise TypeError('Enter the constraint equations in a list')

            o = len(coneqs)

            #Creating the multipliers
            self.lam_vec = Matrix(dynamicsymbols('lam1:' + str(o + 1)))

            #Extracting the coeffecients of the multipliers
            coneqs_mat = Matrix(coneqs)
            qd = self._qdots
            self.lam_coeffs = -coneqs_mat.jacobian(qd)

            #Determining the third term in Lagrange's EOM
            #term3 = ((self.lam_vec).transpose() * self.lam_coeffs).transpose()
            self._term3 = self.lam_coeffs.transpose() * self.lam_vec

            #Taking the time derivative of the constraint equations
            diffconeqs = [diff(i, dynamicsymbols._t) for i in coneqs]

            #Extracting the coeffecients of the qdds from the diff coneqs
            diffconeqs_mat = Matrix(diffconeqs)
            qdd = self._qdoubledots
            self._m_cd = diffconeqs_mat.jacobian(qdd)

            #The remaining terms i.e. the 'forcing' terms in diff coneqs
            qddzero = dict(zip(qdd, [0] * n))
            self._f_cd = -diffconeqs_mat.subs(qddzero)

        else:
            self._term3 = zeros(n, 1)

        if self.forcelist is not None:
            forcelist = self.forcelist
            N = self.inertial
            if not isinstance(N, ReferenceFrame):
                raise TypeError('Enter a valid ReferenceFrame')
            if not isinstance(forcelist, (list, tuple)):
                raise TypeError('Forces must be supplied in a list of: lists'
                        ' or tuples')
            self._term4 = zeros(n, 1)
            for i, v in enumerate(qd):
                for j, w in enumerate(forcelist):
                    if isinstance(w[0], ReferenceFrame):
                        speed = w[0].ang_vel_in(N)
                        self._term4[i] += speed.diff(v, N) & w[1]
                    elif isinstance(w[0], Point):
                        speed = w[0].vel(N)
                        self._term4[i] += speed.diff(v, N) & w[1]
                    else:
                        raise TypeError('First entry in force pair is a point'
                                        ' or frame')

        else:
            self._term4 = zeros(n, 1)

        self.eom = self._term1 - self._term2 - self._term3 - self._term4

        return self.eom

    @property
    def mass_matrix(self):
        """ Returns the mass matrix, which is augmented by the Lagrange
        multipliers, if necessary.

        If the system is described by 'n' generalized coordinates and there are
        no constraint equations then an n X n matrix is returned.

        If there are 'n' generalized coordinates and 'm' constraint equations
        have been supplied during initialization then an n X (n+m) matrix is
        returned. The (n + m - 1)th and (n + m)th columns contain the
        coefficients of the Lagrange multipliers.

        """

        if self.eom is None:
            raise ValueError('Need to compute the equations of motion first')

        #The 'dynamic' mass matrix is generated by the following
        self._m_d = (self.eom).jacobian(self._qdoubledots)

        if len(self.lam_coeffs) != 0:
            return (self._m_d).row_join((self.lam_coeffs).transpose())
        else:
            return self._m_d

    @property
    def mass_matrix_full(self):
        """ Augments the coefficients of qdots to the mass_matrix. """

        n = len(self._q)
        if self.eom is None:
            raise ValueError('Need to compute the equations of motion first')
        #THE FIRST TWO ROWS OF THE MATRIX
        row1 = eye(n).row_join(zeros(n, n))
        row2 = zeros(n, n).row_join(self.mass_matrix)
        if self.coneqs is not None:
            m = len(self.coneqs)
            I = eye(n).row_join(zeros(n, n + m))
            below_eye = zeros(n + m, n)
            A = (self.mass_matrix).col_join((self._m_cd).row_join(zeros(m, m)))
            below_I = below_eye.row_join(A)
            return I.col_join(below_I)
        else:
            A = row1.col_join(row2)
            return A

    @property
    def forcing(self):
        """ Returns the forcing vector from 'lagranges_equations' method. """

        if self.eom is None:
            raise ValueError('Need to compute the equations of motion first')

        qdd = self._qdoubledots
        qddzero = dict(zip(qdd, [0] * len(qdd)))

        if self.coneqs is not None:
            lam = self.lam_vec
            lamzero = dict(zip(lam, [0] * len(lam)))

            #The forcing terms from the eoms
            self._f_d = -((self.eom).subs(qddzero)).subs(lamzero)

        else:
            #The forcing terms from the eoms
            self._f_d = -(self.eom).subs(qddzero)

        return self._f_d

    @property
    def forcing_full(self):
        """ Augments qdots to the forcing vector above. """

        if self.eom is None:
            raise ValueError('Need to compute the equations of motion first')
        if self.coneqs is not None:
            return (Matrix(self._qdots)).col_join((self.forcing).col_join(self._f_cd))
        else:
            return (Matrix(self._qdots)).col_join(self.forcing)

    def rhs(self, inv_method=None, **kwargs):
        """ Returns equations that can be solved numerically

        Parameters
        ==========

        inv_method : str
            The specific sympy inverse matrix calculation method to use. For a
            list of valid methods, see :py:method:
            `~sympy.matrices.matrices.MatrixBase.inv`

        """

        if 'method' in kwargs:
            #The method kwarg is deprecated in favor of inv_method.
            SymPyDeprecationWarning(feature="method kwarg",
                    useinstead="inv_method kwarg",
                    deprecated_since_version="0.7.6").warn()
            #For now accept both
            inv_method = kwargs['method']

        if inv_method is None:
            self._rhs = _mat_inv_mul(self.mass_matrix_full,
                                          self.forcing_full)
        else:
            self._rhs = (self.mass_matrix_full.inv(inv_method,
                         try_block_diag=True) * self.forcing_full)
        return self._rhs
