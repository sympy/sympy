from __future__ import print_function, division

__all__ = ['LagrangesMethod']

from sympy import diff, zeros, Matrix, eye, sympify
from sympy.physics.vector import (dynamicsymbols, ReferenceFrame, Point)
from sympy.physics.mechanics.functions import _mat_inv_mul, \
        _find_dynamicsymbols, _subs_keep_derivs
from sympy.physics.mechanics.linearize import Linearizer
from sympy.utilities import default_sort_key
from sympy.utilities.exceptions import SymPyDeprecationWarning
from sympy.utilities.misc import filldedent
import collections
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

    def __init__(self, Lagrangian, q_list, coneqs=None, forcelist=None,
            frame=None, hol_coneqs=None, nonhol_coneqs=None):
        """Supply the following for the initialization of LagrangesMethod

        Lagrangian : Sympifyable

        q_list : list
            A list of the generalized coordinates

        hol_coneqs: list
            A list of the holonomic constraint equations

        nonhol_coneqs: list
            A list of the nonholonomic constraint equations

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

        # Deal with constraint equations
        if coneqs:
            SymPyDeprecationWarning(filldedent("""The `coneqs` kwarg is
            deprecated in favor of `hol_coneqs` and `nonhol_coneqs`. Please
            update your code""")).warn()
            self.coneqs = coneqs
        else:
            mat_build = lambda x: Matrix(x) if x else Matrix()
            hol_coneqs = mat_build(hol_coneqs)
            nonhol_coneqs = mat_build(nonhol_coneqs)
            self.coneqs = Matrix([hol_coneqs.diff(dynamicsymbols._t),
                    nonhol_coneqs])
            self._hol_coneqs = hol_coneqs

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

        if self.coneqs:
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
        if self.coneqs:
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

        if self.coneqs:
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
        if self.coneqs:
            return (Matrix(self._qdots)).col_join((self.forcing).col_join(self._f_cd))
        else:
            return (Matrix(self._qdots)).col_join(self.forcing)

    def to_linearizer(self, q_ind=None, qd_ind=None, q_dep=None, qd_dep=None):
        """Returns an instance of the Linearizer class, initiated from the
        data in the LagrangesMethod class. This may be more desirable than using
        the linearize class method, as the Linearizer object will allow more
        efficient recalculation (i.e. about varying operating points).

        Parameters
        ----------
        q_ind, qd_ind : array_like, optional
            The independent generalized coordinates and speeds.
        q_dep, qd_dep : array_like, optional
            The dependent generalized coordinates and speeds.
        """

        # Compose vectors
        t = dynamicsymbols._t
        q = Matrix(self._q)
        u = Matrix(self._qdots)
        ud = u.diff(t)
        # Get vector of lagrange multipliers
        lams = self.lam_vec

        mat_build = lambda x: Matrix(x) if x else Matrix()
        q_i = mat_build(q_ind)
        q_d = mat_build(q_dep)
        u_i = mat_build(qd_ind)
        u_d = mat_build(qd_dep)

        # Compose general form equations
        f_c = self._hol_coneqs
        f_v = self.coneqs
        f_a = f_v.diff(t)
        f_0 = u
        f_1 = -u
        f_2 = self._term1
        f_3 = -(self._term2 + self._term4)
        f_4 = -self._term3

        # Check that there are an appropriate number of independent and
        # dependent coordinates
        if len(q_d) != len(f_c) or len(u_d) != len(f_v):
            raise ValueError(("Must supply {:} dependent coordinates, and " +
                    "{:} dependent speeds").format(len(f_c), len(f_v)))
        if set(Matrix([q_i, q_d])) != set(q):
            raise ValueError("Must partition q into q_ind and q_dep, with " +
                    "no extra or missing symbols.")
        if set(Matrix([u_i, u_d])) != set(u):
            raise ValueError("Must partition qd into qd_ind and qd_dep, " +
                    "with no extra or missing symbols.")

        # Find all other dynamic symbols, forming the forcing vector r.
        # Sort r to make it canonical.
        insyms = set(Matrix([q, u, ud, lams]))
        r = list(_find_dynamicsymbols(f_3, insyms))
        r.sort(key=default_sort_key)
        # Check for any derivatives of variables in r that are also found in r.
        for i in r:
            if diff(i, dynamicsymbols._t) in r:
                raise ValueError('Cannot have derivatives of specified \
                                 quantities when linearizing forcing terms.')

        return Linearizer(f_0, f_1, f_2, f_3, f_4, f_c, f_v, f_a, q, u, q_i,
                q_d, u_i, u_d, r, lams)

    def linearize(self, q_ind=None, qd_ind=None, q_dep=None, qd_dep=None,
            **kwargs):
        """ Linearize the equations of motion about a symbolic operating point.

        If kwarg A_and_B is False (default), returns M, A, B, r for the
        linearized form, M*[q', u']^T = A*[q_ind, u_ind]^T + B*r.

        If kwarg A_and_B is True, returns A, B, r for the linearized form
        dx = A*x + B*r, where x = [q_ind, u_ind]^T. Note that this is
        computationally intensive if there are many symbolic parameters. For
        this reason, it may be more desirable to use the default A_and_B=False,
        returning M, A, and B. Values may then be substituted in to these
        matrices, and the state space form found as
        A = P.T*M.inv()*A, B = P.T*M.inv()*B, where P = Linearizer.perm_mat.

        In both cases, r is found as all dynamicsymbols in the equations of
        motion that are not part of q, u, q', or u'. They are sorted in
        canonical form.

        The operating points may be also entered using the `op_point` kwarg.
        This takes a dictionary of {symbol: value}, or a an iterable of such
        dictionaries. The values may be numberic or symbolic. The more values
        you can specify beforehand, the faster this computation will run.

        For more documentation, please see the `Linearizer` class."""

        linearizer = self.to_linearizer(q_ind, qd_ind, q_dep, qd_dep)
        result = linearizer.linearize(**kwargs)
        return result + (linearizer.r,)

    def solve_multipliers(self, op_point=None, sol_type='dict'):
        """Solves for the values of the lagrange multipliers symbolically at
        the specified operating point

        Parameters
        ----------
        op_point : dict or iterable of dicts, optional
            Point at which to solve at. The operating point is specified as
            a dictionary or iterable of dictionaries of {symbol: value}. The
            value may be numeric or symbolic itself.

        sol_type : str, optional
            Solution return type. Valid options are:
            - 'dict': A dict of {symbol : value} (default)
            - 'Matrix': An ordered column matrix of the solution
        """

        # Determine number of multipliers
        k = len(self.lam_vec)
        if k == 0:
            raise ValueError("System has no lagrange multipliers to solve for.")
        # Compose dict of operating conditions
        if isinstance(op_point, dict):
            op_point_dict = op_point
        elif isinstance(op_point, collections.Iterable):
            op_point_dict = {}
            for op in op_point:
                op_point_dict.update(op)
        else:
            op_point_dict = {}
        # Compose the system to be solved
        mass_matrix = self.mass_matrix.col_join((-self.lam_coeffs.row_join(
                zeros(k, k))))
        force_matrix = self.forcing.col_join(self._f_cd)
        # Sub in the operating point
        mass_matrix = _subs_keep_derivs(mass_matrix, op_point_dict)
        force_matrix = _subs_keep_derivs(force_matrix, op_point_dict)
        # Solve for the multipliers
        sol_list = mass_matrix.LUsolve(-force_matrix)[-k:]
        if sol_type == 'dict':
            return dict(zip(self.lam_vec, sol_list))
        elif sol_type == 'Matrix':
            return Matrix(sol_list)
        else:
            raise ValueError("Unknown sol_type {:}.".format(sol_type))

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
