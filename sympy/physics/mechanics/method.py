from __future__ import annotations
from abc import ABC, abstractmethod
from sympy import linear_eq_to_matrix

__all__ = ['Method']


class Method(ABC):
    """Abstract Base Class for all methods for forming the equations of motion
    of multiody systems.

    Minimal coordinate equations of motion take this first order form:

    N(q, t) q' = G(u, q, t)
    M(q, t) u' = F(u, q, t)


    """
    # System: (init, attr) frame
    # KanesMethod: (init) frame, (attr) _inertial & frame
    # LagrangesMethod: (init) frame, (attr) inertial & frame
    # JointsMethod: (init) newtonian, (attr) frame
    # TODO : Deprecate LagrangesMethod.inertial
    @property
    @abstractmethod
    def frame(self):
        """Inertial reference frame that the equations of motion were
        formulated with respect to."""
        pass

    # KanesMethod: (init) q_ind & q_dependent, (attr) q
    # LagrangesMethod: (init) qs, (attr) q
    # JointsMethod: (attr) q
    # System: (attr) q
    @property
    @abstractmethod
    def q(self):
        """Column matrix of N functions of time that represent the multibody
        system's (generalized) coordinates."""
        pass

    # KanesMethod: (init) u_ind & u_dependent, (attr) u
    # LagrangesMethod: (init) NA, (attr) u
    # JointsMethod: (attr) u
    # System: (attr) u
    # TODO : Lagranges should have the option to use defined symbols for u. It
    # does not really put the equations in first order form now.
    @property
    @abstractmethod
    def u(self):
        """Column matrix of N functions of time that represent the multibody
        system's (generalized) speeds."""
        pass

    # KanesMethod: NA
    # LagrangesMethod: (attr) lam_vec
    # JointsMethod: NA
    # System: NA
    # TODO : Add to all methods.
    # TODO : This is not a great name. Better is lam or lagrange_multipliers or
    # lambda or even l.
    #@property
    #@abstractmethod
    def lam_vec(self):
        """Column matrix of Lagrange multipliers."""
        pass

    # KanesMethod: (init) bodies, (attr) bodies & bodylist [deprecated], (kanes_equations) bodies
    # LagrangesMethod: (init) bodies, (attr) bodies
    # JointsMethod: (attr) bodies
    # System: (attr) bodies
    @property
    @abstractmethod
    def bodies(self):
        """List of :py:class:`Particle`, :py:class:`RigidBody`, or
        :py:class:`Body` objects that make up the multibody system."""
        pass

    # KanesMethod: (init) forcelist, (attr) loads & forcelist [deprecated], (kanes_equations) loads
    # LagrangesMethod: (init) forcelist, (attr) loads & forcelist [deprecated], does not include conservative forces
    # JointsMethod: (attr) loads
    # System: (attr) loads
    @property
    @abstractmethod
    def loads(self):
        """List of :py:class:`Force`, :py:class:`Torque`, tuple(Point, Vector),
        tuple(ReferenceFrame, Vector) loads applied to multibody system."""
        pass

    # KanesMethod: (init) configuration_constraints, (attr) _f_h, holonomic_constraints
    # LagrangesMethod: (init) hol_coneqs, (attr) _hol_coneqs, holonomic_constraints
    # JointsMethod: (attr) holonomic_constraints
    # System: (attr) holonomic_constraints
    @property
    @abstractmethod
    def holonomic_constraints(self):
        """Column matrix of shape(M, 1) of configuration constraint residuals f
        where:

        f(q, t) = 0

        """
        pass

    # KanesMethod: (init, attr) nonholonomic_constraints
    # LagrangesMethod: (init) nonhol_coneqs, (attr) coneqs[M:], nonholonomic_constraints
    # JointsMethod: (attr) nonholonomic_constraints
    # System: (attr) nonholonomic_constraints
    @property
    @abstractmethod
    def nonholonomic_constraints(self):
        """Column matrix of shape(m, 1) nonholonomic residuals f where:

        f(q', q, t) = 0

        or

        f(u, q, t) = 0

        """
        pass

    # KanesMethod: (init) velocity_constraints, (attr) _k_nh, _f_nh, velocity_constraints
    # LagrangesMethod: (init) hol_coneqs & nonhol_eqs, (attr) coneqs, velocity_constraints
    # JointsMethod: (attr) velocity_constraints
    # System: (attr) velocity_constraints
    @property
    @abstractmethod
    def velocity_constraints(self):
        """Column matrix of shape(m + M, 1) motion constraint residuals f
        where::

            f(q', q, t) = [fh'(q', q, t)] = 0
                          [fn(q', q, t) ]

        or::

            f(u, q, t) = [fh'(u, q, t)] = 0
                         [fn(u, q, t) ]

        This includes the time differentiated configuration constraints.

        """
        pass

    # KanesMethod: (init) acceleration_constraints, (attr) _f_dnh, _k_dnh, acceleration_constraints
    # LagrangesMethod: (attr) _m_cd, _f_cd, acceleration_constraints
    # JointsMethod: (attr) acceleration_constraints
    # System: (attr) acceleration_constraints
    @property
    @abstractmethod
    def acceleration_constraints(self):
        """Column matrix of shape(m, 1) or shape(m + M, 1) acceleration
        constraint residuals f where:

        f(q'', q', q, t) = 0

        or

        f(u', u, q, t) = 0

        """
        pass

    # KanesMethod: (attr) mass_matrix
    # LagrangesMethod: (attr) mass_matrix, but can also be [Md C.T]
    # JointsMethod: (attr) returns the *Method's method
    # System: (attr) returns the *Method's method
    # TODO : Deprecate the augmented form in LagrangesMethod
    @property
    @abstractmethod
    def mass_matrix(self):
        """Linear coefficient matrix for the second time derivative of the
        coordiantes or the first time derivative of the speeds.

        M in M*q'' = F

        or

        M in M*u' = F

        """
        pass

    # KanesMethod: (attr) forcing
    # LagrangesMethod: (attr) forcing
    # JointsMethod: (attr) returns the *Method's method
    # System: (attr) returns the *Method's method
    @property
    @abstractmethod
    def forcing(self):
        """Terms that are not linear in the the second time derivative of the
        coordiantes or the first time derivative of the speeds.

        F in M*q'' = F

        or

        F in M*u' = F

        """
        pass


    # KanesMethod: (attr) mass_matrix_full
    # LagrangesMethod: (attr) mass_matrix_full, same as Kane's except if there
    # are constraints, then is gives the fully augmented form
    # JointsMethod: (attr) returns the *Method's method
    # System: (attr) returns the *Method's method
    # TODO : Deprecate the augmented form in LagrangesMethod
    @property
    @abstractmethod
    def mass_matrix_full(self):
        """Linear coefficient matrix for the first order form.

        M in M*[q' ] = F
               [q'']

        or

        M in M*[q'] = F
               [u']


        """
        pass

    # KanesMethod: (attr) forcing_full
    # LagrangesMethod: (attr) forcing_full, same as Kane's except if there
    # are constraints, then is gives the fully augmented form
    # JointsMethod: (attr) returns the *Method's method
    # System: (attr) returns the *Method's method
    # TODO : Deprecate the augmented form in LagrangesMethod
    @property
    @abstractmethod
    def forcing_full(self):
        """Linear coefficient matrix for the first order form.

        M in M*[q' ] = F
               [q'']

        or

        M in M*[q'] = F
               [u']
        """
        pass

    @property
    #@abstractmethod
    def _mass_matrix_aug(self):
        """Returns the mass matrix augmented when Lagrange multipliers and/or
        noncontributing forces are included in the system.

        M in M*[q''] = F
               [lam]

        M in M*[u' ] = F
               [lam]

        """
        """
        [Mk 0  0    0 ][qd ] = [Fk]
        [0  Md Ml.T 0 ][ud ]   [Fd]
        [0  Ml 0    0 ][lam]   [Fl]
        [0  Mu 0    Mf][f  ]   [Ff]

        """
        pass

    @property
    #@abstractmethod
    def _forcing_aug(self):
        """Returns the augmented forcing vector when Lagrange multipliers
        and/or noncontributing forces are included in the system.

        F in M*[q''] = F
               [lam]

        F in M*[u' ] = F
               [lam]
        """
        pass

    # KanesMethod: NA
    # LagrangesMethod: (attr) -lam_coeffs
    # JointsMethod: NA
    # System: NA
    def constraints_jacobian(self):
        """Returns a shape(M + m, n) coefficient matrix for the motion level
        constraints.

        C in C*q' + f(q, t) = 0

        or

        C in C*u + f(q, t) = 0

        """
        C, _ = linear_eq_to_matrix(self.velocity_constraints, self.u[:])
        return C

    def rhs(self, inv_method=None, **kwargs):
        """Returns equations that can be solved numerically.

        Parameters
        ==========

        inv_method : str
            The specific sympy inverse matrix calculation method to use. For a
            list of valid methods, see
            :meth:`~sympy.matrices.matrixbase.MatrixBase.inv`

        """

        if inv_method is None:
            self._rhs = self.mass_matrix_full.LUsolve(self.forcing_full)
        else:
            self._rhs = (self.mass_matrix_full.inv(inv_method,
                         try_block_diag=True) * self.forcing_full)
        return self._rhs

    def _form_eoms(self):
        raise NotImplementedError("Subclasses must implement this.")
