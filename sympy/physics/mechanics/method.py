from __future__ import annotations
from abc import ABC, abstractmethod
from sympy import linear_eq_to_matrix

__all__ = ['MethodBase']


class MethodBase(ABC):
    """Abstract Base Class for all methods for forming the equations of motion
    of multiody systems.

    Minimal coordinate equations of motion take this form as first order
    ordinary differential equations.

    Kinematical differential equations:

    .. code:: text

        Mk(q, t) q' = Fk(u, q, t)

    Dynamical differential equations:

    .. code:: text

        Md(q, t) u' = Fd(u, q, t)

    Combined they take this full form:

    .. code:: text

        M = [Mk 0] [q'] = [Fk] = F
            [0 Md] [u']   [Fd]

    If there are additional holonomic or nonholonomic constraints, these
    equations can be augmented with the Lagrange multipliers:

    .. code:: text

        M = [Mk 0  0  ] [q'] = [Fk] = F
            [0  Md MjT] [u']   [Fd]
            [0  Mj 0  ] [j']   [Fj]

    where ``Mj`` is the "Jacobian of the constraints" and ``MjT`` is its
    transpose. ``j'=lambda`` are the Lagrange multipliers representing the the
    constraint forces and then ``j`` are generalized impulses of these forces.

    The equations of motion can also be augmented to reveal any noncontributing
    force and take this form:

    .. code:: text

        M = [Mk 0  0] [q'] = [Fk] = F
            [0  Md 0] [u']   [Fd]
            [0  Mj I] [j']   [Fj]

    With ``j'`` being the measure number of a noncontributing force.

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
    # TODO : Add to all methods or only have on Lagrange/TMT?
    # TODO : This is not a great name. Better is lam or lagrange_multipliers or
    # lambda or even l.
    #@property
    #@abstractmethod
    def _lam_vec(self):
        """Column matrix of Lagrange multipliers."""
        pass

    # KanesMethod: (init) bodies, (attr) bodies & bodylist [deprecated], (kanes_equations) bodies
    # LagrangesMethod: (init) bodies, (attr) bodies
    # JointsMethod: (attr) bodies
    # System: (attr) bodies
    @property
    @abstractmethod
    def bodies(self):
        """List of :py:class:`~sympy.physics.mechanics.particle.Particle`,
        :py:class:`~sympy.physics.mechanics.rigidbody.RigidBody`, or
        :py:class:`~sympy.physics.mechanics.body.Body` objects that make up the
        multibody system."""
        pass

    # KanesMethod: (init) forcelist, (attr) loads & forcelist [deprecated], (kanes_equations) loads
    # LagrangesMethod: (init) forcelist, (attr) loads & forcelist [deprecated], does not include conservative forces
    # JointsMethod: (attr) loads
    # System: (attr) loads
    @property
    @abstractmethod
    def loads(self):
        """List of :py:class:`~sympy.physics.mechanics.loads.Force`,
        :py:class:`~sympy.physics.mechanics.loads.Torque`,
        tuple(:py:class:`~sympy.physics.vector.point.Point`,
        :py:class:`~sympy.physics.vector.vector.Vector`),
        tuple(:py:class:`~sympy.physics.vector.frame.ReferenceFrame`,
        :py:class:`~sympy.physics.vector.vector.Vector`) loads applied to
        multibody system."""
        pass

    # KanesMethod: (init) configuration_constraints, (attr) _f_h, holonomic_constraints
    # LagrangesMethod: (init) hol_coneqs, (attr) _hol_coneqs, holonomic_constraints
    # JointsMethod: (attr) holonomic_constraints
    # System: (attr) holonomic_constraints
    @property
    @abstractmethod
    def holonomic_constraints(self):
        """M x 1 column matrix of holonomic configuration constraint residual
        expressions ``fh`` where:

        .. code:: text

            fh(q, t) = 0

        """
        pass

    # KanesMethod: (init, attr) nonholonomic_constraints
    # LagrangesMethod: (init) nonhol_coneqs, (attr) coneqs[M:], nonholonomic_constraints
    # JointsMethod: (attr) nonholonomic_constraints
    # System: (attr) nonholonomic_constraints
    @property
    @abstractmethod
    def nonholonomic_constraints(self):
        """m x 1 column matrix of nonholonomic residual expressions ``fn``
        where:

        .. code:: text

            fn(q', q, t) = fn(u, q, t) = 0

        """
        pass

    # KanesMethod: (init) velocity_constraints, (attr) _k_nh, _f_nh, velocity_constraints
    # LagrangesMethod: (init) hol_coneqs & nonhol_eqs, (attr) coneqs, velocity_constraints
    # JointsMethod: (attr) velocity_constraints
    # System: (attr) velocity_constraints
    @property
    @abstractmethod
    def velocity_constraints(self):
        """m + M x 1 column matrix of motion/velocity constraint residual
        expressions ``fv`` where:

        .. code:: text

            fv(q', q, t) = [fh'(q', q, t)] = fv(u, q, t) = [fh'(u, q, t)] = 0
                           [fn(q', q, t) ]                 [fn(u, q, t) ]

        The time differentiated holonomic configuration constraints should be
        stacked on top of the nonholonomic constraints.

        """
        pass

    # KanesMethod: (init) acceleration_constraints, (attr) _f_dnh, _k_dnh, acceleration_constraints
    # LagrangesMethod: (attr) _m_cd, _f_cd, acceleration_constraints
    # JointsMethod: (attr) acceleration_constraints
    # System: (attr) acceleration_constraints
    @property
    @abstractmethod
    def acceleration_constraints(self):
        """m + M x 1 column matrix of acceleration constraint residual
        expressions ``fa`` where:

        .. code:: text

            fv' = fa(q'', q', q, t) = fa(u', u, q, t) = 0

        The twice time differentiated holonomic configuration constraints
        should be stacked on top of time differentieated nonholonomic
        constraints.

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
        """Linear coefficient matrix ``Md`` for the second time derivative of
        the coordinates or the first time derivative of the speeds:

        .. code:: text

            Md q'' = Md u' = Fd

        """
        pass

    # KanesMethod: (attr) forcing
    # LagrangesMethod: (attr) forcing
    # JointsMethod: (attr) returns the *Method's method
    # System: (attr) returns the *Method's method
    @property
    @abstractmethod
    def forcing(self):
        """Nonlinear forcing terms ``Fd`` in the dynamical differential
        equations:

        .. code:: text

            Md q'' = Md u' = Fd

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
        """Linear coefficient matrix ``M`` for the full first order form of the
        equations of motion:

        .. code:: text

            M [q' ] = M [q'] =  F
              [q'']     [u']

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
        """Nonlinear forcing terms ``F`` in the full first order form of the
        equations of motion:

        .. code:: text

            M [q' ] = M [q'] =  F
              [q'']     [u']

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
        """Returns an M + m x N coefficient matrix ``C`` which is the Jacobian
        of the constraints.

        .. code:: text

            fv = C*q' + gv(q, t) = C*u + gv(q, t) = 0

        """
        C, _ = linear_eq_to_matrix(self.velocity_constraints, self.u[:])
        return C

    def rhs(self, inv_method=None, **kwargs):
        """Returns the right hand side of the full first order form of the
        equations of motion in explicit form:

        .. code:: text

            rhs(u, q, t) = Inv(M) F

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

    @abstractmethod
    def _form_eoms(self):
        raise NotImplementedError("Subclasses must implement this.")
