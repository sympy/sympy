"""Implementations of musculotendon models.

Musculotendon models are a critical component of biomechanical models, one that
differentiates them from pure multibody systems. Musculotendon models produce a
force dependent on their level of activation, their length, and their
extension velocity. Length- and extension velocity-dependent force production
are governed by force-length and force-velocity characteristics.
These are normalized functions that are dependent on the musculotendon's state
and are specific to a given musculotendon model.

"""

from enum import IntEnum, unique

from sympy.core.numbers import Integer, Rational
from sympy.core.symbol import Symbol
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.trigonometric import sin
from sympy.matrices.dense import MutableDenseMatrix as Matrix, diag, eye, zeros
from sympy.physics._biomechanics.activation import ActivationBase
from sympy.physics._biomechanics.curve import (
    FiberForceLengthActiveDeGroote2016,
    FiberForceLengthPassiveDeGroote2016,
    FiberForceVelocityDeGroote2016,
    TendonForceLengthDeGroote2016,
)
from sympy.physics._biomechanics._mixin import _NamedMixin
from sympy.physics.mechanics.actuator import ForceActuator
from sympy.physics.vector.functions import dynamicsymbols


__all__ = [
    'MusculotendonDeGroote2016',
    'MusculotendonFormulation',
]


@unique
class MusculotendonFormulation(IntEnum):
    """Enumeration of types of musculotendon dynamics formulations.

    Explanation
    ===========

    An (integer) enumeration is used as it allows for clearer selection of the
    different formulations of musculotendon dynamics.

    Members
    =======

    RIGID_TENDON : 0
        Uses a rigid tendon model.
    FIBER_LENGTH_EXPLICIT : 1
        Uses an explicit elastic tendon model with the muscle fiber length
        (l_M) as the state variable.
    TENDON_FORCE_EXPLICIT : 2
        Uses an explicit elastic tendon model with the tendon force (F_T) as
        the state variable.
    FIBER_LENGTH_IMPLICIT : 3
        Uses an implicit elastic tendon model with the muscle fiber length
        (l_M) as the state variable and the muscle fiber velocity as an
        additional input variable.
    TENDON_FORCE_IMPLICIT : 4
        Uses an implicit elastic tendon model with the tendon force (F_T) as
        the state variable as the muscle fiber velocity as an additional input
        variable.

    """

    RIGID_TENDON = 0
    FIBER_LENGTH_EXPLICIT = 1
    TENDON_FORCE_EXPLICIT = 2
    FIBER_LENGTH_IMPLICIT = 3
    TENDON_FORCE_IMPLICIT = 4


_DEFAULT_MUSCULOTENDON_FORMULATION = MusculotendonFormulation.RIGID_TENDON


class MusculotendonDeGroote2016(ForceActuator, _NamedMixin):
    """Abstract base class for all musculotendon classes to inherit from.

    Explanation
    ===========

    A musculotendon generates a contractile force based on its activation,
    length, and shortening velocity. This abstract base class is to be inherited
    by all musculotendon subclasses that implement different characteristic
    musculotendon curves. Characteristic musculotendon curves are required for
    the tendon force-length, passive fiber force-length, active fiber force-
    length, and fiber force-velocity relationships.

    Parameters
    ==========
    name : str
        The name identifier associated with the musculotendon. This name is used
        as a suffix when automatically generated symbols are instantiated. It
        must be a string of nonzero length.
    pathway : PathwayBase
        The pathway that the actuator follows. This must be an instance of a
        concrete subclass of ``PathwayBase``, e.g. ``LinearPathway``.
    activation_dynamics : ActivationBase

    musculotendon_dynamics : MusculotendonFormulation | int

    tendon_slack_length :

    peak_isometric_force :

    optimal_fiber_length :

    maximal_fiber_velocity :

    optimal_pennation_angle :

    fiber_damping_coefficient :

    """

    def __init__(
        self,
        name,
        pathway,
        activation_dynamics,
        *,
        musculotendon_dynamics=_DEFAULT_MUSCULOTENDON_FORMULATION,
        tendon_slack_length=None,
        peak_isometric_force=None,
        optimal_fiber_length=None,
        maximal_fiber_velocity=None,
        optimal_pennation_angle=None,
        fiber_damping_coefficient=None,
    ):
        self.name = name

        # Supply a placeholder force to the super initializer, this will be
        # replaced later
        super().__init__(Symbol('F'), pathway)

        # Activation dynamics
        if not isinstance(activation_dynamics, ActivationBase):
            msg = (
                f'Can\'t set attribute `activation_dynamics` to '
                f'{activation_dynamics} as it must be of type '
                f'`ActivationBase`, not {type(activation_dynamics)}.'
            )
            raise TypeError(msg)
        self._activation_dynamics = activation_dynamics
        self._child_objects = (self._activation_dynamics, )

        # Constants
        if tendon_slack_length is not None:
            self._l_T_slack = tendon_slack_length
        else:
            self._l_T_slack = Symbol(f'l_T_slack_{self.name}')
        if peak_isometric_force is not None:
            self._F_M_max = peak_isometric_force
        else:
            self._F_M_max = Symbol(f'F_M_max_{self.name}')
        if optimal_fiber_length is not None:
            self._l_M_opt = optimal_fiber_length
        else:
            self._l_M_opt = Symbol(f'l_M_opt_{self.name}')
        if maximal_fiber_velocity is not None:
            self._v_M_max = maximal_fiber_velocity
        else:
            self._v_M_max = Symbol(f'v_M_max_{self.name}')
        if optimal_pennation_angle is not None:
            self._alpha_opt = optimal_pennation_angle
        else:
            self._alpha_opt = Symbol(f'alpha_opt_{self.name}')
        if fiber_damping_coefficient is not None:
            self._beta = fiber_damping_coefficient
        else:
            self._beta = Symbol(f'beta_{self.name}')

        # Musculotendon dynamics
        if musculotendon_dynamics == MusculotendonFormulation.RIGID_TENDON:
            self._rigid_tendon_musculotendon_dynamics()
        elif musculotendon_dynamics == MusculotendonFormulation.FIBER_LENGTH_EXPLICIT:
            self._fiber_length_explicit_musculotendon_dynamics()
        elif musculotendon_dynamics == MusculotendonFormulation.TENDON_FORCE_EXPLICIT:
            self._tendon_force_explicit_musculotendon_dynamics()
        elif musculotendon_dynamics == MusculotendonFormulation.FIBER_LENGTH_IMPLICIT:
            self._fiber_length_implicit_musculotendon_dynamics()
        elif musculotendon_dynamics == MusculotendonFormulation.TENDON_FORCE_IMPLICIT:
            self._tendon_force_implicit_musculotendon_dynamics()
        else:
            msg = (
                f'Musculotendon dynamics {repr(musculotendon_dynamics)} '
                f'passed to `musculotendon_dynamics` was of type '
                f'{type(musculotendon_dynamics)}, must be '
                f'{MusculotendonFormulation}.'
            )
            raise TypeError(msg)
        self._musculotendon_dynamics = musculotendon_dynamics

        # Must override the placeholder value in `self._force` now that the
        # actual force has been calculated by
        # `self._<MUSCULOTENDON FORMULATION>_musculotendon_dynamics`.
        # Note that `self._force` assumes forces are expansile, musculotendon
        # forces are contractile hence the minus sign preceeding `self._F_T`
        # (the tendon force).
        self._force = -self._F_T

    @classmethod
    def with_defaults(
        cls,
        name,
        pathway,
        activation_dynamics,
        *,
        musculotendon_dynamics=_DEFAULT_MUSCULOTENDON_FORMULATION,
        tendon_slack_length=None,
        peak_isometric_force=None,
        optimal_fiber_length=None,
    ):
        v_M_max = Integer(10)
        alpha_opt = Integer(0)
        beta = Rational(1, 10)
        return cls(
            name,
            pathway,
            activation_dynamics=activation_dynamics,
            musculotendon_dynamics=musculotendon_dynamics,
            tendon_slack_length=tendon_slack_length,
            peak_isometric_force=peak_isometric_force,
            optimal_fiber_length=optimal_fiber_length,
            maximal_fiber_velocity=v_M_max,
            optimal_pennation_angle=alpha_opt,
            fiber_damping_coefficient=beta,
        )

    @property
    def tendon_slack_length(self):
        return self._l_T_slack

    @property
    def l_T_slack(self):
        return self._l_T_slack

    @property
    def peak_isometric_force(self):
        return self._F_M_max

    @property
    def F_M_max(self):
        return self._F_M_max

    @property
    def optimal_fiber_length(self):
        return self._l_M_opt

    @property
    def l_M_opt(self):
        return self._l_M_opt

    @property
    def maximal_fiber_velocity(self):
        return self._v_M_max

    @property
    def v_M_max(self):
        return self._v_M_max

    @property
    def optimal_pennation_angle(self):
        return self._alpha_opt

    @property
    def alpha_opt(self):
        return self._alpha_opt

    @property
    def fiber_damping_coefficient(self):
        return self._beta

    @property
    def beta(self):
        return self._beta

    @property
    def activation_dynamics(self):
        """Activation dynamics model governing this musculotendon's activation.

        Explanation
        ===========

        Returns the instance of a subclass of ``ActivationBase`` that governs
        the relationship between excitation and activation that is used to
        represent the activation dynamics of this musculotendon.

        """
        return self._activation_dynamics

    @property
    def excitation(self):
        """Dynamic symbol representing excitation.

        Explanation
        ===========

        The alias `e` can also be used to access the same attribute.

        """
        return self._activation_dynamics._e

    @property
    def e(self):
        """Dynamic symbol representing excitation.

        Explanation
        ===========

        The alias `excitation` can also be used to access the same attribute.

        """
        return self._activation_dynamics._e

    @property
    def activation(self):
        """Dynamic symbol representing activation.

        Explanation
        ===========

        The alias `a` can also be used to access the same attribute.

        """
        return self._activation_dynamics._a

    @property
    def a(self):
        """Dynamic symbol representing activation.

        Explanation
        ===========

        The alias `activation` can also be used to access the same attribute.

        """
        return self._activation_dynamics._a

    @property
    def musculotendon_dynamics(self):
        return self._musculotendon_dynamics

    def _rigid_tendon_musculotendon_dynamics(self):
        """Rigid tendon musculotendon."""
        self._l_MT = self.pathway.length
        self._v_MT = self.pathway.extension_velocity
        self._l_T = self._l_T_slack
        self._l_T_tilde = Integer(1)
        self._l_M = sqrt((self._l_MT - self._l_T)**2 + (self._l_M_opt*sin(self._alpha_opt))**2)
        self._l_M_tilde = self._l_M/self._l_M_opt
        self._v_M = self._v_MT*(self._l_MT - self._l_T_slack)/self._l_M
        self._v_M_tilde = self._v_M/self._v_M_max
        self._fl_T = TendonForceLengthDeGroote2016.with_defaults(self._l_T_tilde)
        self._fl_M_pas = FiberForceLengthPassiveDeGroote2016.with_defaults(self._l_M_tilde)
        self._fl_M_act = FiberForceLengthActiveDeGroote2016.with_defaults(self._l_M_tilde)
        self._fv_M = FiberForceVelocityDeGroote2016.with_defaults(self._v_M_tilde)
        self._F_M_tilde = self.a*self._fl_M_act*self._fv_M + self._fl_M_pas + self._beta*self._v_M_tilde
        self._F_T_tilde = self._F_M_tilde
        self._F_M = self._F_M_tilde*self._F_M_max
        self._F_T = self._F_M

        # Containers
        self._state_vars = zeros(0, 1)
        self._input_vars = zeros(0, 1)
        self._state_eqns = zeros(0, 1)

    def _fiber_length_explicit_musculotendon_dynamics(self):
        """Elastic tendon musculotendon using `l_M_tilde` as a state."""
        self._l_M_tilde = dynamicsymbols(f'l_M_tilde_{self.name}')
        self._l_MT = self.pathway.length
        self._v_MT = self.pathway.extension_velocity
        self._l_M = self._l_M_tilde*self._l_M_opt
        self._l_T = self._l_MT - sqrt(self._l_M**2 - (self._l_M_opt*sin(self._alpha_opt))**2)
        self._l_T_tilde = self._l_T/self._l_T_slack
        self._cos_alpha = (self._l_MT - self._l_T)/self._l_M
        self._fl_T = TendonForceLengthDeGroote2016.with_defaults(self._l_T_tilde)
        self._fl_M_pas = FiberForceLengthPassiveDeGroote2016.with_defaults(self._l_M_tilde)
        self._fl_M_act = FiberForceLengthActiveDeGroote2016.with_defaults(self._l_M_tilde)
        self._F_T_tilde = self._fl_T
        self._F_T = self._F_T_tilde*self._F_M_max
        self._F_M = self._F_T/self._cos_alpha
        self._F_M_tilde = self._F_M/self._F_M_max
        self._fv_M = (self._F_M_tilde - self._fl_M_pas)/(self.a*self._fl_M_act)
        self._v_M_tilde = FiberForceVelocityDeGroote2016.with_defaults(self._fv_M)
        self._dl_M_tilde_dt = (self._v_M_max/self._l_M_opt)*self._v_M_tilde

        self._state_vars = Matrix([self._l_M_tilde])
        self._input_vars = zeros(0, 1)
        self._state_eqns = Matrix([self._dl_M_tilde_dt])

    def _tendon_force_explicit_musculotendon_dynamics(self):
        """Elastic tendon musculotendon using `F_T_tilde` as a state."""
        raise NotImplementedError

    def _fiber_length_implicit_musculotendon_dynamics(self):
        raise NotImplementedError

    def _tendon_force_implicit_musculotendon_dynamics(self):
        raise NotImplementedError

    @property
    def state_vars(self):
        """Ordered column matrix of functions of time that represent the state
        variables.

        Explanation
        ===========

        As zeroth-order activation dynamics simply maps excitation to
        activation, this class has no associated state variables and so this
        property return an empty column ``Matrix`` with shape (0, 1).

        The alias ``x`` can also be used to access the same attribute.

        """
        state_vars = [self._state_vars]
        for child in self._child_objects:
            state_vars.append(child.state_vars)
        return Matrix.vstack(*state_vars)

    @property
    def x(self):
        """Ordered column matrix of functions of time that represent the state
        variables.

        Explanation
        ===========

        As zeroth-order activation dynamics simply maps excitation to
        activation, this class has no associated state variables and so this
        property return an empty column ``Matrix`` with shape (0, 1).

        The alias ``state_vars`` can also be used to access the same attribute.

        """
        state_vars = [self._state_vars]
        for child in self._child_objects:
            state_vars.append(child.state_vars)
        return Matrix.vstack(*state_vars)

    @property
    def input_vars(self):
        """Ordered column matrix of functions of time that represent the input
        variables.

        Explanation
        ===========

        Excitation is the only input in zeroth-order activation dynamics and so
        this property returns a column ``Matrix`` with one entry, ``e``, and
        shape (1, 1).

        The alias ``r`` can also be used to access the same attribute.

        """
        input_vars = [self._input_vars]
        for child in self._child_objects:
            input_vars.append(child.input_vars)
        return Matrix.vstack(*input_vars)

    @property
    def r(self):
        """Ordered column matrix of functions of time that represent the input
        variables.

        Explanation
        ===========

        Excitation is the only input in zeroth-order activation dynamics and so
        this property returns a column ``Matrix`` with one entry, ``e``, and
        shape (1, 1).

        The alias ``input_vars`` can also be used to access the same attribute.

        """
        input_vars = [self._input_vars]
        for child in self._child_objects:
            input_vars.append(child.input_vars)
        return Matrix.vstack(*input_vars)

    @property
    def constants(self):
        """Ordered column matrix of non-time varying symbols present in ``M``
        and ``F``.

        Explanation
        ===========

        As zeroth-order activation dynamics simply maps excitation to
        activation, this class has no associated constants and so this property
        return an empty column ``Matrix`` with shape (0, 1).

        The alias ``p`` can also be used to access the same attribute.

        """
        constants = [Matrix([
            self._l_T_slack,
            self._F_M_max,
            self._l_M_opt,
            self._v_M_max,
            self._alpha_opt,
            self._beta,
        ])]
        for child in self._child_objects:
            constants.append(child.constants)
        return Matrix.vstack(*constants)

    @property
    def p(self):
        """Ordered column matrix of non-time varying symbols present in ``M``
        and ``F``.

        Explanation
        ===========

        As zeroth-order activation dynamics simply maps excitation to
        activation, this class has no associated constants and so this property
        return an empty column ``Matrix`` with shape (0, 1).

        The alias ``constants`` can also be used to access the same attribute.

        """
        constants = [Matrix([
            self._l_T_slack,
            self._F_M_max,
            self._l_M_opt,
            self._v_M_max,
            self._alpha_opt,
            self._beta,
        ])]
        for child in self._child_objects:
            constants.append(child.constants)
        return Matrix.vstack(*constants)

    @property
    def M(self):
        """Ordered square matrix of coefficients on the LHS of ``M x' = F``.

        Explanation
        ===========

        The square matrix that forms part of the LHS of the linear system of
        ordinary differential equations governing the activation dynamics:

        ``M(x, r, t, p) x' = F(x, r, t, p)``.

        As zeroth-order activation dynamics have no state variables, this
        linear system has dimension 0 and therefore ``M`` is an empty square
        ``Matrix`` with shape (0, 0).

        """
        M = [eye(len(self._state_vars))]
        for child in self._child_objects:
            M.append(child.M)
        return diag(*M)

    @property
    def F(self):
        """Ordered column matrix of equations on the RHS of ``M x' = F``.

        Explanation
        ===========

        The column matrix that forms the RHS of the linear system of ordinary
        differential equations governing the activation dynamics:

        ``M(x, r, t, p) x' = F(x, r, t, p)``.

        As zeroth-order activation dynamics have no state variables, this
        linear system has dimension 0 and therefore ``F`` is an empty column
        ``Matrix`` with shape (0, 1).

        """
        F = [self._state_eqns]
        for child in self._child_objects:
            F.append(child.F)
        return Matrix.vstack(*F)

    def rhs(self):
        """Ordered column matrix of equations for the solution of ``M x' = F``.

        Explanation
        ===========

        The solution to the linear system of ordinary differential equations
        governing the activation dynamics:

        ``M(x, r, t, p) x' = F(x, r, t, p)``.

        As zeroth-order activation dynamics have no state variables, this
        linear has dimension 0 and therefore this method returns an empty
        column ``Matrix`` with shape (0, 1).

        """
        return self.M.solve(self.F)

    def __repr__(self):
        return (
            f'{self.__class__.__name__}({self.name!r}, '
            f'pathway={self.pathway!r}, '
            f'activation_dynamics={self.activation_dynamics!r}, '
            f'musculotendon_dynamics={self.musculotendon_dynamics}, '
            f'tendon_slack_length={self._l_T_slack!r}, '
            f'peak_isometric_force={self._F_M_max!r}, '
            f'optimal_fiber_length={self._l_M_opt!r}, '
            f'maximal_fiber_velocity={self._v_M_max!r}, '
            f'optimal_pennation_angle={self._alpha_opt!r}, '
            f'fiber_damping_coefficient={self._beta!r})'
        )
