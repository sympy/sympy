"""Implementations of musculotendon models.

Musculotendon models are a critical component of biomechanical models, one that
differentiates them from pure multibody systems. Musculotendon models produce a
force dependent on their level of activation, their length, and their
extension velocity. Length- and extension velocity-dependent force production
are governed by force-length and force-velocity characteristics.
These are normalized functions that are dependent on the musculotendon's state
and are specific to a given musculotendon model.

"""

from abc import abstractmethod
from enum import IntEnum, unique

from sympy.core.numbers import Integer, Rational
from sympy.core.symbol import Symbol, symbols
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.trigonometric import cos, sin
from sympy.matrices.dense import MutableDenseMatrix as Matrix, diag, eye, zeros
from sympy.physics._biomechanics.activation import ActivationBase
from sympy.physics._biomechanics.curve import (
    CharacteristicCurveCollection,
    FiberForceLengthActiveDeGroote2016,
    FiberForceLengthPassiveDeGroote2016,
    FiberForceLengthPassiveInverseDeGroote2016,
    FiberForceVelocityDeGroote2016,
    FiberForceVelocityInverseDeGroote2016,
    TendonForceLengthDeGroote2016,
    TendonForceLengthInverseDeGroote2016,
)
from sympy.physics._biomechanics._mixin import _NamedMixin
from sympy.physics.mechanics.actuator import ForceActuator
from sympy.physics.vector.functions import dynamicsymbols


__all__ = [
    'MusculotendonBase',
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
        A rigid tendon model.
    FIBER_LENGTH_EXPLICIT : 1
        An explicit elastic tendon model with the muscle fiber length (l_M) as
        the state variable.
    TENDON_FORCE_EXPLICIT : 2
        An explicit elastic tendon model with the tendon force (F_T) as the
        state variable.
    FIBER_LENGTH_IMPLICIT : 3
        An implicit elastic tendon model with the muscle fiber length (l_M) as
        the state variable and the muscle fiber velocity as an additional input
        variable.
    TENDON_FORCE_IMPLICIT : 4
        An implicit elastic tendon model with the tendon force (F_T) as the
        state variable as the muscle fiber velocity as an additional input
        variable.

    """

    RIGID_TENDON = 0
    FIBER_LENGTH_EXPLICIT = 1
    TENDON_FORCE_EXPLICIT = 2
    FIBER_LENGTH_IMPLICIT = 3
    TENDON_FORCE_IMPLICIT = 4


_DEFAULT_MUSCULOTENDON_FORMULATION = MusculotendonFormulation.RIGID_TENDON


class MusculotendonBase(ForceActuator, _NamedMixin):
    r"""Abstract base class for all musculotendon classes to inherit from.

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
        The activation dynamics that will be modeled within the musculotendon.
        This must be an instance of a concrete subclass of ``ActivationBase``,
        e.g. ``FirstOrderActivationDeGroote2016``.
    musculotendon_dynamics : MusculotendonFormulation | int
        The formulation of musculotendon dynamics that should be used
        internally, i.e. rigid or elastic tendon model, the choice of
        musculotendon state etc. This must be a member of the integer
        enumeration ``MusculotendonFormulation`` or an integer that can be cast
        to a member. The default is ``MusculotendonFormulation.RIGID_TENDON``
        (or ``0``), which corresponds to a rigid tendon formulation.
    tendon_slack_length : Expr | None
        The length of the tendon when the musculotendon is in its unloaded
        state. In a rigid tendon model the tendon length is the tendon slack
        length. In all musculotendon models, tendon slack length is used to
        normalize tendon length to give
        :math:`\tilde{l}^T = \frac{l^T}{l^T_{slack}}`.
    peak_isometric_force : Expr | None
        The maximum force that the muscle fiber can produce when it is
        undergoing an isometric contraction (no lengthening velocity). In all
        musculotendon models, peak isometric force is used to normalized tendon
        and muscle fiber force to give
        :math:`\tilde{F}^T = \frac{F^T}{F^M_{max}}`.
    optimal_fiber_length : Expr | None
        The muscle fiber length at which the muscle fibers produce no passive
        force and their maximum active force. In all musculotendon models,
        optimal fiber length is used to normalize muscle fiber length to give
        :math:`\tilde{l}^M = \frac{l^M}{l^M_{opt}}`.
    maximal_fiber_velocity : Expr | None
        The fiber velocity at which, during muscle fiber shortening, the muscle
        fibers are unable to produce any active force. In all musculotendon
        models, maximal fiber velocity is used to normalize muscle fiber
        extension velocity to give :math:`\tilde{v}^M = \frac{v^M}{v^M_{max}}`.
    optimal_pennation_angle : Expr | None
        The pennation angle when muscle fiber length equals the optimal fiber
        length.
    fiber_damping_coefficient : Expr | None
        The coefficient of damping to be used in the damping element in the
        muscle fiber model.
    with_defaults : bool
        Whether ``with_defaults`` alternate constructors should be used when
        automatically constructing child classes. Default is ``False``.

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
        with_defaults=False,
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
        self._with_defaults = with_defaults
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
            with_defaults=True,
        )

    @abstractmethod
    def curves(cls):
        """Return a ``CharacteristicCurveCollection`` of the curves related to
        the specific model."""
        pass

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

        The alias ``e`` can also be used to access the same attribute.

        """
        return self._activation_dynamics._e

    @property
    def e(self):
        """Dynamic symbol representing excitation.

        Explanation
        ===========

        The alias ``excitation`` can also be used to access the same attribute.

        """
        return self._activation_dynamics._e

    @property
    def activation(self):
        """Dynamic symbol representing activation.

        Explanation
        ===========

        The alias ``a`` can also be used to access the same attribute.

        """
        return self._activation_dynamics._a

    @property
    def a(self):
        """Dynamic symbol representing activation.

        Explanation
        ===========

        The alias ``activation`` can also be used to access the same attribute.

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
        if self._with_defaults:
            self._fl_T = self.curves.tendon_force_length.with_defaults(self._l_T_tilde)
            self._fl_M_pas = self.curves.fiber_force_length_passive.with_defaults(self._l_M_tilde)
            self._fl_M_act = self.curves.fiber_force_length_active.with_defaults(self._l_M_tilde)
            self._fv_M = self.curves.fiber_force_velocity.with_defaults(self._v_M_tilde)
        else:
            fl_T_constants = symbols('c_0:4_fl_T')
            self._fl_T = self.curves.tendon_force_length(self._l_T_tilde, *fl_T_constants)
            fl_M_pas_constants = symbols('c_0:2_fl_M_pas')
            self._fl_M_pas = self.curves.fiber_force_length_passive(self._l_M_tilde, *fl_M_pas_constants)
            fl_M_act_constants = symbols('c_0:12_fl_M_act')
            self._fl_M_act = self.curves.fiber_force_length_active(self._l_M_tilde, *fl_M_act_constants)
            fv_M_constants = symbols('c_0:4_fv_M')
            self._fv_M = self.curves.fiber_force_velocity(self._v_M_tilde, *fv_M_constants)
        self._F_M_tilde = self.a*self._fl_M_act*self._fv_M + self._fl_M_pas + self._beta*self._v_M_tilde
        self._F_T_tilde = self._F_M_tilde
        self._F_M = self._F_M_tilde*self._F_M_max
        self._cos_alpha = cos(self._alpha_opt)
        self._F_T = self._F_M*self._cos_alpha

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
        if self._with_defaults:
            self._fl_T = self.curves.tendon_force_length.with_defaults(self._l_T_tilde)
            self._fl_M_pas = self.curves.fiber_force_length_passive.with_defaults(self._l_M_tilde)
            self._fl_M_act = self.curves.fiber_force_length_active.with_defaults(self._l_M_tilde)
        else:
            fl_T_constants = symbols('c_0:4_fl_T')
            self._fl_T = self.curves.tendon_force_length(self._l_T_tilde, *fl_T_constants)
            fl_M_pas_constants = symbols('c_0:2_fl_M_pas')
            self._fl_M_pas = self.curves.fiber_force_length_passive(self._l_M_tilde, *fl_M_pas_constants)
            fl_M_act_constants = symbols('c_0:12_fl_M_act')
            self._fl_M_act = self.curves.fiber_force_length_active(self._l_M_tilde, *fl_M_act_constants)
        self._F_T_tilde = self._fl_T
        self._F_T = self._F_T_tilde*self._F_M_max
        self._F_M = self._F_T/self._cos_alpha
        self._F_M_tilde = self._F_M/self._F_M_max
        self._fv_M = (self._F_M_tilde - self._fl_M_pas)/(self.a*self._fl_M_act)
        if self._with_defaults:
            self._v_M_tilde = self.curves.fiber_force_velocity_inverse.with_defaults(self._fv_M)
        else:
            fv_M_constants = symbols('c_0:4_fv_M')
            self._v_M_tilde = self.curves.fiber_force_velocity_inverse(self._fv_M, *fv_M_constants)
        self._dl_M_tilde_dt = (self._v_M_max/self._l_M_opt)*self._v_M_tilde

        self._state_vars = Matrix([self._l_M_tilde])
        self._input_vars = zeros(0, 1)
        self._state_eqns = Matrix([self._dl_M_tilde_dt])

    def _tendon_force_explicit_musculotendon_dynamics(self):
        """Elastic tendon musculotendon using `F_T_tilde` as a state."""
        raise NotImplementedError
        self._F_T_tilde = dynamicsymbols(f'F_T_tilde_{self.name}')
        self._l_MT = self.pathway.length
        self._v_MT = self.pathway.extension_velocity
        self._fl_T = self._F_T_tilde
        if self._with_defaults:
            self._fl_T_inv = self.curves.tendon_force_length_inverse.with_defaults(self._fl_T)
        else:
            fl_T_constants = symbols('c_0:4_fl_T')
            self._fl_T_inv = self.curves.tendon_force_length_inverse(self._fl_T, *fl_T_constants)
        self._l_T_tilde = self._fl_T_inv
        self._l_T = self._l_T_tilde*self._l_T_slack
        self._l_M = sqrt((self._l_MT - self._l_T)**2 + (self._l_M_opt*sin(self._alpha_opt))**2)
        self._l_M_tilde = self._l_M/self._l_M_opt
        if self._with_defaults:
            self._fl_M_pas = self.curves.fiber_force_length_passive.with_defaults(self._l_M_tilde)
            self._fl_M_act = self.curves.fiber_force_length_active.with_defaults(self._l_M_tilde)
        else:
            fl_M_pas_constants = symbols('c_0:2_fl_M_pas')
            self._fl_M_pas = self.curves.fiber_force_length_passive(self._l_M_tilde, *fl_M_pas_constants)
            fl_M_act_constants = symbols('c_0:12_fl_M_act')
            self._fl_M_act = self.curves.fiber_force_length_active(self._l_M_tilde, *fl_M_act_constants)
        self._cos_alpha = (self._l_MT - self._l_T)/self._l_M
        self._F_T = self._F_T_tilde*self._F_M_max
        self._F_M = self._F_T/self._cos_alpha
        self._F_M_tilde = self._F_M/self._F_M_max
        self._fv_M = (self._F_M_tilde - self._fl_M_pas)/(self.a*self._fl_M_act)
        if self._with_defaults:
            self._fv_M_inv = self.curves.fiber_force_velocity_inverse.with_defaults(self._fv_M)
        else:
            fv_M_constants = symbols('c_0:4_fv_M')
            self._fv_M_inv = self.curves.fiber_force_velocity_inverse(self._fv_M, *fv_M_constants)
        self._v_M_tilde = self._fv_M_inv
        self._v_M = self._v_M_tilde*self._v_M_max
        self._v_T = self._v_MT - (self._v_M/self._cos_alpha)
        self._v_T_tilde = self._v_T/self._l_T_slack
        if self._with_defaults:
            self._fl_T = self.curves.tendon_force_length.with_defaults(self._l_T_tilde)
        else:
            fl_T_constants = symbols('c_0:4_fl_T')
            self._fl_T = self.curves.tendon_force_length(self._l_T_tilde, *fl_T_constants)
        self._dF_T_tilde_dt = self._fl_T.diff(dynamicsymbols._t).subs({self._l_T_tilde.diff(dynamicsymbols._t): self._v_T_tilde})

        self._state_vars = Matrix([self._F_T_tilde])
        self._input_vars = zeros(0, 1)
        self._state_eqns = Matrix([self._dF_T_tilde_dt])

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
        """Returns a string representation to reinstantiate the model."""
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

    def __str__(self):
        """Returns a string representation of the expression for musculotendon
        force."""
        return str(self.force)


class MusculotendonDeGroote2016(MusculotendonBase):
    """Musculotendon model using the curves of De Groote et al., 2016 [1].

    References
    ==========

    .. [1] De Groote, F., Kinney, A. L., Rao, A. V., & Fregly, B. J., Evaluation
           of direct collocation optimal control problem formulations for
           solving the muscle redundancy problem, Annals of biomedical
           engineering, 44(10), (2016) pp. 2922-2936

    """

    curves = CharacteristicCurveCollection(
        tendon_force_length=TendonForceLengthDeGroote2016,
        tendon_force_length_inverse=TendonForceLengthInverseDeGroote2016,
        fiber_force_length_passive=FiberForceLengthPassiveDeGroote2016,
        fiber_force_length_passive_inverse=FiberForceLengthPassiveInverseDeGroote2016,
        fiber_force_length_active=FiberForceLengthActiveDeGroote2016,
        fiber_force_velocity=FiberForceVelocityDeGroote2016,
        fiber_force_velocity_inverse=FiberForceVelocityInverseDeGroote2016,
    )
