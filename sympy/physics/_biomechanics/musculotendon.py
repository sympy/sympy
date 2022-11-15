"""Implementations of, and factories for, musculotendon models.

Musculotendon models are probably the most important component of biomechanical
models that differentiate them from pure multibody systems. Musculotendon models
produce a force which is dependent on their level of activate, their length, and
their shortening velocity. Length- and shortening velocity-dependent force
production are governed by force-length and force-velocity characteristics.
These are normalized functions that are dependent on the musculotendon's state
and are specific to a given musculotendon model.

"""


import abc
from typing import Optional

import sympy as sm
import sympy.physics.mechanics as me
from sympy.physics._biomechanics.mixin import _NamedMixin


__all__ = [
    'Brockie2021Musculotendon',
    'DeGroote2016Musculotendon',
    'Millard2013Musculotendon',
    'Musculotendon',
]


class MusculotendonBase(abc.ABC, _NamedMixin):
    """Abstract base class for all musculotendon classes to inherit from.

    Explanation
    ===========

    A musculotendon generates a contractile force based on its activation,
    length, and shortening velocity. This abstract base class is to be inherited
    by all musculotendon subclasses that implement different characteristic
    musculotendon curves. Characteristic musculotendon curves are required for
    the tendon force-length, passive fiber force-length, active fiber force-
    length, and fiber force-velocity relationships.

    Attributes
    ==========
    name : str
        The name identifier associated with the musculotendon. This name is used
        as a suffix when automatically generated symbols are instantiated. It
        must be a string of nonzero length.
    origin : `sympy.physics.mechanics.Point`
        The point from which the musculotendon originates.
    insertion : `sympy.physics.mechanics.Point`
        The point to which the musculotendon inserts.
    optimal_fiber_length : Optional[float]
        The value by which the musculotendon's fiber force-length relationships
        are normalized. It corresponds to the fiber length at which the
        musculotendon can produce its maximal force. This value maps to the
        symbol attribute `l_M_opt`.
    maximal_fiber_velocity : float, default is 10.0
        The value by which the musculotendon's fiber force-velcoity relationship
        is normalized. It corresponds to the fiber shortening velocity at which
        the musculotendon is unable to produce any force whatsoever. This value
        maps to the symbol attribute `v_M_max`.
    peak_isometric_force : Optional[float]
        The force the musculotendon can produce when the fiber length is equal
        to the optimal fiber length, the shortening velocity is zero, and at
        full activation. "Isometric" refers to the fact that the musculotendon
        is undergoing an isometric contraction (i.e. maintaining length) at this
        magnitude of force. This value maps to the symbol attribute `F_M_max`.
    tendon_slack_length : Optional[float]
        The value by which the musculotendon's tendon force-length relationship
        is normalized. It corresponds to the tendon length at which the tendon
        is unable to produce any tensile force due to becoming slack. For a
        rigid tendon model, this value is exactly equal to the tendon length.
        This value maps to the symbol attribute `l_T_slack`.
    optimal_pennation_angle : float, default is 0.0
        The value of pennation angle when the musculotendon's fiber has length
        equal to the optimal fiber length. For unpennated fibers, this value is
        zero. This value maps to the symbol attribute `alpha_opt`.
    fiber_damping_coefficient : float, default is 0.1
        The value of damping in the musculotendon's fiber. While not strictly
        physical, this parameter is important for producing numerically stable
        musculotendon dynamics. For undamped fibers, this value is zero. This
        value maps to the symbol attribute `beta`.
        @property
    l_M_opt : `sympy.Symbol`
        Accessor for the optimal fiber length symbol.
    v_M_max : `sympy.Symbol`
        Accessor for the maximal fiber velocity symbol.
    F_M_max : `sympy.Symbol`
        Accessor for the peak isometric force symbol.
    l_T_slack : `sympy.Symbol`
        Accessor for the tendon slack length symbol.
    alpha_opt : `sympy.Symbol`
        Accessor for the optimal pennation angle symbol.
    beta : `sympy.Symbol`
        Accessor for the fiber damping coefficient symbol.
    l_MT : dynamic symbol
        Accessor for the musculotendon length dynamic symbol.
    v_MT : dynamic symbol
        Accessor for the musculotendon shortening velocity dynamic symbol.
    l_T : dynamic symbol
        Accessor for the tendon length dynamic symbol.
    v_T : dynamic symbol
        Accessor for the tendon shortening velocity dynamic symbol.
    l_M : dynamic symbol
        Accessor for the fiber length dynamic symbol.
    v_M : dynamic symbol
        Accessor for the fiber shortening velocity dynamic symbol.
    l_T_tilde : dynamic symbol
        Accessor for the normalized tendon length dynamic symbol.
    v_T_tilde : dynamic symbol
        Accessor for the normalized tendon shortening velocity dynamic symbol.
    l_M_tilde : dynamic symbol
        Accessor for the normalized fiber length dynamic symbol.
    v_M_tilde : dynamic symbol
        Accessor for the normalized fiber shortening velocity dynamic symbol.
    F_T : dynamic symbol
        Accessor for the tendon force dynamic symbol.
    F_M : dynamic symbol
        Accessor for the fiber force dynamic symbol.
    F_T_tilde : dynamic symbol
        Accessor for the normalized tendon force dynamic symbol.
    F_M_tilde : dynamic symbol
        Accessor for the normalized fiber force dynamic symbol.
    fl_T : dynamic symbol
        Accessor for the tendon force-length dynamic symbol.
    fl_M_pas : dynamic symbol
        Accessor for the passive fiber force-length dynamic symbol.
    fl_M_act : dynamic symbol
        Accessor for the active fiber force-length dynamic symbol.
    fv_M : dynamic symbol
        Accessor for the fiber force-velocity dynamic symbol.
    cos_alpha : dynamic symbol
        Accessor for the cosine of the pennation angle dynamic symbol.
    dF_T_tilde_dt : dynamic symbol
        Accessor for the rate of change in tendon force dynamic symbol.
    dl_M_tilde_dt : dynamic symbol
        Accessor for the rate of change in fiber length dynamic symbol.
    F_orig_x : dynamic symbol
        Accessor for the x-axis origin force dynamic symbol.
    F_orig_y : dynamic symbol
        Accessor for the y-axis origin force dynamic symbol.
    F_orig_z : dynamic symbol
        Accessor for the z-axis origin force dynamic symbol.
    F_insr_x : dynamic symbol
        Accessor for the x-axis insertion force dynamic symbol.
    F_insr_y : dynamic symbol
        Accessor for the y-axis insertion force dynamic symbol.
    F_insr_z : dynamic symbol
        Accessor for the z-axis insertion force dynamic symbol.

    """

    def __init__(
        self,
        name: str,
        *,
        origin: me.Point,
        insertion: me.Point,
        optimal_fiber_length: Optional[float] = None,
        maximal_fiber_velocity: float = 10.0,
        peak_isometric_force: Optional[float] = None,
        tendon_slack_length: Optional[float] = None,
        optimal_pennation_angle: float = 0.0,
        fiber_damping_coefficient: float = 0.1,
    ) -> None:
        """Initializer for a musculotendon's instance attributes.

        Parameters
        ==========
        name : str
            The name identifier associated with the musculotendon. This name is
            used as a suffix when automatically generated symbols are
            instantiated. It must be a string of nonzero length.
        origin : `sympy.physics.mechanics.Point`
            The point from which the musculotendon originates.
        insertion : `sympy.physics.mechanics.Point`
            The point to which the musculotendon inserts.
        optimal_fiber_length : Optional[float]
            The value by which the musculotendon's fiber force-length
            relationships are normalized. It corresponds to the fiber length at
            which the musculotendon can produce its maximal force. This value
            maps to the symbol attribute `l_M_opt`.
        maximal_fiber_velocity : float, default is 10.0
            The value by which the musculotendon's fiber force-velcoity
            relationship is normalized. It corresponds to the fiber shortening
            velocity at which the musculotendon is unable to produce any force
            whatsoever. This value maps to the symbol attribute `v_M_max`.
        peak_isometric_force : Optional[float]
            The force the musculotendon can produce when the fiber length is
            equal to the optimal fiber length, the shortening velocity is zero,
            and at full activation. "Isometric" refers to the fact that the
            musculotendon is undergoing an isometric contraction (i.e.
            maintaining length) at this magnitude of force. This value maps to
            the symbol attribute `F_M_max`.
        tendon_slack_length : Optional[float]
            The value by which the musculotendon's tendon force-length
            relationship is normalized. It corresponds to the tendon length at
            which the tendon is unable to produce any tensile force due to
            becoming slack. For a rigid tendon model, this value is exactly
            equal to the tendon length. This value maps to the symbol attribute
            `l_T_slack`.
        optimal_pennation_angle : float, default is 0.0
            The value of pennation angle when the musculotendon's fiber has
            length equal to the optimal fiber length. For unpennated fibers,
            this value is zero. This value maps to the symbol attribute
            `alpha_opt`.
        fiber_damping_coefficient : float, default is 0.1
            The value of damping in the musculotendon's fiber. While not
            strictly physical, this parameter is important for producing
            numerically stable musculotendon dynamics. For undamped fibers, this
            value is zero. This value maps to the symbol attribute `beta`.

        """
        self.name = name
        self.origin = origin
        self.insertion = insertion

        # Constants
        self.optimal_fiber_length = optimal_fiber_length
        self.maximal_fiber_velocity = maximal_fiber_velocity
        self.peak_isometric_force = peak_isometric_force
        self.tendon_slack_length = tendon_slack_length
        self.optimal_pennation_angle = optimal_pennation_angle
        self.fiber_damping_coefficient = fiber_damping_coefficient

        # Symbols
        self._l_MT = me.dynamicsymbols(f"l_MT_{self.name}")
        self._v_MT = me.dynamicsymbols(f"v_MT_{self.name}")
        self._l_T = me.dynamicsymbols(f"l_T_{self.name}")
        self._v_T = me.dynamicsymbols(f"v_T_{self.name}")
        self._l_M = me.dynamicsymbols(f"l_M_{self.name}")
        self._v_M = me.dynamicsymbols(f"v_M_{self.name}")

        self._l_T_tilde = me.dynamicsymbols(f"l_T_tilde_{self.name}")
        self._v_T_tilde = me.dynamicsymbols(f"v_T_tilde_{self.name}")
        self._l_M_tilde = me.dynamicsymbols(f"l_M_tilde_{self.name}")
        self._v_M_tilde = me.dynamicsymbols(f"v_M_tilde_{self.name}")

        self._F_T = me.dynamicsymbols(f"F_T_{self.name}")
        self._F_M = me.dynamicsymbols(f"F_M_{self.name}")

        self._F_T_tilde = me.dynamicsymbols(f"F_T_tilde_{self.name}")
        self._F_M_tilde = me.dynamicsymbols(f"F_M_tilde_{self.name}")

        self._fl_T = me.dynamicsymbols(f"fl_T_{self.name}")
        self._fl_M_pas = me.dynamicsymbols(f"fl_M_pas_{self.name}")
        self._fl_M_act = me.dynamicsymbols(f"fl_M_act_{self.name}")
        self._fv_M = me.dynamicsymbols(f"fv_M_{self.name}")

        self._cos_alpha = me.dynamicsymbols(f"cos_alpha_{self.name}")
        self._dF_T_tilde_dt = me.dynamicsymbols(f"dF_T_tilde_dt_{self.name}")
        self._dl_M_tilde_dt = me.dynamicsymbols(f"dl_M_tilde_dt_{self.name}")

        self._F_orig_x = me.dynamicsymbols(f"F_orig_x_{self.name}")
        self._F_orig_y = me.dynamicsymbols(f"F_orig_y_{self.name}")
        self._F_orig_z = me.dynamicsymbols(f"F_orig_z_{self.name}")
        self._F_insr_x = me.dynamicsymbols(f"F_insr_x_{self.name}")
        self._F_insr_y = me.dynamicsymbols(f"F_insr_y_{self.name}")
        self._F_insr_z = me.dynamicsymbols(f"F_insr_z_{self.name}")

    @property
    def origin(self) -> me.Point:
        """The `Point` from which the musculotendon originates."""
        return self._origin

    @origin.setter
    def origin(self, origin: me.Point) -> None:
        if not isinstance(origin, me.Point):
            msg = (
                f'Value {repr(origin)} passed to `origin` was of type '
                f'{type(origin)}, must be {type(me.Point)}.'
            )
            raise TypeError(msg)
        self._origin = origin

    @property
    def insertion(self) -> me.Point:
        """The `Point` to which the musculotendon inserts."""
        return self._insertion

    @insertion.setter
    def insertion(self, insertion: me.Point) -> None:
        if not isinstance(insertion, me.Point):
            msg = (
                f'Value {repr(insertion)} passed to `insertion` was of type '
                f'{type(insertion)}, must be {type(me.Point)}.'
            )
            raise TypeError(msg)
        self._insertion = insertion

    @property
    def optimal_fiber_length(self) -> Optional[float]:
        """Numeric value representing the optimal fiber length parameter.

        The value by which the musculotendon's fiber force-length relationships
        are normalized. It corresponds to the fiber length at which the
        musculotendon can produce its maximal force. This value maps to the
        symbol attribute `l_M_opt`.

        See Also
        ========
        `l_M_opt`: The symbol to which this numeric value maps.

        """
        return self._optimal_fiber_length

    @optimal_fiber_length.setter
    def optimal_fiber_length(
        self,
        optimal_fiber_length: Optional[float],
    ) -> None:
        if optimal_fiber_length is None:
            self._optimal_fiber_length = optimal_fiber_length
        else:
            optimal_fiber_length = float(optimal_fiber_length)
            if optimal_fiber_length <= 0.0:
                msg = (
                    f'`optimal_fiber_length` of {optimal_fiber_length} must be '
                    f'greater than {0.0}.'
                )
                raise ValueError(msg)
            self._optimal_fiber_length = optimal_fiber_length

    @property
    def maximal_fiber_velocity(self) -> float:
        """Numeric value representing the maximal fiber velocity parameter.

        The value by which the musculotendon's fiber force-velcoity relationship
        is normalized. It corresponds to the fiber shortening velocity at which
        the musculotendon is unable to produce any force whatsoever. This value
        maps to the symbol attribute `v_M_max`.

        See Also
        ========
        `v_M_max`: The symbol to which this numeric value maps.

        """
        return self._maximal_fiber_velocity

    @maximal_fiber_velocity.setter
    def maximal_fiber_velocity(
        self,
        maximal_fiber_velocity: Optional[float],
    ) -> None:
        maximal_fiber_velocity = float(maximal_fiber_velocity)
        if maximal_fiber_velocity <= 0.0:
            msg = (
                f'`maximal_fiber_velocity` of {maximal_fiber_velocity} '
                f'must be greater than {0.0}.'
            )
            raise ValueError(msg)
        self._maximal_fiber_velocity = maximal_fiber_velocity

    @property
    def peak_isometric_force(self) -> Optional[float]:
        """Numeric value representing the peak isometric force parameter.

        The force the musculotendon can produce when the fiber length is equal
        to the optimal fiber length, the shortening velocity is zero, and at
        full activation. "Isometric" refers to the fact that the musculotendon
        is undergoing an isometric contraction (i.e. maintaining length) at this
        magnitude of force. This value maps to the symbol attribute `F_M_max`.

        See Also
        ========
        `F_M_max`: The symbol to which this numeric value maps.

        """
        return self._peak_isometric_force

    @peak_isometric_force.setter
    def peak_isometric_force(
        self,
        peak_isometric_force: Optional[float],
    ) -> None:
        if peak_isometric_force is None:
            self._peak_isometric_force = peak_isometric_force
        else:
            peak_isometric_force = float(peak_isometric_force)
            if peak_isometric_force <= 0.0:
                msg = (
                    f'`peak_isometric_force` of {peak_isometric_force} must be '
                    f'greater than {0.0}.'
                )
                raise ValueError(msg)
            self._peak_isometric_force = peak_isometric_force

    @property
    def tendon_slack_length(self) -> Optional[float]:
        """Numeric value representing the tendon slack length parameter.

        The value by which the musculotendon's tendon force-length relationship
        is normalized. It corresponds to the tendon length at which the tendon
        is unable to produce any tensile force due to becoming slack. For a
        rigid tendon model, this value is exactly equal to the tendon length.
        This value maps to the symbol attribute `l_T_slack`.

        See Also
        ========
        `l_T_slack`: The symbol to which this numeric value maps.

        """
        return self._tendon_slack_length

    @tendon_slack_length.setter
    def tendon_slack_length(
        self,
        tendon_slack_length: Optional[float],
    ) -> None:
        if tendon_slack_length is None:
            self._tendon_slack_length = tendon_slack_length
        else:
            tendon_slack_length = float(tendon_slack_length)
            if tendon_slack_length <= 0.0:
                msg = (
                    f'`tendon_slack_length` of {tendon_slack_length} must be '
                    f'greater than {0.0}.'
                )
                raise ValueError(msg)
            self._tendon_slack_length = tendon_slack_length

    @property
    def optimal_pennation_angle(self) -> float:
        """Numeric value representing the optimal pennation angle parameter.

        The value of pennation angle when the musculotendon's fiber has length
        equal to the optimal fiber length. For unpennated fibers, this value is
        zero. This value maps to the symbol attribute `alpha_opt`.

        See Also
        ========
        `alpha_opt`: The symbol to which this numeric value maps.

        """
        return self._optimal_pennation_angle

    @optimal_pennation_angle.setter
    def optimal_pennation_angle(
        self,
        optimal_pennation_angle: Optional[float],
    ) -> None:
        optimal_pennation_angle = float(optimal_pennation_angle)
        if optimal_pennation_angle < 0.0:
            msg = (
                f'`optimal_pennation_angle` of {optimal_pennation_angle} must '
                f'be greater than or equal to {0.0}.'
            )
            raise ValueError(msg)
        self._optimal_pennation_angle = optimal_pennation_angle

    @property
    def fiber_damping_coefficient(self) -> float:
        """Numeric value representing the fiber damping coefficient parameter.

        The value of damping in the musculotendon's fiber. While not strictly
        physical, this parameter is important for producing numerically stable
        musculotendon dynamics. For undamped fibers, this value is zero. This
        value maps to the symbol attribute `beta`.

        See Also
        ========
        `beta`: The symbol to which this numeric value maps.

        """
        return self._fiber_damping_coefficient

    @fiber_damping_coefficient.setter
    def fiber_damping_coefficient(
        self,
        fiber_damping_coefficient: float,
    ) -> None:
        fiber_damping_coefficient = float(fiber_damping_coefficient)
        if fiber_damping_coefficient < 0.0:
            msg = (
                f'`fiber_damping_coefficient` of {fiber_damping_coefficient} '
                f'must be greater than or equal to {0.0}.'
            )
            raise ValueError(msg)
        self._fiber_damping_coefficient = fiber_damping_coefficient

    @property
    def l_M_opt(self) -> sm.Symbol:
        """Accessor for the optimal fiber length symbol.

        See Also
        ========

        `optimal_fiber_length`: The optional numeric value to which this symbol
            maps.

        """
        return self._l_M_opt

    @property
    def v_M_max(self) -> sm.Symbol:
        """Accessor for the maximal fiber velocity symbol.

        See Also
        ========

        `maximal_fiber_velocity`: The optional numeric value to which this
            symbol maps.

        """
        return self._v_M_max

    @property
    def F_M_max(self) -> sm.Symbol:
        """Accessor for the peak isometric force symbol.

        See Also
        ========

        `peak_isometric_force`: The optional numeric value to which this symbol
            maps.

        """
        return self._F_M_max

    @property
    def l_T_slack(self) -> sm.Symbol:
        """Accessor for the tendon slack length symbol.

        See Also
        ========

        `tendon_slack_length`: The optional numeric value to which this symbol
            maps.

        """
        return self._l_T_slack

    @property
    def alpha_opt(self) -> sm.Symbol:
        """Accessor for the optimal pennation angle symbol.

        See Also
        ========

        `optimal_pennation_angle`: The optional numeric value to which this
            symbol maps.

        """
        return self._alpha_opt

    @property
    def beta(self) -> sm.Symbol:
        """Accessor for the fiber damping coefficient symbol.

        See Also
        ========

        `fiber_damping_coefficient`: The optional numeric value to which this
            symbol maps.

        """
        return self._beta

    @property
    def l_MT(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the musculotendon length dynamic symbol."""
        return self._l_MT

    @property
    def v_MT(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the musculotendon shortening velocity dynamic symbol."""
        return self._v_MT

    @property
    def l_T(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the tendon length dynamic symbol."""
        return self._l_T

    @property
    def v_T(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the tendon shortening velocity dynamic symbol."""
        return self._v_T

    @property
    def l_M(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the fiber length dynamic symbol."""
        return self._l_M

    @property
    def v_M(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the fiber shortening velocity dynamic symbol."""
        return self._v_M

    @property
    def l_T_tilde(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the normalized tendon length dynamic symbol."""
        return self._l_T_tilde

    @property
    def v_T_tilde(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the normalized tendon shortening velocity dynamic
        symbol.

        """
        return self._v_T_tilde

    @property
    def l_M_tilde(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the normalized fiber length dynamic symbol."""
        return self._l_M_tilde

    @property
    def v_M_tilde(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the normalized fiber shortening velocity dynamic symbol.

        """
        return self._v_M_tilde

    @property
    def F_T(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the tendon force dynamic symbol."""
        return self._F_T

    @property
    def F_M(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the fiber force dynamic symbol."""
        return self._F_M

    @property
    def F_T_tilde(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the normalized tendon force dynamic symbol."""
        return self._F_T_tilde

    @property
    def F_M_tilde(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the normalized fiber force dynamic symbol."""
        return self._F_M_tilde

    @property
    def fl_T(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the tendon force-length dynamic symbol."""
        return self._fl_T

    @property
    def fl_M_pas(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the passive fiber force-length dynamic symbol."""
        return self._fl_M_pas

    @property
    def fl_M_act(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the active fiber force-length dynamic symbol."""
        return self._fl_M_act

    @property
    def fv_M(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the fiber force-velocity dynamic symbol."""
        return self._fv_M

    @property
    def cos_alpha(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the cosine of the pennation angle dynamic symbol."""
        return self._cos_alpha

    @property
    def dF_T_tilde_dt(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the rate of change in tendon force dynamic symbol."""
        return self._dF_T_tilde_dt

    @property
    def dl_M_tilde_dt(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the rate of change in fiber length dynamic symbol."""
        return self._dl_M_tilde_dt

    @property
    def F_orig_x(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the x-axis origin force dynamic symbol."""
        return self._F_orig_x

    @property
    def F_orig_y(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the y-axis origin force dynamic symbol."""
        return self._F_orig_y

    @property
    def F_orig_z(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the z-axis origin force dynamic symbol."""
        return self._F_orig_z

    @property
    def F_insr_x(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the x-axis insertion force dynamic symbol."""
        return self._F_insr_x

    @property
    def F_insr_y(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the y-axis insertion force dynamic symbol."""
        return self._F_insr_y

    @property
    def F_insr_z(self) -> sm.core.backend.AppliedUndef:
        """Accessor for the z-axis insertion force dynamic symbol."""
        return self._F_insr_z

    def __str__(self):
        """String representation of the musculotendon model instance."""
        return f'{self.__class__.__name__}({repr(self.name)})'

    def __repr__(self):
        """Printable representation of the musculotendon model instance."""
        return (
            f'{self.__class__.__name__}({repr(self.name)}, '
            f'origin={repr(self.origin)}, insertion={repr(self.insertion)}), '
            f'optimal_fiber_length={repr(self.optimal_fiber_length)}, '
        )


class Brockie2021Musculotendon(MusculotendonBase):
    """Musculotendon model with curves based on Brockie, 2021 [1].

    References
    ==========

    .. [1] Brockie, S. G., Predictive Simulation of Musculoskeletal Models Using
           Direct Collocation, PhD Thesis, University of Cambridge, (2021) pp.
           234-241

    """


class DeGroote2016Musculotendon(MusculotendonBase):
    """Musculotendon model with curves based on De Groote et al., 2016 [1].

    References
    ==========

    .. [1] De Groote, F., Kinney, A. L., Rao, A. V., & Fregly, B. J., Evaluation
           of direct collocation optimal control problem formulations for
           solving the muscle redundancy problem, Annals of biomedical
           engineering, 44(10), (2016) pp. 2922-2936

    """


class Millard2013Musculotendon(MusculotendonBase):
    """Musculotendon model with curves based on Millard et al., 2013 [1].

    References
    ==========

    .. [1] Millard, M., Uchida, T., Seth, A., Delp, S. L., Flexing computational
           muscle: modeling and simulation of musculotendon dynamics. Journal of
           biomechanical engineering, 135(2), (2013)

    """


def Musculotendon(
    name: str,
    identifier: str,
    *args,
    **kwargs,
) -> MusculotendonBase:
    """Factory function for instantiating a specific musculotendon model by its
    name as an identifier.

    Explanation
    ===========

    Sometimes it can be easier to use a factory function to instantiate a
    specific subclass of a parent class. This function lets users easily change
    which musculotendon model they are instantiating simply by changing the
    string identifier passed to the `identifier` argument. This is useful in
    instances where a user wants to be able to change the type of all
    musculotendon models in their system simply by changing one variable in
    their code.

    Examples
    ========

    Instantiate a `DeGroote2016Musculotendon object` given origin and insertion
    points:

    >>> origin = me.Point('origin')
    >>> insertion = me.Insertion('insertion')
    >>> Musculotendon('muscle', 'DeGroote', origin=origin, insertion=insertion)
    DeGroote2016Musculotendon('muscle', origin=origin, insertion=insertion)

    Parameters
    ==========

    name : str
        The name identifier associated with the musculotendon. This name is used
        as a suffix when automatically generated symbols are instantiated. It
        must be a string of nonzero length.
    identifier : str
        Identifier used to choose which musculotendon class to instantiate. Must
        be one of 'Brockie', 'DeGroote', or 'Millard'.
    *args
        Positional arguments to be passed to the musculotendon class's
        initializer.
    **kwargs
        Keyword arguments to be passed to the musculotendon class's initialiser.

    Returns
    =======

    An instance of either Brockie2021Musculotendon, DeGroote2016Musculotendon,
    or Millard2013Musculotendon.

    Raises
    ======

    ValueError
        If an invalid musculotendon model identifier is passed to
        `identifier`.

    """
    if identifier.lower() in {'brockie', 'brockie2021'}:
        return Brockie2021Musculotendon(name, *args, **kwargs)
    elif identifier.lower() in {'degroote', 'degroote2016'}:
        return DeGroote2016Musculotendon(name, *args, **kwargs)
    elif identifier.lower() in {'millard', 'millard2013'}:
        return Millard2013Musculotendon(name, *args, **kwargs)
    msg = (f'{name} is an invalid musculotendon model identifier, must be one '
           'of: "Brockie", "DeGroote", or "Millard".')
    raise ValueError(msg)
