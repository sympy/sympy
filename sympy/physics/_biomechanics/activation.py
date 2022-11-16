r"""Implementations of, and factories for, activation dynamics models.

Musculotendon models are able to produce active force when they are activated,
that is when a chemical process has taken place within the muscle fibers causing
them to voluntarily contract. Biologically this chemical process (the diffusion
of $Ca^{2+}$ ions) is not the control in the system, electrical signals from the
nervous system are. These are termed excitations. Activation dynamics, which
relate the normalized excitation level to the normalized activation level, can
be modelled by the models present in this module.

"""


from abc import ABC
from typing import Union

from sympy.core.backend import Symbol
from sympy.physics.mechanics import dynamicsymbols
from sympy.physics._biomechanics.mixin import _NamedMixin


__all__ = [
    'ActivationDynamics',
    'DeGroote2016ActivationDynamics',
    'ZerothOrderActivationDynamics',
]


class ActivationDynamicsBase(ABC, _NamedMixin):
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
    order : int
        The order of the ordinary differential equation (ODEs) governing the
        activation dynamics.
    name : str
        The name identifier associated with the activation dynamics. This name
        is used as a suffix when automatically generated symbols are
        instantiated. It must be a string of nonzero length.

    """

    def __init__(self, name):
        """Initializer for all activation dynamics instance attributes.

        Parameters
        ==========
        name : str
            The name identifier associated with the activation dynamics. This
            name is used as a suffix when automatically generated symbols are
            instantiated. It must be a string of nonzero length.

        """
        self.name = name

        # Symbols
        self.a = dynamicsymbols(f"a_{name}")
        self.e = dynamicsymbols(f"e_{name}")

    @property
    def order(self):
        """Order of the ODEs governing the activation dynamics."""
        return self._ORDER

    def symbol_to_constant_mapping(self) -> dict[Symbol, float]:
        """Mapping between symbolic attributes and associated constants."""
        return {}

    def __str__(self):
        """String representation of the activation dynamics instance."""
        return f'{self.__class__.__name__}({repr(self.name)})'

    def __repr__(self):
        """Printable representation of the activation dynamics instance."""
        return f'{self.__class__.__name__}({repr(self.name)})'


class ZerothOrderActivationDynamics(ActivationDynamicsBase):
    """Activation dynamics model that directly maps activation to excitation."""

    _ORDER = 0

    def __init__(self, name):
        """Initializer for zeroth-order activation dynamics instance attributes.

        """
        super().__init__(name)


class DeGroote2016ActivationDynamics(ActivationDynamicsBase):
    """Activation dynamics model with curves based on De Groote et al., 2016
    [1].

    Attributes
    ==========
    activation time constant : float, defaults to 0.015
        The value of the activation time constant governing the delay between
        excitation and activation when excitation exceeds activation. This value
        maps to the symbol attribute `tau_a`.
    deactivation time constant : float, defaults to 0.060
        The value of the deactivation time constant governing the delay between
        excitation and activation when activation exceeds excitation. This value
        maps to the symbol attribute `tau_d`.
    tau_a : `sympy.Symbol`
        Accessor for the activation time constant symbol.
    tau_d : `sympy.Symbol`
        Accessor for the deactivation time constant symbol.

    References
    ==========

    .. [1] De Groote, F., Kinney, A. L., Rao, A. V., & Fregly, B. J., Evaluation
           of direct collocation optimal control problem formulations for
           solving the muscle redundancy problem, Annals of biomedical
           engineering, 44(10), (2016) pp. 2922-2936

    """

    _ORDER = 1

    def __init__(
        self,
        name: str,
        activation_time_constant: float = 0.015,
        deactivation_time_constant: float = 0.060,
    ):
        """Initializer for De Groote 2016 activation dynamics instance
        attributes.

        Parameters
        ==========
        activation_time_constant : float, default is 0.015
            The value of the activation time constant governing the delay
            between excitation and activation when excitation exceeds
            activation. This value maps to the symbol attribute `tau_a`.
        deactivation_time_constant : float, default is 0.060
            The value of the deactivation time constant governing the delay
            between excitation and activation when activation exceeds
            excitation. This value maps to the symbol attribute `tau_d`.

        """
        super().__init__(name)

        # Symbols
        self._tau_a = Symbol(f'tau_a_{self.name}')
        self._tau_d = Symbol(f'tau_d_{self.name}')

        # Constants
        self.activation_time_constant = activation_time_constant
        self.deactivation_time_constant = deactivation_time_constant

    @property
    def activation_time_constant(self) -> float:
        """Numeric value representing the activation time constant parameter.

        The value of the activation time constant governing the delay between
        excitation and activation when excitation exceeds activation. This value
        maps to the symbol attribute `tau_a`.

        See Also
        ========
        `tau_a`: The symbol to which this numeric value maps.

        """
        return self._activation_time_constant

    @activation_time_constant.setter
    def activation_time_constant(
        self,
        activation_time_constant: float,
    ) -> None:
        MINIMUM = 0.0
        activation_time_constant = float(activation_time_constant)
        if activation_time_constant <= MINIMUM:
            msg = (
                f'`activation_time_constant` of {activation_time_constant} '
                f'must be greater than {MINIMUM}.'
            )
            raise ValueError(msg)
        self._activation_time_constant = activation_time_constant

    @property
    def deactivation_time_constant(self) -> float:
        """Numeric value representing the deactivation time constant parameter.

        The value of the deactivation time constant governing the delay between
        excitation and activation when activation exceeds excitation. This value
        maps to the symbol attribute `tau_d`.

        See Also
        ========
        `tau_d`: The symbol to which this numeric value maps.

        """
        return self._deactivation_time_constant

    @deactivation_time_constant.setter
    def deactivation_time_constant(
        self,
        deactivation_time_constant: float,
    ) -> None:
        MINIMUM = 0.0
        deactivation_time_constant = float(deactivation_time_constant)
        if deactivation_time_constant <= MINIMUM:
            msg = (
                f'`deactivation_time_constant` of {deactivation_time_constant} '
                f'must be greater than {MINIMUM}.'
            )
            raise ValueError(msg)
        self._deactivation_time_constant = deactivation_time_constant

    @property
    def tau_a(self) -> Symbol:
        """Accessor for the activation time constant symbol.

        See Also
        ========

        `activation_time_constant`: The numeric value to which this symbol maps.

        """
        return self._tau_a

    @property
    def tau_d(self) -> Symbol:
        """Accessor for the deactivation time constant symbol.

        See Also
        ========

        `deactivation_time_constant`: The numeric value to which this symbol
            maps.

        """
        return self._tau_d

    def symbol_to_constant_mapping(self) -> dict[Symbol, float]:
        """Mapping between symbolic attributes and associated constants."""
        mapping = {
            self.tau_a: self.activation_time_constant,
            self.tau_d: self.deactivation_time_constant,
        }
        return mapping

    def __repr__(self):
        """Printable representation of the activation dynamics instance."""
        return (
            f'{self.__class__.__name__}({repr(self.name)}, '
            f'activation_time_constant={repr(self.activation_time_constant)}), '
            f'deactiation_time_constant={repr(self.deactivation_time_constant)}'
            f')'
        )


def ActivationDynamics(
    name: str,
    identifier: Union[int, str],
    *args,
    **kwargs,
) -> ActivationDynamicsBase:
    """Factory function for instantiating a specific activation dynamics model
    by its name or order as an identifier.

    Explanation
    ===========

    Sometimes it can be easier to use a factory function to instantiate a
    specific subclass of a parent class. This function lets users easily change
    which activation dynamics model they are instantiating simply by changing
    the string or integer identifier passed to the `activation_dynamics`
    argument. This is useful in instances where a user wants to be able to
    change the type of all activation dynamics models in their system simply by
    changing one variable in their code.

    Examples
    ========

    Instantiate a `DeGroote2016ActivationDynamics` object:

    >>> from sympy.physics._biomechanics import ActivationDynamics
    >>> ActivationDynamics('activation', 'DeGroote')
    DeGroote2016ActivationDynamics('activation')

    Parameters
    ==========

    name : str
        The name identifier associated with the musculotendon. This name is used
        as a suffix when automatically generated symbols are instantiated. It
        must be a string of nonzero length.
    identifier : str
        Identifier used to choose which musculotendon class to instantiate. Must
        be one of 'Zeroth', or 'DeGroote'.
    *args
        Positional arguments to be passed to the activation dynamics class's
        initializer.
    **kwargs
        Keyword arguments to be passed to the activation dynamics class's
        initialiser.

    Returns
    =======

    An instance of either `ZerothOrderActivationDynamics`, or
    `DeGroote2016ActivationDynamics`.

    Raises
    ======

    ValueError
        If an invalid activation dynamics model identifier is passed to
        `identifier`.

    """
    ORDER_DEFAULTS: dict[Union[int, str], str] = {0: "zeroth", 1: "degroote"}
    identifier = ORDER_DEFAULTS.get(identifier, identifier)
    if str(identifier).lower() in {"zero", "zeroth"}:
        return ZerothOrderActivationDynamics(name)
    elif str(identifier).lower() in {"degroote", "degroote2016"}:
        return DeGroote2016ActivationDynamics(name, *args, **kwargs)
    msg = (f'{name} is an invalid activation dynamics model identifier, must '
           'be one of: "Zeroth" (alternatively 0), or "DeGroote" '
           '(alternatively 1).')
    raise ValueError(msg)
