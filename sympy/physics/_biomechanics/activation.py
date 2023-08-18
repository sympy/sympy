"""Activation dynamics for musclotendon models.

Musculotendon models are able to produce active force when they are activated,
which is when a chemical process has taken place within the muscle fibers
causing them to voluntarily contract. Biologically this chemical process (the
diffusion of $Ca^{2+}$ ions) is not the input in the system, electrical signals
from the nervous system are. These are termed excitations. Activation dynamics,
which relates the normalized excitation level to the normalized activation
level, can be modeled by the models present in this module.

"""

from abc import ABC, abstractmethod

from sympy.matrices import Matrix
from sympy.matrices.dense import zeros
from sympy.physics._biomechanics._mixin import _NamedMixin
from sympy.physics.mechanics import dynamicsymbols


__all__ = [
    'ActivationBase',
    'ZerothOrderActivation',
]


class ActivationBase(ABC, _NamedMixin):
    """Abstract base class for all activation dynamics classes to inherit from.

    Notes
    =====

    Instances of this class cannot be directly instantiated by users. However,
    it can be used to created custom activation dynamics types through
    subclassing.

    """

    def __init__(self, name):
        """Initializer for ``ActivationBase``."""
        self.name = str(name)

        # Symbols
        self._e = dynamicsymbols(f"e_{name}")
        self._a = dynamicsymbols(f"a_{name}")

    @property
    def excitation(self):
        """Dynamic symbol representing excitation.

        Explanation
        ===========

        The alias `e` can also be used to access the same attribute.

        """
        return self._e

    @property
    def e(self):
        """Dynamic symbol representing excitation.

        Explanation
        ===========

        The alias `excitation` can also be used to access the same attribute.

        """
        return self._e

    @property
    def activation(self):
        """Dynamic symbol representing activation.

        Explanation
        ===========

        The alias `a` can also be used to access the same attribute.

        """
        return self._a

    @property
    def a(self):
        """Dynamic symbol representing activation.

        Explanation
        ===========

        The alias `activation` can also be used to access the same attribute.

        """
        return self._a

    @property
    @abstractmethod
    def order(self):
        """Order of the (differential) equation governing activation."""
        pass

    @property
    @abstractmethod
    def state_vars(self):
        """Ordered column matrix of functions of time that represent the state
        variables.

        Explanation
        ===========

        The alias `x` can also be used to access the same attribute.

        """
        pass

    @property
    @abstractmethod
    def x(self):
        """Ordered column matrix of functions of time that represent the state
        variables.

        Explanation
        ===========

        The alias `state_vars` can also be used to access the same attribute.

        """
        pass

    @property
    @abstractmethod
    def input_vars(self):
        """Ordered column matrix of functions of time that represent the input
        variables.

        Explanation
        ===========

        The alias `r` can also be used to access the same attribute.

        """
        pass

    @property
    @abstractmethod
    def r(self):
        """Ordered column matrix of functions of time that represent the input
        variables.

        Explanation
        ===========

        The alias `input_vars` can also be used to access the same attribute.

        """
        pass

    @property
    @abstractmethod
    def constants(self):
        """Ordered column matrix of non-time varying symbols present in ``M``
        and ``F``.

        Explanation
        ===========

        The alias `p` can also be used to access the same attribute.

        """
        pass

    @property
    @abstractmethod
    def p(self):
        """Ordered column matrix of non-time varying symbols present in ``M``
        and ``F``.

        Explanation
        ===========

        The alias `constants` can also be used to access the same attribute.

        """
        pass

    @property
    @abstractmethod
    def M(self):
        """Ordered square matrix of coefficients on the LHS of ``M x' = F``.

        Explanation
        ===========

        The square matrix that forms part of the LHS of the linear system of
        ordinary differential equations governing the activation dynamics:

        ``M(x, r, t, p) x' = F(x, r, t, p)``.

        """
        pass

    @property
    @abstractmethod
    def F(self):
        """Ordered column matrix of equations on the RHS of ``M x' = F``.

        Explanation
        ===========

        The column matrix that forms the RHS of the linear system of ordinary
        differential equations governing the activation dynamics:

        ``M(x, r, t, p) x' = F(x, r, t, p)``.

        """
        pass

    @abstractmethod
    def rhs(self):
        """

        Explanation
        ===========

        The solution to the linear system of ordinary differential equations
        governing the activation dynamics:

        ``M(x, r, t, p) x' = F(x, r, t, p)``.

        """
        pass

    def __repr__(self):
        """Default representation of activation dynamics."""
        return f'{self.__class__.__name__}({self.name!r})'


class ZerothOrderActivation(ActivationBase):
    """Simple zeroth-order activation dynamics mapping excitation to
    activation.

    Explanation
    ===========

    Zeroth-order activation dynamics are useful in instances where you want to
    reduce the complexity of your musculotendon dynamics as they simple map
    exictation to activation. As a result, no additional state equations are
    introduced to your system. They also remove a potential source of delay
    between the input and dynamics of your system as no (ordinary) differential
    equations are involed.

    """

    def __init__(self, name):
        """Initializer for ``ZerothOrderActivation``.

        Parameters
        ==========

        name : str
            The name identifier associated with the instance. Must be a string
            of length at least 1.

        """
        super().__init__(name)

        # Zeroth-order activation dynamics has activation equal excitation so
        # overwrite the symbol for activation with the excitation symbol.
        self._a = self._e

    @property
    def order(self):
        """Order of the (differential) equation governing activation."""
        return 0

    @property
    def state_vars(self):
        """Ordered column matrix of functions of time that represent the state
        variables.

        Explanation
        ===========

        As zeroth-order activation dynamics simply maps excitation to
        activation, this class has no associated state variables and so this
        property return an empty column ``Matrix`` with shape (0, 1).

        The alias `x` can also be used to access the same attribute.

        """
        return zeros(0, 1)

    @property
    def x(self):
        """Ordered column matrix of functions of time that represent the state
        variables.

        Explanation
        ===========

        As zeroth-order activation dynamics simply maps excitation to
        activation, this class has no associated state variables and so this
        property return an empty column ``Matrix`` with shape (0, 1).

        The alias `state_vars` can also be used to access the same attribute.

        """
        return zeros(0, 1)

    @property
    def input_vars(self):
        """Ordered column matrix of functions of time that represent the input
        variables.

        Explanation
        ===========

        Excitation is the only input in zeroth-order activation dynamics and so
        this property returns a column ``Matrix`` with one entry, ``e``, and
        shape (1, 1).

        The alias `r` can also be used to access the same attribute.

        """
        return Matrix([self._e])

    @property
    def r(self):
        """Ordered column matrix of functions of time that represent the input
        variables.

        Explanation
        ===========

        Excitation is the only input in zeroth-order activation dynamics and so
        this property returns a column ``Matrix`` with one entry, ``e``, and
        shape (1, 1).

        The alias `input_vars` can also be used to access the same attribute.

        """
        return Matrix([self._e])

    @property
    def constants(self):
        """Ordered column matrix of non-time varying symbols present in ``M``
        and ``F``.

        Explanation
        ===========

        As zeroth-order activation dynamics simply maps excitation to
        activation, this class has no associated constants and so this property
        return an empty column ``Matrix`` with shape (0, 1).

        The alias `p` can also be used to access the same attribute.

        """
        return zeros(0, 1)

    @property
    def p(self):
        """Ordered column matrix of non-time varying symbols present in ``M``
        and ``F``.

        Explanation
        ===========

        As zeroth-order activation dynamics simply maps excitation to
        activation, this class has no associated constants and so this property
        return an empty column ``Matrix`` with shape (0, 1).

        The alias `constants` can also be used to access the same attribute.

        """
        return zeros(0, 1)

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
        return Matrix([])

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
        return zeros(0, 1)

    def rhs(self):
        """

        Explanation
        ===========

        The solution to the linear system of ordinary differential equations
        governing the activation dynamics:

        ``M(x, r, t, p) x' = F(x, r, t, p)``.

        As zeroth-order activation dynamics have no state variables, this
        linear has dimension 0 and therefore this method returns an empty
        column ``Matrix`` with shape (0, 1).

        """
        return zeros(0, 1)
