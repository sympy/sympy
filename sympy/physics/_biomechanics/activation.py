"""Activation dynamics for musclotendon models."""

from abc import ABC, abstractmethod

from sympy.physics._biomechanics._mixin import _NamedMixin
from sympy.physics.mechanics import dynamicsymbols


__all__ = ['ActivationBase']


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
