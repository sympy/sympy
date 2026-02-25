from abc import ABC, abstractmethod

class _Methods(ABC):
    """Abstract Base Class for all methods."""

    @property
    @abstractmethod
    def q(self):
        """Generalized coordinates of the system."""
        ...

    @property
    @abstractmethod
    def u(self):
        """Generalized speeds of the system."""
        ...

    @property
    @abstractmethod
    def bodies(self):
        """Bodies in the system."""
        ...

    @property
    @abstractmethod
    def loads(self):
        """Loads applied to the system."""
        ...

    @property
    @abstractmethod
    def mass_matrix(self):
        """The mass matrix of the system."""
        ...

    @property
    @abstractmethod
    def forcing(self):
        """The forcing vector of the system."""
        ...

    @property
    @abstractmethod
    def mass_matrix_full(self):
        """The mass matrix of the system, augmented by the kinematic equations."""
        ...

    @property
    @abstractmethod
    def forcing_full(self):
        """The forcing vector of the system, augmented by the kinematic equations."""
        ...

    def _form_eoms(self):
        raise NotImplementedError("Subclasses must implement this.")
