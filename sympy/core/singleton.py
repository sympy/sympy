"""Singleton mechanism"""

from __future__ import division, print_function

from .core import Registry
from .assumptions import ManagedProperties
from .sympify import sympify


class SingletonRegistry(Registry):
    """
    A map from singleton classes to the corresponding instances.
    E.g. S.Exp == Exp()
    """
    __slots__ = []

    # Also allow things like S(5)
    __call__ = staticmethod(sympify)

    def __init__(self):
        self._classes_to_install = {}
        # Dict of classes that have been registered, but that have not have been
        # installed as an attribute of this SingletonRegistry.
        # Installation automatically happens at the first attempt to access the
        # attribute.
        # The purpose of this is to allow registration during class
        # initialization during import, but not trigger object creation until
        # actual use (which should not happen until after all imports are
        # finished).

    def register(self, cls):
        self._classes_to_install[cls.__name__] = cls

    def __getattr__(self, name):
        """Python calls __getattr__ if no attribute of that name was installed
        yet.

        This __getattr__ checks whether a class with the requested name was
        already registered but not installed; if no, raises an AttributeError.
        Otherwise, retrieves the class, calculates its singleton value, installs
        it as an attribute of the given name, and unregisters the class."""
        if name not in self._classes_to_install:
            raise AttributeError(
                "Attribute '%s' was not installed on SymPy registry %s" % (
                name, self))
        class_to_install = self._classes_to_install[name]
        value_to_install = class_to_install()
        self.__setattr__(name, value_to_install)
        del self._classes_to_install[name]
        return value_to_install

    def __repr__(self):
        return "S"

S = SingletonRegistry()


class Singleton(ManagedProperties):
    """
    Metaclass for singleton classes.

    A singleton class has only one instance which is returned every time the
    class is instantiated. Additionally, this instance can be accessed through
    the global registry object S as S.<class_name>.

    Examples
    ========

        >>> from sympy import S, Basic
        >>> from sympy.core.singleton import Singleton
        >>> from sympy.core.compatibility import with_metaclass
        >>> class MySingleton(with_metaclass(Singleton, Basic)):
        ...     pass
        >>> Basic() is Basic()
        False
        >>> MySingleton() is MySingleton()
        True
        >>> S.MySingleton is MySingleton()
        True

    Notes
    =====

    Instance creation is delayed until the first time the value is accessed.
    (SymPy versions before 0.7.7 would create the instance during class
    creation time, which would be prone to import cycles.)

    This metaclass is a subclass of ManagedProperties because that is the
    metaclass of many classes that need to be Singletons (Python does not allow
    subclasses to have a different metaclass than the superclass, except the
    subclass may use a subclassed metaclass).
    """

    _instances = {}
    "Maps singleton classes to their instances."

    def __new__(cls, *args, **kwargs):
        result = super(Singleton, cls).__new__(cls, *args, **kwargs)
        S.register(result)
        return result

    def __call__(self, *args, **kwargs):
        # Called when application code says SomeClass(), where SomeClass is a
        # class of which Singleton is the metaclas.
        # __call__ is invoked first, before __new__() and __init__().
        if self not in Singleton._instances:
            Singleton._instances[self] = \
                super(Singleton, self).__call__(*args, **kwargs)
                # Invokes the standard constructor of SomeClass.
        return Singleton._instances[self]

        # Inject pickling support.
        def __getnewargs__(self):
            return ()
        self.__getnewargs__ = __getnewargs__
