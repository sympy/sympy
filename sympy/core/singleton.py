"""Singleton mechanism"""

from __future__ import print_function, division

from .core import Registry
from .assumptions import ManagedProperties
from .sympify import sympify


class SingletonRegistry(Registry):
    """
    A map between singleton classes and the corresponding instances.
    E.g. S.Exp == C.Exp()
    """
    __slots__ = []

    __call__ = staticmethod(sympify)

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

    ** Developer notes **
        The class is instantiated immediately at the point where it is defined
        by calling cls.__new__(cls). This instance is cached and cls.__new__ is
        rebound to return it directly.

        The original constructor is also cached to allow subclasses to access it
        and have their own instance.

    """

    def __init__(cls, name, bases, dict_):
        super(Singleton, cls).__init__(cls, name, bases, dict_)

        for ancestor in cls.mro():
            if '__new__' in ancestor.__dict__:
                break
        if isinstance(ancestor, Singleton) and ancestor is not cls:
            ctor = ancestor._new_instance
        else:
            ctor = cls.__new__
        cls._new_instance = staticmethod(ctor)

        the_instance = ctor(cls)

        def __new__(cls):
            return the_instance
        cls.__new__ = staticmethod(__new__)

        setattr(S, name, the_instance)

        # Inject pickling support.
        def __getnewargs__(self):
            return ()
        cls.__getnewargs__ = __getnewargs__
