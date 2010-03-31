"""Singleton mechanism"""

from core import BasicMeta, C
from sympify import sympify
from basic import Basic

class SingletonMeta(BasicMeta):
    """Metaclass for all singletons

       All singleton classes should put this into their __metaclass__, and
       _not_ to define __new__

       example:

       class Zero(Integer):
           __metaclass__ = SingletonMeta

           p = 0
           q = 1
    """

    def __init__(cls, *args, **kw):
        BasicMeta.__init__(cls, *args, **kw)

        # we are going to inject singletonic __new__, here it is:
        def cls_new(cls):
            try:
                obj = getattr(SingletonFactory, cls.__name__)

            except AttributeError:
                obj = Basic.__new__(cls, *(), **{})
                setattr(SingletonFactory, cls.__name__, obj)

            return obj

        cls_new.__name__ = '%s.__new__' % (cls.__name__)

        assert not cls.__dict__.has_key('__new__'), \
                'Singleton classes are not allowed to redefine __new__'

        # inject singletonic __new__
        cls.__new__      = staticmethod(cls_new)

        # Inject pickling support.
        def cls_getnewargs(self):
            return ()
        cls_getnewargs.__name__ = '%s.__getnewargs__' % cls.__name__

        assert not cls.__dict__.has_key('__getnewargs__'), \
                'Singleton classes are not allowed to redefine __getnewargs__'
        cls.__getnewargs__ = cls_getnewargs


        # tag the class appropriately (so we could verify it later when doing
        # S.<something>
        cls.is_Singleton = True


class SingletonFactory:
    """
    A map between singleton classes and the corresponding instances.
    E.g. S.Exp == C.Exp()
    """

    def __getattr__(self, clsname):
        if clsname == "__repr__":
            return lambda: "S"

        cls = getattr(C, clsname)
        assert cls.is_Singleton
        obj = cls()

        # store found object in own __dict__, so the next lookups will be
        # serviced without entering __getattr__, and so will be fast
        setattr(self, clsname, obj)
        return obj

S = SingletonFactory()

S.__call__ = sympify
