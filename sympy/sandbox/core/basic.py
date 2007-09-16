
class BasicType(type):

    classnamespace = dict()
    def __init__(cls,*args,**kws):
        # todo: save only classes that are defined in core.
        n = cls.__name__
        c = BasicType.classnamespace.get(n)
        if c is None:
            BasicType.classnamespace[n] = cls
        else:
            print 'Ignoring redefinition of %s: %s defined earlier than %s' % (n, c, cls)
        type.__init__(cls, *args, **kws)

    def __getattr__(cls, name):
        try: return BasicType.classnamespace[name]
        except KeyError: pass
        raise AttributeError("'%s' object has no attribute '%s'"%
                             (cls.__name__, name))

class Basic(object):

    __metaclass__ = BasicType

    @staticmethod
    def sympify(a):
        if isinstance(a, Basic):
            return a
        elif isinstance(a, bool):
            raise NotImplementedError("bool support")
        elif isinstance(a, (int, long)):
            return Basic.Integer(a)
        raise ValueError("%r is NOT a valid SymPy expression" % a)

    def __new__(cls, *args, **kwds):
        r = cls.canonize(*args, **kwds)
        if isinstance(r, Basic): return r
        if r is not None:
            args, kwds = r
        return cls._new(cls, *args, **kwds)

    @staticmethod
    def _new(cls, *args, **kwds):
        assert not args,`args`
        obj = object.__new__(cls)
        obj.__dict__.update(kwds)
        return obj

    @classmethod
    def canonize(cls, *args, **kwds):
        return args, kwds

    def __repr__(self):
        if isinstance(self, type):
            return self.__class__.torepr(self)
        return self.torepr()


class Atom(Basic):

    def torepr(self):
        return '%s()' % (self.__class__.__name__)


class Composite(Basic, tuple):

    _new = tuple.__new__

    @classmethod
    def canonize(cls, *args, **kwds):
        args = tuple(map(Basic.sympify, args))
        return args, kwds

    def torepr(self):
        return '%s(%s)' % (self.__class__.__name__, ', '.join(map(repr, self)))

    def __hash__(self):
        try:
            return self._hash
        except AttributeError:
            h = hash((self.__class__.__name__,)+tuple(self))
            self._hash = h
            return h

