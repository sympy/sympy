
class BasicType(type):

    classnamespace = dict()
    def __init__(cls,*args,**kws):
        if not cls.undefined_Function:
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
    undefined_Function = False

    @staticmethod
    def sympify(a):
        if isinstance(a, Basic):
            return a
        elif isinstance(a, bool):
            raise NotImplementedError("bool support")
        elif isinstance(a, (int, long)):
            return Basic.Integer(a)
        raise ValueError("%s is NOT a valid SymPy expression" % `a`)

    def __new__(cls, *args, **kwds):
        r = cls.canonize(*args, **kwds)
        if isinstance(r, Basic): return r
        if r is not None:
            args, kwds = r
        if kwds:
            return cls._new(cls, *args, **kwds)
        return cls._new(cls, *args)

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

class Composite(Basic):

    pass

class CompositeTuple(Composite, tuple):

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

class CompositeDict(Composite, dict):

    @staticmethod
    def _new(cls, *args, **kwds):
        obj = dict.__new__(cls)
        obj.update(args[0])
        obj.update(kwds)
        return obj

    def __init__(self, *args, **kwds):
        pass
        
    @classmethod
    def canonize(cls, *args, **kwds):
        if len(args)==1 and isinstance(args[0],dict):
            return (args[0].copy(),), kwds
        d = dict()
        for a in args:
            if not isinstance(a, Basic):
                if isinstance(a, dict):
                    for k,v in a.items():
                        try:
                            d[k] += v
                        except KeyError:
                            d[k] = v
                    continue
                a = Basic.sympify(a)
            try:
                d[a] += 1
            except KeyError:
                d[a] = 1
        if not d: return Basic.Integer(0)
        if len(d)==1:
            k,v = d.popitem()
            if v==1:
                return Basic.sympify(k)
            return Basic.sympify(v*k)
        return (d,), kwds

    def torepr(self):
        return '%s(%s)' % (self.__class__.__name__, dict(self))
