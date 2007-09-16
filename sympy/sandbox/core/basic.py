
class Symbolic(object):
    """
    Defines
    """

    _new = object.__new__

    def __repr__(self):
        if isinstance(self, type):
            return self.__class__.torepr(self)
        return self.torepr()

class Basic(Symbolic):

    def __new__(cls, *args, **kwds):
        r = cls.canonize(*args, **kwds)
        if isinstance(r, Basic): return r
        if r is not None:
            args, kwds = r
        return cls._new(cls, *args, **kwds)

    @classmethod
    def canonize(cls, *args, **kwds):
        return args, kwds

class FunctionSignature:

    def __init__(self, argument_classes = (Basic,), value_classes = (Basic,)):
        self.argument_classes = argument_classes
        self.value_classes = value_classes
        if argument_classes is None:
            # unspecified number of arguments
            self.nof_arguments = None
        else:
            self.nof_arguments = len(argument_classes)
        if value_classes is None:
            # unspecified number of arguments
            self.nof_values = None
        else:
            self.nof_values = len(value_classes)

    def validate(self, args):
        if self.nof_arguments is not None:
            if self.nof_arguments!=len(args):
                # todo: improve exception messages
                raise TypeError('wrong number of arguments')
            for a,cls in zip(args, self.argument_classes):
                if not isinstance(a, cls):
                    raise TypeError('wrong argument type %r, expected %s' % (a, cls))

    def __repr__(self):
        return '%s(%r, %r)' % (self.__class__.__name__,
                               self.argument_classes,
                               self.value_classes)

class FunctionClass(Basic, type):

    _new = type.__new__

    def torepr(cls):
        return cls.__name__
        if cls.nofargs is not None:
            return 'Function(%r, %r)' % (cls.__name__, cls.nofargs)
        return 'Function(%r)' % (cls.__name__)

    @classmethod
    def canonize(cls, *args, **kwds):
        if not kwds:
            if len(args)==2:
                basecls = args[0]
                name = args[1]
                d = basecls.__dict__.copy()
                assert isinstance(name, str),`name`
                return (name, (basecls,), d), {}
            elif len(args)==2:
                basecls = args[0]
                name = args[1]
                signature = args[2]
                d = basecls.__dict__.copy()
                d['signature'] = signature
                assert isinstance(name, str),`name`
                return (name, (basecls,), d), {}
        return args, kwds

class Atom(Basic):
    pass

class Composite(Basic, tuple):

    _new = tuple.__new__

    def torepr(self):
        return '%s(%s)' % (self.__class__.__name__, ', '.join(map(repr, self)))

class Symbol(Atom, str):

    _new = str.__new__

    def __new__(cls, name):
        assert isinstance(name, str), `name`
        return cls._new(cls, name)

    def torepr(self):
        return '%s(%r)' % (self.__class__.__name__, self[:])


class Function(Composite):

    __metaclass__ = FunctionClass
    
    signature = FunctionSignature(None, None)
    def __new__(cls, *args, **kwds):
        if cls is Function or cls is Functional:
            r = cls.canonize(*args, **kwds)
            if isinstance(r, Basic):
                return r
            raise NotImplementedError('%s.__new__(%r,%r)' % (cls, args, kwds))
        else:
            r = cls.canonize(*args, **kwds)
            if isinstance(r, Basic): return r
            args, kwds = r
            cls.signature.validate(*args)
            return Composite._new(cls, *args, **kwds)

    @classmethod
    def canonize(cls, *args, **kwds):
        if cls is Function or cls is Functional:
            return FunctionClass(cls, *args, **kwds)
        return ((args,), kwds)

class Functional(Function):
    """
    Scalar-valued functions.
    """
    signature = FunctionSignature(None, (Basic,))

class sin(Functional):

    signature = FunctionSignature((Basic,), (Basic,))

class Add(Composite):

    pass
