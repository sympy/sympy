
from basic import Atom, CompositeTuple, Basic, BasicType
from methods import ArithMeths

class FunctionSignature:
    """
    Function signature defines valid function arguments
    and its expected return values.
    """

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


class FunctionClass(ArithMeths, Atom, BasicType):

    _new = type.__new__

    def torepr(cls):
        return cls.__name__

    @classmethod
    def canonize(cls, *args, **kwds):
        if not kwds and isinstance(args[0], type):
            basecls = args[0]
            name = args[1]
            d = basecls.__dict__.copy()
            d['undefined_Function'] = True
            if len(args)==3:
                signature = args[2]
                d['signature'] = signature
            return (name, (basecls,), d), {}
        return args, kwds

class Function(CompositeTuple):

    __metaclass__ = FunctionClass
    
    signature = FunctionSignature(None, None)

    @classmethod
    def canonize(cls, *args, **kwds):
        if cls is Function or cls is Functional:
            return FunctionClass(cls, *args, **kwds)
        args = tuple(map(cls.sympify, args))
        cls.signature.validate(args)
        return ((args,), kwds)

class Functional(ArithMeths, Function):
    """
    Scalar-valued functions.
    """
    signature = FunctionSignature(None, (Basic,))

class sin(Functional):

    signature = FunctionSignature((Basic,), (Basic,))
