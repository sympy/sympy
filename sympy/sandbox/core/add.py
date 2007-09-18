
from utils import memoizer_immutable_args
from basic import Basic, Composite
from methods import ArithMeths

class MutableAdd(ArithMeths, Composite, dict):

    def __new__(cls, *args, **options):
        """
        To make MutableAdd immutable, execute
          obj.__class__ = Add
        """
        obj = dict.__new__(MutableAdd)
        [obj.update(a) for a in args]
        obj.__class__ = cls
        return obj

    def __init__(self, *args, **options):
        pass

    def update(self, a):
        if isinstance(a, MutableAdd):
            for k,v in a.items():
                try:
                    self[k] += v
                except KeyError:
                    self[k] = v
            return
        if isinstance(a, dict):
            assert len(self)==0,`len(self)`
            super(MutableAdd, self).update(a)
            return
        a = Basic.sympify(a)
        if isinstance(a, Basic.Number):
            k, v = Basic.Integer(1), a
        else:
            k, v = a, Basic.Integer(1)
        try:
            self[k] += v
        except KeyError:
            self[k] = v

    def torepr(self):
        return '%s(%s)' % (self.__class__.__name__, dict(self))

    def __iadd__(self, other):
        self.update(other)
        return self


class Add(MutableAdd):

    @memoizer_immutable_args
    def __new__(cls, *args, **options):
        obj = dict.__new__(MutableAdd)
        [obj.update(a) for a in args]
        obj.__class__ = cls
        return obj    

    def __setitem__(self, k, v):
        raise TypeError('%s instance is immutable' % (self.__class__.__name__))

    def __delitem__(self, k, v):
        raise TypeError('%s instance is immutable' % (self.__class__.__name__))

    def popitem(self):
        raise TypeError('%s instance is immutable' % (self.__class__.__name__))

    def __iadd__(self, other):
        return Add(self, other)
