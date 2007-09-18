
from utils import memoizer_immutable_args
from basic import Basic, Composite
from methods import ArithMeths, ImmutableMeths

class MutableMul(ArithMeths, Composite, dict):
    """ Represents a product.

    3 * a * b**2 is Mul({1:3, a:1, b:2})
    """

    # constructor methods
    def __new__(cls, *args, **options):
        """
        To make MutableAdd immutable, execute
          obj.__class__ = Add
        """
        obj = dict.__new__(cls)
        [obj.update(a) for a in args]
        return obj

    def __init__(self, *args, **options):
        pass

    def update(self, a):
        if isinstance(a, MutableMul):
            for k,v in a.items():
                try:
                    self[k] += v
                except KeyError:
                    self[k] = v
            return
        if isinstance(a, dict) and not isinstance(a, Basic):
            assert len(self)==0,`len(self)`
            super(MutableMul, self).update(a)
            return
        a = Basic.sympify(a)
        if isinstance(a, Basic.Number):
            k, v = Basic.Integer(1), a
            try:
                self[k] *= v
            except KeyError:
                self[k] = v
        else:
            k, v = a, Basic.Integer(1)
            try:
                self[k] += v
            except KeyError:
                self[k] = v

    # representation methods
    def torepr(self):
        return '%s(%s)' % (self.__class__.__name__, dict(self))

    # arithmetics methods
    def __imul__(self, other):
        self.update(other)
        return self

class Mul(ImmutableMeths, MutableMul):

    # constructor methods
    @memoizer_immutable_args
    def __new__(cls, *args, **options):
        obj = MutableMul(*args, **options)
        for k,v in obj.items():
            if v==0:
                # Mul({a:0}) -> 1
                del obj[k]
        c = obj.pop(1, Basic.Integer(1))
        if len(obj)==0:
            return c
        obj.__class__ = cls
        if len(obj)==1:
            # Mul({a:1}) -> a
            k,v = obj.items()[0]
            if v==1:
                obj = k
        if c is not None:
            obj = Basic.MutableAdd({obj:c})
            obj.__class__ = Basic.Add
        return obj

    # arithmetics methods
    def __imul__(self, other):
        return Mul(self, other)

    # object identity methods
    def __hash__(self):
        try:
            return self.__dict__['_cached_hash']
        except KeyError:
            h = self._cached_hash = sum(map(hash, self.items()))
        return h

class Pow(Basic):

    def __new__(cls, a, b):
        a = Basic.sympify(a)
        b = Basic.sympify(b)
        if b==0: return Integer(1)
        if b==1: return b
        p = a._eval_power(b)
        if p is not None: return p
        m = MutableMul({a:b})
        m.__class__ = Mul
        return m
