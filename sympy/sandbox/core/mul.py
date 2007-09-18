
from utils import memoizer_immutable_args
from basic import Basic, Composite
from methods import ArithMeths, ImmutableMeths

class MutableMul(ArithMeths, Composite, dict):
    """Mutable base class for Mul. This class is used temporarily
    during construction of Mul objects."""

    # constructor methods
    def __new__(cls, *args, **options):
        """
        To make MutableAdd immutable, execute
          obj.__class__ = Add
        """
        if 0 in args:
            return Basic.Integer(0)
        obj = dict.__new__(cls)
        [obj.update(a) for a in args]
        return obj.canonical()

    def __init__(self, *args, **options):
        pass

    def update(self, a):
        if isinstance(a, MutableMul):
            for k,v in a.items():
                if k.is_Add and len(k)==1:
                    # Add({x:3})**2 -> Mul({x:1*2,1:3**2})
                    k1,v1 = k.items()[0]
                    self.update(v1**v)
                    k = k1
                assert not k.is_Mul and not k.is_Add,`k`
                try:
                    self[k] += v
                except KeyError:
                    self[k] = v
            return
        if isinstance(a, dict) and not isinstance(a, Basic):
            # construct Mul instance from a canonical dictionary
            assert len(self)==0,`len(self)`
            super(MutableMul, self).update(a)
            return
        a = Basic.sympify(a)
        if a.is_Number:
            k, v = Basic.Integer(1), a
            try:
                self[k] *= v
            except KeyError:
                self[k] = v
        else:
            k, v = a, Basic.Integer(1)
            if k.is_Add and len(k)==1:
                # Add({x:3}) -> Mul({x:1,1:3})
                k1,v1 = k.items()[0]
                self.update(v1)
                k = k1
            if k.is_Mul:
                for k1,v1 in k.items():
                    assert not k1.is_Mul and not k1.is_Add,`k1`
                    try:
                        self[k1] += v1
                    except KeyError:
                        self[k1] = v1
                return
            assert not k.is_Mul and not k.is_Add,`k`
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

    # canonize methods
    def canonical(self):
        obj = self
        for k,v in obj.items():
            if v==0:
                # Mul({a:0}) -> 1
                del obj[k]
        c = obj.pop(1, Basic.Integer(1))
        obj.__class__ = Mul
        if len(obj)==0:
            return c
        if len(obj)==1:
            # Mul({a:1}) -> a
            k,v = obj.items()[0]
            if v==1:
                obj = k
        if c is not None:
            obj = Basic.MutableAdd({obj:c}).canonical()
        return obj

class Mul(ImmutableMeths, MutableMul):
    """Represents a product, with repeated factors collected into
    powers using a dict representation. The product a**m * b**n
    (where m and n may be arbitrary expressions and not just
    integers) is represented as Mul({a:m, b:n}).

    Note that purely rational multiples are counted using the Add
    class, so 3*x*y**2 --> Add({Mul({x:1, y:2}):3}).
    """
    # constructor methods
    @memoizer_immutable_args
    def __new__(cls, *args, **options):
        return MutableMul(*args, **options).canonical()

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
