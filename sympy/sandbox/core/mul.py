
from utils import memoizer_immutable_args
from basic import Basic, MutableCompositeDict
from methods import ArithMeths, ImmutableMeths

class MutableMul(ArithMeths, MutableCompositeDict):
    """Mutable base class for Mul. This class is used temporarily
    during construction of Mul objects."""

    # canonize methods

    def update(self, a, p=1):
        """
        Mul({}).update(a,p) -> Mul({a:p})
        """
        if a.__class__ is dict:
            # construct Mul instance from a canonical dictionary
            assert len(self)==0,`len(self)` # make sure no data is overwritten
            assert p is 1,`p`
            super(MutableMul, self).update(a)
            return
        a = Basic.sympify(a)
        if a.is_Number:
            if p==1: v = a
            else: v = a ** p
            try:
                self[1] *= v
            except KeyError:
                self[1] = v
        elif a.is_Add and len(a)==1:
            # Mul({x:3,1:4}).update(Add({x:2})) -> Mul({x:3+1,1:4*2})
            k, v = a.items()[0]
            self.update(k, p)
            self.update(v, p)
            return
        elif a.is_MutableMul:
            # Mul({x:3}).update(Mul({x:2}), 4) -> Mul({x:3}).update(x,2*4)
            for k,v in a.items():
                # todo?: make it noncommutative product for (a**2)**(1/2)
                #        (a**z)**w where z,w are complex numbers
                self.update(k, v * p) 
        else:
            try:
                self[a] += p
            except KeyError:
                self[a] = p

    def canonical(self):
        # self will be modified in-place,
        # always return an immutable object
        obj = self
        c = obj.pop(1, Basic.Integer(1))
        for k,v in obj.items():
            if v==0:
                # Mul({a:0}) -> 1
                del obj[k]
        if c==0:
            # todo: handle 0*oo->nan, either here or in Number
            return c
        if len(obj)==0:
            return c
        obj.__class__ = Mul
        if len(obj)==1:
            # Mul({a:1}) -> a
            k,v = obj.items()[0]
            if v==1:
                obj = k
        if c!=1:
            # Mul({1:c,rest:power}) -> Add({Mul({rest:power}):c})
            obj = Basic.Add({obj:c})
        return obj

    # arithmetics methods
    def __imul__(self, other):
        self.update(other)
        return self



class Mul(ImmutableMeths, MutableMul):
    """Represents a product, with repeated factors collected into
    powers using a dict representation. The product a**m * b**n
    (where m and n may be arbitrary expressions and not just
    integers) is represented as Mul({a:m, b:n}).

    Note that purely rational multiples are counted using the Add
    class, so 3*x*y**2 --> Add({Mul({x:1, y:2}):3}).
    """
    # constructor methods
    @memoizer_immutable_args("Mul.__new__")
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

    def canonical(self):
        return self

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
