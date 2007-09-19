
from utils import memoizer_immutable_args
from basic import Basic, MutableCompositeDict
from methods import ArithMeths, ImmutableMeths

class MutableAdd(ArithMeths, MutableCompositeDict):
    """ Represents a sum.

    3 + a + 2*b is Add({1:3, a:1, b:2})
    """

    # canonize methods

    def update(self, a, p=1):
        """
        Add({}).update(a,p) -> Add({a:p})
        """
        if a.__class__ is dict:
            # construct Add instance from a canonical dictionary
            assert len(self)==0,`len(self)` # make sure no data is overwritten
            super(MutableAdd, self).update(a)
            return
        a = Basic.sympify(a)
        if a.is_Number:
            if p==1: v = a
            else: v = a * p
            try:
                self[1] += v
            except KeyError:
                self[1] = v
        elif a.is_Add:
            # (3*x) + ((2*x)*4) -> (3*x) + (8*x)
            # Add({x:3}).update(Add({x:2}), 4) -> Add({x:3}).update(x,2*4)
            for k,v in a.items():
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
        for k,v in obj.items():
            if v==0:
                # Add({a:0}) -> 0
                del obj[k]
        if len(obj)==0:
            return Basic.Integer(0)
        if len(obj)==1:
            try:
                # Add({1:3}) -> 3
                return obj[1]
            except KeyError:
                pass
            # Add({a:1}) -> a
            k,v = obj.items()[0]
            if v==1:
                return k
        obj.__class__ = Add
        return obj

    # arithmetics methods
    def __iadd__(self, other):
        self.update(other)
        return self

class Add(ImmutableMeths, MutableAdd):

    # constructor methods
    @memoizer_immutable_args("Add.__new__")
    def __new__(cls, *args, **options):
        return MutableAdd(*args, **options).canonical()

    # canonize methods
    def canonical(self):
        return self

    # arithmetics methods
    def __iadd__(self, other):
        return Add(self, other)

    # object identity methods
    def __hash__(self):
        try:
            return self.__dict__['_cached_hash']
        except KeyError:
            h = self._cached_hash = sum(map(hash, self.items()))
        return h


