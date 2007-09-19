
from utils import memoizer_immutable_args
from basic import Basic, MutableCompositeDict
from methods import ArithMeths, ImmutableMeths

class MutableAdd(ArithMeths, MutableCompositeDict):
    """ Represents a sum.

    3 + a + 2*b is Add({1:3, a:1, b:2})
    """
    
    def update(self, a):
        if isinstance(a, MutableAdd):
            for k,v in a.items():
                try:
                    self[k] += v
                except KeyError:
                    self[k] = v
            return
        if isinstance(a, dict) and not isinstance(a, Basic):
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

    # canonize methods
    def canonical(self):
        obj = self
        for k,v in obj.items():
            if v==0:
                # Add({a:0}) -> 0
                del obj[k]
        obj.__class__ = Add
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
            if v==1: return k
        return obj

    # arithmetics methods
    def __iadd__(self, other):
        self.update(other)
        return self



class Add(ImmutableMeths, MutableAdd):

    # constructor methods
    @memoizer_immutable_args
    def __new__(cls, *args, **options):
        return MutableAdd(*args, **options).canonical()

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
