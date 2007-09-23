
from utils import memoizer_immutable_args, memoizer_Symbol_new
from basic import Atom
from methods import ArithMeths, RelationalMeths

class Symbol(ArithMeths, RelationalMeths, Atom, str):

    """ Represents a symbol.

    Symbol('x', dummy=True) returns a unique Symbol instance.
    """

    _dummy_count = 0
    is_dummy = False
    
    @memoizer_Symbol_new
    def __new__(cls, name, dummy=False, **options):
        # when changing the Symbol signature then change memoizer_Symbol_new
        # accordingly
        assert isinstance(name, str), `name`
        obj = str.__new__(cls, name)
        obj.is_dummy = dummy
        if dummy:
            Symbol._dummy_count += 1
            obj.dummy_index = Symbol._dummy_count
        return obj

    def torepr(self):
        if self.is_dummy:
            return '%s(%r, dummy=True)' % (self.__class__.__name__, str(self))
        return '%s(%r)' % (self.__class__.__name__, str(self))

    def tostr(self):
        return str(self)

    def compare(self, other):
        if self is other: return 0
        c = cmp(self.__class__, other.__class__)
        if c: return c
        if self.is_dummy or other.is_dummy:
            return cmp(id(self), id(other))
        return cmp(str(self), str(other))

    
