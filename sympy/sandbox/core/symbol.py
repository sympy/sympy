
from basic import Atom
from methods import ArithMeths, RelationalMeths

class Symbol(ArithMeths, RelationalMeths, Atom, str):

    def __new__(cls, name, **options):
        assert isinstance(name, str), `name`
        assert not options, `options`
        return str.__new__(cls, name)

    def torepr(self):
        return '%s(%r)' % (self.__class__.__name__, str(self))
