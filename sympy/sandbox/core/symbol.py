
from basic import Atom

class Symbol(Atom, str):

    _new = str.__new__

    def __new__(cls, name):
        assert isinstance(name, str), `name`
        return cls._new(cls, name)

    def torepr(self):
        return '%s(%r)' % (self.__class__.__name__, str(self))

    def __hash__(self):
        try:
            return self._hash
        except AttributeError:
            h = hash((self.__class__.__name__, str(self)))
            self._hash = h
            return h
