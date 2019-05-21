from sympy import Symbol, Basic
from sympy.stats.rv import RandomSymbol

class IndexedRandomSymbol(RandomSymbol):

    def __new__(cls, symbol, type):
        if not isinstance(symbol, Symbol):
            raise TypeError("symbol should be of type Symbol")
        return Basic.__new__(cls, symbol, type)

    symbol = property(lambda self: self.args[0])
    name = property(lambda self: self.symbol.name)
    type = property(lambda self: self.args[1])
    is_random = True

    def __getattr__(self, attr):
        return getattr(self.type, attr)

    def _hashable_content(self):
        return self.type, self.symbol

    @property
    def free_symbols(self):
        return {self}

    def __getitem__(self, key):
        return self.type.__getitem__(key)
