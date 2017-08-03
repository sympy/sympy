from sympy.external import import_module
matchpy = import_module("matchpy")

if matchpy:
    Symbol = matchpy.Symbol
else:
    class Symbol(object):
        def __init__(self, name, variable_name=None):
            self.name = name
            self.head = self
        @staticmethod
        def dot(x):
            pass

class VariableSymbol(Symbol):
    pass

class matchpyInteger(Symbol):
    def __init__(self, value, variable_name=None):
        super(self.__class__, self).__init__(name=str(value), variable_name=variable_name)