from sympy.external import import_module
matchpy = import_module("matchpy")

if matchpy:
    Symbol, Wildcard = matchpy.Symbol, matchpy.Wildcard
else:
    raise ImportError('MatchPy could not be imported')

class VariableSymbol(Symbol):
    pass

class Integer(Symbol):
    def __init__(self, value, variable_name=None):
        super(self.__class__, self).__init__(name=str(value), variable_name=variable_name)
