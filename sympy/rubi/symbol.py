from matchpy import Symbol, Wildcard
import matchpy

class VariableSymbol(Symbol):
    pass

'''
class ConstantSymbol(Symbol):
    def __init__(self, value):
        super(self.__class__, self).__init__(name=str(value))
        self.value = value
'''

class ConstantSymbol(Symbol):
    def __init__(self, value, variable_name=None):
        super(self.__class__, self).__init__(name=str(value), variable_name=variable_name)
        self.value = value
        
