from matchpy import Symbol

class VariableSymbol(Symbol):
    pass

class ConstantSymbol(Symbol):
    def __init__(self, value):
        super(self.__class__, self).__init__(str(value))
        self.value = value
