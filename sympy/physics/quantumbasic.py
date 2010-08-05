from sympy import Expr
from sympy.core.decorators import call_highest_priority

#-----------------------------------------------------------------------------
# Error handling
#-----------------------------------------------------------------------------
class QuantumError(Exception):
    pass

#-----------------------------------------------------------------------------
# Basic Quantum Object from which all objects descend
#-----------------------------------------------------------------------------
class QuantumBasic(Expr):
    _op_priority = 100.0
    __slots__ = ['evaluates', 'hilbert_space']

    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        from sympy.physics.qmul import QMul
        return QMul(self, other)

    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        from sympy.physics.qmul import QMul
        return QMul(other, self)

    @call_highest_priority('__radd__')
    def __add__(self, other):
        from sympy.physics.qadd import QAdd    
        return QAdd(self, other)

    @call_highest_priority('__add__')
    def __radd__(self, other):
        from sympy.physics.qadd import QAdd       
        return QAdd(other, self)

    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        from sympy.physics.qpow import QPow
        return QPow(self, other)

    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        from sympy.physics.qpow import QPow
        return QPow(other, self)

