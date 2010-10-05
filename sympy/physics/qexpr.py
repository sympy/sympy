from sympy import Expr
from sympy.core.decorators import call_highest_priority
from sympy.core.basic import S, C

#-----------------------------------------------------------------------------
# Error handling
#-----------------------------------------------------------------------------

class QuantumError(Exception):
    pass

#-----------------------------------------------------------------------------
# Basic Quantum Expression from which all objects descend
#-----------------------------------------------------------------------------

class QExpr(Expr):

    # All quantum objects have a _op_priority that is higher than that of
    # Expr so that all operators are overloaded and the quantum versions
    # of Pow, Add, Mul, etc. are used.
    _op_priority = 100.0

    # In sympy, slots are for instance attributes that are computed
    # dynamically by the __new__ method. They are not part of args, but they
    # derive from args.
    # * 'acts_like' tells whether a binary operation acts like a Bra, Ket
    #   or Operator. This help us determine what types of subsequent
    #   operations are possible with that expression. This slot is set
    #   to the class that the object acts like.
    # * 'hilbert_space' tells us to which Hilbert space a quantum Object
    #   belongs. It is an instance of a HilbertSpace subclass.
    __slots__ = ['acts_like', 'hilbert_space']

    def __pos__(self):
        return self

    def __neg__(self):
        from sympy.physics.qmul import QMul
        return QMul(S.NegativeOne, self)

    def __abs__(self):
        return C.abs(self)

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

    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        from sympy.physics.qpow import QPow
        from sympy.physics.qmul import QMul
        return QMul(self, QPow(other, S.NegativeOne))

    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        from sympy.physics.qpow import QPow
        from sympy.physics.qmul import QMul
        return QMul(other, QPow(self, S.NegativeOne))

    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        from sympy.physics.qadd import QAdd
        return QAdd(self, -other)

    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        from sympy.physics.qadd import QAdd
        return QAdd(other, -self)

    __truediv__ = __div__

    __rtruediv__ = __rdiv__

