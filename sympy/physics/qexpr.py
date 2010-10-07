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
    # * 'acts_like' tells whether a binary operation acts like a BraBase,
    #   KetBase or Operator. These are the only possibilities. This help us
    #   determine what types of subsequent operations are possible with that
    #   expression. This slot is set to the class that the object acts like.
    # * 'hilbert_space' tells us to which Hilbert space a quantum Object
    #   belongs. It is an instance of a HilbertSpace subclass.
    __slots__ = ['acts_like', 'hilbert_space']

    @classmethod
    def _new_rawargs(cls, acts_like, hilbert_space, *args):
        """Create new instance of own class with args exactly as provided.

        This class method sets the ``acts_like`` and ``hilbert_space`` slots
        to those specified by the caller, rather than computing them
        dynamically.

        This is handy when we want to optimize things by not having to do
        a careful application of the rules and dynamically computing the
         ``acts_like`` and ``hilbert_space`` slots. If we
        know what it evaluates to and to what hilbert_space it will belong,
        one can hardcode that in. Thus, It's a quicker way of instantiating a
        QAssocOp If we know what type we will get. This was stolen from
        sympy.core.operations.

            >>> from sympy.physics.hilbert import HilbertSpace
            >>> from sympy.physics.quantum import Ket, Operator
            >>> from sympy.physics.qmul import QMul
            >>> QMul._new_rawargs(Ket, HilbertSpace(), Operator('psi'), Ket('a'))
            psi*|a>
            >>> a = _
            >>> a.acts_like
            <class 'sympy.physics.quantum.Ket'>
            >>> a.hilbert_space
            HilbertSpace()
        """

        if len(args) == 1:
            return args[0]
        obj = Expr.__new__(cls, *args)
        obj.acts_like = acts_like
        obj.hilbert_space = hilbert_space
        return obj

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

