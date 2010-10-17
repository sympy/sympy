from sympy import Expr

#-----------------------------------------------------------------------------
# Error handling
#-----------------------------------------------------------------------------

class QuantumError(Exception):
    pass

#-----------------------------------------------------------------------------
# Basic Quantum Expression from which all objects descend
#-----------------------------------------------------------------------------

class QExpr(Expr):

    # In sympy, slots are for instance attributes that are computed
    # dynamically by the __new__ method. They are not part of args, but they
    # derive from args.
    # * 'hilbert_space' tells us to which Hilbert space a quantum Object
    #   belongs. It is an instance of a HilbertSpace subclass.
    __slots__ = ['hilbert_space']

    @classmethod
    def _new_rawargs(cls, hilbert_space, *args):
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

        obj = Expr.__new__(cls, *args)
        obj.hilbert_space = hilbert_space
        return obj


def split_commutative_parts(e):
    """Split into commutative and non-commutative parts."""
    c_part = [p for p in e.args if p.is_commutative]
    nc_part = [p for p in e.args if not p.is_commutative]
    return c_part, nc_part


def split_qexpr_parts(e):
    """Split an expression into Expr and noncommutative QExpr parts."""
    expr_part = []
    qexpr_part = []
    for arg in e.args:
        if not isinstance(arg, QExpr):
            expr_part.append(arg)
        else:
            qexpr_part.append(arg)
    return expr_part, qexpr_part
