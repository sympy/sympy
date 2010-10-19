from sympy import Expr, sympify
from sympy.printing.pretty.stringpict import prettyForm
from sympy.core.containers import Tuple


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

    _label_separator = ''

    def __new__(cls, label, **old_assumptions):
        # First compute args and call Expr.__new__ to create the instance
        label = cls._eval_label(label)
        inst = Expr.__new__(cls, label, **{'commutative':False})
        # Now set the slots on the instance
        inst.hilbert_space = cls._eval_hilbert_space(label)
        return inst

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

    #-------------------------------------------------------------------------
    # Properties
    #-------------------------------------------------------------------------

    @property
    def label(self):
        """The label is the unique set of identifiers for the state.

        The label of a state is what distinguishes it from other states. For
        eigenstates, the label is usually a sympy.core.containers.Tuple of the
        quantum numbers. For an abstract state, it would just be a single
        element Tuple of the symbol of the state.
        """
        return self.args[0]

    @property
    def is_symbolic(self):
        return True

    #-------------------------------------------------------------------------
    # _eval_* methods
    #-------------------------------------------------------------------------

    @classmethod
    def _eval_label(cls, label):
        """Make sure that label is a sympy.core.containers.Tuple.

        The elements of the Tuple must be run through sympify.
        """
        if not isinstance(label, Tuple):
            if isinstance(label, (list, tuple)):
                # Convert a raw tuple or list to a Tuple
                newlabel = Tuple(*label)
            else:
                # Single element label gets wrapped into a tuple.
                newlabel = Tuple(label)
        else:
            newlabel = label
        newlabel = Tuple(*[sympify(item) for item in newlabel])
        return newlabel

    @classmethod
    def _eval_hilbert_space(cls, label):
        """Compute the Hilbert space instance from the label."""
        from sympy.physics.hilbert import HilbertSpace
        return HilbertSpace()

    def _eval_dagger(self):
        """Compute the Dagger of this state."""
        raise NotImplementedError('_eval_dagger not defined on: %r' % self)

    #-------------------------------------------------------------------------
    # Printing
    #-------------------------------------------------------------------------

    def _print_label(self, printer, *args):
        result = []
        for item in self.label:
            result.append(printer._print(item, *args))
        return self._label_separator.join(result)

    def _print_label_repr(self, printer, *args):
        return printer._print(self.label, *args)

    def _print_label_pretty(self, printer, *args):
        pform = printer._print(self.label[0], *args)
        for item in self.label[1:]:
            pform = prettyForm(*pform.right((self._label_separator)))
            nextpform = printer._print(item, *args)
            pform = prettyForm(*pform.right((nextpform)))
        return pform

    def _print_contents(self, printer, *args):
        return self._print_label(printer, *args)

    def _print_contents_repr(self, printer, *args):
        return self._print_label_repr(printer, *args)

    def _print_contents_pretty(self, printer, *args):
        return self._print_label_pretty(printer, *args)

    def _sympystr(self, printer, *args):
        return self._print_contents(printer, *args)

    def _sympyrepr(self, printer, *args):
        return self._print_contents_repr(printer, *args)

    def _pretty(self, printer, *args):
        pform = self._print_contents_pretty(printer, *args)
        return pform

    #-------------------------------------------------------------------------
    # Methods from Basic and Expr
    #-------------------------------------------------------------------------
    
    def doit(self, **kw_args):
        return self

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
