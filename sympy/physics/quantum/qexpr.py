from sympy import Expr, sympify, Symbol
from sympy.printing.pretty.stringpict import prettyForm
from sympy.core.containers import Tuple

__all__ = [
    'QuantumError',
    'QExpr'
]


#-----------------------------------------------------------------------------
# Error handling
#-----------------------------------------------------------------------------

class QuantumError(Exception):
    pass

#-----------------------------------------------------------------------------
# Basic Quantum Expression from which all objects descend
#-----------------------------------------------------------------------------

class QExpr(Expr):
    """A base class for all quantum object like operators and states."""

    # In sympy, slots are for instance attributes that are computed
    # dynamically by the __new__ method. They are not part of args, but they
    # derive from args.

    # * The Hilbert space a quantum Object belongs to. It is an instance of
    # a HilbertSpace subclass.
    __slots__ = ['hilbert_space']

    # The separator used in printing the label
    _label_separator = ''

    def __new__(cls, label, **old_assumptions):
        """Construct a new quantum object.

        Parameters
        ==========
        label : tuple, sympy.core.containers.Tuple
            The list of numbers or parameters that uniquely specify the
            quantum object. For a state, this will be its symbol or its
            set of quantum numbers. The label does not include time for
            time dependent operators or states.

        Examples
        ========

        >>> from sympy.physics.quantum.qexpr import QExpr
        >>> q = QExpr(0)
        >>> q
        0
        >>> q.label
        Tuple(0)
        >>> q.hilbert_space
        H
        >>> q.is_commutative
        False
        """

        # First compute args and call Expr.__new__ to create the instance
        label = cls._eval_label(label)
        inst = Expr.__new__(cls, label, **{'commutative':False})
        # Now set the slots on the instance
        inst.hilbert_space = cls._eval_hilbert_space(label)
        return inst

    @classmethod
    def _new_rawargs(cls, hilbert_space, *args):
        """Create new instance of this class with hilbert_space and args.

        This is used to bypass the more complex logic in the ``__new__``
        method in cases where you already have the exact ``hilbert_space``
        and ``args``. This should be used when you are positive these
        arguments are valid, in their final, proper form and want to optimize
        the creation of the object.
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
        l = []
        for item in newlabel:
            if isinstance(item, basestring):
                i = Symbol(item)
            else:
                i = sympify(item)
            l.append(i)
        newlabel = Tuple(*l)
        return newlabel

    @classmethod
    def _eval_hilbert_space(cls, label):
        """Compute the Hilbert space instance from the label."""
        from sympy.physics.quantum.hilbert import HilbertSpace
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
        return '%s(%s)' % (
            self.__class__.__name__, self._print_contents_repr(printer, *args)
        )

    def _pretty(self, printer, *args):
        pform = self._print_contents_pretty(printer, *args)
        return pform

    def _latex(self, printer, *args):
        return self._print_contents(printer, *args)

    #-------------------------------------------------------------------------
    # Methods from Basic and Expr
    #-------------------------------------------------------------------------
    
    def doit(self, **kw_args):
        return self

    def _eval_rewrite(self, pattern, rule, **hints):
        # TODO: Make Basic.rewrite get the rule using the class name rather
        # than str(). See L1072 of basic.py.
        # This will call self.rule(*self.args) for rewriting.
        if hints.get('deep', False):
            args = [ a._eval_rewrite(pattern, rule, **hints) for a in self ]
        else:
            args = self.args

        if pattern is None or isinstance(self, pattern):
            if hasattr(self, rule):
                rewritten = getattr(self, rule)(*args)

                if rewritten is not None:
                    return rewritten

        return self

    #-------------------------------------------------------------------------
    # Represent
    #-------------------------------------------------------------------------

    def _represent(self, basis, **options):
        """Represent this object in a given basis.

        This method dispatches to the actual methods that perform the
        representation. Subclases of QExpr should define various methods to
        determine how the object will be represented in various bases. The
        format of these methods is::
        
            def _represent_BasisName(self, basis, **options):

        This to define how a quantum object is represented in the basis of
        the operator Position, you would define::

            def _represent_Position(self, basis, **options):

        Usually, Basis object will be instances of Operator subclasses, but
        there is a chance we will relax this in the future to accomodate other
        types of basis sets that are not associated with an operator.

        Parameters
        ==========
        basis : Operator
            The Operator whose basis functions will be used as the basis for
            representation.
        options : dict
            A dictionary of key/value pairs that give options and hints for
            the representation, such as the number of basis functions to
            be used.
        """
        return dispatch_method(self, '_represent', basis, **options)


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


def dispatch_method(self, basename, arg, **options):
    """Dispatch a method to the proper handlers."""
    method_name = '%s_%s' % (basename, arg.__class__.__name__)
    if hasattr(self, method_name):
        f = getattr(self, method_name)
        # This can raise and we will allow it to propagate.
        result = f(arg, **options)
        if result is not None:
            return result
    raise NotImplementedError(
        "%s.%s can't handle: %r" % \
            (self.__class__.__name__, basename, arg)
    )
