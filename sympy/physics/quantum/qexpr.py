
from sympy import Expr, sympify, Symbol, Matrix
from sympy.printing.pretty.stringpict import prettyForm
from sympy.core.containers import Tuple
from sympy.core.compatibility import is_sequence

from sympy.physics.quantum.matrixutils import (
    numpy_ndarray, scipy_sparse_matrix,
    to_sympy, to_numpy, to_scipy_sparse
)

__all__ = [
    'QuantumError',
    'QExpr'
]


#-----------------------------------------------------------------------------
# Error handling
#-----------------------------------------------------------------------------

class QuantumError(Exception):
    pass


def _qsympify_sequence(seq):
    """Convert elements of a sequence to standard form.

    This is like sympify, but it performs special logic for arguments passed
    to QExpr. The following conversions are done:

    * (list, tuple, Tuple) => _qsympify_sequence each element and convert
      sequence to a Tuple.
    * basestring => Symbol
    * Matrix => Matrix
    * other => sympify

    Strings are passed to Symbol, not sympify to make sure that variables like
    'pi' are kept as Symbols, not the Sympy built-in number subclasses.

    Examples
    ========

    >>> from sympy.physics.quantum.qexpr import _qsympify_sequence
    >>> _qsympify_sequence((1,2,[3,4,[1,]]))
    (1, 2, (3, 4, (1,)))

    """

    return tuple(__qsympify_sequence_helper(seq))

def __qsympify_sequence_helper(seq):
    """
       Helper function for _qsympify_sequence
       This function does the actual work.
    """
    #base case. If not a list, do Sympification
    if not is_sequence(seq):
        if isinstance(seq, Matrix):
            return seq
        elif isinstance(seq, basestring):
            return Symbol(seq)
        else:
            return sympify(seq)

    #if list, recurse on each item in the list
    result = [__qsympify_sequence_helper(item) for item in seq]

    return Tuple(*result)


#-----------------------------------------------------------------------------
# Basic Quantum Expression from which all objects descend
#-----------------------------------------------------------------------------

class QExpr(Expr):
    """A base class for all quantum object like operators and states."""

    # In sympy, slots are for instance attributes that are computed
    # dynamically by the __new__ method. They are not part of args, but they
    # derive from args.

    # The Hilbert space a quantum Object belongs to.
    __slots__ = ['hilbert_space']

    # The separator used in printing the label.
    _label_separator = u''

    def __new__(cls, *args, **old_assumptions):
        """Construct a new quantum object.

        Parameters
        ==========
        args : tuple
            The list of numbers or parameters that uniquely specify the
            quantum object. For a state, this will be its symbol or its
            set of quantum numbers.

        Examples
        ========

        >>> from sympy.physics.quantum.qexpr import QExpr
        >>> q = QExpr(0)
        >>> q
        0
        >>> q.label
        (0,)
        >>> q.hilbert_space
        H
        >>> q.args
        (0,)
        >>> q.is_commutative
        False
        """

        # First compute args and call Expr.__new__ to create the instance
        args = cls._eval_args(args)
        inst = Expr.__new__(cls, *args, **{'commutative':False})
        # Now set the slots on the instance
        inst.hilbert_space = cls._eval_hilbert_space(args)
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

        obj = Expr.__new__(cls, *args, **{'commutative':False})
        obj.hilbert_space = hilbert_space
        return obj

    #-------------------------------------------------------------------------
    # Properties
    #-------------------------------------------------------------------------

    @property
    def label(self):
        """The label is the unique set of identifiers for the object.

        Usually, this will include all of the information about the state
        *except* the time (in the case of time-dependent objects).

        This must be a tuple, rather than a Tuple.
        """
        return self.args

    @property
    def is_symbolic(self):
        return True

    #-------------------------------------------------------------------------
    # _eval_* methods
    #-------------------------------------------------------------------------

    @classmethod
    def _eval_args(cls, args):
        """Process the args passed to the __new__ method.

        This simply runs args through _qsympify_sequence.
        """
        return _qsympify_sequence(args)

    @classmethod
    def _eval_hilbert_space(cls, args):
        """Compute the Hilbert space instance from the args.
        """
        from sympy.physics.quantum.hilbert import HilbertSpace
        return HilbertSpace()

    def _eval_dagger(self):
        """Compute the Dagger of this state."""
        raise NotImplementedError('_eval_dagger not defined on: %r' % self)

    #-------------------------------------------------------------------------
    # Printing
    #-------------------------------------------------------------------------

    # Utilities for printing: these operate on raw sympy objects

    def _print_sequence(self, seq, sep, printer, *args):
        result = []
        for item in seq:
            result.append(printer._print(item, *args))
        return sep.join(result)

    def _print_sequence_pretty(self, seq, sep, printer, *args):
        pform = printer._print(seq[0], *args)
        for item in seq[1:]:
            pform = prettyForm(*pform.right((sep)))
            pform = prettyForm(*pform.right((printer._print(item, *args))))
        return pform

    # Utilities for printing: these operate prettyForm objects

    def _print_subscript_pretty(self, a, b):
        top = prettyForm(*b.left(' '*a.width()))
        bot = prettyForm(*a.right(' '*b.width()))
        return prettyForm(binding=prettyForm.POW, *bot.below(top))

    def _print_superscript_pretty(self, a, b):
        return a**b

    def _print_parens_pretty(self, pform, left='(', right=')'):
        return prettyForm(*pform.parens(left=left, right=right))

    # Printing of labels

    def _print_label(self, printer, *args):
        return self._print_sequence(
            self.label, self._label_separator, printer, *args
        )

    def _print_label_repr(self, printer, *args):
        return self._print_sequence(
            self.label, ',', printer, *args
        )

    def _print_label_pretty(self, printer, *args):
        return self._print_sequence_pretty(
            self.label, self._label_separator, printer, *args
        )

    def _print_label_latex(self, printer, *args):
        return self._print_sequence(
            self.label, self._label_separator, printer, *args
        )

    # Printing of contents

    def _print_contents(self, printer, *args):
        return self._print_label(printer, *args)

    def _print_contents_repr(self, printer, *args):
        return self._print_label_repr(printer, *args)

    def _print_contents_pretty(self, printer, *args):
        return self._print_label_pretty(printer, *args)

    def _print_contents_latex(self, printer, *args):
        return self._print_label_latex(printer, *args)

    # Main methods

    def _sympystr(self, printer, *args):
        return self._print_contents(printer, *args)

    def _sympyrepr(self, printer, *args):
        classname = self.__class__.__name__
        contents = self._print_contents_repr(printer, *args)
        return '%s(%s)' % (classname, contents)

    def _pretty(self, printer, *args):
        pform = self._print_contents_pretty(printer, *args)
        return pform

    def _latex(self, printer, *args):
        return self._print_contents_latex(printer, *args)

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

    def _represent_default_basis(self, **options):
        raise NotImplementedError('This object does not have a default basis')

    def _represent(self, **options):
        """Represent this object in a given basis.

        This method dispatches to the actual methods that perform the
        representation. Subclases of QExpr should define various methods to
        determine how the object will be represented in various bases. The
        format of these methods is::

            def _represent_BasisName(self, basis, **options):

        Thus to define how a quantum object is represented in the basis of
        the operator Position, you would define::

            def _represent_Position(self, basis, **options):

        Usually, basis object will be instances of Operator subclasses, but
        there is a chance we will relax this in the future to accomodate other
        types of basis sets that are not associated with an operator.

        If the ``format`` option is given it can be ("sympy", "numpy",
        "scipy.sparse"). This will ensure that any matrices that result from
        representing the object are returned in the appropriate matrix format.

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
        basis = options.pop('basis',None)
        if basis is None:
            result = self._represent_default_basis(**options)
        else:
            result = dispatch_method(self, '_represent', basis, **options)

        # If we get a matrix representation, convert it to the right format.
        format = options.get('format', 'sympy')
        if format == 'sympy' and not isinstance(result, Matrix):
            result = to_sympy(result)
        elif format == 'numpy' and not isinstance(result, numpy_ndarray):
            result = to_numpy(result)
        elif format == 'scipy.sparse' and\
            not isinstance(result, scipy_sparse_matrix):
            result = to_scipy_sparse(result)
        return result


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

