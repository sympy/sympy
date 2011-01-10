"""Dirac notation for Bras, Kets, Operators, and so on.

TODO:

* Fix early 0 in apply_operators.
* Debug and test apply_operators.
* Get cse working with classes in this file.
* Doctests and documentation of special methods for InnerProduct, Commutator,
  AntiCommutator, represent, apply_operators.
* Decide how to handle the label of Operators. Currently, if there is only
  1 element in label, we treat it differently.
"""

import copy

from sympy import S, Expr, sympify, Add, Mul, Function, Matrix, Pow, Integer
from sympy.printing.pretty.stringpict import prettyForm, stringPict
from sympy.core.numbers import NumberSymbol
import sympy.mpmath.libmp as mlib
from sympy.utilities.iterables import all

from sympy.physics.qexpr import (
    QuantumError, QExpr, split_commutative_parts, dispatch_method
)

__all__ = [
    'KetBase',
    'BraBase',
    'StateBase',
    'State',
    'Ket',
    'Bra',
    'TimeDepState',
    'TimeDepBra',
    'TimeDepKet',
    'Operator',
    'HermitianOperator',
    'UnitaryOperator',
    'OuterProduct',
    'InnerProduct',
    'Dagger',
    'KroneckerDelta',
    'Commutator',
    'AntiCommutator',
    'represent',
    'hbar',
    'TensorProduct',
    'tensor_product_simp',
    'apply_operators',
    'apply_operators_Mul',
    'apply_single_op',
    'apply_list_of_ops'
]


#-----------------------------------------------------------------------------
# States, bras and kets.
#-----------------------------------------------------------------------------

class StateBase(QExpr):
    """Abstract base class for general abstract states in quantum mechanics.

    All other state classes defined will need to inherit from this class. It
    carries the basic structure for all other states such as dual, _eval_dagger
    and label.

    This is an abstract base class and you should not instantiate it directly,
    instead use State.
    """

    is_continuous = False
    is_discrete = False

    #-------------------------------------------------------------------------
    # _eval_* methods
    #-------------------------------------------------------------------------

    def _eval_innerproduct(self, other, **hints):
        return dispatch_method(self, '_eval_innerproduct', other, **hints)

    #-------------------------------------------------------------------------
    # Dagger/dual
    #-------------------------------------------------------------------------

    @property
    def dual(self):
        """Return the dual state of this one."""
        return self.dual_class._new_rawargs(self.hilbert_space, *self.args)

    @property
    def dual_class(self):
        """Return the class used to construt the dual."""
        raise NotImplementedError(
            'dual_class must be implemented in a subclass'
        )

    def _eval_dagger(self):
        """Compute the dagger of this state using the dual."""
        return self.dual

    #-------------------------------------------------------------------------
    # Printing
    #-------------------------------------------------------------------------

    def _sympystr(self, printer, *args):
        return '%s%s%s' % (self.lbracket, self._print_contents(printer, *args),
                           self.rbracket)

    def _pretty(self, printer, *args):
        from sympy.printing.pretty.stringpict import prettyForm
        pform = self._print_contents_pretty(printer, *args)
        pform = prettyForm(*pform.left((self.lbracket_pretty)))
        pform = prettyForm(*pform.right((self.rbracket_pretty)))
        return pform

    def _latex(self, printer, *args):
        contents = self._print_contents(printer, *args)
        return '%s%s%s' % (self.lbracket_latex, contents, self.rbracket_latex)


class KetBase(StateBase):
    """Base class for Kets.

    This class defines the dual property and the brackets for printing. This
    is an abstract base class and you should not instantiate it directly,
    instead use Ket.
    """

    lbracket = '|'
    rbracket = '>'
    lbracket_pretty = prettyForm(u'\u2758')
    rbracket_pretty = prettyForm(u'\u27E9')
    lbracket_latex = r'\left|'
    rbracket_latex = r'\right\rangle '

    @property
    def dual_class(self):
        return BraBase

    def __mul__(self, other):
        """KetBase*other"""
        if isinstance(other, BraBase):
            return OuterProduct(self, other)
        else:
            return Expr.__mul__(self, other)

    def __rmul__(self, other):
        """other*KetBase"""
        if isinstance(other, BraBase):
            return InnerProduct(other, self)
        else:
            return Expr.__rmul__(self, other)

    def _apply_operator(self, op, **options):
        """Apply an Operator to this Ket.

        This method will dispatch to methods having the format::

            def _apply_operator_OperatorName(op, **options):

        Subclasses should define these methods (one for each OperatorName) to
        teach the Ket how operators act on it.

        Parameters
        ==========
        op : Operator
            The Operator that is acting on the Ket.
        options : dict
            A dict of key/value pairs that control how the operator is applied
            to the Ket.
        """
        return dispatch_method(self, '_apply_operator', op, **options)


class BraBase(StateBase):
    """Base class for Bras.

    This class defines the dual property and the brackets for printing. This
    is an abstract base class and you should not instantiate it directly,
    instead use Bra.
    """

    lbracket = '<'
    rbracket = '|'
    lbracket_pretty = prettyForm(u'\u27E8')
    rbracket_pretty = prettyForm(u'\u2758')
    lbracket_latex = r'\left\langle '
    rbracket_latex = r'\right|'

    @property
    def dual_class(self):
        return KetBase

    def __mul__(self, other):
        """BraBase*other"""
        if isinstance(other, KetBase):
            return InnerProduct(self, other)
        else:
            return Expr.__mul__(self, other)

    def __rmul__(self, other):
        """other*BraBase"""
        if isinstance(other, KetBase):
            return OuterProduct(other, self)
        else:
            return Expr.__rmul__(self, other)


class State(StateBase):
    """General abstract quantum state used as a base class for Ket and Bra."""
    pass



class Ket(State, KetBase):
    """A general time-independent Ket in quantum mechanics.

    Inherits from State and KetBase. This class should be used as the base
    class for all physical, time-independent Kets in a system. This class
    and its subclasses will be the main classes that users will use for 
    expressing Kets in Dirac notation.

    Parameters
    ==========
    label : tuple, sympy.core.containers.Tuple
        The list of numbers or parameters that uniquely specify the
        ket. This will usually be its symbol or its quantum numbers.

    Examples
    ========

    Create a simple Ket and looking at its properties::

        >>> from sympy.physics.quantum import Ket, Bra
        >>> from sympy import symbols, I
        >>> k = Ket('psi')
        >>> k
        |psi>
        >>> k.hilbert_space
        H
        >>> k.is_commutative
        False
        >>> k.label
        Tuple(psi)

    Ket's know about their associated bra::

        >>> k.dual
        <psi|
        >>> k.dual_class
        <class 'sympy.physics.quantum.Bra'>

    Take a linear combination of two kets::

        >>> k0 = Ket(0)
        >>> k1 = Ket(1)
        >>> 2*I*k0 - 4*k1
        -4*|1> + 2*I*|0>

    Compound labels are passed as tuples::

        >>> n, m = symbols('nm')
        >>> k = Ket((n,m))
        >>> k
        |nm>

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Bra-ket_notation
    """

    @property
    def dual_class(self):
        return Bra


class Bra(State, BraBase):
    """A general time-independent Bra in quantum mechanics.

    Inherits from State and BraBase. A Bra is the dual of a Ket [1]. This
    class and its subclasses will be the main classes that users will use for 
    expressing Bras in Dirac notation.

    Parameters
    ==========
    label : tuple, sympy.core.containers.Tuple
        The list of numbers or parameters that uniquely specify the
        ket. This will usually be its symbol or its quantum numbers.

    Examples
    ========

    Create a simple Bra and looking at its properties::

        >>> from sympy.physics.quantum import Ket, Bra
        >>> from sympy import symbols, I
        >>> b = Bra('psi')
        >>> b
        <psi|
        >>> b.hilbert_space
        H
        >>> b.is_commutative
        False

    Bra's know about their dual Ket's::

        >>> b.dual
        |psi>
        >>> b.dual_class
        <class 'sympy.physics.quantum.Ket'>

    Like Ket's Bras can have compound labels and be manipulated::

        >>> n, m = symbols('nm')
        >>> b = Bra((n,m)) - I*Bra((m,n))
        >>> b
        -I*<mn| + <nm|

    Symbols ina Bra can be substituted using ``.subs``::

        >>> b.subs(n,m)
        -I*<mm| + <mm|

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Bra-ket_notation
    """

    @property
    def dual_class(self):
        return Ket

#-----------------------------------------------------------------------------
# Time dependent states, bras and kets.
#-----------------------------------------------------------------------------

class TimeDepState(StateBase):
    """Base class for a general time-dependent quantum state.

    This class is used as a base class for any time-dependent state. The main
    difference between this class and the time-independent state is that this
    class takes a second argument that is the time in addition to the usual
    label argument.

    Parameters
    ==========
    label : tuple, sympy.core.containers.Tuple
        The list of numbers or parameters that uniquely specify the
        ket. This will usually be its symbol or its quantum numbers.
    time : float
        The time of the state.
    """

    def __new__(cls, label, time, **old_assumptions):
        # First compute args and call Expr.__new__ to create the instance
        label = cls._eval_label(label)
        time = cls._eval_time(time)
        inst = Expr.__new__(cls, label, time, **{'commutative':False})
        # Now set the slots on the instance
        inst.hilbert_space = cls._eval_hilbert_space(label)
        return inst

    @property
    def time(self):
        """The time of the state."""
        return self.args[1]

    @classmethod
    def _eval_time(cls, time):
        return sympify(time)

    def _print_time(self, printer, *args):
        return printer._print(self.time, *args)

    def _print_time_repr(self, printer, *args):
        return printer._print(self.time, *args)

    def _print_time_pretty(self, printer, *args):
        pform = printer._print(self.time, *args)
        return pform

    def _print_contents(self, printer, *args):
        label = self._print_label(printer, *args)
        time = self._print_time(printer, *args)
        return '%s;%s' % (label, time)

    def _print_contents_repr(self, printer, *args):
        label = self._print_label_repr(printer, *args)
        time = self._print_time_repr(printer, *args)
        return '%s,%s' % (label, time)

    def _print_contents_pretty(self, printer, *args):
        pform = self._print_label_pretty(printer, *args)
        pform = prettyForm(*pform.right((';')))
        nextpform = self._print_time_pretty(printer, *args)
        pform = prettyForm(*pform.right((nextpform)))
        return pform


class TimeDepKet(TimeDepState, KetBase):
    """General time-dependent Ket in quantum mechanics.

    This inherits from TimeDepState and KetBase and is the main class that
    should be used for Kets that vary with time. Its dual is a TimeDepBra.

    Parameters
    ==========
    label : tuple, sympy.core.containers.Tuple
        The list of numbers or parameters that uniquely specify the
        ket. This will usually be its symbol or its quantum numbers.
    time : float
        The time of the state.

    Examples
    ========

    Create a TimeDepKet and look at its attributes::

        >>> from sympy.physics.quantum import TimeDepKet
        >>> k = TimeDepKet('psi', 't')
        >>> k
        |psi;t>
        >>> k.time
        t
        >>> k.label
        Tuple(psi)
        >>> k.hilbert_space
        H

    TimeDepKets know about their dual bra::

        >>> k.dual
        <psi;t|
        >>> k.dual_class
        <class 'sympy.physics.quantum.TimeDepBra'>
    """

    @property
    def dual_class(self):
        return TimeDepBra


class TimeDepBra(TimeDepState, BraBase):
    """General time-dependent Bra in quantum mechanics.

    This inherits from TimeDepState and BraBase and is the main class that
    should be used for Bras that vary with time. Its dual is a TimeDepBra.

    Parameters
    ==========
    label : tuple, sympy.core.containers.Tuple
        The list of numbers or parameters that uniquely specify the
        ket. This will usually be its symbol or its quantum numbers.
    time : float
        The time of the state.

    Examples
    ========

        >>> from sympy.physics.quantum import TimeDepBra
        >>> from sympy import symbols, I
        >>> b = TimeDepBra('psi', 't')
        >>> b
        <psi;t|
        >>> b.time
        t
        >>> b.label
        Tuple(psi)
        >>> b.hilbert_space
        H
        >>> b.dual
        |psi;t>
    """

    @property
    def dual_class(self):
        return TimeDepKet

#-----------------------------------------------------------------------------
# Operators and outer products
#-----------------------------------------------------------------------------

class Operator(QExpr):
    """Base class for non-commuting quantum operators.

    An operator maps one ket to another [1]. In quantum mechanics, Hermitian
    operators correspond to observables [2].

    Parameters
    ==========
    label : tuple, sympy.core.containers.Tuple
        The list of numbers or parameters that uniquely specify the
        operator.

    Examples
    ========

    Create an operator and examine its attributes::

        >>> from sympy.physics.quantum import Operator
        >>> from sympy import symbols, I
        >>> A = Operator('A')
        >>> A
        A
        >>> A.hilbert_space
        H
        >>> A.label
        Tuple(A)
        >>> A.is_commutative
        False

    Create another operator and do some arithmetic operations::

        >>> B = Operator('B')
        >>> C = 2*A*A + I*B
        >>> C
        I*B + 2*A**2

    Operators don't commute::

        >>> A.is_commutative
        False
        >>> B.is_commutative
        False
        >>> A*B == B*A
        False

    Polymonials of operators respect the commutation properties::

        >>> e = (A+B)**3
        >>> e.expand()
        A**2*B + B**2*A + A*B**2 + B*A**2 + A**3 + B**3 + A*B*A + B*A*B

    Operator inverses are handle symbolically::

        >>> A.inv()
        1/A
        >>> A*A.inv()
        1

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Operator
    [2] http://en.wikipedia.org/wiki/Observable
    """

    #-------------------------------------------------------------------------
    # Printing
    #-------------------------------------------------------------------------

    _label_separator = ','

    def _print_operator_name(self, printer, *args):
        return printer._print(self.__class__.__name__, *args)

    def _print_operator_name_pretty(self, printer, *args):
        return prettyForm(self.__class__.__name__)

    def _print_contents(self, printer, *args):
        if len(self.label) == 1:
            return self._print_label(printer, *args)
        else:
            return '%s(%s)' % (
                self._print_operator_name(printer, *args),
                self._print_label(printer, *args)
            )

    def _print_contents_pretty(self, printer, *args):
        if len(self.label) == 1:
            return self._print_label_pretty(printer, *args)
        else:
            pform = self._print_operator_name_pretty(printer, *args)
            label_pform = self._print_label_pretty(printer, *args)
            label_pform = prettyForm(
                *label_pform.parens(left='(', right=')')
            )
            pform = prettyForm(*pform.right((label_pform)))
            return pform

    #-------------------------------------------------------------------------
    # _eval_* methods
    #-------------------------------------------------------------------------

    def _eval_commutator(self, other, **options):
        """Evaluate [self, other] if known, return None if not known."""
        return dispatch_method(self, '_eval_commutator', other, **options)

    def _eval_anticommutator(self, other, **options):
        """Evaluate [self, other] if known."""
        return dispatch_method(self, '_eval_anticommutator', other, **options)

    #-------------------------------------------------------------------------
    # Operator application
    #-------------------------------------------------------------------------

    def _apply_operator(self, ket, **options):
        if not isinstance(ket, KetBase):
            raise TypeError('KetBase expected, got: %r' % ket)
        return dispatch_method(self, '_apply_operator', ket, **options)

    def matrix_element(self, *args):
        raise NotImplementedError('matrix_elements is not defined')

    #-------------------------------------------------------------------------
    # Printing
    #-------------------------------------------------------------------------

    def inverse(self):
        return self._eval_inverse()

    inv = inverse

    def _eval_inverse(self):
        # TODO: make non-commutative Exprs print powers using A**-1, not 1/A.
        return self**(-1)


class HermitianOperator(Operator):
    """A Hermitian operator that satisfies H == Dagger(H).

    Parameters
    ==========
    label : tuple, sympy.core.containers.Tuple
        The list of numbers or parameters that uniquely specify the
        operator.

    Examples
    ========

    >>> from sympy.physics.quantum import Dagger, HermitianOperator
    >>> H = HermitianOperator('H')
    >>> Dagger(H)
    H
    """

    def _eval_dagger(self):
        return self


class UnitaryOperator(Operator):
    """A unitary operator that satisfies U*Dagger(U) == 1.

    Parameters
    ==========
    label : tuple, sympy.core.containers.Tuple
        The list of numbers or parameters that uniquely specify the
        operator.

    Examples
    ========

    >>> from sympy.physics.quantum import Dagger, UnitaryOperator
    >>> U = UnitaryOperator('U')
    >>> U*Dagger(U)
    1
    """

    def _eval_dagger(self):
        return self._eval_inverse()


class OuterProduct(Operator):
    """An unevaluated outer product between a ket and kra.

    This constructs an outer product between any subclass of KetBase and
    BraBase as |a><b|. An OuterProduct inherits from Operator as they act as
    operators in quantum expressions.  For reference see [1].

    Parameters
    ==========
    ket : KetBase or subclass
        The ket on the left side of the outer product.
    bar : BraBase or subclass
        The bra on the right side of the outer product.

    Examples
    ========

    Create a simple outer product by hand and take its dagger::

        >>> from sympy.physics.quantum import Ket, Bra, OuterProduct, Dagger
        >>> from sympy.physics.quantum import Operator

        >>> k = Ket('k')
        >>> b = Bra('b')
        >>> op = OuterProduct(k, b)
        >>> op
        |k><b|
        >>> op.hilbert_space
        H
        >>> op.ket
        |k>
        >>> op.bra
        <b|
        >>> Dagger(op)
        |b><k|

    In simple products of kets and bras outer products will be automatically
    identified and created::

        >>> k*b
        |k><b|

    But in more complex expressions, outer products are not automatically
    created::

        >>> A = Operator('A')
        >>> A*k*b
        A*|k>*<b|

    A user can force the creation of an outer product in a complex expression
    by using parentheses to group the ket and bra::

        >>> A*(k*b)
        A*|k><b|

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Outer_product
    """

    def __new__(cls, ket, bra, **old_assumptions):
        if not isinstance(ket, KetBase):
            raise TypeError('KetBase subclass expected, got: %r' % ket)
        if not isinstance(bra, BraBase):
            raise TypeError('BraBase subclass expected, got: %r' % ket)
        if not ket.dual_class == bra.__class__:
            raise TypeError(
                'ket and bra are not dual classes: %r, %r' % \
                (ket.__class__, bra.__class__)
            )
        # TODO: make sure the hilbert spaces of the bra and ket are compatible
        obj = Expr.__new__(cls, *(ket, bra), **{'commutative': False})
        obj.hilbert_space = ket.hilbert_space
        return obj

    @property
    def ket(self):
        """Return the ket on the left side of the outer product."""
        return self.args[0]

    @property
    def bra(self):
        """Return the bra on the right side of the outer product."""
        return self.args[1]

    def _eval_dagger(self):
        return OuterProduct(Dagger(self.bra), Dagger(self.ket))

    def _sympyrepr(self, printer, *args):
        return '%s(%s,%s)' % (self.__class__.__name__, 
            printer._print(self.ket, *args), printer._print(self.bra, *args))

    def _sympystr(self, printer, *args):
        return str(self.ket)+str(self.bra)

    def _pretty(self, printer, *args):
        pform = self.ket._pretty(printer, *args)
        return prettyForm(*pform.right(self.bra._pretty(printer, *args)))

    def _latex(self, printer, *args):
        k = printer._print(self.ket, *args)
        b = printer._print(self.bra, *args)
        return k+b

#-----------------------------------------------------------------------------
# Subclasses of Expr
#-----------------------------------------------------------------------------


class KroneckerDelta(Function):
    """The discrete, or Kronecker, delta function.

    A function that takes in two integers i and j. It returns 0 if i and j are
    not equal or it returns 1 if i and j are equal.

    Parameters
    ==========
    i : Number, Symbol
        The first index of the delta function.
    j : Number, Symbol
        The second index of the delta function.

    Examples
    ========

    A simple example with integer indices::

        >>> from sympy.physics.quantum import KroneckerDelta
        >>> KroneckerDelta(1,2)
        0
        >>> KroneckerDelta(3,3)
        1

    Symbolic indices::

        >>> from sympy import symbols
        >>> i, j, k = symbols('i j k')
        >>> KroneckerDelta(i, j)
        d(i,j)
        >>> KroneckerDelta(i, i)
        1
        >>> KroneckerDelta(i, i+1)
        0
        >>> KroneckerDelta(i, i+1+k)
        d(i,1 + i + k)

    References
    ==========

    http://en.wikipedia.org/wiki/Kronecker_delta
    """

    nargs = 2
    is_commutative=True

    @classmethod
    def eval(cls, i, j):
        """
        Evaluates the discrete delta function.
        """
        if i > j:
            return cls(j,i)
        diff = i-j
        if diff == 0:
            return S.One
        elif diff.is_number:
            return S.Zero

    def _eval_subs(self, old, new):
        r = KroneckerDelta(self.args[0].subs(old, new), self.args[1].subs(old,\
        new))
        return r

    def _eval_dagger(self):
        return self

    def _latex_(self,printer):
        return "\\delta_{%s%s}"% (self.args[0].name,self.args[1].name)

    def _sympyrepr(self, printer, *args):
        return "%s(%s,%s)"% (self.__class__.__name__, self.args[0],\
        self.args[1])

    def _sympystr(self, printer, *args):
        return 'd(%s,%s)'% (self.args[0],self.args[1])

    def _pretty(self, printer, *args):
        pform = printer._print(self.args[0], *args)
        pform = prettyForm(*pform.right((prettyForm(','))))
        pform = prettyForm(*pform.right((printer._print(self.args[1], *args))))        
        a = stringPict(u'\u03b4')
        b = pform
        top = stringPict(*b.left(' '*a.width()))
        bot = stringPict(*a.right(' '*b.width()))
        return prettyForm(binding=prettyForm.POW, *bot.below(top))

    def _latex(self, printer, *args):
        i = printer._print(self.args[0], *args)
        j = printer._print(self.args[1], *args)
        return '\\delta_{%s %s}' % (i,j)


class Dagger(Expr):
    """General Hermitian conjugate operation.

    For matrices this operation is equivalent to transpose and complex
    conjugate [1].

    Parameters
    ==========
    arg : Expr
        The sympy expression that we want to take the dagger of.

    Examples
    ========

    Daggering various quantum objects:

        >>> from sympy.physics.quantum import Dagger, Ket, Bra, Operator
        >>> Dagger(Ket('psi'))
        <psi|
        >>> Dagger(Bra('phi'))
        |phi>
        >>> Dagger(Operator('A'))
        Dagger(A)

    Inner and outer products::

        >>> from sympy.physics.quantum import InnerProduct, OuterProduct
        >>> Dagger(InnerProduct(Bra('a'), Ket('b')))
        <b|a>
        >>> Dagger(OuterProduct(Ket('a'), Bra('b')))
        |b><a|

    Powers, sums and products::

        >>> A = Operator('A')
        >>> B = Operator('B')
        >>> Dagger(A*B)
        Dagger(B)*Dagger(A)
        >>> Dagger(A+B)
        Dagger(A) + Dagger(B)
        >>> Dagger(A**2)
        Dagger(A)**2

    Dagger also seamlessly handles complex numbers and matrices::

        >>> from sympy import Matrix, I
        >>> m = Matrix([[1,I],[2,I]])
        >>> m
        [1, I]
        [2, I]
        >>> Dagger(m)
        [ 1,  2]
        [-I, -I]

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Hermitian_transpose
    """

    def __new__(cls, arg, **old_assumptions):
        # Matrices are not sympify friendly
        if isinstance(arg, Matrix):
            return arg.H
        arg = sympify(arg)
        r = cls.eval(arg)
        if isinstance(r, Expr):
            return r
        #make unevaluated dagger commutative or non-commutative depending on arg
        if arg.is_commutative:
            obj = Expr.__new__(cls, arg, **{'commutative':True})
        else:
            obj = Expr.__new__(cls, arg, **{'commutative':False})
        if isinstance(obj, QExpr):
            obj.hilbert_space = arg.hilbert_space
        return obj

    @classmethod
    def eval(cls, arg):
        """
        Evaluates the Dagger instance.
        """
        try:
            d = arg._eval_dagger()
        except (NotImplementedError, AttributeError):
            if isinstance(arg, Expr):
                if isinstance(arg, Operator):
                    # Operator without _eval_dagger
                    return None
                if arg.is_Add:
                    return Add(*[Dagger(i) for i in arg.args])
                if arg.is_Mul:
                    return Mul(*[Dagger(i) for i in reversed(arg.args)])                    
                if arg.is_Pow:
                    return Pow(Dagger(arg.args[0]),arg.args[1])
                else:
                    return arg.conjugate()
            else:
                return None
        else:
            return d

    def _eval_subs(self, old, new):
        r = Dagger(self.args[0].subs(old, new))
        return r

    def _eval_dagger(self):
        return self.args[0]

    def _sympyrepr(self, printer, *args):
        return '%s(%s)' % (self.__class__.__name__, self.args[0])

    def _sympystr(self, printer, *args):
        return '%s(%s)' % (self.__class__.__name__, self.args[0])

    def _pretty(self, printer, *args):
        from sympy.printing.pretty.stringpict import prettyForm
        pform = printer._print(self.args[0], *args)
        pform = pform**prettyForm(u'\u2020')
        return pform

    def _latex(self, printer, *args):
        arg = printer._print(self.args[0])
        return '%s^{\\dag}' % arg


# InnerProduct is not an QExpr because it is really just a regular comutative
# number. We have gone back and forth about this, but we gain a lot by having
# it subclass Expr. The main challenges were getting Dagger to work
# (we use _eval_conjugate) and represent (we can use atoms and subs). Having
# it be an Expr, mean that there are no commutative QExpr subclasses, 
# which simplifies the design of everything.

class InnerProduct(Expr):
    """An unevaluated inner product between a Bra and a Ket.

    Parameters
    ==========
    bra : BraBase or subclass
        The bra on the left side of the inner product.
    ket : KetBase or subclass
        The ket on the right side of the inner product.

    Examples
    ========

    Create an InnerProduct and check its properties:

        >>> from sympy.physics.quantum import Bra, Ket, InnerProduct
        >>> b = Bra('b')
        >>> k = Ket('k')
        >>> ip = b*k
        >>> ip
        <b|k>
        >>> ip.bra
        <b|
        >>> ip.ket
        |k>

    In simple products of kets and bras inner products will be automatically
    identified and created::

        >>> b*k
        <b|k>

    But in more complex expressions, there is ambiguity in whether inner or
    outer products should be created::

        >>> k*b*k*b
        |k><b|*|k>*<b|

    A user can force the creation of a inner products in a complex expression
    by using parentheses to group the bra and ket::

        >>> k*(b*k)*b
        <b|k>*|k>*<b|

    Notice how the inner product <b|k> moved to the left of the expression
    because inner products are commutative complex numbers.

    References
    ==========

    http://en.wikipedia.org/wiki/Inner_product
    """

    def __new__(cls, bra, ket, **old_assumptions):
        if not isinstance(ket, KetBase):
            raise TypeError('KetBase subclass expected, got: %r' % ket)
        if not isinstance(bra, BraBase):
            raise TypeError('BraBase subclass expected, got: %r' % ket)
        obj = Expr.__new__(cls, *(bra, ket), **{'commutative':True})
        return obj

    @property
    def bra(self):
        return self.args[0]

    @property
    def ket(self):
        return self.args[1]

    def _eval_dagger(self):
        return InnerProduct(Dagger(self.ket), Dagger(self.bra))

    def _eval_conjugate(self):
        return self._eval_dagger()

    def _sympyrepr(self, printer, *args):
        return '%s(%s,%s)' % (self.__class__.__name__, 
            printer._print(self.bra, *args), printer._print(self.ket, *args))

    def _sympystr(self, printer, *args):
        sbra = str(self.bra)
        sket = str(self.ket)
        return '%s|%s' % (sbra[:-1], sket[1:])

    def _pretty(self, printer, *args):
        pform = prettyForm(u'\u276C')
        pform = prettyForm(*pform.right(self.bra._print_contents_pretty(printer, *args)))
        return prettyForm(*pform.right(self.ket._pretty(printer, *args)))

    def _latex(self, printer, *args):
        bra_label = self.bra._print_contents(printer, *args)
        ket = printer._print(self.ket, *args)
        return r'\left\langle %s \right. %s' % (bra_label, ket)

    def doit(self, **hints):
        try:
            r = self.ket._eval_innerproduct(self.bra, **hints)
        except NotImplementedError:
            try:
                r = self.bra._eval_innerproduct(self.ket, **hints)
            except NotImplementedError:
                r = None
        if r is not None:
            return r
        return self


class Commutator(Expr):
    """The standard commutator, in an unevaluated state.

    The commutator is defined [1] as: [A, B] = A*B - B*A, but in this class
    the commutator is initially unevaluated. To muliple the commutator out,
    use the ``doit`` method.

    The arguments of the commutator are put into canonical order using
    ``__cmp__``, so that [B,A] becomes -[A,B].

    Parameters
    ==========
    A : Expr
        The first argument of the commutator [A,B].
    B : Expr
        The second argument of the commutator [A,B].

    Examples
    ========

        >>> from sympy import symbols
        >>> from sympy.physics.quantum import Commutator, Operator, Dagger
        >>> x, y = symbols('xy')
        >>> A = Operator('A')
        >>> B = Operator('B')
        >>> C = Operator('C')

    Create some commutators and use ``doit`` to multiply them out.

        >>> comm = Commutator(A,B); comm
        [A,B]
        >>> comm.doit()
        A*B - B*A

    The commutator orders it arguments in canonical order::

        >>> comm = Commutator(B,A); comm
        -[A,B]

    Scalar constants are factored out::

        >>> Commutator(3*x*A,x*y*B)
        3*y*x**2*[A,B]

    Using ``expand(commutator=True)``, the standard commutator expansion rules
    can be applied::

        >>> Commutator(A+B,C).expand(commutator=True)
        [A,C] + [B,C]
        >>> Commutator(A,B+C).expand(commutator=True)
        [A,B] + [A,C]
        >>> Commutator(A*B,C).expand(commutator=True)
        [A,C]*B + A*[B,C]
        >>> Commutator(A,B*C).expand(commutator=True)
        [A,B]*C + B*[A,C]

    Commutator works with Dagger::

        >>> Dagger(Commutator(A,B))
        -[Dagger(A),Dagger(B)]

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Commutator
    """

    def __new__(cls, A, B, **old_assumptions):
        r = cls.eval(A, B)
        if r is not None:
            return r
        obj = Expr.__new__(cls, *(A, B), **{'commutative': False})
        return obj

    @classmethod
    def eval(cls, a, b):
        """The Commutator [A,B] is on canonical form if A < B.
        """
        if not (a and b): return S.Zero
        if a == b: return S.Zero
        if a.is_commutative or b.is_commutative:
            return S.Zero

        # [xA,yB]  ->  xy*[A,B]
        # from sympy.physics.qmul import QMul
        c_part = c_part2 = []
        nc_part = nc_part2 = []
        if isinstance(a, Mul):
            c_part, nc_part = split_commutative_parts(a)
        if isinstance(b, Mul):
            c_part2, nc_part2 = split_commutative_parts(b)
            c_part.extend(c_part2)
        if c_part:
            a = nc_part or [a]
            b = nc_part2 or [b]
            return Mul(*c_part)*cls(Mul(*a),Mul(*b))

        # Canonical ordering of arguments
        if a.compare(b) == 1:
            return S.NegativeOne*cls(b,a)

    def _eval_expand_commutator(self, **hints):
        A = self.args[0].expand(**hints)
        B = self.args[1].expand(**hints)

        result = None

        if isinstance(A, Add):
            # [A+B,C]  ->  [A,C] + [B,C]
            result = Add(
                *[Commutator(term,B).expand(**hints)\
                  for term in A.args]
            )
        elif isinstance(B, Add):
            # [A,B+C]  ->  [A,B] + [A,C]
            result = Add(
                *[Commutator(A,term).expand(**hints)\
                  for term in B.args]
            )
        elif isinstance(A, Mul):
            # [A*B,C] -> A*[B,C] + [A,C]*B
            a = A.args[0]
            b = Mul(*A.args[1:])
            c = B
            comm1 = Commutator(b,c).expand(**hints)
            comm2 = Commutator(a,c).expand(**hints)
            first = Mul(a, comm1)
            second = Mul(comm2, b)
            result = Add(first, second)
        elif isinstance(B, Mul):
            # [A,B*C] -> [A,B]*C + B*[A,C]
            a = A
            b = B.args[0]
            c = Mul(*B.args[1:])
            comm1 = Commutator(a,b).expand(**hints)
            comm2 = Commutator(a,c).expand(**hints)
            first = Mul(comm1, c)
            second = Mul(b, comm2)
            result = Add(first, second)

        if result is None:
            # No changes, so return self
            return self
        else:
            return result

    def doit(self, **hints):
        A = self.args[0]
        B = self.args[1]
        if isinstance(A, Operator) and isinstance(B, Operator):
            try:
                comm = A._eval_commutator(B, **hints)
            except NotImplementedError:
                try:
                    comm = -1*B._eval_commutator(A, **hints)
                except NotImplementedError:
                    comm = None
            if comm is not None:
                return comm.doit(**hints)
        return (A*B - B*A).doit(**hints)

    def represent(self, basis, **options):
        # TODO: should Commutator know how to represent?
        rep_method = '_represent_%s' % basis.__class__.__name__
        if hasattr(self, rep_method):
            f = getattr(self, rep_method)
            rep = f(basis, **options)
            if rep is not None:
                return rep
        raise NotImplementedError("Can't represent %r in basis: %r" % (
            self, basis
        ))

    def _eval_dagger(self):
        return Commutator(Dagger(self.args[1]), Dagger(self.args[0]))

    def _sympyrepr(self, printer, *args):
        return "%s(%s,%s)" % (self.__class__.__name__, self.args[0],\
        self.args[1])

    def _sympystr(self, printer, *args):
        return "[%s,%s]" % (self.args[0], self.args[1])

    def _pretty(self, printer, *args):
        pform = printer._print(self.args[0], *args)
        pform = prettyForm(*pform.right((prettyForm(u','))))
        pform = prettyForm(*pform.right((printer._print(self.args[1], *args))))
        pform = prettyForm(*pform.parens(left='[', right=']'))
        return pform

    def _latex(self, printer, *args):
        return "\\left[%s,%s\\right]" % tuple([
            printer._print(arg, *args) for arg in self.args])


class AntiCommutator(Expr):
    """The standard anicommutator, in an unevaluated state.

    The commutator is defined [1] as: {A, B} = A*B + B*A, but in this class
    the anticommutator is initially unevaluated. To muliple the anticommutator
    out, use the ``doit`` method.

    The arguments of the anticommutator are put into canonical order using
    ``__cmp__``, so that {B,A} becomes {A,B}.

    Parameters
    ==========
    A : Expr
        The first argument of the anticommutator {A,B}.
    B : Expr
        The second argument of the anticommutator {A,B}.

    Examples
    ========

        >>> from sympy import symbols
        >>> from sympy.physics.quantum import AntiCommutator, Operator, Dagger
        >>> x, y = symbols('xy')
        >>> A = Operator('A')
        >>> B = Operator('B')

    Create an anticommutator and use ``doit`` to multiply them out.

        >>> ac = AntiCommutator(A,B); ac
        {A,B}
        >>> ac.doit()
        A*B + B*A

    The commutator orders it arguments in canonical order::

        >>> ac = AntiCommutator(B,A); ac
        {A,B}

    Scalar constants are factored out::

        >>> AntiCommutator(3*x*A,x*y*B)
        3*y*x**2*{A,B}

    Dagger is alto handled::

        >>> Dagger(AntiCommutator(A,B))
        {Dagger(A),Dagger(B)}

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Commutator
    """

    def __new__(cls, A, B, **old_assumptions):
        r = cls.eval(A, B)
        if r is not None:
            return r
        obj = Expr.__new__(cls, *(A, B), **{'commutative': False})
        return obj

    @classmethod
    def eval(cls, a, b):
        """The Commutator [A,B] is on canonical form if A < B.
        """
        if not (a and b): return S.Zero
        if a == b: return Integer(2)*a**2
        if a.is_commutative or b.is_commutative:
            return Integer(2)*a*b

        # [xA,yB]  ->  xy*[A,B]
        # from sympy.physics.qmul import QMul
        c_part = []
        nc_part = []
        nc_part2 = []
        if isinstance(a, Mul):
            c_part, nc_part = split_commutative_parts(a)
        if isinstance(b, Mul):
            c_part2, nc_part2 = split_commutative_parts(b)
            c_part.extend(c_part2)
        if c_part:
            a = nc_part or [a]
            b = nc_part2 or [b]
            return Mul(Mul(*c_part), cls(Mul(*a), Mul(*b)))

        # Canonical ordering of arguments
        if a.compare(b) == 1:
            return cls(b,a)

    def _eval_expand_anticommutator(self, **hints):
        # No changes, so return self
        return self

    def doit(self, **hints):
        A = self.args[0]
        B = self.args[1]
        if isinstance(A, Operator) and isinstance(B, Operator):
            try:
                comm = A._eval_anticommutator(B, **hints)
            except NotImplementedError:
                try:
                    comm = -1*B._eval_anticommutator(A, **hints)
                except NotImplementedError:
                    comm = None
            if comm is not None:
                return comm.doit(**hints)
        return (A*B + B*A).doit(**hints)

    def _eval_dagger(self):
        return AntiCommutator(Dagger(self.args[0]), Dagger(self.args[1]))

    def _sympyrepr(self, printer, *args):
        return "%s(%s,%s)" % (self.__class__.__name__, self.args[0],\
        self.args[1])

    def _sympystr(self, printer, *args):
        return "{%s,%s}" % (self.args[0], self.args[1])

    def _pretty(self, printer, *args):
        pform = printer._print(self.args[0], *args)
        pform = prettyForm(*pform.right((prettyForm(u','))))
        pform = prettyForm(*pform.right((printer._print(self.args[1], *args))))
        pform = prettyForm(*pform.parens(left='{', right='}'))
        return pform

    def _latex(self, printer, *args):
        return "\\left{}%s,%s\\right}" % tuple([
            printer._print(arg, *args) for arg in self.args])


#-----------------------------------------------------------------------------
# Tensor product
#-----------------------------------------------------------------------------

# TODO: move this to the main sympy.matrices module.
def matrix_tensor_product(*matrices):
    """Compute the tensor product of a sequence of sympy Matrices.

    This is the standard Kronecker product of matrices [1].

    Parameters
    ==========
    matrices : tuple of Matrix instances
        The matrices to take the tensor product of.

    Returns
    =======
    matrix : Matrix
        The tensor product matrix.

    Examples
    ========

        >>> from sympy import I, Matrix, symbols
        >>> from sympy.physics.quantum import matrix_tensor_product

        >>> m1 = Matrix([[1,2],[3,4]])
        >>> m2 = Matrix([[1,0],[0,1]])
        >>> matrix_tensor_product(m1, m2)
        [1, 0, 2, 0]
        [0, 1, 0, 2]
        [3, 0, 4, 0]
        [0, 3, 0, 4]
        >>> matrix_tensor_product(m2, m1)
        [1, 2, 0, 0]
        [3, 4, 0, 0]
        [0, 0, 1, 2]
        [0, 0, 3, 4]

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Kronecker_product
    """
    # Make sure we have a sequence of Matrices
    testmat = [isinstance(m, Matrix) for m in matrices]
    if not all(testmat):
        raise TypeError('Sequence of Matrices expected, got: %r' % matrices)

    # Pull out the first element in the product.
    matrix_expansion  = matrices[-1]
    # Do the tensor product working from right to left.
    for mat in reversed(matrices[:-1]):
        rows = mat.rows
        cols = mat.cols
        # Go through each row appending tensor product to.
        # running matrix_expansion.
        for i in range(rows):
            start = matrix_expansion*mat[i*cols]
            # Go through each column joining each item
            for j in range(cols-1):
                start = start.row_join(
                    matrix_expansion*mat[i*cols+j+1]
                )
            # If this is the first element, make it the start of the
            # new row.
            if i == 0:
                next = start
            else:
                next = next.col_join(start)
        matrix_expansion = next
    return matrix_expansion


class TensorProduct(Expr):
    """The tensor product of two or more arguments.

    For matrices, this uses ``matrix_tensor_product`` to compute the 
    Kronecker or tensor product matrix. For other objects a symbolic
    ``TensorProduct`` instance is returned. The tensor product is a
    non-commutative multiplication that is used primarily with operators
    and states in quantum mechanics.

    Current, the tensor product distinguishes between commutative and non-
    commutative arguments.  Commutative arguments are assumed to be scalars
    and are pulled out in front of the ``TensorProduct``. Non-commutative
    arguments remain in the resulting ``TensorProduct``.

    Parameters
    ==========
    args : tuple
        The objects to take the tensor product of.

    Examples
    ========

    Start with a simple tensor product of sympy matrices::

        >>> from sympy import I, Matrix, symbols
        >>> from sympy.physics.quantum import TensorProduct

        >>> m1 = Matrix([[1,2],[3,4]])
        >>> m2 = Matrix([[1,0],[0,1]])
        >>> TensorProduct(m1, m2)
        [1, 0, 2, 0]
        [0, 1, 0, 2]
        [3, 0, 4, 0]
        [0, 3, 0, 4]
        >>> TensorProduct(m2, m1)
        [1, 2, 0, 0]
        [3, 4, 0, 0]
        [0, 0, 1, 2]
        [0, 0, 3, 4]

    We can also construct tensor products of non-commutative symbols::

        >>> from sympy import Symbol
        >>> A = Symbol('A',commutative=False)
        >>> B = Symbol('B',commutative=False)
        >>> tp = TensorProduct(A, B)
        >>> tp
        (AxB)

    We can take the dagger of a tensor product::

        >>> from sympy.physics.quantum import Dagger
        >>> Dagger(tp)
        (conjugate(A)xconjugate(B))

    Expand can be used to distribute a tensor product across addition::

        >>> C = Symbol('C',commutative=False)
        >>> tp = TensorProduct(A+B,C)
        >>> tp
        ((A + B)xC)
        >>> tp.expand(tensorproduct=True)
        (AxC) + (BxC)
    """

    def __new__(cls, *args, **assumptions):
        if isinstance(args[0], Matrix):
            return matrix_tensor_product(*args)
        c_part, new_args = cls.flatten(args)
        c_part = Mul(*c_part)
        if len(new_args) == 0:
            return c_part
        elif len(new_args) == 1:
            return c_part*new_args[0]
        else:
            tp = Expr.__new__(cls, *new_args, **{'commutative': False})
            return c_part*tp

    @classmethod
    def flatten(cls, args):
        # TODO: disallow nested TensorProducts.
        c_part = []
        nc_parts = []
        for arg in args:
            if isinstance(arg, Mul):
                cp, ncp = split_commutative_parts(arg)
                ncp = Mul(*ncp)
            else:
                if arg.is_commutative:
                    cp = [arg]; ncp = 1
                else:
                    cp = []; ncp = arg
            c_part.extend(cp)
            nc_parts.append(ncp)
        return c_part, nc_parts

    def _eval_dagger(self):
        return TensorProduct(*[Dagger(i) for i in self.args])

    def _sympystr(self, printer, *args):
        from sympy.printing.str import sstr
        length = len(self.args)
        s = ''
        for i in range(length):
            if isinstance(self.args[i], (Add, Pow, Mul)):
                s = s + '('
            s = s + sstr(self.args[i])
            if isinstance(self.args[i], (Add, Pow, Mul)):
                s = s + ')'
            if i != length-1:
                s = s + 'x'
        return s

    def _pretty(self, printer, *args):
        length = len(self.args)
        pform = printer._print('', *args)
        for i in range(length):
            next_pform = printer._print(self.args[i], *args)
            if isinstance(self.args[i], (Add, Mul)):
                next_pform = prettyForm(
                    *next_pform.parens(left='(', right=')')
                )
            pform = prettyForm(*pform.right(next_pform))
            if i != length-1:
                pform = prettyForm(*pform.right(u'\u2a02' + u' '))
        return pform

    def _latex(self, printer, *args):
        length = len(self.args)
        s = ''
        for i in range(length):
            if isinstance(self.args[i], (Add, Mul)):
                s = s + '\\left('
            s = s + printer._print(self.args[i], *args)
            if isinstance(self.args[i], (Add, Mul)):
                s = s + '\\right)'
            if i != length-1:
                s = s + '\\otimes '
        return s

    def doit(self, **hints):
        return TensorProduct(*[item.doit(**hints) for item in self.args])

    def _eval_expand_tensorproduct(self, **hints):
        """Distribute TensorProducts across addition."""
        args = self.args
        add_args = []
        stop = False
        for i in range(len(args)):
            if isinstance(args[i], Add):
                for aa in args[i].args:
                    add_args.append(TensorProduct(*args[:i]+(aa,)+args[i+1:]))
                stop = True
            if stop: break
        if add_args:
            return Add(*add_args).expand(**hints)
        else:
            return self

    def expand(self, **hints):
        tp = TensorProduct(*[item.expand(**hints) for item in self.args])
        return Expr.expand(tp, **hints)


def tensor_product_simp_Mul(e):
    """Simplify a Mul with TensorProducts.

    Current the main use of this is to simplify a ``Mul`` of
    ``TensorProduct``s to a ``TensorProduct`` of ``Muls``. It currently only
    works for relatively simple cases where the initial ``Mul`` only has
    scalars and raw ``TensorProduct``s, not ``Add``, ``Pow``, ``Commutator``s
    of ``TensorProduct``s.

    Parameters
    ==========
    e : Expr
        A ``Mul`` containing ``TensorProduct``s to be simplified.

    Returns
    =======
    e : Expr
        A ``TensorProduct`` of ``Mul``s.

    Examples
    ========

    This is an example of the type of simplification that this function
    performs::

        >>> from sympy.physics.quantum import tensor_product_simp_Mul, TensorProduct
        >>> from sympy import Symbol
        >>> A = Symbol('A',commutative=False)
        >>> B = Symbol('B',commutative=False)
        >>> C = Symbol('C',commutative=False)
        >>> D = Symbol('D',commutative=False)
        >>> e = TensorProduct(A,B)*TensorProduct(C,D)
        >>> e
        (AxB)*(CxD)
        >>> tensor_product_simp_Mul(e)
        ((A*C)x(B*D))

    """
    # TODO: This won't work with Muls that have other composites of 
    # TensorProducts, like an Add, Pow, Commutator, etc.
    # TODO: This only works for the equivalent of single Qbit gates.
    if not isinstance(e, Mul):
        return e
    c_part, nc_part = split_commutative_parts(e)
    n_nc = len(nc_part)
    if n_nc == 0 or n_nc == 1:
        return e
    elif e.has(TensorProduct):
        current = nc_part[0]
        if not isinstance(current, TensorProduct):
            raise TypeError('TensorProduct expected, got: %r' % current)
        n_terms = len(current.args)
        new_args = list(current.args)
        for next in nc_part[1:]:
            # TODO: check the hilbert spaces of next and current here.
            if isinstance(next, TensorProduct):
                if n_terms != len(next.args):
                    raise QuantumError(
                        'TensorProducts of different lengths: %r and %r' % \
                        (current, next)
                    )
                for i in range(len(new_args)):
                    new_args[i] = new_args[i]*next.args[i]
            else:
                # this won't quite work as we don't want next in the TensorProduct
                for i in range(len(new_args)):
                    new_args[i] = new_args[i]*next
            current = next
        return Mul(*c_part)*TensorProduct(*new_args)
    else:
        return e


def tensor_product_simp(e, **hints):
    """Try to simplify and combine TensorProducts.

    In general this will try to pull expressions inside of ``TensorProducts``.
    It is best to see what it does by showing examples.

    Examples
    ========

        >>> from sympy.physics.quantum import tensor_product_simp
        >>> from sympy.physics.quantum import TensorProduct
        >>> from sympy import Symbol
        >>> A = Symbol('A',commutative=False)
        >>> B = Symbol('B',commutative=False)
        >>> C = Symbol('C',commutative=False)
        >>> D = Symbol('D',commutative=False)

    First see what happens to products of tensor products::

        >>> e = TensorProduct(A,B)*TensorProduct(C,D)
        >>> e
        (AxB)*(CxD)
        >>> tensor_product_simp(e)
        ((A*C)x(B*D))

    This is the core logic of this function, and it works inside, powers,
    sums, commutators and anticommutators as well::

        >>> tensor_product_simp(e**2)
        ((A*C)x(B*D))**2    
    """
    if isinstance(e, Add):
        return Add(*[tensor_product_simp(arg) for arg in e.args])
    elif isinstance(e, Pow):
        return tensor_product_simp(e.base)**e.exp
    elif isinstance(e, Mul):
        return tensor_product_simp_Mul(e)
    elif isinstance(e, Commutator):
        return Commutator(*[tensor_product_simp(arg) for arg in e.args])
    elif isinstance(e, AntiCommutator):
        return AntiCommutator(*[tensor_product_simp(arg) for arg in e.args])
    else:
        return e

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------


def represent(expr, basis, **options):
    """Represent the quantum expression in the given basis.

    In quantum mechanics abstract states and operators can be represented in
    various basis sets. Under this operator, states become vectors or
    functions and operators become matrices or differential operators. This
    function is the top-level interface for this action.

    This function walks the sympy expression tree looking for ``QExpr``
    instances that have a ``_represent`` method. This method is then called
    and the object is replaced by the representation returned by this method.
    By default, the ``_represent`` method will dispatch to other methods
    that handle the representation logic for a particular basis set. The
    naming convention for these methods is the following::

        def _represent_FooBasis(self, e, basis, **options)

    This function will have the logic for representing instances of its class
    in the basis set having a class named ``FooBasis``.

    Parameters
    ==========
    expr  : Expr
        The expression to represent.
    basis : Operator, basis set
        An object that contains the information about the basis set. If an
        operator is used, the basis is assumed to be the orthonormal
        eigenvectors of that operator. In general though, the basis argument
        can be any object that contains the basis set information.
    options : dict
        Key/value pairs of options that are passed to the underlying method
        that does finds the representation. These options can be used to
        control how the representation is done. For example, this is where
        the size of the basis set would be set.

    Returns
    =======
    e : Expr
        The sympy expression of the represented quantum expression.

    Examples
    ========

    Here we subclass ``Operator`` and ``Ket`` to create the z-spin operator
    and its spin 1/2 up eigenstate. By definining the ``_represent_SzOp``
    method, the ket can be represented in the z-spin basis.

        >>> from sympy.physics.quantum import Operator, represent, Ket
        >>> from sympy import Matrix

        >>> class SzUpKet(Ket):
        ...     def _represent_SzOp(self, basis, **options):
        ...         return Matrix([1,0])
        ...     
        >>> class SzOp(Operator):
        ...     pass
        ... 
        >>> sz = SzOp('Sz')
        >>> up = SzUpKet('up')
        >>> represent(up, sz)
        [1]
        [0]
    """
    if isinstance(expr, QExpr):
        return expr._represent(basis, **options)
    elif isinstance(expr, Add):
        result = S.Zero
        for args in expr.args:
            if not result:
                result = represent(args, basis, **options)
            else:
                result += represent(args, basis, **options)
        return result
    elif isinstance(expr, Pow):
        return represent(expr.base, basis, **options)**expr.exp
    elif isinstance(expr, TensorProduct):
        new_args = [represent(arg, basis, **options) for arg in expr.args]
        return TensorProduct(*new_args)
    elif not isinstance(expr, Mul):
        return expr

    if not isinstance(expr, Mul):
        raise TypeError('Mul expected, got: %r' % expr)

    result = S.One
    for arg in reversed(expr.args):
        result = represent(arg, basis, **options)*result
    if isinstance(result, Matrix):
        if result.shape == (1,1):
            result = result[0]
    return result


def apply_operators(e, **options):
    """Apply operators to states in a quantum expression.

    By default, operators acting on states (O|psi>) are left in symbolic,
    unevaluated form. This function uses various special methods to attempt
    and apply operators to states. When this happens there are two possible
    outcomes: i) it is not known how the operator acts on the state, which
    will result in the expression being unchanged and ii) it is known how the
    operator acts on the sate, which will result in the action being carried
    out.

    Parameters
    ==========
    e : Expr
        The expression containing operators and states. This expression tree
        will be walked to find operators acting on states symbolically.
    options : dict
        A dict of key/value pairs that determine how the operator actions
        are carried out.

    Returns
    =======
    e : Expr
        The original expression, but with the operators applied to states.
    """

    # TODO: Fix early 0 in apply_operators.

    # This may be a bit aggressive but ensures that everything gets expanded
    # to its simplest form before trying to apply operators. This includes
    # things like (A+B+C)*|a> and A*(|a>+|b>) and all Commutators and
    # TensorProducts. The only problem with this is that if we can't apply
    # all the Operators, we have just expanded everything.
    # TODO: don't expand the scalars in front of each Mul.
    e = e.expand(commutator=True, tensorproduct=True).doit()

    # If we just have a raw ket, return it.
    # TODO: make this acts_like_KetBase
    if isinstance(e, KetBase):
        return e

    # We have an Add(a, b, c, ...) and compute 
    # Add(apply_operators(a), apply_operators(b), ...)
    elif isinstance(e, Add):
        result = 0
        for arg in e.args:
            result += apply_operators(arg, **options)
        return result

    # We have a Mul where there might be actual operators to apply to kets.
    elif isinstance(e, Mul):
        return apply_operators_Mul(e, **options)

    # In all other cases (State, Operator, Pow, Commutator, InnerProduct,
    # OuterProduct) we won't ever have operators to apply to kets.
    else:
        return e


def apply_operators_Mul(e, **options):

    # The general form of e at this point is:
    #
    # ZeroOrMore(scalars)*ZeroOrMore(Bra)*ZeroOrMore(Operators)*Ket*
    # ZeroOrMore(Bra*ZeroMore(Operators)*Ket)
    # 
    # We want to pick out pieces that look like:
    #
    # Bra*ZeroOrMore(Operators)*Ket

    if not isinstance(e, Mul):
        return e

    # Split the Mul into a Mul of scalars and a list of non-commutative parts.
    c_part, nc_part = split_commutative_parts(e)
    c_part = Mul(*c_part)
    n_nc = len(nc_part)
    if n_nc == 0 or n_nc == 1:
        return c_part

    terms = []
    last_ket = 0
    for i in range(n_nc):
        # Make this acts_like_KetBase
        if isinstance(nc_part[i], KetBase):
            terms.append(nc_part[last_ket:i+1])
            last_ket = i+1
    last_part = nc_part[last_ket:]
    if last_part:
        terms.append(last_part)

    result = 1
    for term in terms:
        # print 'term', term

        bra = 1
        # TODO: make this acts_like_BraBase
        if isinstance(term[0], BraBase):
            bra = term.pop(0)

        ket = 1
        # TODO: make this acts_like(term[-1], KetBase)
        if isinstance(term[-1], KetBase):
            ket = term.pop(-1)

        # print 'bra', bra
        # print 'ket', ket

        if bra == 1 and ket == 1:
            tresult = Mul(*term)

        elif len(term) == 0:
            # This cover the cases where either bra or ket is 1 and will also
            # automatically create an InnerProduct if neither are 1.
            tresult = bra*ket

        # The <bra|*A*B case
        elif ket == 1:
            try:
                new_term = [Dagger(i) for i in reversed(term)]
            except NotImplementedError:
                tresult = bra*Mul(*term)
            else:
                # TODO: what is the best way of handling copy.copy in all this
                tresult = Dagger(
                    apply_list_of_ops(1, copy.copy(new_term), Dagger(bra))
                )

        # The full <bra|*A*B*|ket> case
        # TODO: Get scalar<bra|*|ket> into an InnerProduct.
        else:
            try:
                tresult = bra*apply_list_of_ops(1, copy.copy(term), ket)
            except NotImplementedError:
                tresult = bra*Mul(*term)*ket

        # print 'tresult', tresult
        # print
        result *= tresult

    return c_part*result


def apply_list_of_ops(scalar, ops, ket):
    # scalar is a number
    # ops is a list of Operators or Powers of Operators.
    # ket can be a*|alpha> + b*|beta> + ... but cannot include any Operators.

    if scalar == 0 or ket == 0:
        return S.Zero

    # print scalar, ops, ket

    if not ops:
        return scalar*ket

    if isinstance(ket, Add):
        result = 0
        for arg in ket.args:
            if isinstance(arg, KetBase):
                result += apply_list_of_ops(scalar, ops, arg)
            elif isinstance(arg, Mul):
                k = arg.args[-1]
                s = Mul(scalar, *arg.args[:-1])
                result += apply_list_of_ops(s, copy.copy(ops), k)
            else:
                raise TypeError('Ket or Mul expected, got: %r' % arg)
        return result

    if isinstance(ket, Mul):
        scalar = Mul(scalar, *ket.args[:-1])
        ket = ket.args[-1]

    next_op = ops.pop(-1)
    if isinstance(next_op, Pow):
        ops.append(next_op.base**(next_op.exp-1))
        next_op = next_op.base

    # If unsuccessful, this will raise NotImplementedError which we let
    # propagate.
    new_ket = apply_single_op(next_op, ket)

    if new_ket == 0:
        return 0

    # TODO: Try is without this expand, I don't think we will need it.
    # We definitely need it!!!
    new_ket = new_ket.expand()

    return apply_list_of_ops(scalar, ops, new_ket)


def apply_single_op(op, ket):
    # Apply a single Operator to a Ket, at this point, we must have only an
    # Operator and a Ket and nothing else!
    if not isinstance(ket, KetBase):
        raise TypeError('Ket expected, got: %r' % ket)
    if not isinstance(op, Operator):
        raise TypeError('Operator expected, got: %r' % op)
    try:
        result = op._apply_operator(ket)
    except NotImplementedError:
        try:
            result = ket._apply_operator(op)
        except NotImplementedError:
            raise
    return result


#-----------------------------------------------------------------------------
# Constants
#-----------------------------------------------------------------------------

class HBar(NumberSymbol):
    """Reduced Plank's constant in numerical and symbolic form [1].

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Planck_constant
    """

    is_real = True
    is_positive = True
    is_negative = False
    is_irrational = True

    __slots__ = []

    def _as_mpf_val(self, prec):
        return mlib.from_float(1.05457162e-34, prec)

    def _sympyrepr(self, printer, *args):
        return 'HBar()'

    def _sympystr(self, printer, *args):
        return 'hbar'

    def _pretty(self, printer, *args):
        return prettyForm(u'\u210f')

# Create an instance for everyone to use.
hbar = HBar()

