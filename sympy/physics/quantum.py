"""Dirac notation for Bras, Kets, Operators, and so on.

TODO:

* tpsimp.
* Fix early 0 in apply_operators.
* Debug and test apply_operators.
* Get cse working with classes in this file.
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
    'tpsimp',
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
        raise NotImplementedError(
            'dual_class must be implemented in a subclass'
        )

    def _eval_dagger(self):
        """Compute the Dagger of this state."""
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
    """General abstract quantum state."""
    pass



class Ket(State, KetBase):
    """A Ket state for quantum mechanics.

    Inherits from State and KetBase. In a state represented by a ray in Hilbert
    space, a Ket is a vector that points along that ray [1].

    A Ket takes in a label as its argument in order to be differentiated from
    other Kets.

    Examples
    ========

    Creating and using a Ket:

        >>> from sympy.physics.quantum import Ket
        >>> psi = Ket('psi')
        >>> psi
        |psi>
        >>> psi.label
        psi
        >>> psi.dual
        <psi|


    References
    ==========

    [1] http://en.wikipedia.org/wiki/Bra-ket_notation
    """

    @property
    def dual_class(self):
        return Bra


class Bra(State, BraBase):
    """A Bra state for quantum mechanics.

    Inherits from State and BraBase. A Bra is the dual of a Ket [1].

    A Bra takes in a label as its argument in order to be differentiated from
    other Bras.

    Examples
    ========

    Creating and using a Bra:

        >>> from sympy.physics.quantum import Bra
        >>> b = Bra('bus')
        >>> b
        <bus|
        >>> b.label
        bus
        >>> b.dual
        |bus>


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
    """General abstract time dependent quantum state.

    Used for sub-classing time dependent states in quantum mechanics.

    Takes in label and time arguments.
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
    """A time dependent Ket state for quantum mechanics.

    Inherits from TimeDepState and KetBase. Its dual is a time dependent
    Bra state.

    Takes in label and time arguments.
    """

    @property
    def dual_class(self):
        return TimeDepBra


class TimeDepBra(TimeDepState, BraBase):
    """A time dependent Bra state for quantum mechanics.

    Inherits from TimeDepState and BraBase. Its dual is a time dependent
    Ket state.

    Takes in label and time arguments.
    """

    @property
    def dual_class(self):
        return TimeDepKet

#-----------------------------------------------------------------------------
# Operators and outer products
#-----------------------------------------------------------------------------

class Operator(QExpr):
    """Base class for non-commuting quantum operators.

    An operator is a map from one vector space to another [1]. In quantum
    mechanics, operators correspond to observables [2].

    An operator takes in a name argument to differentiate it from other
    operators.

    Examples
    ========

    Creating and using an Operator.

        >>> from sympy import sympify
        >>> from sympy.physics.quantum import Operator
        >>> a = Operator('A')
        >>> a
        A
        >>> a.name
        A

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
    """A Hermitian operator"""

    def _eval_dagger(self):
        return self


class UnitaryOperator(Operator):
    """A unitary operator."""

    def _eval_dagger(self):
        return self._eval_inverse()


class OuterProduct(Operator):
    """An unevaluated outer product between a Ket and Bra.

    Because a Ket is essentially a column vector and a Bra is essentially a row
    vector, the outer product evaluates (acts_like) to an operator
    (i.e. a matrix) [1].

    An InnerProduct takes first a Ket and then a Bra as its arguments.

    Examples
    ========

    Create an OuterProduct and check its properties:

        >>> from sympy.physics.quantum import Ket, Bra, OuterProduct
        >>> op = OuterProduct(Ket('a'), Bra('b'))
        >>> op
        |a><b|
        >>> op.ket
        |a>
        >>> op.bra
        <b|

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
        # TODO: fix this, to allow different symbolic dimensions
        # if not ket.hilbert_space == bra.hilbert_space:
        #     raise HilbertSpaceError(
        #         'Incompatible hilbert spaces: %r and %r' % (ket, bra)
        #     )
        obj = Expr.__new__(cls, *(ket, bra), **{'commutative': False})
        obj.hilbert_space = ket.hilbert_space
        return obj

    @property
    def ket(self):
        return self.args[0]

    @property
    def bra(self):
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


#-----------------------------------------------------------------------------
# Subclasses of Expr
#-----------------------------------------------------------------------------


class KroneckerDelta(Function):
    """The discrete delta function.

    A function that takes in two integers i and j. It returns 0 if i and j are
    not equal or it returns 1 if i and j are equal.

    Examples
    ========

        >>> from sympy.physics.quantum import KroneckerDelta
        >>> KroneckerDelta(1,2)
        0
        >>> KroneckerDelta(3,3)
        1

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


class Dagger(Expr):
    """General Hermitian conjugate operation.

    When an object is daggered it is transposed and then the complex conjugate
    is taken [1].

    A Dagger can take in any quantum object as its argument.

    Examples
    ========

    Daggering various quantum objects:

        >>> from sympy.physics.quantum import Dagger, Ket, Bra, InnerProduct,\
        OuterProduct, Operator
        >>> Dagger(Ket('psi'))
        <psi|
        >>> Dagger(Bra('bus'))
        |bus>
        >>> Dagger(InnerProduct(Bra('a'), Ket('b')))
        <b|a>
        >>> Dagger(OuterProduct(Ket('a'), Bra('b')))
        |b><a|
        >>> Dagger(Operator('A'))
        Dagger(A)

    Notice how daggering the operator 'A' results in isympy returning an
    unevaluated Dagger object (because A had no _eval_dagger method).

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


# InnerProduct is not an QExpr because it is really just a regular comutative
# number. We have gone back and forth about this, but we gain a lot by having
# this. The main challenges were getting Dagger to work (we use _eval_conjugate)
# and represent (we can use atoms and subs). Having it be an Expr, mean that
# there are no commutative QExpr subclasses, which simplifies the logic in QMul
# a lot.

class InnerProduct(Expr):
    """An unevaluated inner product between a Bra and a Ket.

    Because a Bra is essentially a row vector and a Ket is essentially a column
    vector, the inner product evaluates (acts_like) to a complex number.

    An InnerProduct takes first a Bra and then a Ket as its arguments.

    Examples
    ========

    Create an InnerProduct and check its properties:

        >>> from sympy.physics.quantum import Bra, Ket, InnerProduct
        >>> ip = InnerProduct(Bra('a'), Ket('b'))
        >>> ip
        <a|b>
        >>> ip.bra
        <a|
        >>> ip.ket
        |b>

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
        pform = prettyForm(*pform.right(self.bra._print_label_pretty(printer, *args)))
        return prettyForm(*pform.right(self.ket._pretty(printer, *args)))

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
    """The commutator function for quantum mechanics.

    This function behaves as such: [A, B] = A*B - B*A

    The arguments are ordered according to .__cmp__()

    Examples
    ========

    >>> from sympy import symbols
    >>> from sympy.physics.quantum import Commutator
    >>> A, B = symbols('A B', **{'commutative':False})
    >>> Commutator(B, A)
    -1*Commutator(A, B)

    Evaluate the commutator with .doit()

    >>> comm = Commutator(A,B); comm
    Commutator(A, B)
    >>> comm.doit()
    A*B - B*A

    References
    ==========

    http://en.wikipedia.org/wiki/Commutator
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

    def _latex_(self,printer):
        return "\\left[%s,%s\\right]" % tuple([
            printer._print(arg) for arg in self.args])


class AntiCommutator(Expr):
    """The commutator function for quantum mechanics.

    This function behaves as such: [A, B] = A*B - B*A

    The arguments are ordered according to .__cmp__()

    Examples
    ========

    >>> from sympy import symbols
    >>> from sympy.physics.quantum import Commutator
    >>> A, B = symbols('A B', **{'commutative':False})
    >>> Commutator(B, A)
    -1*Commutator(A, B)

    Evaluate the commutator with .doit()

    >>> comm = Commutator(A,B); comm
    Commutator(A, B)
    >>> comm.doit()
    A*B - B*A

    References
    ==========

    http://en.wikipedia.org/wiki/Commutator
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

    def _latex_(self,printer):
        return "\\left{%s,%s\\right}"%tuple([
            printer._print(arg) for arg in self.args])


#-----------------------------------------------------------------------------
# Tensor product
#-----------------------------------------------------------------------------

def matrix_tensor_product(*matrices):
    """Compute the tensor product of a sequence of sympy Matrices."""
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
        if length > 1:
            s = '(%s)' % s
        return s

    def _pretty(self, printer, *args):
        length = len(self.args)
        pform = printer._print('', *args)
        for i in range(length):
            next_pform = printer._print(self.args[i], *args)
            if isinstance(self.args[i], (Add, Pow, Mul)):
                next_pform = prettyForm(
                    *next_pform.parens(left='(', right=')')
                )
            pform = prettyForm(*pform.right(next_pform))
            if i != length-1:
                pform = prettyForm(*pform.right(u'\u2a02' + u' '))
        if length > 1:
            pform = prettyForm(*pform.parens(left='(', right=')'))
        return pform

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


def tpsimp_Mul(e):
    # TODO: This won't work with Muls that have other composites of 
    # TensorProducts, like an Add, Pow, Commutator, etc. We need to move
    # to the full parallel subs approach.
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


def tpsimp(e, **hints):
    """Try to simplify and combine TensorProducts."""
    if isinstance(e, Add):
        return Add(*[tpsimp(arg) for arg in e.args])
    elif isinstance(e, Pow):
        return tpsimp(e.base)**e.exp
    elif isinstance(e, Mul):
        return tpsimp_Mul(e)
    elif isinstance(e, (Commutator,)):
        return Commutator(*[tpsimp(arg) for arg in e.args])
    else:
        return e

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------


def represent(expr, basis, **options):
    """Represent the quantum expression in the given basis.
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

hbar = HBar()

