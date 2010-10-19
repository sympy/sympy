from sympy import Expr, sympify, Add, Mul, Function, S, Matrix, Pow, Integer
from sympy.printing.pretty.stringpict import prettyForm
from sympy.core.numbers import NumberSymbol
import sympy.mpmath.libmp as mlib

from sympy.physics.qexpr import (
    QuantumError, QExpr, split_commutative_parts, split_qexpr_parts
)
from sympy.physics.hilbert import HilbertSpace, HilbertSpaceError

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
    'tpsimp'
]

#-----------------------------------------------------------------------------
# Abstract base classes
#-----------------------------------------------------------------------------

class Representable(object):
    """An object that can be represented.
    """

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

#-----------------------------------------------------------------------------
# States, bras and kets.
#-----------------------------------------------------------------------------

class StateBase(QExpr, Representable):
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
        return None

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

    def _sympyrepr(self, printer, *args):
        return '%s(%s)' % (
            self.__class__.__name__, self._print_contents_repr(printer, *args)
        )

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
        if isinstance(other, self.dual_class):
            return OuterProduct(self, other)
        else:
            return Expr.__mul__(self, other)

    def __rmul__(self, other):
        """other*KetBase"""
        if isinstance(other, self.dual_class):
            return InnerProduct(other, self)
        else:
            return Expr.__rmul__(self, other)

    def apply_operator(self, op):
        if not isinstance(op, Operator):
            raise TypeError('Operator expected, got: %r' % op)
        apply_method = '_apply_operator_%s' % op.__class__.__name__
        if hasattr(self, apply_method):
            f = getattr(self, apply_method)
            result = f(op)
            if result is not None:
                return result
        raise NotImplementedError(
            "Don't know how apply operator %r to Ket %r" % (op, self)
        )


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
        if isinstance(other, self.dual_class):
            return InnerProduct(self, other)
        else:
            return Expr.__mul__(self, other)

    def __rmul__(self, other):
        """other*BraBase"""
        if isinstance(other, self.dual_class):
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

class Operator(QExpr, Representable):
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

    def _eval_commutator(self, other):
        """Evaluate [self, other] if known, return None if not known."""
        return None

    def _eval_rcommutator(self, other):
        """Evaluate [other, self] if known, return None if not known."""
        return None

    def _eval_anticommutator(self, other):
        """Evaluate [self, other] if known."""
        return None

    def _eval_ranticommutator(self, other):
        """Evaluate [other, self] if known."""
        return None

    #-------------------------------------------------------------------------
    # Operator application
    #-------------------------------------------------------------------------

    def apply_to_ket(self, ket):
        if not isinstance(ket, KetBase):
            raise TypeError('KetBase expected, got: %r' % ket)
        apply_method = '_apply_to_ket_%s' % ket.__class__.__name__
        if hasattr(self, apply_method):
            f = getattr(self, apply_method)
            result = f(ket)
            if result is not None:
                return result
        raise NotImplementedError(
            "Don't know how apply operator %r to Ket %r" % (self, ket)
        )

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
        if not ket.hilbert_space == bra.hilbert_space:
            raise HilbertSpaceError(
                'Incompatible hilbert spaces: %r and %r' % (ket, bra)
            )
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

    def doit(self,**kw_args):
        return self

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

class InnerProduct(Expr, Representable):
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
        if not ket.dual_class == bra.__class__:
            raise TypeError(
                'bra and ket are not dual classes: %r, %r' % \
                (bra.__class__, ket.__class__)
            )
        if not ket.hilbert_space == bra.hilbert_space:
            raise HilbertSpaceError(
                'Incompatible hilbert spaces: %r and %r' % (ket, bra)
        )
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
        r = self.bra._eval_innerproduct(self.ket, **hints)
        if r is None:
            r = self.ket._eval_innerproduct(self.bra, **hints)
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
        # TODO: apply this recursively to handle mutiple rules at once.

        A = self.args[0].expand(**hints)
        B = self.args[1].expand(**hints)

        # [A+B,C]  ->  [A,C] + [B,C]
        if isinstance(A, Add):
            return Add(*[Commutator(term,B) for term in A.args])
        if isinstance(B, Add):
            return Add(*[Commutator(A,term) for term in B.args])

        if isinstance(A, Mul):
            # [a*b,c] -> a*[b,c] + [a,c]*b
            a = A.args[0]
            b = Mul(*A.args[1:])
            c = B
            first = Mul(a, Commutator(b, c))
            second = Mul(Commutator(a, c), b)
            return Add(first, second)
        if isinstance(B, Mul):
            # [a,b*c] -> [a,b]*c + b*[a,c]
            a = A
            b = B.args[0]
            c = Mul(*B.args[1:])
            first = Mul(Commutator(a, b), c)
            second = Mul(b, Commutator(a, c))
            return Add(first, second)

        # No changes, so return self
        return self

    def doit(self, **hints):
        A = self.args[0]
        B = self.args[1]
        if isinstance(A, Operator) and isinstance(B, Operator):
            comm = A._eval_commutator(B)
            if comm is None:
                comm = B._eval_rcommutator(A)
            if comm is not None:
                return comm.doit(**hints)
        return (A*B - B*A).doit(**hints)

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
            comm = A._eval_anticommutator(B)
            if comm is None:
                comm = B._eval_ranticommutator(A)
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

class TensorProduct(Expr):

    def __new__(cls, *args):
        matrix = False
        new_args = []
        for arg in args:
            if isinstance(arg, Matrix):
                matrix = True
                new_args.append(arg)
            else:
                if matrix:
                    raise TypeError('Matrix expected, got: %r' % arg)
                else:
                    new_args.append(sympify(arg))
        if not matrix:
            c_part, new_args = cls.flatten(new_args)
            c_part = Mul(*c_part)
        else:
            c_part = 1
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
        string = ''
        for i in range(length):
            if isinstance(self.args[i], (Add, Pow)):
                string = string + '('
            string = string + sstr(self.args[i])
            if isinstance(self.args[i], (Add, Pow)):
                string = string + ')'
            if i != length-1:
                string = string + 'x'
        return string

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
        # TODO: doit doesn't work with a Mul of TensorProducts of Matrices.
        # Possibly want to switch to using expand instead?
        new_args = self.args
        if isinstance(new_args[0], Matrix):
            # Pull out the first element in the product.
            matrix_expansion  = new_args[-1]
            # Do the tensor product working from right to left.
            for mat in reversed(new_args[:-1]):
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
        else:
            new_args = [arg.doit(**hints) for arg in self.args]
            return TensorProduct(*new_args)

def tpsimp_Mul(e):
    # TODO: This won't work with Muls that have other composites of 
    # TensorProducts, like an Add, Pow, Commutator, etc. We need to move
    # to the full parallel subs approach.
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


def tpsimp(e):
    # e = e.expand()
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

    muls = e.atoms(Mul)
    subs_list = [(m,tpsimp_Mul(m)) for m in iter(muls)]
    return e.subs(subs_list)

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------

def represent(expr, basis, **options):
    """Represent the quantum expression in the given basis.
    """
    # from sympy.physics.qadd import QAdd
    # from sympy.physics.qmul import QMul
    # from sympy.physics.qpow import QPow
    if isinstance(expr, Representable):
        return expr.represent(basis, **options)
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
    elif not isinstance(expr, Mul):
        return expr

    if not isinstance(expr, Mul):
        raise TypeError('Mul expected, got: %r' % expr)

    result = S.One
    for arg in reversed(expr.args):
        result = represent(arg, basis, **options)*result
    return result


# def apply_Mul(m):
#     """
#     Take a Mul instance with operators and apply them to states.
# 
#     This method applies all operators with integer state labels
#     to the actual states.  For symbolic state labels, nothing is done.
#     When inner products of FockStates are encountered (like <a|b>),
#     the are converted to instances of InnerProduct.
# 
#     This does not currently work on double inner products like,
#     <a|b><c|d>.
# 
#     If the argument is not a Mul, it is simply returned as is.
#     """
#     if not isinstance(m, Mul):
#         return m
#     c_part, nc_part = split_commutative_parts(m)
#     n_nc = len(nc_part)
#     if n_nc == 0 or n_nc == 1:
#         return m
#     else:
#         last = nc_part[-1]
#         next_to_last = nc_part[-2]
#         if isinstance(last, KetBase):
#             if isinstance(next_to_last, Operator):
#                 try:
#                     result = next_to_last.apply_to_ket(last)
#                 except NotImplementedError:
#                     try:
#                         result = last.apply_operator(next_to_last)
#                     except NotImplementedError:
#                         return m
#                 if result == 0:
#                     return 0
#                 else:
#                     return apply_Mul(Mul(*(c_part+nc_part[:-2]+[result])))
#             elif isinstance(next_to_last, Pow):
#                 if isinstance(next_to_last.base, Operator) and \
#                     next_to_last.exp.is_Integer:
#                     result = last
#                     op = next_to_last.base
#                     for i in range(next_to_last.exp):
#                         try:
#                             result = op.apply_to_ket(result)
#                         except NotImplementedError:
#                             try:
#                                 result = result.apply_operator(op)
#                             except NotImplementedError:
#                                 return m
#                         if result == 0: break
#                     if result == 0:
#                         return 0
#                     else:
#                         return apply_Mul(Mul(*(c_part+nc_part[:-2]+[result])))
#                 else:
#                     return m
#             elif isinstance(next_to_last, last.dual_class):
#                 result = InnerProduct(next_to_last, last)
#                 if result == 0:
#                     return 0
#                 else:
#                     return apply_Mul(Mul(*(c_part+nc_part[:-2]+[result])))
#             else:
#                 return m
#         else:
#             return m
# 
# 
# def apply_operators(e):
#     """
#     Take a sympy expression with operators and states and apply the operators.
#     """
#     e = e.expand()
#     muls = e.atoms(Mul)
#     subs_list = [(m,apply_Mul(m)) for m in iter(muls)]
#     return e.subs(subs_list)
#     
#     
# def apply_operators(e):
# 
#     # if all we have is a Qbit without any gates, return the qbit
#     if isinstance(circuit, Qbit):
#         return circuit
# 
#     #if we have a Mul object, get the state of the system
#     elif isinstance(circuit, QMul):
#         states = circuit.args[len(circuit.args)-1]
#         states = states.expand()
# 
#     #if we have an add object with gates mixed in, apply_gates recursively
#     elif isinstance(circuit, QAdd):
#         result = 0
#         for i in circuit.args:
#             result = result + apply_gates(i, basis, floatingPoint)
#         if floatingPoint:
#             result = result.evalf()
#         return result
# 
#     state_coeff = 1
#     #pick out each object that multiplies the state
#     for multiplier in reversed(circuit.args[:len(circuit.args)-1]):
# 
#         #if the object that mutliplies is a Gate, we will apply it once
#         if isinstance(multiplier, Gate):
#             gate = multiplier
#             number_of_applications = 1
# 
#         #if the object that multiplies is a Pow who's base is a Gate, we will
#         # apply Pow.exp times
#         elif isinstance(multiplier, QPow) and isinstance(multiplier.base, Gate):
#             gate = multiplier.base
#             number_of_applications = multiplier.exp
# 
#         #if the object that multiplies is not a gate of any sort,
#         # we apply it by multiplying
#         else:
#             state_coeff = multiplier*state_coeff
#             continue
# 
#         #if states is in superposition of states (a sum of qbits states),
#         # applyGates to each state contined within
#         if isinstance(states, QAdd):
#             #special check for non-distributivity, do all at once
#             if isinstance(gate, NondistributiveGate):
#                 states = gate.measure(states)
#             else:
#                 result = 0
#                 for state in states.args:
#                     result = result + \
#                     apply_gates(gate**number_of_applications*state, basis,\
#                     floatingPoint)
#                 if floatingPoint:
#                     result = result.evalf()
#                 states = result
#                 states = states.expand()
# 
#         #if we have a mul, apply gate to each register and multiply result
#         elif isinstance(states, QMul):
#             #find the Qbits in the Mul
#             for i in range(len(states.args)):
#                 if isinstance(states.args[i],Qbit):
#                     break
#             #if we didn't find one, something is wrong
#             if not isinstance(states.args[i],Qbit):
#                 print states
#                 raise QuantumError()
# 
#             #apply the gate the right number of times to this state
#             coefficient = Mul(*(states.args[:i]+states.args[i+1:]))#TODO
#             states = apply_gates(gate**(number_of_applications)*states.args[i],\
#             basis, floatingPoint)
#             states = coefficient*states
#             states = states.expand()
# 
#         #If we have a single Qbit, apply to this Qbit
#         elif isinstance(states, Qbit):
#             if isinstance(gate, NondistributiveGate):
#                 states = gate.measure(states)
#             else:
#                 basis_name = basis.__class__.__name__
#                 apply_method_name = '_apply_%s' % basis_name
#                 apply_method = getattr(gate, apply_method_name)
#                 states = apply_method(states)
#                 states = states.expand()
#                 number_of_applications -= 1
#                 while number_of_applications > 0:
#                     states = apply_gates(gate*states, basis, floatingPoint)
#                     number_of_applications -= 1
# 
#         #if it's not one of those, there is something wrong
#         else:
#             raise QuantumError()
# 
#     #tack on any coefficients that were there before and simplify
#     states = state_coeff*states
#     if isinstance(states, (Add,Pow, Mul)):
#         states = states.expand()
#     return states


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

