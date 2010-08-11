from sympy import (
    Expr, Basic, sympify, Add, Mul, Pow, 
    I, Function, Integer, S, sympify, Matrix, oo
)
from sympy.core.sympify import _sympify
from sympy.core.numbers import Number
from sympy.core.symbol import Symbol, symbols
from sympy.physics.quantumbasic import QuantumError, QuantumBasic
from sympy.printing.pretty.stringpict import prettyForm
from sympy.physics.hilbert import HilbertSpace

"""
Questions:

* What is an appropriate base class for Operator?
* What does doit do and should be use it to do things like apply operators?
* How should we handle assumptions?  New or old system?
* How should we handle different types of operators like symmetric, hermitian,
  etc? Assumptions?
* How do we handle the null state?  Do we use Integer(0) or an actual
  ZeroState Singleton?
* Should we allow states to have names that are composite from the ground up?
* Should __repr__ return a nice representation of things or the more pythonic
  thing.  |a> or Ket('a').
* What should a basis set consist of?  Somehow it needs to know about a set
  of states and maybe the operator they are eigenvectors of. We need a way
  of labeling these states by name and index: |a[0]>, |a[1]>, etc.
* Should Operator.name and State.name be a string or symbol?
"""

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
    'InnerProduct',
    'OuterProduct',
    'Operator',
    'Dagger',
    'KroneckerDelta',
    'Commutator',
    # 'AntiCommutator',
]

#-----------------------------------------------------------------------------
# Main objects
#-----------------------------------------------------------------------------

class Representable(object):
    """
    An object that can be represented.
    """

    def represent(self, basis, **options):
        rep_method = '_represent_%s' % basis.__class__.__name__
        if hasattr(self, rep_method):
            f = getattr(self, rep_method)
            rep = f(basis, **options)
            return rep
        else:
            raise QuantumError(
                'Object %r does not know how to represent itself in basis: %r' % (self, basis)
            )

class StateBase(QuantumBasic, Representable):
    """
    Base class for general abstract states in quantum mechanics.
    """

    # Slots are for instance variables that can always be computed dynamically
    # from self.args, but that we don't want users to have to pass each time.
    # Class level attributes that are the same for all instances go here.
    # Because the hilbert_space can change, it can't go here. All instance
    # level things must go in self.args.
    is_continuous = False
    is_discrete = False
    # I am a little worried that the basis set might need to be in self.args.
    # Or maybe this needs to be a slot if it can be computed from self.args.
    basis_set = None

    @property
    def name(self):
        raise NotImplementedError('name must be implemented in a subclass')

    @property
    def dual(self):
        raise NotImplementedError('dual must be implemented in a subclass')

    def _eval_dagger(self):
        return self.dual

    @property
    def is_symbolic(self):
        return True

    @property
    def evaluates(self):
        return self.__class__

    def _print_name(self, printer, *args):
        return printer._print(self.name, *args)

    def _print_name_pretty(self, printer, *args):
        pform = printer._print(self.name, *args)
        return pform

    def _sympyrepr(self, printer, *args):
        return '%s(%s)' % (self.__class__.__name__, self._print_name(printer, *args))

    def _sympystr(self, printer, *args):
        return '%s%s%s' % (self.lbracket, self._print_name(printer, *args), self.rbracket)

    def _pretty(self, printer, *args):
        from sympy.printing.pretty.stringpict import prettyForm
        pform = self._print_name_pretty(printer, *args)
        pform = prettyForm(*pform.left((self.lbracketPretty)))
        pform = prettyForm(*pform.right((self.rbracketPretty)))
        return pform

class KetBase(StateBase):
    """
    Base class for Ket states.
    """

    lbracket = '|'
    rbracket = '>'
    lbracketPretty = prettyForm(u'\u2758')
    rbracketPretty = prettyForm(u'\u276D')
 

    @property
    def dual(self):
        return BraBase(*self.args)

class BraBase(StateBase):
    """
    Base class for Bra states.
    """

    lbracket = '<'
    rbracket = '|'
    lbracketPretty = prettyForm(u'\u276C')
    rbracketPretty = prettyForm(u'\u2758') 

    @property
    def dual(self):
        return KetBase(*self.args)

class State(StateBase):
    """
    General abstract quantum state.

    Anywhere you can have a State, you can also have Integer(0).
    All code must check for this!
    """

    def __new__(cls, name):
        # First compute args and call Expr.__new__ to create the instance
        name = cls._eval_name(name)
        inst = Expr.__new__(cls, name, **{'commutative':False})
        # Now set the slots on the instance
        inst.hilbert_space = cls._eval_hilbert_space(name)
        return inst

    @classmethod
    def _eval_name(cls, name):
        return sympify(name)

    @property
    def name(self):
        return self.args[0]

    @classmethod
    def _eval_hilbert_space(cls, name):
        return HilbertSpace()

    def doit(self, **kw_args):
        return self

class Ket(State, KetBase):
    """
    A class for creating Ket objects in sympy. Inherits from State and KetBase.
    In a state represented by a ray in Hilbert space, a Ket is a vector that
    that points along that ray [1].

    Examples
    ========

    Creating and using a Ket:
        
        >>> from sympy.physics.quantum import Ket
        >>> psi = Ket('psi')
        >>> psi
        |psi>
        >>> psi.name
        psi
        >>> psi.dual
        <psi|


    References
    ==========

    [1] http://en.wikipedia.org/wiki/Bra-ket_notation
    """

    @property
    def dual(self):
        return Bra(*self.args)

class Bra(State, BraBase):
    """
    A class for creating Bra objects in sympy. Inherits from State and BraBase.
    A Bra is the dual of a Ket [1].

    Examples
    ========

    Creating and using a Bra:
        
        >>> from sympy.physics.quantum import Bra
        >>> b = Bra('bus')
        >>> b
        <bus|
        >>> b.name
        bus
        >>> b.dual
        |bus>


    References
    ==========

    [1] http://en.wikipedia.org/wiki/Bra-ket_notation
    """

    @property
    def dual(self):
        return Ket(*self.args)

class TimeDepState(StateBase):
    """
    General abstract time dependent quantum state.
    """
    def __new__(cls, name, time):
        # First compute args and call Expr.__new__ to create the instance
        name = cls._eval_name(name, time)
        time = cls._eval_time(name, time)
        inst = Expr.__new__(cls, name, time, **{'commutative':False})
        # Now set the slots on the instance
        inst.hilbert_space = cls._eval_hilbert_space(name, time)
        return inst

    @property
    def name(self):
        return self.args[0]

    @classmethod
    def _eval_name(cls, name, time):
        return sympify(name)

    @classmethod
    def _eval_time(cls, name, time):
        return sympify(time)

    @classmethod
    def _eval_hilbert_space(cls, name, time):
        return HilbertSpace()

class TimeDepKet(TimeDepState, KetBase):

    @property
    def dual(self):
        return TimeDepBra(*self.args)

class TimeDepBra(TimeDepState, BraBase):

    @property
    def dual(self):
        return TimeDepKet(*self.args)

class BasisSet(Basic):
    """
    A basis set for a Hilbert space.
    """

    def __new__(cls, dimension):
        dimension = sympify(dimension)
        return Basic.__new__(cls, dimension)

    @property
    def dimension(self):
        return self.args[0]

    def __len__(self):
        return self.dimension

    @property
    def unity(self):
        """Return the unity operator for this basis set."""
        raise NotImplementedError('Not implemented')

class Operator(QuantumBasic, Representable):
    """
    Base class for non-commuting quantum operators. An operator is a map from
    one vector space to another [1]. In quantum mechanics, operators correspond
    to observables [2].

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

    def __new__(cls, name):
        name = sympify(name)
        obj = Expr.__new__(cls, name, **{'commutative': False})
        obj.hilbert_space = cls._eval_hilbert_space(name)        
        return obj

    @classmethod
    def _eval_hilbert_space(cls, name):
        return HilbertSpace()

    @property
    def name(self):
        return self.args[0]

    @property
    def evaluates(self):
        return self.__class__

    def doit(self,**kw_args):
        return self

    def _sympyrepr(self, printer, *args):
        return '%s(%s)' % (self.__class__.__name__, printer._print(self.name, *args))

    def _sympystr(self, printer, *args):
        return printer._print(self.name, *args)

    def _pretty(self, printer, *args):
        return printer._print(self.name, *args)

class InnerProduct(QuantumBasic):
    """
    An unevaluated inner product between a Bra and a Ket. Because a Bra is 
    essentially a row vector and a Ket is essentially a column vector, the
    inner product evaluates to a complex number.

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

    def __new__(cls, bra, ket):
        #What about innerProd(1,1), should it auto simplify?
        if not (bra and ket):
            raise Exception('InnerProduct requires a leading Bra and a trailing Ket')
        assert issubclass(bra.evaluates, Bra), 'First argument must be a Bra'
        assert issubclass(ket.evaluates, Ket), 'Second argument must be a Ket'
        assert bra.hilbert_space == ket.hilbert_space
        r = cls.eval(bra, ket)
        if isinstance(r, Expr):
            return r
        obj = Expr.__new__(cls, bra, ket)
        obj.hilbert_space = bra.hilbert_space
        return obj

    @classmethod
    def eval(cls, bra, ket):
        # We need to decide what to do here. We probably will ask if the
        # bra and ket know how to do the inner product.
        return None

    @property
    def bra(self):
        return self.args[0]

    @property
    def ket(self):
        return self.args[1]

    @property
    def evaluates(self):
        return self.__class__

    def _eval_dagger(self):
        return InnerProduct(Dagger(self.ket), Dagger(self.bra))

    def _sympyrepr(self, printer, *args):
        return '%s(%s,%s)' % (self.__class__.__name__, printer._print(self.bra, *args), printer._print(self.ket, *args))

    def _sympystr(self, printer, *args):
        sbra = str(self.bra)
        sket = str(self.ket)
        return '%s|%s' % (sbra[:-1], sket[1:])

    def _pretty(self, printer, *args):
        pform = prettyForm(u'\u276C')
        pform = prettyForm(*pform.right(self.bra._print_name(printer, *args)))
        return prettyForm(*pform.right(self.ket._pretty(printer, *args)))

class OuterProduct(Operator):
    """
    An unevaluated outer product between a Ket and Bra. Because a Ket is
    essentially a column vector and a Bra is essentially a row vector, the
    outer product evaluates to an operator (i.e. matrix) [1].

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

    def __new__(cls, ket, bra):
        if not (ket and bra):
            raise Exception('OuterProduct requires a leading Ket and a trailing Bra')
        assert issubclass(ket.evaluates, Ket), 'First argument must be a Ket'
        assert issubclass(bra.evaluates, Bra), 'Second argument must be a Bra'
        assert ket.hilbert_space == bra.hilbert_space
        r = cls.eval(ket, bra)
        if isinstance(r, Expr):
            return r
        obj = Expr.__new__(cls, *(ket, bra), **{'commutative': False})
        obj.hilbert_space = ket.hilbert_space
        return obj

    @classmethod
    def eval(cls, ket, bra):
        # We need to decide what to do here. We probably will ask if the
        # bra and ket know how to do the outer product.
        return None

    @property
    def ket(self):
        return self.args[0]

    @property
    def bra(self):
        return self.args[1]

    @property
    def evaluates(self):
        return self.__class__

    def _eval_dagger(self):
        return OuterProduct(Dagger(self.bra), Dagger(self.ket))

    def _sympyrepr(self, printer, *args):
        return '%s(%s,%s)' % (self.__class__.__name__, printer._print(self.ket, *args), printer._print(self.bra, *args))

    def _sympystr(self, printer, *args):
        return str(self.ket)+str(self.bra)

    def _pretty(self, printer, *args):
        pform = self.ket._pretty(printer, *args)
        return prettyForm(*pform.right(self.bra._pretty(printer, *args)))

class Dagger(QuantumBasic):
    """
    General Hermitian conjugate operation. When an object is daggered it is
    transposed and then the complex conjugate is taken [1].

    Examples
    ========

    Daggering various quantum objects:

        >>> from sympy.physics.quantum import Dagger, Ket, Bra, InnerProduct, OuterProduct, Operator
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

    def __new__(cls, arg):
        if isinstance(arg, Matrix):
            return arg.H
        arg = sympify(arg)
        r = cls.eval(arg)
        if isinstance(r, Expr):
            return r
        #make unevaluated dagger commutative or non-commutative depending on arg
        if arg.is_commutative:
            obj = Expr.__new__(cls, arg)
        else:
            obj = Expr.__new__(cls, arg, **{'commutative':False})
        return obj

    @classmethod
    def eval(cls, arg):
        """
        Evaluates the Dagger instance.
        """
        if hasattr(arg, '_eval_dagger'):
            return arg._eval_dagger()
        elif not isinstance(arg, (InnerProduct, OuterProduct, Operator, StateBase)):
            return arg.conjugate()
        elif isinstance(arg, Matrix):
            arg = arg.T
            for i in range(arg.rows*arg.cols):
                arg[i] = Dagger(arg[i])
            return arg
        else:
            return None

    @property
    def evaluates(self):
        return self.args[0].evaluates

    @property
    def hilbert_space(self):
        return self.args[0].hilbert_space

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

class KroneckerDelta(Function):
    """
    Discrete delta function. A function that takes in two integers i and j.
    It returns 0 if i and j are not equal or it returns 1 if i and j are equal.

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
            return Integer(1)
        elif diff.is_number:
            return S.Zero

    def _eval_subs(self, old, new):
        r = KroneckerDelta(self.args[0].subs(old, new), self.args[1].subs(old, new))
        return r

    def _eval_dagger(self):
        return self

    def _latex_(self,printer):
        return "\\delta_{%s%s}"% (self.args[0].name,self.args[1].name)

    def __repr__(self):
        return "KroneckerDelta(%s,%s)"% (self.args[0],self.args[1])

    def __str__(self):
        return 'd(%s,%s)'% (self.args[0],self.args[1])


class Commutator(Function):
    """
    The Commutator:  [A, B] = A*B - B*A

    The arguments are ordered according to .__cmp__()

    >>> from sympy import symbols
    >>> from sympy.physics.quantum import Commutator
    >>> A, B = symbols('A B', **{'commutative':False})
    >>> Commutator(B, A)
    -Commutator(A, B)

    Evaluate the commutator with .doit()

    >>> comm = Commutator(A,B); comm
    Commutator(A, B)
    >>> comm.doit()
    A*B - B*A
    """

    is_commutative = False
    nargs = 2

    @classmethod
    def eval(cls, a, b):
        """
        The Commutator [A,B] is on canonical form if A < B


        """
        if not (a and b): return S.Zero
        if a == b: return S.Zero
        if a.is_commutative or b.is_commutative:
            return S.Zero

        #
        # [A+B,C]  ->  [A,C] + [B,C]
        #
        a = a.expand()
        if isinstance(a,Add):
            return Add(*[cls(term,b) for term in a.args])
        b = b.expand()
        if isinstance(b,Add):
            return Add(*[cls(a,term) for term in b.args])

        #
        # [xA,yB]  ->  xy*[A,B]
        #
        c_part = []
        nc_part = []
        nc_part2 = []
        if isinstance(a,Mul):
            c_part,nc_part = split_commutative_parts(a)
        if isinstance(b,Mul):
            c_part2,nc_part2 = split_commutative_parts(b)
            c_part.extend(c_part2)
        if c_part:
            a = nc_part or [a]
            b = nc_part2 or [b]
            return Mul(*c_part)*cls(Mul(*a),Mul(*b))

        #
        # Canonical ordering of arguments
        #
        if a > b:
            return S.NegativeOne*cls(b,a)

    def doit(self, **hints):
        a = self.args[0]
        b = self.args[1]
        return (a*b - b*a).doit(**hints)

    def _eval_dagger(self):
        return Commutator(Dagger(self.args[1]), Dagger(self.args[0]))

    def __repr__(self):
        return "Commutator(%s,%s)" %(self.args[0], self.args[1])

    def __str__(self):
        return "[%s,%s]" %(self.args[0], self.args[1])

    def _latex_(self,printer):
        return "\\left[%s,%s\\right]"%tuple([
            printer._print(arg) for arg in self.args])


#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------

def represent(expr, basis, **options):
    """
    Represent the quantum expression in the given basis.
    """
    from sympy.physics.qadd import QAdd
    from sympy.physics.qmul import QMul
    from sympy.physics.qpow import QPow
    if isinstance(expr, Representable):
        return expr.represent(basis, **options)
    elif isinstance(expr, QAdd):
        result = S.Zero
        for args in expr.args:
            if not result:
                result = represent(args, basis, **options)
            else:
                result += represent(args, basis, **options)
        return result
    elif isinstance(expr, QPow):
        return represent(expr.base, basis, **options)**expr.exp
    elif not isinstance(expr, QMul):
        return expr

    if not isinstance(expr, QMul):
        raise TypeError('Mul expected, got: %r' % expr)

    result = S.One
    for arg in reversed(expr.args):
        result = represent(arg, basis, **options)*result
    return result

def split_product(expr):
    """
    Separates a (valid) Mul of quantum objects and inner/outer products into a 
    Mul.

    * Only works for simple Muls (e.g. a*b*c*d) right now.
    """
    new_expr = Mul()
    for arg in expr.args:
        if isinstance(arg, InnerProduct):
            new_expr = new_expr*Mul(*arg.args)
        elif isinstance(arg, OuterProduct):
            new_expr = new_expr*Mul(*arg.args)
        else:
            new_expr = new_expr*arg
    return new_expr

def split_commutative_parts(m):
    c_part = [p for p in m.args if p.is_commutative]
    nc_part = [p for p in m.args if not p.is_commutative]
    return c_part, nc_part
