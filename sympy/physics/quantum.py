from sympy import (
    Expr, Basic, sympify, Add, Mul, Pow, 
    I, Function, Integer, S, sympify, Matrix, oo
)
from sympy.core.sympify import _sympify
from sympy.core.decorators import call_highest_priority
from sympy.physics.hilbert import *
from sympy.core.numbers import Number
from sympy.core.symbol import Symbol, symbols

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
    'StateBase',
    'State',
    'Ket',
    'Bra',
    'InnerProduct',
    'OuterProduct',
    'Operator',
    'Dagger',
    'KroneckerDelta',
    'Commutator',
    # 'AntiCommutator',
]


#-----------------------------------------------------------------------------
# Error handling
#-----------------------------------------------------------------------------

class QuantumError(Exception):
    pass


#-----------------------------------------------------------------------------
# Main objects
#-----------------------------------------------------------------------------

class Representable(object):
    """An object that can be represented."""

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

class StateBase(Expr, Representable):

    # Slots are for instance variables that can always be computed dynamically
    # from self.args, but that we don't want users to have to pass each time.
    __slots__ = ['hilbert_space']

    # Class level attributes that are the same for all instances go here.
    # Because the hilbert_space can change, it can't go here. All instance
    # level things must go in self.args.
    is_continuous = False
    is_discrete = False
    # I am a little worried that the basis set might need to be in self.args.
    # Or maybe this needs to be a slot if it can be computed from self.args.
    basis_set = None

    _op_priority = 100.0

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

    def _sympyrepr(self, printer, *args):
        return '%s(%s)' % (self.__class__.__name__, self._print_name(printer, *args))

    def _sympystr(self, printer, *args):
        return '%s%s%s' % (self.lbracket, self._print_name(printer, *args), self.rbracket)

    def _pretty(self, printer, *args):
        from sympy.printing.pretty.stringpict import prettyForm
        pform = self._print_name_pretty(printer, *args)
        pform = prettyForm(*pform.left(prettyForm(self.lbracket)))
        pform = prettyForm(*pform.right(prettyForm(self.rbracket)))
        return pform

    def _print_name(self, printer, *args):
        return printer._print(self.name, *args)

    def _print_name_pretty(self, printer, *args):
        pform = printer._print(self.name, *args)
        return pform

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
        compare_hilbert(self, other)
        _validate_add(self, other)
        return Add(self, other)

    @call_highest_priority('__add__')
    def __radd__(self, other):
        compare_hilbert(other, self)
        _validate_add(other, self)
        return Add(other, self)

    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        raise NotImplementedError("Can't do Tensor Products of States Yet")

    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        raise QuantumError("Can't raise %s to a State" % (other.__class__.__name__,))

    def _sympyrepr(self, printer, *args):
        return '%s(%s)' % (self.__class__.__name__, self._print_name(printer, *args))

    def _sympystr(self, printer, *args):
        return '%s%s%s' % (self.lbracket, self._print_name(printer, *args), self.rbracket)

    def _pretty(self, printer, *args):
        from sympy.printing.pretty.stringpict import prettyForm
        pform = self._print_name_pretty(printer, *args)
        pform = prettyForm(*pform.left(prettyForm(self.lbracket)))
        pform = prettyForm(*pform.right(prettyForm(self.rbracket)))
        return pform

class KetBase(StateBase):

    lbracket = '|'
    rbracket = '>'

    @property
    def dual(self):
        return BraBase(*self.args)

class BraBase(StateBase):

    lbracket = '<'
    rbracket = '|'

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

    @property
    def dual(self):
        return Bra(*self.args)

class Bra(State, BraBase):

    @property
    def dual(self):
        return Ket(*self.args)

class TimeDepState(StateBase):

    def __new__(cls, name, time):
        # First compute args and call Expr.__new__ to create the instance
        name = cls._eval_name(name, time)
        time = cls._eval_time(name, time)
        inst = Expr.__new__(cls, name, time, **{'commutative':False})
        # Now set the slots on the instance
        inst.hilbert_space = cls._eval_hilbert_space(name, time)
        return inst

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
    """A basis set for a Hilbert Space"""

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

class Operator(Expr, Representable):
    """
    Base class for non-commuting Quantum operators.

    >>> from sympy import simplify
    >>> from sympy.physics.quantum import Operator
    >>> A = Operator('A')
    >>> print A
    A
    """
    _op_priority = 100.0

    __slots__ = ['hilbert_space']

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
        #if not isinstance(other, (Mul, Add, Pow, Number, Symbol)):
        return Pow(self, other)

    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        #??operator**operator??
        return Pow(other, self)

    def doit(self,**kw_args):
        return self

    def _sympyrepr(self, printer, *args):
        return '%s(%s)' % (self.__class__.__name__, printer._print(self.name, *args))

    def _sympystr(self, printer, *args):
        return printer._print(self.name, *args)

    def _pretty(self, printer, *args):
        return printer._print(self.name, *args)

class InnerProduct(Expr):
    """
    An unevaluated inner product between a Bra and Ket.
    """

    def __new__(cls, bra, ket):
        #What about innerProd(1,1), should it auto simplify?
        if not (bra and ket):
            raise Exception('InnerProduct requires a leading Bra and a trailing Ket')
        assert issubclass(bra.evaluates, Bra), 'First argument must be a Bra'
        assert issubclass(ket.evaluates, Ket), 'Second argument must be a Ket'
        r = cls.eval(bra, ket)
        if isinstance(r, Expr):
            return r
        obj = Expr.__new__(cls, bra, ket)
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
        return Number

    def _eval_dagger(self):
        return InnerProduct(Dagger(self.ket), Dagger(self.bra))

    def _sympyrepr(self, printer, *args):
        return '%s(%s,%s)' % (self.__class__.__name__, printer._print(self.bra, *args), printer._print(self.ket, *args))

    def _sympystr(self, printer, *args):
        sbra = str(self.bra)
        sket = str(self.ket)
        return '%s|%s' % (sbra[:-1], sket[1:])

    def _pretty(self, printer, *args):
        return self.bra._pretty(printer, *args)*self.ket._pretty(printer, *args)


class OuterProduct(Expr):
    """
    An unevaluated outer product between a Ket and Bra.
    """

    __slots__ = ['hilbert_space']

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
        return self.ket._pretty(printer, *args)*self.bra._pretty(printer, *args)


class Dagger(Expr):
    """
    General hermitian conjugate operation.
    """

    def __new__(cls, arg):
        if isinstance(arg, Matrix):
            return cls.eval(arg)
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
        try:
            d = arg._eval_dagger()
        except:
            if isinstance(arg, Expr):
                if arg.is_Add:
                    return Add(*tuple(map(Dagger, arg.args)))
                if arg.is_Mul:
                    return Mul(*tuple(map(Dagger, reversed(arg.args))))
                if arg.is_Number:
                    return arg
                if arg.is_Pow:
                    return Pow(Dagger(arg.args[0]), Dagger(arg.args[1]))
                if arg == I:
                    return -arg
            elif isinstance(arg, Matrix):
                arg = arg.T
                for i in range(arg.rows*arg.cols):
                    arg[i] = Dagger(arg[i])
                return arg
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


class KroneckerDelta(Function):
    """
    Discrete delta function.
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
    >>> from sympy.physics.secondquant import Commutator
    >>> A, B = symbols('A B', commutative=False)
    >>> Commutator(B, A)
    Commutator(B, A)

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
    """Represent the quantum expression in the given basis."""
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

def _evaluate_type(sympy_binop):
    #determines if a Mul or Add evaluates to a Ket, Bra, Operator, or Number
    #Willbe Used in _validate_add; this is absurd, We should just have an extra thing in Mul that says "I am a ket"
    #Mul is not extendable at all. When you can't implement Matrix operations in your data and bin-op model, your model is wrong.
    #But I digress
    if isinstance(sympy_binop, Add):
        return _evaluate_type(sympy_binop.args[0])
    elif isinstance(sympy_binop, Mul):
        result = [Number, '']
        for item in sympy_binop.args:
            if issubclass(result[0], Number):
                result = [item.__class__, ''] #FIXME 
            elif issubclass(result[0], Operator):
                if isinstance(item, (Number, Operator)):
                    result = [Operator, '']
                elif isinstance(item, KetBase):
                    return [State, 'ket']
                else:
                    raise QuantumError("Can't multiply %s and %s" % (result.__name__, item.__class__.__name__))
            elif issubclass(result[0], State):
                if result[1] == 'ket':
                    if isinstance(item, BraBase):
                        return [Number, '']
                    elif isinstance(item, Number):
                        return 
                else:
                    pass                
    #base case
    elif isinstance(sympy_binop, (OuterProduct, Operator)):
        return Operator
    elif isinstance(sympy_binop, (InnerProduct, Number)):
        return Number  
    elif isinstance(sympy_binop, State):
        return sympy_binop.__class__
    else:
        raise QuantumError("Don't Know how you got here. Contact your system administrator?")

def _validate_add(expr1, expr2):
    if isinstance(expr1, Add):
        _validate_add(expr1.args[-1], expr2)
        return 
    elif isinstance(expr2, Add):
        _validate_add(expr1, expr2.args[0])
        return
    else:
        if isinstance(expr1, StateBase) and isinstance(expr2, StateBase):
            if expr1.__class__ == expr2.__class__:
                return
        elif isinstance(expr1, Operator) and isinstance(expr2, Operator):
            return
        elif expr1 == 0 or expr2 == 0:
            return
        else:
            raise QuantumError("Can't add %s and %s" % (expr1.__class__.__name__, expr2.__class__.__name__))

def _validate_mul(expr1, expr2):
    if isinstance(expr1, Mul):
        _validate_mul(expr1.args[-1], expr2)
    elif isinstance(expr2, Mul):
        _validate_mul(expr1, expr2.args[0])
    elif isinstance(expr1, Add):
        _validate_mul(expr1.args[-1], expr2)
    elif isinstance(expr2, Add):
        _validate_mul(expr1, expr2.args[-1])
    elif isinstance(expr1, State) and isinstance(expr2, StateBase):
        if isinstance(expr1, Ket) and isinstance(expr2, KetBase):
            raise NotImplementedError("TensorProducts of ket%ket not implemented")
        if isinstance(expr1, Bra) and isinstance(expr2, BraBase):
            raise NotImplementedError("TensorProducts of bra%bra not implemented")
    elif isinstance(expr1, (Operator, OuterProduct)) and isinstance(expr2, BraBase):
        raise QuantumError('(Operator or OuterProduct)*Bra is invalid.\n(Try using parentheses to form inner and outer products)')
    elif isinstance(expr2, (Operator, OuterProduct)) and isinstance(expr1, Ket):
        raise QuantumError('Ket*(Operator or OuterProduct) is invalid.\n(Try using parentheses to form inner and outer products)')

def _qmul(expr1, expr2):
    """
    Check to see if arg1 or arg2 can combine to be a valid Mul.
    """
    if isinstance(expr1, Bra) and isinstance(expr2, Ket):
        return InnerProduct(expr1, expr2)
    elif isinstance(expr1, Ket) and isinstance(expr2, Bra):
        return OuterProduct(expr1, expr2)
    else:
        return Mul(expr1, expr2)

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
