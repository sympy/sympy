from sympy import (
    Expr, Basic, sympify, Add, Mul, Pow, 
    I, Function, Integer, S, sympify, Matrix
)
from sympy.core.decorators import call_highest_priority
from sympy.physics.hilbert import *

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


class BasisSet(Basic):
    """A basis set for a Hilbert Space"""

    hilbert_space = HilbertSpace()

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


class State(Expr, Representable):
    """
    General abstract quantum state.

    Anywhere you can have a State, you can also have Integer(0).
    All code must check for this!
    """

    _op_priority = 100.0

    hilbert_space = HilbertSpace()

    def __new__(cls, name, kind='ket'):
        if not (kind=='ket' or kind=='bra'):
            raise ValueError("kind must be either 'ket' or 'bra', got: %r" % kind)
        name = sympify(name)
        obj = Expr.__new__(cls, name, kind, **{'commutative': False})
        return obj

    @property
    def name(self):
        return self.args[0]

    @property
    def kind(self):
        return self.args[1]

    @property
    def is_ket(self):
        if self.kind == 'ket':
            return True
        elif self.kind == 'bra':
            return False

    @property
    def is_bra(self):
        if self.kind == 'bra':
            return True
        elif self.kind == 'ket':
            return False

    @property
    def is_symbolic(self):
        return True

    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return _qmul(self, other)

    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return _qmul(other, self)

    @call_highest_priority('__radd__')
    def __add__(self, other):
        _validate_add(self, other)
        return Add(self, other)

    @call_highest_priority('__add__')
    def __radd__(self, other):
        _validate_add(self, other)
        return Add(other, self)

    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        raise NotImplementedError("Can't do Tensor Products of States Yet")

    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        raise QuantumError("Can't raise %s to a State" % (other.__class__.__name__,))

    def _print_name(self, printer, *args):
        return printer._print(self.args[0], *args)

    def _print_name_pretty(self, printer, *args):
        pform = printer._print(self.args[0], *args)
        return pform

    def doit(self, **kw_args):
        return self

    def _eval_dagger(self):
        return self.dual

    @property
    def dual(self):
        if self.is_ket:
            return State(self.args[0], 'bra')
        elif self.is_bra:
            return State(self.args[0], 'ket')

    @property
    def lbracket(self):
        if self.is_ket:
            return '|'
        else:
            return '<'

    @property
    def rbracket(self):
        if self.is_ket:
            return '>'
        else:
            return '|'

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

def Ket(name):
    return State(name, kind='ket')

def Bra(name):
    return State(name, kind='bra')


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

    hilbert_space = HilbertSpace()

    def __new__(cls, name):
        name = sympify(name)
        obj = Expr.__new__(cls, name, **{'commutative': False})
        return obj

    @property
    def name(self):
        return self.args[0]

    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return _qmul(self, other)

    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return _qmul(other, self)

    @call_highest_priority('__radd__')
    def __add__(self, other):
        _validate_add(self, other)
        return Add(self, other)

    @call_highest_priority('__add__')
    def __radd__(self, other):
        _validate_add(self, other)
        return Add(other, self)

    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        if not isinstance(other, (Mul, Add, Pow, Number, Symbol)):
            raise QuantumError("Can't raise Operator to %s" % (other.__class__.__name__,))
        return Pow(self, other)

    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        #??operator**operator??
        pass

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
        if not (bra and ket):
            raise Exception('InnerProduct requires a leading Bra and a trailing Ket')
        assert (isinstance(bra, State) and bra.kind == 'bra'), 'First argument must be a Bra'
        assert (isinstance(ket, State) and ket.kind == 'ket'), 'Second argument must be a Ket'
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

    def __new__(cls, ket, bra):
        if not (ket and bra):
            raise Exception('OuterProduct requires a leading Ket and a trailing Bra')
        assert (isinstance(ket, State) and ket.kind == 'ket'), 'First argument must be a Ket'
        assert (isinstance(bra, State) and bra.kind == 'bra'), 'Second argument must be a Bra'
        r = cls.eval(ket, bra)
        if isinstance(r, Expr):
            return r
        obj = Expr.__new__(cls, *(ket, bra), **{'commutative': False})
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

def is_bra(expr):
    if isinstance(expr, State):
        if expr.is_bra:
            return True
    return False

def is_ket(expr):
    if isinstance(expr, State):
        if expr.is_bra:
            return True
    return False


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

def _validate_add(expr1, expr2):
    if isinstance(expr1, Add):
        _validate_add(expr1.args[-1], expr2)
        return 
    elif isinstance(expr2, Add):
        _validate_add(expr1, expr2.args[0])
        return
    else:
        if isinstance(expr1, State) and isinstance(expr2, State):
            if expr1.kind == expr2.kind:
                return
        elif isinstance(expr1, Operator) and isinstance(expr2, Operator):
            return
        else:
            raise QuantumError("Can't add %s and %s" % (expr1.__class__.__name__, expr2.__class__.__name__))

def _qmul(expr1, expr2):
    """
    Check to see if arg1 or arg2 can combine to be a valid Mul.
    """
    def _validate(expr, mul, i):
        if isinstance(mul.args[i], State) and isinstance(expr, State):
            if mul.args[i].kind == 'bra' and expr.kind == 'bra': 
                raise NotImplementedError
            elif mul.args[i].kind == 'ket' and expr.kind == 'ket':
                raise NotImplementedError
        if (isinstance(mul.args[i], Operator) or isinstance(mul.args[i], OuterProduct)) and (isinstance(expr, State) and expr.kind == 'bra'):
            raise QuantumError('(Operator or OuterProduct)*Bra is invalid in quantum mechanics.')
        if (isinstance(mul.args[i], State) and mul.args[i].kind == 'ket') and (isinstance(expr, Operator) or isinstance(expr, OuterProduct)):
            raise QuantumError('Ket*(Operator or OuterProduct) is invalid in quantum mechanics.')
        if i == -1:
            new_args = mul.args[:i] + (_qmul(mul.args[i], expr), )
            return Mul(*new_args)
        else:
            new_args = (_qmul(expr, mul.args[i]), ) + mul.args[1:]
            return Mul(*new_args)

    if isinstance(expr1, Mul):
        return _validate(expr2, expr1, -1)     
    elif isinstance(expr2, Mul):
        return _validate(expr1, expr2, 0)
    else:
        if (isinstance(expr1, State) and expr1.kind == 'bra') and (isinstance(expr2, State) and expr2.kind == 'ket'):
            return InnerProduct(expr1, expr2)
        elif (isinstance(expr1, State) and expr1.kind == 'ket') and (isinstance(expr2, State) and expr2.kind == 'bra'):
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
