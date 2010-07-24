from sympy import (
    Expr, Basic, sympify, Add, Mul, Pow, 
    I, Function, Integer, S, sympify, Matrix
)
from sympy.physics.hilbert import *

"""
Notes:

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
* How are the names of an Operator different than the class?
* Should operators be singletons?
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

    def __mul__(self, other):
        qmul = Mul(self, other)
        if validate_mul(qmul):
            return qmul

    def __rmul__(self, other):
        qrmul = Mul(other, self)
        if validate_mul(qrmul):
            return qrmul

    def _print_name(self, printer, *args):
        return printer._print(self.args[0], *args)

    def _print_name_pretty(self, printer, *args):
        pform = printer._print(self.args[0], *args)
        return pform

    def doit(self,**kw_args):
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

    hilbert_space = HilbertSpace()

    def __new__(cls, name):
        name = sympify(name)
        obj = Expr.__new__(cls, name, **{'commutative': False})
        return obj

    @property
    def name(self):
        return self.args[0]

    def __mul__(self, other):
        qmul = Mul(self, other)
        if validate_mul(qmul):
            return qmul

    def __rmul__(self, other):
        qrmul = Mul(other, self)
        if validate_mul(qrmul):
            return qrmul

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
        r = cls.eval(*args)
        if isinstance(r, Expr):
            return r
        obj = Expr.__new__(cls, bra, ket)
        return obj

    @property
    def bra(self):
        return self.args[0]

    @property
    def ket(self):
        return self.args[1]

    def _eval_dagger(self):
        return InnerProduct(Dagger(self.ket), Dagger(self.bra))

    def _innerproduct_printer(self, printer, *args):
        print_elements = []
        for arg in self.args:
            print_elements.append(printer._print(arg, *args))
        return print_elements

    def _sympyrepr(self, printer, *args):
        innerproduct_reprs = self._innerproduct_printer(printer, *args)
        return '%s(%s)' % (self.__class__.__name__, ', '.join(innerproduct_reprs))

    def _sympystr(self, printer, *args):
        if len(self.args) == 2:
            sbra = str(self.bra)
            sket = str(self.ket)
            return '%s|%s' % (sbra[:-1], sket[1:])
        else:
            innerproduct_strs = self._innerproduct_printer(printer, *args)
            return '%s' % ''.join(innerproduct_strs)

    def _pretty(self, printer, *args):
#        if len(self.args) == 2:
#            return '<%s|%s>' % (self.bra._pretty(printer, *args), self.ket._pretty(printer, *args))
#        else:
            pretty_args = None
            for arg in self.args:
                if pretty_args != None:
                    pretty_args = pretty_args*arg._pretty(printer, *args)
                else:
                    pretty_args = arg._pretty(printer, *args)
            return pretty_args


class OuterProduct(Expr):
    """
    An unevaluated outer product between a Ket and Bra.
    """

    def __new__(cls, ket, bra):
        if not (ket and bra):
            raise Exception('OuterProduct requires a leading Ket and a trailing Bra')
        assert isinstance(ket, Ket), 'First argument must be a Ket'
        assert isinstance(bra, Bra), 'Second argument must be a Bra'
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
            result += represent(basis, **options)
        return result
    elif isinstance(expr, Pow):
        return represent(expr.base, basis, **options)**expr.exp
    else:
        return expr

    if not isinstance(expr, Mul):
        raise TypeError('Mul expected, got: %r' % expr)

    for arg in reversed(expr.args):
        result = S.One
        result *= represent(arg, basis, **options)
    return result


def validate_mul(expr):
    """
    Check to see if the Mul containing quantum objects is valid and return True.

    * Only works for simple Muls (e.g. a*b*c*d) right now.
    """
    expr = split_commutative_parts(expr)[1] #obtain noncommutative parts of mul
    old_arg = None
    for arg in expr:
        if old_arg != None:
            if isinstance(old_arg, Ket) and isinstance(arg, Ket):
                raise NotImplementedError
            if isinstance(old_arg, Bra) and isinstance(arg, Bra):
                raise NotImplementedError
            if (isinstance(old_arg, Operator) or isinstance(old_arg, OuterProduct)) and isinstance(arg, Bra):
                raise Exception('(Operator or OuterProduct)*Bra is invalid in quantum mechanics.')
            if isinstance(old_arg, Ket) and (isinstance(arg, Operator) or isinstance(arg, OuterProduct)):
                raise Exception('Ket*(Operator or OuterProduct) is invalid in quantum mechanics.')
        old_arg = arg
    return True


def combine_innerproduct(expr):
    """
    Combines a (valid) Mul of quantum objects into inner products (if possible).

    * Only works for simple Muls (e.g. a*b*c*d) right now.
    """
    if validate_mul(expr):
        new_expr = Mul()
        inner_product = []
        left_of_bra = False
        for arg in expr.args:
            if isinstance(arg, Bra):
                left_of_bra = True
                inner_product.append(arg)
            elif (isinstance(arg, Ket) and left_of_bra == True):
                inner_product.append(arg)
                new_expr = new_expr*InnerProduct(*inner_product)
                inner_product = []
                left_of_bra = False
            elif left_of_bra == True:
                inner_product.append(arg)
            else:
                new_expr = new_expr*arg
        return new_expr*Mul(*inner_product)


def combine_outerproduct(expr):
    """
    Combines a (valid) Mul of quantum objects into outer products (if possible).

    * Only works for simple Muls (e.g. a*b*c*d) right now.
    """
    if validate_mul(expr):
        old_arg = None
        new_expr = Mul()
        left_of_bra = False
        for arg in expr.args:
            if old_arg != None:
                if (isinstance(old_arg, Ket) and isinstance(arg, Bra)):
                    new_expr = new_expr*OuterProduct(old_arg, arg)
                    left_of_bra = True
                elif left_of_bra == True:
                    left_of_bra = False
                else:
                    new_expr = new_expr*old_arg
            old_arg = arg
        if left_of_bra == True:
            return new_expr
        else:
            return new_expr*old_arg


def split_product(expr):
    """
    Separates a (valid) Mul of quantum objects and inner/outer products into a 
    Mul.

    * Only works for simple Muls (e.g. a*b*c*d) right now.
    """
    if validate_mul(expr):
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


