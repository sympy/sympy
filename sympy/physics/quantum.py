from sympy import Expr, sympify, Add, Mul, Pow, I, Function, Integer, S, sympify, Matrix
from sympy.core.basic import Atom
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

class State(Expr):
    """
    General abstract quantum state.

    Anywhere you can have a State, you can also have Integer(0).
    All code must check for this!
    """

    lbracket = 'State('
    rbracket = ')'

    hilbert_space = HilbertSpace()

    def __new__(cls, name):
        name = sympify(name)
        obj = Expr.__new__(cls, name, **{'commutative': False})
        return obj

    @property
    def name(self):
        return self.args[0]

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

    def _sympyrepr(self, printer, *args):
        return '%s(%s)' % (self.__class__.__name__, self._print_name(printer, *args))

    def _sympystr(self, printer, *args):
        return '%s%s%s' % (self.lbracket, self._print_name(printer, *args), self.rbracket)

    def _pretty(self, printer, *args):
        from sympy.printing.pretty.stringpict import prettyForm, stringPict
        pform = self._print_name_pretty(printer, *args)
        pform = prettyForm(*pform.left(prettyForm(self.lbracket)))
        pform = prettyForm(*pform.right(prettyForm(self.rbracket)))
        return pform

class Ket(State):

    lbracket = '|'
    rbracket = '>'

    @property
    def dual(self):
        return Bra(*self.args)


class Bra(State):

    lbracket = '<'
    rbracket = '|'

    @property
    def dual(self):
        return Ket(*self.args)


class InnerProduct(Expr):
    """
    An unevaluated inner product between a Bra and Ket.
    """

    def __new__(cls, bra, ket):
        assert isinstance(bra, Bra), 'First argument must be a Bra'
        assert isinstance(ket, Ket), 'Second argument must be a Ket'
        r = cls.eval(bra, ket)
        if isinstance(r, Expr):
            return r
        obj = Expr.__new__(cls, *(bra, ket), **{'commutative': True})
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
        return "%s|%s" % (sbra[:-1], sket[1:])

class OuterProduct(Expr):
    """
    An unevaluated inner product between a Bra and Ket.
    """

    def __new__(cls, ket, bra):
        assert isinstance(ket, Ket), 'First argument must be a Ket'
        assert isinstance(bra, Bra), 'Second argument must be a Bra'
        r = cls.eval(ket, bra)
        if isinstance(r, Expr):
            return r
        obj = Expr.__new__(cls, *(ket, bra), **{'commutative': True})
        return obj

    @classmethod
    def eval(cls, ket, bra):
        # We need to decide what to do here. We probably will ask if the
        # bra and ket know how to do the outer product.
        return None

    @property
    def bra(self):
        return self.args[1]

    @property
    def ket(self):
        return self.args[0]

    def _eval_dagger(self):
        return OuterProduct(Dagger(self.bra), Dagger(self.ket))

    def _sympyrepr(self, printer, *args):
        return '%s(%s,%s)' % (self.__class__.__name__, printer._print(self.ket, *args), printer._print(self.bra, *args))

    def _sympystr(self, printer, *args):
        return str(self.ket)+str(self.bra)

class Operator(Expr):
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

def validate_mul(expr):
    """
    Check to see if the Mul containing quantum objects is valid and return True.

    * Only works for simple Muls (e.g. a*b*c*d) right now.
    """
    expr = split_commutative_parts(expr)[1] #obtain noncommutative parts of mul
    old_arg = None
    for arg in expr:
        if old_arg:
            if isinstance(old_arg, Ket) and isinstance(arg, Ket):
                raise NotImplementedError
            if isinstance(old_arg, Bra) and isinstance(arg, Bra):
                raise NotImplementedError
            if isinstance(old_arg, Operator) and isinstance(arg, Bra):
                raise Exception('Operator*Bra is invalid in quantum mechanics.')
            if isinstance(old_arg, Ket) and isinstance(arg, Operator):
                raise Exception('Ket*Operator is invalid in quantum mechanics.')
        old_arg = arg
    return True

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
        from sympy.printing.pretty.stringpict import prettyForm, stringPict
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

def split_commutative_parts(m):
    c_part = [p for p in m.args if p.is_commutative]
    nc_part = [p for p in m.args if not p.is_commutative]
    return c_part, nc_part

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
