from sympy import Expr, sympify, Add, Mul, Pow, I, Function, Integer, S, sympify
from sympy.core.basic import Atom

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
    'Dagger',
    'KroneckerDelta',
    'Operator',
    'Commutator',
    # 'AntiCommutator',
    'State',
    'Ket',
    'Bra',
    'InnerProduct'
]

class Dagger(Expr):
    """
    General hermitian conjugate operation.
    """

    def __new__(cls, arg):
        arg = sympify(arg)
        r = cls.eval(arg)
        if isinstance(r, Expr):
            return r
        obj = Expr.__new__(cls, arg)
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
                    return Pow(Dagger(arg.args[0]),arg.args[1])
                if arg == I:
                    return -arg
            else:
                return None
        else:
            return d

    def _eval_subs(self, old, new):
        r = Dagger(self.args[0].subs(old, new))
        return r

    def _eval_dagger(self):
        return self.args[0]

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

class Operator(Expr):
    """
    Base class for non-commuting Quantum operators.

    >>> from sympy import simplify
    >>> from sympy.physics.quantum import Operator
    >>> A = Operator('A')
    >>> print A
    A
    """

    def __new__(cls, name):
        name = sympify(name)
        obj = Expr.__new__(cls, name, commutative=False)
        return obj

    @property
    def name(self):
        return self.args[0]

    def doit(self,**kw_args):
        return self

    def _sympyrepr_(self, printer, *args):
        return "%s(%s)" % (self.__class__.__name__, printer._print(self.name, *args))

    def _sympystr_(self, printer, *args):
        return printer._print(self.name, *args)

    def _pretty_(self, printer, *args):
        return printer._print(self.name, *args)

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

class State(Expr):
    """
    General abstract quantum state.

    Anywhere you can have a State, you can also have Integer(0).
    All code must check for this!
    """

    lbracket = 'State('
    rbracket = ')'

    def __new__(cls, name):
        obj = Expr.__new__(cls, name, commutative=False)
        return obj

    @property
    def name(self):
        return self.args[0]

    @property
    def is_symbolic(self):
        return True

    def doit(self,**kw_args):
        return self

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "%s%s%s" % (self.lbracket, self.name, self.rbracket)

class Ket(State):

    lbracket = '|'
    rbracket = '>'

    def _eval_dagger(self):
        return Bra(*self.args)

class Bra(State):

    lbracket = '<'
    rbracket = '|'

    def _eval_dagger(self):
        return Ket(*self.args)

    def __mul__(self, other):
        if isinstance(other, Ket):
            return InnerProduct(self, other)
        else:
            return Expr.__mul__(self, other)

class InnerProduct(Expr):
    """
    An unevaluated inner product between a Bra and Ket.
    """

    def __new__(cls, bra, ket):
        assert isinstance(bra, Bra), 'must be a Bra'
        assert isinstance(ket, Ket), 'must be a Ket'
        r = cls.eval(bra, ket)
        if isinstance(r, Expr):
            return r
        obj = Expr.__new__(cls, *(bra, ket), **dict(commutative=True))
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

    def _eval_subs(self, old, new):
        r = self.__class__(self.bra.subs(old,new), self.ket.subs(old,new))
        return r

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        sbra = str(self.bra)
        sket = str(self.ket)
        return "%s|%s" % (sbra[:-1], sket[1:])

