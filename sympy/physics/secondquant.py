"""
Second quantization operators and states for bosons.

This follow the formulation of Fetter and Welecka, "Quantum Theory
of Many-Particle Systems."
"""

from sympy import (
    Basic, Function, var, Mul, sympify, Integer, Add, sqrt,
    Number, Matrix, zeros, Pow, I, S,Symbol, latex
)

from sympy.utilities.decorator import deprecated

__all__ = [
    'Dagger',
    'KroneckerDelta',
    'BosonicOperator',
    'AnnihilateBoson',
    'CreateBoson',
    'FockState',
    'FockStateBra',
    'FockStateKet',
    'BBra',
    'BKet',
    # 'FBra',
    # 'FKet',
    'F',
    'Fd',
    'B',
    'Bd',
    'apply_operators',
    'InnerProduct',
    'BosonicBasis',
    'VarBosonicBasis',
    'FixedBosonicBasis',
    'commutator',
    'matrix_rep'
]


class Dagger(Basic):
    """
    Hermitian conjugate of creation/annihilation operators.
    """

    def __new__(cls, arg):
        arg = sympify(arg)
        r = cls.eval(arg)
        if isinstance(r, Basic):
            return r
        obj = Basic.__new__(cls, arg)
        return obj

    @classmethod
    @deprecated
    def canonize(cls, arg):
        return cls.eval(arg)

    @classmethod
    def eval(cls, arg):
        try:
            d = arg._dagger_()
        except:
            if isinstance(arg, Basic):
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

    def _dagger_(self):
        return self.args[0]

class KroneckerDelta(Function):
    """
    Discrete delta function.
    """

    nargs = 2
    is_commutative=True

    @classmethod
    def eval(cls, i, j):
        if i > j:
            return cls(j,i)
        diff = i-j
        if diff == 0:
            return Integer(1)
        elif diff.is_number:
            return S.Zero

        if i.assumptions0.get("below_fermi") and j.assumptions0.get("above_fermi"):
            return S.Zero
        if j.assumptions0.get("below_fermi") and i.assumptions0.get("above_fermi"):
            return S.Zero

    def _eval_subs(self, old, new):
        r = KroneckerDelta(self.args[0].subs(old, new), self.args[1].subs(old, new))
        return r

    @property
    def is_above_fermi(self):
        """
        True if Delta can be non-zero above fermi
        """
        if self.args[0].assumptions0.get("below_fermi"):
            return False
        if self.args[1].assumptions0.get("below_fermi"):
            return False
        return True
    @property
    def is_below_fermi(self):
        """
        True if Delta can be non-zero below fermi
        """
        if self.args[0].assumptions0.get("above_fermi"):
            return False
        if self.args[1].assumptions0.get("above_fermi"):
            return False
        return True

    @property
    def indices_contain_equal_information(self):
        if (self.args[0].assumptions0.get("below_fermi") and
                self.args[1].assumptions0.get("below_fermi")):
            return True
        if (self.args[0].assumptions0.get("above_fermi")
                and self.args[1].assumptions0.get("above_fermi")):
            return True

        # if both indices are general we are True, else false
        return self.is_below_fermi and self.is_above_fermi


    @property
    def preffered_index(self):
        if self._get_preffered_index():
            return self.args[1]
        else:
            return self.args[0]

    @property
    def killable_index(self):
        if self._get_preffered_index():
            return self.args[0]
        else:
            return self.args[1]

    def _get_preffered_index(self):
        """
        Returns the index which is preffered to keep in the final expression.

        The preffered index is the index with more information regarding fermi
        level.  If indices contain same information, index 0 is returned.
        """
        if not self.is_above_fermi:
            if self.args[0].assumptions0.get("below_fermi"):
                return 0
            else:
                return 1
        elif not self.is_below_fermi:
            if self.args[0].assumptions0.get("above_fermi"):
                return 0
            else:
                return 1
        else:
            return 0

    def _dagger_(self):
        return self

    def _latex_(self,printer):
        return "\\delta_{%s%s}"% (self.args[0].name,self.args[1].name)

    def __repr__(self):
        return "KroneckerDelta(%s,%s)"% (self.args[0],self.args[1])

    def __str__(self):
        if not self.is_above_fermi:
            return 'd<(%s,%s)'% (self.args[0],self.args[1])
        elif not self.is_below_fermi:
            return 'd>(%s,%s)'% (self.args[0],self.args[1])
        else:
            return 'd(%s,%s)'% (self.args[0],self.args[1])



class SqOperator(Basic):
    """
    Base class for Second Quantization operators.
    """

    op_symbol = 'sq'

    def __new__(cls, k):
        obj = Basic.__new__(cls, sympify(k), commutative=False)
        return obj

    def _eval_subs(self, old, new):
        r = self.__class__(self.args[0].subs(old, new))
        return r

    @property
    def state(self):
        return self.args[0]

    @property
    def is_symbolic(self):
        if self.state.is_Integer:
            return False
        else:
            return True

    def doit(self,**kw_args):
        """
        FIXME: hack to prevent crash further up...
        """
        return self

    def __repr__(self):
        return NotImplemented

    def __str__(self):
        return "%s(%r)" % (self.op_symbol, self.state)

    def apply_operator(self, state):
        raise NotImplementedError('implement apply_operator in a subclass')

class BosonicOperator(SqOperator):
    pass

class Annihilator(SqOperator):
    pass

class Creator(SqOperator):
    pass


class AnnihilateBoson(BosonicOperator, Annihilator):
    """
    Bosonic annihilation operator
    """

    op_symbol = 'b'

    def _dagger_(self):
        return CreateBoson(self.state)

    def apply_operator(self, state):
        if not self.is_symbolic and isinstance(state, FockStateKet):
                element = self.state
                amp = sqrt(state[element])
                return amp*state.down(element)
        else:
            return Mul(self,state)

    def __repr__(self):
        return "AnnihilateBoson(%s)"%self.state

class CreateBoson(BosonicOperator, Creator):
    """
    Bosonic creation operator
    """

    op_symbol = 'b+'

    def _dagger_(self):
        return AnnihilateBoson(self.state)

    def apply_operator(self, state):
        if not self.is_symbolic and isinstance(state, FockStateKet):
                element = self.state
                amp = sqrt(state[element] + 1)
                return amp*state.up(element)
        else:
            return Mul(self,state)

    def __repr__(self):
        return "CreateBoson(%s)"%self.state

B = AnnihilateBoson
Bd = CreateBoson


class FockState(Basic):
    """
    Many particle Fock state with a sequence of occupation numbers.

    Anywhere you can have a FockState, you can also have Integer(0).
    All code must check for this!
    """

    def __new__(cls, occupations):
        """
        occupations is a list with two possible meanings:

        - For bosons it is a list of occupation numbers.
          Element i is the number of particles in state i.

        - For fermions it is a list of occupied orbits.
          Element 0 is the state that was occupied first, element i
          is the i'th occupied state.
        """
        o = map(sympify, occupations)
        obj = Basic.__new__(cls, tuple(o), commutative=False)
        return obj

    def _eval_subs(self, old, new):
        r = self.__class__([o.subs(old, new) for o in self.args[0]])
        return r


    def __getitem__(self, i):
        i = int(i)
        return self.args[0][i]

    def __repr__(self):
        return ("FockState(%r)") % (self.args)

    def __str__(self):
        return "%s%r%s" % (self.lbracket,self._labels(),self.rbracket)

    def _labels(self):
        return self.args[0]

    def __len__(self):
        return len(self.args[0])

class BosonState(FockState):
    """
    Many particle Fock state with a sequence of occupation numbers.

    occupation numbers can be any integer >= 0
    """
    def up(self, i):
        i = int(i)
        new_occs = list(self.args[0])
        new_occs[i] = new_occs[i]+Integer(1)
        return self.__class__(new_occs)

    def down(self, i):
        i = int(i)
        new_occs = list(self.args[0])
        if new_occs[i]==Integer(0):
            return Integer(0)
        else:
            new_occs[i] = new_occs[i]-Integer(1)
            return self.__class__(new_occs)



class FockStateKet(FockState):

    lbracket = '|'
    rbracket = '>'


class FockStateBra(FockState):


    lbracket = '<'
    rbracket = '|'


    def __mul__(self, other):
        if isinstance(other, FockStateKet):
            return InnerProduct(self, other)
        else:
            return Basic.__mul__(self, other)

class FockStateBosonKet(BosonState,FockStateKet):
    def _dagger_(self):
        return FockStateBosonBra(*self.args)

class FockStateBosonBra(BosonState,FockStateBra):
    def _dagger_(self):
        return FockStateBosonKet(*self.args)


BBra = FockStateBosonBra
BKet = FockStateBosonKet

def split_commutative_parts(m):
    c_part = [p for p in m.args if p.is_commutative]
    nc_part = [p for p in m.args if not p.is_commutative]
    return c_part, nc_part


def apply_Mul(m):
    """
    Take a Mul instance with operators and apply them to states.

    This method applies all operators with integer state labels
    to the actual states.  For symbolic state labels, nothing is done.
    When inner products of FockStates are encountered (like <a|b>),
    the are converted to instances of InnerProduct.

    This does not currently work on double inner products like,
    <a|b><c|d>.

    If the argument is not a Mul, it is simply returned as is.
    """
    if not isinstance(m, Mul):
        return m
    c_part, nc_part = split_commutative_parts(m)
    n_nc = len(nc_part)
    if n_nc == 0 or n_nc == 1:
        return m
    else:
        last = nc_part[-1]
        next_to_last = nc_part[-2]
        if isinstance(last, FockStateKet):
            if isinstance(next_to_last, SqOperator):
                if next_to_last.is_symbolic:
                    return m
                else:
                    result = next_to_last.apply_operator(last)
                    if result == 0:
                        return 0
                    else:
                        return apply_Mul(Mul(*(c_part+nc_part[:-2]+[result])))
            elif isinstance(next_to_last, Pow):
                if isinstance(next_to_last.base, SqOperator) and \
                    next_to_last.exp.is_Integer:
                    if next_to_last.base.is_symbolic:
                        return m
                    else:
                        result = last
                        for i in range(next_to_last.exp):
                            result = next_to_last.base.apply_operator(result)
                            if result == 0: break
                        if result == 0:
                            return 0
                        else:
                            return apply_Mul(Mul(*(c_part+nc_part[:-2]+[result])))
                else:
                    return m
            elif isinstance(next_to_last, FockStateBra):
                result = InnerProduct(next_to_last, last)
                if result == 0:
                    return 0
                else:
                    return apply_Mul(Mul(*(c_part+nc_part[:-2]+[result])))
            else:
                return m
        else:
            return m


def apply_operators(e):
    """
    Take a sympy expression with operators and states and apply the operators.
    """
    e = e.expand()
    muls = e.atoms(Mul)
    subs_list = [(m,apply_Mul(m)) for m in iter(muls)]
    return e.subs(subs_list)


class InnerProduct(Basic):
    """
    An unevaluated inner product between a bra and ket.

    Currently this class just reduces things to a prouct of
    KroneckerDeltas.  In the future, we could introduce abstract
    states like |a> and |b>, and leave the inner product unevaluated as
    <a|b>.
    """
    def __new__(cls, bra, ket):
        assert isinstance(bra, FockStateBra), 'must be a bra'
        assert isinstance(ket, FockStateKet), 'must be a key'
        r = cls.eval(bra, ket)
        if isinstance(r, Basic):
            return r
        obj = Basic.__new__(cls, *(bra, ket), **dict(commutative=True))
        return obj

    @classmethod
    @deprecated
    def canonize(cls, bra, ket):
        return cls.eval(bra, ket)

    @classmethod
    def eval(cls, bra, ket):
        result = Integer(1)
        for i,j in zip(bra.args[0], ket.args[0]):
            result *= KroneckerDelta(i,j)
            if result == 0: break
        return result

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
        sbra = repr(self.bra)
        sket = repr(self.ket)
        return "%s|%s" % (sbra[:-1], sket[1:])

    def __str__(self):
        return self.__repr__()


def matrix_rep(op, basis):
    """
    Find the representation of an operator in a basis.
    """
    a = zeros((len(basis), len(basis)))
    for i in range(len(basis)):
        for j in range(len(basis)):
            a[i,j] = apply_operators(Dagger(basis[i])*op*basis[j])
    return a


class BosonicBasis(object):
    """
    Base class for a basis set of bosonic Fock states.
    """
    pass


class VarBosonicBasis(object):
    """
    A single state, variable particle number basis set.
    """

    def __init__(self, n_max):
        self.n_max = n_max
        self._build_states()

    def _build_states(self):
        self.basis = []
        for i in range(self.n_max):
            self.basis.append(FockStateBosonKet([i]))
        self.n_basis = len(self.basis)

    def index(self, state):
        return self.basis.index(state)

    def state(self, i):
        return self.basis[i]

    def __getitem__(self, i):
        return self.state(i)

    def __len__(self):
        return len(self.basis)

    def __repr__(self):
        return repr(self.basis)


class FixedBosonicBasis(BosonicBasis):
    """
    Fixed particle number basis set.
    """
    def __init__(self, n_particles, n_levels):
        self.n_particles = n_particles
        self.n_levels = n_levels
        self._build_particle_locations()
        self._build_states()

    def _build_particle_locations(self):
        tup = ["i"+str(i) for i in range(self.n_particles)]
        first_loop = "for i0 in range(%i)" % self.n_levels
        other_loops = ''
        for i in range(len(tup)-1):
            temp = "for %s in range(%s + 1) " % (tup[i+1],tup[i])
            other_loops = other_loops + temp
        var = "("
        for i in tup[:-1]:
            var = var + i + ","
        var = var + tup[-1] + ")"
        cmd = "result = [%s %s %s]" % (var, first_loop, other_loops)
        exec cmd
        if self.n_particles==1:
            result = [(item,) for item in result]
        self.particle_locations = result

    def _build_states(self):
        self.basis = []
        for tuple_of_indices in self.particle_locations:
            occ_numbers = self.n_levels*[0]
            for level in tuple_of_indices:
                occ_numbers[level] += 1
            self.basis.append(FockStateBosonKet(occ_numbers))
        self.n_basis = len(self.basis)

    def index(self, state):
        return self.basis.index(state)

    def state(self, i):
        return self.basis[i]

    def __getitem__(self, i):
        return self.state(i)

    def __len__(self):
        return len(self.basis)

    def __repr__(self):
        return repr(self.basis)


# def move(e, i, d):
#     """
#     Takes the expression "e" and moves the operator at the position i by "d".
#     """
#     if e.is_Mul:
#         if d == 1:
#             # e = a*b*c*d
#             a = Mul(*e.args[:i])
#             b = e.args[i]
#             c = e.args[i+1]
#             d = Mul(*e.args[i+2:])
#             if isinstance(b, Dagger) and not isinstance(c, Dagger):
#                 i, j = b.args[0].args[0], c.args[0]
#                 return a*c*b*d-a*KroneckerDelta(i, j)*d
#             elif not isinstance(b, Dagger) and isinstance(c, Dagger):
#                 i, j = b.args[0], c.args[0].args[0]
#                 return a*c*b*d-a*KroneckerDelta(i, j)*d
#             else:
#                 return a*c*b*d
#         elif d == -1:
#             # e = a*b*c*d
#             a = Mul(*e.args[:i-1])
#             b = e.args[i-1]
#             c = e.args[i]
#             d = Mul(*e.args[i+1:])
#             if isinstance(b, Dagger) and not isinstance(c, Dagger):
#                 i, j = b.args[0].args[0], c.args[0]
#                 return a*c*b*d-a*KroneckerDelta(i, j)*d
#             elif not isinstance(b, Dagger) and isinstance(c, Dagger):
#                 i, j = b.args[0], c.args[0].args[0]
#                 return a*c*b*d-a*KroneckerDelta(i, j)*d
#             else:
#                 return a*c*b*d
#         else:
#             if d > 1:
#                 while d >= 1:
#                     e = move(e, i, 1)
#                     d -= 1
#                     i += 1
#                 return e
#             elif d < -1:
#                 while d <= -1:
#                     e = move(e, i, -1)
#                     d += 1
#                     i -= 1
#                 return e
#     elif isinstance(e, Add):
#         a, b = e.as_two_terms()
#         return move(a, i, d) + move(b, i, d)
#     raise NotImplementedError()

def commutator(a, b):
    """
    Return the commutator:  [a, b] = a*b - b*a
    """
    return a*b - b*a



