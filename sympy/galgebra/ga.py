# sympy/galgebra/ga.py

"""
ga.py implements the symbolic geometric algebra of an n-dimensional
vector space with constant metric (future versions will allow for a
metric that is a function of coordinates) with an arbitrary set of
basis vectors (whether they are orthogonal or not depends on the
metric the use inputs).

All the products of geometric algebra -

    geometric
    outer (wedge)
    inner (dots)
    left contraction
    right contaction

and the geometric derivative applied to all products from both sides.

For more information on the details of geometric algebra please look
at the documentation for the module.
"""

from __future__ import print_function

from sympy.core.compatibility import reduce
from itertools import combinations

import copy
import operator

from sympy import Symbol, Expr, expand, Mul, Add, S, collect, \
    Function, simplify, diff, trigsimp, sqrt, Number, \
    factor_terms, sin, cos, sinh, cosh
from sympy import N as Nsympy
from sympy.core.compatibility import range

from sympy.galgebra.printing import GA_Printer, GA_LatexPrinter, enhance_print, latex
from sympy.galgebra.vector import Vector
from sympy.galgebra.debug import oprint
from sympy.galgebra.stringarrays import fct_sym_array, str_combinations
from sympy.galgebra.ncutil import linear_expand, bilinear_product, nc_substitue, \
    get_commutative_coef, ONE_NC


def diagpq(p, q=0):
    """
    Returns string equivalent metric tensor for signature (p, q).
    """
    n = p + q
    D = []
    for i in range(p):
        D.append((i*'0 ' +'1 '+ (n-i-1)*'0 ')[:-1])
    for i in range(p,n):
        D.append((i*'0 ' +'-1 '+ (n-i-1)*'0 ')[:-1])
    return ','.join(D)


def arbitrary_metric(n):
    """
    Returns string equivalent metric tensor for arbitrary signature.
    """
    return ','.join(n*[(n*'# ')[:-1]])


def arbitrary_metric_conformal(n):
    """
    Returns string equivalent metric tensor for arbitrary signature (n+1,1).
    """
    str1 = ','.join(n*[n*'# '+'0 0'])
    return ','.join([str1, n*'0 '+'1 0', n*'0 '+'0 -1'])


def make_coef(self, coef_str):
    if self.fct:
        if self.vars is not None:
            return Function(coef_str)(*self.vars)
        elif MV.coords is not None:
                return Function(coef_str)(*MV.coords)
        else:
            return Symbol(coef_str)
    else:
        return Symbol(coef_str)


class MV(object):
    """
    'MV' class wraps sympy expressions of the form

        s = s_0 + s_1*b_1 + ... + s_K*b_K

    where the s_i are real sympy scalars (commutative expressions) and
    the b_i are non-commutative sympy symbols.  For an N-dimensional
    vector space K = 2**N - 1.

    The linear combination of scalar (commutative) sympy quatities and the
    basis multivectors form the multivector space.  If the number of basis
    vectors is 'n' the dimension of the multivector space is 2**n. If the
    basis of the underlying vector space is (a_1,...,a_n) then the bases
    of the multivector space are the noncommunicative geometric products
    of the basis vectors of the form a_i1*a_i2*...*a_ir where i1 < i2 < ... < ir
    (normal order) and the scalar 1.  A multivector space is the vector
    space with these bases over the sympy scalars.  A basic assumption of
    the geometric product, '*', is that it is associative and that the
    geometric product of a vector with itself is a scalar.  Thus we define
    for any two vectors -

        a.b = (a*b + b*a)/2 [1] (D&L 4.7)

    noting then that a.a = a*a, a.b = b.a, and that a.b is a scalar. The
    order of the geometric product of any two vectors can be reversed
    with -

        b*a = 2*(a.b) - a*b [2] (D&L 4.30)

    This is all that is required to reduce the geometric product of any
    number of basis vectors in any order to a linear combination of
    normal order basis vectors.  Note that a dot product for these bases
    has not yet been defined and when it is the bases will not be orthogonal
    unless the basis vectors are orthogonal.

    The outer product of two vectors is defined to be -

        a^b = (a*b - b*a)/2 [3] (D&L 4.8)

    This is generalized by the formula

        a^R_k = (a*R_k + (-1)**k*R_k*a)/2 [4] (D&L 4.38)

    where R_k is the outer product of k vectors (k-blade) and equation
    [4] recursively defines the outer product of k + 1 vectors in terms of
    the linear combination of geometric products of terms with k + 1 and
    fewer vectors.

    D&L is "Geometric Algebra for Physicists" by Chris Doran and
    Anthony Lasenby, Cambridge University Press.
    """

    ##########Methods for products (*,^,|) of orthogonal blades#########
    """
    No multiplication tables (*,^,|) are calculated if the basis vectors
    are orthogonal.  All types of products are calculated on the fly and
    the basis bases and blades are identical.
    """

    latex_flg = False
    dot_mode = 's'  # 's' - symmetric, 'l' - left contraction, 'r' - right contraction

    @staticmethod
    def product_orthogonal_blades(blade1, blade2):
        blade_index = list(MV.blade_to_index[blade1] + MV.blade_to_index[blade2])
        repeats = []
        sgn = 1
        for i in range(1, len(blade_index)):
            save = blade_index[i]
            j = i
            while j > 0 and blade_index[j - 1] > save:
                sgn = -sgn
                blade_index[j] = blade_index[j - 1]
                j -= 1
            blade_index[j] = save
            if blade_index[j] == blade_index[j - 1]:
                repeats.append(save)
        result = S(sgn)
        for i in repeats:
            blade_index.remove(i)
            blade_index.remove(i)
            result *= MV.metric[i]
        result *= MV.index_to_blade[tuple(blade_index)]
        return result

    @staticmethod
    def dot_orthogonal_blades(blade1, blade2):

        index1 = MV.blade_to_index[blade1]
        index2 = MV.blade_to_index[blade2]
        index = list(index1 + index2)
        grade1 = len(index1)
        grade2 = len(index2)

        if MV.dot_mode == 's':
            if grade1 == 0:
                return S.Zero
            elif grade2 == 0:
                return S.Zero
            else:
                grade = abs(grade1 - grade2)
        elif MV.dot_mode == 'l':
            grade = grade2 - grade1
            if grade < 0:
                return S.Zero
            if grade1 == 0:
                return blade2
        elif MV.dot_mode == 'r':
            grade = grade1 - grade2
            if grade < 0:
                return S.Zero
            if grade2 == 0:
                return blade1
        n = len(index)
        sgn = 1
        result = S.One
        ordered = False
        while n > grade:
            ordered = True
            i2 = 1
            while i2 < n:
                i1 = i2 - 1
                index1 = index[i1]
                index2 = index[i2]
                if index1 == index2:
                    n -= 2
                    if n < grade:
                        return S.Zero
                    result *= MV.metric[index1]
                    index = index[:i1] + index[i2 + 1:]
                elif index1 > index2:
                    ordered = False
                    index[i1] = index2
                    index[i2] = index1
                    sgn = -sgn
                    i2 += 1
                else:
                    i2 += 1
            if ordered:
                break
        if n > grade:
            return S.Zero
        else:
            return sgn * result * MV. index_to_blade[tuple(index)]

    ######################Multivector Constructors######################

    def __init__(self, base=None, mvtype=None, fct=False, blade_rep=False):

        """
        Initialization of multivector X. Inputs are as follows

        mvtype           base                       result

        default         default                  Zero multivector
        'basisvector'   int i                    ith basis vector
        'basisbivector' int i                    ith basis bivector
        'scalar'         x                       scalar of value x
                        's'
        'grade'         [A]                      X.grade(i) = A
                        's,i'
        'vector'        [A]                      X.grade(1) = [A]
                        's'
 'grade2' or 'bivector' [A]                      X.grade(2) = A
                        's'
        'pseudo'         x                       X.grade(n) = x
                        's'
        'spinor'        's'                      spinor with coefficients
                                                 s__indices and name s
        'mv'            's'                      general multivector with
                                                 s__indices and name s

        If fct is 'True' and MV.coords is defined in MV.setup then a
        multivector field of MV.coords is instantiated.

        Multivector data members are:

            obj       - a sympy expression consisting of a linear
                        combination of sympy scalars and bases/blades.

            blade_rep - 'True' if 'MV' representation is a blade expansion,
                        'False' if 'MV' representation is a base expansion.
        """

        def make_scalar(self, base):  # make a scalar (grade 0)
            if isinstance(base, str):
                if self.fct:
                    self.obj = Function(base)(*MV.coords) * MV.ONE
                else:
                    self.obj = make_coef(self, base) * MV.ONE
            else:
                self.obj = base * MV.ONE
            self.igrade = 0
            self.blade_rep = True
            return self

        def make_vector(self, base):  # make a vector (grade 1)
            if isinstance(base, str):
                if self.fct:
                    base_lst = str_combinations(base, MV.coords, rank=1, mode='__')
                    fct_lst = fct_sym_array(base_lst, MV.coords)
                    self.obj = reduce(operator.add, tuple(x0 * x1 for x0, x1 in zip(fct_lst, MV.blades[1])))
                else:
                    if MV.coords is not None:
                        base_lst = str_combinations(base, MV.coords, rank=1, mode='__')
                    else:
                        base_lst = str_combinations(base, MV.subscripts, rank=1, mode='__')
                    fct_lst = fct_sym_array(base_lst, None)
                    self.obj = reduce(operator.add, tuple(x0 * x1 for x0, x1 in zip(fct_lst, MV.blades[1])))
            else:
                result = S.Zero
                for (coef, base) in zip(base, MV.blades[1]):
                    result += coef * base
                self.obj = result
            self.igrade = 1
            self.blade_rep = True
            return self

        def make_basisvector(self, base):
            raise NotImplementedError("Don't know how to compute basis vectors of class %" % self.__class__)

        def make_basisbivector(self, base):
            raise NotImplementedError("Don't know how to compute basis bivectors of class %" % self.__class__)

        def make_grade(self, base):  # if base is 'A,n' then make a grade n multivector
            if isinstance(base, str):
                base_lst = base.split(',')
                base = base_lst[0]
                n = int(base_lst[1])
                if self.fct:
                    base_lst = str_combinations(base, MV.coords, rank=n, mode='__')
                    fct_lst = fct_sym_array(base_lst, MV.coords)
                    self.obj = reduce(operator.add, tuple(x0 * x1 for x0, x1 in zip(fct_lst, MV.blades[n])))
                else:
                    if MV.coords is not None:
                        base_lst = str_combinations(base, MV.coords, rank=n, mode='__')
                    else:
                        base_lst = str_combinations(base, MV.subscripts, rank=n, mode='__')
                    fct_lst = fct_sym_array(base_lst, None)
                    self.obj = reduce(operator.add, tuple(x0 * x1 for x0, x1 in zip(fct_lst, MV.blades[n])))
            else:
                raise TypeError('Cannot make_grade for base = %s' % base)
            self.igrade = n
            self.blade_rep = True
            return self

        def make_grade2(self, base):  # grade 2 multivector
            if isinstance(base, str):
                if self.fct:
                    base_lst = str_combinations(base, MV.coords, rank=2, mode='__')
                    fct_lst = fct_sym_array(base_lst, MV.coords)
                    self.obj = reduce(operator.add, tuple(x0 * x1 for x0, x1 in zip(fct_lst, MV.blades[2])))
                else:
                    if MV.coords is not None:
                        base_lst = str_combinations(base, MV.coords, rank=2, mode='__')
                    else:
                        base_lst = str_combinations(base, MV.subscripts, rank=2, mode='__')
                    fct_lst = fct_sym_array(base_lst, None)
                    self.obj = reduce(operator.add, tuple(x0 * x1 for x0, x1 in zip(fct_lst, MV.blades[2])))
            else:
                raise TypeError('!!!!Cannot make_grade2 for base = ' + str(base) + '!!!!\n')
            self.igrade = 2
            self.blade_rep = True
            return self

        def make_pseudo(self, base):  # multivector of grade MV.dim
            if isinstance(base, str):
                if self.fct:
                    base_lst = str_combinations(base, MV.coords, rank=MV.dim, mode='__')
                    fct_lst = fct_sym_array(base_lst, MV.coords)
                    self.obj = reduce(operator.add, tuple(x0 * x1 for x0, x1 in zip(fct_lst, MV.blades[MV.dim])))
                else:
                    if MV.coords is not None:
                        base_lst = str_combinations(base, MV.coords, rank=MV.dim, mode='__')
                    else:
                        base_lst = str_combinations(base, MV.subscripts, rank=MV.dim, mode='__')
                    fct_lst = fct_sym_array(base_lst, None)
                    self.obj = reduce(operator.add, tuple(x0 * x1 for x0, x1 in zip(fct_lst, MV.blades[MV.dim])))
            else:
                raise TypeError('!!!!Cannot make_pseudo for base = ' + str(base) + '!!!!\n')
            self.igrade = MV.dim
            self.blade_rep = True
            return self

        def make_spinor(self, base):  # multivector with all even grades
            if isinstance(base, str):
                if self.fct:
                    self.obj = Function(base)(*MV.coords) * MV.ONE
                else:
                    self.obj = Symbol(base) * MV.ONE
                for rank in range(2, MV.dim1, 2):
                    if self.fct:
                        base_lst = str_combinations(base, MV.coords, rank=rank, mode='__')
                        fct_lst = fct_sym_array(base_lst, MV.coords)
                        self.obj += reduce(operator.add, tuple(x0 * x1 for x0, x1 in zip(fct_lst, MV.blades[rank])))
                    else:
                        if MV.coords is not None:
                            base_lst = str_combinations(base, MV.coords, rank=rank, mode='__')
                        else:
                            base_lst = str_combinations(base, MV.subscripts, rank=rank, mode='__')
                        fct_lst = fct_sym_array(base_lst, None)
                        self.obj += reduce(operator.add, tuple(x0 * x1 for x0, x1 in zip(fct_lst, MV.blades[rank])))
            else:
                raise TypeError('Cannot make_mv for base = %s' % base)
            self.igrade = -1
            self.blade_rep = True
            return self

        def make_mv(self, base):
            if isinstance(base, str):
                if self.fct:
                    self.obj = Function(base)(*MV.coords) * MV.ONE
                else:
                    self.obj = Symbol(base) * MV.ONE
                for rank in range(1, MV.dim1):
                    if self.fct:
                        base_lst = str_combinations(base, MV.coords, rank=rank, mode='__')
                        fct_lst = fct_sym_array(base_lst, MV.coords)
                        self.obj += reduce(operator.add, tuple(x0 * x1 for x0, x1 in zip(fct_lst, MV.blades[rank])))
                    else:
                        if MV.coords is not None:
                            base_lst = str_combinations(base, MV.coords, rank=rank, mode='__')
                        else:
                            base_lst = str_combinations(base, MV.subscripts, rank=rank, mode='__')
                        fct_lst = fct_sym_array(base_lst, None)
                        self.obj += reduce(operator.add, tuple(x0 * x1 for x0, x1 in zip(fct_lst, MV.blades[rank])))
            else:
                raise TypeError('!!!!Cannot make_mv for base = ' + str(base) + '!!!!\n')
            self.igrade = -1
            self.blade_rep = True
            return self

        MVtypes = {'scalar': make_scalar,
                   'vector': make_vector,
                   'basisvector': make_basisvector,
                   'basisbivector': make_basisbivector,
                   'grade': make_grade,
                   'grade2': make_grade2,
                   'bivector': make_grade2,
                   'pseudo': make_pseudo,
                   'spinor': make_spinor,
                   'mv': make_mv}

        self.fct = fct
        self.is_base = False
        self.is_grad = False
        self.print_blades = MV.print_blades
        self.fmt = 1

        if mvtype is None:
            if base in (None, S.Zero):  # Default is zero multivector
                self.blade_rep = True
                self.obj = S.Zero
                self.igrade = 0
            elif isinstance(base, str):  # Base or blade basis multivector
                self.is_base = True
                if '*' in base:
                    self.blade_rep = False
                    self.igrade = -1
                else:
                    if '^' in base:
                        self.blade_rep = True
                        self.igrade = base.count('^') + 1
                    else:
                        self.blade_rep = blade_rep
                        self.igrade = 1
                self.obj = Symbol(base, commutative=False)
            elif isinstance(base, MV):  # Copy constructor
                self.blade_rep = base.blade_rep
                self.obj = base.obj
                self.igrade = base.igrade
                self.fct = base.fct
                self.is_base = base.is_base
                self.is_grad = base.is_grad
            elif isinstance(base, (Expr, Symbol)):  # Gets properties of multivector from Expr
                if base.is_commutative:
                    self.obj = base * MV.ONE
                    self.blade_rep = True
                    self.igrade = 0
                else:
                    if isinstance(base, (Add, Mul)):  # Complex expression
                        MV.characterize_expression(self, base)
                    elif isinstance(base, Symbol):
                        if not base.is_commutative:
                            if base == MV.ONE:
                                self.obj = base
                                self.blade_rep = True
                                self.igrade = 0
                            elif base in MV.blades_flat:  # basis blade
                                self.obj = base
                                self.blade_rep = True
                                self.igrade = MV.blade_grades[base]
                            elif base in MV.bases_flat:  # basis base
                                self.obj = base
                                self.blade_rep = False
                                self.igrade = -1
                            else:
                                raise ValueError('MV(' + str(base) + ') is not allowed in constructor\n' +
                                                 'non-commutative argument is not a base\n')
                        else:  # scalar sympy symbol
                            self.obj = base * MV.ONE
                            self.igrade = 0
                            self.blade_rep = True
                    elif isinstance(base, Number):
                        self.obj = base * MV.ONE
                        self.igrade = 0
                        self.blade_rep = True
        else:  # Preconfigured multivector types
            MVtypes[mvtype](self, base)

    def Fmt(self, fmt=1, title=None):
        self.fmt = fmt
        if title is not None:
            print(title + ' = ' + str(self))
            return
        return self

    def __str__(self):
        if GA_LatexPrinter.latex_flg:
            Printer = GA_LatexPrinter
        else:
            Printer = GA_Printer
        self.discover_and_set_grade()
        self.obj = expand(self.obj).collect(MV.blades_flat)
        return Printer().doprint(self)

    ######################## Operator Definitions#######################

    def coef(self, base):
        (coefs, bases) = linear_expand(self.obj)
        if base.obj in bases:
            icoef = bases.index(base.obj)
            return coefs[icoef]
        else:
            return S.Zero

    def func(self, fct):
        (coefs, bases) = linear_expand(self.obj)
        result = S.Zero
        for (coef, base) in zip(coefs, bases):
            result += fct(coef) * base
        fself = MV(self)
        fself.obj = result
        return fself

    def __eq__(self, mv):
        if not isinstance(mv, MV):
            mv = MV(mv)
        if self.obj == mv.obj:
            return True
        else:
            return False

    def __neg__(self):  # -self
        nself = MV(self)
        nself.obj = -self.obj
        return nself

    def __pos__(self):  # +self
        return self

    def __add__(self, b):  # self + b
        if isinstance(b, MV):
            self.base_to_blade()
            b.base_to_blade()
            self_add_b = MV.basic_add(self, b)
            self_add_b.discover_and_set_grade()
        else:
            self_add_b = MV(self)
            self_add_b.obj = self.obj + b * MV.ONE
            if self.igrade != 0:
                self_add_b.igrade = -1
        return self_add_b

    def __radd__(self, b):  # b + self
        b_add_self = MV(self)
        b_add_self.obj = b * MV.ONE + self.obj
        if self.igrade != 0:
            b_add_self.igrade = -1
        return b_add_self

    def __add_ab__(self, b):  # self += b
        selfb = MV()
        selfb.obj += b.obj
        return selfb

    def __sub__(self, b):  # self - b
        if isinstance(b, MV):
            self.base_to_blade()
            b.base_to_blade()
            self_sub_b = MV.basic_sub(self, b)
            self_sub_b.discover_and_set_grade()
        else:
            self_sub_b = MV(self)
            self_sub_b.obj = self.obj - b * MV.ONE
            if self.igrade != 0:
                self_sub_b.igrade = -1
        return self_sub_b

    def __sub_ab__(self, b):  # self -= b
        selfb = MV()
        selfb.obj -= b.obj
        return selfb

    def __rsub__(self, b):  # b - self
        b_sub_self = MV(self)
        b_sub_self.obj = b * MV.ONE - self.obj
        if self.igrade != 0:
            b_sub_self.igrade = -1
        return b_sub_self

    def __mul__(self, b):  # self*b
        if isinstance(b, MV):
            if self.is_grad:  # left geometric derivative
                result = MV()
                if self.connection:
                    for (coord, brecp, bnorm) in zip(self.coords, self.rcpr_bases_MV, self.tangent_norm):
                        result += (brecp * b.diff(coord)) / bnorm
                    if b.igrade > 0:
                        result.obj += nc_substitue(b.obj, self.connection['left'])
                else:
                    for (coord, brecp) in zip(self.coords, self.rcpr_bases_MV):
                        result += brecp * b.diff(coord)
                return result
            elif b.is_grad:  # right geometric derivative
                result = MV()
                if b.connection:
                    for (coord, brecp, bnorm) in zip(b.coords, b.rcpr_bases_MV, b.tangent_norm):
                        result += (self.diff(coord) * brecp) / bnorm
                    if self.igrade > 0:
                        result.obj += nc_substitue(self.obj, b.connection['right'])
                else:
                    for (coord, brecp) in zip(b.coords, b.rcpr_bases_MV):
                        result += self.diff(coord) * brecp
                return MV(result)
            else:
                if not MV.is_orthogonal:
                    self.blade_to_base()
                    b.blade_to_base()
                obj = bilinear_product(self.obj * b.obj, MV.geometric_product)
                self_mul_b = MV(obj)
                if not MV.is_orthogonal:
                    self_mul_b.igrade = -1
                    self_mul_b.blade_rep = False
                self_mul_b.base_to_blade()
                self_mul_b.discover_and_set_grade()
                return self_mul_b
        else:
            if self.is_grad:
                result = MV()
                for (coord, brecp) in zip(self.coords, self.rcpr_bases_MV):
                    result += brecp * diff(b, coord)
                return result
            else:
                self_mul_b = MV(self)
                self_mul_b.obj *= b
        return self_mul_b

    def __mul_ab__(self, b):  # self *= b
        selfb = MV(self)
        selfb.obj *= b.obj
        return selfb

    def __rmul__(self, b):  # b * self
        b_mul_self = MV(self)
        b_mul_self.obj = b * self.obj
        return b_mul_self

    def __div__(self, b):  # self / b
        if not isinstance(b, MV):
            self_div_b = MV(self)
            self_div_b.obj = self.obj / b
            return self_div_b
        else:
            raise TypeError('No multivector division for divisor = ' + str(b) + '\n')

    __truediv__ = __div__

    def __or__(self, b):  # self | b
        if isinstance(b, MV):
            if self.is_grad:  # left dot/div (inner) derivative
                result = MV()
                if b.igrade == 0:
                    return result
                if self.connection:
                    for (coord, brecp, bnorm) in zip(self.coords, self.rcpr_bases_MV, self.tangent_norm):
                            result += (brecp | b.diff(coord)) / bnorm
                    result.obj += nc_substitue(b.obj, self.connection['left_dot'])
                else:
                    for (coord, brecp) in zip(self.coords, self.rcpr_bases_MV):
                            result += brecp | b.diff(coord)
                return result
            elif b.is_grad:  # right dot/div (inner) derivative
                result = MV()
                if b.connection:
                    for (coord, brecp, bnorm) in zip(b.coords, b.rcpr_bases_MV, b.tangent_norm):
                            result += (self.diff(coord) | brecp) / bnorm
                    result.obj += nc_substitue(self.obj, b.connection['right_dot'])
                else:
                    for (coord, brecp) in zip(b.coords, b.rcpr_bases_MV):
                            result += self.diff(coord) | brecp
                return MV(result)
            else:
                if MV.is_orthogonal:
                    MV.dot_mode = 's'
                    result = bilinear_product(self.obj * b.obj, MV.dot_orthogonal_blades)
                    return MV(result)
                else:
                    return MV.non_orthogonal_products(self, b, mode='s')
        else:  # dot product returns zero for r.h.s. scalar multiplier
            return MV()

    def __ror__(self, b):  # b | self
        b_dot_self = MV()
        return b_dot_self

    def __lt__(self, b):  # left contraction
        if isinstance(b, MV):
            if self.is_grad:  # left derivative for left contraction
                result = MV()
                if self.connection:
                    for (coord, brecp, bnorm) in zip(self.coords, self.rcpr_bases_MV, self.tangent_norm):
                        result += (brecp < b.diff(coord)) / bnorm
                    result.obj += nc_substitue(b.obj, self.connection['left_dot'])
                else:
                    for (coord, brecp) in zip(self.coords, self.rcpr_bases_MV):
                        result += brecp < b.diff(coord)
                return result
            elif b.is_grad:  # right derivative for left contraction
                result = MV()
                if b.connection:
                    for (coord, brecp, bnorm) in zip(b.coords, b.rcpr_bases_MV, b.tangent_norm):
                        result += (self.diff(coord) < brecp) / bnorm
                    result.obj += nc_substitue(self.obj, b.connection['right_dot'])
                else:
                    for (coord, brecp) in zip(b.coords, b.rcpr_bases_MV):
                        result += self.diff(coord) < brecp
                return MV(result)
            else:
                if MV.is_orthogonal:
                    MV.dot_mode = 'l'
                    result = bilinear_product(self.obj * b.obj, MV.dot_orthogonal_blades)
                    return MV(result)
                else:
                    return MV.non_orthogonal_products(self, b, mode='l')

    def __gt__(self, b):  # right contraction
        if isinstance(b, MV):
            if self.is_grad:  # left derivative for right contraction
                result = MV()
                if self.connection:
                    for (coord, brecp, bnorm) in zip(self.coords, self.rcpr_bases_MV, self.tangent_norm):
                        result += (brecp > b.diff(coord)) / bnorm
                    result.obj += nc_substitue(b.obj, self.connection['left_dot'])
                else:
                    for (coord, brecp) in zip(self.coords, self.rcpr_bases_MV):
                        result += brecp > b.diff(coord)
                return result
            elif b.is_grad:  # right derivative for right contraction
                result = MV()
                if b.connection:
                    for (coord, brecp, bnorm) in zip(b.coords, b.rcpr_bases_MV, b.tangent_norm):
                        result += (self.diff(coord) > brecp) / bnorm
                    result.obj += nc_substitue(self.obj, b.connection['right_dot'])
                else:
                    for (coord, brecp) in zip(b.coords, b.rcpr_bases_MV):
                        result += self.diff(coord) > brecp
                return MV(result)
            else:
                if MV.is_orthogonal:
                    MV.dot_mode = 'r'
                    result = bilinear_product(self.obj * b.obj, MV.dot_orthogonal_blades)
                    return MV(result)
                else:
                    return MV.non_orthogonal_products(self, b, mode='r')

    def __xor__(self, b):  # self ^ b
        if isinstance(b, MV):
            if self.is_grad:  # left wedge/curl (outer) derivative
                result = MV()
                if self.connection:
                    for (coord, brecp, bnorm) in zip(self.coords, self.rcpr_bases_MV, self.tangent_norm):
                        result += (brecp ^ b.diff(coord)) / bnorm
                    result.obj += nc_substitue(b.obj, self.connection['left_wedge'])
                else:
                    for (coord, brecp) in zip(self.coords, self.rcpr_bases_MV):
                        result += brecp ^ b.diff(coord)
                return result
            elif b.is_grad:  # right wedge/curl (outer) derivative
                result = MV()
                if b.connection:
                    for (coord, brecp, bnorm) in zip(b.coords, b.rcpr_bases_MV, b.tangent_norm):
                        result += (self.diff(coord) ^ brecp) / bnorm
                    result.obj += nc_substitue(self.obj, b.connection['right_wedge'])
                else:
                    for (coord, brecp) in zip(b.coords, b.rcpr_bases_MV):
                        result += self.diff(coord) ^ brecp
                return MV(result)
            else:
                if MV.is_orthogonal:
                    result = bilinear_product(self.obj * b.obj, MV.wedge_product)
                    return MV(result)
                else:
                    return MV.non_orthogonal_products(self, b, mode='w')
        else:
            if self.is_grad:
                result = MV()
                for (coord, brecp) in zip(self.coords, self.rcpr_bases_MV):
                    result += brecp * diff(b, coord)
                return result
            else:
                return self * b

    def __rxor__(self, b):  # b ^ self
        b_W_self = MV(self)
        b_W_self.obj = b * self.obj
        return b_W_self

    def scalar(self):
        (coefs, blades) = linear_expand(self.obj)
        result = S.Zero
        for (coef, blade) in zip(coefs, blades):
            if MV.blade_grades[blade] == 0:
                result += coef
        return result

    def set_coef(self, igrade, ibase, value):
        if self.blade_rep:
            base = MV.blades[igrade][ibase]
        else:
            base = MV.bases[igrade][ibase]
        (coefs, bases) = linear_expand(self.obj)
        bases_lst = list(bases)  # python 2.5
        if base in bases:
            self.obj += (value - coefs[bases_lst.index(base)]) * base
        else:
            self.obj += value * base
        return

    def grade(self, igrade=0):
        if igrade > MV.dim:
            return MV()
        if self.igrade > -1:
            if self.igrade == igrade:
                return self
            else:
                return MV()
        else:
            (coefs, blades) = linear_expand(self.obj)
            result = S.Zero
            for (coef, blade) in zip(coefs, blades):
                if MV.blade_grades[blade] == igrade:
                    result += coef * blade
            self_igrade = MV(result)
            self_igrade.igrade = igrade
            return self_igrade

    def get_grades(self):  # grade decomposition of multivector
        self.base_to_blade()
        (coefs, bases) = linear_expand(self.obj)
        grades = {}
        for (coef, base) in zip(coefs, bases):
            igrade = MV.blade_grades[base]
            if igrade in grades:
                grades[igrade] += coef * base
            else:
                grades[igrade] = coef * base
        for key in grades:  # convert sympy expression to multivector
            grade = MV(grades[key])
            grade.blad_rep = True
            grade.igrade = key
            grades[key] = grade
        return grades

    def discover_and_set_grade(self):
        self.base_to_blade()
        (coefs, bases) = linear_expand(self.obj)
        old_grade = -1
        first_flg = True
        for (coef, base) in zip(coefs, bases):
            igrade = MV.blade_grades[base]
            if igrade != old_grade and first_flg:
                first_flg = False
                old_grade = igrade
            elif igrade != old_grade:
                self.igrade = -1
                return
        self.igrade = old_grade
        return

    def get_normal_order_str(self):
        self.obj = expand(self.obj)
        if self.blade_rep:
            self.obj = self.obj.collect(MV.blades_flat1)
        else:
            self.obj = self.obj.collect(MV.bases_flat1)
        terms = zip(*linear_expand(self.obj))
        if self.blade_rep:
            terms = sorted(terms, key=lambda x: MV.blades_flat1_lst.index(x[1]))  # Python 2.5
        else:
            terms = sorted(terms, key=lambda x: MV.bases_flat1_lst.index(x[1]))  # Python 2.5
        o_str = ''
        first = True
        if self.fmt == 2:
            if terms[0][1] == MV.ONE:
                grade = 0
            else:
                s = str(factor_terms(terms[0][0]))
                grade = max(s.count('^'), s.count('*')) + 1
        for term in terms:
            if term[1] == MV.ONE:
                tmp = str(factor_terms(term[0]))
            else:
                if isinstance(term[0], Add):
                    tmp = str(factor_terms(term[0]))
                    tmp = '(' + tmp + ')*' + enhance_print.enhance_base(str(term[1]))
                else:
                    coef_str = str(factor_terms(term[0]))
                    if coef_str == '1':
                        coef_str = ''
                    elif coef_str == '-1':
                        coef_str = '-'
                    else:
                        coef_str += '*'
                    tmp = coef_str + enhance_print.enhance_base(str(term[1]))
            if first:
                first = False
                o_str += tmp
            else:
                nl = ''
                if self.fmt == 2:
                    s = str(term[1])
                    new_grade = max(s.count('^'), s.count('*')) + 1
                    if new_grade > grade:
                        nl = '\n'
                        grade = new_grade
                if tmp[0] == '-':
                    o_str += nl + ' - ' + tmp[1:]
                else:
                    o_str += nl + ' + ' + tmp
            if self.fmt == 3:
                o_str += '\n'

        if o_str[-1] == '\n':
            o_str = o_str[:-1]
        return o_str

    def get_latex_normal_order_str(self):

        latex_sep = {'^': r'\W ', '*': ' '}

        def base_string(base_obj):
            base_str = GA_LatexPrinter.Basic__str__(base_obj)
            sep = '^'
            if '*' in base_str:
                sep = '*'
            base_lst = base_str.split(sep)

            lstr = r'\bm{' + latex(Symbol(base_lst[0]))
            for base in base_lst[1:]:
                lstr += latex_sep[sep] + latex(Symbol(base))
            lstr += '}'
            return lstr

        self.obj = expand(self.obj)
        if self.blade_rep:
            self.obj = self.obj.collect(MV.blades_flat)
        else:
            self.obj = self.obj.collect(MV.bases_flat)
        terms = zip(*linear_expand(self.obj))
        if self.blade_rep:
            bgrades = MV.blade_grades
            terms = sorted(terms, key=lambda x: MV.blades_flat1_lst.index(x[1]))  # Python 2.5
        else:
            bgrades = MV.base_grades
            terms = sorted(terms, key=lambda x: MV.bases_flat1_lst.index(x[1]))  # Python 2.5
        grades = []
        bases = []
        old_grade = -1
        for term in terms:
            new_grade = bgrades[term[1]]
            if old_grade != new_grade:
                if old_grade > -1:
                    grades.append(bases)
                    bases = []
                old_grade = new_grade
            bases.append(term)
        if len(bases) > 0:
            grades.append(bases)

        o_str = ''
        grade_strs = []
        nbases = 0
        for grade in grades:  # create [grade[base]] list of base strings
            base_strs = []
            for base in grade:
                nbases += 1
                if base[1] == MV.ONE:
                    base_str = latex(simplify(base[0]))
                else:
                    if isinstance(base[0], Add):
                        base_str = r'\left ( ' + latex(simplify(base[0])) + r'\right ) ' + base_string(base[1])
                    else:
                        coef_str = latex(simplify(base[0]))
                        if coef_str == '1':
                            coef_str = ''
                        elif coef_str == '-1':
                            coef_str = '-'
                        base_str = coef_str + base_string(base[1])
                    if base_str[0] != '-':
                        base_str = '+' + base_str
                base_strs.append(base_str)
            grade_strs.append(base_strs)
        if grade_strs[0][0][0] == '+':
            grade_strs[0][0] = grade_strs[0][0][1:]

        o_str = ''
        ngrades = len(grade_strs)

        if (self.fmt == 2 and ngrades > 1) or (self.fmt == 3 and nbases > 1):
            o_str += '\\begin{align*} '

        for base_strs in grade_strs:
            if self.fmt == 2 and ngrades > 1:
                o_str += ' & '
            for base_str in base_strs:
                if self.fmt == 3 and nbases > 1:
                    o_str += ' & ' + base_str + ' \\\\ '
                else:
                    o_str += base_str
            if self.fmt == 2 and ngrades > 1:
                o_str += ' \\\\ '

        if (self.fmt == 2 and ngrades > 1) or (self.fmt == 3 and nbases > 1):
            o_str += '\\end{align*} \n'
        else:
            pass
        return o_str

    @staticmethod
    def characterize_expression(self, expr):
        (coefs, bases) = linear_expand(expr)
        self.obj = expr
        if not MV.is_orthogonal:
            if len(set(bases) & MV.bases_set) != 0:
                self.blade_rep = False
                self.igrade = -1
            else:
                self.blade_rep = True
                self.igrade = 0
                return self
        else:
            self.blade_rep = True
            self.igrade = -1

        if self.blade_rep:
            self.igrade = MV.blade_grades[bases[0]]
            for base in bases[1:]:
                igrade = MV.blade_grades[base]
                if self.igrade != igrade:
                    self.igrade = -1
                    break
            return self
        return self

    def db(self):
        print('(blade_rep,igrade,obj) =', self.blade_rep, self.igrade, self.obj)
        return

    ##########################Member Functions##########################

    def dd(self, v):
        (coefs, bases) = linear_expand(v.obj)
        dderiv = MV()
        for (coef, base) in zip(coefs, bases):
            dderiv += coef * self.diff(MV.dd_dict[base])
        return dderiv

    def diff(self, var):
        dself = MV(self)
        dself.obj = diff(self.obj, var)
        return dself

    def simplify(self):
        (coefs, bases) = linear_expand(self.obj)
        obj = 0
        for (coef, base) in zip(coefs, bases):
            coef = simplify(coef)
            obj += coef * base
        sself = MV(self)
        sself.obj = obj
        return sself

    def trigsimp(self, **kwargs):
        (coefs, bases) = linear_expand(self.obj)
        obj = 0
        for (coef, base) in zip(coefs, bases):
            coef = trigsimp(coef, **kwargs)
            obj += coef * base
        ts_self = MV(self)
        ts_self.obj = obj
        return ts_self

    def exp(self, alpha=1, norm=0, mode='T'):
        if self.is_blade():
            self_sq = (self * self).scalar()
            if mode == 'T':
                if norm == 0:
                    norm = sqrt(-self_sq)
                return cos(alpha * norm) + sin(alpha * norm) * self / norm
            else:
                if norm == 0:
                    norm = sqrt(self_sq)
                return cosh(alpha * norm) + sinh(alpha * norm) * self / norm
        else:
            raise TypeError('!!! ' + str(self) + ' is not a blade in member function "exp" !!!\n')

    def expand(self):
        xself = MV(self)
        xself.obj = expand(self.obj)
        return xself

    def factor(self):
        fself = MV(self)
        fself.obj = factor_terms(self.obj)
        return fself

    def subs(self, x):
        xsubs = self.obj.subs(x)
        return MV(xsubs)

    def collect(self, x):
        (coefs, bases) = linear_expand(self.obj)
        result = S.Zero
        for (coef, base) in zip(coefs, bases):
            result += collect(coef, x) * base
        return MV(result)

    def is_scalar(self):
        self.discover_and_set_grade()
        if self.igrade == 0:
            return True
        else:
            return False

    def is_blade(self):
        self_sq = self * self
        if self_sq.is_scalar():
            return True
        return False

    def dual(self):
        dself = MV.I * self
        dself.discover_and_set_grade()
        return dself

    def even(self):
        if self.igrade > -1:
            if self.igrade % 2 == 0:
                return self
            else:
                return MV()
        else:
            (coefs, blades) = linear_expand(self.obj)
            result = S.Zero
            for (coef, blade) in zip(coefs, blades):
                if MV.blade_grades[blade] % 2 == 0:
                    result += coef * blade
            return MV(result)

    def odd(self):
        if self.igrade > -1:
            if self.igrade % 2 == 1:
                return self
            else:
                return MV()
        else:
            (coefs, blades) = linear_expand(self.obj)
            result = S.Zero
            for (coef, blade) in zip(coefs, blades):
                if MV.blade_grades[blade] % 2 == 1:
                    result += coef * blade
            return MV(result)

    def norm(self):
        norm_sq = self * self.rev()
        norm_sq.discover_and_set_grade()
        if norm_sq.igrade == 0:
            norm_sq = norm_sq.scalar()
            if isinstance(norm_sq, Number):
                norm = sqrt(abs(norm_sq))
            else:
                norm = sqrt(abs(norm_sq))
            return norm
        else:
            raise ValueError('In norm self*self.rev() = ' + str(norm_sq) +
                             ' is not a scalar!\n')

    def norm2(self):
        norm_sq = self * self.rev()
        norm_sq.discover_and_set_grade()
        if norm_sq.igrade == 0:
            norm_sq = norm_sq.scalar()
            return norm_sq
        else:
            raise ValueError('In norm self*self.rev() = ' + str(norm_sq) +
                             ' is not a scalar!\n')

    def rev(self):
        self.base_to_blade()
        (coefs, bases) = linear_expand(self.obj)
        result = S.Zero
        for (coef, base) in zip(coefs, bases):
            grade = MV.blade_grades[base]
            if grade < 2:
                result += coef * base
            else:
                sgn_pow = (grade * (grade - 1)) / 2 % 2
                if sgn_pow == 1:
                    result -= coef * base
                else:
                    result += coef * base
        self_rev = MV()
        self_rev.obj = result
        self_rev.igrade = self.igrade
        self_rev.blade_rep = self.blade_rep
        self_rev.fct = self.fct
        self_rev.is_grad = self.is_grad
        self_rev.print_blades = MV.print_blades
        self_rev.obj = simplify(self_rev.obj)
        return self_rev

    def inv(self):
        self_rev = self.rev()
        norm = self * self_rev
        norm.obj = expand(norm.obj)
        norm.discover_and_set_grade()
        if norm.igrade == 0:
            return self_rev / norm.obj
        else:
            raise ValueError('Cannot take inv(A) since A*rev(A) = ' + str(norm) +
                             ' is not a scalar.\n')

    #######################Reduce Combined Indexes######################

    @staticmethod
    def reduce_basis_loop(blst):
        """
        blst is a list of integers [i_{1},...,i_{r}] representing the geometric
        product of r basis vectors a_{{i_1}}*...*a_{{i_r}}.  reduce_basis_loop
        searches along the list [i_{1},...,i_{r}] untill it finds i_{j} == i_{j+1}
        and in this case contracts the list, or if i_{j} > i_{j+1} it revises
        the list (~i_{j} means remove i_{j} from the list)

        Case 1: If i_{j} == i_{j+1}, return a_{i_{j}}**2 and
                [i_{1},..,~i_{j},~i_{j+1},...,i_{r}]

        Case 2: If i_{j} > i_{j+1}, return a_{i_{j}}.a_{i_{j+1}},
                [i_{1},..,~i_{j},~i_{j+1},...,i_{r}], and
                [i_{1},..,i_{j+1},i_{j},...,i_{r}]
        """
        nblst = len(blst)  # number of basis vectors
        if nblst <= 1:
            return True  # a scalar or vector is already reduced
        jstep = 1
        while jstep < nblst:
            istep = jstep - 1
            if blst[istep] == blst[jstep]:  # basis vectorindex is repeated
                i = blst[istep]  # save basis vector index
                if len(blst) > 2:
                    blst = blst[:istep] + blst[jstep + 1:]  # contract blst
                else:
                    blst = []
                if len(blst) <= 1 or jstep == nblst - 1:
                    blst_flg = True  # revision of blst is complete
                else:
                    blst_flg = False  # more revision needed
                return MV.metric[i, i], blst, blst_flg
            if blst[istep] > blst[jstep]:  # blst not in normal order
                blst1 = blst[:istep] + blst[jstep + 1:]  # contract blst
                a1 = MV.metric2[blst[jstep], blst[istep]]  # coef of contraction
                blst = blst[:istep] + [blst[jstep]] + [blst[istep]] + blst[jstep + 1:]  # revise blst
                if len(blst1) <= 1:
                    blst1_flg = True  # revision of blst is complete
                else:
                    blst1_flg = False  # more revision needed
                return a1, blst1, blst1_flg, blst
            jstep += 1
        return True  # revision complete, blst in normal order

    @staticmethod
    def reduce_basis(blst):
        """
        Repetitively applies reduce_basis_loop to blst
        product representation until normal form is realized.
        """
        if blst == []:  # blst represents scalar
            blst_coef = [S.One]
            blst_expand = [[]]
            return blst_coef, blst_expand
        blst_expand = [blst]
        blst_coef = [S.One]
        blst_flg = [False]
        # reduce untill all blst revise flgs are True
        while not reduce(operator.and_, blst_flg):
            for i in range(len(blst_flg)):
                if not blst_flg[i]:  # keep revising if revise flg is False
                    tmp = MV.reduce_basis_loop(blst_expand[i])
                    if isinstance(tmp, bool):
                        blst_flg[i] = tmp  # revision of blst_expand[i] complete
                    elif len(tmp) == 3:  # blst_expand[i] contracted
                        blst_coef[i] = tmp[0] * blst_coef[i]
                        blst_expand[i] = tmp[1]
                        blst_flg[i] = tmp[2]
                    else:  # blst_expand[i] revised
                        blst_coef[i] = -blst_coef[i]
                        # if revision force one more pass in case revision
                        # causes repeated index previous to revised pair of
                        # indexes
                        blst_flg[i] = False
                        blst_expand[i] = tmp[3]
                        blst_coef.append(-blst_coef[i] * tmp[0])
                        blst_expand.append(tmp[1])
                        blst_flg.append(tmp[2])
        new_blst_coef = []
        new_blst_expand = []
        for (coef, expand) in zip(blst_coef, blst_expand):
            if expand in new_blst_expand:
                i = new_blst_expand.index(expand)
                new_blst_coef[i] += coef
            else:
                new_blst_expand.append(expand)
                new_blst_coef.append(coef)
        return new_blst_coef, new_blst_expand

    ##########################Bases Construction########################

    @staticmethod
    def symbol_product_bases(i1, i2):
        if i1 == ():
            if i2 == ():
                return S.One
            else:
                return MV.index_to_base[i2]
        else:
            if i2 == ():
                return MV.index_to_base[i1]

        index = list(i1 + i2)
        result = S.Zero
        (coefs, indexes) = MV.reduce_basis(index)
        for (coef, index) in zip(coefs, indexes):
            result += coef * MV.index_to_base[tuple(index)]
        return result

    @staticmethod
    def make_base_blade_symbol(ibase):
        if len(ibase) == 1:
            base_str = MV.basis_names[ibase[0]]
            return Symbol(base_str, commutative=False), base_str, \
                Symbol(base_str, commutative=False), base_str
        else:
            base_str = ''
            blade_str = ''
            for index in ibase:
                vector_str = MV.basis_names[index]
                base_str += vector_str + '*'
                blade_str += vector_str + '^'
            base_str = base_str[:-1]
            blade_str = blade_str[:-1]
            return Symbol(base_str, commutative=False), base_str, \
                Symbol(blade_str, commutative=False), blade_str

    ################Geometric, Wedge, and Dot Products##################

    @staticmethod
    def basic_geometric_product(obj1, obj2):
        """
        basic_geometric_product assumes that mv1 and mv2 are both
        mulitvectors, not scalars and both are in the base and not the
        blade representation.  No multivector flags are checked.
        This function is used to construct the blades from the bases.
        """
        def mul_table(b1, b2):
            return MV.base_mul_table[(b1, b2)]

        obj12 = bilinear_product(obj1 * obj2, mul_table)

        return obj12

    @staticmethod
    def geometric_product(b1, b2):
        """
        geometric_product(b1, b2) calculates the geometric
        product of the multivectors b1 and b2 (b1*b2).
        """
        if MV.is_orthogonal:
            return MV.product_orthogonal_blades(b1, b2)
        else:
            result = MV.base_mul_table[(b1, b2)]
            return result

    @staticmethod
    def basic_add(mv1, mv2):
        """
        basic_add assummes that mv1 and mv2 are multivectors both in the
        base or blade representation. It sets no flags for the output
        and forms mv1.obj+mv2.obj.  It is used to form the base expansion
        of the blades.
        """
        obj = expand(mv1.obj + mv2.obj)
        return MV(obj)

    @staticmethod
    def basic_sub(mv1, mv2):
        """
        basic_sub assummes that mv1 and mv2 are multivectors both in the
        base or blade representation. It sets no flags for the output
        and forms mv1.obj-mv2.obj.  It is used to form the base expansion
        of the blades.
        """
        obj = expand(mv1.obj - mv2.obj)
        return MV(obj)

    @staticmethod
    def dot_product(b1, b2):
        if MV.is_orthogonal:
            return MV.dot_orthogonal_blades(b1, b2)
        else:
            grade1 = MV.blade_grades[b1]
            grade2 = MV.blade_grades[b2]
            if MV.dot_mode == 's':  # dot product
                return MV.blade_dot_table[(b1, b2)]
            elif MV.dot_mode == 'l':  # left contraction
                grade = grade2 - grade1
            elif MV.dot_mode == 'r':  # right contraction
                grade = grade1 - grade2
            if grade < 0:
                return MV()
            else:
                return MV.blade_dot_table[(b1, b2)]

    @staticmethod
    def wedge_product(b1, b2):
        i1 = MV.blade_to_index[b1]
        i2 = MV.blade_to_index[b2]
        i1_plus_i2 = list(i1 + i2)
        if len(i1_plus_i2) > MV.dim:
            return S.Zero
        (sgn, i1_W_i2) = MV.blade_reduce(i1_plus_i2)
        if sgn != 0:
            return sgn * MV.index_to_blade[tuple(i1_W_i2)]
        else:
            return S.Zero

    @staticmethod
    def blade_reduce(lst):
        sgn = 1
        for i in range(1, len(lst)):
            save = lst[i]
            j = i
            while j > 0 and lst[j - 1] > save:
                sgn = -sgn
                lst[j] = lst[j - 1]
                j -= 1
            lst[j] = save
            if lst[j] == lst[j - 1]:
                return 0, None
        return sgn, lst

    @staticmethod
    def non_orthogonal_products(mv1, mv2, mode='w'):
        if isinstance(mv1, MV) and isinstance(mv2, MV):  # both sides are mv
            mv1_grades = mv1.get_grades()
            mv2_grades = mv2.get_grades()
            result = MV()
            for grade1 in mv1_grades:
                for grade2 in mv2_grades:
                    if mode == 'w':  # wedge product
                        grade = grade1 + grade2
                    elif mode == 's':  # dot product
                        if grade1 == 0:
                            grade = -1
                        elif grade2 == 0:
                            grade = -1
                        else:
                            grade = abs(grade1 - grade2)
                    elif mode == 'l':  # left contraction
                        grade = grade2 - grade1
                    elif mode == 'r':  # right contraction
                        grade = grade1 - grade2
                    if grade >= 0 and grade <= MV.dim:
                        mv1mv2 = mv1_grades[grade1] * mv2_grades[grade2]
                        mv1mv2_grades = MV(mv1mv2).get_grades()
                        if grade in mv1mv2_grades:
                            result += mv1mv2_grades[grade]
            return result
        elif isinstance(mv1, MV):  # rhs is sympy scalar
            if mode == 'w':  # wedge product
                mv1mv2 = MV(mv1)
                mv1mv2.obj = mv2 * mv1.obj
                return mv1mv2
            else:  # dot product or contractions
                return MV()
        elif isinstance(mv2, MV):  # lhs is sympy scalar
                mv1mv2 = MV(mv1)
                mv1mv2.obj = mv2 * mv1.obj
                return mv1mv2
        else:  # both sides are sympy scalars
            if mode == 'w':
                return MV(mv1 * mv2)
            else:
                return MV()

    ###################Blade Base conversion functions##################

    def blade_to_base(self):
        if MV.is_orthogonal:
            return
        if self.igrade == 0 or self.igrade == 1:
            return self
        if self.blade_rep:
            self.blade_rep = False
            self.obj = expand(self.obj)
            self.obj = self.obj.subs({S.One**2: S.One})
            self.obj = simplify(self.obj.subs(MV.blade_expand))

            return

    def base_to_blade(self):
        if MV.is_orthogonal:
            return
        if self.igrade == 0 or self.igrade == 1:
            return
        if not self.blade_rep:
            self.blade_rep = True
            self.obj = expand(self.obj)
            self.obj = self.obj.subs(MV.base_expand)
            self.obj = expand(self.obj)
            self.obj = simplify(self.obj)
            return

    @staticmethod
    def build_base_blade_arrays(debug):
        indexes = tuple(range(MV.dim))
        MV.index = [()]
        for i in indexes:
            MV.index.append(tuple(combinations(indexes, i + 1)))
        MV.index = tuple(MV.index)

        # Set up base and blade and index arrays

        if not MV.is_orthogonal:
            MV.bases_flat = []
            MV.bases = [MV.ONE]
            MV.base_to_index = {MV.ONE: ()}
            MV.index_to_base = {(): MV.ONE}
            MV.base_grades = {MV.ONE: 0}
            MV.base_grades[S.One] = 0

        MV.blades = [MV.ONE]
        MV.blades_flat = []
        MV.blade_grades = {MV.ONE: 0}
        MV.blade_grades[S.One] = 0
        MV.blade_to_index = {MV.ONE: ()}
        MV.index_to_blade = {(): MV.ONE}

        ig = 1  # pseudo grade and grade index
        for igrade in MV.index[1:]:
            if not MV.is_orthogonal:
                bases = []  # base symbol array within pseudo grade
            blades = []  # blade symbol array within grade
            ib = 0  # base index within grade
            for ibase in igrade:
                # build base name string
                (base_sym, base_str, blade_sym, blade_str) = MV.make_base_blade_symbol(ibase)

                if not MV.is_orthogonal:
                    bases.append(base_sym)
                    MV.bases_flat.append(base_sym)

                blades.append(blade_sym)
                MV.blades_flat.append(blade_sym)
                base_index = MV.index[ig][ib]

                # Add to dictionarys relating symbols and indexes
                if not MV.is_orthogonal:
                    MV.base_to_index[base_sym] = base_index
                    MV.index_to_base[base_index] = base_sym
                    MV.base_grades[base_sym] = ig

                MV.blade_to_index[blade_sym] = base_index
                MV.index_to_blade[base_index] = blade_sym
                MV.blade_grades[blade_sym] = ig

                ib += 1
            ig += 1

            if not MV.is_orthogonal:
                MV.bases.append(tuple(bases))

            MV.blades.append(tuple(blades))

        if not MV.is_orthogonal:
            MV.bases = tuple(MV.bases)
            MV.bases_flat = tuple(MV.bases_flat)
            MV.bases_flat1 = (MV.ONE, ) + MV.bases_flat
            MV.bases_set = set(MV.bases_flat[MV.dim:])

            MV.bases_flat1_lst = list(MV.bases_flat1)  # Python 2.5

        MV.blades = tuple(MV.blades)
        MV.blades_flat = tuple(MV.blades_flat)
        MV.blades_flat1 = (MV.ONE, ) + MV.blades_flat
        MV.blades_set = set(MV.blades_flat[MV.dim:])

        MV.blades_flat1_lst = list(MV.blades_flat1)  # Python 2.5

        if debug:
            if not MV.is_orthogonal:
                oprint('MV Class Global Objects:', None,
                       'index(tuple)', MV.index,
                       'bases(Symbol)', MV.bases,
                       'base_to_index(Symbol->tuple)', MV.base_to_index,
                       'index_to_base(tuple->Symbol)', MV.index_to_base,
                       'bases flat', MV.bases_flat,
                       'bases set', MV.bases_set,
                       'blades(Symbol)', MV.blades,
                       'blade_grades(int)', MV.blade_grades,
                       'blade_to_index(Symbol->tuple)', MV.blade_to_index,
                       'index_to_blade(tuple->Symbol)', MV.index_to_blade,
                       'blades flat', MV.blades_flat,
                       'blades set', MV.blades_set, dict_mode=True)
            else:
                oprint('MV Class Global Objects:', None,
                       'index(tuple)', MV.index,
                       'blades(Symbol)', MV.blades,
                       'blades flat', MV.blades_flat,
                       'blades set', MV.blades_set,
                       'blade_grades(int)', MV.blade_grades,
                       'blade_to_index(Symbol->tuple)', MV.blade_to_index,
                       'index_to_blade(tuple->Symbol)', MV.index_to_blade, dict_mode=True)
        return

    @staticmethod
    def build_base_mul_table(debug):

        # Calculate geometric product multiplication table for bases

        MV.base_mul_table = {(MV.ONE, MV.ONE): MV.ONE}

        for ig1 in MV.index[1:]:
            for ib1 in ig1:
                b1 = MV.index_to_base[ib1]
                MV.base_mul_table[(MV.ONE, b1)] = b1
                MV.base_mul_table[(b1, MV.ONE)] = b1
                for ig2 in MV.index[1:]:
                    for ib2 in ig2:
                        b2 = MV.index_to_base[ib2]
                        b1b2 = MV.symbol_product_bases(ib1, ib2)
                        key = (b1, b2)
                        MV.base_mul_table[key] = simplify(b1b2)

        if debug:
            oprint('Geometric Product (*) Table for Bases', MV.base_mul_table, dict_mode=True)
        return

    @staticmethod
    def build_base_blade_expansion_tables(debug):

        # Expand blades in terms of bases

        MV.blade_expand = {}
        for (blade, base) in zip(MV.blades[1], MV.bases[1]):
            MV.blade_expand[blade] = base

        sgn = -S.One
        igrade = 2
        while igrade <= MV.dim:
            for ibase in MV.index[igrade]:
                pre_index = (ibase[0], )
                post_index = ibase[1:]
                a = MV.index_to_blade[pre_index]
                B = MV.index_to_blade[post_index]
                B_expand = MV.blade_expand[B]
                # a^B = (a*B+sgn*B*a)/2
                if sgn == 1:
                    result = MV.basic_geometric_product(a, B_expand) + MV.basic_geometric_product(B_expand, a)
                else:
                    result = MV.basic_geometric_product(a, B_expand) - MV.basic_geometric_product(B_expand, a)
                MV.blade_expand[MV.index_to_blade[ibase]] = simplify(expand(result / S(2)))
            igrade += 1
            sgn = -sgn

        if debug:
            oprint('Blade Expansion Table', MV.blade_expand, dict_mode=True)

        # Expand bases in terms of blades

        MV.base_expand = {}
        for (blade, base) in zip(MV.blades[1], MV.bases[1]):
            MV.base_expand[base] = blade

        ig = 2
        while ig <= MV.dim:
            tmp_dict = {}
            for ib in MV.index[ig]:
                base = MV.index_to_base[ib]
                blade = MV.index_to_blade[ib]
                tmp = MV.blade_expand[blade]
                tmp = tmp.subs({base: -blade})
                tmp = tmp.subs(MV.base_expand)
                tmp_dict[base] = simplify(expand(-tmp))
            MV.base_expand.update(tmp_dict)
            ig += 1

        if debug:
            oprint('Base Expansion Table', MV.base_expand, dict_mode=True)

            print('Test Blade Expansion:')
            for key in MV.blade_expand:
                test = MV.blade_expand[key].subs(MV.base_expand)
                print(str(key) + ' = ' + str(test))
            print('Test Base Expansion:')
            for key in MV.base_expand:
                test = MV.base_expand[key].subs(MV.blade_expand)
                print(str(key) + ' = ' + str(test))

        return

    @staticmethod
    def build_reciprocal_basis(debug):
        MV.I = MV(MV.blades_flat[-1])
        MV.rcpr_norm = get_commutative_coef(simplify((MV.I * MV.I).obj))
        duals = list(MV.blades_flat[-(MV.dim + 1): -1])
        duals.reverse()
        sgn = 1
        MV.rcpr_bases_MV = []
        for dual in duals:
            recpv = sgn * MV(dual) * MV.I
            MV.rcpr_bases_MV.append(recpv)
            sgn = -sgn

        if debug:
            print('Reciprocal Norm =', MV.rcpr_norm)
            oprint('Reciprocal Basis', MV.rcpr_bases_MV)

        if MV.coords is not None:
            rcpr_bases_MV = []
            MV.grad = MV()

            result = S.Zero
            for (coef, rbase) in zip(MV.coords, MV.rcpr_bases_MV):
                nbase_obj = rbase.obj / MV.rcpr_norm
                term = coef * rbase
                result += term
                rcpr_bases_MV.append(MV(nbase_obj))

            MV.dd_dict = {}

            if MV.is_orthogonal:
                bases = MV.blades[1]
            else:
                bases = MV.bases[1]

            for (coord, base) in zip(MV.coords, bases):
                MV.dd_dict[base] = coord

            MV.rcpr_bases_MV = rcpr_bases_MV

            MV.grad.is_grad = True
            MV.grad.blade_rep = True
            MV.grad.igrade = 1
            MV.grad.rcpr_bases_MV = tuple(rcpr_bases_MV)
            MV.grad.coords = MV.coords
            MV.grad.norm = MV.rcpr_norm
            MV.grad.connection = {}

            if debug:
                print('grad =', MV.grad)
                oprint('reciprocal bases', MV.rcpr_bases_MV)

        if debug and MV.coords is not None and not MV.is_orthogonal:
            print('Reciprocal Vector Test:')
            for v1 in MV.blades_MV:
                for (v2, rv2) in zip(MV.blades_MV, MV.rcpr_bases_MV):
                    print(str(v1) + '|Reciprocal(' + str(v2) + ') = ' + str(simplify(expand((v1 | rv2).obj)) / MV.rcpr_norm))
            print('I**2 =', MV.rcpr_norm)
            print('Grad Vector:', MV.grad)

        return

    @staticmethod
    def build_curvilinear_connection(debug):
        """
        Vector.dtau_dict[basis vector symbol,coordinate symbol] = derivative of basis vector as sympy expression
        """

        MV.connection = True
        MV.tangent_norm = Vector.norm
        MV.tangent_derivatives_MV = {}

        rcpr_bases_MV = []

        for (rbase, norm) in zip(MV.rcpr_bases_MV, MV.tangent_norm):
            rcpr_bases_MV.append(rbase / norm)

        MV.rcpr_bases_MV = rcpr_bases_MV

        for key in Vector.dtau_dict.keys():
            MV.tangent_derivatives_MV[key] = MV(Vector.dtau_dict[key])

        if debug:
            oprint('Tangent Vector Derivatives', MV.tangent_derivatives_MV, dict_mode=True)

        MV.left_connection = {}
        MV.right_connection = {}
        MV.left_wedge_connection = {}
        MV.right_wedge_connection = {}
        MV.left_dot_connection = {}
        MV.right_dot_connection = {}

        for base in MV.blades[1]:
            right_result = MV()
            left_result = MV()
            for (coord, rblade) in zip(MV.coords, MV.rcpr_bases_MV):
                left_result += rblade * MV.tangent_derivatives_MV[(base, coord)]
                right_result += MV.tangent_derivatives_MV[(base, coord)] * rblade
            left_result.obj = expand(left_result.obj)
            right_result.obj = expand(right_result.obj)
            left_result.discover_and_set_grade()
            right_result.discover_and_set_grade()
            MV.left_connection[base] = left_result.obj
            MV.right_connection[base] = right_result.obj
            MV.left_wedge_connection[base] = left_result.grade(2).obj
            MV.right_wedge_connection[base] = right_result.grade(2).obj
            MV.left_dot_connection[base] = left_result.grade(0).obj
            MV.right_dot_connection[base] = right_result.grade(0).obj

        for grade in MV.blades[2:]:
            for blade in grade:
                index = MV.blade_to_index[blade]
                N = len(index)
                left_result = MV()
                right_result = MV()
                for (coord, rblade) in zip(MV.coords, MV.rcpr_bases_MV):
                    tmp = MV()
                    for i in range(N):
                        i_pre = index[:i]
                        i_dtan = index[i]
                        i_post = index[i + 1:]
                        base = MV.blades[1][i_dtan]
                        tmp += (MV(MV.index_to_blade[i_pre]) ^ MV.tangent_derivatives_MV[(base, coord)]) ^ MV(MV.index_to_blade[i_post])
                    left_result += rblade * tmp
                    right_result += tmp * rblade
                left_result.discover_and_set_grade()
                right_result.discover_and_set_grade()
                MV.left_connection[blade] = left_result.obj
                MV.right_connection[blade] = right_result.obj
                MV.left_wedge_connection[blade] = left_result.grade(N + 1).obj
                MV.right_wedge_connection[blade] = right_result.grade(N + 1).obj
                MV.left_dot_connection[blade] = left_result.grade(abs(N - 1)).obj
                MV.right_dot_connection[blade] = right_result.grade(abs(N - 1)).obj

        MV.connection = {'left': MV.left_connection, 'right': MV.right_connection,
                         'left_wedge': MV.left_wedge_connection, 'right_wedge': MV.right_wedge_connection,
                         'left_dot': MV.left_dot_connection, 'right_dot': MV.right_dot_connection}
        MV.grad.connection = MV.connection

        if debug:
            oprint('Left Mutlivector Connection', MV.left_connection,
                   'Right Mutlivector Connection', MV.right_connection,
                   'Left Wedge Mutlivector Connection', MV.left_wedge_connection,
                   'Right Wedge Mutlivector Connection', MV.right_wedge_connection,
                   'Left Dot Mutlivector Connection', MV.left_dot_connection,
                   'Right Dot Mutlivector Connection', MV.right_dot_connection, dict_mode=True)
        return

    @staticmethod
    def setup(basis, metric=None, coords=None, rframe=False, debug=False, curv=(None, None)):
        """
        MV.setup() creates all the arrays and dictionaries required to construct, multiply, add,
        and differentiate multivectors in linear and curvilinear coordinate systems.  The inputs
        to MV.setup() are as follows -

            basis:  A string that defines the noncommutative symbols that represent the basis
                    vectors of the underlying vector space of the multivector space.  If the
                    string consists of substrings separated by spaces or commas each substring
                    will be the name of the basis vector symbol for example basis='e_x e_y e_z'
                    or basis='i j k'.  Another way to enter the basis symbols is to specify a
                    base string with a list of subscript strings.  This is done with the following
                    notation so that 'e*x|y|z' is equivalent to 'e_x e_y e_z'.

            metric:

            rframe:

            coords:

            debug:

            curv:

        """
        MV.print_blades = False
        MV.connection = False

        MV.ONE = ONE_NC

        MV.basis_vectors = Vector.setup(basis, metric=metric, coords=coords, curv=curv, debug=debug)
        MV.curv_norm = curv[1]
        MV.metric = Vector.metric
        MV.subscripts = Vector.subscripts
        MV.coords = Vector.coords
        MV.metric2 = 2 * Vector.metric
        MV.is_orthogonal = Vector.is_orthogonal

        MV.basis_names = []
        for base in MV.basis_vectors:
            MV.basis_names.append(str(base))

        if debug:
            oprint('Basis Names', MV.basis_names)

        MV.dim = len(MV.basis_vectors)
        MV.dim1 = MV.dim + 1

        MV.build_base_blade_arrays(debug)

        if not MV.is_orthogonal:
            MV.build_base_mul_table(debug)
            MV.build_base_blade_expansion_tables(debug)

        MV.blades_MV = []
        for b in MV.blades[1]:
            mv = MV()
            mv.obj = b
            mv.blade_rep = True
            mv.igrade = 1
            MV.blades_MV.append(mv)
        MV.build_reciprocal_basis(debug)

        MV.blades_MV = tuple(MV.blades_MV)

        if curv != (None, None):
            MV.build_curvilinear_connection(debug)

        MV.print_blades = True

        MV.I = MV(MV.blades_flat[-1])
        MV.Isq = simplify((MV.I * MV.I).scalar())
        MV.Iinv = MV.I / MV.Isq

        if coords is not None:
            return MV.blades_MV + (MV.grad, )
        else:
            return MV.blades_MV


def Format(Fmode=True, Dmode=True, ipy=False):
    "Initialize the LaTeX printer with the given mode."
    GA_LatexPrinter.Dmode = Dmode
    GA_LatexPrinter.Fmode = Fmode
    GA_LatexPrinter.ipy = ipy
    MV.latex_flg = True
    GA_LatexPrinter.redirect(ipy)
    return


def ga_print_on():
    """
    Turn on the galgebra-aware string printer.

    This function is intended for interactive use only.
    Use

      with GA_Printer():
        xxx

    instead of

      ga_print_on()
      xxx
      ga_print_off()
    """
    GA_Printer._on()
    return


def ga_print_off():
    """
    Turn off the galgebra-aware string printer.

    This function is intended for interactive use only.
    See ga_print_on for the noninteractive technique.
    """
    GA_Printer._off()
    return


def DD(v, f):
    if isinstance(f, MV):
        return f.dd(v)
    sf = MV(f, 'scalar')
    return sf.dd(v)


def Nga(x, prec=5):
    if isinstance(x, MV):
        Px = MV(x)
        Px.obj = Nsympy(x.obj, prec)
        return Px
    else:
        return Nsympy(x, prec)


def Com(A, B):
    "Commutator of A and B divided by 2."
    return (A * B - B * A) / S(2)


def inv(B):
    "Invert B if B*B.rev() is scalar."
    Bnorm = B * B.rev()
    if Bnorm.is_scalar():
        invB = B.rev() / Bnorm.obj
        return invB
    else:
        raise TypeError('Cannot calculate inverse of ' + str(B) + ' since \n'
                        + 'B*Brev() = ' + str(Bnorm) + ' is not a scalar.\n')


def proj(B, A):
    "Project blade A on blade B."
    result = (A < B) * inv(B)
    result.trigsimp()
    return result


def rotor(theta, n):
    n_sq = (n * n).obj
    if n_sq != S.One:
        n /= sqrt(n_sq)
    N = n.dual()
    R = cos(theta) + sin(theta) * N
    return R


def rot(itheta, A):
    "Rotate blade A by angle itheta."
    theta = itheta.norm()
    i = itheta / theta
    result = (cos(theta / 2) - i * sin(theta / 2)) * A * (cos(theta / 2) + i * sin(theta / 2))
    # result.trigsimp(recursive=True) #trigsimp doesn't work for half angle formulas
    return result


def refl(B, A):
    "Reflect blade A in blade B."
    j = B.is_blade()
    k = A.is_blade()
    if j > -1 and k > -1:
        result = (-1)**(j * (k + 1)) * B * A * inv(B)
        result.trigsimp()
        return result
    else:
        raise ValueError('Can only reflect blades')


def dual(M):
    return M * MV.Iinv


def cross(M1, M2):
    return -MV.I * (M1 ^ M2)


def ScalarFunction(TheFunction):
    return MV() + TheFunction


def ReciprocalFrame(basis, mode='norm'):
    dim = len(basis)

    indexes = tuple(range(dim))
    index = [()]

    for i in indexes[-2:]:
        index.append(tuple(combinations(indexes, i + 1)))

    MFbasis = []

    for igrade in index[-2:]:
        grade = []
        for iblade in igrade:
            blade = MV(1, 'scalar')
            for ibasis in iblade:
                blade ^= basis[ibasis]
            blade = blade.trigsimp(deep=True, recursive=True)
            grade.append(blade)
        MFbasis.append(grade)
    E = MFbasis[-1][0]
    E_sq = trigsimp((E * E).scalar(), deep=True, recursive=True)

    duals = copy.copy(MFbasis[-2])

    duals.reverse()
    sgn = 1
    rbasis = []
    for dual in duals:
        recpv = (sgn * dual * E).trigsimp(deep=True, recursive=True)
        rbasis.append(recpv)
        sgn = -sgn

    if mode != 'norm':
        rbasis.append(E_sq)
    else:
        for i in range(dim):
            rbasis[i] = rbasis[i] / E_sq

    return tuple(rbasis)
