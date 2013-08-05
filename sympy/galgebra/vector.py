# sympy/galgebra/vector.py

"""
vector.py is a helper class for the MV class that defines the basis
vectors and metric and calulates derivatives of the basis vectors for
the MV class.
"""

import itertools
import copy

from sympy import Symbol, S, Matrix, trigsimp, diff, expand

from sympy.galgebra.printing import GA_Printer
from sympy.galgebra.stringarrays import str_array
from sympy.galgebra.ncutil import linear_derivation, bilinear_product
from sympy.galgebra.debug import oprint


def flatten(lst):
    return list(itertools.chain(*lst))


def TrigSimp(x):
    return trigsimp(x, recursive=True)


class Vector(object):
    """
    Vector class.

    Setup is done by defining a set of basis vectors in static function 'Bases'.
    The linear combination of scalar (commutative) sympy quatities and the
    basis vectors form the vector space.  If the number of basis vectors
    is 'n' the metric tensor is formed as an n by n sympy matrix of scalar
    symbols and represents the dot products of pairs of basis vectors.
    """

    is_orthogonal = False

    @staticmethod
    def setup(base, n=None, metric=None, coords=None, curv=(None, None), debug=False):
        """
        Generate basis of vector space as tuple of vectors and
        associated metric tensor as Matrix.  See str_array(base,n) for
        usage of base and n and str_array(metric) for usage of metric.

        To overide elements in the default metric use the character '#'
        in the metric string.  For example if one wishes the diagonal
        elements of the metric tensor to be zero enter metric = '0 #,# 0'.

        If the basis vectors are e1 and e2 then the default metric -

            Vector.metric = ((dot(e1,e1),dot(e1,e2)),dot(e2,e1),dot(e2,e2))

        becomes -

            Vector.metric = ((0,dot(e1,e2)),(dot(e2,e1),0)).

        The function dot returns a Symbol and is symmetric.

        The functions 'Bases' calculates the global quantities: -

            Vector.basis
                tuple of basis vectors
            Vector.base_to_index
                dictionary to convert base to base inded
            Vector.metric
                metric tensor represented as a matrix of symbols and numbers

        """
        Vector.is_orthogonal = False
        Vector.coords = coords
        Vector.subscripts = []
        base_name_lst = base.split(' ')

        # Define basis vectors

        if '*' in base:
            base_lst = base.split('*')
            base = base_lst[0]
            Vector.subscripts = base_lst[1].split('|')
            base_name_lst = []
            for subscript in Vector.subscripts:
                base_name_lst.append(base + '_' + subscript)
        else:
            if len(base_name_lst) > 1:
                Vector.subscripts = []
                for base_name in base_name_lst:
                    tmp = base_name.split('_')
                    Vector.subscripts.append(tmp[-1])
            elif len(base_name_lst) == 1 and Vector.coords is not None:
                base_name_lst = []
                for coord in Vector.coords:
                    Vector.subscripts.append(str(coord))
                    base_name_lst.append(base + '_' + str(coord))
            else:
                raise TypeError("'%s' does not define basis vectors" % base)

        basis = []
        base_to_index = {}
        index = 0
        for base_name in base_name_lst:
            basis_vec = Vector(base_name)
            basis.append(basis_vec)
            base_to_index[basis_vec.obj] = index
            index += 1

        Vector.base_to_index = base_to_index
        Vector.basis = tuple(basis)

        # define metric tensor

        default_metric = []
        for bv1 in Vector.basis:
            row = []
            for bv2 in Vector.basis:
                row.append(Vector.basic_dot(bv1, bv2))
            default_metric.append(row)
        Vector.metric = Matrix(default_metric)
        if metric is not None:
            if metric[0] == '[' and metric[-1] == ']':
                Vector.is_orthogonal = True
                metric_str_lst = metric[1:-1].split(',')
                Vector.metric = []
                for g_ii in metric_str_lst:
                    Vector.metric.append(S(g_ii))
                Vector.metric = Matrix(Vector.metric)
            else:
                metric_str_lst = flatten(str_array(metric))
                for index in range(len(metric_str_lst)):
                    if metric_str_lst[index] != '#':
                        Vector.metric[index] = S(metric_str_lst[index])

        Vector.metric_dict = {}  # Used to calculate dot product
        N = range(len(Vector.basis))
        if Vector.is_orthogonal:
            for ii in N:
                    Vector.metric_dict[Vector.basis[ii].obj] = Vector.metric[ii]
        else:
            for irow in N:
                for icol in N:
                    Vector.metric_dict[(Vector.basis[irow].obj, Vector.basis[icol].obj)] = Vector.metric[irow, icol]

        # calculate tangent vectors and metric for curvilinear basis

        if curv != (None, None):
            X = S.Zero
            for (coef, base) in zip(curv[0], Vector.basis):
                X += coef * base.obj
            Vector.tangents = []
            for (coord, norm) in zip(Vector.coords, curv[1]):
                tau = diff(X, coord)
                tau = trigsimp(tau)
                tau /= norm
                tau = expand(tau)
                Vtau = Vector()
                Vtau.obj = tau
                Vector.tangents.append(Vtau)
            metric = []
            for tv1 in Vector.tangents:
                row = []
                for tv2 in Vector.tangents:
                    row.append(tv1 * tv2)
                metric.append(row)
            metric = Matrix(metric)
            metric = metric.applyfunc(TrigSimp)
            Vector.metric_dict = {}
            if metric.is_diagonal:
                Vector.is_orthogonal = True
                tmp_metric = []
                for ii in N:
                    tmp_metric.append(metric[ii, ii])
                    Vector.metric_dict[Vector.basis[ii].obj] = metric[ii, ii]
                Vector.metric = Matrix(tmp_metric)
            else:
                Vector.is_orthogonal = False
                Vector.metric = metric
                for irow in N:
                    for icol in N:
                        Vector.metric_dict[(Vector.basis[irow].obj, Vector.basis[icol].obj)] = Vector.metric[irow, icol]
            Vector.norm = curv[1]

            if debug:
                oprint('Tangent Vectors', Vector.tangents,
                       'Metric', Vector.metric,
                       'Metric Dictionary', Vector.metric_dict,
                       'Normalization', Vector.norm, dict_mode=True)

            # calculate derivatives of tangent vectors

            Vector.dtau_dict = None
            dtau_dict = {}

            for x in Vector.coords:
                for (tau, base) in zip(Vector.tangents, Vector.basis):
                    dtau = tau.diff(x).applyfunc(TrigSimp)
                    result = S.Zero
                    for (t, b) in zip(Vector.tangents, Vector.basis):
                        t_dtau = TrigSimp(t * dtau)
                        result += t_dtau * b.obj
                    dtau_dict[(base.obj, x)] = result

            Vector.dtau_dict = dtau_dict

            if debug:
                oprint('Basis Derivatives', Vector.dtau_dict, dict_mode=True)

        return tuple(Vector.basis)

    def __init__(self, basis_str=None):
        if isinstance(basis_str, Vector):
            self.obj = basis_str
        else:
            if basis_str is None or basis_str == '0':
                self.obj = S(0)
            else:
                self.obj = Symbol(basis_str, commutative=False)

    """
    def diff(self, x):
        (coefs, bases) = linear_expand(self.obj)
        result = S.Zero
        for (coef, base) in zip(coefs, bases):
            result += diff(coef, x) * base
        return result
    """

    def diff(self, x):
        Dself = Vector()
        if isinstance(Vector.dtau_dict, dict):
            Dself.obj = linear_derivation(self.obj, Vector.Diff, x)
        else:
            Dself.obj = diff(self.obj, x)
        return Dself

    @staticmethod
    def basic_dot(v1, v2):
        """
        Dot product of two basis vectors returns a Symbol
        """
        i1 = list(Vector.basis).index(v1)  # Python 2.5
        i2 = list(Vector.basis).index(v2)  # Python 2.5
        if i1 < i2:
            dot_str = '(' + str(Vector.basis[i1]) + '.' + str(Vector.basis[i2]) + ')'
        else:
            dot_str = '(' + str(Vector.basis[i2]) + '.' + str(Vector.basis[i1]) + ')'
        return Symbol(dot_str)

    @staticmethod
    def dot(b1, b2):
        if Vector.is_orthogonal:
            if b1 != b2:
                return S.Zero
            else:
                return Vector.metric_dict[b1]
        else:
            return Vector.metric_dict[(b1, b2)]

    @staticmethod
    def Diff(b, x):
        return Vector.dtau_dict[(b, x)]

    ######################## Operator Definitions#######################

    def __str__(self):
        return GA_Printer().doprint(self)

    def __mul__(self, v):
        if not isinstance(v, Vector):
            self_x_v = Vector()
            self_x_v.obj = self.obj * v
            return self_x_v
        else:
            result = expand(self.obj * v.obj)
            result = bilinear_product(result, Vector.dot)
            return result

    def __rmul__(self, s):
        s_x_self = Vector()
        s_x_self.obj = s * self.obj
        return s_x_self

    def __add__(self, v):
        self_p_v = Vector()
        self_p_v.obj = self.obj + v.obj
        return self_p_v

    def __add_ab__(self, v):
        self.obj += v.obj
        return

    def __sub__(self, v):
        self_m_v = Vector()
        self_m_v.obj = self.obj - v.obj
        return self_m_v

    def __sub_ab__(self, v):
        self.obj -= v.obj
        return

    def __pos__(self):
        return self

    def __neg__(self):
        n_self = copy.deepcopy(self)
        n_self.obj = -self.obj
        return n_self

    def applyfunc(self, fct):
        fct_self = Vector()
        fct_self.obj = fct(self.obj)
        return fct_self
