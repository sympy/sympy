from sympy.core.assumptions import StdFactKB
from sympy.matrices import Matrix, eye
from sympy.core import (Basic, Expr, Dummy, Symbol, symbols, AtomicExpr,
                        sympify, diff, Pow, Mul, Add, S)
from sympy.core.numbers import Zero
from sympy.solvers import solve
from sympy.simplify import simplify, trigsimp
from sympy.core.decorators import call_highest_priority, _sympifyit
from sympy.core.compatibility import as_int
from sympy.functions import sin, cos, tan, sqrt, atan

import copy

def base_scalars(names, coord_sys, ind=1):
    """
    Creates BaseScalars (coordinates) in a tuple for convenience.

    Parameters
    ==========

    names : string
        A space separated string of names for base_scalars.
    coord_sys : CoordSys
        The coordinate system in which the base scalars are to be creted.
    ind : integer
        Index to account for the position of the base scalar.

    Not to be used directly.

    Examples
    ========

    >>> from sympy.vector import CoordSysRect, base_scalars
    >>> c0 = CoordSysRect('c')
    >>> x, y, z = base_scalars('x y z', c0)
    >>> type(x)
    sympy.vector.vector.BaseScalar

    """
    base_scalar_list = []
    try:
        names = names.split(' ')
    except:
        pass

    for i, name in enumerate(names):
        base_scalar_list.append(BaseScalar(name, coord_sys, i + ind))
    return tuple(base_scalar_list)

def vectors(names, coord_sys, ind=1):
    """
    Creates BaseVectors in a tuple for convenience.

    Parameters
    ==========

    names : string
        A space separated string of names for base_vectors.
    coord_sys : CoordSys
        The coordinate system in which the base vectors are to be creted.
    ind : integer
        Index to account for the position of the base vectors.

    Not to be used directly.

    Examples
    ========

    >>> from sympy.vector import CoordSysRect, base_scalars
    >>> c0 = CoordSysRect('c')
    >>> ex, ey, ez = base_vectors('ex ey ez', c0)
    >>> type(ex)
    sympy.vector.vector.BaseVector

    """
    vector_list = []
    try:
        names = names.split(' ')
    except:
        pass

    for i, name in enumerate(names):
        vector_list.append(
            BaseVector(name, coord_sys=coord_sys, position=i+ind)
        )
    return tuple(vector_list)

def is_const_vect(vect):
    """
    Check if a vector field is a constant.

    Parameters
    ==========

    vect : Vector
        The vector to test.

    Examples
    ========

    >>> from sympy.vector import CoordSysRect, is_const_vect
    >>> c0 = CoordSysRect('c0')
    >>> v = 5 * c0.x + 7 * c0.y
    >>> is_const_vect(v)
    True

    """
    # First check if it is a vector at all
    if not vect.is_Vector:
        raise TypeError("argument of function needs to be a vector")
    vect = vect.expand()

    for atom in vect.atoms():
        if _has_base_scalar(atom):
            return False

        if (isinstance(atom, BaseVector) and
            (isinstance(atom.coord_sys, CoordSysCyl) or
             isinstance(atom.coord_sys, CoordSysSph))):
            return False

    return True

def _diff_scalar(scalar, s):
    """
    Differentiate a scalar.
    s : a BaseScalar
    """
    # Express scalar in rect coord_sys aligned with s.coord_sys
    c_rect = s.coord_sys._change_coord_sys(CoordSysRect)
    scalar = _express_scalar(scalar, c_rect)
    return diff(scalar, s)

def _dot_same(vect_a, vect_b):
    """
    The dot product of two vectors - both in same coordinate system.
    """
    a_comp = vect_a.components
    b_comp = vect_b.components
    length = len(a_comp)
    r = S.Zero
    for i in range(length):
        r = r + a_comp[i] * b_comp[i]
    return r

def dot(vect_a, vect_b, coord_sys=None):
    """
    Generalized dot product of two vectors.

    Parameters
    ==========

    vect_a : Vector
        First vector in the dot product.
    vect_b : Vector
        Second vector in the dot product.
    coord_sys : CoordSys
        A coordinate system to express the result in. If not provided,
        result will be expressed in all the coordinate systems present.

    Examples
    ========

    >>> from sympy.vector import CoordSysRect, CoordSysSph, dot
    >>> c0 = CoordSysRect('c0')
    >>> cs = CoordSysSph('cs')
    >>> v1 = c0.x * c0.e_x + c0.y * c0.e_y + c0.z * c0.e_z
    >>> v2 = cs.r * cs.e_r
    >>> dot(v1, v2)
    cs.r*sqrt(c0.x**2 + c0.y**2 + c0.z**2)

    """
    ret = S.Zero
    # Get two lists - each having vectors separated by coordinate system

    vect_a_coord_sys = _all_coordinate_systems(vect_a)
    vect_b_coord_sys = _all_coordinate_systems(vect_b)
    if not coord_sys:
        if(len(vect_a_coord_sys) == len(vect_b_coord_sys) == 1 and
           vect_a_coord_sys[0] == vect_b_coord_sys[0]):
            coord_sys = vect_a_coord_sys[0]

    if coord_sys:
        # Express each vector in coord_sys and call _dot_same
        vect_a = vect_a.express(coord_sys)
        vect_b = vect_b.express(coord_sys)
        return _dot_same(vect_a, vect_b)
    # We don't have a coord_sys so, proceed to return the output in all
    # the base scalars

    vect_a = vect_a.factor_vect()
    vect_b = vect_b.factor_vect()

    # Decompose the vectors
    vect_a_dec = []
    for arg in vect_a._all_args:
        t = (arg.scalar, arg.vector)
        vect_a_dec.append(t)

    vect_b_dec = []
    for arg in vect_b._all_args:
        t = (arg.scalar, arg.vector)
        vect_b_dec.append(t)

    for i in vect_a_dec:
        for j in vect_b_dec:
            ret = ret + i[0] * j[0] * dot(i[1], j[1], i[1].coord_sys)

    ret = ret.simplify()
    return ret

def grad(expr, coord_sys=None):
    """
    Calculate the gradient of a vector field.

    Parameters
    ==========

    expr : SymPy expression
        The SymPy expression whose gradient needs to be calculated.
        Needs to have at least one base scalar for non zero gradient.
    coord_sys : CoordSys
        A coordinate system to express the gradient in. If expr is in
        only one coordinate system, then, this parameter can be left
        out. Otherwise, raises ValueError.

    Examples
    ========

    >>> from sympy.vector import CoordSysRect, CoordSysSph, dot
    >>> c0 = CoordSysRect('c0')
    >>> cs = CoordSysSph('cs')
    >>> expr1 = c0.x**2 + exp(c0.z)*c0.y
    >>> grad(expr1)
    exp(c0.z) * c0.e_y + c0.y*exp(c0.z) * c0.e_z + 2*c0.x * c0.e_x
    >>> expr2 = 1/cs.r
    >>> grad(expr2)
    -1/cs.r**2 * cs.e_r

    # TODO : Add  more examples (once it is tested and works).
    """
    if not _has_base_scalar(expr):
        return ZeroVector

    coord_list = _coord_sys_scalar_list(expr)
    if not coord_sys:
        if len(coord_list) == 1:
            coord_sys = coord_list[0]
        else:
            raise TypeError("Coordinate system to express the vector in not known")


    # We are using the _express_scalar method to just express the given
    # sympy expression in the given coord_sys
    expr = _express_scalar(expr, coord_sys)

    # Now the expr has been converted completely to given coord_sys
    # We differentiate now
    ret = ZeroVector
    h_list = coord_sys.h_list
    base_scalars = coord_sys.base_scalars
    base_vectors = coord_sys.base_vectors

    for i, h in enumerate(h_list):
        l = (1/h) * diff(expr, base_scalars[i])
        ret = ret + l * base_vectors[i]
    return ret

def cross(vect_a, vect_b, coord_sys=None):
    """
    Generalized cross product of two vectors.

    Parameters
    ==========

    vect_a : Vector
        First vector in the cross product.
    vect_b : Vector
        Second vector in the cross product.
    coord_sys : CoordSys
        A coordinate system to express the result in. If not provided,
        result will be expressed in all the coordinate systems present.

    Examples
    ========

    >>> from sympy.vector import CoordSysRect, CoordSysSph, cross
    >>> c0 = CoordSysRect('c0')
    >>> v1 = c0.x * c0.e_y + 2*c0.y * c0.e_z
    >>> v2 = c0.z * c0.e_x
    >>> cross(v1, v2)
    -c0.x*c0.z * c0.e_z + 2*c0.y*c0.z * c0.e_y

    """
    ret = ZeroVector
    # Get two lists - each having vectors separated by coordinate system

    vect_a_coord_sys = _all_coordinate_systems(vect_a)
    vect_b_coord_sys = _all_coordinate_systems(vect_b)
    if not coord_sys:
        if (len(vect_a_coord_sys) == 1 and len(vect_b_coord_sys) == 1 and
           vect_a_coord_sys[0] == vect_b_coord_sys[0]):
            coord_sys = vect_a_coord_sys[0]

    if coord_sys:
        vect_a = vect_a.express(coord_sys)
        vect_b = vect_b.express(coord_sys)
        return _cross_same(vect_a, vect_b)

    ret = ZeroVector

    vect_a = vect_a.factor_vect()
    vect_b = vect_b.factor_vect()

    vect_a_dec = []
    for arg in vect_a._all_args:
        t = (arg.scalar, arg.vector)
        vect_a_dec.append(t)

    vect_b_dec = []
    for arg in vect_b._all_args:
        t = (arg.scalar, arg.vector)
        vect_b_dec.append(t)

    for i in vect_a_dec:
        for j in vect_b_dec:
            ret = ret + i[0] * j[0] * cross(i[1], j[1], i[1].coord_sys)

    ret = ret.expand().factor_vect()
    return ret

def _cross_same(vect_a, vect_b):
    """Cross product between two vectors of same coordinate system."""
    coord_a = _all_coordinate_systems(vect_a)
    coord_b = _all_coordinate_systems(vect_b)

    try:
        assert len(coord_a) == len(coord_b) == 1 and coord_a[0] == coord_b[0]
    except AssertionError:
        raise ValueError("Both vectors need to be defined in the same \
                          coordinate system")

    base_vectors = coord_a[0].base_vectors
    comp_a = vect_a.components
    comp_b = vect_b.components

    ret = ((comp_a[1] * comp_b[2] - comp_a[2] * comp_b[1]) * base_vectors[0] +
          (comp_a[2] * comp_b[0] - comp_a[0] * comp_b[2]) * base_vectors[1] +
          (comp_a[0] * comp_b[1] - comp_a[1] * comp_b[0]) * base_vectors[2])

    return ret

def div(vect, coord_sys=None):
    """
    Calculate the divergence of a vector field.

    Parameters
    ==========

    vect : Vector
        The vector field to whose divergence is to be calculated.
    coord_sys : CoordSys
        The coordinate system to express the results in.

    Examples
    ========

    >>> from sympy.vector import CoordSysRect, CoordSysSph, div
    >>> c0 = CoordSysRect('c0')
    >>> cs = CoordSysSph('cs')
    >>> v = c0.y**2*c0.z*c0.e_x + c0.y**3*c0.e_y + c0.x*c0.z*c0.e_z
    >>> div(v)
    c0.x + 3*c0.y**2
    >>> v = cs.r**2 * cs.e_r
    >>> div(v)
    4*cs.r

    """
    coord_list = _all_coordinate_systems(vect)
    if not coord_sys:
        if len(coord_list) == 1:
            coord_sys = coord_list[0]
        else:
            raise TypeError("Coordinate system to express the result in not known")

    if coord_sys.dim > 3:
        raise NotImplementedError

    vect = vect.express(coord_sys)
    base_scalars = coord_sys.base_scalars
    vect_comp = vect.components

    indices = [(1,2), (2, 0), (0, 1)] # Indices
    h_list = coord_sys.h_list
    ret = S.Zero

    for i, index in enumerate(indices):
        ret = ret + diff(vect_comp[i] * h_list[index[0]] * \
                         h_list[index[1]], base_scalars[i]).factor()
    ret = (1/(h_list[0] * h_list[1] * h_list[2])) * ret

    return ret.factor()


def curl(vect, coord_sys=None):
    """
    Calculate the curl of a vector field.

    Parameters
    ==========

    vect : Vector
        The vector field to whose curl is to be calculated.
    coord_sys : CoordSys
        The coordinate system to express the results in.

    Examples
    ========

    >>> from sympy.vector import CoordSysRect, CoordSysSph, div
    >>> c0 = CoordSysRect('c0')
    >>> cs = CoordSysSph('cs')
    >>> v = c0.y*c0.e_x - c0.x * c0.e_y
    >>> div(v)
    -2 * c0.e_z
    >>> v = (1/cs.r**2) * cs.e_r
    >>> res = div(v)
    >>> res
    0
    >>> type(res)
    sympy.vector.vector.ZeroVectorClass

    """

    coord_list = _all_coordinate_systems(vect)
    if not coord_sys:
        if len(coord_list) == 1:
            coord_sys = coord_list[0]
        else:
            raise TypeError("Coordinate system to express the vector in not known")

    if coord_sys.dim > 3:
        raise NotImplementedError
    vect = vect.express(coord_sys)
    bv = coord_sys.base_vectors
    bs = coord_sys.base_scalars
    comp = vect.components

    indices = [(2, 1), (0, 2), (1, 0)] # Coordinate Indices
    h = coord_sys.h_list
    ret = ZeroVector

    for i, index in enumerate(indices):
        diff_1 = h[index[0]] * comp[index[0]]
        diff_2 = h[index[1]] * comp[index[1]]
        l = ((diff(diff_1, bs[index[1]]) - diff(diff_2, bs[index[0]])) *
            (h[i] * bv[i]))
        # add if l in non zero
        if l:
            ret = ret + l.expand().factor_vect()
    return ret.factor_vect()

def laplacian(expr, coord_sys=None):
    """
    Calculate the laplacian of a vector field.

    Parameters
    ==========

    vect : Vector
        The vector field to whose laplacian is to be calculated.
    coord_sys : CoordSys
        The coordinate system to express the results in.

    # TODO : Add examples once testing is done.
    """

    return div(grad(expr, coord_sys))

def _has_base_scalar(expr):
    """
    Check whether a SymPy expression has a BaseScalar in it.

    Parameters
    ==========

    expr : a SymPy expression
        The expression to check.

    Examples
    ========

    >>> from sympy.vector.vector import _has_base_scalar
    >>> from sympy.vector import CoordSysRect
    >>> from sympy import Symbol
    >>> c0 = CoordSysRect('c0')
    >>> expr = 5 * (6 + c0.x) + 2/(3 + (2/c0.y))
    >>> _has_base_scalar(expr)
    True
    >>> x = Symbol('x')
    >>> expr = 2 * (x + (3/x))
    >>> _has_base_scalar(expr)
    False

    """
    try:
        expr = expr.expand()
    except AttributeError:
        # Some object which does not have an expand method
        return False

    for arg in expr.args:
        if isinstance(arg, BaseScalar):
            return True
        elif _has_base_scalar(arg):
            return True
    return False

def _all_base_scalars(scalar):
    """
    Returns all BaseScalars contained in a SymPy expression.

    Parameters
    ==========

    scalar : a SymPy expression.
         The scalar from which BaseScalars are to be extracted.

    Examples
    ========

    >>> from sympy.vector.vector import _all_base_scalars
    >>> from sympy.vector import CoordSysRect
    >>> from sympy import Symbol
    >>> x = Symbol('x')
    >>> c0 = CoordSysRect('c0')
    >>> expr = 5 * (6 + c0.x) + x/(3 + (2/c0.y))
    >>> _all_base_scalars(expr)
    [c0.y, c0.x]

    """
    if isinstance(scalar, BaseScalar):
        return [scalar]

    ret = []
    scalar = scalar.expand()
    for arg in scalar.args:
        if isinstance(arg, BaseScalar):
            ret.append(arg)
        else:
            l = _all_base_scalars(arg)
            if l:
                ret = ret + l
    return _remove_duplicates(ret)

def _separate_to_vectors(vect):
    # TODO: Used where?
    coord_dict = vect.separate()
    r = [v for v in coord_dict.values()]
    return r

def _all_coordinate_systems(vector):
    """
    Return all coordinate systems used to express a vector.

    Paramters
    =========

    vector : Vector
        The vector for which coordinate systems need to be extracted.

    Examples
    ========

    >>> from sympy.vector.vector import _all_cooordinate_systems
    >>> from sympy.vector import CoordSysRect, CoordSysSph
    >>> c0 = CoordSysRect('c0')
    >>> cs = CoordSysSph('c0')
    >>> v = c0.x * c0.e_y + cs.r * (c0.e_x + cs.e_r)
    >>> _all_base_scalars(expr)
    [cs, c0]

    """
    vector = vector.expand()
    coord_list = []
    # all_args is a separate method that return only the vector args
    # arg is either BaseVector or VectMul - becuase we have expanded
    # the vector
    for arg in vector._all_args:
        if isinstance(arg, BaseVector):
            coord_list.append(arg.coord_sys)
        elif isinstance(arg, VectMul):
            coord_list = coord_list + _coord_sys_scalar_list(arg.scalar)
            coord_list.append(arg.vector.coord_sys)
        else:
            # Shouldn't happen
            raise ValueError("Couldn't expand vector")

    return _remove_duplicates(coord_list)

def _all_coordinate_systems_scalar(scalar):
    """
    Return all coordinate syetms in a scalar.

    Parameters
    ==========

    scalar : a SymPy expression.
        The scalar from which coordinate systems have to be extracted.
    """
    coord_list = []
    base_scalars = _all_base_scalars(scalar)
    for b in base_scalars:
        if b.coord_sys not in coord_list:
            coord_list.append(b.coord_sys)
    return coord_list

def _coord_sys_scalar_list(scalar):
    """
    Return a list of coordinate systems contained within scalar.

    Parameters
    ==========

    expr : a SymPy expression
        The expression to extract coordinate systems from.

    See also
    ========

    _all_coordinate_systems : Find all coordinate systems used to
    construct a vector.
    """
    coord_list = []
    for atom in scalar.atoms():
        if isinstance(atom, BaseScalar):
            coord_list.append(atom.coord_sys)
    return _remove_duplicates(coord_list)


class BaseScalar(AtomicExpr):
    """
    BaseScalar instances are used to express coordinate variables for field.
    Not to be instantiated by the user directly.
    """
    is_comparable = False

    @property
    def _diff_wrt(self):
        """Allow derivatives wrt to BaseScalars"""
        return True

    def __init__(self, name, coord_sys, position, **assumptions):
        """
        Intialize a BaseScalar.
        Not to be used directly by the user.
        """
        if (position not in [1, 2, 3]
            and not isinstance(coord_sys, CoordSysRect)):
            raise ValueError("Position of scalar not specified. \
                            See `position` in docstring of `BaseScalar`")
        is_commutative = assumptions.get('commutative', True)
        assumptions['commutative'] = is_commutative

        if name is None:
            name = 'BaseSc_' + str(Dummy._count)
            Dummy._count += 1
        super(BaseScalar, self).__init__()
        self.coord_sys = coord_sys
        self.name = name
        self.position = S(position)
        self._assumptions = StdFactKB(assumptions)
        # self._args = (self.coord_sys, self.position)
        self._args = ()

    def __str__(self, printer=None):
        return self.coord_sys.name + "." + self.name

    __repr__ = __str__
    _sympystr = __str__

    def grad(self, coord_sys):
        __doc__ = grad.__doc__
        return grad(self, coord_sys)

    @staticmethod
    def coord_rel_rect(base_scalar, coord_sys):
        """
        Returns a relation between the given base_scalar base scalars
        in the given coordinate system which must be rectangular.

        Parameters
        ==========

        base_scalar : BaseScalar
            A base scalar whose relation is desired.
        coord_sys : CoordSys
            The coordinate system to which the relation of base_scalar
            is desired.
        # TODO : Where was this method used? If at all, add examples.
        """
        if base_scalar.coord_sys == coord_sys:
            return base_scalar

        if isinstance(base_scalar, CoordSysSph):
            x, y, z = coord_sys.base_scalars
            r, theta, phi = base_scalar.coord_sys.base_scalars

            subs_dict = { r: sqrt(x**2 + y**2 + z**2),
                          theta: atan(sqrt(x**2 + y**2)/z),
                          phi: atan(y/x) }
            return subs_dict[base_scalar]

        if isinstance(base_scalar, CoordSysCyl):
            x, y, z = coord_sys.base_scalars
            rho, phi, _z = base_scalar.coord_sys.base_scalars

            subs_dict = { rho: sqrt(x**2 + y**2),
                          phi: atan(y/x),
                          z: _z }

            return subs_dict[base_scalar]

    def _convert_to_rect(self, coord_sys, ind=1):
        """
        Takes a BaseScalar and converts it to combination of BaseScalars
        in coord_sys.

        Parameters
        ==========

        coord_sys: CoordSys
            The coordinate system to convert to.

        """
        # First check that both coordinate systems have same position and
        # orientation
        flag = False
        if self.coord_sys.position != coord_sys.position:
            flag = True
        if self.coord_sys._dcm_global != coord_sys._dcm_global:
            flag = True

        if flag:
            raise ValueError("coord_sys doesn't have same position and \
                              orientaion as coordinate system of self")

        # Check whether coord_sys.__class__ is same as self's coordinate
        # system
        if type(self.coord_sys) == type(coord_sys):
            return coord_sys.base_scalars[self.position - 1]
        rv = CoordSys._pos_to_rect(coord_sys.base_scalars,
                                   coord_conv[self.coord_sys.__class__])
        return rv[self.position - ind]

    def _convert_to_sph(self, coord_sys, ind=1):
        """
        Takes a BaseScalar and converts it to combination of BaseScalars
        in coord_sys.

        Parameters
        ==========

        coord_sys: CoordSys
            The coordinate system to convert to.

        """
        # First check that both coordinate systems have same position and
        # orientation
        flag = False
        if self.coord_sys.position != coord_sys.position:
            flag = True
        if self.coord_sys._dcm_global != coord_sys._dcm_global:
            flag = True

        if flag:
            raise ValueError("coord_sys doesn't have same position and \
                              orientaion as coordinate system of self")

        # Check whether coord_sys.__class__ is same as self's coordinate
        # system
        if type(self.coord_sys) == type(coord_sys):
            return coord_sys.base_scalars[self.position - 1]
        rv = CoordSys._pos_to_sph(coord_sys.base_scalars,
                                   coord_conv[self.coord_sys.__class__])
        return rv[self.position - ind]

    def _convert_to_cyl(self, coord_sys, ind=1):
        """
        Takes a BaseScalar and converts it to combination of BaseScalars
        in coord_sys.

        Parameters
        ==========

        coord_sys: CoordSys
            The coordinate system to convert to.

        """
        # First check that both coordinate systems have same position and
        # orientation
        flag = False
        if self.coord_sys.position != coord_sys.position:
            flag = True
        if self.coord_sys._dcm_global != coord_sys._dcm_global:
            flag = True

        if flag:
            raise ValueError("coord_sys doesn't have same position and \
                              orientaion as coordinate system of self")

        # Check whether coord_sys.__class__ is same as self's coordinate
        # system
        if type(self.coord_sys) == type(coord_sys):
            return coord_sys.base_scalars[self.position - 1]
        rv = CoordSys._pos_to_cyl(coord_sys.base_scalars,
                                   coord_conv[self.coord_sys.__class__])
        return rv[self.position - ind]

    def _hashable_content(self):
        return (repr(self),) + tuple(sorted(self.assumptions0.items()))

class CoordSys(Basic):
    """
    Superclass for all CoordSys<type> classes. Not to be intialized
    directly.
    """
    def __new__(cls, name=None, dim=None, position=None, pos_flag=None,
                orient_type=None, orient_amount=None, rot_order=None,
                parent=None, basis_vectors=None, coordinates=None):
        if orient_amount:
            orient_amount = tuple(orient_amount)
        return Basic.__new__(cls, name, dim, position, pos_flag,
                             orient_type, orient_amount, rot_order,
                             parent, basis_vectors, coordinates)

    def __init__(
            self, name=None, dim=3, position=None, pos_flag=True,
            orient_type=None, orient_amount=None, rot_order=None,
            parent=None, basis_vectors=None, coordinates=None
            ):
        """
        Special __init___ method for CoordSys class.

        Parameters
        ==========

        name : string
            The name of the coordinate system.
        dim : integer
            The dimension of the coordinate system. Can take values != 3
            for only rectangular coordinates.
        position : Vector
            The vector position of the coordinate system.
        pos_flag : Boolean
            A bool to determine whether the position is to be interpreted
            as being wrt the parent, if one exists. Defaults to True.
            If parent exists, the position is taken as if it was expressed
            wrt the parents by default.
        orient_type : string
            The type of orientation to be used.
            Values:
            ------
            'Axis' : Rotate the coordinate system about an axis.
            'Body' : Rotate the coordinate system about the axes,
                     one at a time in a specific order.
             #TODO : Implement two remaining type of orientations.
        orient_amount : tuple
             The tuple that tells the amount by which the rotation is to
             be performed.
             If orient_type is 'Axis':
                 tuple if of length 2.
                 First element is the angle in radians.
                 Second element (Vector) is the axis of rotation.
             If orient_type is 'Body':
                 tuple of length 3.
                 All elements are angles in radians that will be applied
                 to the rotation in the specific order (see rot_order below).
        rot_order : string (length 3)
             The string to represent rotation order in case of 'Body' type
             of orientation. Not to be provided for 'Axis' orientation.
        parent : CoordSys
             The parent coordinate system of this coordinate system.
        basis_vectors : string
             If the default naming of the base vectors is to be changed, this
             string (space separated names of vectors) should be provided.
        coordinates : string
             If the default naming of the base scalars is to be changed, this
             string (space separated names of scalars) should be provided.

        Not to be initialized directly.
        """
        if name is None:
            name = 'CoordSys_' + str(Dummy._count)
            Dummy._count += 1

        self.name = name
        self.dim = int(dim)
        # Moved setting of parent to here
        self.parent = parent

        self._coord_names = coordinates
        self._basis_names = basis_vectors
        self._coordinates = base_scalars(self._coord_names, self)
        for i, v in enumerate(self._coord_names):
            self.__setattr__(v, self._coordinates[i])
        self._basis = vectors(self._basis_names, self)
        for i, v in enumerate(self._basis_names):
            self.__setattr__(v, self._basis[i])

        if position:
            # check whether position is a vector
            if not position.is_Vector:
                raise TypeError("vector expected for position")
            if not is_const_vect(position):
                raise ValueError("Position vector needs to be a constant")

            # The position vector needs to be constant, hence, it cannot be in
            # any other coordinate system other than rectangular.
            if parent and pos_flag:
                self.position = _vect_add_const(parent.position, position,
                                                parent)
            elif parent and not pos_flag:
                # Take the position to be expressed as global coordinates
                self.position = position
            else:
                self.position = position
        else:
            self.position = ZeroVector

        # _dcm_parent stores matrix to go from self -> parent
        # _dcm_global stores matrix to go from global -> self
        self._dcm_global = eye(3)
        self._dcm_parent = eye(3)

        # TODO : Remove the string interface; use functions
        self.orient_type = orient_type
        self.orient_amount = orient_amount
        self.rot_order = rot_order
        if orient_type:
            self._check_orient_raise(orient_type, orient_amount, rot_order)
            if parent:
                self._dcm_parent = self._dcm_parent_method(orient_type,
                                                         orient_amount,
                                                         rot_order)

            self._dcm_global = self._dcm_global_method(orient_type,
                                                       orient_amount,
                                                       rot_order)

    def __str__(self, printer=None):
        return self.name
    __repr__ = __str__
    _sympystr = __str__
    _sympyrepr = __str__

    """
    def _hashable_content(self):
        return (repr(self),) + tuple(sorted(self.assumptions0.items()))
    """

    def _dcm_parent_method(self, orient_type, amounts, rot_order=''):
        """Intialize the self._dcm_parent attribute."""
        orient_type = orient_type.capitalize()
        if orient_type == 'Axis':
            if rot_order:
                raise TypeError('Axis orientation takes no rotation order')
            if not (isinstance(amounts, (list, tuple)) & (len(amounts) == 2)):
                raise TypeError('Amounts are a list or tuple of length 2')
            theta = amounts[0]
            axis = amounts[1]
            axis = axis.express(self.parent).normalize().as_mat().T
            parent_orient = ((eye(3) - axis * axis.T) * cos(theta) +
                    Matrix([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]],
                        [-axis[1], axis[0], 0]]) * sin(theta) + axis * axis.T)
            return parent_orient

        elif orient_type == 'Body':
            if not (len(amounts) == 3 & len(rot_order) == 3):
                raise TypeError('Body orientation takes 3 values & 3 orders')
            a1 = int(rot_order[0])
            a2 = int(rot_order[1])
            a3 = int(rot_order[2])
            parent_orient = (self._rot(a1, amounts[0]) *
                            self._rot(a2, amounts[1]) *
                            self._rot(a3, amounts[2]))
            return parent_orient


    def _dcm_global_method(self, orient_type, amounts, rot_order=''):
        """Initialize the self._dcm_global attribute."""
        orient_type = orient_type.capitalize()
        if self.parent:
            # Parent given therefore the given angle is wrt parent.
            # Go from self -> parent -> global
            parent = self.parent
            r = parent.dcm().T * parent.dcm(self)
            return r.T
        else:
            # Parent not set. So, directly initialize the dcm wrt global
            # Using sympy.physics.mechanics logic here.
            if orient_type == 'Axis':
                if rot_order:
                    raise ValueError('Axis orientation \
                                      takes no rotation order')
                if not (isinstance(amounts, (list, tuple))
                        & (len(amounts) == 2)):
                    raise TypeError('Amounts are a list or tuple of length 2')
                theta = amounts[0]
                axis = amounts[1]
                # TODO : No such method. Do it via express.
                c_rect = CoordSysRect # Dummy coord_sys for expressing axis
                axis = axis.express(c_rect)
                axis = axis.normalize().as_mat()
                global_orient = ((eye(3) - axis * axis.T) * cos(theta) +
                        Matrix([[0, -axis[2], axis[1]],
                                [axis[2], 0, -axis[0]],
                                [-axis[1], axis[0], 0]]) * \
                                 sin(theta) + axis * axis.T)
                return global_orient
            if orient_type == 'Body':
                if not (len(amounts) == 3 & len(rot_order) == 3):
                    raise TypeError('Body orientation takes \
                                     3 values & 3 orders')
                a1 = int(rot_order[0])
                a2 = int(rot_order[1])
                a3 = int(rot_order[2])
                global_orient = (self._rot(a1, amounts[0]) *
                                self._rot(a2, amounts[1]) *
                                self._rot(a3, amounts[2]))
                global_orient.applyfunc(lambda i: trigsimp(i, method='fu'))
                return global_orient

    def _rot(self, axis, angle):
        """DCM about one of the 3 axes."""
        if axis == 1:
            return Matrix([[1, 0, 0],
                [0, cos(angle), -sin(angle)],
                [0, sin(angle), cos(angle)]])
        elif axis == 2:
            return Matrix([[cos(angle), 0, sin(angle)],
                [0, 1, 0],
                [-sin(angle), 0, cos(angle)]])
        elif axis == 3:
            return Matrix([[cos(angle), -sin(angle), 0],
                [sin(angle), cos(angle), 0],
                [0, 0, 1]])

    def _check_orient_raise(self, orient_type, orient_amount, rot_order):
        orient_type = orient_type.capitalize()
        if orient_type == 'Axis':
            if len(orient_amount) != 2:
                raise ValueError("orient_amount should be a \
                                list of lenght 2")
            if rot_order:
                raise ValueError("orient_axis works only with 'Body' \
                                type orientaion")
            if not is_const_vect(orient_amount[1]):
                raise ValueError("Orient axis needs to be constant vector")

        elif orient_type == 'Body':
            if len(orient_amount) != 3:
                raise ValueError("orient_amount must be a list of \
                                length 3")

            possible_orders = (
                        '123', '132', '213', '312', '312', '321',
                        '121', '131', '212', '313', '313', '323'
                        )
            if rot_order not in possible_orders:
                raise ValueError("Wrong rotation order")

        else:
            raise ValueError("orient_type not known")

    @property
    def base_scalars(self):
        """Return a list of base scalars."""
        return self._coordinates

    @property
    def base_vectors(self):
        """Return a list of base vectors."""
        return self._basis

    def dcm(self, other='global'):
        """
        Returns the direction cosine matrix for converting from `other`
        to `self`, i.e.,
        self.xyz = self.dcm(other) * other.xyz

        Parameters
        ==========

        other : CoordSys
            The coordinate system from which we want the DCM to self.

        Examples
        ========

        >>> from sympy.vector import CoordSysRect
        >>> c0 = CoordSysRect('c0')
        >>> c1 = CoordSysRect
        >>> cr = CoordSysRect('cr', orient_type='Axis', \
        orient_amount=(q1, c0.e_z), parent=c0)
        >>> cr.dcm(c0)
        Matrix([
        [ cos(q1), sin(q1), 0],
        [-sin(q1), cos(q1), 0],
        [       0,       0, 1]])

        """
        if other == 'global':
            return self._dcm_global

        # First check self and other are parent/child frames
        if self.parent == other:
            return self._dcm_parent.T
        elif other.parent == self:
            return other._dcm_parent
        else:
            # Now we can either traverse the frames to find the nearest
            # common frame. This will require multiplying many matrices.
            # So, we go via the global frame
            # TODO : Implement path A->B->C->D rather than A->global->D
            r = other.dcm() * self.dcm().T
            r.applyfunc(lambda i: trigsimp(i, method='fu'))
            return r

    def _dcm_path(self, other=None):
        # self -> child, other -> parent
        # First handle other
        other_coord_list = [other]

        other_parent = other
        while 1:
            other_parent = other.parent
            if other_parent:
                other_coord_list.append(p)
            else:
                break

        common_ancestor = None
        self_coord_list = [self]
        self_parent = self
        while 1:
            self_parent = self.parent
            if self.parent:
                if self_parent in other_coord_list:
                    common_ancestor = self_parent
                else:
                    self_coord_list.append(self_parent)
            else:
                break

        # Construct a list going from other to self
        new_list = []
        for c in other_coord_list:
            new_list.append(c)
        for c in reversed(self_coord_list):
            new_list.append(c)

        return new_list, common_ancestor

    def orientnew(self, name=None, orient_type=None, orient_amount=None,
                 rot_order=None):
        """
        Create a coordinate system with a new orientation from an exisitng
        coordinate system.

        Parameters
        ==========

        See the parameters in the docstring of CoordSys.__init__

        Examples
        ========

        >>> from sympy import Symbol
        >>> from sympy.vector import CoordSysRect
        >>> q0 = Symbol('q0')
        >>> c0 = CoordSysRect('c0')
        >>> c1 = c0.orientnew('c1', 'Axis', [q0, c0.e_z])
        >>> c1.dcm(c0)
        Matrix([
        [ cos(q0), sin(q0), 0],
        [-sin(q0), cos(q0), 0],
        [       0,       0, 1]])

        """
        if name == None:
            name = 'CoordSys_' + str(Dummy._count)
            Dummy._count += 1

        newframe = self.__class__(
            name=name, dim=self.dim, position=self.position, pos_flag=False,
            orient_type=orient_type, orient_amount=orient_amount,
            rot_order=rot_order, parent=self, basis_vectors=self._basis_names,
            coordinates=self._coord_names)

        return newframe

    def posnew(self, position, name=None):

        """
        Create a coordinate system with a new position from an exisitng
        coordinate system.

        Parameters
        ==========

        position : Vector
            The position vector of the new coordinate system.
        name : string
            Name of the coordinate system.

        Examples
        ========

        >>> from sympy import Symbol
        >>> from sympy.vector import CoordSysRect
        >>> c0 = CoordSysRect('c0')
        >>> v = 2 * c0.e_x + 3 * c0.e_y
        >>> c1 = posnew(v)

        """
        if not position.is_Vector:
            raise TypeError("vector expected for position")
        if not position.is_constant:
            raise TypeError("Position variables need to be constant.")

        if not name:
            name = 'CoordSys_' + str(Dummy._count)
            Dummy._count += 1

        newframe = self.__class__(
            name=name, dim=self.dim, position=position,
            orient_type=self.orient_type, orient_amount=self.orient_amount,
            rot_order=self.rot_order, parent=self,
            basis_vectors=self._basis_names, coordinates=self._coord_names)

        return newframe

    def _change_coord_sys(self, coord_sys_class, name=None, flag=False):
        """
        Return a coordinate system after changing the type of coordinate
        system used, keeping other things same.

        Parameters
        ==========

        coord_sys_class : class object
             The class in which the resultant coordinate system is
             desired.

        name : string
             The name of the returned coordinate system.

        flag : bool
             Value that enforces the use of same names of bases as the
             original coordinate system, if true.

        """
        if not name:
            name = "coord_sys_" + str(Dummy._count)
            Dummy._count += 1
        CoordSysClass = coord_sys_class
        if flag:
            ret_coord_sys = CoordSysClass(name=name,
                                    position=self.position, pos_flag=False,
                                    orient_type=self.orient_type,
                                    orient_amount=self.orient_amount,
                                    rot_order=self.rot_order,
                                    parent=self.parent,
                                    basis_vectors=self._basis_names,
                                    coordinates=self._coord_names)
        else:
            ret_coord_sys = CoordSysClass(name=name,
                                    position=self.position, pos_flag=False,
                                    orient_type=self.orient_type,
                                    orient_amount=self.orient_amount,
                                    rot_order=self.rot_order,
                                    parent=self.parent)

        return ret_coord_sys

    @staticmethod
    def _pos_to_rect(coord, conv_from):
        if conv_from == 'cyl':
            rho = sqrt(coord[0]**2 + coord[1]**2)
            phi = coord[2]
            z = coord[2]
            return (rho, phi, z)

        if conv_from == 'sph':
            r = sqrt(coord[0]**2 + coord[1]**2 + coord[2]**2)
            theta = atan(sqrt(coord[0]**2 + coord[1]**2)/coord[2])
            phi = atan(coord[1]/coord[0])
            return (r, theta, phi)


    @staticmethod
    def _pos_to_cyl(coord, conv_from):
        if conv_from == 'rect':
            x = coord[0] * cos(coord[1])
            y = coord[0] * sin(coord[1])
            z = coord[2]
            return (x, y, z)


        if conv_from == 'sph':
            r = sqrt(coord[0]**2 + coord[2]**2)
            theta = atan(coord[0]/coord[2])
            phi = coord[1]
            return (r, theta, phi)


    @staticmethod
    def _pos_to_sph(coord, conv_from):
        if conv_from == 'rect':
            x = coord[0] * sin(coord[1]) * cos(coord[2])
            y = coord[0] * sin(coord[1]) * sin(coord[2])
            z = coord[0] * cos(coord[1])
            return (x, y, z)


        if conv_from == 'cyl':
            rho = coord[0] * sin(coord[1])
            phi = coord[2]
            z = coord[0] * cos(coord[1])
            return (rho, phi, z)


    @staticmethod
    def _orient_along(vect, coord_sys):
        # TODO : Used where?
        if vect.coord_sys.dim > 3 or coord_sys.dim > 3:
            raise NotImplementedError
        Ax0, Ay0, Az0 = vect.components
        mat = coord_sys.dcm(vect.coord_sys)
        mat_components = mat * Matrix([[Ax0], [Ay0], [Az0]])

        ex0, ey0, ez0 = vect.base_vectors
        mat_base_vectors = mat * Matrix([[ex0], [ey0], [ez0]])

        base_vectors_list = mat_base_vectors._mat

        ret = []
        for i, comp in enumerate(mat_components._mat):
            ret.append(VectMul(comp, base_vectors_list[i]))
        return VectAdd(*ret)

    @staticmethod
    def _convert_base_vect_rect(vect, coord_sys):
        """
        Returns base vectors in rectangular coordinates converted
        to another set of rectangular coordinates while changing
        the orientation.

        Parameters
        ==========

        vect : Vector (in rectangular coordinates)
            The vector to be converted.
        coord_sys: A CoordSys object.
            The coordinate system to express the result in.
            Needs to be rectangular.

        """
        if vect.coord_sys == coord_sys:
            return vect

        Ax0, Ay0, Az0 = vect.components

        mat = vect.coord_sys.dcm(coord_sys)
        mat_components = mat * Matrix([[Ax0], [Ay0], [Az0]])

        # construct the vector to return the value
        ret = []
        base_vectors = coord_sys.base_vectors

        for i, comp in enumerate(mat_components._mat):
            ret.append(VectMul(comp, base_vectors[i]))
        vect = VectAdd(*ret)
        return vect

    @staticmethod
    def _convert_base_sclr_rect(sclr, coord_sys):
        """
        Returns base scalars in rectangular coordinates converted
        to another set of rectangular coordinates while changing
        the orientation.

        Parameters
        ==========

        sclr : a BaseScalar (in rectangular coordinates)
            The scalar to be converted.
        coord_sys: A CoordSys object.
            The coordinate system to express the result in.
            Needs to be rectangular.

        """
        if sclr.coord_sys == coord_sys:
            return sclr

        # Checking for constant terms (for example, in a constant vector)
        if not isinstance(sclr, BaseScalar):
            return sclr

        # Let the sclr.coord_sys == B, coord_sys == A
        # B.xyz = pos_A_wrt_B + B.dcm(A) * A.xyz

        # A.x, A.y, A.z
        x0, y0, z0 = coord_sys.base_scalars

        mat = sclr.coord_sys.dcm(coord_sys)
        mat_components = mat * Matrix([[x0], [y0], [z0]])
        mat_components = mat_components._mat

        # Now, position vectors are just constant vectors. So, we can
        # call express on them without causing a recursion error.
        csA = coord_sys._change_coord_sys(CoordSysRect, 'csA')
        csB = sclr.coord_sys._change_coord_sys(CoordSysRect, 'csB')
        rel_pos_vect = coord_sys.position.express(csA) - \
                       sclr.coord_sys.position.express(csB)
        if not rel_pos_vect:
            rel_pos_vect = ZeroVector

        rel_pos_vect_comp = rel_pos_vect.expand().factor_vect().components
        res = []

        for i, comp in enumerate(rel_pos_vect_comp):
            res.append(S(mat_components[i] + comp))

        return res[sclr.position - 1]


class CoordSysRect(CoordSys):
    """The class to represent rectangular coordinate system."""
    def __init__(self, *args, **kwargs):
        """
        Intialize a rectangular coordinate system.

        See the docstring of CoordSys class for more.
        """
        if 'coordinates' not in kwargs.keys():
            try:
                kwargs['coordinates'] = \
                        ['x' + str(i + 1) for i in range(kwargs['dim'])]
            except:
                kwargs['coordinates'] = ('x', 'y', 'z')
        if 'basis_vectors' not in kwargs.keys():
            try:
                kwargs['basis_vectors'] = \
                    ['e_' + str(i + 1) for i in range(kwargs['dim'])]
            except:
                kwargs['basis_vectors'] = ('e_x', 'e_y', 'e_z')

        # Convert to tuple
        kwargs['coordinates'] = tuple(kwargs['coordinates'])
        kwargs['basis_vectors'] = tuple(kwargs['basis_vectors'])

        super(CoordSysRect, self).__init__(*args, **kwargs)
        self.h_list = tuple([S.One] * self.dim)


    @staticmethod
    def _convert_to_sph(vector, coord_sys):
        """
        Return a vector in this coordinate system converted to
        spherical coordinates.

        Parameters
        ==========

        vector : Vector
            The vector to convert.
        coord_sys : CoordSys
            The coordinate system to convert to.

        """
        # TODO : Implement a coord_sys check method to check same orientation
        # and position of vector.coord_sys and coord_sys
        Ax, Ay, Az = vector.components
        x, y, z = vector.coord_sys.base_scalars
        r, theta, phi = coord_sys.base_scalars
        subs_dict = { x: r * sin(theta) * cos(phi),
                      y: r * sin(theta) * sin(phi),
                      z: r * cos(theta) }
        Ax = Ax.subs(subs_dict)
        Ay = Ay.subs(subs_dict)
        Az = Az.subs(subs_dict)
        mat =  Matrix([
                       [sin(theta)*cos(phi), sin(theta)*sin(phi),  cos(theta)],
                       [cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)],
                       [          -sin(phi),            cos(phi),           0]
                      ])
        r = mat*Matrix([ [Ax], [Ay], [Az] ])
        res = []
        for i, vect in enumerate(r._mat):
            res.append(VectMul(vect, coord_sys.base_vectors[i]))
        return VectAdd(*res)

    @staticmethod
    def _convert_to_cyl(vector, coord_sys):
        """
        Return a vector in this coordinate system converted to
        cylindrical coordinates.

        Parameters
        ==========

        vector : Vector
            The vector to convert.
        coord_sys : CoordSys
            The coordinate system to convert to.

        """
        Ax, Ay, Az = vector.components
        x, y, z = vector.coord_sys.base_scalars
        rho, phi, _z = coord_sys.base_scalars
        subs_dict = { x: rho * cos(phi),
                      y: rho * sin(phi),
                      z: _z }
        Ax = Ax.subs(subs_dict)
        Ay = Ay.subs(subs_dict)
        Az = Az.subs(subs_dict)
        mat =  Matrix([ [ cos(phi), sin(phi), 0],
                        [-sin(phi), cos(phi), 0],
                        [        0,        0, 1] ])
        r = mat*Matrix([ [Ax], [Ay], [Az] ])
        res = []
        for i, vect in enumerate(r._mat):
            res.append(VectMul(vect, coord_sys.base_vectors[i]))
        return VectAdd(*res)

    def __str__(self):
        return self.name

    __repr__ = __str__
    _sympystr = __str__
    _sympyrepr = __str__


class CoordSysSph(CoordSys):
    """The class to represent spherical polar coordinate system."""
    def __init__(self, *args, **kwargs):
        """
        Initialize a spherical coordinate system.

        See the docstring of CoordSys for more.
        """
        if 'coordinates' not in kwargs.keys():
            kwargs['coordinates'] = ('r', 'theta', 'phi')
        if 'basis_vectors' not in kwargs.keys():
            kwargs['basis_vectors'] = ('e_r', 'e_theta', 'e_phi')
        super(CoordSysSph, self).__init__(*args, **kwargs)

        self.dim = S(3)

        self.h_list = (S.One,
                       self.base_scalars[0],
                       self.base_scalars[0]*sin(self.base_scalars[1]))

    @staticmethod
    def _convert_to_rect(vector, coord_sys):
        """
        Return a vector in this coordinate system converted to
        rectangular coordinates.

        Parameters
        ==========

        vector : Vector
            The vector to convert.
        coord_sys : CoordSys
            The coordinate system to convert to.

        """
        Ar, At, Ap = vector.components
        r, theta, phi = vector.coord_sys.base_scalars
        x, y, z = coord_sys.base_scalars
        subs_dict = { r: sqrt(x**2 + y**2 + z**2),
                      theta: atan(sqrt(x**2 + y**2)/z),
                      phi: atan(y/x) }
        Ar = Ar.subs(subs_dict)
        At = At.subs(subs_dict)
        Ap = Ap.subs(subs_dict)

        mat = Matrix([
            [x/sqrt(x**2 + y**2 + z**2),
             x*z/(sqrt(x**2 + y**2)*sqrt(x**2 + y**2 + z**2)),
             -y/sqrt(x**2 + y**2)],
            [y/sqrt(x**2 + y**2 + z**2),
             y*z/(sqrt(x**2 + y**2)*sqrt(x**2 + y**2 + z**2)),
             x/sqrt(x**2 + y**2)],
            [z/sqrt(x**2 + y**2 + z**2),
             -sqrt(x**2 + y**2)/sqrt(x**2 + y**2 + z**2), 0]])

        r = mat*Matrix([ [Ar], [At], [Ap] ])
        res = []
        for i, vect in enumerate(r._mat):
            res.append(VectMul(vect, coord_sys.base_vectors[i]))
        return VectAdd(*res)

    @staticmethod
    def _convert_to_cyl(vector, coord_sys):
        """
        Return a vector in this coordinate system converted to
        cylindrical coordinates.

        Parameters
        ==========

        vector : Vector
            The vector to convert.
        coord_sys : CoordSys
            The coordinate system to convert to.

        """
        Ar, At, Ap = vector.components
        r, theta, phi = vector.coord_sys.base_scalars
        rho, _phi, z = coord_sys.base_scalars
        subs_dict = { r: sqrt(rho**2 + z**2),
                      theta: atan(rho/z),
                      phi: _phi }

        Ar = Ar.subs(subs_dict)
        At = At.subs(subs_dict)
        Ap = Ap.subs(subs_dict)
        mat = Matrix([
                        [rho/(z*sqrt(rho**2/z**2 + 1)),
                        1/sqrt(rho**2/z**2 + 1), 0 ],
                        [0, 0, 1],
                        [ 1/sqrt(rho**2/z**2 + 1),
                        -rho/(z*sqrt(rho**2/z**2 + 1)), 0]
                    ])
        r = mat*Matrix([ [Ar], [At], [Ap] ])
        res = []
        for i, vect in enumerate(r._mat):
            res.append(VectMul(vect, coord_sys.base_vectors[i]))
        return VectAdd(*res)

    def __str__(self):
        return self.name

    __repr__ = __str__
    _sympystr = __str__
    _sympyrepr = __str__


class CoordSysCyl(CoordSys):
    """The class to represent cylindrical coordinate system."""
    def __init__(self, *args, **kwargs):
        """
        Initialize a cylindrical coordinate system.

        See the docstring of CoordSys for more.
        """
        if 'coordinates' not in kwargs.keys():
            kwargs['coordinates'] = ('rho', 'phi', 'z')
        if 'basis_vectors' not in kwargs.keys():
            kwargs['basis_vectors'] = ('e_rho', 'e_phi', 'e_z')
        super(CoordSysCyl, self).__init__(*args, **kwargs)

        self.dim = S(3)

        self.h_list = (S.One,
                       self.base_scalars[0],
                       S.One)

    @staticmethod
    def _convert_to_rect(vector, coord_sys):
        """
        Return a vector in this coordinate system converted to
        rectangular coordinates.

        Parameters
        ==========

        vector : Vector
            The vector to convert.
        coord_sys : CoordSys
            The coordinate system to convert to.

        """
        Ar, Ap, Az = vector.components
        rho, phi, z = vector.base_scalars
        x, y, _z = coord_sys.base_scalars
        subs_dict = { rho: sqrt(x**2 + y**2),
                      phi: atan(y/x),
                      z: _z }

        Ar = Ar.subs(subs_dict)
        Ap = Ap.subs(subs_dict)
        Az = Az.subs(subs_dict)
        mat = Matrix([
            [    1/sqrt(1 + y**2/x**2), -y/(x*sqrt(1 + y**2/x**2)), 0],
            [y/(x*sqrt(1 + y**2/x**2)),      1/sqrt(1 + y**2/x**2), 0],
            [                        0,                          0, 1]
        ])

        r = mat*Matrix([ [Ar], [Ap], [Az] ])
        res = []
        for i, vect in enumerate(r._mat):
            res.append(VectMul(vect, coord_sys.base_vectors[i]))
        return VectAdd(*res)

    @staticmethod
    def _convert_to_sph(vector, coord_sys):
        """
        Return a vector in this coordinate system converted to
        spherical coordinates.

        Parameters
        ==========

        vector : Vector
            The vector to convert.
        coord_sys : CoordSys
            The coordinate system to convert to.

        """
        Ar, Ap, Az = vector.components
        rho, phi, z = vector.coord_sys.base_scalars
        r, theta, _phi = coord_sys.base_scalars
        subs_dict = { rho: r*sin(theta),
                      phi: _phi,
                      z: r*cos(theta) }

        Ar = Ar.subs(subs_dict)
        Ap = Ap.subs(subs_dict)
        Az = Az.subs(subs_dict)
        mat = Matrix([
            [sin(theta), 0,  cos(theta)],
            [cos(theta), 0, -sin(theta)],
            [         0, 1,           0]
        ])

        r = mat*Matrix([ [Ar], [Ap], [Az] ])
        res = []
        for i, vect in enumerate(r._mat):
            res.append(VectMul(vect, coord_sys.base_vectors[i]))
        return VectAdd(*res)

    def __str__(self):
        return self.name

    __repr__ = __str__
    _sympystr = __str__
    _sympyrepr = __str__


class Vector(AtomicExpr):
    """Base class of all vector classes."""

    is_Vector = True
    _op_priority = 11.0

    def dot(self, other, coord_sys=None):
        __doc__ = dot.__doc__
        return dot(self, other, coord_sys)

    def cross(self, other, coord_sys=None):
        __doc__ = cross.__doc__
        return dot(self, other, coord_sys)

    def div(self, coord_sys=None):
        __doc__ = div.__doc__
        return div(self, coord_sys)

    def curl(self, coord_sys=None):
        __doc__ = curl.__doc__
        return curl(self, coord_sys)

    def __neg__(self):
        return VectMul(S.NegativeOne, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__radd__')
    def __add__(self, other):
        return _vect_add(self, other)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__add__')
    def __radd__(self, other):
        return _vect_add(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return _vect_add(self, -other)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return _vect_add(other, -self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return _vect_mul(self, other)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return _vect_mul(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        return _vect_div(self, other)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        raise TypeError("Cannot divide by vector")

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def __getitem__(self, i):
        if self.coord_sys.dim > 3:
            raise NotImplementedError
        comp = self.components
        if i > 0 and i <= 3:
            return comp[i]
        else:
            raise IndexError

    def simplify(self):
        """Simplify a vector."""
        ret = ZeroVector
        vect = self.expand().factor_vect()
        for arg in vect._all_args:
            scalar = arg.scalar
            vector = arg.vector
            scalar = scalar.simplify()
            ret = ret + (scalar * vector)
        return ret

    def normalize(self):
        """Normalize the vector. Returns vector divided by its magnitude."""
        return self/sqrt(dot(self, self))

    def as_mat(self):
        """
        Returns the components of vector in matrix form.
        Returns a column matrix.
        The vector should be entirely in one coordinate system.
        """
        coord_list = _all_coordinate_systems(self)
        if len(coord_list) > 1:
            raise ValueError("Vector must be in a single coordinate system")
        return Matrix([self.components])

class BaseVector(Vector, Symbol):
    """Class to represent basis vectors."""
    # __new__ required because it the __new__ of super class takes coord_sys,
    # position as a part of assuptions
    def __new__(cls, name, coord_sys, position, **assumptions):
        # Call to Vector's __new__ method so that BaseVectors with same names
        # don't get cached
        obj = Vector.__new__(BaseVector, name, **assumptions)
        return obj

    def __init__(self, name, coord_sys, position):
        """
        Intialize a BaseVector.
        Not to be initialized directly.

        Parameters
        ==========

        name : string
            The name of the basis vector.

        coord_sys : CoordSys
            The coordinate system in which the basis vector has been
            defined in.

        position : integer
            Integer to represent which coordinate does the vector represent.
            Possible values : 1, 2 or 3

        """
        self.name = name
        if not isinstance(coord_sys, CoordSys):
            raise TypeError("coord_sys must be isinstance(CoordSys)")
        self.coord_sys = coord_sys
        self.position = S(position)
        if position > coord_sys.dim:
            raise ValueError("Vector outside coordinate system dimensionality")

        is_commutative = True
        assumptions = {}
        assumptions['commutative'] = is_commutative
        self._assumptions = StdFactKB(assumptions)
        #self._args = (self.coord_sys, self.position)
        self._args = ()

    def separate(self):
        # We just have a Vector - just return it
        coord_sys_dict = {self.coord_sys: [[], [self]]}
        return coord_sys_dict

    def express(self, coord_sys):
        __doc__ = express.__doc__
        return express(self, coord_sys)

    def expand(self):
        return self

    def factor_vect(self):
        return self

    @property
    def _all_args(self):
        return [self]

    @property
    def scalar(self):
        return S.One

    @property
    def vector(self):
        return self

    @property
    def components(self):
        """
        Returns a list of the components of a Vector.

        In this case, it is a basis vector, so the list is all zeros except for
        1 entry of unity; e.g. for the first basis vector in a rect coordinate
        system, the returned list is [1, 0, 0]
        """
        # Since it is a base vector, so return a list
        # of len == dim(coord_sys) with the pos element
        # as unity.
        r = [S.Zero] * self.coord_sys.dim
        r[int(self.position) - 1] = S.One
        return r

    def _convert_to_rect(self, coord_sys):
        """
        coord_sys: an instance of subclass of CoordSys class
        converts a BaseVector in a given coordinate system into rectangular
        coordinates.
        """
        # First check that both coordinate systems have same position and
        # orientation
        flag = False
        if self.coord_sys.position != coord_sys.position:
            flag = True
        if self.coord_sys._dcm_global != coord_sys._dcm_global:
            flag = True

        if flag:
            raise ValueError("coord_sys doesn't have same position and \
                              orientation as coordinate system of self")

        if type(self.coord_sys) == type(coord_sys):
            return coord_sys.base_vectors[self.position - 1]
        rv = self.coord_sys._convert_to_rect(self, coord_sys)
        return rv

    def _convert_to_sph(self, coord_sys):
        """
        coord_sys: an instance of subclass of CoordSys class
        converts a BaseVector in a given coordinate system into rectangular
        coordinates.
        """
        # First check that both coordinate systems have same position and
        # orientation
        flag = False
        if self.coord_sys.position != coord_sys.position:
            flag = True
        if self.coord_sys._dcm_global != coord_sys._dcm_global:
            flag = True

        if flag:
            raise ValueError("coord_sys doesn't have same position and \
                              orientaion as coordinate system of self")

        if type(self.coord_sys) == type(coord_sys):
            return coord_sys.base_vectors[self.position - 1]
        rv = self.coord_sys._convert_to_sph(self, coord_sys)
        return rv

    def _convert_to_cyl(self, coord_sys):
        """
        coord_sys: an instance of subclass of CoordSys class
        converts a BaseVector in a given coordinate system into rectangular
        coordinates.
        """
        # First check that both coordinate systems have same position and
        # orientation
        flag = False
        if self.coord_sys.position != coord_sys.position:
            flag = True
        if self.coord_sys._dcm_global != coord_sys._dcm_global:
            flag = True

        if flag:
            raise ValueError("coord_sys doesn't have same position and \
                              orientaion as coordinate system of self")

        if type(self.coord_sys) == type(coord_sys):
            return coord_sys.base_vectors[self.position - 1]
        rv = self.coord_sys._convert_to_cyl(self, coord_sys)
        return rv

    def diff(self, s):
        """
        Differentiate a vector.
        s : a BaseScalar to differentiate against
        """
        if not isinstance(s, BaseScalar):
            raise TypeError("`s` should be a BaseScalar")
        if is_const_vect(self):
            return ZeroVector
        # If self.coord_sys if in rectangular coordinates, is_cons_vect will
        # return zero. So, we can be certain that  base vector is in some
        # coordinate system

        c_rect = s.coord_sys._change_coord_sys(CoordSysRect)
        vect = self.express(c_rect)

        return vect.diff(s)

    def __str__(self, printer=None):
        return self.coord_sys.name + "." + self.name

    __repr__ = __str__
    _sympystr = __str__
    _sympyrepr = __str__

    def _hashable_content(self):
        return (repr(self),) + tuple(sorted(self.assumptions0.items()))


class VectAdd(Add, Vector):
    """
    Container to hold added BaseVectors/VectMuls
    """

    def __new__(cls, *args, **options):
        for arg in args:
            try:
                arg.is_Vector
            except:
                raise TypeError(str(arg) + " is not a vector.")
        obj = super(VectAdd, cls).__new__(cls, *args, **options)

        is_commutative = True
        assumptions = {}
        assumptions['commutative'] = is_commutative

        obj._assumptions = StdFactKB(assumptions)
        return obj

    __init__ = Add.__init__

    def separate(self):
        # Flatten the VectMul so that there are no nested VectAdds
        vect = self.expand()

        # vect.args are either BaseVector or VectMul
        coord_sys_dict = {}
        coord_list = _all_coordinate_systems(vect)
        for c in coord_list:
            coord_sys_dict[c] = [[],[]]

        for arg in vect._all_args:
            # arg is either BaseVector or VectMul
            if isinstance(arg, BaseVector):
                coord_sys_dict[arg.coord_sys] = [[], [arg]]
            elif isinstance(arg, VectMul):
                arg_coord_dict = arg.separate()
                # Now, merge these dicts
                for c, v in arg_coord_dict.iteritems():
                    coord_sys_dict[c][0] += v[0]
                    coord_sys_dict[c][1] += v[1]
            else:
                # Shouldn't happen
                raise ValueError

        for c, v in coord_sys_dict.iteritems():
            coord_sys_dict[c][0] = _remove_duplicates(v[0])
            coord_sys_dict[c][1] = _remove_duplicates(v[1])

        return coord_sys_dict

    def expand(self):
        res = ZeroVector
        for arg in self._all_args:
            res = res + arg.expand()
        return res

    def factor_vect(self):
        # self can contain BaseVector or VectMul
        factor_dict = {}
        for arg in self._all_args:
            scalar = arg.scalar
            vector = arg.vector

            if not factor_dict.has_key(vector):
                factor_dict[arg.vector] = []
            if isinstance(arg, BaseVector) or isinstance(arg, VectMul):
                factor_dict[vector].append(scalar)
            else:
                # Should never get here
                raise ValueError

        res = []
        for arg, val in factor_dict.iteritems():
            res.append(VectMul(Add(*val), arg))
        return VectAdd(*res)

    @property
    def vector(self):
        vect = self.expand().factor_vect()
        if not isinstance(vect, VectAdd):
            return vect.vector
        else:
            raise TypeError("Cannot separate vector from scalar for VectAdd")

    @property
    def scalar(self):
        # First we expand
        vect = self.expand()
        # Now we factor the base vectors
        vect = vect.factor_vect()
        # Now, we have either a BaseVector, VectMul or VectAdd
        if not isinstance(vect, VectAdd):
            return vect.scalar
        else:
            raise TypeError("Cannot separate vector from scalar for VectAdd")

    @property
    def _all_args(self):
        return list(self.args)

    @property
    def components(self):
        if len(_all_coordinate_systems(self)) != 1:
            raise TypeError("The vector isn't in a single coordinate system.\
                             Cannot compute components")

        vect = self.expand().factor_vect()
        if not isinstance(vect, VectAdd):
            return vect.components

        # Fix this to accomodate n-dmesions
        r = [S.Zero] * 3
        for arg in vect._all_args:
            if isinstance(arg, BaseVector):
               # Component corresponding to it will be unity
               r[int(arg.position) - 1] = S.One
            elif isinstance(arg, VectMul):
               # The scalar part would be the component
               r[int(arg.vector.position) - 1] = arg.scalar

        return r

    @property
    def coord_sys(self):
        coord_list = _all_coordinate_systems(self)
        if not len(coord_list) == 1:
            raise ValueError("More than one coordinate systems")
        else:
            return coord_list.pop()

    def express(self, coord_sys):
        __doc__ = express.__doc__
        return express(self, coord_sys)

    def diff(self, s):
        """
        Differentiate a vector.
        s : a BaseScalar
        """
        # check the base scalar
        if not isinstance(s, BaseScalar):
            raise TypeError("`s` must be a BaseScalar")

        # Now, self is a VectAdd - so it is a combination of BaseVectors
        # and VectMul (which might have nested VectAdds)
        ret = ZeroVector
        for arg in self._all_args:
            # arg is BaseVectot or VectMul
            ret = ret + arg.diff(s)

        return ret

    def __str__(self, printer=None):
        # VectAdd can contain either BaseVectors or VectMuls
        # TODO : This method works but can be improved.
        ret_str = ''
        for arg in self._all_args:
            # arg must be BaseVector or VectMul
            assert isinstance(arg, BaseVector) or isinstance(arg, VectMul)
            ret_str += " + " + repr(arg)

        # Remove the leading " + "
        ret_str = ret_str[3:]
        return ret_str

    __repr__ = __str__
    _sympystr = __str__
    _sympyrepr = __str__


class VectMul(Mul, Vector):
    """
    Container to hold Vectors/scalar products
    """
    def __new__(cls, *args, **options):
        counter = 0
        for arg in args:
            if arg.is_Vector:
                counter += 1
        if counter > 1:
            raise TypeError("Cannot multiply two or more vectors.")
        obj = super(VectMul, cls).__new__(cls, *args, **options)

        is_commutative = True
        assumptions = {}
        assumptions['commutative'] = is_commutative
        obj._assumptions = StdFactKB(assumptions)
        return obj

    __init__ = Mul.__init__

    def separate(self):
        # First we flatten things out - so that there are no nested VectAdd
        vect = self.expand()
        if isinstance(vect, VectAdd) or isinstance(vect, BaseVector):
            #This causes max recursion depth
            return vect.separate()

        # Now we are sure that vect is just VectMul - no nesting
        scalar = vect.scalar
        vector = vect.vector

        coord_sys_dict = {}
        coord_list = _all_coordinate_systems(vect)
        for c in coord_list:
            coord_sys_dict[c] = [[], []]

        # First process the scalar
        if isinstance(scalar, BaseScalar):
            coord_sys_dict[scalar.coord_sys][0].append(scalar)
        else:
            for arg_scalar in scalar.args:
                if isinstance(arg_scalar, BaseScalar):
                    coord_sys_dict[scalar.coord_sys][0].append(scalar)

        # Now process the vector part
        # The vector part is necessarily BaseVector becauase we have
        # used expand before
        assert isinstance(vector, BaseVector)
        coord_sys_dict[vector.coord_sys][1].append(vector)

        return coord_sys_dict

    @property
    def coord_sys(self):
        """
        If the object is only defined in one coordinate system
        then return that cooordinate system.
        Cases:
        1. scalar*BaseVector
        2. scalar*VectAdd
        """
        vect = self.vector
        if isinstance(vect, BaseVector):
            return vect.coord_sys

        elif isinstance(vect, VectAdd):
            coord_list = _all_coordinate_systems(vect)
            if len(coord_list) == 1:
                return coord_list[0]

            else:
                raise ValueError("Cannot recognize " + str(vect))

    def expand(self):
        # VectMul can be either BaseVector*scalar
        # or VectAdd*scalar or BaseVector
        # We process all cases
        scalar = self.scalar
        vector = self.vector
        ret = ZeroVector

        # Only if the scalar is of type Add can it be distributed
        if not isinstance(scalar, Add):
            for arg in vector._all_args:
                if isinstance(arg, BaseVector):
                    v = VectMul(arg, scalar)
                elif isinstance(arg, VectMul):
                    v = VectMul(arg, scalar).expand()
                else:
                    # Shouldn't happen
                    raise TypeError
                ret = ret + v
            return ret

        # Now we know that scalar must be Add. So it has to distribute
        if isinstance(scalar, BaseScalar):
            scalar_args = [scalar]
        else:
            scalar_args = scalar.args

        if not scalar_args:
            # scalar.args was empty. Get atoms
            scalar_args = list(scalar.atoms())

        if isinstance(vector, BaseVector):
            ret = ZeroVector
            for arg in scalar_args:
                ret = ret + VectMul(arg, vector)
            return ret

        if isinstance(vector, VectAdd):
            ret = ZeroVector
            vector = vector.expand()
            # Now, vect.args should be either BaseVector or VectMul
            for arg in vector._all_args:
                ret = ret +  (arg * scalar).expand()
            return ret

    def factor_vect(self):
        # This is of type scalar * BaseVector and scalar * VectAdd
        scalar = self.scalar
        vector = self.vector
        if isinstance(vector, BaseVector):
            # Just factor out the scalar and return
            return VectMul(scalar.factor(), vector)
        elif isinstance(vector, VectAdd):
            return VectMul(scalar.factor(), vector.factor_vect())
        else:
            # Shouldn't happen
            raise ValueError

    def as_coeff_Mul(self):
        __doc__ = super(VectMul, self).as_coeff_mul.__doc__
        # This is a temporary hack.
        # When flatten method encounters VectMul such as 3 * c0.e_x, it
        # creates a dict where c0.e_x is stored as being repeated thrice.
        # This wouldn't happen if as_coeff_mul returns itself (as happens
        # in case of x * c0.e_x).
        # This can also be corrected by writing flatten methods for vector
        # classes.
        return S.One, self

    @property
    def vector(self):
        """
        Return the vector contained in VectMul.
        Returns a BaseVector or a VectAdd
        """
        for arg in self.args:
            if isinstance(arg, BaseVector) or isinstance(arg, VectAdd):
                return arg
        # If we haven't found a vector, raise error
        raise TypeError(str(self) + " doesn't have a BaseVector/VectAdd \
                        its args")

    @property
    def _all_args(self):
        return [self]

    @property
    def scalar(self):
        """
        Returns the scalar contained in a VectMul object
        """
        # It can be either scalar*BaseVector or scalar*VectAdd
        # We treat both cases
        vector = self.vector
        scalar = S.One
        for arg in self.args:
            if arg.is_Vector:
                continue
            scalar = scalar * arg

        # To avoid confusion between attribute vector and var vector
        v = vector
        if isinstance(v, VectAdd):
            # vector is a VectAdd.
            v = v.factor_vect()
            # We get either scalar1 * BaseVector or scalar1 * VectAdd as result
            # In second case, we just raise an error
            # In the first case, we get the scalar part from it and multiply
            # that to already find scalat from above
            # The vector part of it is a BaseVector or a VectAdd

            if isinstance(v.vector, VectAdd):
                raise TypeError("More than one component present in vector")
            # Now, v is just a VectMul and its vector component
            # must be BaseVector
            else:
                scalar = scalar * v.scalar

        return scalar

    @property
    def components(self):
        vect = self.expand().factor_vect()
        if not isinstance(vect, VectMul):
            return vect.components

        # Since vect is just VectMul, so return the scalar part
        r = [S.Zero]*vect.coord_sys.dim
        r[int(vect.vector.position) - 1] = vect.scalar
        return r

    def express(self, coord_sys):
        __doc__ = express.__doc__
        return express(self, coord_sys)

    def diff(self, s):
        """
        Differentiate a vector.
        s : a BaseScalar
        """
        # (fg)' = f'g + fg'
        scalar = self.scalar
        vector = self.vector

        ret = ZeroVector
        ret = ret + _diff_scalar(scalar, s) * vector + vector.diff(s)
        return ret

    @property
    def func(self):
        # VectMul(scalar, VectAdd) does not work (see comments in _vect_mul helper).
        # So, we need to pass it through _vect_mul instead
        return VectMul._func_vect_mul

    @classmethod
    def _func_vect_mul(cls, *args):
        # args can contain a large number of scalars but only one vector
        ret = S.One
        flag = False
        for arg in args:
            if arg.is_Vector:
                flag = False
            if not flag:
                ret = ret * arg
            else:
                raise TypeError("Cannot multiply two or more vectors")
        return ret

    def __str__(self, printer=None):
        # VectMul can contain another VectAdd
        # TODO : This method works but can be improved.
        ret_str = ''
        scalar = self.scalar
        vector = self.vector
        if isinstance(scalar, Add):
            ret_str = "(" + repr(scalar) + ") * "
        else:
            ret_str = repr(scalar) + " * "

        if isinstance(vector, BaseVector):
            ret_str += repr(vector)
        elif isinstance(vector, VectAdd):
            ret_str += "(" + repr(vector) + ")"
        else:
            # Shouldn't happen
            raise TypeError

        return ret_str

    __repr__ = __str__
    _sympystr = __str__
    _sympyrepr = __str__


class ZeroVectorClass(Zero, Vector):
    """
    A class to represent the zero vector
    """
    is_Vector = True
    # Note that _op_priority is higher than other vector classes
    # This is for the usage of _zero_vector decorator
    _op_priority = 11.1

    def __neg__(self):
        return VectMul(S.NegativeOne, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__radd__')
    def __add__(self, other):
        return _vect_add(S.Zero, other)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__add__')
    def __radd__(self, other):
        return _vect_add(other, S.Zero)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rsub__')
    def _sub__(self, other):
        return _vect_add(S.Zero, -other)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return _vect_add(other, -S.Zero)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return _vect_mul(S.Zero, other)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return _vect_mul(other, S.Zero)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        return _vect_div(S.Zero, other)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        raise TypeError("Cannot divide by vector")

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    @property
    def components(self):
        return (S.Zero, S.Zero, S.Zero)

    def express(self, coord_sys):
        return ZeroVector

    def expand(self):
        return self

    def factor_vect(self):
        return self

def _vect_add(one, other):
    # We are adding this method to check for cases involving S.Zero
    if ((one.is_Vector or one == S.Zero) and
        (other.is_Vector or other == S.Zero)):
        return VectAdd(one, other)
    else:
        raise TypeError("Cannot add a scalar with a vector")

def _vect_mul(one, other):
    # TODO: We cannot yet initialized objects of type scalar * VectAdd
    # because of the flatten method being called in the AssocOp.__new__
    # method. This results in an Add being created internally. The fix
    # consist of writing a flatten method for all vector classes.

    # For the time begin, we manually 'expand out' the VectAdd, if there
    # is one.
    if isinstance(one, BaseVector) or isinstance(other, BaseVector):
        ret = VectMul(one, other)
    elif isinstance(one, VectAdd) or isinstance(other, VectAdd):
        ret = ZeroVector
        # First get the vector and scalar
        vector = one if isinstance(one, VectAdd) else other
        scalar = one if vector == other else other
        for arg in vector._all_args:
            t = scalar * arg
            ret = ret + t
    elif isinstance(one, VectMul) or isinstance(other, VectMul):
        ret = VectMul(one, other)
    else:
        # scalar * scalar
        ret = Mul(one, other)
    return ret

def _vect_div(one, other):
    if one.is_Vector and other.is_Vector:
        raise TypeError("Cannot divide two vectors")
    if not one.is_Vector and not other.is_Vector:
        # This shouldn't happen unless _vect_div is called directly
        return Mul(one, Pow(other, S.NegativeOne))
    # We cannot have scalar/vector
    if other.is_Vector:
        raise ValueError("Cannot divide by vector")
    else:
        # Now we know that either one is a vector. Remaining is scalar
        return _vect_mul(one, Pow(other, S.NegativeOne))

def _vect_add_const(one, other, coord_sys):
    # Constant vectors are different from other vectors in the sense
    # that they cannot be moved in space. They, thus, need to be handled
    # differently.
    assert one.is_Vector and other.is_Vector

    # Check that both the vectors are non zero
    if not one or not other:
        return _vect_add(one, other)

    if one.coord_sys == other.coord_sys:
        return _vect_add(one, other).factor_vect()
    other = other.express(coord_sys)
    one_comp = one.components
    other_comp = other.components
    ret = [S.Zero] * one.coord_sys.dim

    for i in ret:
        ret[i] = one_comp[i] + other_comp[i]
    return to_vector(ret, coord_sys)

def to_vector(components, coord_sys):
    """
    Create a vector from given components.
    components : a tuple/list of components
    coord_sys : an instance of subclass of CoordSys class. This object
                will be used get the base vectors for the vector
    Returns an object with is_Vector == True
    """
    ret = ZeroVector
    for i, comp in enumerate(components):
        ret = ret + comp * coord_sys.base_vectors[i]
    return ret

def express(vect, coord_sys):
    """
    coord_sys: Instance of subclasses of CoordSys
    Express a vector in a given coordinate system.
    """
    coord_list = _all_coordinate_systems(vect)
    subs_dict = {}
    # We proceed to build a dictionary that, upon substituting in the
    # original vector, will give us a vector in desired coord_sys

    if len(coord_list) == 1 and vect.coord_sys == coord_sys:
        return vect
    coord_sys_t = coord_sys._change_coord_sys(CoordSysRect, 'coord_sys_t')
    vect = vect.expand()
    coord_sys_dict = _coord_sys_dict(vect)
    # First express everything in rect coordinates
    for arg in vect._all_args:
        # arg is necessarily a VectMul or a BaseVector
        assert isinstance(arg, BaseVector) or isinstance(arg, VectMul)
        scalar = arg.scalar
        vector = arg.vector

        # First process the scalar
        if isinstance(scalar, BaseScalar):
            # Because BaseScalar.args returns coord_sys as well
            c_rect = coord_sys_dict[scalar.coord_sys]
            subs_dict[scalar] = scalar._convert_to_rect(c_rect)

        else:
            for arg_scalar in _all_base_scalars(scalar):
                if not subs_dict.has_key(arg_scalar):
                    c_rect = coord_sys_dict[arg_scalar.coord_sys]
                    subs_dict[arg_scalar] = arg_scalar._convert_to_rect(c_rect)
        # Now process the vector
        # vector is necessarily BaseVector
        assert isinstance(vector, BaseVector)

        c_rect = coord_sys_dict[vector.coord_sys]
        subs_dict[vector] = vector._convert_to_rect(c_rect)

    # Now performig the substitution, we have changed all involved
    # variables into rectangular coordinates
    vect = vect.subs(subs_dict)
    subs_dict = {}

    vect = vect.expand()

    # First phase complete
    # Now, we find subs for base vectors and substitute it in
    for arg in vect._all_args:
        assert isinstance(arg, BaseVector) or isinstance(arg, VectMul)
        vector = arg.vector
        subs_dict[vector] = CoordSys._convert_base_vect_rect(vector,
                                                             coord_sys_t)
    # Now subs in the vectors
    vect = vect.subs(subs_dict)

    subs_dict = {}
    vect = vect.expand()

    # Now we need to find the subs_dict for each BaseScalar
    for arg in vect._all_args:
        assert isinstance(arg, BaseVector) or isinstance(arg, VectMul)
        scalar = arg.scalar

        # Process the scalar
        if isinstance(scalar, BaseScalar):
            # Because BaseScalar.args returns coord_sys as well
            subs_dict[scalar] = CoordSys._convert_base_sclr_rect(scalar,
                                                            coord_sys_t)
        else:
            for arg_scalar in _all_base_scalars(scalar):
                if not subs_dict.has_key(arg_scalar):
                    subs_dict[arg_scalar] = \
                        CoordSys._convert_base_sclr_rect(arg_scalar, coord_sys_t)

    # Now performig the substitution, we have changed all involved
    # variables into rectangular coordinates
    vect = vect.subs(subs_dict)

    # Now we need to convert back from coord_sys_t to coord_sys
    subs_dict = {}
    vect = vect.expand()

    for arg in vect._all_args:
        # arg is necessarily a VectMul or a BaseVector
        assert isinstance(arg, BaseVector) or isinstance(arg, VectMul)
        scalar = arg.scalar
        vector = arg.vector

        # First process the scalar
        if isinstance(scalar, BaseScalar):
            exec "conv_func = scalar._convert_to_" + coord_conv[type(coord_sys)]
            # Because BaseScalar.args returns coord_sys as well
            subs_dict[scalar] = conv_func(coord_sys)

        else:
            for arg_scalar in _all_base_scalars(scalar):
                if not subs_dict.has_key(arg_scalar):
                    exec "conv_func = arg_scalar._convert_to_" + \
                        coord_conv[type(coord_sys)]
                    subs_dict[arg_scalar] = \
                                            arg_scalar._convert_to_rect(coord_sys)

        # Now process the vector
        # vector is necessarily BaseVector
        assert isinstance(vector, BaseVector)
        exec "conv_func = vector._convert_to_" + coord_conv[type(coord_sys)]
        subs_dict[vector] = conv_func(coord_sys)

    # Now performig the substitution, we have changed all involved
    # variables into rectangular coordinates
    vect = vect.subs(subs_dict)
    subs_dict = {}

    return vect.factor_vect()

def _express_scalar(expr, coord_sys):
    """
    Express a scalar in the given coord_sys
    expr : a SymPy expression
    coord_sys : an instance of subclass of CoordSys class
    """
    # First check if expr is entirely in coord_sys
    coord_list = _all_coordinate_systems_scalar(expr)
    if len(coord_list) == 1 and coord_sys == coord_list[0]:
        return expr

    expr = expr.expand()
    subs_dict = {}
    all_base_scalars = _all_base_scalars(expr)
    coord_dict = _coord_sys_dict_scalar(expr)
    for arg in all_base_scalars:
        # Convert base scalars to c_rect
        c_rect = coord_dict[arg.coord_sys]
        subs_dict[arg] = arg._convert_to_rect(c_rect)
    expr = expr.subs(subs_dict)
    subs_dict = {}

    # Create a coordinate from coord_sys - just that it is rectangular
    coord_sys_t = coord_sys._change_coord_sys(CoordSysRect, 'coord_sys_t')
    all_base_scalars = _all_base_scalars(expr)

    for arg in all_base_scalars:
        subs_dict[arg] = CoordSys._convert_base_sclr_rect(arg, coord_sys_t)

    expr = expr.subs(subs_dict)
    subs_dict = {}
    all_base_scalars = _all_base_scalars(expr)

    # Now convert back to coord_sys
    for arg in all_base_scalars:
        subs_dict[arg] = arg._convert_to_rect(coord_sys)

    expr = expr.subs(subs_dict)
    return expr

def _coord_sys_dict(vect):
    coord_list = _all_coordinate_systems(vect)
    ret = {}
    for c in coord_list:
        ret[c] = c._change_coord_sys(CoordSysRect)
    return ret

def _coord_sys_dict_scalar(scalar):
    coord_list = _all_coordinate_systems_scalar(scalar)
    ret = {}
    for c in coord_list:
        ret[c] = c._change_coord_sys(CoordSysRect)
    return ret

def _remove_duplicates(l):
    ret = []
    for i in l:
        if i not in ret:
            ret.append(i)
    return ret

ZeroVector = ZeroVectorClass()

coord_conv = {
    CoordSysRect: 'rect',
    CoordSysCyl: 'cyl',
    CoordSysSph: 'sph'
    }
