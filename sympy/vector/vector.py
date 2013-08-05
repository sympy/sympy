from sympy.core.assumptions import StdFactKB
from sympy.matrices import Matrix, eye
from sympy.core import (Basic, Expr, Dummy, Function, Symbol, symbols,
                        sympify, diff, Pow, Mul, Add, S, AtomicExpr)
from sympy.core.numbers import Zero
from sympy.solvers import solve
from sympy.simplify import simplify, trigsimp
from sympy.core.decorators import call_highest_priority, _sympifyit
from sympy.core.compatibility import as_int
from functools import wraps

import copy

def base_scalars(names, coord_sys, ind=1):
    """
    Creates BaseScalars (coordinates) in a tuple for convenience.

    Takes the names of the BaseScalars, the coordinate system the are
    defined in, and an initial index.
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
    Creates Vectors in a tuple for convenience.

    Takes the names of the Vectors, the coordinate system the are defined
    in, and an initial index.
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
    Generalized dot product.
    """
    ret = S.Zero
    # Get two lists - each having vectors separated by coordinate system
    a_vectors = _separate_to_vectors(vect_a)
    b_vectors = _separate_to_vectors(vect_b)

    if not coord_sys:
        if (len(a_vectors) == 1 and len(b_vectors) == 1
            and a_vectors[0].coord_sys == b_vectors[0].coord_sys):
           coord_sys = a_vectors[0].coord_sys
        else:
           raise ValueError("Coordinate system not provided.")

    r_vect_a = None
    for vector in a_vectors:
        r_vect_a += vector.express(coord_sys)
    r_vect_b = None
    for vector in b_vectors:
        r_vect_b += vector.express(coord_sys)
    # Now we have two vectors, each in the provided coord_sys
    return _dot_same(r_vect_a, r_vect_b)

def _separate_to_vectors(vect):
    coord_dict = vect.separate()
    r = [v for v in coord_dict.values()]
    return r

def _all_coordinate_systems(vector):
    vector = vector.expand()
    coord_list = []
    # all_args is a separate method that return only the vector args
    # arg is either BaseVector or VectMul - becuase we have expanded the vector
    for arg in vector._all_args:
        if isinstance(arg, BaseVector):
            coord_list.append(arg.coord_sys)
        elif isinstance(arg, VectMul):
            coord_list = coord_list + _coord_sys_scalar_list(arg.scalar)
            coord_list.append(arg.vector.coord_sys)
        else:
            # Shouldn't happen
            raise ValueError("Couldn't expand vector")

    return list(set(coord_list))

def _coord_sys_scalar_list(scalar):
    # scalar can be any sympy basic type
    # return a list of coordinate systems contained within scalar
    coord_list = []
    if isinstance(scalar, BaseScalar):
        coord_list.append(scalar.coord_sys)
        return coord_list
    for arg in scalar.args:
        if isinstance(arg, BaseScalar):
            coord_list.append(arg.coord_sys)
        return coord_list


class BaseScalar(Expr):
    """
    BaseScalar instances are used to express coordinate variables for field.
    Not to be instantiated by the user.
    """
    is_comparable = False

    @property
    def _diff_wrt(self):
        """Allow derivatives wrt to BaseScalars"""
        return True

    def __init__(self, name, coord_sys, position, **assumptions):
        if (position not in [1, 2, 3]
            and not isinstance(coord_sys, CoordSysRect)):
            raise ValueError("Position of scalar not specified. \
                            See `position` in docstring of `BaseScalar`")
        is_commutative = assumptions.get('commutative', True)
        assumptions['commutative'] = is_commutative

        if name is None:
            name = 'Dummy_' + str(Dummy._count)
            Dummy._count += 1
        super(BaseScalar, self).__init__()
        self.coord_sys = coord_sys
        self.name = name
        self.position = position
        self._assumptions = StdFactKB(assumptions)

    def __str__(self):
        return self.coord_sys.name + "." + self.name

    __repr__ = __str__

    @staticmethod
    def coord_rel_rect(base_scalar, coord_sys):
        """
        base_scalar: A BaseScalar object in CoordSysRect
        coord_sys: A CoordSys object
        Returns a relation between the given base_scalar and the
        given coord_sys
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


class CoordSys(Basic):
    """
    Superclass for all CoordSys<type> classes. Not to be intialized
    directly.
    """
    def __new__(cls, name=None, dim=None, position=None, position_coord=None,
                orient_type=None, orient_amount=None, rot_order=None,
                parent=None, basis_vectors=None, coordinates=None):
        return Basic.__new__(cls, name, dim, position, position_coord,
                             orient_type, orient_amount, rot_order, parent,
                             basis_vectors, coordinates)

    def __init__(
            self, name=None, dim=3, position=None, position_coord='rect',
            orient_type=None, orient_amount=None, rot_order=None,
            parent=None, basis_vectors=None, coordinates=None
            ):

        if name is None:
            name = 'CoordSys_' + str(Dummy._count)
            Dummy._count += 1

        self.name = name
        self.dim = int(dim)
        #Moved setting of parent to here
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
            if not position.is_constant:
                raise TypeError("Position variables cannot be BaseScalars")

            if parent:
                c_rect = CoordSysRect('c_rect')
                parent_pos = parent.position.express(c_rect)
                self_pos = position.express(c_rect)
                self_pos = parent_pos + self_pos
                self.position = self_pos.factor()
            else:
                self.position = position

        self._dcm_global = eye(3)
        self._dcm_parent = eye(3)

        if orient_type:
            self._check_orient_raise(orient_type, orient_amount, rot_order)
            if parent:
                self._dcm_parent = self._dcm_parent_func(
                    orient_type, orient_amount, rot_order
                )

            self._dcm_global = self._dcm_global_method(orient_type,
                                                       orient_amount)

    #def __getattr__(self, i):
    #    return None
    #    return self._coordinates[int(i[1:]) - 1]

    def __str__(self, printer=None):
        return self.name
    __repr__ = __str__
    _sympystr = __str__
    _sympyrepr = __str__

    def _dcm_parent_method(self, orient_type, amounts, rot_order=''):
        orient_type = orient_type.capitalize()
        if orient_type == 'Axis':
            if not rot_order == '':
                raise TypeError('Axis orientation takes no rotation order')
            if not (isinstance(amounts, (list, tuple)) & (len(amounts) == 2)):
                raise TypeError('Amounts are a list or tuple of length 2')
            theta = amounts[0]
            axis = amounts[1]
            axis = axis.express(self.parent).normalize().as_mat()
            axis = axis.args[0][0]
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
        if hasattr(self, 'parent'):
            # Parent given therefore the given angle is wrt parent.
            # DCM(global<-parent)*DCM(parent<-self) == DCM(global <- self)
            parent = self.parent
            r = parent.dcm().T * parent.dcm(self)
            return r
        else:
            # Parent not set. So, directly initialize the dcm wrt global
            # Using sympy.physics.mechanics logic here.
            if orient_type == 'AXIS':
                if not rot_order == '':
                    raise ValueError('Axis orientation \
                                      takes no rotation order')
                if not (isinstance(amounts, (list, tuple))
                        & (len(amounts) == 2)):
                    raise TypeError('Amounts are a list or tuple of length 2')
                theta = amounts[0]
                axis = amounts[1]
                axis = axis.in_global()
                axis = axis.normalize().as_mat()
                axis = axis.args[0][0]
                global_orient = ((eye(3) - axis * axis.T) * cos(theta) +
                        Matrix([[0, -axis[2], axis[1]],
                                [axis[2], 0, -axis[0]],
                                [-axis[1], axis[0], 0]]) * \
                                 sin(theta) + axis * axis.T)
                return global_orient
            if orient_type == 'BODY':
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
        """DCM about one of the 3 axes"""
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
        """
        Return a list of base scalars
        """
        return self._coordinates

    @property
    def base_vectors(self):
        """
        Return a list of base scalars
        """
        return self._basis

    def dcm(self, other='global'):
        """
        Returns the direction cosine matrix for converting from `other`
        to self, i.e.,
        self.xyz = self.dcm(other)*other.xyz
        """
        if other == 'global':
            return self._dcm_global

        # First check self and other are parent/child frames
        if self.parent == other:
            return self._dcm_parent
        elif other.parent == self:
            return self._dcm_parent.T
        # TODO : Implement path A->B->C->D rather than A->global->D
        else:
            # Now we can either traverse the frames to find the nearest
            # common frame. This will require multiplying many matrices.
            # So, we go via the global frame
            r = other.dcm() * self.dcm().T
            r.applyfunc(lambda i: trigsimp(method='fu'))
            return r

    def orientnew(self, name=None, orient_type=None, orient_amount=None,
                 rot_order=None):
        """
        Create a coordinate system with a new orientation from an exisitng
        coordinate system.
        """
        newframe = copy.copy(self)

        if name == None:
            name = 'CoordSys_' + str(Dummy._count)
            Dummy._count += 1
        newframe.name = name
        newframe._check_orient_raise(orient_type, orient_amount, rot_order)
        newframe.parent = self
        newframe._dcm_parent = newframe._dcm_parent(orient_type,
                                                    orient_amount, rot_order)
        newframe._dcm_global = newframe._dcm_global(orient_type, orient_amount)
        return newframe

    def posnew(self, position):
        if not position.is_Vector:
            raise TypeError("vector expected for position")
        if not position.is_constant:
            raise TypeError("Position variables cannot be BaseScalars")

        parent = self
        newframe = copy.copy(self)
        c_rect = CoordSysRect('c_rect')

        parent_pos = parent.position.express(c_rect)
        newframe_pos = position.express(c_rect)
        newframe_pos = parent_pos + newframe_pos

        newframe.position = newframe_pos.factor()
        return newframe

    def _convert_coord_sys(self, coord_sys, name=None):
        if not name:
            name = "coord_sys_" + str(Dummy._count)
            Dummy._count += 1
        CoordSysClass = coord_sys
        ret_coord_sys = CoordSysClass(name=name, position=self.positon,
                                    orient_type=self.orient_type,
                                    orient_amount=self.orient_amount,
                                    rot_order=self.rot_order,
                                    parent=self.parent,
                                    basis_vectors=self.basis_vectors,
                                    coordinates=self.coordinates)
        return ret_coord_sys

    @staticmethod
    def _pos_to_rect(coord, conv_from):
        if conv_from == 'cyl':
            x = coord[0] * cos(coord[1])
            y = coord[0] * sin(coord[1])
            z = coord[2]
            return (x, y, z)

        if conv_from == 'sph':
            x = coord[0] * sin(coord[1]) * cos(coord[2])
            y = coord[0] * sin(coord[1]) * sin(coord[2])
            z = coord[0] * cos(coord[1])
            return (x, y, z)

    @staticmethod
    def _pos_to_cyl(coord, conv_from):
        if conv_from == 'rect':
            rho = sqrt(coord[0]**2 + coord[1]**2)
            phi = coord[2]
            z = coord[2]
            return (rho, phi, z)

        if conv_from == 'sph':
            rho = coord[0] * sin(coord[1])
            phi = coord[2]
            z = coord[0] * cos(coord[1])
            return (rho, phi, z)

    @staticmethod
    def _pos_to_sph(coord, conv_from):
        if conv_from == 'rect':
            r = sqrt(coord[0]**2 + coord[1]**2 + coord[2]**2)
            theta = atan(sqrt(coord[0]**2 + coord[1]**2)/coord[2])
            phi = atan(coord[1]/coord[0])
            return (r, theta, phi)

        if conv_from == 'cyl':
            r = sqrt(coord[0]**2 + coord[2]**2)
            theta = atan(coord[0]/coord[2])
            phi = coord[1]
            return (r, theta, phi)

    @staticmethod
    def _pos_at(vect, new_coord_sys):
        # 1. Convert location of both the points in space to rect
        old_coord_sys = vect.coord_sys
        pos_old = old_coord_sys.position
        pos_new = new_coord_sys.position
        if pos_old == S.Zero:
            pos_old = (0, 0, 0)
        if pos_new == S.Zero:
            pos_new = (0, 0, 0)
        pos_old = CoordSys._pos_to_rect(pos_old,
                                        coord_conv[type(old_coord_sys)])
        pos_new = CoordSys._pos_to_rect(pos_new,
                                        coord_conv[type(new_coord_sys)])
        # 2. Determine shifting relations
        # pos1 ->x, y, z and pos2 ->X, Y, Z
        # We need to go from pos1 to pos2
        if old_coord_sys.dim > 3 or new_coord_sys.dim > 3:
            raise NotImplementedError
        x0, y0, z0 = old_coord_sys.base_scalars
        X, Y, Z = new_coord_sys.base_scalars

        ex0, ey0, ez0 = old_coord_sys.base_scalars
        eX, eY, eZ = new_coord_sys.base_scalars

        x = X + (pos_old[0] - pos_new[0])
        y = Y + (pos_old[1] - pos_new[1])
        z = Z + (pos_old[2] - pos_new[2])

        # 3. Convert vect to rectangular coordinates
        dummy_system_rect = CoordSysRect('DummyRect')
        vect = vect._convert_to_rect(vect, dummy_system_rect)
        # 4. Subs for x, y, z in terms of X, Y, Z in vect
        subs_dict = { x0: x,
                      y0: y,
                      z0: z,
                      ex0: eX,
                      ey0: eY,
                      ez0: eZ }
        vect = vect.subs(subs_dict)
        # Now we have the vector field in the new coordinates - just
        # that it is in rectangular system.
        return vect

    @staticmethod
    def _orient_along(vect, coord_sys):
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
    def _change_sys(vect, coord_sys, func):
        """
        Change the coordinate systems of the the vector.
        """
        return func(vect, coord_sys)


class CoordSysRect(CoordSys):
    """
    The rectangular coordinate system.
    """
    def __init__(self, *args, **kwargs):
        if 'coordinates' not in kwargs.keys():
            try:
                kwargs['coordinates'] = \
                        ['x' + str(i + 1) for i in range(kwargs['dim'])]
            except:
                kwargs['coordinates'] = ['x', 'y', 'z']
        if 'basis_vectors' not in kwargs.keys():
            try:
                kwargs['basis_vectors'] = \
                    ['e_' + str(i + 1) for i in range(kwargs['dim'])]
            except:
                kwargs['basis_vectors'] = ['e_x', 'e_y', 'e_z']

        super(CoordSysRect, self).__init__(*args, **kwargs)
        self.h_list = tuple([S.One] * self.dim)


    @staticmethod
    def _convert_to_sph(vector, coord_sys):
        """
        vector : a VectAdd
        coord_sys : a CoordSys object - the coordinate system to convert to
        returns an object with is_Vector == True
        """
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
        __doc__ = CoordSysRect._convert_to_sph.__doc__

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

class CoordSysSph(CoordSys):
    """
    The spherical polar coordinate system.
    """
    def __init__(self, *args, **kwargs):
        if 'coordinates' not in kwargs.keys():
            kwargs['coordinates'] = ['r', 'theta', 'phi']
        if 'basis_vectors' not in kwargs.keys():
            kwargs['basis_vectors'] = ['e_r', 'e_theta', 'e_phi']
        super(CoordSysSph, self).__init__(*args, **kwargs)

        self.dim = S(3)

        if name is None:
            name = 'CoordSysSph_' + str(Dummy._count)
            Dummy._count += 1

        self.name = name

        self.one, self.two, self.three = symbols('_1 _2 _3')
        self.h_list = (S.One,
                       self.one,
                       self.one*sin(self.two))

    @staticmethod
    def _convert_to_rect(vector, coord_sys):
        __doc__ = CoordSysRect._convert_to_sph.__doc__
        Ar, At, Ap = vector.components
        r, theta, phi = vector.coord_sys.base_scalars
        x, y, z = coord_sys.base_scalars
        subs_dict = { r: sqrt(x**2 + y**2 + z**2),
                      theta: atan(sqrt(x**2 + y**2)/z),
                      phi: atan(y/x) }

        Ar = Ar.subs(subs_dict)
        At = At.subs(subs_dict)
        Ap = Ap.subs(subs_dict)
        # Rect to Sph conv matrix. So we need to take tranpose.
        mat =  Matrix([
                      [sin(theta)*cos(phi), sin(theta)*sin(phi),  cos(theta)],
                      [cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)],
                      [          -sin(phi),            cos(phi),           0]
                     ])
        mat = mat.T
        mat = mat.subs(subs_dict)
        mat.simplify()
        r = mat*Matrix([ [Ar], [At], [Ap] ])
        res = []
        for i, vect in enumerate(r._mat):
            res.append(VectMul(vect, coord_sys.base_vectors[i]))
        return VectAdd(*res)

    @staticmethod
    def _convert_to_cyl(vector, coord_sys):
        __doc__ = CoordSysRect._convert_to_sph.__doc__
        Ar, At, Ap = vector.components
        r, theta, phi = vector.coord_sys.base_scalars
        rho, _phi, z = base_scalars
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


class CoordSysCyl(CoordSys):
    """
    The cylindrical coordinate system.
    """
    def __init__(self, *args, **kwargs):
        if 'coordinates' not in kwargs.keys():
            kwargs['coordinates'] = ['rho', 'phi', 'z']
        if 'basis_vectors' not in kwargs.keys():
            kwargs['basis_vectors'] = ['e_rho', 'e_phi', 'e_z']
        super(CoordSysCyl, self).__init__(*args, **kwargs)

        self.dim = S(3)

        if name is None:
            name = 'CoordSysCyl_' + str(Dummy._count)
            Dummy._count += 1

        self.name = name

        self.one, self.two, self.three = symbols('_1 _2 _3')
        self.h_list = (S.One,
                       self.one,
                       S.One)


    @staticmethod
    def _convert_to_rect(vector, coord_sys):
        """
        vector : a VectAdd
        coord_sys : an instance of subclass of CoordSys class
        returns an object with is_Vector == True
        """
        Ar, Ap, Az = vector.components
        rho, phi, z = vector.base_scalars
        x, y, _z = base_scalars
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
            res.append(VectMul(vect, base_vectors[i]))
        return VectAdd(*res)

    @staticmethod
    def _convert_to_sph(vector, coord_sys):
        __doc__ = CoordSysRect._convert_to_sph.__doc__
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


class Vector(Expr):
    """
    Class instances will represent base vectors
    """
    is_Vector = True
    _op_priority = 11.0


    #def __init__(self, name, coord_sys, position):
    #    self.name = name
    #    if not isinstance(coord_sys, CoordSys):
    #        raise TypeError("coord_sys must be isinstance(CoordSys)")
    #    self.coord_sys = coord_sys
    #     self.position = position
    #     if position > coord_sys.dim:
    #         raise ValueError("Vector outside coordinate system\
    #                           dimensionality")

    #def __str__(self, printer=None):
    #    return self.name

    #_sympystr = __str__
    #_sympyrepr = __str__


    def __neg__(self):
        return VectMul(S.NegativeOne, self)

    @call_highest_priority('__radd__')
    def __add__(self, other):
        return _vect_add(self, other)

    @call_highest_priority('__add__')
    def __radd__(self, other):
        return _vect_add(other, self)

    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return _vect_add(self, -other)

    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return _vect_add(other, -self)

    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return _vect_mul(self, other)

    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return _vect_mul(other, self)

    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        return _vect_div(self, other)

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


class BaseVector(Vector, Symbol):
    def __init__(self, name, coord_sys, position):
        self.name = name
        if not isinstance(coord_sys, CoordSys):
            raise TypeError("coord_sys must be isinstance(CoordSys)")
        self.coord_sys = coord_sys
        self.position = position
        if position > coord_sys.dim:
            raise ValueError("Vector outside coordinate system dimensionality")

        is_commutative = True
        assumptions = {}
        assumptions['commutative'] = is_commutative
        self._assumptions = StdFactKB(assumptions)

    def __str__(self):
        return self.coord_sys.name + "." + self.name

    __repr__ = __str__

    def dot(self, other, coord_sys=None):
        return dot(self, other, coord_sys)

    def separate(self):
        # We just have a Vector - just return it
        coord_sys_dict = {self.coord_sys: [[], [self]]}
        return coord_sys_dict

    def express(self, coord_sys):
        __doc__ = express.__doc__
        return express(self, coord_sys)

    def expand(self):
        return self

    def factor(self):
        return self

    @property
    def _all_args(self):
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

    def _convert_to_rect(self, coord_sys, ind=1):
        """
        coord_sys: An instance of subclass of CoordSys class
        Takes a BaseScalar and converts it to combination of BaseScalars
        in coord_sys.
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

        rv = CoordSys._pos_to_rect(self.coord_sys.base_scalars,
                                   coord_conv[coord_sys.__class__])
        return rv[self.position - ind]


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
        import ipdb;ipdb.set_trace()
        vect = self.expand()

        # vect.args are either BaseVector or VectMul
        coord_sys_dict = {}
        coord_list = _all_coordinate_systems(vect)
        for c in coord_list:
            coord_sys_dict[c] = [[],[]]

        for arg in vect.args:
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
            coord_sys_dict[c][0] = list(set(v[0]))
            coord_sys_dict[c][1] = list(set(v[1]))

        return coord_sys_dict

    def expand(self):
        res = ZeroVector
        for arg in self.args:
            res = res + arg.expand()
        return res

    def factor(self):
        # self can contain BaseVector or VectMul
        factor_dict = {}
        for arg in self.args:
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
        vect = self.expand().factor()
        if not isinstance(vect, VectAdd):
            return vect.vector
        else:
            raise TypeError("Cannot separate vector from scalar for VectAdd")

    @property
    def scalar(self):
        # First we expand
        vect = self.expand()
        # Now we factor the base vectors
        vect = vect.factor()
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
            raise TypeError("The vector isn't in a single coordinate system")

        vect = self.expand()
        if not isinstance(vect, VectAdd):
            return vect.components

        r = [S.Zero]*self.max_dim
        for arg in vect.args:
            if isinstance(arg, BaseVector):
               # Component corresponding to it will be unity
               r[int(arg.position) + 1] = S.One
            elif isinstance(arg, VectMul):
               # The scalar part would be the component
               r[int(arg.vector.position) + 1] = arg.scalar

        return r

    @property
    def coord_sys(self):
        coord_list = _all_coordinate_systems(self)
        if not len(coord_sys) == 1:
            raise ValueError("More than one coordinate systems")
        else:
            return coord_list.pop()

    def express(self, coord_sys):
        __doc__ = express.__doc__
        return express(self, coord_sys)


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
        scalar = S.One
        for arg in self.args:
            if arg.is_Vector:
                continue
            scalar = scalar * arg

        if isinstance(scalar, BaseScalar):
            scalar_args = [scalar]
        else:
            scalar_args = scalar.args
        vect = self.vector

        if isinstance(vect, BaseVector):
            ret = ZeroVector
            for arg in scalar_args:
                ret = ret + VectMul(arg, vect)
            return ret

        if isinstance(vect, VectAdd):
            ret = ZeroVector
            exp_vect = vect.expand()
            # Now, exp_vect.args should be either BaseVector or VectMul
            for arg in exp_vect.args:
                ret = ret +  (arg * scalar).expand()
            return ret

    def factor(self):
        # This is of type scalar * BaseVector and scalar * VectAdd
        scalar = self.scalar
        vector = self.vector
        if isinstance(vector, BaseVector):
            # Just factor out the scalar and return
            return VectMul(scalar.factor(), vector)
        elif isinstance(vector, VectAdd):
            return VectMul(scalar.factor(), vector.factor())
        else:
            # Shouldn't happen
            raise ValueError

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
            v = v.factor()
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
        vect = self.expnad()
        if not isinstance(vect, VectMul):
            return vect.components

        # Since vect is just VectMul, so return the scalar part
        r = [S.One]*vect.coord_sys.dim
        r[int(vect.vector.position) + 1] = vect.vector.scalar
        return r

    def express(self, coord_sys):
        __doc__ = express.__doc__
        return express(self, coord_sys)


class ZeroVectorClass(Zero):
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

def _vect_add(one, other):
    # We are adding this method to check for cases involving S.Zero
    if ((one.is_Vector or one == S.Zero) and
        (other.is_Vector or other == S.Zero)):
        return VectAdd(one, other)
    else:
        raise TypeError("Cannot add a scalar with a vector")

def _vect_mul(one, other):
    #if one.is_Vector and other.is_Vector:
    #    raise TypeError("Cannot multiply two vectors")
    #if not one.is_Vector and not other.is_Vector:
    #    raise TypeError("At least one argument should be a vector")
    # Now we know that either one or other is a vector. Remaining is scalar
    #if one.is_Vector:
    #    return VectMul(other, one)
    #else:
    return VectMul(one, other)

def _vect_div(one, other):
    if one.is_Vector and other.is_Vector:
        raise TypeError("Cannot divide two vectors")
    if not one.is_Vector and not other.is_Vector:
        raise TypeError("At least one argument should be a vector")
    # Now we know that either one or other is a vector. Remaining is scalar
    if one.is_Vector:
        return VectMul(one, Pow(other, S.NegativeOne))
    else:
        raise TypeError("Cannot divide by vector")


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

    # First express everything in rect coordinates
    for arg in vect.args:
        # arg is necessarily a VectMul or a BaseVector
        assert isinstance(arg, BaseVector) or isinstance(arg, VectMul)
        scalar = arg.scalar
        vector = arg.vector

        # First process the scalar
        if isinstance(scalar, BaseScalar):
            # Because BaseScalar.args returns coord_sys as well
            c_rect = scalar.coord_sys._change_coord_sys(CoordSysRect, 'c_rect')
            subs_dict[scalar] = scalar.convert_to_rect(c_rect)

        else:
            for arg_scalar in scalar.args:
                if (isinstance(arg_scalar, BaseScalar) and not
                    subs_dict.has_key(arg_scalar)):
                    c_rect = arg_scalar.coord_sys._change_coord_sys(
                                                        CoordSysRect, 'c_rect')
                    subs_dict[arg_scalar] = BaseScalar.convert_to_rect(vect,
                                                                        c_rect)
        # Now process the vector
        # vector is necessarily BaseVector
        assert isinstance(vector, BaseVector)

        c_rect = vector.coord_sys._change_coord_sys(CoordSysRect, 'c_rect')
        subs_dict[vector] = vector.convert_to_rect(c_rect)

    # Now performig the substitution, we have changed all involved
    # variables into rectangular coordinates
    vect = vect.subs(subs_dict)
    subs_dict = {}

    vect = vect.expand()

    # First phase complete
    # Now, we find subs for base vectors and substitute it in
    for arg in vect.args:
        assert isinstance(arg, BaseVector) or isinstance(arg, VectMul)
        vector = arg.vector
        subs_dict[vector] = CoordSys._convert_base_vect_rect(vector,
                                                             coord_sys_t)
    # Now subs in the vectors
    vect = vect.subs(subs_dict)

    subs_dict = {}
    vect = vect.expand()

    # Now we need to find the subs_dict for each BaseScalar
    for arg in vect.args:
        assert isinstance(arg, BaseVector) or isinstance(arg, VectMul)
        scalar = arg.scalar

        # Process the scalar
        if isinstance(scalar, BaseVector):
            # Because BaseScalar.args returns coord_sys as well
            subs_dict[v] = CoordSys._convert_base_sclr_rect(scalar,
                                                            coord_sys_t)
        else:
            for arg_scalar in scalar.args:
                if (isinstance(arg_scalar, BaseScalar) and not
                    subs_dict.has_key(arg_scalar)):
                    subs_dict[arg_scalar] = \
                        CoordSys._convert_base_sclr_rect(vect, coord_sys_t)

    # Now performig the substitution, we have changed all involved
    # variables into rectangular coordinates
    vect = vect.subs(subs_dict)

    # Now we need to convert back from coord_sys_t to coord_sys
    subs_dict = {}
    vect = vect.expand()

    for arg in vect.args:
        # arg is necessarily a VectMul or a BaseVector
        assert isinstance(arg, BaseVector) or isinstance(arg, VectMul)
        scalar = arg.scalar
        vector = arg.vector

        # First process the scalar
        if isinstance(scalar, BaseScalar):
            # Because BaseScalar.args returns coord_sys as well
            subs_dict[scalar] = scalar.convert_to_rect(coord_sys)

        else:
            for arg_scalar in scalar.args:
                if (isinstance(arg_scalar) and not
                    subs_dict.has_key(arg_scalar)):
                    subs_dict[v] = arg_scalar.convert_to_rect(coord_sys)

        # Now process the vector
        # vector is necessarily BaseVector
        assert isinstance(vector, BaseVector)
        subs_dict[vector] = vector.convert_to_rect(coord_sys)

        # Now performig the substitution, we have changed all involved
        # variables into rectangular coordinates
        vect = vect.subs(subs_dict)
        subs_dict = {}

        return vect.factor()

"""
def express_in_global(vect, coord_sys=None):
"""
    #coord_sys: Instance of subclasses of CoordSys.
    #Express a vector in global coodinates.
"""
    if not coord_sys:
        coord_sys = CoordSysRect('rect_coord_sys')
    all_vectors_dict = vect.separate()
    coord_list = [c for c in all_vectors_dict.itervalues()]

    # Check whether given coord_sys is in global
    if coord_sys.position != S.Zero or not coord_sys._dcm_global.is_Identity:
        raise ValueError(str(coord_sys) + " is either psoitioned or oriented \
                         w.r.t. the global frame.")

    if (len(coord_list) == 1 and not vect.coord_sys.position
       and vect.coord_sys._dcm_global.is_Identity):
            # Vector is already in global coordnates
            # No position, no orientation
            # Return as is
            return vect
    else:
        ret_vector = ZeroVector

        for c, v in all_vectors_dict.iteritems():
            if c == coord_sys:
                ret_vector  = ret_vector + v
                continue

            exec "func = self._convert_to_" +  coord_conv[type(coord_sys)]

            # Tranforms we have to apply here:
            # 1. Coordinate system - rect, cyl or sph
            # 2. Orientation
            # 3. Position
            vect_i = CoordSys._pos_at(v, coord_sys)
            # vect_i is in rectangular coordinates -expressed with
            # changed position
            vect_i = CoordSys._orient_along(v, coord_sys)
            # Construct a vector, pass it to func
            vect_i = CoordSys._change_sys(v, coord_sys, func)
            ret_vector = ret_vector + vect_i
        return ret_vector.factor()
"""

ZeroVector = ZeroVectorClass()

coord_conv = {
    CoordSysRect: 'rect',
    CoordSysCyl: 'cyl',
    CoordSysSph: 'sph'
    }
