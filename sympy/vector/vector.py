from sympy.matrices import Matrix, eye
from sympy.core import (Basic, Expr, Dummy, Function, Symbol, symbols,
     sympify, diff, Pow, Mul, Add)
from sympy.core.numbers import Zero
from sympy.solvers import solve
from sympy.simplify import simplify, trigsimp

import copy

def base_scalars(names, coord_sys):
    base_scalar_list = []
    for i, name in enumerate(names.split(' ')):
        base_scalar_list.append(BaseScalar(name, coord_sys, i+1))
    return base_scalar_list

def vectors(names, coord_sys):
    vector_list = []
    for i, name in enumerate(names.split(' ')):
        vector_list.append(BaseScalar(name, coord_sys, i+1))
    return vector_list

def _dot_same(vect_a, vect_b):
    """
    The dot product of two vectors - both in same coordinate system.
    """
    a_comp = vect_a.components
    b_comp = vect_b.components
    length = len(a_comp)
    r = S.Zero
    for i in range(length):
        r += self_comp[i] * other_comp[i]
    return r

def dot(vect_a, vect_b, coord_sys=None):
    """
    Generalized dot product.
    """
    a_vectors = _separate_to_vectors(vect_a)
    b_vectors = _separate_to_vectors(vect_b)

    if not coord_sys and not len(_all_coord_system(a_vectors, b_vectors)) == 1:
        raise ValueError("Coordinate system not provided.")
    else:
        # coord_sys not given and all vectors have same coord_sys
        coord_sys = a_vectors[0].coord_sys

    r_vect_a = None
    for vector in a_vectors:
        r_vect_a += vector.express_in(coord_sys)
    r_vect_b = None
    for vector in b_vectors:
        r_vect_b += vector.express_in(coord_sys)
    # Now we have two vectors, each in the provided coord_sys
    return _dot_same(r_vect_a, r_vect_b)

def _separate_to_vector(vector):
    coord_dict = vect.separate()
    r = [vect for vect in coord_dict.itervalues()]

def _all_coordinate_sysetms(vector):
    vector = vector.expand()
    coord_list = []
    # all_args is a separate method that return only the vector args
    for arg in vector._all_args:
        if isinstance(arg, Vector):
            coord_list.append(arg.coord_sys)
        if isinstance(arg, VectMul):
            try:
                return [arg.coord_sys]
            except Exception as e:
                raise ValueError("Could not separate " + str(vector))


class BaseScalar(Symbol):
    """
    BaseScalar instances are used to express coordinate variables for field
    """
    def __new__(cls, name, coord_sys, position, **assumptions):
        if position not in ['1', '2', '3']:
            raise ValueError("Position of scalar not specified. \
                            See `position` in docstring of `BaseScalar`")
        if name is None:
            name = 'Dummy_' + str(Dummy._count)
        obj = Symbol.__new__(cls, name, **assumptions)
        obj.coord_sys = coord_sys
        return obj


class CoordSys(Basic):
    """
    Superclass for all CoordSys<type> classes. Not to be intialized
    directly.
    """
    def __init__(
            self, name, dim=S(3), position=None, position_coord='rect',
            orient_type=None, orient_amount=None, rot_order=None,
            parent=None
            ):
        if name is None:
            name = 'CoordSys_' + str(Dummy._count)

        self.name = name
        self.dim = as_int(dim)

        if position:
            for p in position:
                if isinstance(p, BaseScalar):
                    raise TypeError("Position variables cannot be BaseScalars")
            self.position_coord = position_coord
            if parent:
                pos_list = []
                exec "pos_func = _pos_to_" + self.position_coord
                parent_pos = pos_func(parent.position, parent.position_coord)
                for i in len(position):
                    pos_list.append(parent_pos[i] + position[i])
                self.position = pos_list
            else:
                self.position = position

        self._dcm_global = None
        self._dcm_parent = None

        if orient_type:
            self._check_orient_raise(orient_type, orient_amount, rot_order)
            if parent:
                self.parent = parent
                self._dcm_parent = self._dcm_parent(orient_type,
                                                    orient_amount, rot_order)

            self._dcm_global = self._dcm_global(orient_type, orient_amount)


    def _dcm_parent(self, orient_type, amounts, rot_order=''):
        orient_type = orient_type.capitalize()
        if orient_type == 'AXIS':
            if not rot_order == '':
                raise TypeError('Axis orientation takes no rotation order')
            if not (isinstance(amounts, (list, tuple)) & (len(amounts) == 2)):
                raise TypeError('Amounts are a list or tuple of length 2')
            theta = amounts[0]
            axis = amounts[1]
            if not axis.dt() == 0:
                raise ValueError('Axis cannot be time-varying')
            axis = axis.express(parent).normalize().as_mat()
            axis = axis.args[0][0]
            parent_orient = ((eye(3) - axis * axis.T) * cos(theta) +
                    Matrix([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]],
                        [-axis[1], axis[0], 0]]) * sin(theta) + axis * axis.T)
            return parent_orient

        elif rot_type == 'BODY':
            if not (len(amounts) == 3 & len(rot_order) == 3):
                raise TypeError('Body orientation takes 3 values & 3 orders')
            a1 = int(rot_order[0])
            a2 = int(rot_order[1])
            a3 = int(rot_order[2])
            parent_orient = (self._rot(a1, amounts[0]) *
                            self._rot(a2, amounts[1]) *
                            self._rot(a3, amounts[2]))
            parent_orient.applyfunc(lambda i: trigsimp(i, method='fu'))
            return parent_orient


    def _dcm_global(self, orient_type, orient_amount, rot_type):
        if self.parent:
            # Parent given therefore the given angle is wrt parent.
            # DCM(global<-parent)*DCM(parent<-self) == DCM(global <- self)
            r = parent.dcm().T * parent.dcm(self)
            r.applyfunc(lambda i: trigsimp(i, method='fu'))
            return r
        else:
            # Parent not set. So, directly initialize the dcm wrt global
            # Using sympy.physics.mechanics logic here.
            if rot_type == 'AXIS':
                if not rot_order == '':
                    raise TypeError('Axis orientation takes no rotation order')
                if not (isinstance(amounts, (list, tuple)) & (len(amounts) == 2)):
                    raise TypeError('Amounts are a list or tuple of length 2')
                theta = amounts[0]
                axis = amounts[1]
                axis = axis.in_global()
                if not axis.dt() == 0:
                    raise ValueError('Axis cannot be time-varying')
                axis = axis.normalize().as_mat()
                axis = axis.args[0][0]
                global_orient = ((eye(3) - axis * axis.T) * cos(theta) +
                        Matrix([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]],
                            [-axis[1], axis[0], 0]]) * sin(theta) + axis * axis.T)
                return global_orient
            if rot_type == 'BODY':
                if not (len(amounts) == 3 & len(rot_order) == 3):
                    raise TypeError('Body orientation takes 3 values & 3 orders')
                a1 = int(rot_order[0])
                a2 = int(rot_order[1])
                a3 = int(rot_order[2])
                global_orient = (self._rot(a1, amounts[0]) *
                                self._rot(a2, amounts[1]) *
                                self._rot(a3, amounts[2]))
                global_orient.applyfunc(lambda i: trigsimp(i, method='fu'))
                return parent_orient

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
        if orient_type == 'AXIS':
            if len(orient_amount) != 2:
                raise ValueError("orient_amount should be a \
                                list of lenght 2")
            if rot_order:
                raise ValueError("orient_axis works only with 'Body' \
                                type orientaion")

        elif orient_type == 'BODY':
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
        Create a new frame by positioning/orienting an exisitng frame.
        """
        newframe = copy.copy(self)

        if name == None:
            name = 'CoordSys_' + str(Dummy._count)
        newframe.name = name
        newframe._check_orient_raise(orient_type, orient_amount, rot_order)
        newframe.parent = self
        newframe._dcm_parent = newframe._dcm_parent(orient_type,
                                                    orient_amount, rot_order)
        newframe._dcm_global = newframe._dcm_global(orient_type, orient_amount)
        return newframe

    def posnew(self, position, position_coord='rect'):
        newframe = copy.copy(self)
        newframe.position_coord = position_coord
        pos_list = []
        exec "pos_func = _pos_to_" + self.position_coord
        parent_pos = pos_func(parent.position, parent.position_coord)
        for i in len(position):
            pos_list.append(parent_pos[i] + position[i])
        self.position = pos_list

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
            return (x, y, z)

    @staticmethod
    def _pos_to_sph(coord, con_from):
        if conv_from == 'rect':
            r = sqrt(coord[0]**2 + coord[1]**2 + coord[2]**2)
            theta = atan(sqrt(coord[0]**2 + coord[1]**2)/coord[2])
            phi = atan(coord[1]/coord[0])
            return (r, theta, phi)

        if conv_from == 'cyl':
            r = sqrt(coord[0]**2 + coord[2]**2)
            theta = atan(coord[0]/coord[2])
            phi = coord[1]
            return (t, theta, phi)

class CoordSysRect(CoordSys):
    """
    The rectangular coordinate system.
    """
    def __init__(self, *args, **kwargs):
        super(CoordSysRect, self).__init__(*args, **kwargs)

        self.h_list = [S.One]*self.dim

    @staticmethod
    def _convert_to_sph(vector, base_scalrs, base_vectors):
        """
        vector : a VectAdd
        base_scalrs : list or tuple of base scalars
        base_vectors : list or tuple of base vectos
        returns an object with is_Vector == True
        """
        Ax, Ay, Az = vector.components
        x, y, z = vector.base_scalars
        r, theta, phi = base_scalars
        subs_dict = {
                        x : r*sin(theta)*cos(phi),
                        y : r*sin(tehta)*sin(phi),
                        z : r*cos(theta)
                    }
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
            res.append(VectMul(vect, base_vectors[i]))
        return VectAdd(*res)

    @staticmethod
    def _convert_to_cyl(vector, base_scalars, base_vectors):
        """
        vector : a VectAdd
        base_scalrs : list or tuple of base scalars
        base_vectors : list or tuple of base vectos
        returns an object with is_Vector == True
        """
        Ax, Ay, Az = vector.components
        x, y, z = vector.base_scalars
        rho, phi, _z = base_scalars
        subs_dict = {
                        x : rho*cos(phi),
                        y : rho*sin(phi)
                        z : _z
                    }
        Ax = Ax.subs(subs_dict)
        Ay = Ay.subs(subs_dict)
        Az = Az.subs(subs_dict)
        mat =  Matrix([
                        [ cos(phi), sin(phi), 0],
                        [-sin(phi), cos(phi), 0],
                        [        0,        0, 1]
                     ])
        r = mat*Matrix([ [Ax], [Ay], [Az] ])
        res = []
        for i, vect in enumerate(r._mat):
            res.append(VectMul(vect, base_vectors[i]))
        return VectAdd(*res)

class CoordSysSph(CoordSys):
    """
    The spherical polar coordinate system.
    """
    def __init__(self, *args, **kwargs):
        super(CoordSysSph, self).__init__(*args, **kwargs)

        self.one, self.two, self.three = symbols('_1 _2 _3')
        self.h_list = [
                        S.One,
                        self.one,
                        self.one*sin(self.two)
                      ]

    @staticmethod
    def _convert_to_rect(vector, base_scalrs, base_vectors):
        """
        vector : a VectAdd
        base_scalrs : list or tuple of base scalars
        base_vectors : list or tuple of base vectos
        returns an object with is_Vector == True
        """
        Ar, At, Ap = vector.components
        r, theta, phi = vector.base_scalars
        x, y, z = base_scalars
        subs_dict = {
                        r : sqrt(x**2 + y**2 + z**2),
                        theta : atan(sqrt(x**2 + y**2)/z),
                        phi : atan(y/x)
                    }
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
            res.append(VectMul(vect, base_vectors[i]))
        return VectAdd(*res)

    @staticmethod
    def _convert_to_cyl(vector, base_scalrs, base_vectors):
        """
        vector : a VectAdd
        base_scalrs : list or tuple of base scalars
        base_vectors : list or tuple of base vectos
        returns an object with is_Vector == True
        """
        Ar, At, Ap = vector.components
        r, theta, phi = vector.base_scalars
        rho, _phi, z = base_scalars
        subs_dict = {
                        r : sqrt(rho**2 + z**2),
                        theta : atan(rho/z),
                        phi : _phi

                    }
        Ar = Ar..subs(subs_dict)
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
            res.append(VectMul(vect, base_vectors[i]))
        return VectAdd(*res)


class CoordSysCyl(CoordSys):
    """
    The spherical polar coordinate system.
    """
    def __init__(self, *args, **kwargs):
        super(CoordSysCyl, self).__init__(*args, **kwargs)

        self.one, self.two, self.three = symbols('_1 _2 _3')
        self.h_list = [
                        S.One,
                        self.one,
                        S.One
                      ]

    @staticmethod
    def _convert_to_rect(vector, base_scalrs, base_vectors):
        """
        vector : a VectAdd
        base_scalrs : list or tuple of base scalars
        base_vectors : list or tuple of base vectos
        returns an object with is_Vector == True
        """
        Ar, Ap, Az = vector.components
        rho, phi, z = vector.base_scalars
        x, y, _z = base_scalars
        subs_dict = {
                        rho : sqrt(x**2 + y**2),
                        phi : atan(y/x),
                        z   : _z
                    }
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
    def _convert_to_sph(vector, base_scalrs, base_vectors):
        """
        vector : a VectAdd
        base_scalrs : list or tuple of base scalars
        base_vectors : list or tuple of base vectos
        returns an object with is_Vector == True
        """
        Ar, Ap, Az = vector.components
        rho, phi, z = vector.base_scalars
        r, theta, _phi = base_scalars
        subs_dict = {
                        rho : r*sin(theta),
                        phi : _phi,
                        z   : r*cos(theta)
                    }
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
            res.append(VectMul(vect, base_vectors[i]))
        return VectAdd(*res)


class Vector(Basic):
    """
    Class instances will represent base vectors
    """
    is_Vector = True

    def __init__(self, name, coord_sys, position):
        if not name:
            name = 'Vector_' + str(Dummy._count)

        self.name = name
        if not isinstance(coord_sys, CoordSys):
            raise TypeError("coord_sys must be isinstance(CoordSys)")
        self.coord_sys = coord_sys

        self.position = position

    def __add__(self, other):
        return VectAdd(self, other)

    def __mul__(self, other):
        return VectMul(self, other)

    def separate(self):
        # We just have a Vector - just return it
        coord_sys_dict = {self.coord_sys: coord_sys}
        ret coord_sys_dict

    @property
    def _all_args(self):
        return self

    @property
    def components(self):
        # Since it is a base vector, so return a list
        # of len == dim(coord_sys) with the pos element
        # as unity.
        r = [S.One]*self.coord_sys.dim


class VectAdd(Add):
    """
    Container to hold added Vectors/VectMuls
    """
    is_Vector = True
    def __new__(cls, *args, **options):
        for arg in args:
            if not arg.is_Vector:
                except TypeError(str(arg) + " is not a vector.")
        obj = super(VectAdd, cls).__new__(cls, *args, **options)
        return obj

    def __add__(self, other):
        return VectAdd(self, other)

    def __mul__(self, other):
        return VectMul(self, other)

    def dot(self, other, coord_sys=None):
        return dot(self, other, coord_sys)

    def separate(self):
        # Flatten the VectMul so that there are no nested VectAdds
        vect = self.expand().factor()
        args = vect.args

        coord_sys_dict = {}
        for arg in args:
            # arg is either Vector or VectMul

            if coord_sys_dict.has_key(v.coord_sys):
                if isinstance(arg, Vector):
                    coord_sys_list[arg.coord_sys].append(arg)
                elif isinstance(arg, VectMul):
                    try:
                        coord_sys_list[arg.coord_sys].append(arg)
                    except Exception as e:
                        raise ValueError("Could not expand " + str(vect))

            else:
                coord_sys_dict[v.coord_sys] = []

        for c, v in coord_sys_list.iteritems():
            coord_sys_dict[c] = VectAdd(*v)

        return coord_sys_dict

    @property
    def _all_args(self):
        return list(self.args)


class VectMul(Add):
    """
    Container to hold added Vectors/VectMuls
    """
    is_Vector = True
    def __new__(cls, *args, **options):
        counter = 0
        for arg in args:
            if arg.is_Vector:
                counter += 1
        if counter > 1:
            raise TypeError("Cannot multiply two or more vectors.")
        obj = super(VectMul, cls).__new__(cls, *args, **options)
        return obj

    def __add__(self, other):
        return VectAdd(self, other)

    def __mul__(self, other):
        return VectMul(self, other)

    def separate(self):
        # First we flatten things out - so that there are no nested VectAdd
        vect = self.expand().factor()
        if isinstance(vect, VectAdd) or isinstance(vect, Vector):
            return vect.separate()

        # Now we are sure that vect is just VectMul - no nesting
        coord_list_dict = {vect.vector.coord_sys: vect}
        return coord_list_dict

    @property
    def coord_sys(self):
        """
        If the object is only defined in one coordinate system
        then return that cooordinate system.
        Cases:
        1. scalar*Vector
        2. scalar*VectAdd
        """
        vect = self.vector
        if isinstance(vect, Vector):
            return vect.coord_sys

        elif isinstance(vect, Vector):
            coord_list = _all_coordinate_systems(vect)
            if len(coord_list) == 1:
                return coord_list[0]

        else:
            raise ValueError("Cannot recognize " + str(vector))

    @property
    def vector(self):
        """
        Return the vector contained in VectMul.
        Returns a Vector or a VectAdd
        """
        for arg in self.args:
            if isinstance(arg, Vector) or isinstance(arg, VectAdd):
                return arg
        # If we haven't found a vector, raise error
        raise TypeError(str(self) + " doesn't have a Vector/VectAdd \
                        its args")

    @property
    def _all_args(self):
        return [self]
