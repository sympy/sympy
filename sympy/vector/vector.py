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

def dot(vect_a, vect_b):
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
            self.position_coord = position_coord
            if parent:
                pos_list = []
                exec "pos_func = _pos_to_" + self.position_coord
                for i, p in enumerate(position):
                    pos_list.append(pos_func(parent.position, i+1) + p)
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
        for i, p in enumerate(position):
            pos_list.append(pos_func(parent.position) + p)
        self.position = pos_list

    @staticmethod
    def _pos_to_rect(coord, pos):



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
    def _convert_to_cyl(vector, base_scalrs, base_vectors):
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
    def _convert_to_cyl(vector, base_scalrs, base_vectors):
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

