from sympy import *
from copy import deepcopy

class Vector:
    """
    This is the class used to define vectors.  It along with reference frame are the building blocks of pydy.  
    Class attributes include: 
        subscript_indices - a 3 character string used for printing
    """

    subscript_indices = "xyz"

    def __init__(self,mat,frame):
        """
        This is the constructor for the Vector class.  
        It should only be used in construction of the basis vectors,
        which is part of the ReferenceFrame construction.  
        It takes in a SymPy matrix and a ReferenceFrame.  
        """
        self.args=[[mat,frame]]

def __str__(self):
        """
        Printing method.  Uses Vector Attribute subscript_indices to choose how
        to show basis vector indices.
        """
        ar = self.args # just to shorten things
        ol = [] # output list, to be concatenated to a string
        for i in range(len(ar)):
            for j in 0,1,2:
                if ar[i][0][j] == 1: # if the coef of the basis vector is 1, we skip the 1
                    if len(ol) != 0: 
                        ol.append(' + ')
                    ol.append( ar[i][1].name.lower() + self.subscript_indices[j] + '>' )
                elif ar[i][0][j] == -1: # if the coef of the basis vector is -1, we skip the 1
                    if len(ol) != 0:
                        ol.append(' ')
                    ol.append( '- ' + ar[i][1].name.lower() + self.subscript_indices[j] + '>' )
                elif ar[i][0][j] != 0: 
                    # If the coefficient of the basis vector is not 1 or -1, 
                    # we wrap it in parentheses, for readability.
                    if len(ol) != 0:
                        ol.append(' + ')
                    ol.append('(' +  `ar[i][0][j]` + ')*' + ar[i][1].name.lower() + self.subscript_indices[j] + '>' )
        return ''.join(ol)

    def __repr__(self):
        """
        Wraps __str__
        """
        return self.__str__()

    def __add__(self, other):
        """
        The add operator for Vector. 
        It checks that other is a Vector, otherwise it throws an error.
        """
        if not(isinstance(other,Vector)): # Rejects adding a scalar to a vector
            raise TypeError('You can only add two Vectors')
        self2 = self.copy()
        other2 = other.copy()
        for i in range(len(self2.args)):
            for j in range(len(other2.args)):
                if self2.args[i][1] == other2.args[j][1]:
                    self2.args[i][0] += other2.args[j][0]
                    other2.args.remove(j)
        self2.args += other.args[j]
        
        return self2

    def __and__(self, other):
        """
        Dot product of two vectors.  
        """
        out = 0
        for i in range(len(self.args)):
            for j in range(len(other.args)):
                out += ((other.args[j][0].T)\
                        *(self.args[i][1].dcm(other.args[j][1]))\
                        *(self.args[i][0]))[0]
        return out

    def __div__(self, other):
        """
        This uses mul and inputs self and 1 divided by other.  
        """
        return self.__mul__(1/other)

    def __mul__(self, other):
        """
        Multiplies the Vector by a scalar. 
        Throws an error if another Vector is entered.  
        """
        if isinstance(other, Vector): # Rejects scalar multiplication of two vectors
            raise TypeError('Why u try to mul vecs?')
        self2 = self.copy()
        for i in range(len(self2.args)):
            self2.args[i][0] *= other
        return self2

    def __rmul__(self, other):
        """
        This wraps mul. 
        """
        return self.__mul__(other)

    def __sub__(self, other):
        """
        The subraction operator. 
        Reuses add and multiplication operations.  
        """
        return self.__add__(other*-1)

    def copy(self):
        newvec = Vector(Matrix([0,0,0]),None)
        for i in range(len(self.args)-1):
            newvec.args.append([Matrix([0,0,0]),None])
        for i in range(len(self.args)):
            newvec.args[i][0][0:] = self.args[i][0][0:]
            newvec.args[i][1] = self.args[i][1]
        return newvec

    def dot(self, other):
        """
        Wraps around & operator, which is the designated operator for dot.
        """
        return self&other

    def express(self, otherframe):
        """
        Returns the vector, expressed in the other frame.
        Uses a DCM from each basis vector triplet's current frame to the other
        frame.
        Takes in a frame.
        """
        assert isinstance(otherframe, ReferenceFrame), 'Needs a frame to express in'
        out = Vector(Matrix([0,0,0], ReferenceFrame))
        for i in range(len(self.args)):
            if self.args[i][1] == otherframe:
                out.args[i][0] += self.args[i][0]
            else:
                out.args[i][0] += self.args[i][1].dcm(otherframe) *\
                self.args[i][0]
        return out


class ReferenceFrame:
    """
    ReferenceFrame is a reference frame.
    It will store its basis vectors as attributes, 
    and orientation information to a parent frame,
    or it will be at the top of a tree.
    """

    def __init__(self, name=''):
        """
        Constructor for ReferenceFrame
        """
        self.name = name
        self.parent = None
        self.x = Vector(Matrix([1,0,0]), self)
        self.y = Vector(Matrix([0,1,0]), self)
        self.z = Vector(Matrix([0,0,1]), self)

    def __str__(self):
        """
        Returns the name of the frame.
        """
        return self.name
    
    def __repr__(self):
        """
        Wraps __str__
        """
        return self.name

    def dcm(self,other):
        """
        This will return the direction cosine matrix from self to other;
        other's basis = DCM * self's basis is the definition.
        All it takes in is the other frame.  
        nxyz = N.dcm(B)*bxyz
        """
        assert isinstance(other, ReferenceFrame), 'You have to use 2\
        reference frames'
        leg1 = [self]
        ptr = self
        while ptr.parent != None:
            ptr = ptr.parent
            leg1.append(ptr)
        leg2 = [other]
        ptr = other
        while ptr.parent != None:
            ptr = ptr.parent
            leg2.append(ptr)
        try:
            commonframe = (set(leg1) & set(leg2)).pop()
        except:
            raise ValueError('No Common Parent')
        # form DCM from self to first common frame
        leg1 = eye(3)
        ptr = self
        while ptr != commonframe:
            leg1 *= ptr.parent_orient
            ptr = ptr.parent
        # form DCM from other to first common frame
        leg2 = eye(3)
        ptr = other
        while ptr != commonframe:
            leg2 *= ptr.parent_orient
            ptr = ptr.parent
        return leg1.T * leg2

    def orientnew(self, newname, rot_type, amounts, rot_order):
        newframe = ReferenceFrame(newname)
        newframe.orient(self, rot_type, amounts, rot_order)
        return newframe

    def orient(self, parent, rot_type, amounts, rot_order):
        """
        This function will be used to define the orientation of a
        ReferenceFrame relative to a parent.  It takes in the parent frame,
        type of rotation, amount(s) of rotation, and if applicable, order of
        rotations.  
        The format for this needs to be spelled out and made explicit, with
        example.  
        """

        def _rot(axis, angle): 
            """
            Returns direction cosine matrix for simple 1,2,3 rotations
            """
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

        self.parent = parent
        approved_orders = ('123','231','312','132',\
                '213','321','121','131','212','232',\
                '313','323','1','2','3')
        rot_order = str(rot_order).upper() # Now we need to make sure XYZ = 123
        rot_type  = rot_type.upper()
        rot_order = [i.replace('X', '1') for i in rot_order]
        rot_order = [i.replace('Y', '2') for i in rot_order]
        rot_order = [i.replace('Z', '3') for i in rot_order]
        rot_order = ''.join(rot_order)
        assert rot_order in approved_orders, 'Not approved order'

        if rot_type == 'AXIS':
            raise NotImplementedError('Axis rotation not yet implemented')
        elif rot_type == 'EULER':
            assert len(amounts)==4, 'Euler orietation requires 4 values'
            assert rot_order=='', 'Euler orientation take no rotation order'
            # TODO need to finish this up...
        elif rot_type == 'BODY':
            assert len(amounts)==3, 'Body orientation requires 3 values'
            assert len(rot_order)==3, 'Body orientation requires 3 orders'
            a1 = int(rot_order[0])
            a2 = int(rot_order[1])
            a3 = int(rot_order[2])
            self.parent_orient = _rot(a1, amounts[0]) * _rot(a2, amounts[1])\
                    * _rot(a3, amounts[2])
        elif rot_type == 'SPACE':
            assert len(amounts)==3, 'Space orientation requires 3 values'
            assert len(rot_order)==3, 'Space orientation requires 3 orders'
            a1 = int(rot_order[0])
            a2 = int(rot_order[1])
            a3 = int(rot_order[2])
            self.parent_orient = _rot(a3, amounts[2]) * _rot(a2, amounts[1])\
                    * _rot(a1, amounts[0])
        elif rot_type == 'SIMPLE':
            assert not(isinstance(amounts,(list,tuple))), 'Simple takes in a\
            single value for amount'
            assert not(isinstance(rot_order,(list,tuple))), 'Simple takes in a\
            single value for rotation order'
            a = int(rot_order)
            self.parent_orient = _rot(a, amounts)
        else:
            raise NotImplementedError('That is not an implemented rotation')


