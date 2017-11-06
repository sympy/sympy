from sympy.physics.mechanics import Point, ReferenceFrame, dynamicsymbols,\
Body, inertia, Vector, mprint
from sympy import symbols, pretty_print
from sympy.matrices import Matrix, zeros, eye

# The purpose of this script is to derive the equations of motion of a 1-link
# pendulum using the spatial operator framework. Much of this example
# is developed using Abhi's text and his paper in IJRR with KKD and GR.

# Section 1.1
def RotationMatrix(toFrameF, fromFrameG):
    return toFrameF.dcm(fromFrameG)

def HomogeneousTransform( toFrameF, fromFrameG, posVecFromFtoG ):
    RotMat = RotationMatrix(toFrameF, fromFrameG)
    posVecCol = posVecFromFtoG.args[0][0]
    return Matrix([ [RotMat, posVecCol], [zeros(1,3), Matrix([1])] ])

# Create two frames F and G
q1 = dynamicsymbols('q1')
F = ReferenceFrame('F')
G = F.orientnew('G', 'Axis', [q1, F.z])
# Now create two points to serve as origins of each frame
f = Point('f')
l_fg = symbols('l_fg')
g = f.locatenew('g', l_fg * F.x)
# And create a third point o using g
p_go = symbols('p_go')
o = g.locatenew('o', p_go * G.y)

# TEST: HomogeneousTransform and Rotation Matrix
pgoCol4 = Matrix( [[o.pos_from(g).args[0][0]],[ Matrix([1])]] )
lFtoG = g.pos_from(f)
HomTranGtoF = HomogeneousTransform(F, G, lFtoG)
o.pos_from(f).express(F).args[0][0]==(HomTranGtoF*pgoCol4)[0]

# Section 1.3 Spatial Vectors
""" Create the SpatialVector class. This serves as the superclass for
SpatialVelocity and, later on, will be used to develop a SpatialForce
"""
class SpatialVector(object):
    def __init__(self):
        self.w = Matrix([0, 0, 0])
        self.v = Matrix([0, 0, 0])
        self.frame = None

    def __call__(self):
        return ( Matrix( [ [self.w], [self.v] ] ), self.frame)

#    def tilde(self):
#        z = zeros(3,3)
#        return Matrix( [self.w.tilde(), z], [self.v.tilde(), self.w.tilde()] )

class SpatialVelocity(SpatialVector):
    """ This defines a class to generate a spatial velocity at a point on
    a body.

    """

    def __init__(self, point=None, bodyPointIsOn=None, observerFrame=None):
        if point == None and bodyPointIsOn == None and observerFrame == None:
            SpatialVector.__init__(self)
        else:
            self.point = point
            self.body = bodyPointIsOn
            bodyFrame = self.body.frame
            if bodyFrame.ang_vel_in(observerFrame).express(bodyFrame).args != []:
                self.w = bodyFrame.ang_vel_in(observerFrame).express(bodyFrame).args[0][0]
            if point.vel(observerFrame).express(bodyFrame).args != []:
                self.v = point.vel(observerFrame).express(bodyPointIsOn.frame).args[0][0]
            if bodyPointIsOn != None:
                self.frame = bodyPointIsOn.frame
            else:
                self.frame = None

    def express(self, inFrame):
        fromFrame = self.__call__()[1]
        zeroPosVec = Vector(0)
        return phiStar(inFrame, fromFrame, zeroPosVec)*self.__call__()[0]

# TEST: SpatialVelocity class
I = ReferenceFrame( 'I' )
l = symbols( 'l' )
q = dynamicsymbols( 'q' )
qd = dynamicsymbols('q', 1)
B = I.orientnew('B', 'Axis', [q, I.z] )
# o is the origin of the inertial frame
# x will eventually be a point on a body
x = o.locatenew( 'x', I.x + I.y)
x.set_vel( I, 10 * B.y )
y = x.locatenew( 'y', l * B.x )

# Now we will create a body
m = symbols('m')
mass_center = x.locatenew('mass_center',l/2*B.x + l/3*B.y + l/4*B.z)
ixx, iyy, izz = symbols('ixx iyy izz')
inertia = inertia(B, ixx, iyy, izz)
body = Body('body', mass_center, m, B, inertia)
# Constructing the spatial velocity of x will make it a point of the Body
# object 'body'
Vx = SpatialVelocity( x, body, I)
Vx()

# Section 1.4
# TEST: Computing the spatial velocity of another point 'y' also on the
# same rigid body B using the Phi operator, also known as the rigid body
# transformation matrix (Robot and Multibody Dynamics by Abhi Jain).

def spatialVelocity2PtThm(ofY, usingX, onBody):
    """
    This is the 2 point theorem to compute SPATIAL velocity of the point y
    on a rigid body B using the spatial velocity of point x which is also on
    the smae rigid body B.
        
    Personal Note: I can see this being implemented as a method of the
    Point Class. If so, the implementation would compute the spatial
    velocity of x as x.spatialVelocity2PtThm(y, B)

    """
    if not isinstance(ofY, Point):
        raise TypeError('Supply Point for first argument')
    if not isinstance(usingX, Point):
        raise TypeError('Supply Point for second argument')
    Vy = SpatialVelocity()
    Vx = SpatialVelocity(usingX, onBody, I)
    posVecXtoY = ofY.pos_from(usingX)
    phi_Star = phiStar(onBody.frame, onBody.frame, posVecXtoY )
    spVel_y = phi_Star * Vx()[0]
    Vy.w = spVel_y[0]
    Vy.v = spVel_y[1]
    Vy.frame = Vx.frame
    return Vy

def phi(x, y, posVecXY):
    """
    Computes the rigid body transformation matrix given on page 11, section 1.4.
    The phi matrix is a 6X6 Matrix. It can however also be organized as a
    matrix of four 3X3 matrices:

        Phi = [I tilde_p(x,y); Z I]
    
    where I is the 3X3 identity matrix, Z is the 3X3 zero matrix, and
    tilde_p(x,y) converts the position vector from x to y into a 3X3 matrix.
    This conersion is through an operator called the tilde operator.
    """
    if not isinstance(x, ReferenceFrame) or not isinstance(y,\
                                                                ReferenceFrame):
        raise TypeError('Please provide ReferenceFrame objects as the first'\
                        ' and second arguments')
    else:
        RotMat = RotationMatrix(x, y) #should be a Matrix object
    if not isinstance(posVecXY, Vector):
        raise TypeError('Please supply a Vector object for the last argument.')
    else:
        l_XY = posVecXY
    Phi11 = eye(3) * RotMat
    Phi22 = Phi11 
    Phi12 = l_XY.express(x).tilde() * RotMat
    Phi21 = zeros(3)
    return Matrix( [ [Phi11, Phi12], [ Phi21, Phi22 ] ] )

def phiStar(x, y, posVecXY):
    if not isinstance(x, ReferenceFrame) or not isinstance(y,\
                                                                ReferenceFrame):
        raise TypeError('Please provide ReferenceFrame objects as the first'\
                        ' and second arguments')
    else:
        RotMat = RotationMatrix(x, y) #should be a Matrix object
    if not isinstance(posVecXY, Vector):
        raise TypeError('Please supply a Vector object for the last argument.')
    else:
        l_XY = posVecXY
    M11 = eye(3) * RotMat
    M22 = M11 
    M21 = -l_XY.express(x).tilde() * RotMat
    M12 = zeros(3)
    return Matrix( [ [M11, M12], [ M21, M22 ] ] )

# TEST: Compute the spatial velocity of y using the spatialVelocty2PtThm which
# is on the same body as x
Vy = spatialVelocity2PtThm(y, x, body)

## THUS ENDETH SINGLE BODY KINEMATICS

# Chapter 2
class SpatialInertia(object):
    def __init__(self, bodyName, pointName):
        self._point = pointName
        self._body = bodyName
        self._mass = bodyName.mass
        #if pointName.pos_from(self._body.masscenter) == []:
        #    self._posvec =  
        #else:
        self._posvec = pointName.pos_from(self._body.masscenter)
        #if self._point == bodyName.masscenter:
        #    self.inertia = bodyname.central_inertia
        #else:
        self.inertia = bodyName.central_inertia.to_matrix(self._body.frame) - self._mass *\
                self._posvec.tilde() * self._posvec.tilde() #Equation 2.11 on page 19
        self.spI = self._compute()
        
    def __call__(self, pointName):
        """ Equation 2.12, page 20- Jain textbook
        This is the parallel axis theorem to compute spatial inertia.
        Returns a Matrix object."""
        posVec = pointName.pos_from(self._point)
        bodyFrame = self._body.frame
        return phi(bodyFrame, bodyFrame,\
                   posVec)*self.spI*phiStar(bodyFrame, bodyFrame, posVec)

    def _compute(self):
        """Equation 2.13, page 20- Jain textbook
        This is the parallel axis theorem to compute the moment of inertia.
        Returns a Matrix object."""
        offDiagonalMatrix = self._mass * self._posvec.tilde()
        massMatrixify = self._mass * eye( 3 )
        return Matrix( [ [self.inertia, offDiagonalMatrix ], [\
            -offDiagonalMatrix, massMatrixify ] ])

# This tests the SpatialInertia class.
M = SpatialInertia(body, mass_center)
pretty_print( M(mass_center) )
pretty_print( M(x))
#

# IRRELEVANT TO THE PULL REQUEST
# Section 1.5 Spatial Force
#class SpatialForce(SpatialVector):
#    """
#    The spatial force is defined for a frame attached to a point x.
#    """
#
#    def __init__(self, pointName=None, bodyName=None):
#        if pointName == None and bodyName == None:
#            SpatialVector.__init__(self)
#        else:
#
#    def __call__(pointName):
#        return Matrix([self.w, self.v]) 
#
##torques = body.applyTorque()
##forces = body.applyForce()
##Together they can be accessed by:
#body.loads
#
## how do we compute the equations of motion?
## Jain says there are many didefine a spatial momentum
