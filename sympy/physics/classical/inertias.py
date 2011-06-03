from kinematics import *

class Particle(Point):
    """
    This is for a Particle; it has position, velocity, acceleration, and
    a mass.  Forces can also be applied.  
    """

    # Needs a new init w/ mass attribute? or just mass attribute outside?
    # Or maybe a mass property?
    # Needs a way to store forces: Vector attribute? or getter/setter?

class RigidBody(Particle, ReferenceFrame):
    """
    This is for a RigidBody.  It extends Particle and ReferenceFrame, and 
    adds inertia tensor and torque. 
    """

    # Inertia terms: Property makes most sense
    # Also, torques; store same as Forces
