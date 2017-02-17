from sympy.physics.mechanics import ReferenceFrame, dynamicsymbols
from sympy import pi

N = ReferenceFrame('N')
A = ReferenceFrame('A')
q = dynamicsymbols('q')

#currently orient aligns A.z with N.z
#A.orient(N,'Axis', [q, N.z])

#but maybe what we want is to actually align A.x with N.z and then rotate
a1 = N.orientnew('a1', 'Axis', [pi/2, N.z]) #makes a1.x = N.y
a2 = a1.orientnew('a2', 'Axis', [pi/2, a1.x]) #makes a2.y = a1.z
Af = a2.orientnew('Af', 'Axis', [q,a2.y]) #simple rotation where A.z and N.z
                                          # aren't aligned.

def alignyz(childframe,parentframe):
    """Aligns the childframe's y-axis with parentframe's z-axis..

    Parameters
    ==========

    childframe : Frame to be rotated

    parentframe : Frame that is fixed
    """
    if not isinstance(childframe,ReferenceFrame):
        raise TypeError('You suck, child')
    if not isinstance(parentframe,ReferenceFrame):
        raise TypeError('Parent is sucky')
    # The funciton will automatically generate two new frames, b1, b2 and
    # perform two orthogonal rotations
    else:
        b1 = parentframe.orientnew('b1', 'Axis', [pi/2, parentframe.z]) #makes b1.x = N.y
        b2 = b1.orientnew('b2', 'Axis', [pi/2, b1.x]) #makes b2.y = b1.z
        childframe.orient(b2, 'Axis', [q, b2.y]) #simple rotation where A.y = N.z
    return childframe.dcm(parentframe)

def align(childframedirection,parentframedirection, q=0):
    """ Supposed to be a general align function.

    Parameters
    ==========

    rotframeaxis : Axis to be rotated

    parentframeaxis : Axis that rotframeaxis should finally be aligned with

    """
    if not isinstance(childframedirection.__dict__['args'][0][1],
                      ReferenceFrame):
        raise TypeError('You suck, child')
    if not isinstance(parentframedirection.__dict__['args'][0][1],
                      ReferenceFrame):
        raise TypeError('You are a sucky parent')
    childframe = childframedirection.__dict__['args'][0][1]
    parentframe = parentframedirection.__dict__['args'][0][1]
    childframe.orient(parentframe, 'Body', [pi/2, pi/2, q], 'zxy')
    return childframedirection == childframe.y
