__all__ = ['RigidBody']

from sympy import sympify
from sympy.physics.classical.point import Point
from sympy.physics.classical.essential import ReferenceFrame
from sympy.physics.classical.inertiadyadic import InertiaDyadic

class RigidBody(object):
    """An idealized rigid body.

    >>> B = RigidBody()

    """

    def __init__(self):
        self._frame = None
        self._mc = None
        self._mass = None
        self._inertia = None

    @property
    def frame(self):
        return self._frame

    @frame.setter
    def frame(self, F):
        if not isinstance(F, ReferenceFrame):
            raise TypeError("RigdBody frame must be a ReferenceFrame object.")
        self._frame = F

    @property
    def mc(self):
        return self._mc

    @mc.setter
    def mc(self, p):
        if not isinstance(p, Point):
            raise TypeError("RigidBody mass center must be a Point object")
        self._mc = p

    @property
    def mass(self):
        return self._mass

    @mass.setter
    def mass(self, m):
        self._mass = sympify(m)

    @property
    def inertia(self):
        return self._inertia

    @inertia.setter
    def inertia(self, I):
        if not isinstance(I, InertiaDyadic):
            raise TypeError("RigidBody central inertia must be a InertiaDyadic "
                            "object.")
        self._inertia = I

if __name__ == "__main__":
    import doctest
    doctest.testmod()
