__all__ = ['cross',
           'dot',
           'dynamicsymbols',
           'express',
           'outer',
           'inertia']

from sympy.physics.mechanics.essential import Vector, Dyad, ReferenceFrame
from sympy.physics.mechanics.dynamicsymbol import DynamicSymbol
from sympy import sympify, symbols

def cross(vec1, vec2):
    """Cross product convenience wrapper for Vector.cross(): \n"""
    if not isinstance(vec1, (Vector, Dyad)):
        raise TypeError('Cross product is between two vectors')
    return vec1 ^ vec2
cross.__doc__ += Vector.cross.__doc__

def dot(vec1, vec2):
    """Dot product convenience wrapper for Vector.dot(): \n"""
    if not isinstance(vec1, (Vector, Dyad)):
        raise TypeError('Dot product is between two vectors')
    return vec1 & vec2
dot.__doc__ += Vector.dot.__doc__

def dynamicsymbols(names):
    """Wraps sympy.symbols to use DynamicSymbol. """
    return symbols(names, cls=DynamicSymbol)

def express(vec, frame, frame2=None):
    """Express convenience wrapper for Vector.express(): \n"""
    if not isinstance(vec, (Vector, Dyad)):
        raise TypeError('Can only express Vectors')
    if isinstance(vec, Vector):
        return vec.express(frame)
    else:
        return vec.express(frame, frame2)

express.__doc__ += Vector.express.__doc__

def outer(vec1, vec2):
    """Outer prodcut convenience wrapper for Vector.outer():\n"""
    if not isinstance(vec1, Vector):
        raise TypeError('Outer product is between two Vectors')
    return vec1 | vec2
outer.__doc__ += Vector.express.__doc__

def inertia(frame, ixx, iyy, izz, ixy=0, iyz=0, izx=0):
    """Simple way to create inertia Dyad object.

    If you don't know what a Dyad is, just treat this like the inertia tensor.
    Then, do the easy thing and define it in a body-fixed frame.

    Parameters
    ==========
    frame : ReferenceFrame
        The frame the inertia is defined in
    ixx : Sympifyable
        the xx element in the inertia dyad
    iyy : Sympifyable
        the yy element in the inertia dyad
    izz : Sympifyable
        the zz element in the inertia dyad
    ixy : Sympifyable
        the xy element in the inertia dyad
    iyz : Sympifyable
        the yz element in the inertia dyad
    izx : Sympifyable
        the zx element in the inertia dyad

    Examples
    ========

    >>> from sympy.physics.mechanics import ReferenceFrame, inertia
    >>> N = ReferenceFrame('N')
    >>> inertia(N, 1, 2, 3)
    nx>nx> + (2)*ny>ny> + (3)*nz>nz>

    """

    if not isinstance(frame, ReferenceFrame):
        raise TypeError('Need to define the inertia in a frame')
    ol = sympify(ixx) * (frame.x | frame.x)
    ol += sympify(ixy) * (frame.x | frame.y)
    ol += sympify(izx) * (frame.x | frame.z)
    ol += sympify(ixy) * (frame.y | frame.x)
    ol += sympify(iyy) * (frame.y | frame.y)
    ol += sympify(iyz) * (frame.y | frame.z)
    ol += sympify(izx) * (frame.z | frame.x)
    ol += sympify(iyz) * (frame.z | frame.y)
    ol += sympify(izz) * (frame.z | frame.z)
    return ol
