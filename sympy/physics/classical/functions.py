__all__ = ['cross',
           'dot',
           'dynamicsymbols',
           'express',
           'outer']

from sympy.physics.classical.essential import Vector, Dyad
from sympy.physics.classical.dynamicsymbol import DynamicSymbol
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
    if not isinstance(vec1, (Vector, Dyad)):
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
