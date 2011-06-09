__all__ = ['cross',
           'dot',
           'dynamicsymbols',
           'express']

from sympy.physics.classical.essential import Vector
from sympy.physics.classical.dynamicsymbol import DynamicSymbol
from sympy import sympify

def cross(vec1, vec2):
    """Cross product convenience wrapper for Vector.cross(): \n"""
    if not isinstance(vec1, Vector):
        raise TypeError('Cross product is between two vectors')
    return vec1.cross(vec2)
cross.__doc__ += Vector.cross.__doc__

def dot(vec1, vec2):
    """Dot product convenience wrapper for Vector.dot(): \n"""
    if not isinstance(vec1, Vector):
        raise TypeError('Doc product is between two vectors')
    return vec1.dot(vec2)
dot.__doc__ += Vector.dot.__doc__

def dynamicsymbols(basename, count, diffno=0):
    """Returns a list of DynamicSymbols, and a number of their time derivatives.
    Needs a base name supplied, number of DynamicSymbols, and number of
    time derivatives of the Dynamic Symbols.

    """

    sympify(basename)
    outlist = []
    for i in range(diffno + 1):
        innerlist = []
        tag = ''
        for j in range(i):
            tag += 'd'
        for j in range(count):
            innerlist.append(DynamicSymbol(basename + str(j) + tag))
        outlist.append(innerlist)
    try:
        outlist[1]
        return outlist
    except:
        return outlist[0]

def express(vec, frame):
    """Express convenience wrapper for Vector.express(): \n"""
    if not isinstance(vec, Vector):
        raise TypeError('Can only express Vectors')
    return vec.express(frame)
express.__doc__ += Vector.express.__doc__
