__all__ = ['cross',
           'dot',
           'dynamicsymbols',
           'express']

from sympy.physics.classical.essential import Vector

def cross(vec1, vec2):
    """
    Returns the dot product of the two vectors.
    """
    assert isinstance(vec1, Vector), 'Cross product is between two vectors'
    return vec1.cross(vec2)

def dot(vec1, vec2):
    """
    Returns the dot product of the two vectors
    """
    assert isinstance(vec1, Vector), 'Dot product is between two vectors'
    return vec1.dot(vec2)

def dynamicsymbols(basename, count, diffno=0):
    """
    Returns a list of DynamicSymbols, and a number of their time derivatives.
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
    """
    Returns the input Vector in the input ReferenceFrame.
    """
    assert isinstance(vec, Vector), 'Can only express Vectors in a frame'
    return vec.express(frame)

