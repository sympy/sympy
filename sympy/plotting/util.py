from pyglet.gl import *

def get_matrix():
    """
    Returns the current modelview matrix.
    """
    m = (c_float*16)()
    glGetFloatv(GL_MODELVIEW_MATRIX, m)
    return m

def get_projection_matrix():
    """
    Returns the current modelview matrix.
    """
    m = (c_float*16)()
    glGetFloatv(GL_PROJECTION_MATRIX, m)
    return m

def billboard_matrix():
    """
    Removes rotational components of
    current matrix so that primitives
    are always drawn facing the viewer.

    |1|0|0|x|
    |0|1|0|x|
    |0|0|1|x| (x means left unchanged)
    |x|x|x|x|
    """
    m = get_matrix()
    m[0] =1;m[1] =0;m[2] =0
    m[4] =0;m[5] =1;m[6] =0
    m[8] =0;m[9] =0;m[10]=1
    glLoadMatrixf(m)

def get_direction_vectors():
    m = get_matrix()
    return ((m[0], m[4], m[8]),
            (m[1], m[5], m[9]),
            (m[2], m[6], m[10]))

def get_view_direction_vectors():
    m = get_matrix()
    return ((m[0], m[1], m[2]),
            (m[4], m[5], m[6]),
            (m[8], m[9], m[10]))

def get_basis_vectors():
    return ((1,0,0), (0,1,0), (0,0,1))

def invert_vector(v):
    return (-v[0], -v[1], -v[2])

def get_matrix_d():
    """
    Returns the current modelview matrix (double).
    """
    m = (c_double*16)()
    glGetDoublev(GL_MODELVIEW_MATRIX, m)
    return m

def get_projection_matrix_d():
    """
    Returns the current modelview matrix (double).
    """
    m = (c_double*16)()
    glGetDoublev(GL_PROJECTION_MATRIX, m)
    return m

def get_viewport():
    v = (c_int*4)()
    glGetIntegerv(GL_VIEWPORT, v)
    return v

def model_to_screen(c):
    m = get_matrix_d()
    p = get_projection_matrix_d()
    v = get_viewport()
    x,y,z = c_double(),c_double(),c_double()
    if GL_FALSE == gluProject(c[0], c[1], c[2],
                              m, p, v, x, y, z):
        raise Exception("gluProject failed.")
    return (x,y,z)

def model_to_screen_ratio(c):
    m = get_matrix_d()
    p = get_projection_matrix_d()
    v = get_viewport()
    x,y,z = c_double(),c_double(),c_double()
    if GL_FALSE == gluProject(c[0], c[1], c[2],
                              m, p, v, x, y, z):
        raise Exception("gluProject failed.")
    return (x.value/float(v[2])-0.5,y.value/float(v[3])-0.5)

def inner_product(v1, v2):
    s = 0
    for i in xrange(len(v1)):
        s += v1[i] * v2[i]
    return s

def mat_mult(m, v):
    if len(v) != 4:
        v = (v[0],v[1],v[2],1)
    return [inner_product((m[i], m[i+1], m[i+2], m[i+3]), v)
            for i in range(0,len(m),4)]

def frange(start, end=None, inc=None):
    """
    A range function, that does accept float increments...
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/66472
    """

    if end is None:
        end = start + 0.0
        start = 0.0

    if inc is None:
        inc = 1.0

    L = []
    while 1:
        next = start + len(L) * inc
        if inc > 0 and next >= end:
            break
        elif inc < 0 and next <= end:
            break
        L.append(next)
        
    return L
