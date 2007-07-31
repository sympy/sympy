from pyglet.gl import *

def get_model_matrix(array_type=c_float, glGetMethod=glGetFloatv):
    """
    Returns the current modelview matrix.
    """
    m = (array_type*16)()
    glGetMethod(GL_MODELVIEW_MATRIX, m)
    return m

def get_projection_matrix(array_type=c_float, glGetMethod=glGetFloatv):
    """
    Returns the current modelview matrix.
    """
    m = (array_type*16)()
    glGetMethod(GL_PROJECTION_MATRIX, m)
    return m

def get_viewport():
    """
    Returns the current viewport.
    """
    m = (c_int*4)()
    glGetIntegerv(GL_VIEWPORT, m)
    return m

def get_direction_vectors():
    m = get_model_matrix()
    return ((m[0], m[4], m[8]),
            (m[1], m[5], m[9]),
            (m[2], m[6], m[10]))

def get_view_direction_vectors():
    m = get_model_matrix()
    return ((m[0], m[1], m[2]),
            (m[4], m[5], m[6]),
            (m[8], m[9], m[10]))

def get_basis_vectors():
    return ((1,0,0), (0,1,0), (0,0,1))

def screen_to_model(x,y,z):
    m = get_model_matrix(c_double, glGetDoublev)
    p = get_projection_matrix(c_double, glGetDoublev)
    w = get_viewport()
    mx,my,mz = c_double(),c_double(),c_double()
    gluUnProject(x,y,z,m,p,w,mx,my,mz)
    return mx,my,mz

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
    m = get_model_matrix()
    m[0] =1;m[1] =0;m[2] =0
    m[4] =0;m[5] =1;m[6] =0
    m[8] =0;m[9] =0;m[10]=1
    glLoadMatrixf(m)

def find_bounds_2d(vertices):
    bounds = [[None,None], [None,None], [None,None]]
    for v in vertices:
        v = v[0]
        if v is None: continue
        for axis in range(len(bounds)):
            if v[axis] is not None:
                if bounds[axis][0] is None: bounds[axis][0] = v[axis]
                else: bounds[axis][0] = min( [v[axis], bounds[axis][0]] )
                if bounds[axis][1] is None: bounds[axis][1] = v[axis]
                else: bounds[axis][1] = max( [v[axis], bounds[axis][1]] )
    return bounds

def find_bounds_3d(vertices):
    bounds = [[None,None], [None,None], [None,None]]
    for w in vertices:
        for v in w:
            v = v[0]
            if v is None: continue
            for axis in range(len(bounds)):
                if v[axis] is not None:
                    if bounds[axis][0] is None: bounds[axis][0] = v[axis]
                    else: bounds[axis][0] = min( [v[axis], bounds[axis][0]] )
                    if bounds[axis][1] is None: bounds[axis][1] = v[axis]
                    else: bounds[axis][1] = max( [v[axis], bounds[axis][1]] )
    return bounds

def interpolate(a_min, a_max, a_ratio):
    return a_min + a_ratio * (a_max - a_min)

def rinterpolate(a_min, a_max, a_value):
    a_range = a_max-a_min
    if a_range == 0:
        a_range = 1.0
    return (a_value - a_min) / float(a_range)

def interpolate_color(color1, color2, ratio):
    return [interpolate(color1[i], color2[i], ratio) for i in range(3)]

def calc_color_1(color_function, p, pbounds, i, ibounds):
    if not p:
        return (0, 0, 0)
    k = ( rinterpolate(pbounds[0][0], pbounds[0][1], p[0]),
          rinterpolate(pbounds[1][0], pbounds[1][1], p[1]),
          rinterpolate(pbounds[2][0], pbounds[2][1], p[2]),
          rinterpolate(ibounds[0][0], ibounds[0][1], i[0]) )
    return color_function(*k)

def calc_color_2(color_function, p, pbounds, i, ibounds):
    if not p:
        return (0, 0, 0)
    k = ( rinterpolate(pbounds[0][0], pbounds[0][1], p[0]),
          rinterpolate(pbounds[1][0], pbounds[1][1], p[1]),
          rinterpolate(pbounds[2][0], pbounds[2][1], p[2]),
          rinterpolate(ibounds[0][0], ibounds[0][1], i[0]),
          rinterpolate(ibounds[1][0], ibounds[1][1], i[1]) )
    return color_function(*k)
