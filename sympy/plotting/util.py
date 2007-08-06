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
    return float(mx.value),float(my.value),float(mz.value)

def model_to_screen(x,y,z):
    m = get_model_matrix(c_double, glGetDoublev)
    p = get_projection_matrix(c_double, glGetDoublev)
    w = get_viewport()
    mx,my,mz = c_double(),c_double(),c_double()
    gluProject(x,y,z,m,p,w,mx,my,mz)
    return float(mx.value),float(my.value),float(mz.value)

def vec_subs(a,b):
    return tuple(a[i]-b[i] for i in xrange(len(a)))

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

def scale_value(v, v_min, v_len):
    return (v-v_min)/v_len

def scale_value_list(flist):
        v_min, v_max = min(flist), max(flist)
        v_len = v_max-v_min
        return list(scale_value(f,v_min,v_len) for f in flist)

def strided_range(r_min, r_max, stride):
    if abs(r_min-r_max) < 0.001: return []
    try: xrange(int(r_min-r_max))
    except: return []
    assert r_min < r_max
    r_min_s = r_min % stride
    r_max_s = r_max % stride
    if r_min_s and r_max_s == 0.0:
        r_max_s += stride
    elif r_max_s and r_min_s == 0.0:
        r_min_s += stride
    r_min -= r_min_s
    r_max += r_max_s
    r = list()
    r_steps = int( (r_max-r_min) / stride )
    r = list(r_min+e*stride for e in xrange(r_steps+1))
    #print "%s-%s: %s" % (r_min, r_max, r)
    return r

def parse_option_string(s):
    if not isinstance(s, str):
        return None
    options = {}
    for token in s.split(';'):
        pieces = token.split('=')
        if len(pieces) == 1:
            option, value = pieces[0], ""
        elif len(pieces) == 2:
            option, value = pieces
        else:
            raise ValueError("Plot option string '%s' is malformed." % (s))
        options[option.strip()] = value.strip()
    return options
