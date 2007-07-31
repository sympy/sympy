from sympy import Basic, symbols, sin, cos, pi
from plot_interval import PlotInterval

from math import sin as psin
from math import cos as pcos

plot_modes = {}

plot_mode_aliases = {
    0: [], # position 0 holds every alias string
    1: {1: {},
        2: {}},
    2: {1: {},
        2: {}},
    3: {1: {},
        2: {}},
}

plot_mode_defaults = {
    1: {1: '',
        2: ''},
    2: {1: '',
        2: ''},
    3: {1: '',
        2: ''},
}

plot_mode_default_i_vars = {}
plot_mode_default_intervals = {}

default_mode_error = ("Don't know how to plot %i independent " +
                      "and %i dependent variables.")
unknown_mode = ("Didn't recognize plot mode string '%s'.")
alias_error = ("No appropriate %s mode found for %i independent " +
               "and %i dependent variables.")

def get_plot_mode(mode_str, d_var_c, i_var_c):
    if not mode_str:
        def find_mode_from_defaults(i):            
            try: return plot_mode_defaults[d_var_c][i]
            except:
                if i < 2: return find_mode_from_defaults(i+1)
                else: raise ValueError(default_mode_error % (i_var_c, d_var_c))
        mode_str = find_mode_from_defaults(i_var_c)

    if mode_str in plot_mode_aliases[0]:
        def find_mode_from_alias(mode_str, i, i_var_c):
            try: return plot_mode_aliases[d_var_c][i][mode_str]
            except:
                if i < 2: return find_mode_from_alias(mode_str, i+1, i_var_c)
                else: raise ValueError(alias_error % (mode_str, i_var_c, d_var_c))
        mode_str = find_mode_from_alias(mode_str, i_var_c, i_var_c)

    try: mode = plot_modes[mode_str]
    except: raise ValueError(unknown_mode % (mode_str))

    return mode, mode_str

def fill_i_vars(mode_str, i_vars):
    default_i_vars = plot_mode_default_i_vars[mode_str]
    for i in xrange(len(i_vars), len(default_i_vars)):
        i_vars.append(default_i_vars[i])
    return i_vars

def fill_intervals(mode_str, i_vars, intervals):
    default_intervals = plot_mode_default_intervals[mode_str]
    for i in xrange(len(intervals), len(default_intervals)):
        intervals.append(PlotInterval(i_vars[i], *default_intervals[i]))
    return intervals

def plot_mode(i_var_str, d_var_str, default_intervals, aliases=None, default=False):
    default_i_vars = symbols(i_var_str)
    default_d_vars = symbols(d_var_str)

    d = len(default_d_vars)
    i = len(default_i_vars)

    if aliases is None: aliases = []

    def plot_mode_decorator(mode_function):
        name = mode_function.__name__
        plot_modes[name] = mode_function

        for a in aliases:
            plot_mode_aliases[d][i][a] = name
            if a not in plot_mode_aliases[0]:
                plot_mode_aliases[0].append(a)
        if default:
            plot_mode_defaults[d][i] = name

        plot_mode_default_i_vars[name] = default_i_vars
        plot_mode_default_intervals[name] = default_intervals

        def base_plot_mode(*args):
            d_vars = list(args[:d])
            i_vars = list(args[d:])
            if len(d_vars) != d:
                raise ValueError(("%s takes %i dependent function(s) " +
                                  " (%i given).") % (name, d, len(d_vars)))
            if len(i_vars) > i:
                raise ValueError(("%s takes at most %i independent" +
                                  " variables (%i given).") % (name, i, len(i_vars)) )
            for j in range(len(i_vars), i):
                i_vars.append(default_i_vars[j])
            return mode_function(*(d_vars+i_vars))

        return base_plot_mode

    return plot_mode_decorator

def float_vector(x, y, z):
    try:
        #print (x,y,z)
        float(x.evalf())
        float(y.evalf())
        float(z.evalf())
    except:
        return None
    return x, y, z

@plot_mode('t', 'xyz', [[0,2*pi,60]], aliases=['parametric'], default=True)
def parametric_curve3d(fx, fy, fz, t):
    fx = Basic.sympify(fx)
    fy = Basic.sympify(fy)
    fz = Basic.sympify(fz)
    def _f(_t):
        try:
            x = fx.subs(t,_t)
            y = fy.subs(t,_t)
            z = fz.subs(t,_t)
        except:
            return None
        return float_vector(x,y,z)
    return _f

@plot_mode('t', 'xy', [[0,2*pi,60]], aliases=['parametric'], default=True)
def parametric_curve2d(fx, fy, t):
    return parametric_curve3d(fx, fy, 0.0, t)

@plot_mode('uv', 'xyz', [[-1,1,10], [-1,1,10]], aliases=['parametric'], default=True)
def parametric_surface(fx, fy, fz, u, v):
    fx = Basic.sympify(fx)
    fy = Basic.sympify(fy)
    fz = Basic.sympify(fz)
    def _f(_u, _v):
        try:
            x = fx.subs(u,_u).subs(v,_v)
            y = fy.subs(u,_u).subs(v,_v)
            z = fz.subs(u,_u).subs(v,_v)
        except:
            return None
        return float_vector(x,y,z)
    return _f

@plot_mode('x', 'y', [[-5,5,60]], aliases=['cartesian'], default=True)
def cartesian_curve(fy, x):
    return parametric_curve2d(x, fy, x)

@plot_mode('xy', 'z', [[-1,1,10], [-1,1,10]], aliases=['cartesian'], default=True)
def cartesian_surface(fz, x, y):
    return parametric_surface(x, y, fz, x, y)

@plot_mode('t', 'r', [[0,2*pi,60]], aliases=['polar'])
def polar_curve(fr, t):
    #return parametric_curve2d(fr * cos(t), fr * sin(t), t)
    fr = Basic.sympify(fr)
    def _f(_t):
        try:
            _t = float(_t)
            _r = float(fr.subs(t,_t).evalf())
        except:
            return None
        return (_r*pcos(_t), _r*psin(_t), 0.0)
    return _f

@plot_mode('th', 'r', [[0,2*pi,20], [-0.5,0.5,6]], aliases=['cylindrical','polar'])
def cylindrical_surface(fr, t, h):
    #return parametric_surface(fr * cos(t), fr * sin(t), h, t, h)
    fr = Basic.sympify(fr)
    def _f(_t, _h):
        try:
            _r = float(fr.subs(t,_t).subs(h,_h).evalf())
        except:
            return None
        _t = float(_t)
        _h = float(_h)
        return (_r*pcos(_t), _r*psin(_t), _h)
    return _f

@plot_mode('pt', 'r', [[0,pi,12], [0,2*pi,16]], aliases=['spherical'])
def spherical_surface(fr, p, t):
    #return parametric_surface(fr * sin(p) * cos(t), fr * sin(p) * sin(t), fr * cos(p), t, h)
    fr = Basic.sympify(fr)
    def _f(_p, _t):
        try:
            _r = float(fr.subs(p,_p).subs(t,_t).evalf())
        except:
            return None
        _p = float(_p)
        _t = float(_t)
        return (_r*psin(_p)*pcos(_t), _r*psin(_p)*psin(_t), _r*pcos(_p))
    return _f

def print_mode_aliases():
    for d in plot_mode_aliases:
        for i in plot_mode_aliases[d]:
            print "[%i][%i]: %s" % (d,i,plot_mode_aliases[d][i])

def print_mode_defaults():
    for d in plot_mode_defaults:
        for i in plot_mode_defaults[d]:
            print "[%i][%i]: %s" % (d,i,plot_mode_defaults[d][i])

def torus(radius=2.5, thickness=1.0):
    a=radius
    b=thickness
    return (a+b*cos(x))*cos(y), (a+b*cos(x))*sin(y), b*sin(x)
