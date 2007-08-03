from plot_curve import PlotCurve
from plot_surface import PlotSurface
from util import scale_value

from sympy import Pi
from math import sin as p_sin
from math import cos as p_cos

def float_vec3(f):
    def inner(*args):
        try:
            v = f(*args)
            return float(v[0]), float(v[1]), float(v[2])
        except:
            return None
    return inner

class Cartesian2D(PlotCurve):
    i_vars, d_vars = 'x', 'y'
    intervals = [[-5,5,60]]
    aliases = ['cartesian']
    is_default = True
    
    def _get_evaluator(self):
        fy = self.d_vars[0]
        x  = self.t_interval.v
        @float_vec3
        def e(_x):
            return ( _x, fy.subs(x, _x), 0.0 )
        return e

    def calculate_one_cvert(self, t):
        return self.color(self.st_set[t],
                          scale_value(self.verts[t][1], self.bounds[1][0], self.bounds[1][2]),
                          scale_value(self.verts[t][2], self.bounds[2][0], self.bounds[2][2]),
                          self.st_set[t], None)

class Cartesian3D(PlotSurface):
    i_vars, d_vars = 'xy', 'z'
    intervals = [[-1,1,20], [-1,1,20]]
    aliases = ['cartesian', 'monge']
    is_default = True
    
    def _get_evaluator(self):
        fz = self.d_vars[0]
        x  = self.u_interval.v
        y  = self.v_interval.v
        @float_vec3
        def e(_x, _y):
            return ( _x, _y, fz.subs(x, _x).subs(y, _y) )
        return e

    def calculate_one_cvert(self, u, v):
        vert = self.verts[u][v]
        return self.color(self.su_set[u], self.sv_set[v],
                          scale_value(vert[2], self.bounds[2][0], self.bounds[2][2]),
                          self.su_set[u], self.sv_set[v])


class ParametricCurve2D(PlotCurve):
    i_vars, d_vars = 't', 'xy'
    intervals = [[0,1,60]]
    aliases = ['parametric']
    is_default = True

    def _get_evaluator(self):
        fx, fy = self.d_vars
        t  = self.t_interval.v
        @float_vec3
        def e(_t):
            return ( fx.subs(t, _t),
                     fy.subs(t, _t),
                     0.0 )
        return e

class ParametricCurve3D(PlotCurve):
    i_vars, d_vars = 't', 'xyz'
    intervals = [[0,1,60]]
    aliases = ['parametric']
    is_default = True

    def _get_evaluator(self):
        fx, fy, fz = self.d_vars
        t  = self.t_interval.v
        @float_vec3
        def e(_t):
            return ( fx.subs(t, _t),
                     fy.subs(t, _t),
                     fz.subs(t, _t) )
        return e

class ParametricSurface(PlotSurface):
    i_vars, d_vars = 'uv', 'xyz'
    intervals = [[-1,1,20], [-1,1,20]]
    aliases = ['parametric']
    is_default = True

    def _get_evaluator(self):
        fx, fy, fz = self.d_vars
        u  = self.u_interval.v
        v  = self.v_interval.v
        @float_vec3
        def e(_u, _v):
            return ( fx.subs(u, _u).subs(v, _v),
                     fy.subs(u, _u).subs(v, _v),
                     fz.subs(u, _u).subs(v, _v) )
        return e

class Polar(PlotCurve):
    i_vars, d_vars = 't', 'r'
    intervals = [[0,2*Pi,64]]
    aliases = ['polar']
    is_default = False

    def _get_evaluator(self):
        fr = self.d_vars[0]
        t  = self.t_interval.v
        def e(_t):
            try:
                _r = float( fr.subs(t, _t) )
                return ( _r*p_cos(_t), _r*p_sin(_t), 0.0 )
            except: return None
        return e

class Cylindrical(PlotSurface):
    i_vars, d_vars = 'th', 'r'
    intervals = [[0,2*Pi,24], [-1,1,10]]
    aliases = ['cylindrical', 'polar']
    is_default = False

    def _get_evaluator(self):
        fr = self.d_vars[0]
        t  = self.u_interval.v
        h  = self.v_interval.v
        def e(_t, _h):
            try:
                _r = float( fr.subs(t, _t).subs(h, _h) )
                return ( _r*p_cos(_t), _r*p_sin(_t), _h )
            except: return None
        return e

class Spherical(PlotSurface):
    i_vars, d_vars = 'tp', 'r'
    intervals = [[0,2*Pi,24], [0,Pi,12]]
    aliases = ['spherical']
    is_default = False

    def _get_evaluator(self):
        fr = self.d_vars[0]
        t  = self.u_interval.v
        p  = self.v_interval.v
        def e(_t, _p):
            try:
                _r = float( fr.subs(t, _t).subs(p, _p) )
                return ( _r*p_cos(_t)*p_sin(_p),
                         _r*p_sin(_t)*p_sin(_p),
                         _r*p_cos(_p) )
            except: return None
        return e

Cartesian2D._register()
Cartesian3D._register()
ParametricCurve2D._register()
ParametricCurve3D._register()
ParametricSurface._register()
Polar._register()
Cylindrical._register()
Spherical._register()
