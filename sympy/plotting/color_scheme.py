from sympy import Basic, Symbol, symbols, lambdify
from util import interpolate, rinterpolate, interpolate_color, create_bounds, update_bounds

default_color_schemes = {}

class ColorScheme(object):

    def __init__(self, *args, **kwargs):
        self.f = None
        #self.gradient = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]
        self.gradient = [[0.4, 0.4, 0.4], [0.9, 0.9, 0.9]]

        e = ValueError("A custom color function must take "
                       "five arguments x,y,z,u,v (t=u for "
                       "curves) and return a tuple (r,g,b). "
                       "(Given args %s)" % str(args))

        if len(args) == 1:
            arg = args[0]
            # ColorScheme(some_existing_color_scheme)
            if isinstance(arg, ColorScheme):
                self.f = arg.f
            # ColorScheme(lambda x,y,z,u,v: (z,y,x))
            elif callable(arg):
                self.f = arg
            elif isinstance(arg, str):
                # ColorScheme('rainbow')
                if arg in default_color_schemes:
                    self.f = default_color_schemes[arg].f
                    self.gradient = default_color_schemes[arg].gradient
                # ColorScheme('z,y,x')
                else: self.f = lambdify(arg, 'x,y,z,u,v')
            else: raise e

        elif len(args) in [5, 6, 7]:
            try:
                if len(args) == 7: raise Exception() # just to skip to the next try below
                # ColorScheme( z, y, x, (0.45,0.45,0.5), (0.95,0.9,0.98) )
                fr, fg, fb, (r_min, g_min, b_min), (r_max, g_max, b_max) = args[:5]
                s = self._fill_in_vars(*args[5:])
                self.f = lambdify([fr, fg, fb], s)
                self.gradient = [[r_min, g_min, b_min], [r_max, g_max, b_max]]
            except:
                try:
                    # ColorScheme( z, (0.45,0.95), y, (0.45,0.9), x, (0.5,0.98) )
                    fr, (r_min, r_max), fg, (g_min, g_max), fb, (b_min, b_max) = args[:6]
                    s = self._fill_in_vars(*args[6:])
                    self.f = lambdify([fr, fg, fb], s)
                    self.gradient = [[r_min, g_min, b_min], [r_max, g_max, b_max]]
                except: raise e

        elif len(args) in [3, 4]:
            s = self._fill_in_vars(*args[3:])
            try:
                # ColorScheme( z, (0,0,1), (1,0,0) )
                fz, (r_min, g_min, b_min), (r_max, g_max, b_max) = args[:3]
                self.f = lambdify([fz, fz, fz], s)
                self.gradient = [[r_min, g_min, b_min], [r_max, g_max, b_max]]
            except:
                try:
                    # ColorScheme( z,y,x )
                    fr, fg, fb = args[:3]
                    self.f = lambdify([fr, fg, fb], s)
                except: raise e

        self._test_color_function(e)
        self.args = args

    def _fill_in_vars(self, *args):
        v_error = ValueError('The last optional argument '
                             'to a ColorScheme must be a '
                             'list or tuple of symbols. '
                             'To specify that one of x,y, '
                             'or z is not used, you can '
                             'specify None for that var. '
                             'Given args: %s' % str(args))
        assert len(args) in [0,1]
        defaults = symbols('xyzuv')
        if len(args) == 0: return defaults
        args = args[0]
        if not isinstance(args, (tuple, list)):
            raise v_error
        if len(args) == 0: return defaults
        for s in args:
            if s is not None and not isinstance(s, Symbol):
                raise v_error
        # when vars are given explicitly, any vars
        # not given are marked 'unbound' as to not
        # be accidentally used in an expression
        vars = [Symbol('unbound%i'%(i)) for i in xrange(1,6)]
        # interpret as t
        if len(args) == 1:
            vars[3] = args[0]
        # interpret as u,v
        elif len(args) == 2:
            vars[3] = args[0]
            vars[4] = args[1]
        # interpret as x,y,z
        elif len(args) >= 3:
            # allow some of x,y,z to be
            # left unbound if not given
            if args[0] is not None: vars[0] = args[0]
            if args[1] is not None: vars[1] = args[1]
            if args[2] is not None: vars[2] = args[2]
            # interpret the rest as t
            if len(args) >= 4:
                vars[3] = args[3]
                # ...or u,v
                if len(args) >= 5:
                    vars[4] = args[4]
        return vars

    def _test_color_function(self, e):
        try:
            assert callable(self.f)
            result = self.f(0,0,0,0,0)
            assert len(result) == 3
        except: raise e

    def __call__(self, x,y,z,u,v):
        try:    return self.f(x,y,z,u,v)
        except: return None

    def apply_to_curve(self, verts, u_set, set_len=None, inc_pos=None):
        """
        Apply this color scheme to a
        set of vertices over a single
        independent variable u.
        """
        bounds = create_bounds()
        cverts = list()
        if callable(set_len): set_len(len(u_set)*2)
        # calculate f() = r,g,b for each vert
        # and find the min and max for r,g,b
        for _u in xrange(len(u_set)):
            if verts[_u] is None:
                cverts.append(None)
            else:
                x,y,z = verts[_u]
                u,v = u_set[_u], None
                c = list(self(x,y,z,u,v))
                update_bounds(bounds, c)
                cverts.append(c)
            if callable(inc_pos): inc_pos()
        # scale the calculated values to the
        # color gradient
        for _u in xrange(len(u_set)):
            if cverts[_u] is not None:
                c = cverts[_u] # modify it in place using c as an alias
                for _c in xrange(3): # _c points to one of r,g,b
                    # scale from [f_min, f_max] to [0,1]
                    i = rinterpolate(bounds[_c][0], bounds[_c][1], c[_c])
                    # scale from [0,1] to [gradient_min, gradient_max]
                    c[_c] = interpolate(self.gradient[0][_c], self.gradient[1][_c], i)
            if callable(inc_pos): inc_pos()
        return cverts

    def apply_to_surface(self, verts, u_set, v_set, set_len=None, inc_pos=None):
        """
        Apply this color scheme to a
        set of vertices over two
        independent variables u and v.
        """
        bounds = create_bounds()
        cverts = list()
        if callable(set_len): set_len(len(u_set)*len(v_set)*2)
        # calculate f() = r,g,b for each vert
        # and find the min and max for r,g,b
        for _u in xrange(len(u_set)):
            column = list()
            for _v in xrange(len(v_set)):
                if verts[_u][_v] is None:
                    column.append(None)
                else:
                    x,y,z = verts[_u][_v]
                    u,v = u_set[_u], v_set[_v]
                    c = list(self(x,y,z,u,v))
                    update_bounds(bounds, c)
                    column.append(c)
                if callable(inc_pos): inc_pos()
            cverts.append(column)
        # scale the calculated values to the
        # color gradient
        for _u in xrange(len(u_set)):
            for _v in xrange(len(v_set)):
                if cverts[_u][_v] is not None:
                    c = cverts[_u][_v] # modify it in place using c as an alias
                    for _c in xrange(3): # _c points to one of r,g,b
                        # scale from [f_min, f_max] to [0,1]
                        i = rinterpolate(bounds[_c][0], bounds[_c][1], c[_c])
                        # scale from [0,1] to [gradient_min, gradient_max]
                        c[_c] = interpolate(self.gradient[0][_c], self.gradient[1][_c], i)
                if callable(inc_pos): inc_pos()
        return cverts

    def str_base(self): return ", ".join(str(a) for a in self.args)
    def __repr__(self): return "%s" % (self.str_base())

x,y,z,t,u,v = symbols('xyztuv')
default_color_schemes['rainbow'] = ColorScheme( z, y, x, (0.45,0.45,0.5), (0.95,0.9,0.98) )
default_color_schemes['zfade'] = ColorScheme( z, (0.5,0.5,0.9), (0.9,0.5,0.5) )

