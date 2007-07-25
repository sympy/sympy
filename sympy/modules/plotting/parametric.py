from pyglet.gl import *
from plot_function import PlotFunction, PlotFunctionRegistry, fsubs, get_vars, count_vars, vrange

class ParametricFunction(PlotFunction):

    def __new__(cls, functions, intervals, options):
        if len(intervals) == 2:
            raise NotImplementedError("Parametric surfaces not supported at this time.")
        if len(intervals) > 2:
            raise Exception("A parametric function cannot use more than 2 parameters.")
        d = len(functions)
        v = None
        # Check that we only have one parameter
        for f in functions:
            for iv in get_vars(f):
                if v is None:
                    v = iv
                elif iv != v:
                    raise Exception("A parametric function cannot use more than 1 parameter (multiple Symbols found in functions).")
        if d in [2,3] and len(intervals) == 1: return object.__new__(ParametricCurve, functions, intervals, options)
        elif d == 1:
            raise ValueError("Cannot plot a parametric function with only 1 dimension.")
        else:
            raise ValueError("Cannot plot a parametric function with %i dimensions." % d)

class ParametricCurve(PlotFunction):

    def __init__(self, functions, intervals, options):
        self.f = functions
        self.intervals = intervals
        self.options = options
        self.calculate_vertices()
        self.calculate_color_vertices()

    def calculate_vertices(self):
        if len(self.intervals) != 1:
            raise NotImplementedError("Automatic interval not implemented.")

        t, t_min, t_max, t_steps = self.intervals[0]
        if t_steps is None: t_steps = 60

        t_set = vrange(t_min, t_max, t_steps)

        if len(self.f) == 2:
            self.vertices = list( (fsubs(self.f[0], t, t_e), fsubs(self.f[1], t, t_e), 0.0) for t_e in t_set )
        elif len(self.f) == 3:
            self.vertices = list( (fsubs(self.f[0], t, t_e), fsubs(self.f[1], t, t_e), fsubs(self.f[2], t, t_e)) for t_e in t_set )

    def calculate_color_vertices(self):
        self.color_vertices = []
        for x, y, z in self.vertices:
            self.color_vertices.append( (0.7, 0.75, 0.8) )

    def render(self):
        glBegin(GL_LINE_STRIP)
        for x in range(0, len(self.vertices)):
            if self.vertices[x][1] is None:
                glEnd()
                glBegin(GL_LINE_STRIP)
                continue
            glColor3f(*self.color_vertices[x])
            glVertex3f(*self.vertices[x])
        glEnd()

PlotFunctionRegistry.register('parametric', ParametricFunction)
