from pyglet.gl import *
from plot_function import PlotFunction, vrange, fsubs
from math import sin, cos

class PolarFunction2d(PlotFunction):

    def __init__(self, f, intervals, options):
        self.f = f
        self.intervals = intervals
        self.options = options
        self.calculate_vertices()
        self.calculate_color_vertices()

    def calculate_vertices(self):
        if len(self.intervals) != 1:
            raise NotImplementedError("Automatic intervals not implemented.")

        t, t_min, t_max, t_steps = self.intervals[0]
        if t_steps is None: t_steps = 60

        t_set = vrange(t_min, t_max, t_steps)

        def polar_fsubs(f, t, _t):
            r = fsubs(f, t, _t)
            if r is None:
                return (None, None, None)
            return (r*cos(_t), r*sin(_t), 0.0)

        self.vertices = list( polar_fsubs(self.f, t, _t) for _t in t_set )

    def calculate_color_vertices(self):
        self.color_vertices = []
        for x, y, z in self.vertices:
            self.color_vertices.append( (0.5,0.5,1.0) )

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
