from pyglet.gl import *

from plot_object import PlotObject
from plot_interval import PlotInterval

from util import find_bounds_2d, calc_color
from util import interpolate as intrp

class PlotCurve(PlotObject):
    """
    """
    def __init__(self, function, t_interval, **kwargs):
        self.function        =  function
        self.line_color      =  kwargs.pop('line_color', self.default_line_color)
        self.t_interval      =  PlotInterval.try_parse(t_interval)
        self.display_list    =  -1

        self.vertices   = list([self.function(t), t] for t in self.t_interval.vrange())

        self.pbounds    = find_bounds_2d(self.vertices)
        self.ibounds    = [[self.t_interval.v_min, self.t_interval.v_max]]

        self.c_vertices = list(calc_color(self.line_color,
                                          self.vertices[t][0], self.pbounds,
                                          self.vertices[t][1:], self.ibounds)
                                          for t in xrange(self.t_interval.v_len))

    def default_line_color(self, x, y, z, t):
        r = intrp(0.6, 0.97, x)
        g = intrp(0.6, 0.95, y)
        b = intrp(0.6, 0.97, z)
        return r, g, b

    def compile(self):
        self.display_list = glGenLists(1)
        glNewList(self.display_list, GL_COMPILE)
        glBegin(GL_LINE_STRIP)
        for t in xrange(self.t_interval.v_len):
            p = self.vertices[t][0]
            c = self.c_vertices[t]
            if p is None:
                glEnd()
                glBegin(GL_LINE_STRIP)
                continue
            glColor3f(*c)
            glVertex3f(*p)
        glEnd()
        glEndList()

    def draw(self):
        if GL_FALSE == glIsList(self.display_list):
            self.compile()

        glCallList(self.display_list)
