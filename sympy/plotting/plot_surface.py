from pyglet.gl import *

from plot_object import PlotObject
from plot_interval import PlotInterval

from util import find_bounds_3d, calc_color
from util import interpolate as intrp

WIREFRAME = 1
SOLID = 2
SOLID_AND_WIREFRAME = WIREFRAME | SOLID

class PlotSurface(PlotObject):
    """
    """
    def __init__(self, function, u_interval, v_interval, **kwargs):
        self.function           =  function
        self.u_interval         =  PlotInterval.try_parse(u_interval)
        self.v_interval         =  PlotInterval.try_parse(v_interval)
        self.surface_color      =  kwargs.pop('surface_color', self.default_surface_color)
        self.line_color         =  kwargs.pop('line_color', self.default_line_color)
        self.render_style       =  kwargs.pop('render_style', SOLID_AND_WIREFRAME)
        self.display_list       =  -1

        if isinstance(self.render_style, str):
            try:
                self.render_style = int(eval(self.render_style))
            except:
                pass

        v_set = list(float(v) for v in self.v_interval.vrange())
        u_set = list(float(u) for u in self.u_interval.vrange())
        self.vertices   = list(list([self.function(u,v), u, v]
                                    for v in v_set)
                                    for u in u_set)

        self.pbounds    = find_bounds_3d(self.vertices)
        self.ibounds    = [[self.u_interval.v_min, self.u_interval.v_max],
                           [self.v_interval.v_min, self.v_interval.v_max]]

        self.c_vertices = list(list((calc_color(self.line_color,
                                                self.vertices[u][v][0], self.pbounds,
                                                self.vertices[u][v][1:3], self.ibounds),

                                     calc_color(self.surface_color,
                                                self.vertices[u][v][0], self.pbounds,
                                                self.vertices[u][v][1:3], self.ibounds))

                                          for v in xrange(self.v_interval.v_len))
                                          for u in xrange(self.u_interval.v_len))

    def default_surface_color(self, x, y, z, u, v):
        r = intrp(0.6, 0.97, z)
        g = intrp(0.6, 0.95, y)
        b = intrp(0.6, 0.97, x)
        return r, g, b

    def default_line_color(self, x, y, z, u, v):
        return (0.5, 0.5, 0.5)

    def gl_surface(self, fill):
        for u in xrange(1, self.u_interval.v_len):
            glBegin(GL_QUAD_STRIP)
            for v in xrange(self.v_interval.v_len):
                pa = self.vertices[u-1][v][0]
                pb = self.vertices[u][v][0]
                if fill:
                    ca = self.c_vertices[u-1][v][1]
                    cb = self.c_vertices[u][v][1]
                else:
                    ca = self.c_vertices[u-1][v][0]
                    cb = self.c_vertices[u][v][0]
                if pa is None or pb is None:
                    glEnd()
                    glBegin(GL_QUAD_STRIP)
                    continue
                glColor3f(*ca)
                glVertex3f(*pa)
                glColor3f(*cb)
                glVertex3f(*pb)
            glEnd()

    def compile(self):
        self.display_list = glGenLists(2)

        glNewList(self.display_list, GL_COMPILE)
        self.gl_surface(True)
        glEndList()

        glNewList(self.display_list+1, GL_COMPILE)
        self.gl_surface(False)
        glEndList()

    def draw(self):
        if glIsList(self.display_list) == GL_FALSE:
            self.compile()

        glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT)
        
        if self.render_style & SOLID:
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
            glCallList(self.display_list)

        if self.render_style & WIREFRAME:
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
            glEnable(GL_POLYGON_OFFSET_LINE)
            glPolygonOffset(-0.005, -50.0)
            glCallList(self.display_list+1)
            glDisable(GL_POLYGON_OFFSET_LINE)

        glPopAttrib()
