from pyglet.gl import *
from plot_function import PlotFunction, PlotFunctionRegistry, vrange, fsubs, rinterpolate, interpolate
from math import sin, cos

class SphericalFunction(PlotFunction):

    def __new__(cls, f, intervals, options):
        return object.__new__(cls, f, intervals, options)

    def __init__(self, f, intervals, options):
        self.f = f
        self.intervals = intervals
        self.options = options
        self.calculate_vertices()
        self.calculate_bounding_box()
        self.calculate_color_vertices()
        
    def calculate_vertices(self):
        if len(self.intervals) != 2:
            raise NotImplementedError("Automatic intervals not implemented.")

        t, t_min, t_max, t_steps = self.intervals[0]
        if t_steps == None: t_steps = 40

        p, p_min, p_max, p_steps = self.intervals[1]
        if p_steps == None: p_steps = 20

        t_set = vrange(t_min, t_max, t_steps)
        p_set = vrange(p_min, p_max, p_steps)

        def cyl_fsubs(f, t, _t, p, _p):
            r = fsubs(f, t, _t, p, _p)
            if r == None:
                return (None, None, None)
            return r*sin(_p)*cos(_t), r*cos(_p), r*sin(_p)*sin(_t)

        self.vertices = list( list( cyl_fsubs(self.f, t, _t, p, _p)
                                    for _p in p_set ) for _t in t_set )

    def calculate_bounding_box(self):
        x_len = len(self.vertices)
        y_len = len(self.vertices[0])

        for x in range(x_len):
            for y in range(y_len):
                if self.vertices[x][y][0] != None:
                    self.x_min = min([self.x_min, self.vertices[x][y][0]])
                    self.x_max = max([self.x_max, self.vertices[x][y][0]])

                if self.vertices[x][y][1] != None:
                    self.y_min = min([self.y_min, self.vertices[x][y][1]])
                    self.y_max = max([self.y_max, self.vertices[x][y][1]])

                if self.vertices[x][y][2] != None:
                    self.z_min = min([self.z_min, self.vertices[x][y][2]])
                    self.z_max = max([self.z_max, self.vertices[x][y][2]])

    def calculate_color_vertices(self, min_brightness=0.55, max_brightness=0.95):
        x_len = len(self.vertices)
        y_len = len(self.vertices[0])

        self.color_vertices = list(list(None for y in range(y_len)) for x in range(x_len))

        for x in range(x_len):
            for y in range(y_len):
                if self.vertices[x][y][2] == None:
                    self.color_vertices[x][y] = (0.0, 0.0, 0.0)
                else:
                    cx = interpolate( min_brightness, max_brightness,
                                      rinterpolate(self.x_min,
                                                   self.x_max,
                                                   self.vertices[x][y][0]) )
                    cy = interpolate( min_brightness, max_brightness,
                                      rinterpolate(self.y_min,
                                                   self.y_max,
                                                   self.vertices[x][y][1]) )
                    cz = interpolate( min_brightness, max_brightness,
                                      rinterpolate(self.z_min,
                                                   self.z_max,
                                                   self.vertices[x][y][2]) )
                    self.color_vertices[x][y] = (cz, cx, cy)

    def render(self):
        t_len = len(self.vertices)
        h_len = len(self.vertices[0])

        for t in range(1, t_len):
            glBegin(GL_TRIANGLE_STRIP)
            for h in range(0, h_len):

                if (self.vertices[t][h][1] == None) or (self.vertices[t-1][h][1] == None):
                    glEnd()
                    glBegin(GL_TRIANGLE_STRIP)
                    continue
                
                glColor3f(*self.color_vertices[t][h])
                glVertex3f(*self.vertices[t][h]);

                glColor3f(*self.color_vertices[t-1][h])
                glVertex3f(*self.vertices[t-1][h])

            glEnd()

PlotFunctionRegistry.register('spherical', SphericalFunction)
