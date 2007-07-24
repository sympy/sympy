from pyglet.gl import *
from plot_function import PlotFunction, vrange, fsubs, rinterpolate, interpolate

from sympy import Basic
from math import tan

class CartesianFunction3d(PlotFunction):
    
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

        x, self.x_min, self.x_max, x_steps = self.intervals[0]
        if x_steps == None: x_steps = 12
        
        y, self.y_min, self.y_max, y_steps = self.intervals[1]
        if y_steps == None: y_steps = 12
        
        x_set = vrange(self.x_min, self.x_max, x_steps)
        y_set = vrange(self.y_min, self.y_max, y_steps)

        self.vertices = list( list( (_x, _y, fsubs(self.f, x, _x, y, _y))
                                    for _y in y_set ) for _x in x_set )

    def calculate_bounding_box(self):
        x_len = len(self.vertices)
        y_len = len(self.vertices[0])

        for x in range(x_len):
            for y in range(y_len):
                self.x_min = min([self.x_min, self.vertices[x][y][0]])
                self.x_max = max([self.x_max, self.vertices[x][y][0]])

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
        x_len = len(self.vertices)
        y_len = len(self.vertices[0])

        for x in range(1, x_len):
            glBegin(GL_TRIANGLE_STRIP)
            for y in range(0, y_len):

                if (self.vertices[x][y][2] == None) or (self.vertices[x-1][y][2] == None):
                    glEnd()
                    glBegin(GL_TRIANGLE_STRIP)
                    continue

                glColor3f(*self.color_vertices[x][y])
                glVertex3f(*self.vertices[x][y])

                glColor3f(*self.color_vertices[x-1][y])
                glVertex3f(*self.vertices[x-1][y])

            glEnd()
