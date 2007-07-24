from pyglet.gl import *
from plot_function import PlotFunction, fsubs, vrange

class CartesianFunction2d(PlotFunction):
    
    def __init__(self, f, intervals, options):
        self.f = f
        self.intervals = intervals
        self.options = options
        self.calculate_vertices()
        self.calculate_color_vertices()
        
    def calculate_vertices(self):
        if len(self.intervals) != 1:
            raise NotImplementedError("Automatic intervals not implemented.")

        x, x_min, x_max, x_steps = self.intervals[0]
        if x_steps == None: x_steps = 60

        x_set = vrange(x_min, x_max, x_steps)
        self.vertices = list( (x_e, fsubs(self.f, x, x_e), 0.0) for x_e in x_set )

    def calculate_color_vertices(self):
        self.color_vertices = []
        for x, y, z in self.vertices:
            self.color_vertices.append( (0.7,0.75,0.8) )

    def render(self):
        glBegin(GL_LINE_STRIP)
        for x in range(0, len(self.vertices)):
            if self.vertices[x][1] == None:
                glEnd()
                glBegin(GL_LINE_STRIP)
                continue
            glColor3f(*self.color_vertices[x])
            glVertex3f(*self.vertices[x])
        glEnd()
