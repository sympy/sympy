from pyglet.gl import *

from plot_object import PlotObject
from plot_interval import PlotInterval

from util import find_bounds_2d, calc_color_1
from util import interpolate as intrp

from threading import Thread, Event

def uniform(r=0.3,g=0.3,b=0.7):
    return (lambda x,y,z,t: (r,g,b))

rainbow = (lambda x,y,z,t: (x*0.5+0.25, y*0.5+0.25, z*0.5+0.25))

class PlotCurve(PlotObject):
    """
    """
    def __init__(self, vertex_function, t_interval, **kwargs):
        self.vertex_function = vertex_function
        self.t_interval      = PlotInterval.try_parse(t_interval)

        self.vertices_ready = Event()
        self.colors_ready = Event()

        color = kwargs.pop('color', rainbow) # alias for color_function
        self.set_color_function(kwargs.pop('color_function', color))
        self.set_line_color(kwargs.pop('line_color', (0.7,0.7,0.7)))

        self.vertices = None
        self.colors = None

        self.vertices_display_list = -1
        self.colors_display_list = -1

        self.calculate()

    def convert_color(self, cs):
        if isinstance(cs, str):
            cs = eval(cs)
        if isinstance(cs, (list, tuple)) and len(cs) == 3:
            cs = list(cs)
            for e in range(len(cs)):
                cs[e] = float(cs[e])
            return tuple(cs)
        raise ValueError("Color must be in the form (r, g, b).")

    def set_color_function(self, cf):
        if isinstance(cf, str):
            cfstr = cf
            try:
                cf = eval(cfstr)
            except:
                try:
                    cf = eval("lambda x,y,z,t: (%s)" % (cfstr))
                    assert len(cf(0,0,0,0)) == 3
                except:
                    raise ValueError(("Color function string '%s' could not be "+
                                      "converted into a callable which returns "+
                                      "an RGB triplet.") % (cfstr))
        if cf is None or callable(cf):
            self.color_function = cf
        self.colors_ready.clear()
        if self.vertices_ready.isSet():
            Thread(target=self.calculate_colors).start()

    def set_line_color(self, cs):
        self.line_color = self.convert_color(cs)

    def set_style(self, rs):
        pass

    def get_calc_state(self):
        return (not self.vertices_ready.isSet(),
                not self.colors_ready.isSet())

    def calculate(self):
        Thread(target=self.calculation_thread).start()

    def calculation_thread(self):
        self.calculate_vertices()
        self.calculate_colors()

    def calculate_vertices(self):
        self.vertices = list((self.vertex_function(t), float(t)) for t in self.t_interval.vrange())
        self.vertices_ready.set()

    def calculate_colors(self):
        if self.color_function is not None and callable(self.color_function):
            self.pbounds = find_bounds_2d(self.vertices)
            self.ibounds = [[self.t_interval.v_min, self.t_interval.v_max]]
            self.colors  = list(calc_color_1(self.color_function,
                                             self.vertices[t][0], self.pbounds,
                                             self.vertices[t][1:], self.ibounds)
                                             for t in xrange(self.t_interval.v_len))
        self.colors_ready.set()

    def compile_display_list(self, use_custom_colors):
        display_list = glGenLists(1)

        glNewList(display_list, GL_COMPILE)
        glBegin(GL_LINE_STRIP)
        for t in xrange(self.t_interval.v_len):
            p = self.vertices[t][0]
            if p is None:
                glEnd()
                glBegin(GL_LINE_STRIP)
                continue
            if use_custom_colors:
                glColor3f(*self.colors[t])
            else:
                glColor3f(*self.line_color)
            glVertex3f(*p)
        glEnd()

        glEndList()
        return display_list

    def draw(self):
        if not self.colors_ready.isSet() and GL_TRUE == glIsList(self.colors_display_list):
            glDeleteLists(self.colors_display_list, 1)

        # Draw with colorization if possible.
        if self.color_function != None and self.colors_ready.isSet():
            if GL_FALSE == glIsList(self.colors_display_list):
                self.colors_display_list = self.compile_display_list(True)
            glCallList(self.colors_display_list)

        # Draw with default colors while waiting for
        # color calculations to finish.
        elif self.vertices_ready.isSet():
            if GL_FALSE == glIsList(self.vertices_display_list):
                self.vertices_display_list = self.compile_display_list(False)
            glCallList(self.vertices_display_list)
        
        else: pass # Can't draw anything yet
