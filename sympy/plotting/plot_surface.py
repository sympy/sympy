from pyglet.gl import *

from plot_object import PlotObject
from plot_interval import PlotInterval

from util import find_bounds_3d, calc_color_2
from util import interpolate as intrp

from threading import Thread, Event

WIREFRAME = 1
SOLID = 2
SOLID_AND_WIREFRAME = WIREFRAME | SOLID

rainbow = (lambda x,y,z,u,v: (intrp(0.6, 0.97, x), intrp(0.6, 0.97, y), intrp(0.6, 0.97, z)))

class PlotSurface(PlotObject):
    """
    """
    def __init__(self, vertex_function, u_interval, v_interval, **kwargs):
        self.vertex_function = vertex_function
        self.u_interval      = PlotInterval.try_parse(u_interval)
        self.v_interval      = PlotInterval.try_parse(v_interval)

        self.vertices_ready = Event()
        self.colors_ready = Event()

        color = kwargs.pop('color', rainbow) # alias for color_function
        self.set_color_function(kwargs.pop('color_function', color))

        self.set_line_color(kwargs.pop('line_color', (0.7,0.7,0.7)))
        self.set_fill_color(kwargs.pop('fill_color', (0.5,0.5,0.9)))
        self.set_style(kwargs.pop('style', SOLID_AND_WIREFRAME))

        self.vertices = None
        self.colors = None

        self.vertices_display_list_fill = -1
        self.vertices_display_list_line = -1
        self.colors_display_list_fill = -1

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

    def set_line_color(self, cs):
        self.line_color = self.convert_color(cs)
        
    def set_fill_color(self, cs):
        self.fill_color = self.convert_color(cs)

    def set_style(self, rs):
        if isinstance(rs, int):
            self.style = rs
        elif isinstance(rs, str):
            rs = rs.lower()
            wireframe = ['wire', 'wireframe', 'line']
            solid = ['solid', 'fill']
            both = ['both']
            both += ["%s and %s" % (w,s) for w in wireframe for s in solid]
            both += ["%s and %s" % (s,w) for w in wireframe for s in solid]
            both += ["%s_and_%s" % (w,s) for w in wireframe for s in solid]
            both += ["%s_and_%s" % (s,w) for w in wireframe for s in solid]

            if rs in wireframe:
                self.style = WIREFRAME
            elif rs in solid:
                self.style = SOLID
            elif rs in both:
                self.style = SOLID_AND_WIREFRAME
            else:
                try: self.style = int(eval(rs))
                except: pass

    def set_color_function(self, cf):
        if isinstance(cf, str):
            cfstr = cf
            try:
                cf = eval(cfstr)
            except:
                try:
                    cf = eval("lambda x,y,z,u,v: (%s)" % (cfstr))
                    assert len(cf(0,0,0,0,0)) == 3
                except:
                    raise ValueError(("Color function string '%s' could not be "+
                                      "converted into a callable which returns "+
                                      "an RGB triplet.") % (cfstr))
        if cf is None or callable(cf):
            self.color_function = cf
        self.colors_ready.clear()
        if self.vertices_ready.isSet():
            Thread(target=self.calculate_colors).start()

    def get_calc_state(self):
        return (not self.vertices_ready.isSet(),
                not self.colors_ready.isSet())

    def calculate(self):
        Thread(target=self.calculation_thread).start()

    def calculation_thread(self):
        self.calculate_vertices()
        self.calculate_colors()

    def calculate_vertices(self):
        v_set = list(float(v) for v in self.v_interval.vrange())
        u_set = list(float(u) for u in self.u_interval.vrange())
        self.vertices = list(list((self.vertex_function(u,v), float(u), float(v)) for v in v_set) for u in u_set)
        self.vertices_ready.set()

    def calculate_colors(self):
        if self.color_function is not None and callable(self.color_function):
            self.pbounds = find_bounds_3d(self.vertices)
            self.ibounds = [[self.u_interval.v_min, self.u_interval.v_max],
                            [self.v_interval.v_min, self.v_interval.v_max]]
            self.colors  = list(list(calc_color_2(self.color_function,
                                                  self.vertices[u][v][0], self.pbounds,
                                                  self.vertices[u][v][1:], self.ibounds)
                                                  for v in xrange(self.v_interval.v_len))
                                                  for u in xrange(self.u_interval.v_len))
        self.colors_ready.set()

    def compile_display_list(self, use_custom_colors, fill):
        display_list = glGenLists(1)
        glNewList(display_list, GL_COMPILE)

        for u in xrange(1, self.u_interval.v_len):
            glBegin(GL_QUAD_STRIP)
            for v in xrange(self.v_interval.v_len):
                pa = self.vertices[u-1][v][0]
                pb = self.vertices[u][v][0]
                if pa is None or pb is None:
                    glEnd()
                    glBegin(GL_QUAD_STRIP)
                    continue
                if use_custom_colors:
                    ca = self.colors[u-1][v]
                    cb = self.colors[u][v]
                else:
                    if fill: ca = cb = self.fill_color
                    else:    ca = cb = self.line_color
                glColor3f(*ca)
                glVertex3f(*pa)
                glColor3f(*cb)
                glVertex3f(*pb)
            glEnd()

        glEndList()
        return display_list

    def draw_line_color(self):
        # Always uses line color
        if self.vertices_ready.isSet():
            if GL_FALSE == glIsList(self.vertices_display_list_line):
                self.vertices_display_list_line = self.compile_display_list(False, False)
            glCallList(self.vertices_display_list_line)

    def draw_fill_color(self):
        # Draw with colorization if possible.
        if not self.colors_ready.isSet() and GL_TRUE == glIsList(self.colors_display_list_fill):
            glDeleteLists(self.colors_display_list_fill, 1)

        if self.color_function != None and self.colors_ready.isSet():
            if GL_FALSE == glIsList(self.colors_display_list_fill):
                self.colors_display_list_fill = self.compile_display_list(True, True)
            glCallList(self.colors_display_list_fill)

        # Draw with default colors while waiting for
        # color calculations to finish.
        elif self.vertices_ready.isSet():
            if GL_FALSE == glIsList(self.vertices_display_list_fill):
                self.vertices_display_list_fill = self.compile_display_list(False, True)
            glCallList(self.vertices_display_list_fill)
        
        else: pass # Can't draw anything yet

    def draw(self):
        glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT)
        
        if self.style == SOLID:
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
            self.draw_fill_color()

        elif self.style == WIREFRAME:
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
            self.draw_fill_color()

        elif self.style == SOLID_AND_WIREFRAME:
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
            self.draw_fill_color()
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
            glEnable(GL_POLYGON_OFFSET_LINE)
            glPolygonOffset(-0.005, -50.0)
            self.draw_line_color()
            #glDisable(GL_POLYGON_OFFSET_LINE)

        glPopAttrib()
