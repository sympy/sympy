from pyglet.gl import *
from pyglet import font

from plot_object import PlotObject
from util import strided_range, billboard_matrix
from util import get_direction_vectors
from util import dot_product, vec_sub, vec_mag
from sympy import oo

class PlotAxes(PlotObject):

    def __init__(self, *args, **kwargs):
        # initialize style parameter
        style = kwargs.pop('style', '').lower()
        # allow alias kwargs to override style kwarg
        if kwargs.pop('none', None) is not None: style = 'none'
        if kwargs.pop('frame', None) is not None: style = 'frame'
        if kwargs.pop('box', None) is not None: style = 'box'
        if kwargs.pop('ordinate', None) is not None: style = 'ordinate'
        if style in ['', 'ordinate']:
            self._render_object = PlotAxesOrdinate(self)
        elif style in ['frame', 'box']:
            self._render_object = PlotAxesFrame(self)
        elif style in ['none']:
            self._render_object = None
        else: raise ValueError(("Unrecognized axes "
                                "style %s.") % (style))

        # initialize stride parameter
        stride = kwargs.pop('stride', 0.25)
        try: stride = eval(stride)
        except: pass
        if isinstance(stride, (list, tuple)):
            assert len(stride) == 3
            self._stride = stride
        else:
            self._stride = [stride, stride, stride]
        self._tick_length = kwargs.pop('tick_length', 0.1)

        # setup bounding box and ticks
        self._origin = [0,0,0]
        self.reset_bounding_box()

        def flexible_boolean(input, default):
            if input in [True, False]:
                return input
            if input in ['f','F','false','False']: return False
            if input in ['t','T','true','True']: return True
            return default

        # initialize remaining parameters
        self._overlay = flexible_boolean(kwargs.pop('overlay',''), True)
        self._colored = flexible_boolean(kwargs.pop('colored',''), False)
        self._label_axes = flexible_boolean(kwargs.pop('label_axes', ''), False)
        self._label_ticks = flexible_boolean(kwargs.pop('label_ticks', ''), True)

        # setup label font
        self.font_face = kwargs.pop('font_face', 'Arial')
        self.font_size = kwargs.pop('font_size', 28)

        self.reset_resources()

    def reset_resources(self):
        self.label_font = None

    def reset_bounding_box(self):
        self._bounding_box = [[None,None], [None,None], [None,None]]
        self._axis_ticks = [[],[],[]]

    def draw(self):
        if self._render_object:
            glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT | GL_DEPTH_BUFFER_BIT)
            if self._overlay:
                glDisable(GL_DEPTH_TEST)
            self._render_object.draw()
            glPopAttrib()

    def adjust_bounds(self, child_bounds):
        b = self._bounding_box
        c = child_bounds
        for i in [0,1,2]:
            if abs(c[i][0]) == oo or abs(c[i][1]) == oo: continue
            b[i][0] = [ min([b[i][0], c[i][0]]), c[i][0] ][ b[i][0] is None ]
            b[i][1] = [ max([b[i][1], c[i][1]]), c[i][1] ][ b[i][1] is None ]
            self._adjust_axis_ticks(i)

    def _adjust_axis_ticks(self, axis):
        b = self._bounding_box
        if b[axis][0] is None or b[axis][1] is None:
            self._axis_ticks[axis] = []
        else:
            self._axis_ticks[axis] = strided_range(b[axis][0], b[axis][1], self._stride[axis])

class PlotAxesBase(PlotObject):

    def draw_text(self, text, position, c, scale=1.0):
        if len(c) == 3: c = [c[0], c[1], c[2], 1.0]
        p = self._parent_axes
        if p.label_font is None:
            p.label_font = font.load(p.font_face, p.font_size, bold=True, italic=False)
        label = font.Text(p.label_font, text, color=c, valign=font.Text.BASELINE, halign=font.Text.CENTER)

        glPushMatrix()
        glTranslatef(*position)
        billboard_matrix()
        scalef = 0.005*scale
        glScalef(scalef, scalef, scalef)
        glColor4f(0, 0, 0, 0)
        label.draw()
        glPopMatrix()

    def draw_line(self, v, c):
        p = self._parent_axes
        o = p._origin
        glBegin(GL_LINES)
        glColor3f(*c)
        glVertex3f(v[0][0] + o[0],
                   v[0][1] + o[1],
                   v[0][2] + o[2])
        glVertex3f(v[1][0] + o[0],
                   v[1][1] + o[1],
                   v[1][2] + o[2])
        glEnd()

class PlotAxesOrdinate(PlotAxesBase):

    def __init__(self, parent_axes):
        self._parent_axes = parent_axes

    def draw(self):
        p = self._parent_axes
        #c = [ [0.2,0.1,0.3], ([0.9,0.3,0.5], [0.5,0.8,0.5], [0.3,0.3,0.9]) ][p._colored]
        c = [ [0.2,0.1,0.3], ([0.9,0.3,0.5], [0.5,1.0,0.5], [0.3,0.3,0.9]) ][p._colored]

        self.draw_axe(p, 0, c)
        self.draw_axe(p, 1, c)
        self.draw_axe(p, 2, c)

    def draw_axe(self, p, axis, c):
        if p._colored: c = c[axis]
        b = p._bounding_box
        a = p._axis_ticks[axis]
        r = p._tick_length / 2.0
        if not a or len(a) < 2: return
        v = [[0,0,0], [0,0,0]]
        v[0][axis] = a[0]
        v[1][axis] = a[-1]
        v = vec_sub(v[1], v[0])
        d = abs(dot_product(v, get_direction_vectors()[2]))
        labels_visible = abs(d/vec_mag(v) - 1) > 0.02
        self.draw_axe_line(axis, c, a[0], a[-1], labels_visible)
        for t in a:
            self.draw_tick_line(axis, c, r, t, labels_visible)

    def draw_axe_line(self, axis, c, a_min, a_max, labels_visible):
        v = [[0,0,0], [0,0,0]]
        v[0][axis] = a_min
        v[1][axis] = a_max
        self.draw_line(v, c)
        if labels_visible:
            self.draw_axe_line_labels(axis, c, a_min, a_max)

    def draw_axe_line_labels(self, axis, c, a_min, a_max):
        if not self._parent_axes._label_axes: return
        v = [[0,0,0], [0,0,0]]
        v[0][axis] = a_min-0.35
        v[1][axis] = a_max+0.35
        a_str = ['X', 'Y', 'Z'][axis]
        max_str = "+" + a_str
        min_str = "-" + a_str
        self.draw_text(min_str, v[0], c)
        self.draw_text(max_str, v[1], c)

    tick_axis = {0: 1, 1: 0, 2: 1}
    def draw_tick_line(self, axis, c, r, t, labels_visible, two_pronged=False):
        def d(a):
            v = [[0,0,0], [0,0,0]]
            v[0][axis] = v[1][axis] = t
            v[0][a], v[1][a] = -r, r
            self.draw_line(v, c)
        if two_pronged:
            alist = [0,1,2]
            alist.remove(axis)
            for a in alist: d(a)
        else: d(PlotAxesOrdinate.tick_axis[axis])
        if labels_visible:
            self.draw_tick_line_label(axis, c, r, t)
            
    def draw_tick_line_label(self, axis, c, r, t):
        if not self._parent_axes._label_axes: return
        v = [0,0,0]
        v[axis] = t
        v[PlotAxesOrdinate.tick_axis[axis]] = [-1,1,1][axis]*r*3.5
        self.draw_text(str(t), v, c, scale=0.5)

class PlotAxesFrame(PlotAxesBase):

    def __init__(self, parent_axes):
        self._parent_axes = parent_axes

    def draw(self):
        raise NotImplementedError()

