from pyglet.gl import *
from plot_object import PlotObject
from util import strided_range
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

class PlotAxesOrdinate(PlotObject):

    def __init__(self, parent_axes):
        self._parent_axes = parent_axes

    def draw(self):
        p = self._parent_axes
        c = [ [0.2,0.1,0.3], [0.1,0.1,0.1] ][p._colored]

        self.draw_axe(p, 0, c)
        self.draw_axe(p, 1, c)
        self.draw_axe(p, 2, c)

    def draw_axe(self, p, axis, c):
        c = c[::]
        if p._colored: c[axis] = 1.0
        b = p._bounding_box
        a = p._axis_ticks[axis]
        r = p._tick_length / 2.0
        if not a or len(a) < 2: return
        self.draw_axe_line(axis, c, a[0], a[-1])
        if len(a)%2: del a[len(a)//2]
        for t in a: self.draw_tick_line(axis, c, r, t)

    def draw_axe_line(self, axis, c, a_min, a_max):
        v = [[0,0,0], [0,0,0]]
        v[0][axis] = a_min
        v[1][axis] = a_max
        self.draw_line(v, c)

    tick_axis = {0: 1, 1: 0, 2: 1}
    def draw_tick_line(self, axis, c, r, t, two_pronged=False):
        def d(a):
            v = [[0,0,0], [0,0,0]]
            v[0][axis] = v[1][axis] = t
            v[0][a], v[1][a] = -r, r
            self.draw_line(v, c)
        if two_pronged:
            alist = [0,1,2]
            alist.remove(axis)
            for a in alist: d(a)
        else:
            d(PlotAxesOrdinate.tick_axis[axis])

    def draw_line(self, v, c):
        o = self._parent_axes._origin
        glBegin(GL_LINES)
        glColor3f(*c)
        glVertex3f(v[0][0] + o[0],
                   v[0][1] + o[1],
                   v[0][2] + o[2])
        glVertex3f(v[1][0] + o[0],
                   v[1][1] + o[1],
                   v[1][2] + o[2])
        glEnd()

class PlotAxesFrame(PlotObject):

    def __init__(self, parent_axes):
        self._parent_axes = parent_axes

    def draw(self):
        raise NotImplementedError()

