from pyglet.gl import *
from plot_mode_base import PlotModeBase
from sympy import oo
from util import scale_value, scale_value_list

class PlotCurve(PlotModeBase):

    style_override = 'wireframe'

    def _on_calculate_verts(self):
        self.t_interval = self.intervals[0]
        self.t_set = list(self.t_interval.frange())
        self.bounds = [ [oo,-oo,0],[oo,-oo,0],[oo,-oo,0] ]
        evaluate = self._get_evaluator()

        self._calculating_verts_pos = 0.0
        self._calculating_verts_len = float(self.t_interval.v_len)

        self.verts = list()
        b = self.bounds
        for t in self.t_set:
            try: _e = evaluate(t)   # calculate vertex
            except: _e = None
            if _e is not None: # update bounding box
                for axis in xrange(3):
                    b[axis][0] = min([b[axis][0], _e[axis]])
                    b[axis][1] = max([b[axis][1], _e[axis]])
            self.verts.append(_e)
            self._calculating_verts_pos += 1.0

        for axis in xrange(3):
            b[axis][2] = b[axis][1] - b[axis][0]
            if b[axis][2] == 0.0: b[axis][2] = 1.0

        self.push_wireframe(self.draw_verts(False))

    def _on_calculate_cverts(self):
        if not self.verts or not self.color: return
        self.st_set = scale_value_list(self.t_set)

        self._calculating_cverts_pos = 0.0
        self._calculating_cverts_len = float(self.t_interval.v_len)

        self.cverts = list()
        for t in xrange(self.t_interval.v_len):
            if self.verts[t] is None: _c = (0,0,0)
            else: _c = self.calculate_one_cvert(t)
            self.cverts.append(_c)
            self._calculating_cverts_pos += 1.0

        self.push_wireframe(self.draw_verts(True))

    def calculate_one_cvert(self, t):
        vert = self.verts[t]
        return self.color(scale_value(vert[0], self.bounds[0][0], self.bounds[0][2]),
                          scale_value(vert[1], self.bounds[1][0], self.bounds[1][2]),
                          scale_value(vert[2], self.bounds[2][0], self.bounds[2][2]),
                          self.st_set[t], None)

    def draw_verts(self, use_cverts):
        def f():
            glBegin(GL_LINE_STRIP)
            for t in xrange( len(self.t_set) ):
                p = self.verts[t]
                if p is None:
                    glEnd()
                    glBegin(GL_LINE_STRIP)
                    continue
                if use_cverts: glColor3f(*self.cverts[t])
                else: glColor3f(*self.default_wireframe_color)
                glVertex3f(*p)
            glEnd()
        return f
