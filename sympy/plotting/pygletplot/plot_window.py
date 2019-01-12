from __future__ import print_function, division

from pyglet.gl import *
from sympy.plotting.pygletplot.managed_window import ManagedWindow

from sympy.plotting.pygletplot.plot_camera import PlotCamera
from sympy.plotting.pygletplot.plot_controller import PlotController

from time import clock


class PlotWindow(ManagedWindow):

    def __init__(self, plot, antialiasing=True, ortho=False,
                 invert_mouse_zoom=False, linewidth=1.5, caption="SymPy Plot",
                 **kwargs):
        """
        Named Arguments
        ===============

        antialiasing = True
            True OR False
        ortho = False
            True OR False
        invert_mouse_zoom = False
            True OR False
        """
        self.plot = plot

        self.camera = None
        self._calculating = False

        self.antialiasing = antialiasing
        self.ortho = ortho
        self.invert_mouse_zoom = invert_mouse_zoom
        self.linewidth = linewidth
        self.title = caption
        self.last_caption_update = 0
        self.caption_update_interval = 0.2
        self.drawing_first_object = True

        super(PlotWindow, self).__init__(**kwargs)

    def setup(self):
        self.camera = PlotCamera(self, ortho=self.ortho)
        self.controller = PlotController(self,
                invert_mouse_zoom=self.invert_mouse_zoom)
        self.push_handlers(self.controller)

        glClearColor(1.0, 1.0, 1.0, 0.0)
        glClearDepth(1.0)

        glDepthFunc(GL_LESS)
        glEnable(GL_DEPTH_TEST)

        glEnable(GL_LINE_SMOOTH)
        glShadeModel(GL_SMOOTH)
        glLineWidth(self.linewidth)

        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

        if self.antialiasing:
            glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
            glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST)

        self.camera.setup_projection()

    def on_resize(self, w, h):
        super(PlotWindow, self).on_resize(w, h)
        if self.camera is not None:
            self.camera.setup_projection()

    def update(self, dt):
        self.controller.update(dt)

    def draw(self):
        self.plot._render_lock.acquire()
        self.camera.apply_transformation()

        calc_verts_pos, calc_verts_len = 0, 0
        calc_cverts_pos, calc_cverts_len = 0, 0

        should_update_caption = (clock() - self.last_caption_update >
                                 self.caption_update_interval)

        if len(self.plot._functions.values()) == 0:
            self.drawing_first_object = True

        try:
            dict.iteritems
        except AttributeError:
            # Python 3
            iterfunctions = iter(self.plot._functions.values())
        else:
            # Python 2
            iterfunctions = self.plot._functions.itervalues()

        for r in iterfunctions:
            if self.drawing_first_object:
                self.camera.set_rot_preset(r.default_rot_preset)
                self.drawing_first_object = False

            glPushMatrix()
            r._draw()
            glPopMatrix()

            # might as well do this while we are
            # iterating and have the lock rather
            # than locking and iterating twice
            # per frame:

            if should_update_caption:
                try:
                    if r.calculating_verts:
                        calc_verts_pos += r.calculating_verts_pos
                        calc_verts_len += r.calculating_verts_len
                    if r.calculating_cverts:
                        calc_cverts_pos += r.calculating_cverts_pos
                        calc_cverts_len += r.calculating_cverts_len
                except ValueError:
                    pass

        for r in self.plot._pobjects:
            glPushMatrix()
            r._draw()
            glPopMatrix()

        if should_update_caption:
            self.update_caption(calc_verts_pos, calc_verts_len,
                                calc_cverts_pos, calc_cverts_len)
            self.last_caption_update = clock()

        if self.plot._screenshot:
            self.plot._screenshot._execute_saving()

        self.plot._render_lock.release()

    def update_caption(self, calc_verts_pos, calc_verts_len,
            calc_cverts_pos, calc_cverts_len):
        caption = self.title
        if calc_verts_len or calc_cverts_len:
            caption += " (calculating"
            if calc_verts_len > 0:
                p = (calc_verts_pos / calc_verts_len) * 100
                caption += " vertices %i%%" % (p)
            if calc_cverts_len > 0:
                p = (calc_cverts_pos / calc_cverts_len) * 100
                caption += " colors %i%%" % (p)
            caption += ")"
        if self.caption != caption:
            self.set_caption(caption)
