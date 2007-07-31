from pyglet.gl import *
from managed_window import ManagedWindow

from plot_camera import PlotCamera
from plot_controller import PlotController

class PlotWindow(ManagedWindow):

    def __init__(self, plot, **kwargs):
        self.plot = plot

        self.camera = None
        self._calculating = False

        self.wireframe    = kwargs.pop('wireframe', False)
        self.antialiasing = kwargs.pop('antialiasing', True)
        self.ortho        = kwargs.pop('ortho', False)
        self.title = kwargs.setdefault('caption', "SymPy Plot")

        super(PlotWindow, self).__init__(**kwargs)

    def setup(self):
        self.camera = PlotCamera(self, ortho = self.ortho)
        self.controller = PlotController(self)
        self.push_handlers(self.controller)

        glClearColor(1.0, 1.0, 1.0, 0.0)
        glClearDepth(1.0)

        glDepthFunc(GL_LESS)
        glEnable(GL_DEPTH_TEST)

        glShadeModel(GL_SMOOTH)

        # now done at the function level
        # glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)

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
        
        calc_verts, calc_colors = 0, 0

        for r in self.plot._pobjects:
            glPushMatrix()
            r._draw()
            glPopMatrix()

        for r in self.plot._functions.itervalues():
            glPushMatrix()
            r._draw()
            try:
                r_calc_verts, r_calc_colors = r.get_calc_state()
                if r_calc_verts: calc_verts += 1
                if r_calc_colors: calc_colors += 1
            except: pass
            glPopMatrix()

        self.update_caption(calc_verts, calc_colors)
        self.plot._render_lock.release()

    def update_caption(self, calc_verts, calc_colors):
        caption = self.title

        if calc_verts > 0:
            caption = self.title + " (calculating vertices...)"
        elif calc_colors > 0:
            caption = self.title + " (calculating colors...)"

        if self.caption != caption:        
            self.set_caption(caption)
