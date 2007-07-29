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

        self.setup_polygon_mode()

        if self.antialiasing:
            glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
            glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST)

        self.camera.setup_projection()

    def on_resize(self, w, h):
        super(PlotWindow, self).on_resize(w, h)
        if self.camera is not None:
            self.camera.setup_projection()

    def setup_polygon_mode(self):
        if self.wireframe:
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
        else:
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)

    def update(self, dt):
        self.update_caption()
        self.controller.update(dt)

    def draw(self):
        self.plot._render_lock.acquire()

        self.camera.apply_transformation()

        for r in self.plot._pobjects:
            glPushMatrix()
            r._draw()
            glPopMatrix()

        for r in self.plot._functions.itervalues():
            glPushMatrix()
            r._draw()
            glPopMatrix()

        self.plot._render_lock.release()

    def update_caption(self):
        if not self._calculating and self.plot._calculations_in_progress > 0:
            self.set_caption(self.title + " (calculating...)")
            self._calculating = True
        elif self._calculating and self.plot._calculations_in_progress == 0:
            self.set_caption(self.title)
            self._calculating = False
