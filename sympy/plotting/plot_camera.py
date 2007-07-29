from pyglet.gl import *
from plot_rotation import get_spherical_rotatation
from util import get_model_matrix

class PlotCamera(object):

    min_dist = 1.0
    max_dist = 500.0

    min_ortho_dist = 100.0
    max_ortho_dist = 10000.0

    _default_dist = 12.0
    _default_ortho_dist = 1200.0

    def __init__(self, window, ortho = False):
        self._dist = 0.0
        self._x, self._y = 0.0, 0.0
        self._rot = None

        self.window = window
        self.ortho = ortho
        if self.ortho:
            self._dist = self._default_ortho_dist
        else:
            self._dist = self._default_dist
        self.init_rot_matrix()

    def init_rot_matrix(self):
        glPushMatrix()
        glLoadIdentity()
        self._rot = get_model_matrix()
        glPopMatrix()

    def mult_rot_matrix(self, rot):
        glPushMatrix()
        glLoadMatrixf(rot)
        glMultMatrixf(self._rot)
        self._rot = get_model_matrix()
        glPopMatrix()

    def setup_projection(self):
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        if self.ortho:
            # yep, this is pseudo ortho (don't tell anyone)
            gluPerspective(0.3, float(self.window.width) / float(self.window.height), self.min_ortho_dist-0.5, self.max_ortho_dist+0.5)
        else:
            gluPerspective(30.0, float(self.window.width) / float(self.window.height), self.min_dist-0.5, self.max_dist+0.5)
        glMatrixMode(GL_MODELVIEW)

    def apply_transformation(self):
        glLoadIdentity()
        glTranslatef(self._x, self._y, -self._dist)
        if self._rot is not None:
            glMultMatrixf(self._rot)

    def spherical_rotate(self, p1, p2, sensitivity=1.0):
        mat = get_spherical_rotatation(p1, p2, self.window.width, self.window.height, sensitivity)
        if mat is not None: self.mult_rot_matrix(mat)

    def euler_rotate(self, angle, x, y, z):
        glPushMatrix()
        glLoadMatrixf(self._rot)
        glRotatef(angle, x, y, z)
        self._rot = get_model_matrix()
        glPopMatrix()

    def zoom_relative(self, clicks, sensitivity):

        if self.ortho:
            dist_d = clicks * sensitivity * 50.0
            min_dist = self.min_ortho_dist
            max_dist = self.max_ortho_dist
        else:
            dist_d = clicks * sensitivity
            min_dist = self.min_dist
            max_dist = self.max_dist
        
        new_dist = (self._dist - dist_d)
        if (clicks < 0 and new_dist < max_dist) or new_dist > min_dist:
            self._dist = new_dist

    def translate(self, dx, dy, sensitivity):
        self._x += dx*sensitivity/200.0
        self._y += dy*sensitivity/200.0
