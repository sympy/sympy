from pyglet.gl import *
from plot_object import PlotObject
from plot_function import PlotFunction
#from pyglet import font
from util import frange, billboard_matrix, get_basis_vectors, get_view_direction_vectors, invert_vector, model_to_screen_ratio
from math import sqrt

class BoundingBox(PlotObject):

    def __init__(self,
                 line_color=(0.7,0.75,0.8),
                 fill_color=(0.97,0.97,0.98),
                 show_axes=False,
                 x_stride=0.5,
                 y_stride=0.5,
                 z_stride=0.5,
                 #font_face="Courier New",
                 #font_size=14,
                 cull_line_front=True,
                 default_radius=0.0):

        d = default_radius/2.0
        self.x_min, self.x_max = -d, d
        self.y_min, self.y_max = -d, d
        self.z_min, self.z_max = -d, d

        self.line_color = line_color
        self.fill_color = fill_color
        self.cull_line_front = cull_line_front

        self.show_axes = show_axes
        self.stride = [x_stride, y_stride, z_stride]
        #self.font_face = font_face
        #self.font_size = font_size
        #self.label_font = None

    def consider_function(self, f):
        def pad_stride(v, axis, isMin):
            #s = self.stride[axis]
            #r = (v % s)
            #if isMin:
            #    return v - r - s
            #return v + 2*s - r
            return v
        self.x_min = min([self.x_min, pad_stride(f.x_min, 0, True)])
        self.x_max = max([self.x_max, pad_stride(f.x_max, 0, False)])
        self.y_min = min([self.y_min, pad_stride(f.y_min, 1, True)])
        self.y_max = max([self.y_max, pad_stride(f.y_max, 1, False)])
        self.z_min = min([self.z_min, pad_stride(f.z_min, 2, True)])
        self.z_max = max([self.z_max, pad_stride(f.z_max, 2, False)])

    def render(self):
        if self.x_min == self.x_max and self.y_min == self.y_max and self.z_min == self.z_max:
            return
        #if self.show_axes:
        #    self.render_axes()
        if self.fill_color is not None:
            self.render_box(line=False)
        if self.line_color is not None:
            self.render_box(line=True)

    def render_direction_vectors(self):
        x, y, z = get_direction_vectors()

        glBegin(GL_LINES)

        glColor3f(1,0,0)
        glVertex3f(0,0,0)
        glVertex3f(*x)

        glColor3f(0,1,0)
        glVertex3f(0,0,0)
        glVertex3f(*y)

        glColor3f(0,0,1)
        glVertex3f(0,0,0)
        glVertex3f(*z)

        glEnd()

    def render_basis_vectors(self):
        x, y, z = get_basis_vectors()

        glBegin(GL_LINES)

        glColor3f(1,0,0)
        glVertex3f(0,0,0)
        glVertex3f(*x)

        glColor3f(0,1,0)
        glVertex3f(0,0,0)
        glVertex3f(*y)

        glColor3f(0,0,1)
        glVertex3f(0,0,0)
        glVertex3f(*z)

        glEnd()

    """def render_axes(self):
        glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT | GL_DEPTH_BUFFER_BIT)
        glDisable(GL_DEPTH_TEST)

        p = 0.2 # padding
        dx, dy, dz = get_view_direction_vectors()
        ax = [self.x_min, self.x_max]
        ay = [self.y_min, self.y_max]
        az = [self.z_min, self.z_max]
        bx = [self.x_min-p, self.x_max+p]
        by = [self.y_min-p, self.y_max+p]
        bz = [self.z_min-p, self.z_max+p]
        self.render_x_axis(ax,ay,az, bx,by,bz, dx,dy,dz)
        self.render_y_axis(ax,ay,az, bx,by,bz, dx,dy,dz)
        self.render_z_axis(ax,ay,az, bx,by,bz, dx,dy,dz)

        glPopAttrib()

    def render_x_axis(self, ax,ay,az, bx,by,bz, dx,dy,dz):
        cy = by[int(dy[2]<0)]
        cz = bz[int(dz[2]>0)]
        self.render_axis_labels(0, (ax[0], cy, cz), (ax[1], cy, cz))

    def render_y_axis(self, ax,ay,az, bx,by,bz, dx,dy,dz):
        cx = bx[int(dx[2]>0)]
        cz = bz[int(dz[2]<0)]
        self.render_axis_labels(1, (cx, ay[0], cz), (cx, ay[1], cz))

    def render_z_axis(self, ax,ay,az, bx,by,bz, dx,dy,dz):
        cx = bx[int(dx[2]>0)]
        cy = by[int(dy[2]<0)]
        self.render_axis_labels(2, (cx, cy, az[0]), (cx, cy, az[1]))

    def render_axis_labels(self, axis, start, end):
        color = [0.25, 0.25, 0.25, 1.0]
        color[axis] = 0.75
        color = tuple(color)
        mag = end[axis]-start[axis]
        steps = mag/self.stride[axis]

        def vstr(val):
            vs = str(val)
            if len(vs) > 4:
                vs = vs[0:5]
            return vs

        #self.draw_text(vstr(start[axis]), start, color)
        for val in frange(start[axis]+self.stride[axis],
                          end[axis], self.stride[axis]):
            vec = list(start)
            vec[axis] = val
            vec = tuple(vec)
            self.draw_text(vstr(val), vec, color)
        #self.draw_text(vstr(end[axis]), end, color)

    def draw_text(self, text, position, color):
        if self.label_font is None:
            self.label_font = font.load(self.font_face, self.font_size,
                                        bold=True, italic=False)
        label = font.Text(self.label_font, text, color=color,
                          valign=font.Text.TOP, halign=font.Text.CENTER)

        glPushMatrix()
        glTranslatef(*position)
        billboard_matrix()
        glScalef(0.007, 0.007, 0.007)
        glColor4f(0, 0, 0, 0)
        label.draw()
        glPopMatrix()"""

    def render_box(self, line=True):
        glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT | GL_DEPTH_BUFFER_BIT)

        if not line or self.cull_line_front:
            glCullFace(GL_FRONT)
            glEnable(GL_CULL_FACE)
            glDisable(GL_DEPTH_TEST)
        if line:
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
            glColor3f(*self.line_color)
        else:
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
            glColor3f(*self.fill_color)

        glBegin(GL_QUADS)

        # Top
        glVertex3f(self.x_max, self.y_max, self.z_max)
        glVertex3f(self.x_max, self.y_max, self.z_min)
        glVertex3f(self.x_min, self.y_max, self.z_min)
        glVertex3f(self.x_min, self.y_max, self.z_max)

        # Bottom
        glVertex3f(self.x_min, self.y_min, self.z_max)
        glVertex3f(self.x_min, self.y_min, self.z_min)
        glVertex3f(self.x_max, self.y_min, self.z_min)
        glVertex3f(self.x_max, self.y_min, self.z_max)

        # Left
        glVertex3f(self.x_min, self.y_max, self.z_min)
        glVertex3f(self.x_min, self.y_min, self.z_min)
        glVertex3f(self.x_min, self.y_min, self.z_max)
        glVertex3f(self.x_min, self.y_max, self.z_max)

        # Right
        glVertex3f(self.x_max, self.y_max, self.z_max)
        glVertex3f(self.x_max, self.y_min, self.z_max)
        glVertex3f(self.x_max, self.y_min, self.z_min)
        glVertex3f(self.x_max, self.y_max, self.z_min)

        # Front
        glVertex3f(self.x_min, self.y_max, self.z_max)
        glVertex3f(self.x_min, self.y_min, self.z_max)
        glVertex3f(self.x_max, self.y_min, self.z_max)
        glVertex3f(self.x_max, self.y_max, self.z_max)

        # Back
        glVertex3f(self.x_max, self.y_max, self.z_min)
        glVertex3f(self.x_max, self.y_min, self.z_min)
        glVertex3f(self.x_min, self.y_min, self.z_min)
        glVertex3f(self.x_min, self.y_max, self.z_min)

        glEnd()

        glPopAttrib()
