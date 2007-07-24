from pyglet.gl import *
from plot_object import PlotObject

class GridPlane(PlotObject):

    def __new__(cls, plane, pitch1=1.0, pitch2=1.0, res1=50, res2=50, origin=(0.0,0.0,-0.01), color=(0.9,0.9,0.9)):
        assert isinstance(plane, str)
        plane = plane.lower()

        if plane not in grid_plane_mapping:
            raise ValueError("'plane' kwarg must be one of %s." % (", ".join(["'"+str(i)+"'" for i in grid_plane_mapping])))

        subcls = grid_plane_mapping[plane]

        self = object.__new__(subcls)

        self.pitch1 = pitch1
        self.pitch2 = pitch2
        self.res1 = res1
        self.res2 = res2
        self.origin = origin
        self.color = color

        self.width1 = res1 * pitch1
        self.min1 = -(self.width1 / 2.0)
        self.max1 =  (self.width1 / 2.0)

        self.width2 = res2 * pitch2
        self.min2 = -(self.width2 / 2.0)
        self.max2 =  (self.width2 / 2.0)

        self.set1 = list( self.min1 + t * self.pitch1 for t in range(self.res1 + 1) )
        self.set2 = list( self.min2 + t * self.pitch2 for t in range(self.res2 + 1) )

        return self

    def order_verts(self, v1, v2):
        raise Exception("This method is abstract.")

    def transform(self, v1, v2):
        va = self.order_verts(v1, v2)
        vb = (va[0]+self.origin[0],
              va[1]+self.origin[1],
              va[2]+self.origin[2])
        return vb

    def render(self):
        glColor3f(*self.color)
        glBegin(GL_LINES)
        for v1 in self.set1:
            glVertex3f(*self.transform(v1, self.min2))
            glVertex3f(*self.transform(v1, self.max2))
        for v2 in self.set2:
            glVertex3f(*self.transform(self.min1, v2))
            glVertex3f(*self.transform(self.max1, v2))
        glEnd()

class GridPlaneXY(GridPlane):
    def order_verts(self, v1, v2):
        return (v1, v2, 0.0)

class GridPlaneXZ(GridPlane):
    def order_verts(self, v1, v2):
        return (v1, 0.0, v2)

class GridPlaneYZ(GridPlane):
    def order_verts(self, v1, v2):
        return (0.0, v1, v2)

grid_plane_mapping = dict(xy = GridPlaneXY, yx = GridPlaneXY,
                          xz = GridPlaneXZ, zx = GridPlaneXZ,
                          yz = GridPlaneYZ, zy = GridPlaneYZ)
