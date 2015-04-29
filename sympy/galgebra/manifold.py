# sympy/galgebra/manifold.py

"""
manifold.py defines the Manifold class which allows one to create a
vector manifold (manifold defined by vector field of coordinates in
embedding vector space) calculate the tangent vectors and derivatives
of tangent vectors.

Once manifold is created multivector fields can be constructed in the
tangent space and all the geometric algebra products and derivatives
of the multivector fields calculated.

Note that all calculations are done in the embedding space.  Future
versions of the code will allow manifolds defined purely in terms of
a metric.
"""

from __future__ import print_function

from itertools import combinations
from os import system
import copy

from sympy import trigsimp, simplify
from sympy.core.compatibility import range

from sympy.galgebra.ga import MV
from sympy.galgebra.debug import oprint
from sympy.galgebra.ncutil import linear_expand
from sympy.galgebra.printing import find_executable


def fct_to_str(fct_names):
    import sys
    current_file = open(sys.argv[0], 'r')
    file_str = current_file.read()
    current_file.close()

    if isinstance(fct_names, str):
        return fct_names

    fcts_str = ''

    for fct_name in fct_names:
        start_def = file_str.find('\ndef ' + fct_name)
        end_def = file_str.find('\ndef ', start_def + 5)
        start_class = file_str.find('\nclass ', start_def + 5)
        end_def = min(end_def, start_class)
        fcts_str += file_str[start_def:end_def]
    return fcts_str


def VectorComponents(X, basis):
    (coefs, bases) = linear_expand(X.obj)
    cdict = {}
    for (coef, base) in zip(coefs, bases):
        cdict[str(base)] = coef
    comp = []
    for base in basis:
        if base in cdict:
            comp.append(cdict[base])
        else:
            comp.append(0)
    return comp


def FillTemplate(self, template):
    Nd = 0
    var = []
    id_old = 0
    while True:
        id_new = template.find('$', id_old + 1)
        if id_new == -1:
            break
        Nd += 1
        if Nd % 2 == 0:
            var.append(template[id_old + 1:id_new])
        id_old = id_new

    var.sort(reverse=True)

    for v in var:
        template = template.replace('$' + v + '$', str(eval('self.' + v)))

    return template


class Manifold:

    def __init__(self, x, coords, debug=False, I=None):
        """
        coords: list of coordinate variables
        x: vector fuction of coordinate variables (parametric surface)
        """
        self.I = I
        self.x = x
        self.coords = coords

        self.basis = []
        self.basis_str = []
        self.embedded_basis = []
        for u in coords:
            tv = x.diff(u)
            self.basis.append(tv)
            (coefs, bases) = linear_expand(tv.obj)
            tc = {}
            for (coef, base) in zip(coefs, bases):
                str_base = str(base)
                tc[str_base] = coef
                if str_base not in self.embedded_basis:
                    self.embedded_basis.append(str_base)
            self.basis_str.append(tc)

        self.gij = []

        for base1 in self.basis:
            tmp = []
            for base2 in self.basis:
                tmp.append(simplify(trigsimp((base1 | base2).scalar())))
            self.gij.append(tmp)

        for tv in self.basis_str:
            for base in self.embedded_basis:
                if base not in tv:
                    tv[base] = 0

        self.dim = len(self.basis)

        indexes = tuple(range(self.dim))
        self.index = [()]
        for i in indexes:
            self.index.append(tuple(combinations(indexes, i + 1)))
        self.index = tuple(self.index)

        self.MFbasis = [[MV.ONE], self.basis]

        for igrade in self.index[2:]:
            grade = []
            for iblade in igrade:
                blade = MV(1, 'scalar')
                for ibasis in iblade:
                    blade ^= self.basis[ibasis]
                blade = blade.trigsimp(deep=True, recursive=True)
                grade.append(blade)
            self.MFbasis.append(grade)
        self.E = self.MFbasis[-1][0]
        self.E_sq = trigsimp((self.E * self.E).scalar(), deep=True, recursive=True)

        duals = copy.copy(self.MFbasis[-2])

        duals.reverse()
        sgn = 1
        self.rbasis = []
        for dual in duals:
            recpv = (sgn * dual * self.E).trigsimp(deep=True, recursive=True)
            self.rbasis.append(recpv)
            sgn = -sgn

        self.dbasis = []

        for base in self.basis:
            dbase = []
            for coord in self.coords:
                d = base.diff(coord).trigsimp(deep=True, recursive=True)
                dbase.append(d)
            self.dbasis.append(dbase)

        self.surface = {}
        (coefs, bases) = linear_expand(self.x.obj)

        for (coef, base) in zip(coefs, bases):
            self.surface[str(base)] = coef

        self.grad = MV()
        self.grad.is_grad = True
        self.grad.blade_rep = True
        self.grad.igrade = 1
        self.grad.rcpr_bases_MV = []
        for rbase in self.rbasis:
            self.grad.rcpr_bases_MV.append(rbase / self.E_sq)
        self.grad.rcpr_bases_MV = tuple(self.grad.rcpr_bases_MV)
        self.grad.coords = self.coords
        self.grad.norm = self.E_sq
        self.grad.connection = {}

        if debug:
            oprint('x', self.x,
                   'coords', self.coords,
                   'basis vectors', self.basis,
                   'index', self.index,
                   'basis blades', self.MFbasis,
                   'E', self.E,
                   'E**2', self.E_sq,
                   '*basis', duals,
                   'rbasis', self.rbasis,
                   'basis derivatives', self.dbasis,
                   'surface', self.surface,
                   'basis strings', self.basis_str,
                   'embedding basis', self.embedded_basis,
                   'metric tensor', self.gij)

    def Basis(self):
        return tuple(self.basis)

    def Grad(self, F):  # Intrisic Derivative
        dF = 0
        for (rbase, coord) in zip(self.rbasis, self.coords):
            dF += rbase * F.diff(coord)
        dF = dF.simplify()
        dF = dF / self.E_sq
        return dF

    def D(self, F):  # Covariant Derivative
        dF = self.Grad(F)
        return self.Proj(dF)

    def S(self, a):  # Shape Tensor

        return

    def Proj(self, F):
        PF = (F < self.E) * self.E
        PF = PF.simplify()
        PF = PF.trigsimp(deep=True, recursive=True)
        return (PF / self.E_sq).simplify()

    def Reject(self, F):
        return (F - self.Proj(F)).simplify()

    def DD(self, v, f, opstr=False):
        mf_comp = []
        for e in self.rbasis:
            mf_comp.append((v | e).scalar() / self.E_sq)
        result = MV()
        op = ''
        for (coord, comp) in zip(self.coords, mf_comp):
            result += comp * (f.diff(coord))
            if opstr:
                op += '(' + str(comp) + ')D{' + str(coord) + '}+'
        if opstr:
            return str(result), op[:-1]
        return result

    def Plot2DSurface(self, u_range, v_range, surf=True, grid=True, tan=1.0, scalar_field=None, skip=[1, 1], fct_def=None):

        plot_template = \
"""
from numpy import mgrid,shape,swapaxes,zeros,log,exp,sin,cos,tan
$fct_def$
eps = 1.0e-6
u_r = $u_range$
v_r = $v_range$
$coords$ = mgrid[u_r[0]:u_r[1]+eps:(u_r[1]-u_r[0])/float(u_r[2]-1),\\
                 v_r[0]:v_r[1]+eps:(v_r[1]-v_r[0])/float(v_r[2]-1)]
X = $surface$
scal_tan = $tan$
x = X['ex']
y = X['ey']
z = X['ez']
du = $basis_str[0]$
dv = $basis_str[1]$
Zero = zeros(shape(x))
if scal_tan > 0.0:
    du_x = Zero+du['ex']
    du_y = Zero+du['ey']
    du_z = Zero+du['ez']
    dv_x = Zero+dv['ex']
    dv_y = Zero+dv['ey']
    dv_z = Zero+dv['ez']

f = $scalar_field$
n = $n$
skip = $skip$
su = skip[0]
sv = skip[1]
if f[0] != None:
    dn_x = f[0]*n[0]
    dn_y = f[0]*n[1]
    dn_z = f[0]*n[2]

from mayavi.mlab import plot3d,quiver3d,mesh,figure
figure(bgcolor=(1.0,1.0,1.0))
if $surf$:
    mesh(x,y,z,colormap="gist_earth")
if $grid$:
    for i in range(shape(u)[0]):
        plot3d(x[i,],y[i,],z[i,],line_width=1.0,color=(0.0,0.0,0.0),tube_radius=None)

    xr = swapaxes(x,0,1)
    yr = swapaxes(y,0,1)
    zr = swapaxes(z,0,1)

    for i in range(shape(u)[1]):
        plot3d(xr[i,],yr[i,],zr[i,],line_width=1.0,color=(0.0,0.0,0.0),tube_radius=None)
if scal_tan > 0.0:
    quiver3d(x[::su,::sv],y[::su,::sv],z[::su,::sv],\\
             du_x[::su,::sv],du_y[::su,::sv],du_z[::su,::sv],scale_factor=scal_tan,\\
             line_width=1.0,color=(0.0,0.0,0.0),scale_mode='vector',mode='arrow',resolution=16)
    quiver3d(x[::su,::sv],y[::su,::sv],z[::su,::sv],\\
             dv_x[::su,::sv],dv_y[::su,::sv],dv_z[::su,::sv],scale_factor=scal_tan,\\
             line_width=1.0,color=(0.0,0.0,0.0),scale_mode='vector',mode='arrow',resolution=16)
if f[0] != None:
    quiver3d(x[::su,::sv],y[::su,::sv],z[::su,::sv],\\
             dn_x[::su,::sv],dn_y[::su,::sv],dn_z[::su,::sv],\\
             line_width=1.0,color=(0.0,0.0,0.0),scale_mode='none',mode='cone',\\
             resolution=16,opacity=0.5)

"""
        if len(self.coords) != 2:
            return

        self.skip = skip
        self.surf = surf
        self.grid = grid
        self.tan = tan
        if fct_def is None:
            self.fct_def = ' '
        else:
            self.fct_def = fct_to_str(fct_def)

        self.u_range = u_range
        self.v_range = v_range
        self.scalar_field = [scalar_field]

        print(self.I, '\n', self.basis[0], '\n', self.basis[1])

        self.normal = -self.I * (self.basis[0] ^ self.basis[1])
        self.n = VectorComponents(self.normal, ['ex', 'ey', 'ez'])

        msurf = open('manifold_surf.py', 'w')
        plot_template = FillTemplate(self, plot_template)
        msurf.write(plot_template)
        msurf.close()
        mayavi2 = find_executable('mayavi2')
        if mayavi2 is None:
            return
        system(mayavi2 + ' manifold_surf.py &')
        return
