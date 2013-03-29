#Manifold.py

import ga_dir
from os import system
import copy
from itertools import izip,islice,combinations,imap,product,ifilter
from sympy import trigsimp,simplify

if ga_dir.GA =='GA':
    from ga import MV
    from ga_sympy import linear_expand
    from ga_debug import oprint
else:
    from sympy.ga.ga import MV
    from sympy.ga.ga_sympy import linear_expand
    from sympy.ga.ga_debug import oprint

def fct_to_str(fct_names):
    import sys
    current_file = open(sys.argv[0],'r')
    file_str = current_file.read()
    current_file.close()

    if isinstance(fct_names,str):
        return(fct_names)

    fcts_str = ''

    for fct_name in fct_names:
        start_def   = file_str.find('\ndef '+fct_name)
        end_def     = file_str.find('\ndef ',start_def+5)
        start_class = file_str.find('\nclass ',start_def+5)
        end_def = min(end_def,start_class)
        fcts_str += file_str[start_def:end_def]
    return(fcts_str)

def VectorComponents(X,basis):
    (coefs,bases) = linear_expand(X.obj)
    cdict = {}
    for (coef,base) in zip(coefs,bases):
        cdict[str(base)] = coef
    comp = []
    for base in basis:
        if base in cdict:
            comp.append(cdict[base])
        else:
            comp.append(0)
    return(comp)

def FillTemplate(self,template):
    Nd = 0
    var = []
    id_old = 0
    while True:
        id_new = template.find('$',id_old+1)
        if id_new == -1:
            break
        Nd += 1
        if Nd%2 == 0:
            var.append(template[id_old+1:id_new])
        id_old = id_new

    var.sort(reverse=True)

    for v in var:
        template = template.replace('$'+v+'$',str(eval('self.'+v)))

    return(template)

class Manifold:

    def __init__(self,x,coords,debug=False,I=None):
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
            (coefs,bases) = linear_expand(tv.obj)
            tc = {}
            for (coef,base) in zip(coefs,bases):
                str_base = str(base)
                tc[str_base] = coef
                if str_base not in self.embedded_basis:
                    self.embedded_basis.append(str_base)
            self.basis_str.append(tc)

        for tv in self.basis_str:
            for base in self.embedded_basis:
                if base not in tv:
                    tv[base] = 0

        self.dim = len(self.basis)

        indexes = tuple(range(self.dim))
        self.index = [()]
        for i in indexes:
            self.index.append(tuple(combinations(indexes,i+1)))
        self.index = tuple(self.index)

        self.MFbasis = [[MV.ONE],self.basis]

        for igrade in self.index[2:]:
            grade = []
            for iblade in igrade:
                blade = MV(1,'scalar')
                for ibasis in iblade:
                    blade ^= self.basis[ibasis]
                blade = blade.trigsimp(deep=True,recursive=True)
                grade.append(blade)
            self.MFbasis.append(grade)
        self.E = self.MFbasis[-1][0]
        self.E_sq = trigsimp((self.E*self.E).scalar(),deep=True,recursive=True)

        duals = copy.copy(self.MFbasis[-2])

        duals.reverse()
        sgn = 1
        self.rbasis = []
        for dual in duals:
            recpv = (sgn*dual*self.E).trigsimp(deep=True,recursive=True)
            self.rbasis.append(recpv)
            sgn = -sgn

        self.dbasis = []

        for base in self.basis:
            dbase = []
            for coord in self.coords:
                d = base.diff(coord).trigsimp(deep=True,recursive=True)
                dbase.append(d)
            self.dbasis.append(dbase)

        self.surface = {}
        (coefs,bases) = linear_expand(self.x.obj)

        for (coef,base) in zip(coefs,bases):
            self.surface[str(base)] = coef

        self.grad = MV()
        self.grad.is_grad   = True
        self.grad.blade_rep = True
        self.grad.igrade    = 1
        self.grad.rcpr_bases_MV  = []
        for rbase in self.rbasis:
            self.grad.rcpr_bases_MV.append(rbase/self.E_sq)
        self.grad.rcpr_bases_MV =tuple(self.grad.rcpr_bases_MV)
        self.grad.coords = self.coords
        self.grad.norm  = self.E_sq
        self.grad.connection = {}

        if debug:
            oprint('x',self.x,\
                   'coords',self.coords,\
                   'basis vectors',self.basis,\
                   'index',self.index,\
                   'basis blades',self.MFbasis,\
                   'E',self.E,\
                   'E**2',self.E_sq,\
                   '*basis',duals,\
                   'rbasis',self.rbasis,\
                   'basis derivatives',self.dbasis,\
                   'surface',self.surface,\
                   'basis strings',self.basis_str,\
                   'embedding basis',self.embedded_basis)


    def Basis(self):
        return(tuple(self.basis))

    def Grad(self,F):
        dF = 0
        for (rbase,coord) in zip(self.rbasis,self.coords):
            dF += rbase*F.diff(coord)
        dF = dF.simplify()
        dF = dF/self.E_sq
        return(dF)

    def Proj(self,F):
        PF = (F<self.E)*self.E
        PF = PF.simplify()
        PF = PF.trigsimp(deep=True,recursive=True)
        return(PF/self.E_sq)

    def DD(self,v,f,opstr=False):
        mf_comp = []
        for e in self.rbasis:
            mf_comp.append((v|e).scalar()/self.E_sq)
        result = MV()
        op = ''
        for (coord,comp) in zip(self.coords,mf_comp):
            result += comp*(f.diff(coord))
            if opstr:
                op += '('+str(comp)+')D{'+str(coord)+'}+'
        if opstr:
            return(str(result),op[:-1])
        return(result)

    def Plot2DSurface(self,u_range,v_range,surf=True,grid=True,tan=1.0,scalar_field=None,skip=[1,1],fct_def=None):

        plot_template = \

        if len(self.coords) != 2:
            return

        self.skip = skip
        self.surf = surf
        self.grid = grid
        self.tan  = tan
        if fct_def == None:
            self.fct_def = ' '
        else:
            self.fct_def = fct_to_str(fct_def)

        self.u_range = u_range
        self.v_range = v_range
        self.scalar_field = [scalar_field]

        self.normal = -self.I*(self.basis[0]^self.basis[1])
        self.n = VectorComponents(self.normal,['ex','ey','ez'])

        msurf = open('manifold_surf.py','w')
        plot_template = FillTemplate(self,plot_template)
        msurf.write(plot_template)
        msurf.close()
        #system('geany manifold_surf.py &')
        system('mayavi2 manifold_surf.py &')
        return
