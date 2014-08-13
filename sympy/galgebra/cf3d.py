#cf3d.py

from sympy.galgebra.ga import Ga
from sympy.galgebra.mv import Mv
from sympy import symbols, S, Symbol

class Cf3d(object):

    ga_flg = False

    @staticmethod
    def make_cf3d():
        if Cf3d.ga_flg:
            return
        Cf3d.ga_flg = True
        Cf3d.basis_sups = ['__x','__y','__z']
        Cf3d.basis = 'e_x e_y e_z e ebar'
        (Cf3d.cf3d,Cf3d.ex,Cf3d.ey,Cf3d.ez,Cf3d.e,Cf3d.ebar) = Ga.build(Cf3d.basis,g=[1,1,1,1,-1])
        Cf3d.xyz = (Cf3d.ex,Cf3d.ey,Cf3d.ez)
        Cf3d.n = Cf3d.e + Cf3d.ebar
        Cf3d.nbar = Cf3d.e - Cf3d.ebar

        return

    @staticmethod
    def basis():
        Cf3d.make_cf3d()
        return (Cf3d.ex,Cf3d.ey,Cf3d.ez)

    def __init__(self,vec,mode='point'):

        Cf3d.make_cf3d()

        if isinstance(vec,Mv) and vec.is_vector() and mode == 'cf3d':
            self.mode = 'cf3d'
            self.x = vec

        if isinstance(vec,Mv) and vec.is_vector() and mode == 'point':  # define conformal 3 point from 3 vector
            x_sq = vec | vec
            self.x = x_sq * Cf3d.n + S(2) * vec - Cf3d.nbar
            self.mode = 'point'

        if isinstance(vec,Mv) and vec.is_vector() and mode == 'trans':
            self.x = S(1) + (Cf3d.n * vec) / S(2)
            self.mode = 'trans'

        if isinstance(vec,Mv) and not vec.is_vector():
            self.x = vec
            self.mode = 'cf3d'

    def X3d(self):
        return ((self.x | Cf3d.ex) * Cf3d.ex + (self.x | Cf3d.ey) * Cf3d.ey + (self.x | Cf3d.ez) * Cf3d.ez) / S(2)

    def __str__(self):

        pt = self.X3d()
        return str(pt)

    def __mul__(self, X):
        if isinstance(X, Cf3d):
            return Cf3d(self.x * X.x, mode='cf3d')
        else:
            return Cf3d(self.x * X, mode='cf3d')

    def __rmul__(self,X):
        if isinstance(X,Cf3d):
            return Cf3d(X.x * self.x, mode='cf3d')
        else:
            return  Cf3d(X * self.x, mode='cf3d')

    def rev(self):
        return Cf3d(self.x.rev(), mode=self.mode)

    def __add__(self, X):
        if isinstance(self, Cf3d) and isinstance(X, Cf3d):
            x = self.X3d() + X.X3d()
            return Cf3d(x,mode=self.mode)

    def __sub__(self, X):
        if isinstance(self, Cf3d) and isinstance(X, Cf3d):
            x = self.X3d() - X.X3d()
            return Cf3d(x,mode=self.mode)

    @staticmethod
    def make_vector(base_str):
        coefs = []
        vec = S(0)
        for (coord,base) in zip(Cf3d.basis_sups,Cf3d.xyz):
            vec += Symbol(base_str + coord, real=True) * base

        return (vec,Cf3d(vec,mode='point'))

if __name__ == "__main__":

    a = symbols('a__x a__y a__z', real=True)

    basis = (ex,ey,ez) = Cf3d.basis()

    """
    a1 = -8649.51*ex+4688.51*ey+600.0*ez
    b1 = 4557.58*ex+679.734*ey+302.5*ez
    a2 = -8625.71*ex+4720.65*ey+600.0*ez
    b2 = 4545.3*ex+641.665*ey+302.5*ez

    A1 = Cf3d(a1)
    B1 = Cf3d(b1)
    A2 = Cf3d(a2)
    B2 = Cf3d(b2)
    """

    (a1,A1) = Cf3d.make_vector('a1')
    (b1,B1) = Cf3d.make_vector('b1')
    (a2,A2) = Cf3d.make_vector('a2')
    (b2,B2) = Cf3d.make_vector('b2')

    print 'a1 =', A1
    print 'b1 =', B1
    print 'a2 =', A2
    print 'b2 =', B2

    b1_to_a1 = Cf3d(a1-b1,mode='trans')

    b1_p = b1_to_a1*B1*b1_to_a1.rev()

    print 'b1 translated to a_1 =',b1_p
    print 'b1_p - b1 =', b1_p - A1

    a1a2 = a1 - a2

    print a1a2

