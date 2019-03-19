# %load Contour Integration.py
from sympy import *
import matplotlib.pyplot as plt
from numpy import degrees

class ComplexIntegral:

    def __init__(self, eq, shape="semi-circle", params=dict(z0=0+0*I, radius=oo, theta1=0, theta2=pi)):

        self.eq = eq
        self.shape = shape

        if self.shape == "semi-circle" :
            self.solveit = self.solveit_semi_cir()
            self.z0 = params["z0"]
            self.radius = params["radius"]
            self.theta1 = params["theta1"]
            self.theta2 = params["theta2"]
            
        else:
            pass            #Enter in anticlockwise direction
            self.solveit = solveit_oth()
            self.points = points

    def contour_cir(self):
        self.radius=1
        plot(im(self.z0) + sqrt(self.radius**2 - (x - re(self.z0))**2),(x, re(self.z0) +  
                     self.radius * cos(self.theta2), im(self.z0) + self.radius * cos(self.theta1) ),
                     xlabel="Real", ylabel="Imaginary"  )
    

    def solveit_semi_cir(self):

        eqn_residues = self.eq.as_numer_denom()
        residues = list(roots(eqn_residues[1]).keys()) + list(roots((eqn_residues[0])**-1).keys())
        sum_residues = 0
        x = list(self.eq.atoms(Symbol))[0]
        for i in residues:   
            if im(i) > 0:    
                sum_residues += residue(self.eq, x, i)
        return simplify(2*pi*I*sum_residues)  
