from sympy import *
import matplotlib.pyplot as plt

class ComplexIntegral:

    def __init__(self, eq, params=dict(shape="semi-circle", origin=(0,0), radius=oo, angle=pi)):

        self.eq = eq
        self.shape = params["shape"]
        if "circle" in self.shape:
            self.solveit = self.solveit_cir()
            self.origin = params["origin"]
            self.rad = params["radius"]
            self.angle = params["angle"]

        else:
            pass 
            #Enter in anticlockwise direction
            self.solveit = solveit_oth()
            self.points = points

    def contour(self):

        #Plot the contour
        pass   

    def residue_calc(self, rooot, mult, x ):
        
        return Rational(1/factorial(mult-1)) *limit(diff(self.eq*(x - rooot)**mult, x, mult - 1), x, complex(rooot))
    
    def solveit_cir(self):

        eqn_residues = solve(self.eq.as_numer_denom()[1]) 
        residues = roots(eqn_residues)
        sum_residues = 0
        x = list(self.eq.atoms())[0]
        for i, j in residues.items(): 
            sum_residues += self.residue_calc(i, j, x)
    
        return simplify(2*pi*I*sum_residues)  