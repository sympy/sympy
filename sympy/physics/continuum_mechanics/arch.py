"""
This module can be used to solve probelsm related to 2D parabolic arches
"""
from sympy.core.sympify import sympify
from sympy.core.symbol import Symbol,symbols
from sympy import diff, sqrt, cos , sin, rad
from sympy.core.relational import Eq
from sympy.solvers.solvers import solve

class Arch:
    """
    This class is used to solve problems related to a three hinged arch(determinate) structure.\n
    An arch is a curved vertical structure spanning an open space underneath it.\n
    Arches can be used to reduce the bending moments in long-span structures.\n

    Arches are used in structural engineering(over windows, door and even bridges)\n
    because they can support a very large mass placed on top of them.

    Example
    ========
    >>> from sympy.physics.continuum_mechanics.arch import Arch
    >>> a = Arch((0,0),(10,0),crown_x=5,crown_y=5)
    >>> a.get_parabola_eqn
    5 - (x - 5)**2/5

    >>> from sympy.physics.continuum_mechanics.arch import Arch
    >>> a = Arch((0,0),(10,1),crown_x=6)
    >>> a.get_parabola_eqn
    9/5 - (x - 6)**2/20
    """
    def __init__(self,left_support,right_support,**kwargs):
        self._shape_eqn = None
        self._left_support  = (sympify(left_support[0]),sympify(left_support[1]))
        self._right_support  = (sympify(right_support[0]),sympify(right_support[1]))
        self._crown_x = None
        self._crown_y = None
        if 'crown_x' in kwargs:
            self._crown_x = sympify(kwargs['crown_x'])
        if 'crown_y' in kwargs:
            self._crown_y = sympify(kwargs['crown_y'])
        self._shape_eqn = self.get_parabola_eqn
        self._conc_loads = {}
        self._distributed_loads = {}
        self._loads = {'concentrated': self._conc_loads, 'distributed':self._distributed_loads}
        self._supports = {'left':'hinged', 'right':'hinge'}
        self._rope = None
        self._reaction_force = {Symbol('R_A_x'):0, Symbol('R_A_y'):0, Symbol('R_B_x'):0, Symbol('R_B_y'):0}
        self._points_disc = []
        # self._crown = (sympify(crown[0]),sympify(crown[1]))

    @property
    def get_parabola_eqn(self):
        "returns the equation of the shape of arch developed"
        if self._shape_eqn:
            return self._shape_eqn

        x,y,c = symbols('x y c')
        a = Symbol('a',positive=False)
        if self._crown_x and self._crown_y:
            x0  = self._crown_x
            y0 = self._crown_y
            parabola_eqn = a*(x-x0)**2 + y0 - y
            eq1 = parabola_eqn.subs({x:self._left_support[0], y:self._left_support[1]})
            solution = solve((eq1),(a))
            parabola_eqn = solution[0]*(x-x0)**2 + y0
            if(parabola_eqn.subs({x:self._right_support[0]}) != self._right_support[1]):
                raise ValueError("provided coordinates of crown and supports are not consistent with parabolic arch")

        elif self._crown_x:
            x0  = self._crown_x
            parabola_eqn = a*(x-x0)**2 + c - y
            eq1 = parabola_eqn.subs({x:self._left_support[0], y:self._left_support[1]})
            eq2 = parabola_eqn.subs({x:self._right_support[0], y:self._right_support[1]})
            solution = solve((eq1,eq2),(a,c))
            if len(solution) <2 or solution[a] == 0:
                raise ValueError("parabolic arch cannot be constructed with the provided coordinates, try providing crown_y")
            parabola_eqn = solution[a]*(x-x0)**2+ solution[c]
            self._crown_y = solution[c]

        else:
            raise KeyError("please provide crown_x to contruct arch")

        return parabola_eqn

    @property
    def get_loads(self):
        """
        return the position of the applied load and angle (for concentrated loads)
        """
        return self._loads

    @property
    def supports(self):
        """
        Returns the type of support
        """
        return self._supports

    @property
    def left_support(self):
        """
        Returns the position of the left support.
        """
        return self._left_support

    @property
    def right_support(self):
        """
        Returns the position of the right support.
        """
        return self._right_support

    @property
    def reaction_force(self):
        """
        return the reaction forces generated
        """
        return self._reaction_force

    def apply_load(self,order,label,x1,mag,x2=None,angle=None):
        """
        This method adds load to the Arch.

        Parameters
        ==========

            order : Integer
                Order of the applied load.

                    - For point/concentrated loads, order = -1
                    - For distributed load, order = 0

            label : String or Symbol
                The label of the load

            x1 : Sympifyable

                    - For concentrated/point loads, x1 is the x coordinate
                    - For distributed loads, x1 is the starting position of distributed load

            mag : Sympifyable
                Magnitude of the appliead load. Must be positive

            x2 : Sympifyable
                Required for distributed loads

                    - For concentrated/point load , x2 is None(may not be given)
                    - For distributed loads, x2 is the end position of distributed load

            angle: Sympifyable
                The angle in degrees, the load vector makes with the horizontal
                in the counter-clockwise direction.

        Examples
        ========
        For applying distributed load

        >>> from sympy.physics.continuum_mechanics.arch import Arch
        >>> a = Arch((0,0),(10,0),crown_x=5,crown_y=5)
        >>> a.apply_load(0,'A',x1=3,x2=5,mag=10)

        For applying point/concentrated_loads

        >>> from sympy.physics.continuum_mechanics.arch import Arch
        >>> a = Arch((0,0),(10,0),crown_x=5,crown_y=5)
        >>> a.apply_load(-1,'B',x1=2,mag=15,angle=45)

        """
        if label in self._loads:
            raise ValueError("load with the given label already exists")

        if order == 0:
            if not x2 or x2<x1:
                raise KeyError("provide x2 greater than x1")

            if x1>self._right_support[0] or x2<self._left_support[0]:
                raise ValueError(f"loads must be applied between {self._left_support[0]} and {self._right_support[0]}")
            self._distributed_loads[label] = {'start':x1, 'end':x2, 'f_y': -mag}
            self._points_disc.append(x1)
            self._points_disc.append(x2)

        if order == -1:
            if not angle:
                raise TypeError("please provide direction of force")
            y = self._shape_eqn.subs({'x':x1})
            if x1>self._right_support[0] or x1<self._left_support[0]:
                raise ValueError(f"loads must be applied between x = {self._left_support[0]} and x = {self._right_support[0]}")
            self._conc_loads[label] = {'x':x1, 'y':y, 'f_x':mag*cos(rad(angle)), 'f_y': mag*sin(rad(angle)), 'magnitude':mag, 'angle':angle}
            self._points_disc.append(x1)

    def remove_load(self,label):
        """
        This methods removes the load applied to the arch

        Parameters
        ==========

        label : String or Symbol
            The label of the applied load
        """
        if label in self._distributed_loads :
            val = self._distributed_loads.pop(label)
            print(f"removed load {label}: {val}")
        elif label in self._conc_loads :
            val = self._conc_loads.pop(label)
            print(f"removed load {label}: {val}")
        else :
            raise ValueError("label not found")

    def add_support(self,left_support,right_support):
        """
        Add the type for support at each end.
        Can use roller or hinge support at each end.
        """
        support_types = ['roller','hinge']
        if left_support not in support_types or right_support not in support_types:
            raise ValueError("supports must only be roller or hinge")
        self._supports['left'] = left_support
        self._supports['right'] = right_support

    def add_rope(self,y):
        if y>=self._crown_y or y<=min(self._left_support[1],  self._right_support[1]):
            raise ValueError(f"position of rope must be between y={min(self._left_support[1],  self._right_support[1])} and y={self._crown_y}")
        x = Symbol('x')
        a = diff(self._shape_eqn,x).subs(x,self._crown_x+1)/2
        x_diff = sqrt((y - self._crown_y)/a)
        x1 = self._crown_x + x_diff
        x2 = self._crown_x - x_diff
        self._rope = (x1,x2,y)

    def solve(self):
        """
        This method solves for the reaction forces generated at the supports,\n
        bending moment and shear force generated in the arch and tension produced in the rope if used.
        """
        # discontinuity_points = sorted(self._points_disc)

        # for reaction forces
        net_x = 0
        net_y = 0
        moment_A = 0
        moment_hinge_right = 0
        for label in self._conc_loads:
            net_x += self._conc_loads[label]['f_x']
            net_y += self._conc_loads[label]['f_y']
            moment_A += self._conc_loads[label]['f_y']*(self._conc_loads[label]['x']-self._left_support[0]) - \
                        self._conc_loads[label]['f_x']*(self._conc_loads[label]['y']-self._left_support[1])

            if self._conc_loads[label]['x']> self._crown_x:
                moment_hinge_right += self._conc_loads[label]['f_y']*(self._conc_loads[label]['x']-self._crown_x) - \
                        self._conc_loads[label]['f_x']*(self._conc_loads[label]['y']-self._crown_y)

        for label in self._distributed_loads:
            start = self._distributed_loads[label]['start']
            end = self._distributed_loads[label]['end']
            tot_force = self._distributed_loads[label]['f_y']*( end - start)
            net_y += tot_force
            moment_A += tot_force*((end+start)/2 - self._left_support[0])

            if self._distributed_loads[label]['end']>self._crown_x:
                st = max(start,self._crown_x)
                force_right = self._distributed_loads[label]['f_y']*(end-st)
                moment_hinge_right += force_right*((end+st)/2 - self._crown_x)

        R_A_x, R_A_y, R_B_x, R_B_y, T = symbols('R_A_x R_A_y R_B_x R_B_y T')

        if self._supports['left'] == 'roller' and self._supports['right'] == 'roller':
            if self._rope[2]>max(self._left_support[1],self._right_support[1]):
                if net_x!=0:
                    raise ValueError("net force in x direction not possible under the specified conditions")
                else:
                    eq1 = Eq(R_A_x ,0)
                    eq2 = Eq(R_B_x, 0)
                    eq3 = Eq(R_A_y + R_B_y + net_y,0)
                    eq4 = Eq(R_B_y*(self._right_support[0]-self._left_support[0])-R_B_x*(self._right_support[1]-self._left_support[1])+moment_A,0)
                    eq5 = Eq(moment_hinge_right + R_B_y*(self._right_support[0]-self._crown_x) + T*(self._rope[2]-self._crown_y),0)
                    solution = solve((eq1,eq2,eq3,eq4,eq5),(R_A_x,R_A_y,R_B_x,R_B_y,T))
            elif self._rope[2]>self._left_support[1]:
                if net_x>0:
                    raise ValueError("net positive force in x direction is not possible under the specified conditions")
                else:
                    eq1 = Eq(R_A_x ,0)
                    eq2 = Eq(R_B_x, 0)
                    eq3 = Eq(R_A_y + R_B_y + net_y,0)
                    eq4 = Eq(R_B_y*(self._right_support[0]-self._left_support[0])-T*(self._rope[2]-self._left_support[1])+moment_A,0)
                    eq5 = Eq(T+net_x,0)
                    solution = solve((eq1,eq2,eq3,eq4,eq5),(R_A_x,R_A_y,R_B_x,R_B_y,T))
            elif self._rope[2]>self._right_support[1]:
                if net_x>0:
                    raise ValueError("net negative force in x direction is not possible under the specified conditions")
                else:
                    eq1 = Eq(R_A_x ,0)
                    eq2 = Eq(R_B_x, 0)
                    eq3 = Eq(R_A_y + R_B_y + net_y,0)
                    eq4 = Eq(R_B_y*(self._right_support[0]-self._left_support[0])+T*(self._rope[2]-self._left_support[1])+moment_A,0)
                    eq5 = Eq(T-net_x,0)
                    solution = solve((eq1,eq2,eq3,eq4,eq5),(R_A_x,R_A_y,R_B_x,R_B_y,T))

        elif self._supports['left'] == 'roller':
            if self._rope[2]>max(self._left_support[1], self._right_support[1]):
                eq1 = Eq(R_A_x ,0)
                eq2 = Eq(R_B_x+net_x,0)
                eq3 = Eq(R_A_y + R_B_y + net_y,0)
                eq4 = Eq(R_B_y*(self._right_support[0]-self._left_support[0])-R_B_x*(self._right_support[1]-self._left_support[1])+moment_A,0)
                eq5 = Eq(R_A_y*(self._left_support[0]-self._crown_x)-T*(self._rope[2]-self._crown_y),0)
                solution = solve((eq1,eq2,eq3,eq4,eq5),(R_A_x,R_A_y,R_B_x,R_B_y,T))

            elif self._rope[2]>self._left_support[0]:
                if net_x>0:
                    eq1 = Eq(R_A_x,0)
                    eq2 = Eq(T,0)
                    eq3 = Eq(R_A_y + R_B_y + net_y, 0)
                    eq4 = Eq(R_B_x+net_x,0)
                    eq5 = Eq(R_B_y*(self._right_support[0]-self._left_support[0])-R_B_x*(self._right_support[1]-self._left_support[1])+moment_A,0)
                    solution = solve((eq1,eq2,eq3,eq4,eq5),(R_A_x,R_B_x,R_B_y,R_A_x,T))
                else:
                    eq1 = Eq(R_A_x ,0)
                    eq2 = Eq(R_B_x+ T +net_x,0)
                    eq3 = Eq(R_A_y + R_B_y + net_y,0)
                    eq4 = Eq(R_B_y*(self._right_support[0]-self._left_support[0])-R_B_x*(self._right_support[1]-self._left_support[1])+moment_A,0)
                    eq5 = Eq(R_A_y*(self._left_support[0]-self._crown_x)-T*(self._rope[2]-self._crown_y),0)
                    solution = solve((eq1,eq2,eq3,eq4,eq5),(R_A_x,R_A_y,R_B_x,R_B_y,T))

            elif self._rope[2]>self._right_support[0]:
                pass

        elif self._supports['right'] == 'roller':
            pass
        else:
            eq1 = Eq(R_A_x + R_B_x + net_x,0)
            eq2 = Eq(R_A_y + R_B_y + net_y,0)
            eq3 = Eq(R_B_y*(self._right_support[0]-self._left_support[0])-R_B_x*(self._right_support[1]-self._left_support[1])+moment_A,0)
            eq4 = Eq(moment_hinge_right + R_B_y*(self._right_support[0]-self._crown_x) + R_B_x*(self._right_support[1]-self._crown_y),0)

            solution = solve((eq1,eq2,eq3,eq4),(R_A_x,R_A_y,R_B_x,R_B_y))

        print(solution)
