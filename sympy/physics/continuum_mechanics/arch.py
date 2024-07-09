"""
This module can be used to solve probelsm related to 2D parabolic arches
"""
from sympy import sympify, symbols, solve, Symbol

class Arch:
    """
    An arch is a curved vertical structure spanning an open space underneath it.
    Arches can be used to reduce the bending moments in long-span structures.

    Arches are used in structural engineering(over windows, door and even bridges)
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
        self._loads = {}
        self._conc_loads = {}
        self._distributed_loads = {}
        self._supports = {'left':None, 'right':None}
        self._rope = None
        self._reaction_force = {Symbol('R_A_x'):0, Symbol('R_A_y'):0, Symbol('R_B_x'):0, Symbol('R_B_y'):0}
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
        loads = {'distributed':self._distributed_loads, 'concentrated':self._conc_loads}
        return loads

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

        """
        if label in self._loads:
            raise ValueError("load with the given label already exists")

        self._loads[label] = mag
        if order == 0:
            if not x2 or x2<x1:
                raise KeyError("provide x2 greater than x1")

            if x1>self._right_support[0] or x2<self._left_support[0]:
                raise ValueError(f"loads must be applied between {self._left_support[0]} and {self._right_support[0]}")
            self._distributed_loads[label] = (x1,x2)
        if order == 1:
            if not angle:
                raise TypeError("please provide direction of force")

            y = self._shape_eqn.subs({'x':x1})
            if x1>self._right_support[0] or x1<self._left_support[0]:
                raise ValueError(f"loads must be applied between x = {self._left_support[0]} and x = {self._right_support[0]}")
            self._conc_loads[label] = (x1,y,angle)

    def remove_load(self,order,label):
        """
        This methods removes the load applied to the arch

        Parameters
        ==========

        order : Integer
            The order of the appplied load.

                - For point loads, order = -1
                - For distributed load, order = 0

        label : String or Symbol
            The label of the applied load
        """
        if label in self._loads:
            mag = self._loads[label]
            self._loads.pop(label)
        else:
            raise KeyError("no such load applied")

        if order==0 and label in self._distributed_loads:
            self._distributed_loads.pop(label)
        elif order==1 and label in self._conc_loads:
            self._conc_loads.pop(label)
        else:
            self._loads[label] = mag
            raise KeyError("no such load in the provided load type or load type does not exist")

    def add_support(self,left_support,right_support):
        """
        Add the type for support at each end.
        Can use roller or hinge support at each end.
        """
        support_types = ['roller','hinge']
        if left_support not in support_types or right_support not in support_types:
            raise ValueError("supports must only be roller or hinged")
        self._supports['left'] = left_support
        self._supports['right'] = right_support

    def add_rope(self,start,end):
        if start<self._left_support[0] or end >self._right_support[0]:
            raise ValueError(f"start and end point of rope must be between {self._left_support[0]} and {self._right_support[0]}")

        x = symbols('x')
        y0 = self._shape_eqn.subs({'x': start})
        y1 = self._shape_eqn.subs({'x': end})
        a = (y1-y0)/(end-start)
        y = a*(x-start) + y0
        self._rope = y
