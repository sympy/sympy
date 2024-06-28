from sympy import sympify, symbols, solve, Symbol

class Arches:
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
        self.get_parabola_eqn
        self._loads = {}
        self._conc_loads = {}
        self._distributed_loads = {}
        self._supports = {'left':None,'right':None}
        self._rope = None
        self._reaction_force = {Symbol('R_A_x'):0, Symbol('R_A_y'):0, Symbol('R_B_x'):0, Symbol('R_B_y'):0}
        # self._crown = (sympify(crown[0]),sympify(crown[1]))

    @property
    def get_parabola_eqn(self):
        if self._shape_eqn:
            return self._shape_eqn

        x,y,a,c = symbols('x y a c')

        if self._crown_x and self._crown_y:
            x0  = self._crown_x
            y0 = self._crown_y
            parabola_eqn = a*(x-x0)**2 + y0 - y
            eq1 = parabola_eqn.subs({x:self._left_support[0], y:self._left_support[1]})
            solution = solve((eq1),(a))
            parabola_eqn = solution[0]*(x-x0)**2 + y0
            if(parabola_eqn.subs({x:self._right_support[0]}) != self._right_support[1]) or solution[0]>=0:
                raise ValueError("provided coordinates of crown and supports are not consistent with parabolic arch")
            self._shape_eqn = parabola_eqn

        elif self._crown_x:
            x0  = self._crown_x
            parabola_eqn = a*(x-x0)**2 + c - y
            eq1 = parabola_eqn.subs({x:self._left_support[0], y:self._left_support[1]})
            eq2 = parabola_eqn.subs({x:self._right_support[0], y:self._right_support[1]})
            solution = solve((eq1,eq2),(a,c))
            if len(solution) <2:
                raise ValueError("parabolic arch cannot be constructed with the provided coordinates")
            parabola_eqn = solution[a]*(x-x0)**2+ solution[c]
            self._crown_y = solution[c]
            self._shape_eqn = parabola_eqn

        else:
            raise KeyError("please provide crown_x to contruct arch")

    @property
    def get_loads(self):
        loads = {'distributed':self._distributed_loads, 'concentrated':self._conc_loads}
        return loads, self._loads

    @property
    def get_reaction_force(self):
        return self._reaction_force

    def apply_load(self,order,label,x1,mag,x2=None,angle=None):
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
            self._conc_loads[label] = (x1,y)

    def remove_loads(self,order,label):
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
        support_types = ['roller','hinge']
        if left_support not in support_types or right_support not in support_types:
            raise ValueError("supports must only be roller or hinged")
        self._supports['left'] = left_support
        self._right_support['right'] = right_support

    def add_rope(self,start,end):
        if start<self._left_support[0] or end >self._right_support[0]:
            raise ValueError(f"start and end point of rope must be between {self._left_support[0]} and {self._right_support[0]}")

        x = symbols('x')
        y0 = self._shape_eqn.subs({'x': start})
        y1 = self._shape_eqn.subs({'x': end})
        a = (y1-y0)/(end-start)
        y = a*(x-start) + y0
        self._rope = y
