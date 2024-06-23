from sympy import sympify, symbols, solve

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

    def apply_load(self,order,label,coord1,coord2,mag,*args):
        if label in self._loads:
            raise ValueError("load with the given label already exists")

        self._loads[label] = mag
        if order == 0:
            x0 = min(coord1,coord2)
            x1 = max(coord1,coord2)
            if x0>self._right_support[0] or x1<self._left_support[0]:
                raise ValueError(f"loads must be applied between {self._left_support[0]} and {self._right_support[0]}")
            self._distributed_loads[label] = (x0,x1)
        if order == 1:
            y0 = min(self._right_support[1],self._left_support[1])
            y1 = max(self._right_support[1],self._left_support[1])

            if coord1>self._right_support[0] or coord1<self._left_support[0]:
                raise ValueError(f"loads must be applied between x = {self._left_support[0]} and x = {self._right_support[0]}")
            if coord2>y1 or coord2<y0:
                raise ValueError(f"loads must be applied between y = {y0} and y = {y1}")
            self._conc_loads[label] = (coord1,coord2,args[0])

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
