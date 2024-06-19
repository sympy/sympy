from sympy import sympify, symbols, solve

class Arches:
    def __init__(self,left_support,right_support,**kwargs):
        self._shape_eqn = None
        self._left_support  = (sympify(left_support[0]),sympify(left_support[1]))
        self._right_support  = (sympify(right_support[0]),sympify(right_support[1]))
        self._crown_x = sympify(kwargs['crown_x'])
        self._crown_y = sympify(kwargs['crown_y'])
        self.get_parabola_eqn
        # self._crown = (sympify(crown[0]),sympify(crown[1]))

    @property
    def get_parabola_eqn(self):
        if self._shape_eqn:
            return self._shape_eqn

        x,y,a,c = symbols('x y a c')
        x0 = self._crown_x
        y0 = self._crown_y

        if self._crown_x and self._crown_y:
            parabola_eqn = a*(x-x0)**2 + y0 - y
            eq1 = parabola_eqn.subs({x:self._left_support[0], y:self._left_support[1]})
            solution = solve((eq1),(a))
            parabola_eqn = solution[0]*(x-x0)**2 + y0
            if(parabola_eqn.subs({x:self._right_support[0]}) != self._right_support[1]) or solution[0]>=0:
                raise ValueError("provided coordinates of crown and supports are not consistent with parabolic arch")
            self._shape_eqn = parabola_eqn
