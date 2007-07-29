from sympy import Symbol, Basic, Integer

class PlotInterval(object):
    """
    """
    def __init__(self, *args):
        if len(args) == 1:
            if isinstance(args[0], str):
                try:
                    args = eval(args[0])
                except:
                    s_eval_error = "Could not interpret string %s."
                    raise ValueError(s_eval_error % (args[0]))
            elif isinstance(args[0], (tuple, list)):
                args = args[0]
        if not isinstance(args, (tuple, list)) or len(args) not in [3,4]:
            f_error = "Interval must be in the format var, min, max, [steps]."
            raise ValueError(f_error)

        if len(args) == 4:
            v, v_min, v_max, v_steps = args
        else:
            v, v_min, v_max = args
            v_steps = 20

        if not isinstance(v, Symbol):
            raise ValueError("v must be a sympy Symbol.")
        self.v = v

        try:
            self.v_min = Basic.sympify(v_min)
        except:
            raise ValueError("v_min could not be interpreted as a number.")

        try:
            self.v_max = Basic.sympify(v_max)
        except:
            raise ValueError("v_max could not be interpreted as a number.")

        if isinstance(v_steps, int):
            self.v_steps = Integer(v_steps)
        elif isinstance(v_steps, Integer):
            self.v_steps = v_steps
        else:
            raise ValueError("v_steps must be an int or sympy Integer.")

        if self.v_steps <= Integer(0):
            raise ValueError("v_steps must be positive.")

        self.v_len = self.v_steps + 1

    @staticmethod
    def try_parse(*args):
        """
        Returns a PlotInterval if args can be interpreted
        as such, otherwise None.
        """
        if len(args) == 1 and isinstance(args[0], PlotInterval):
            return args[0]
        try:
            return PlotInterval(*args)
        except:
            return None

    def _str_base(self):
        return "%s, %s, %s, %i" % (str(self.v),self.v_min,self.v_max,self.v_steps)

    def __repr__(self):
        """
        A string representing the interval in class constructor form.
        """
        return "PlotInterval(%s)" % (self._str_base())

    def __str__(self):
        """
        A string representing the interval in list form.
        """
        return "[%s]" % (self._str_base())

    def vrange(self):
        """
        Yields v_steps+1 sympy numbers ranging from
        v_min to v_max.
        """
        d = (self.v_max - self.v_min) / self.v_steps
        for i in xrange(self.v_steps+1):
            a = self.v_min + (d * Integer(i))
            yield a

    def vrange2(self):
        """
        Yields v_steps pairs of sympy numbers ranging from
        (v_min, v_min + step) to (v_max - step, v_max).
        """
        d = (self.v_max - self.v_min) / self.v_steps
        a = self.v_min + (d * Integer(0))
        for i in xrange(self.v_steps):
            b = self.v_min + (d * Integer(i+1))
            yield a, b
            a = b

