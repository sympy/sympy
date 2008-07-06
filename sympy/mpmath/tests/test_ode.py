from sympy.mpmath.calculus import ODE_step_euler, ODE_step_rk4, odeint, arange

solvers = [ODE_step_euler, ODE_step_rk4]

def test_ode1():
    """
    Let's solve:

    x'' + w**2 * x = 0

    i.e. x1 = x, x2 = x1':

    x1' =  x2
    x2' = -x1
    """
    def derivs((x1, x2), t):
        return x2, -x1

    for solver in solvers:
        t = arange(0, 3.1415926, 0.005)
        sol = odeint(derivs, (0., 1.), t, solver)
        x1 = [a[0] for a in sol]
        x2 = [a[1] for a in sol]
        # the result is x1 = sin(t), x2 = cos(t)
        # let's just check the end points for t = pi
        assert abs(x1[-1]) < 1e-2
        assert abs(x2[-1] - (-1)) < 1e-2

def test_ode2():
    """
    Let's solve:

    x' - x = 0

    i.e. x = exp(x)

    """
    def derivs((x), t):
        return x

    for solver in solvers:
        t = arange(0, 1, 1e-3)
        sol = odeint(derivs, (1.,), t, solver)
        x = [a[0] for a in sol]
        # the result is x = exp(t)
        # let's just check the end point for t = 1, i.e. x = e
        assert abs(x[-1] - 2.718281828) < 1e-2
