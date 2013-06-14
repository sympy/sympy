from sympy import degree_list
from sympy import igcd
from sympy import symbols
from sympy import Add
from sympy import Integer


def diop_solve(eq):

    coeff, var, t = classify_diop(eq)

    if t == "linear":
        return diop_linear(var, coeff)


def classify_diop(eq):

    var = list(eq.free_symbols)
    var.sort()

    if max(degree_list(eq)) == 1:
        diop_type = "linear"
    else:
        raise NotImplementedError("Still not implemented")


    coeff = dict([reversed(t.as_independent(*var)) for t in eq.args])
    for v in coeff:
        if not isinstance(coeff[v], Integer):
            raise TypeError("Coefficients should be Integers")

    return coeff, var, diop_type,


def diop_linear(var, coeff):

    x = var[0]; y = var[1]
    a = coeff[x]; b = coeff[y]

    if len(var) == len(coeff):
        c = 0
    else:
        c = -coeff[Integer(1)]

    if len(var) == 2:
        sol_x, sol_y = base_solution(c, a, b, True)
        return {x: sol_x, y: sol_y}

    elif len(var) > 2:
        X = []; Y = []

        for v in var[2:]:
            sol_x, sol_y  = base_solution(-coeff[v], a, b)
            X.append(sol_x*v); Y.append(sol_y*v)

        sol_x, sol_y = base_solution(c, a, b, True)
        X.append(sol_x); Y.append(sol_y)

        l = []
        if None not in X and None not in Y:
            l.append((x, Add(*X))); l.append((y, Add(*Y)))
            for v in var[2:]:
                l.append((v, v))
        else:
            for v in var:
                l.append((v, None))

        return dict(l)


def base_solution(c, a, b, param=False):

    t = symbols("t", type = Integer)

    if c == 0:
        if param:
            return (b*t , -a*t)
        else:
            return (0, 0)
    else:
        x0, y0, d = extended_euclid(abs(a), abs(b))

        x0 = x0 * (a / abs(a))
        y0 = y0 * (b / abs(b))

        if d == igcd(c, d):
            a = a / d; b = b / d; c = c / d;
            if param:
                return (c*(x0 + b*t), c*(y0 - a*t))
            else:
                return (c * x0, c * y0)
        else:
            return (None, None)


def extended_euclid(a, b):

    if b == 0:
        return (1, 0, a)

    x0, y0, d = extended_euclid(b, a%b)
    x, y = y0, x0 - (a//b) * y0

    return x, y, d
