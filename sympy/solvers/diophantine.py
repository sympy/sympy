from sympy import degree_list
from sympy import igcd
from sympy import symbols, Wild


def diop_solve(eq):

    coeff, var, t = classify_diop(eq)

    if t == "linear":
        return diop_linear(coeff, var)


def classify_diop(eq):

    if max(degree_list(eq)) == 1:
        var = list(eq.free_symbols)
        diop_type = "linear"

        pattern = 0

        for i, v in enumerate(var):
            a = Wild("a" + str(i+1), exclude=var)
            pattern = pattern + a * v

        b = Wild('b', exclude=var)
        pattern = pattern - b

    return eq.match(pattern), var, diop_type


def diop_linear(coeff, var):

    if len(var) == 2:
        a1 = Wild('a1', exclude=var)
        a2 = Wild('a2', exclude=var)
        b = Wild('b', exclude=var)

        t = symbols('t', type = int)
        d, x0, y0 = extended_euclid(coeff[a1], coeff[a2])

        if d == igcd(coeff[b], d):
            return {var[0]:(coeff[b]/d) * (x0 + coeff[a2]*t), var[1]:(coeff[b]/d) * (y0 - coeff[a1]*t)}
        else:
            print "Not Solvable"


def extended_euclid(a, b):

    if b == 0:
        return (a, 1, 0)

    d, x0, y0 = extended_euclid(b, a%b)
    x, y = y0, x0 - (a//b) * y0

    return d, x, y
