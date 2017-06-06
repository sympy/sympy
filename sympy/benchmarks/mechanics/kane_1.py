# This is taken from the example in KanesMethod docstring
from sympy import symbols
import sympy.physics.mechanics as me

q, u = me.dynamicsymbols('q u')
qd, ud = me.dynamicsymbols('q u', 1)
m, c, k = symbols('m c k')
N = me.ReferenceFrame('N')
P = me.Point('P')
P.set_vel(N, u * N.x)

kd = [qd - u]
FL = [(P, (-k * q - c * u) * N.x)]
pa = me.Particle('pa', P, m)
BL = [pa]


def kmethod():
    KM = me.KanesMethod(N, q_ind=[q], u_ind=[u], kd_eqs=kd)
    (fr, frstar) = KM.kanes_equations(BL, FL)
