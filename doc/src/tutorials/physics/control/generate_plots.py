from sympy import Matrix, laplace_transform, inverse_laplace_transform, exp, cos, sqrt, sin
from sympy.abc import s, t
from sympy.physics.control import *

def main_q3():
    g =  Matrix([[exp(-t)*(1 - t), exp(-2*t)], [5*exp((-2*t))-exp((-t)), (cos((sqrt(3)*t)/2) - 3*sqrt(3)*sin((sqrt(3)*t)/2))*exp(-t/2)]])
    G = g.applyfunc(lambda a: laplace_transform(a, t, s)[0])
    G = TransferFunctionMatrix.from_Matrix(G, s)
    return G

def q3_3():
    G = main_q3()
    pole_zero_plot(G[0, 0])

def q3_4():
    G = main_q3()
    tf1 = G[0, 0]
    step_response_plot(tf1)

def q3_5_1():
    G = main_q3()
    tf2 = G[0, 1]
    bode_magnitude_plot(tf2)

def q3_5_2():
    G = main_q3()
    tf2 = G[0, 1]
    bode_magnitude_plot(tf2)

def q5():
    G1 = TransferFunction(1, 10 + s, s)
    G2 = TransferFunction(1, 1 + s, s)
    G3 = TransferFunction(1 + s**2, 4 + 4*s + s**2, s)
    G4 = TransferFunction(1 + s, 6 + s, s)
    H1 = TransferFunction(1 + s, 2 + s, s)
    H2 = TransferFunction(2*(6 + s), 1 + s, s)
    H3 = TransferFunction(1, 1, s)
    sys1 = Series(G3, G4)
    sys2 = Feedback(sys1, H1, 1).doit()
    sys3 = Series(G2, sys2)
    sys4 = Feedback(sys3, H2).doit()
    sys5 = Series(G1, sys4)
    sys6 = Feedback(sys5, H3)
    sys6 = sys6.doit(cancel=True, expand=True)
    pole_zero_plot(sys6)
