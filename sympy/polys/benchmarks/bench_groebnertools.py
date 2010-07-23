"""Benchmark for the Groebner bases algorithm. """

from sympy import Symbol, groebner

V = range(1, 12+1)
E = [(1,2),(2,3),(1,4),(1,6),(1,12),(2,5),(2,7),(3,8),
(3,10),(4,11),(4,9),(5,6),(6,7),(7,8),(8,9),(9,10),
(10,11),(11,12),(5,12),(5,9),(6,10),(7,11),(8,12)]

Vx = [ Symbol('x' + str(i)) for i in V ]
Ex = [ (Vx[i-1], Vx[j-1]) for i, j in E ]

F3 = [ x**3 - 1 for x in Vx ]
Fg = [ x**2 + x*y + y**2 for x, y in Ex ]

x3, x4 = Vx[2], Vx[3]

F_1 = F3 + Fg
F_2 = F3 + Fg + [x3**2 + x3*x4 + x4**2]

def time_vertex_color_12_vertices_23_edges():
    assert groebner(F_1, Vx) != [1]

def time_vertex_color_12_vertices_24_edges():
    assert groebner(F_2, Vx) == [1]

