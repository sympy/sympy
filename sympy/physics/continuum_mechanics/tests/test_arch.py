from sympy.physics.continuum_mechanics.arch import Arch
from sympy import Symbol, simplify

x = Symbol('x')

def test_arch():
    a = Arch((0,0),(10,0),crown_x=5,crown_y=5)
    assert a.get_loads == {'distributed': {}, 'concentrated': {}}
    assert a.reaction_force == {Symbol('R_A_x'):0, Symbol('R_A_y'):0, Symbol('R_B_x'):0, Symbol('R_B_y'):0}
    assert a.supports == {'left':None, 'right':None}
    assert a.left_support == (0,0)
    assert a.right_support == (10,0)
    assert a.get_parabola_eqn == 5 - ((x-5)**2)/5

    a = Arch((0,0),(10,1),crown_x=6)
    assert simplify(a.get_parabola_eqn) == simplify(9/5 - (x - 6)**2/20)
