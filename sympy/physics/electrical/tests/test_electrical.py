from sympy import symbols, simplify, oo, sqrt
from sympy.physics.electrical import EField, ParticleCharge, check_conservative
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics.essential import _check_vector

N = ReferenceFrame('N')

def test_EField():
    """
    Test working of EField.
    """
    
    x, y, z = symbols('x y z')
    field = EField(x * N.x + y * N.y + z * N.z)
    assert check_conservative(field.get_vector())
    assert simplify(field.magnitude() - sqrt(x ** 2 + y ** 2 + z ** 2)) == 0
    assert (field + EField(N.x + N.y + N.z)).get_vector() == (x + 1) * N.x + (y + 1) * N.y + (z + 1) * N.z
    assert (field - EField(N.x + N.y + N.z)).get_vector() == (x - 1) * N.x + (y - 1) * N.y + (z - 1) * N.z
    assert (field * 3).__class__ == EField
    assert (field / 2.0).__class__ == EField
    field1 = EField(- x * N.x - y * N.y - z * N.z)
    assert field == -field1
    field = EField(3*N.x)
    o = Point('o')
    p = Point('p')
    q = Point('q')
    p.set_pos(o,3*N.x + 4*N.y)
    q.set_pos(o,2*N.x + 7*N.y)
    assert field.potential_difference(p, q) == 3
    assert field.potential_difference(p, p) == 0

def test_ParticleCharge():
    """
    Test working of ParticleCharge.
    """
    
    o = Point('o')
    p = Point('p')
    q = Point('q')
    r = Point('r')
    p.set_pos(o, N.x)
    q.set_pos(o, -N.x)
    r.set_pos(o, N.y)
    P1 = ParticleCharge('P1', p, 1)
    Q1 = ParticleCharge('Q1', q, 1)
    assert P1.get_charge() == 1
    #Test whether the vectorial sum of fields at r due to p and q
    #is parallel to N.y
    assert (P1.field_at(r) + Q1.field_at(r)).get_vector().dot(N.x) == 0
    assert P1.field_at(p) == oo
    assert P1.potential_at(p) == oo
    P1.set_charge(-1)
    assert P1.potential_at(o) == -Q1.potential_at(o)

def test_check_conservative():
    """
    Test check_conservative
    """
    
    assert check_conservative(N.x)
    x, y, z = symbols('x y z')
    v = - y /(x ** 2 + y ** 2) * N.x + x /(x ** 2 + y ** 2) * N.y
    assert check_conservative(v)
    assert check_conservative(y * z * N.x) == False


    
    
