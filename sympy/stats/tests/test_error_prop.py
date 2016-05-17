from sympy import symbols, exp
from sympy.stats.symbolic_probability import RandomSymbol, Variance
from sympy.stats.error_prop import variance_prop


def test_variance_prop():
    x, y, z = symbols('x y z')
    phi, t = consts = symbols('phi t')
    var_x = Variance(RandomSymbol(x))
    var_y = Variance(RandomSymbol(y))
    var_z = Variance(RandomSymbol(z))
    cases = {
        x + y: var_x + var_y,
        x + y + z: var_x + var_y + var_z,
        2*x: 4*var_x,
        x*y: var_x*y**2 + var_y*x**2,
        1/x: var_x/x**4,
        x/y: (var_x*y**2 + var_y*x**2)/y**4,
        exp(x): var_x*exp(2*x),
        exp(2*x): 4*var_x*exp(4*x),
        exp(-x*t): t**2*var_x*exp(-2*t*x),
        }
    for inp, out in cases.items():
        obs = variance_prop(inp, consts=consts)
        assert out == obs
