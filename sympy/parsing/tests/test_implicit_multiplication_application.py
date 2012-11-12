from sympy.parsing.sympy_parser import (
    parse_expr,
    standard_transformations,
    convert_xor,
    implicit_multiplication_application
)

def test_implicit_multiplication_application():
    d = {
        '2x': '2*x',  # implicit multiplication
        'x y': 'x*y',
        'xy': 'x*y',
        'sin x': 'sin(x)',  # add parentheses
        '2sin x': '2*sin(x)',
        'x y z': 'x*y*z',
        'sin(2 * 3x)': 'sin(2 * 3 * x)',
        'sin(x) (1 + cos(x))': 'sin(x) * (1 + cos(x))',
        '(x + 2) sin(x)': '(x + 2) * sin(x)',
        '(x + 2) sin x': '(x + 2) * sin(x)',
        'sin(sin x)': 'sin(sin(x))',
        'sin x!': 'sin(factorial(x))',
        'sin x!!': 'sin(factorial2(x))',
        'factorial': 'factorial',  # don't apply a bare function
        'x sin x': 'x * sin(x)',  # both application and multiplication
        'xy sin x': 'x * y * sin(x)',
        '(x+2)(x+3)': '(x + 2) * (x+3)',
        'x**2 + 2xy + y**2': 'x**2 + 2 * x * y + y**2',  # split the xy
        'pi': 'pi',  # don't mess with constants
        'None': 'None',
        'ln sin x': 'ln(sin(x))',  # multiple implicit function applications
        'factorial': 'factorial',  # don't add parentheses
        'sin x**2': 'sin(x**2)',  # implicit application to an exponential
        'alpha': 'Symbol("alpha")',  # don't split Greek letters/subscripts
        'x_2': 'Symbol("x_2")',
        'sin^2 x**2': 'sin(x**2)**2',  # function raised to a power
        'sin**3(x)': 'sin(x)**3',
        '(factorial)': 'factorial'
    }
    transformations = standard_transformations + (convert_xor,)
    transformations2 = transformations + (implicit_multiplication_application,)
    for e in d:
        implicit = parse_expr(e, transformations=transformations2)
        normal = parse_expr(d[e], transformations=transformations)
        assert(implicit == normal)
