from sympy.parsing.sympy_parser import (
    parse_expr,
    standard_transformations,
    convert_xor,
    implicit_multiplication_application,
    implicit_multiplication,
    implicit_application,
    function_exponentiation,
    split_symbols
)


def test_implicit_multiplication():
    d = {
        '5x': '5*x',
        'abc': 'a*b*c',
        '3sin(x)': '3*sin(x)',
        '(x+1)(x+2)': '(x+1)*(x+2)',
        '(5 x**2)sin(x)': '(5*x**2)*sin(x)',
        '2 sin(x) cos(x)': '2*sin(x)*cos(x)'
    }
    transformations = standard_transformations + (convert_xor,)
    transformations2 = transformations + (split_symbols,
                                          implicit_multiplication)
    for e in d:
        implicit = parse_expr(e, transformations=transformations2)
        normal = parse_expr(d[e], transformations=transformations)
        assert(implicit == normal)


def test_implicit_application():
    d = {
        'factorial': 'factorial',
        'sin x': 'sin(x)',
        'tan y**3': 'tan(y**3)',
        'cos 2*x': 'cos(2*x)',
        '(cot)': 'cot'
    }
    transformations = standard_transformations + (convert_xor,)
    transformations2 = transformations + (implicit_application,)
    for e in d:
        implicit = parse_expr(e, transformations=transformations2)
        normal = parse_expr(d[e], transformations=transformations)
        assert(implicit == normal)


def test_function_exponentiation():
    d = {
        'sin**2(x)': 'sin(x)**2',
        'exp^y(z)': 'exp(z)^y'
    }
    transformations = standard_transformations + (convert_xor,)
    transformations2 = transformations + (function_exponentiation,)
    for e in d:
        implicit = parse_expr(e, transformations=transformations2)
        normal = parse_expr(d[e], transformations=transformations)
        assert(implicit == normal)


def test_all_implicit_steps():
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
        '(factorial)': 'factorial',
        'tan 3x': 'tan(3*x)'
    }
    transformations = standard_transformations + (convert_xor,)
    transformations2 = transformations + (implicit_multiplication_application,)
    for e in d:
        implicit = parse_expr(e, transformations=transformations2)
        normal = parse_expr(d[e], transformations=transformations)
        assert(implicit == normal)
