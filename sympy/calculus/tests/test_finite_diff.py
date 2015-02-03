from sympy.core.compatibility import range
from sympy import S, symbols, Function
from sympy.calculus.finite_diff import (
    apply_finite_diff, finite_diff_weights, as_finite_diff
)


def test_apply_finite_diff():
    x, h = symbols('x h')
    f = Function('f')
    assert (apply_finite_diff(1, [x-h, x+h], [f(x-h), f(x+h)], x) -
            (f(x+h)-f(x-h))/(2*h)).simplify() == 0

    assert (apply_finite_diff(1, [5, 6, 7], [f(5), f(6), f(7)], 5) -
            (-S(3)/2*f(5) + 2*f(6) - S(1)/2*f(7))).simplify() == 0


def test_finite_diff_weights():

    d = finite_diff_weights(1, [5, 6, 7], 5)
    assert d[1][2] == [-S(3)/2, 2, -S(1)/2]

    # Table 1, p. 702 in doi:10.1090/S0025-5718-1988-0935077-0
    # --------------------------------------------------------
    # x = [[0], [-1, 0, 1], ...]
    xl = [[j for j in range(-i, i+1)] for i in range(0, 5)]

    # d holds all coefficients
    d = [finite_diff_weights({0: 0, 1: 2, 2: 4, 3: 4, 4: 4}[i],
                             xl[i], 0) for i in range(5)]

    # Zeroeth derivative
    assert d[0][0][0] == [S(1)]

    # First derivative
    assert d[1][1][2] == [-S(1)/2, S(0), S(1)/2]
    assert d[2][1][4] == [S(1)/12, -S(2)/3, S(0), S(2)/3, -S(1)/12]
    assert d[3][1][6] == [-S(1)/60, S(3)/20, -S(3)/4, S(0), S(3)/4, -S(3)/20,
                          S(1)/60]
    assert d[4][1][8] == [S(1)/280, -S(4)/105, S(1)/5, -S(4)/5, S(0), S(4)/5,
                          -S(1)/5, S(4)/105, -S(1)/280]

    # Second derivative
    assert d[1][2][2] == [S(1), -S(2), S(1)]
    assert d[2][2][4] == [-S(1)/12, S(4)/3, -S(5)/2, S(4)/3, -S(1)/12]
    assert d[3][2][6] == [S(1)/90, -S(3)/20, S(3)/2, -S(49)/18, S(3)/2,
                          -S(3)/20, S(1)/90]
    assert d[4][2][8] == [-S(1)/560, S(8)/315, -S(1)/5, S(8)/5, -S(205)/72,
                          S(8)/5, -S(1)/5, S(8)/315, -S(1)/560]

    # Third derivative
    assert d[2][3][4] == [-S(1)/2, S(1), S(0), -S(1), S(1)/2]
    assert d[3][3][6] == [S(1)/8, -S(1), S(13)/8, S(0), -S(13)/8, S(1),
                          -S(1)/8]
    assert d[4][3][8] == [-S(7)/240, S(3)/10, -S(169)/120, S(61)/30, S(0),
                          -S(61)/30, S(169)/120, -S(3)/10, S(7)/240]

    # Fourth derivative
    assert d[2][4][4] == [S(1), -S(4), S(6), -S(4), S(1)]
    assert d[3][4][6] == [-S(1)/6, S(2), -S(13)/2, S(28)/3, -S(13)/2, S(2),
                          -S(1)/6]
    assert d[4][4][8] == [S(7)/240, -S(2)/5, S(169)/60, -S(122)/15, S(91)/8,
                          -S(122)/15, S(169)/60, -S(2)/5, S(7)/240]

    # Table 2, p. 703 in doi:10.1090/S0025-5718-1988-0935077-0
    # --------------------------------------------------------
    xl = [[j/S(2) for j in list(range(-i*2+1, 0, 2))+list(range(1, i*2+1, 2))]
          for i in range(1, 5)]

    # d holds all coefficients
    d = [finite_diff_weights({0: 1, 1: 2, 2: 4, 3: 4}[i], xl[i], 0) for
         i in range(4)]

    # Zeroth derivative
    assert d[0][0][1] == [S(1)/2, S(1)/2]
    assert d[1][0][3] == [-S(1)/16, S(9)/16, S(9)/16, -S(1)/16]
    assert d[2][0][5] == [S(3)/256, -S(25)/256, S(75)/128, S(75)/128,
                          -S(25)/256, S(3)/256]
    assert d[3][0][7] == [-S(5)/2048, S(49)/2048, -S(245)/2048, S(1225)/2048,
                          S(1225)/2048, -S(245)/2048, S(49)/2048, -S(5)/2048]

    # First derivative
    assert d[0][1][1] == [-S(1), S(1)]
    assert d[1][1][3] == [S(1)/24, -S(9)/8, S(9)/8, -S(1)/24]
    assert d[2][1][5] == [-S(3)/640, S(25)/384, -S(75)/64, S(75)/64,
                          -S(25)/384, S(3)/640]
    assert d[3][1][7] == [S(5)/7168, -S(49)/5120,  S(245)/3072, S(-1225)/1024,
                          S(1225)/1024, -S(245)/3072, S(49)/5120, -S(5)/7168]

    # Reasonably the rest of the table is also correct... (testing of that
    # deemed excessive at the moment)


def test_as_finite_diff():
    x, h = symbols('x h')
    f = Function('f')

    # Central 1st derivative at gridpoint
    assert (as_finite_diff(f(x).diff(x), [x-2, x-1, x, x+1, x+2]) -
            (S(1)/12*(f(x-2)-f(x+2)) + S(2)/3*(f(x+1)-f(x-1)))).simplify() == 0

    # Central 1st derivative "half-way"
    assert (as_finite_diff(f(x).diff(x)) -
            (f(x + S(1)/2)-f(x - S(1)/2))).simplify() == 0
    assert (as_finite_diff(f(x).diff(x), h) -
            (f(x + h/S(2))-f(x - h/S(2)))/h).simplify() == 0
    assert (as_finite_diff(f(x).diff(x), [x - 3*h, x-h, x+h, x + 3*h]) -
            (S(9)/(8*2*h)*(f(x+h) - f(x-h)) +
             S(1)/(24*2*h)*(f(x - 3*h) - f(x + 3*h)))).simplify() == 0

    # One sided 1st derivative at gridpoint
    assert (as_finite_diff(f(x).diff(x), [0, 1, 2], 0) -
            (-S(3)/2*f(0) + 2*f(1) - f(2)/2)).simplify() == 0
    assert (as_finite_diff(f(x).diff(x), [x, x+h], x) -
            (f(x+h) - f(x))/h).simplify() == 0
    assert (as_finite_diff(f(x).diff(x), [x-h, x, x+h], x-h) -
            (-S(3)/(2*h)*f(x-h) + 2/h*f(x) -
             S(1)/(2*h)*f(x+h))).simplify() == 0

    # One sided 1st derivative "half-way"
    assert (as_finite_diff(f(x).diff(x), [x-h, x+h, x + 3*h, x + 5*h, x + 7*h])
            - 1/(2*h)*(-S(11)/(12)*f(x-h) + S(17)/(24)*f(x+h)
                       + S(3)/8*f(x + 3*h) - S(5)/24*f(x + 5*h)
                       + S(1)/24*f(x + 7*h))).simplify() == 0

    # Central 2nd derivative at gridpoint
    assert (as_finite_diff(f(x).diff(x, 2), [x-h, x, x+h]) -
            h**-2 * (f(x-h) + f(x+h) - 2*f(x))).simplify() == 0

    assert (as_finite_diff(f(x).diff(x, 2), [x - 2*h, x-h, x, x+h, x + 2*h]) -
            h**-2 * (-S(1)/12*(f(x - 2*h) + f(x + 2*h)) +
                     S(4)/3*(f(x+h) + f(x-h)) - S(5)/2*f(x))).simplify() == 0

    # Central 2nd derivative "half-way"
    assert (as_finite_diff(f(x).diff(x, 2), [x - 3*h, x-h, x+h, x + 3*h]) -
            (2*h)**-2 * (S(1)/2*(f(x - 3*h) + f(x + 3*h)) -
                         S(1)/2*(f(x+h) + f(x-h)))).simplify() == 0

    # One sided 2nd derivative at gridpoint
    assert (as_finite_diff(f(x).diff(x, 2), [x, x+h, x + 2*h, x + 3*h]) -
            h**-2 * (2*f(x) - 5*f(x+h) +
                     4*f(x+2*h) - f(x+3*h))).simplify() == 0

    # One sided 2nd derivative at "half-way"
    assert (as_finite_diff(f(x).diff(x, 2), [x-h, x+h, x + 3*h, x + 5*h]) -
            (2*h)**-2 * (S(3)/2*f(x-h) - S(7)/2*f(x+h) + S(5)/2*f(x + 3*h) -
                         S(1)/2*f(x + 5*h))).simplify() == 0

    # Central 3rd derivative at gridpoint
    assert (as_finite_diff(f(x).diff(x, 3)) -
            (-f(x - 3/S(2)) + 3*f(x - 1/S(2)) -
             3*f(x + 1/S(2)) + f(x + 3/S(2)))).simplify() == 0

    assert (as_finite_diff(
        f(x).diff(x, 3), [x - 3*h, x - 2*h, x-h, x, x+h, x + 2*h, x + 3*h]) -
        h**-3 * (S(1)/8*(f(x - 3*h) - f(x + 3*h)) - f(x - 2*h) +
                 f(x + 2*h) + S(13)/8*(f(x-h) - f(x+h)))).simplify() == 0

    # Central 3rd derivative at "half-way"
    assert (as_finite_diff(f(x).diff(x, 3), [x - 3*h, x-h, x+h, x + 3*h]) -
            (2*h)**-3 * (f(x + 3*h)-f(x - 3*h) +
                         3*(f(x-h)-f(x+h)))).simplify() == 0

    # One sided 3rd derivative at gridpoint
    assert (as_finite_diff(f(x).diff(x, 3), [x, x+h, x + 2*h, x + 3*h]) -
            h**-3 * (f(x + 3*h)-f(x) + 3*(f(x+h)-f(x + 2*h)))).simplify() == 0

    # One sided 3rd derivative at "half-way"
    assert (as_finite_diff(f(x).diff(x, 3), [x-h, x+h, x + 3*h, x + 5*h]) -
            (2*h)**-3 * (f(x + 5*h)-f(x-h) +
                         3*(f(x+h)-f(x + 3*h)))).simplify() == 0
