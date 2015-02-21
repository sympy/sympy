from sympy import integrate, pi, sin, cos
from sympy.series.fourier import SeqFormula
from sympy.abc import x, y, z
from sympy.utilities.pytest import raises
from sympy.core.compatibility import range

def test_SeqFormula_generator():
    seq = SeqFormula(x*(x+1), x, (0, 5))
    assert [val for val in seq] == [n*(n+1) for n in range(6)]

    seq2 = SeqFormula(integrate(x * cos(y*x), (x, -pi, pi)), y, (0, 5))
    assert [val for val in seq2] == [integrate(x * cos(n*x), (x, -pi, pi)) for n in range(6)]

    seq3 = SeqFormula(integrate(x * sin(y*x), (x, -pi, pi)), y, (0, 5))
    assert [val for val in seq3] == [integrate(x * sin(n*x), (x, -pi, pi)) for n in range(6)]

def test_SeqFormula_coeff():
    seq = SeqFormula(sin(x), x, (0, 5))
    assert seq.coeff(3) == sin(3)

    raises(IndexError, lambda: seq.coeff(10))
    raises(IndexError, lambda: seq.coeff(-1))
