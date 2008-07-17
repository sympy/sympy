
from sympy.utilities.decorator import threaded

from sympy import symbols, Eq, Matrix

def test_threaded():
    x, y = symbols('xy')

    @threaded()
    def function(expr, *args):
        return 2*expr + sum(args)

    assert function(Matrix([[x, y], [1, x]]), 1, 2) == \
        Matrix([[2*x+3, 2*y+3], [5, 2*x+3]])

    assert function(Eq(x, y), 1, 2) == Eq(2*x+3, 2*y+3)

    assert function([x, y], 1, 2) == [2*x+3, 2*y+3]
    assert function((x, y), 1, 2) == (2*x+3, 2*y+3)

    assert function(set([x, y]), 1, 2) == set([2*x+3, 2*y+3])

    @threaded()
    def function(expr, n):
        return expr**n

    assert function(x + y, 2) == x**2 + y**2
    assert function(x, 2) == x**2

    @threaded(use_add=False)
    def function(expr, n):
        return expr**n

    assert function(x + y, 2) == (x + y)**2
