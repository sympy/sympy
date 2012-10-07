from sympy.rr.rl import rmid, glom, flatten, unpack, sort, distribute
from sympy import Basic

def test_rmid():
    rmzeros = rmid(lambda x: x == 0)
    assert rmzeros(Basic(0, 1)) == Basic(1)
    assert rmzeros(Basic(0, 0)) == Basic(0)
    assert rmzeros(Basic(2, 1)) == Basic(2, 1)

def test_glom():
    conglomerate = glom(lambda num, x: num * x)
    assert conglomerate(Basic(1, 2, 2)) == Basic(1, 4)
    conglomerate = glom(lambda num, x: x ** num)
    assert conglomerate(Basic(1, 3, 3)) == Basic(1, 9)

def test_flatten():
    assert flatten(Basic(1, 2, Basic(3, 4))) == Basic(1, 2, 3, 4)

def test_unpack():
    assert unpack(Basic(2)) == 2

def test_sort():
    assert sort(str)(Basic(3,1,2)) == Basic(1,2,3)

def test_distribute():
    class T1(Basic):        pass
    class T2(Basic):        pass

    distribute_t12 = distribute(T1, T2)
    assert distribute_t12(T1(1, 2, T2(3, 4), 5)) == \
            T2(T1(1, 2, 3, 5),
               T1(1, 2, 4, 5))

def test_distribute_add_mul():
    from sympy import Add, Mul, symbols
    x, y = symbols('x, y')
    expr = Mul(2, Add(x, y), evaluate=False)
    expected = Add(Mul(2, x), Mul(2, y))
    distribute_mul = distribute(Mul, Add)
    assert distribute_mul(expr) == expected

