from sympy.matrices.dense import Matrix
from sympy.polys.polymatrix import PolyMatrix
from sympy.polys import Poly

from sympy import S, ZZ, QQ, EX, Symbol, symbols

from sympy.abc import x, y, z


def test_polymatrix():
    pm1 = PolyMatrix([[Poly(x**2, x), Poly(-x, x)], [Poly(x**3, x), Poly(-1 + x, x)]])
    v1 = PolyMatrix([[1, 0], [-1, 0]], ring='ZZ[x]')
    m1 = Matrix([[1, 0], [-1, 0]], ring='ZZ[x]')
    A = PolyMatrix([[Poly(x**2 + x, x), Poly(0, x)], \
                    [Poly(x**3 - x + 1, x), Poly(0, x)]])
    B = PolyMatrix([[Poly(x**2, x), Poly(-x, x)], [Poly(-x**2, x), Poly(x, x)]])

    assert A.ring == ZZ[x]
    assert v1.ring == ZZ[x]
    assert isinstance(pm1*v1, PolyMatrix)
    assert pm1*v1 == A
    assert pm1*m1 == A
    assert v1*pm1 == B

    pm2 = PolyMatrix([[Poly(x**2, x, domain='QQ'), Poly(0, x, domain='QQ'), Poly(-x**2, x, domain='QQ'), \
                    Poly(x**3, x, domain='QQ'), Poly(0, x, domain='QQ'), Poly(-x**3, x, domain='QQ')]])
    assert pm2.ring == QQ[x]
    v2 = PolyMatrix([1, 0, 0, 0, 0, 0], ring='ZZ[x]')
    m2 = Matrix([1, 0, 0, 0, 0, 0], ring='ZZ[x]')
    C = PolyMatrix([[Poly(x**2, x, domain='QQ')]])
    assert pm2*v2 == C
    assert pm2*m2 == C

    pm3 = PolyMatrix([[Poly(x**2, x), S.One]], ring='ZZ[x]')
    v3 = S.Half*pm3
    assert v3 == PolyMatrix([[Poly(S.Half*x**2, x, domain='QQ'), S.Half]], ring='EX')
    assert pm3*S.Half == v3
    assert v3.ring == EX

    pm4 = PolyMatrix([[Poly(x**2, x, domain='ZZ'), Poly(-x**2, x, domain='ZZ')]])
    v4 = Matrix([1, -1], ring='ZZ[x]')
    assert pm4*v4 == PolyMatrix([[Poly(2*x**2, x, domain='ZZ')]])

    assert len(PolyMatrix()) == 0
    assert PolyMatrix([1, 0, 0, 1])/(-1) == PolyMatrix([-1, 0, 0, -1])

    a = PolyMatrix([[1, 2], [3, 4]], ring=ZZ[x])
    b = PolyMatrix([[1, 2], [3, 4]], ring=ZZ[y])
    assert (a*b).ring == ZZ[x, y]


def test_matrix_add():
    a = PolyMatrix([[1, 2], [3, 4]], ring=ZZ[x])
    b = PolyMatrix([[1, 2], [3, 4]], ring=ZZ[y])
    c = a + b
    assert c == PolyMatrix([[2, 4], [6, 8]], ring=ZZ[x, y])
    assert c.ring == ZZ[x, y]


def test_scalar_mul():
    a = PolyMatrix([[1, 2], [3, 4]], ring=ZZ[x])
    assert (a * 2).ring == (2 * a).ring == ZZ[x]

    a = PolyMatrix([[Poly(1, x), Poly(2, x)], [Poly(3, x), Poly(4, x)]])
    b = PolyMatrix([[Poly(2, x), Poly(4, x)], [Poly(6, x), Poly(8, x)]])
    c = a * Poly(2, x)
    assert c == b
    assert c.ring == ZZ[x]
    c = Poly(2, x) * a
    assert c == b
    assert c.ring == ZZ[x]


def test_matrix_pow():
    a = PolyMatrix([[1, 2], [3, 4]], ring=ZZ)
    b = a**2
    assert b == PolyMatrix([[7, 10], [15, 22]], ring=ZZ)
    assert b.ring == ZZ

    a = PolyMatrix([[1, 2], [3, 4]], ring=ZZ[x])
    b = a**10
    assert b == PolyMatrix(
        [[4783807, 6972050], [10458075, 15241882]], ring=ZZ[x])
    assert b.ring == ZZ[x]


def test_hadamard_product():
    a = PolyMatrix([[1, 2], [3, 4]], ring='ZZ[x]')
    b = PolyMatrix([[1, 2], [3, 4]], ring='ZZ[y]')
    c = a.multiply_elementwise(b)
    assert c == PolyMatrix([[1, 4], [9, 16]], ring='ZZ[y]')
    assert c.ring == ZZ[x, y]


def test_matrix_slice_extract():
    a = PolyMatrix([
        [Poly(1, x), Poly(2, x), Poly(3, x)],
        [Poly(4, x), Poly(5, x), Poly(6, x)],
        [Poly(7, x), Poly(8, x), Poly(9, x)]],
        ring=ZZ[x, y])
    b = a[:1, :1]
    assert b == PolyMatrix([[Poly(1, x)]], ring=ZZ[x, y])
    assert b.ring == ZZ[x, y]

    b = a.extract([0], [0])
    assert b == PolyMatrix([[Poly(1, x)]], ring=ZZ[x, y])
    assert b.ring == ZZ[x, y]

    b = a.row(0)
    assert b == PolyMatrix(
        [[Poly(1, x), Poly(2, x), Poly(3, x)]], ring=ZZ[x, y])
    assert b.ring == ZZ[x, y]
    b = a.col(0)
    assert b == PolyMatrix(
        [[Poly(1, x)], [Poly(4, x)], [Poly(7, x)]], ring=ZZ[x, y])
    assert b.ring == ZZ[x, y]


def test_matrix_copy():
    a = PolyMatrix([
        [Poly(1, x), Poly(2, x), Poly(3, x)],
        [Poly(4, x), Poly(5, x), Poly(6, x)],
        [Poly(7, x), Poly(8, x), Poly(9, x)]],
        ring=ZZ[x, y])
    b = a.copy()
    assert b == PolyMatrix([
        [Poly(1, x), Poly(2, x), Poly(3, x)],
        [Poly(4, x), Poly(5, x), Poly(6, x)],
        [Poly(7, x), Poly(8, x), Poly(9, x)]],
        ring=ZZ[x, y])
    assert b.ring == ZZ[x, y]


def test_matrix_row_col():
    a = PolyMatrix([[Poly(1, x), Poly(2, x)]], ring=ZZ[x, y])
    b = PolyMatrix([[Poly(3, x), Poly(4, x)]], ring=ZZ[x, z])
    c = a.row_insert(1, b)
    assert c == PolyMatrix(
            [[Poly(1, x), Poly(2, x)], [Poly(3, x), Poly(4, x)]],
            ring=ZZ[x, y, z])
    assert c.ring == ZZ[x, y, z]

    c = a.col_join(b)
    assert c == PolyMatrix(
            [[Poly(1, x), Poly(2, x)], [Poly(3, x), Poly(4, x)]],
            ring=ZZ[x, y, z])
    assert c.ring == ZZ[x, y, z]

    c = PolyMatrix.vstack(a, b)
    assert c == PolyMatrix(
            [[Poly(1, x), Poly(2, x)], [Poly(3, x), Poly(4, x)]],
            ring=ZZ[x, y, z])
    assert c.ring == ZZ[x, y, z]

    a = PolyMatrix([[Poly(1, x)], [Poly(3, x)]], ring=ZZ[x, y])
    b = PolyMatrix([[Poly(2, x)], [Poly(4, x)]], ring=ZZ[x, z])
    c = a.col_insert(1, b)
    assert c == PolyMatrix(
            [[Poly(1, x), Poly(2, x)], [Poly(3, x), Poly(4, x)]],
            ring=ZZ[x, y, z])
    assert c.ring == ZZ[x, y, z]

    c = a.row_join(b)
    assert c == PolyMatrix(
            [[Poly(1, x), Poly(2, x)], [Poly(3, x), Poly(4, x)]],
            ring=ZZ[x, y, z])
    assert c.ring == ZZ[x, y, z]

    c = PolyMatrix.hstack(a, b)
    assert c == PolyMatrix(
            [[Poly(1, x), Poly(2, x)], [Poly(3, x), Poly(4, x)]],
            ring=ZZ[x, y, z])
    assert c.ring == ZZ[x, y, z]


def test_diagonal():
    a = PolyMatrix([
        [Poly(1, x), Poly(2, x), Poly(3, x)],
        [Poly(4, x), Poly(5, x), Poly(6, x)],
        [Poly(7, x), Poly(8, x), Poly(9, x)]],
        ring=ZZ[x, y])
    b = a.diagonal(0)
    assert b == PolyMatrix(
        [[Poly(1, x), Poly(5, x), Poly(9, x)]], ring=ZZ[x, y])
    b = a.diagonal(1)
    assert b == PolyMatrix(
        [[Poly(2, x), Poly(6, x)]], ring=ZZ[x, y])
    b = a.diagonal(-2)
    assert b == PolyMatrix(
        [[Poly(7, x)]], ring=ZZ[x, y])


def test_get_diag_blocks():
    a = PolyMatrix([
        [Poly(1, x), Poly(0, x), Poly(0, x)],
        [Poly(0, x), Poly(5, x), Poly(6, x)],
        [Poly(0, x), Poly(8, x), Poly(9, x)]],
        ring=ZZ[x, y])
    b, c = a.get_diag_blocks()
    assert b == PolyMatrix([[Poly(1, x)]], ring=ZZ[x, y])
    assert b.ring == ZZ[x, y]
    assert c == PolyMatrix(
        [[Poly(5, x), Poly(6, x)], [Poly(8, x), Poly(9, x)]], ring=ZZ[x, y])
    assert c.ring == ZZ[x, y]


def test_reshape():
    a = PolyMatrix([
        [Poly(1, x), Poly(2, x), Poly(3, x)],
        [Poly(4, x), Poly(5, x), Poly(6, x)]],
        ring=ZZ[x, y])
    b = a.reshape(3, 2)
    assert b == PolyMatrix([
        [Poly(1, x), Poly(2, x)],
        [Poly(3, x), Poly(4, x)],
        [Poly(5, x), Poly(6, x)]],
        ring=ZZ[x, y])
    assert b.ring == ZZ[x, y]


def test_charpoly():
    N = lambda n: Poly(n, x, domain=ZZ)
    l = Symbol('l')
    m = PolyMatrix(
        [[N(64 + 11*x), N(1 + 7*x), N(82 + 83*x)],
        [N(75 + 46*x), N(35 + 54*x), N(63 + 8*x)],
        [N(88 + 19*x), N(3 + 99*x), N(73 + 79*x)]],
    )
    assert m.charpoly(l) == \
        Poly(
            l**3 +
            (-144*x - 172)*l**2 +
            (3038*x**2 + 713*x + 1987)*l -
            306664*x**3 - 639500*x**2 + 13322*x + 82617,
            l,
            domain=ZZ[x])


def test_det():
    N = lambda n: Poly(n, x, domain=ZZ)
    m = PolyMatrix(
        [[N(64 + 11*x), N(1 + 7*x), N(82 + 83*x), N(55 + 99*x)],
        [N(75 + 46*x), N(35 + 54*x), N(63 + 8*x), N(23 + 40*x)],
        [N(88 + 19*x), N(3 + 99*x), N(73 + 79*x), N(21 + 34*x)],
        [N(48 + 83*x), N(42 + 82*x), N(25 + 69*x), N(56 + 54*x)]]
    )
    assert m.det(method='bareiss') == \
        Poly(
            -26849754*x**4 + 34317848*x**3 + 57571980*x**2 + 4558659*x -
            4955490, x, domain=ZZ)
    assert m.det(method='berkowitz') == \
        Poly(
            -26849754*x**4 + 34317848*x**3 + 57571980*x**2 + 4558659*x -
            4955490, x, domain=ZZ)


def test_diff():
    a, b, c, d, e, f, g, h = symbols('a:h')
    x = Symbol('x')
    domain = ZZ[a, b, c, d, e, f, g, h]
    m = PolyMatrix(
        [[Poly(a*x + b, x, domain=domain), Poly(c*x + d, x, domain=domain)],
        [Poly(e*x + f, x, domain=domain), Poly(g*x + h, x, domain=domain)]],
    )
    assert m.diff(x) == PolyMatrix(
        [[Poly(a, x, domain=domain), Poly(c, x, domain=domain)],
        [Poly(e, x, domain=domain), Poly(g, x, domain=domain)]],
    )


def test_integrate():
    a, b, c, d, e, f, g, h = symbols('a:h')
    x = Symbol('x')
    domain = ZZ[a, b, c, d, e, f, g, h]
    m = PolyMatrix([
        [Poly(a*x + b, x, domain=domain), Poly(c*x + d, x, domain=domain)],
        [Poly(e*x + f, x, domain=domain), Poly(g*x + h, x, domain=domain)]],
    )
    domain = ZZ.frac_field(a, b, c, d, e, f, g, h)
    assert m.integrate(x) == PolyMatrix([
        [Poly(a/2*x**2 + b*x, x, domain=domain),
         Poly(c/2*x**2 + d*x, x, domain=domain)],
        [Poly(e/2*x**2 + f*x, x, domain=domain),
         Poly(g/2*x**2 + h*x, x, domain=domain)]],
    )


def test_vec():
    a = PolyMatrix([
        [Poly(1, x), Poly(2, x)],
        [Poly(3, x), Poly(4, x)]],
        ring=ZZ[x, y])
    b = a.vec()
    assert b == PolyMatrix(
        [[Poly(1, x)], [Poly(3, x)], [Poly(2, x)], [Poly(4, x)]],
        ring=ZZ[x, y])
    assert b.ring == ZZ[x, y]


def test_jordan_block():
    m = PolyMatrix.jordan_block(3, Poly(-1, x))
    assert m == PolyMatrix([
        [Poly(-1, x), Poly(1, x), Poly(0, x)],
        [Poly(0, x), Poly(-1, x), Poly(1, x)],
        [Poly(0, x), Poly(0, x), Poly(-1, x)]])
    assert m.ring == ZZ[x]


def test_companion():
    m = PolyMatrix.companion(Poly([1, 2, 3, 4], x))
    assert m == PolyMatrix([
        [Poly(0, x), Poly(0, x), Poly(-4, x)],
        [Poly(1, x), Poly(0, x), Poly(-3, x)],
        [Poly(0, x), Poly(1, x), Poly(-2, x)]], ring=ZZ[x])
    assert m.ring == ZZ[x]


def test_transpose():
    a = PolyMatrix(
        [[Poly(1, x), Poly(2, x)], [Poly(3, x), Poly(4, x)]],
        ring=ZZ[x, y])
    b = a.T
    assert b == PolyMatrix(
        [[Poly(1, x), Poly(3, x)], [Poly(2, x), Poly(4, x)]],
        ring=ZZ[x, y])
    assert b.ring == ZZ[x, y]


def test_applyfunc():
    a = PolyMatrix(
        [[Poly(1, x), Poly(2, x)], [Poly(3, x), Poly(4, x)]],
        ring=ZZ[x, y]
    )
    b = a.applyfunc(lambda x: x+1)
    assert b == PolyMatrix(
        [[Poly(2, x), Poly(3, x)], [Poly(4, x), Poly(5, x)]],
        ring=ZZ[x, y])
    assert b.ring == ZZ[x, y]


def test_symmetry_properties():
    m = PolyMatrix([[Poly(0, x), Poly(1, x)], [Poly(1, x), Poly(0, x)]])
    assert m.is_symmetric() is True
    assert m.is_anti_symmetric() is False

    m = PolyMatrix([[Poly(0, x), Poly(1, x)], [Poly(-1, x), Poly(0, x)]])
    assert m.is_symmetric() is False
    assert m.is_anti_symmetric() is True


def test_diagonal_properties():
    m = PolyMatrix([[Poly(1, x), Poly(0, x)], [Poly(0, x), Poly(1, x)]])
    assert m.is_diagonal() is True
    assert m.is_Identity is True

    m = PolyMatrix([[Poly(1, x), Poly(0, x)], [Poly(0, x), Poly(0, x)]])
    assert m.is_diagonal() is True
    assert m.is_Identity is False

    m = PolyMatrix([[Poly(1, x), Poly(1, x)], [Poly(0, x), Poly(1, x)]])
    assert m.is_diagonal() is False
    assert m.is_Identity is False


def test_triangular_properties():
    m = PolyMatrix([
        [Poly(1, x), Poly(1, x), Poly(1, x)],
        [Poly(0, x), Poly(1, x), Poly(1, x)],
        [Poly(0, x), Poly(0, x), Poly(1, x)]])
    assert m.is_upper is True
    assert m.is_lower is False
    assert m.is_upper_hessenberg is True
    assert m.is_lower_hessenberg is False

    m = PolyMatrix([
        [Poly(1, x), Poly(1, x), Poly(1, x)],
        [Poly(1, x), Poly(1, x), Poly(1, x)],
        [Poly(0, x), Poly(1, x), Poly(1, x)]])
    assert m.is_upper is False
    assert m.is_lower is False
    assert m.is_upper_hessenberg is True
    assert m.is_lower_hessenberg is False

    m = PolyMatrix([
        [Poly(1, x), Poly(0, x), Poly(0, x)],
        [Poly(1, x), Poly(1, x), Poly(0, x)],
        [Poly(1, x), Poly(1, x), Poly(1, x)]])
    assert m.is_upper is False
    assert m.is_lower is True
    assert m.is_upper_hessenberg is False
    assert m.is_lower_hessenberg is True

    m = PolyMatrix([
        [Poly(1, x), Poly(1, x), Poly(0, x)],
        [Poly(1, x), Poly(1, x), Poly(1, x)],
        [Poly(1, x), Poly(1, x), Poly(1, x)]])
    assert m.is_upper is False
    assert m.is_lower is False
    assert m.is_upper_hessenberg is False
    assert m.is_lower_hessenberg is True


def test_is_zero_matrix():
    m = PolyMatrix([[Poly(0, x), Poly(0, x)], [Poly(0, x), Poly(0, x)]])
    assert m.is_zero_matrix is True
    m = PolyMatrix([[Poly(0, x), Poly(0, x)], [Poly(0, x), Poly(1, x)]])
    assert m.is_zero_matrix is False


def test_permute_rows_cols():
    a = PolyMatrix(
        [[Poly(1, x), Poly(2, x)], [Poly(3, x), Poly(4, x)]], ring=ZZ[x, y])
    b = a.permute_rows([1, 0])
    assert b == PolyMatrix(
        [[Poly(3, x), Poly(4, x)], [Poly(1, x), Poly(2, x)]], ring=ZZ[x, y])
    assert b.ring == ZZ[x, y]
    b = a.permute_cols([1, 0])
    assert b == PolyMatrix(
        [[Poly(2, x), Poly(1, x)], [Poly(4, x), Poly(3, x)]], ring=ZZ[x, y])
    assert b.ring == ZZ[x, y]


def test_echelon_form():
    X = Symbol('X')
    m = PolyMatrix(8, 8, [Poly(x + i, x) for i in range(8*8)])

    N = lambda n: Poly(n, X, domain=m.ring.get_field())
    assert m._to_field(dummy=X).echelon_form() == \
        Matrix([
            [N(x), N(x+1), N(x+2), N(x+3), N(x+4), N(x+5), N(x+6), N(x+7)],
            [N(0), N(-8), N(-16), N(-24), N(-32), N(-40), N(-48), N(-56)],
            [N(0), N(0), N(0), N(0), N(0), N(0), N(0), N(0)],
            [N(0), N(0), N(0), N(0), N(0), N(0), N(0), N(0)],
            [N(0), N(0), N(0), N(0), N(0), N(0), N(0), N(0)],
            [N(0), N(0), N(0), N(0), N(0), N(0), N(0), N(0)],
            [N(0), N(0), N(0), N(0), N(0), N(0), N(0), N(0)],
            [N(0), N(0), N(0), N(0), N(0), N(0), N(0), N(0)]])


def test_rref():
    X = Symbol('X')
    m = PolyMatrix(8, 8, [Poly(x + i, x) for i in range(8*8)])

    N = lambda n: Poly(n, X, domain=m.ring.get_field())
    assert m._to_field(dummy=X).rref() == (
        Matrix([
            [N(1), N(0), N(-1), N(-2), N(-3), N(-4), N(-5), N(-6)],
            [N(0), N(1), N(2), N(3), N(4), N(5), N(6), N(7)],
            [N(0), N(0), N(0), N(0), N(0), N(0), N(0), N(0)],
            [N(0), N(0), N(0), N(0), N(0), N(0), N(0), N(0)],
            [N(0), N(0), N(0), N(0), N(0), N(0), N(0), N(0)],
            [N(0), N(0), N(0), N(0), N(0), N(0), N(0), N(0)],
            [N(0), N(0), N(0), N(0), N(0), N(0), N(0), N(0)],
            [N(0), N(0), N(0), N(0), N(0), N(0), N(0), N(0)]]),
        (0, 1))


def test_rank():
    m = PolyMatrix([
        [Poly(9850*x**2 + 24970*x + 15845, x, domain='QQ'),
         Poly(7140*x**2 + 9870*x + 863, x, domain='QQ'),
         Poly(6850*x**2 + 14069*x + 6595, x, domain='QQ')],
        [Poly(7140*x**2 + 9870*x + 863, x, domain='QQ'),
         Poly(6370*x**2 + 1274*x + 65, x, domain='QQ'),
         Poly(4746*x**2 + 2289*x + 97, x, domain='QQ')],
        [Poly(6850*x**2 + 14069*x + 6595, x, domain='QQ'),
         Poly(4746*x**2 + 2289*x + 97, x, domain='QQ'),
         Poly(4804*x**2 + 8292*x + 6565, x, domain='QQ')]])
    assert m._to_field().rank() == 2
