"""Test sparse polynomials. """

from sympy.polys.rings import ring, xring, PolyRing, PolyElement
from sympy.polys.fields import field
from sympy.polys.domains import ZZ, QQ, RR
from sympy.polys.monomialtools import lex, grlex
from sympy.polys.polyerrors import GeneratorsError, GeneratorsNeeded

from sympy.utilities.pytest import raises
from sympy.core import Symbol, symbols
from sympy import sqrt

def test_PolyRing___init__():
    x, y, z = map(Symbol, "xyz")

    assert len(PolyRing("x,y,z", ZZ, lex).gens) == 3
    assert len(PolyRing(x, ZZ, lex).gens) == 1
    assert len(PolyRing(("x", "y", "z"), ZZ, lex).gens) == 3
    assert len(PolyRing((x, y, z), ZZ, lex).gens) == 3

    raises(GeneratorsNeeded, lambda: PolyRing([], ZZ, lex))
    raises(GeneratorsError, lambda: PolyRing(0, ZZ, lex))

def test_PolyRing___hash__():
    R, x, y, z = ring("x,y,z", QQ)
    assert hash(R)

def test_PolyRing___eq__():
    assert ring("x,y,z", QQ)[0] == ring("x,y,z", QQ)[0]
    assert ring("x,y,z", QQ)[0] != ring("x,y,z", ZZ)[0]
    assert ring("x,y,z", ZZ)[0] != ring("x,y,z", QQ)[0]
    assert ring("x,y,z", QQ)[0] != ring("x,y", QQ)[0]
    assert ring("x,y", QQ)[0] != ring("x,y,z", QQ)[0]

def test_PolyElement___hash__():
    R, x, y, z = ring("x,y,z", QQ)
    assert hash(x*y*z)

def test_PolyElement___eq__():
    R, x, y = ring("x,y", ZZ, lex)

    assert ((x*y + 5*x*y) == 6) == False
    assert ((x*y + 5*x*y) == 6*x*y) == True
    assert (6 == (x*y + 5*x*y)) == False
    assert (6*x*y == (x*y + 5*x*y)) == True

    assert ((x*y - x*y) == 0) == True
    assert (0 == (x*y - x*y)) == True

    assert ((x*y - x*y) == 1) == False
    assert (1 == (x*y - x*y)) == False

    assert ((x*y - x*y) == 1) == False
    assert (1 == (x*y - x*y)) == False

    assert ((x*y + 5*x*y) != 6) == True
    assert ((x*y + 5*x*y) != 6*x*y) == False
    assert (6 != (x*y + 5*x*y)) == True
    assert (6*x*y != (x*y + 5*x*y)) == False

    assert ((x*y - x*y) != 0) == False
    assert (0 != (x*y - x*y)) == False

    assert ((x*y - x*y) != 1) == True
    assert (1 != (x*y - x*y)) == True

    Rt, t = ring("t", ZZ)
    R, x, y = ring("x,y", Rt.to_domain())

    assert (t**3*x/x == t**3) == True
    assert (t**3*x/x == t**4) == False

def test_PolyElement_copy():
    R, x, y, z = ring("x,y,z", ZZ)

    f = x*y + 3*z
    g = f.copy()

    assert f == g
    g[(1, 1, 1)] = 7
    assert f != g

def test_PolyElement_as_expr():
    R, x, y, z = ring("x,y,z", ZZ)
    f = 3*x**2*y - x*y*z + 7*z**3 + 1

    X, Y, Z = R.symbols
    g = 3*X**2*Y - X*Y*Z + 7*Z**3 + 1

    assert f != g
    assert f.as_expr() == g

    X, Y, Z = symbols("x,y,z")
    g = 3*X**2*Y - X*Y*Z + 7*Z**3 + 1

    assert f != g
    assert f.as_expr(X, Y, Z) == g

    raises(ValueError, lambda: f.as_expr(X))

def test_PolyElement_from_expr():
    x, y, z = symbols("x,y,z")
    R, X, Y, Z = ring((x, y, z), ZZ)

    f = R.from_expr(1)
    assert f == 1 and isinstance(f, R.dtype)

    f = R.from_expr(x)
    assert f == X and isinstance(f, R.dtype)

    f = R.from_expr(x*y*z)
    assert f == X*Y*Z and isinstance(f, R.dtype)

    f = R.from_expr(x*y*z + x*y + x)
    assert f == X*Y*Z + X*Y + X and isinstance(f, R.dtype)

    f = R.from_expr(x**3*y*z + x**2*y**7 + 1)
    assert f == X**3*Y*Z + X**2*Y**7 + 1 and isinstance(f, R.dtype)

    raises(ValueError, lambda: R.from_expr(1/x))
    raises(ValueError, lambda: R.from_expr(2**x))
    raises(ValueError, lambda: R.from_expr(7*x + sqrt(2)))

def test_PolyElement_coeff():
    R, x, y, z = ring("x,y,z", ZZ, lex)
    f = 3*x**2*y - x*y*z + 7*z**3 + 23

    assert f.coeff(1) == 23
    raises(ValueError, lambda: f.coeff(3))

    assert f.coeff(x) == 0
    assert f.coeff(y) == 0
    assert f.coeff(z) == 0

    assert f.coeff(x**2*y) == 3
    assert f.coeff(x*y*z) == -1
    assert f.coeff(z**3) == 7

    raises(ValueError, lambda: f.coeff(3*x**2*y))
    raises(ValueError, lambda: f.coeff(-x*y*z))
    raises(ValueError, lambda: f.coeff(7*z**3))

def test_PolyElement_LC():
    R, x, y = ring("x,y", QQ, lex)
    assert R(0).LC == QQ(0)
    assert (QQ(1,2)*x).LC == QQ(1, 2)
    assert (QQ(1,4)*x*y + QQ(1,2)*x).LC == QQ(1, 4)

def test_PolyElement_LM():
    R, x, y = ring("x,y", QQ, lex)
    assert R(0).LM == (0, 0)
    assert (QQ(1,2)*x).LM == (1, 0)
    assert (QQ(1,4)*x*y + QQ(1,2)*x).LM == (1, 1)

def test_PolyElement_LT():
    R, x, y = ring("x,y", QQ, lex)
    assert R(0).LT == ((0, 0), QQ(0))
    assert (QQ(1,2)*x).LT == ((1, 0), QQ(1, 2))
    assert (QQ(1,4)*x*y + QQ(1,2)*x).LT == ((1, 1), QQ(1, 4))

def test_PolyElement_leading_monom():
    R, x, y = ring("x,y", QQ, lex)
    assert R(0).leading_monom == 0
    assert (QQ(1,2)*x).leading_monom == x
    assert (QQ(1,4)*x*y + QQ(1,2)*x).leading_monom == x*y

def test_PolyElement_leading_term():
    R, x, y = ring("x,y", QQ, lex)
    assert R(0).leading_term == 0
    assert (QQ(1,2)*x).leading_term == QQ(1,2)*x
    assert (QQ(1,4)*x*y + QQ(1,2)*x).leading_term == QQ(1,4)*x*y

def test_PolyElement_terms():
    R, x,y,z = ring("x,y,z", QQ)
    terms = list((x**2/3 + y**3/4 + z**4/5).terms())
    assert set(terms) == set([((2,0,0), QQ(1,3)), ((0,3,0), QQ(1,4)), ((0,0,4), QQ(1,5))])

def test_PolyElement_monoms():
    R, x,y,z = ring("x,y,z", QQ)
    monoms = list((x**2/3 + y**3/4 + z**4/5).monoms())
    assert set(monoms) == set([(2,0,0), (0,3,0), (0,0,4)])

def test_PolyElement_coeffs():
    R, x,y,z = ring("x,y,z", QQ)
    coeffs = list((x**2/3 + y**3/4 + z**4/5).coeffs())
    assert set(coeffs) == set([QQ(1,3), QQ(1,4), QQ(1,5)])

def test_PolyElement___add__():
    Rt, t = ring("t", ZZ)
    Ruv, u,v = ring("u,v", ZZ)
    Rxyz, x,y,z = ring("x,y,z", Ruv.to_domain())

    assert dict(x + 3*y) == {(1, 0, 0): 1, (0, 1, 0): 3}

    assert dict(u + x) == dict(x + u) == {(1, 0, 0): 1, (0, 0, 0): u}
    assert dict(u + x*y) == dict(x*y + u) == {(1, 1, 0): 1, (0, 0, 0): u}
    assert dict(u + x*y + z) == dict(x*y + z + u) == {(1, 1, 0): 1, (0, 0, 1): 1, (0, 0, 0): u}

    assert dict(u*x + x) == dict(x + u*x) == {(1, 0, 0): u + 1}
    assert dict(u*x + x*y) == dict(x*y + u*x) == {(1, 1, 0): 1, (1, 0, 0): u}
    assert dict(u*x + x*y + z) == dict(x*y + z + u*x) == {(1, 1, 0): 1, (0, 0, 1): 1, (1, 0, 0): u}

    raises(TypeError, lambda: t + x)
    raises(TypeError, lambda: x + t)
    raises(TypeError, lambda: t + u)
    raises(TypeError, lambda: u + t)

    Fuv, u,v = field("u,v", ZZ)
    Rxyz, x,y,z = ring("x,y,z", Fuv.to_domain())

    assert dict(u + x) == dict(x + u) == {(1, 0, 0): 1, (0, 0, 0): u}

def test_PolyElement___sub__():
    Rt, t = ring("t", ZZ)
    Ruv, u,v = ring("u,v", ZZ)
    Rxyz, x,y,z = ring("x,y,z", Ruv.to_domain())

    assert dict(x - 3*y) == {(1, 0, 0): 1, (0, 1, 0): -3}

    assert dict(-u + x) == dict(x - u) == {(1, 0, 0): 1, (0, 0, 0): -u}
    assert dict(-u + x*y) == dict(x*y - u) == {(1, 1, 0): 1, (0, 0, 0): -u}
    assert dict(-u + x*y + z) == dict(x*y + z - u) == {(1, 1, 0): 1, (0, 0, 1): 1, (0, 0, 0): -u}

    assert dict(-u*x + x) == dict(x - u*x) == {(1, 0, 0): -u + 1}
    assert dict(-u*x + x*y) == dict(x*y - u*x) == {(1, 1, 0): 1, (1, 0, 0): -u}
    assert dict(-u*x + x*y + z) == dict(x*y + z - u*x) == {(1, 1, 0): 1, (0, 0, 1): 1, (1, 0, 0): -u}

    raises(TypeError, lambda: t - x)
    raises(TypeError, lambda: x - t)
    raises(TypeError, lambda: t - u)
    raises(TypeError, lambda: u - t)

    Fuv, u,v = field("u,v", ZZ)
    Rxyz, x,y,z = ring("x,y,z", Fuv.to_domain())

    assert dict(-u + x) == dict(x - u) == {(1, 0, 0): 1, (0, 0, 0): -u}

def test_PolyElement___mul__():
    Rt, t = ring("t", ZZ)
    Ruv, u,v = ring("u,v", ZZ)
    Rxyz, x,y,z = ring("x,y,z", Ruv.to_domain())

    assert dict(u*x) == dict(x*u) == {(1, 0, 0): u}

    assert dict(2*u*x + z) == dict(x*2*u + z) == {(1, 0, 0): 2*u, (0, 0, 1): 1}
    assert dict(u*2*x + z) == dict(2*x*u + z) == {(1, 0, 0): 2*u, (0, 0, 1): 1}
    assert dict(2*u*x + z) == dict(x*2*u + z) == {(1, 0, 0): 2*u, (0, 0, 1): 1}
    assert dict(u*x*2 + z) == dict(x*u*2 + z) == {(1, 0, 0): 2*u, (0, 0, 1): 1}

    assert dict(2*u*x*y + z) == dict(x*y*2*u + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}
    assert dict(u*2*x*y + z) == dict(2*x*y*u + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}
    assert dict(2*u*x*y + z) == dict(x*y*2*u + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}
    assert dict(u*x*y*2 + z) == dict(x*y*u*2 + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}

    assert dict(2*u*y*x + z) == dict(y*x*2*u + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}
    assert dict(u*2*y*x + z) == dict(2*y*x*u + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}
    assert dict(2*u*y*x + z) == dict(y*x*2*u + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}
    assert dict(u*y*x*2 + z) == dict(y*x*u*2 + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}

    assert dict(3*u*(x + y) + z) == dict((x + y)*3*u + z) == {(1, 0, 0): 3*u, (0, 1, 0): 3*u, (0, 0, 1): 1}

    raises(TypeError, lambda: t*x + z)
    raises(TypeError, lambda: x*t + z)
    raises(TypeError, lambda: t*u + z)
    raises(TypeError, lambda: u*t + z)

    Fuv, u,v = field("u,v", ZZ)
    Rxyz, x,y,z = ring("x,y,z", Fuv.to_domain())

    assert dict(u*x) == dict(x*u) == {(1, 0, 0): u}

def test_PolyElement___div__():
    R, x,y,z = ring("x,y,z", ZZ)

    assert (2*x**2 - 4)/2 == x**2 - 2
    assert (2*x**2 - 3)/2 == x**2

    assert (x**2 - 1)/x == (x**2 - 1).quo(x) == x
    assert (x**2 - x)/x == (x**2 - x).quo(x) == x - 1

    assert (x**2 - 1)/(2*x) == (x**2 - 1).quo(2*x)== 0
    assert (x**2 - x)/(x - 1) == (x**2 - x).quo(x - 1) == x

    R, x,y,z = ring("x,y,z", ZZ)
    assert len(list((x**2/3 + y**3/4 + z**4/5).terms())) == 0

    R, x,y,z = ring("x,y,z", QQ)
    assert len(list((x**2/3 + y**3/4 + z**4/5).terms())) == 3

    Rt, t = ring("t", ZZ)
    Ruv, u,v = ring("u,v", ZZ)
    Rxyz, x,y,z = ring("x,y,z", Ruv.to_domain())

    assert dict((u**2*x + u)/u) == {(1, 0, 0): u, (0, 0, 0): 1}
    raises(TypeError, lambda: u/(u**2*x + u))

    raises(TypeError, lambda: t/x)
    # TODO: raises(TypeError, lambda: x/t)
    raises(TypeError, lambda: t/u)
    raises(TypeError, lambda: u/t)

def test_PolyElement_pow():
    R, x = ring("x", ZZ, grlex)
    f = 2*x + 3

    assert f**0 == 1
    assert f**1 == f

    assert f**2 == f._pow_generic(2) == f._pow_multinomial(2) == 4*x**2 + 12*x + 9
    assert f**3 == f._pow_generic(3) == f._pow_multinomial(3) == 8*x**3 + 36*x**2 + 54*x + 27
    assert f**4 == f._pow_generic(4) == f._pow_multinomial(4) == 16*x**4 + 96*x**3 + 216*x**2 + 216*x + 81
    assert f**5 == f._pow_generic(5) == f._pow_multinomial(5) == 32*x**5 + 240*x**4 + 720*x**3 + 1080*x**2 + 810*x + 243

    R, x,y,z = ring("x,y,z", ZZ, grlex)
    f = x**3*y - 2*x*y**2 - 3*z + 1
    g = x**6*y**2 - 4*x**4*y**3 - 6*x**3*y*z + 2*x**3*y + 4*x**2*y**4 + 12*x*y**2*z - 4*x*y**2 + 9*z**2 - 6*z + 1

    assert f**2 == f._pow_generic(2) == f._pow_multinomial(2) == g

    raises(ValueError, lambda: f**-2)

def test_PolyElement_div():
    R, x = ring("x", ZZ, grlex)

    f = x**3 - 12*x**2 - 42
    g = x - 3

    q = x**2 - 9*x - 27
    r = -123

    assert f.div([g]) == ([q], r)

    R, x = ring("x", ZZ, grlex)
    f = x**2 + 2*x + 2
    assert f.div([R(1)]) == ([f], 0)

    R, x = ring("x", QQ, grlex)
    f = x**2 + 2*x + 2
    assert f.div([R(2)]) == ([QQ(1,2)*x**2 + x + 1], 0)

    R, x,y = ring("x,y", ZZ, grlex)
    f = 4*x**2*y - 2*x*y + 4*x - 2*y + 8

    assert f.div([R(2)]) == ([2*x**2*y - x*y + 2*x - y + 4], 0)
    assert f.div([2*y]) == ([2*x**2 - x - 1], 4*x + 8)

    f = x - 1
    g = y - 1

    assert f.div([g]) == ([0], f)

    f = x*y**2 + 1
    G = [x*y + 1, y + 1]

    Q = [y, -1]
    r = 2

    assert f.div(G) == (Q, r)

    f = x**2*y + x*y**2 + y**2
    G = [x*y - 1, y**2 - 1]

    Q = [x + y, 1]
    r = x + y + 1

    assert f.div(G) == (Q, r)

    G = [y**2 - 1, x*y - 1]

    Q = [x + 1, x]
    r = 2*x + 1

    assert f.div(G) == (Q, r)

def test_PolyElement_rem():
    R, x = ring("x", ZZ, grlex)

    f = x**3 - 12*x**2 - 42
    g = x - 3
    r = -123

    assert f.rem([g]) == f.div([g])[1] == r

    R, x,y = ring("x,y", ZZ, grlex)

    f = 4*x**2*y - 2*x*y + 4*x - 2*y + 8

    assert f.rem([R(2)]) == f.div([R(2)])[1] == 0
    assert f.rem([2*y]) == f.div([2*y])[1] == 4*x + 8

    f = x - 1
    g = y - 1

    assert f.rem([g]) == f.div([g])[1] == f

    f = x*y**2 + 1
    G = [x*y + 1, y + 1]
    r = 2

    assert f.rem(G) == f.div(G)[1] == r

    f = x**2*y + x*y**2 + y**2
    G = [x*y - 1, y**2 - 1]
    r = x + y + 1

    assert f.rem(G) == f.div(G)[1] == r

    G = [y**2 - 1, x*y - 1]
    r = 2*x + 1

    assert f.rem(G) == f.div(G)[1] == r

def test_PolyElement_deflate():
    R, x = ring("x", ZZ)

    assert (2*x**2).deflate(x**4 + 4*x**2 + 1) == ((2,), [2*x, x**2 + 4*x + 1])

    R, x,y = ring("x,y", ZZ)

    assert R(0).deflate(R(0)) == ((1, 1), [0, 0])
    assert R(1).deflate(R(0)) == ((1, 1), [1, 0])
    assert R(1).deflate(R(2)) == ((1, 1), [1, 2])
    assert R(1).deflate(2*y) == ((1, 1), [1, 2*y])
    assert (2*y).deflate(2*y) == ((1, 1), [2*y, 2*y])
    assert R(2).deflate(2*y**2) == ((1, 2), [2, 2*y])
    assert (2*y**2).deflate(2*y**2) == ((1, 2), [2*y, 2*y])

    f = x**4*y**2 + x**2*y + 1
    g = x**2*y**3 + x**2*y + 1

    assert f.deflate(g) == ((2, 1), [x**2*y**2 + x*y + 1, x*y**3 + x*y + 1])

def test_PolyElement_diff():
    R, X = xring("x:11", QQ)

    f = QQ(288,5)*X[0]**8*X[1]**6*X[4]**3*X[10]**2 + 8*X[0]**2*X[2]**3*X[4]**3 +2*X[0]**2 - 2*X[1]**2

    assert f.diff(X[0]) == QQ(2304,5)*X[0]**7*X[1]**6*X[4]**3*X[10]**2 + 16*X[0]*X[2]**3*X[4]**3 + 4*X[0]
    assert f.diff(X[4]) == QQ(864,5)*X[0]**8*X[1]**6*X[4]**2*X[10]**2 + 24*X[0]**2*X[2]**3*X[4]**2
    assert f.diff(X[10]) == QQ(576,5)*X[0]**8*X[1]**6*X[4]**3*X[10]

def test_PolyElement_clear_denoms():
    R, x,y = ring("x,y", QQ)

    assert R(1).clear_denoms() == (ZZ(1), 1)
    assert R(7).clear_denoms() == (ZZ(1), 7)

    assert R(QQ(7,3)).clear_denoms() == (3, 7)
    assert R(QQ(7,3)).clear_denoms() == (3, 7)

    assert (3*x**2 + x).clear_denoms() == (1, 3*x**2 + x)
    assert (x**2 + QQ(1,2)*x).clear_denoms() == (2, 2*x**2 + x)

    rQQ, x,t = ring("x,t", QQ, lex)
    rZZ, X,T = ring("x,t", ZZ, lex)

    F = [x - QQ(17824537287975195925064602467992950991718052713078834557692023531499318507213727406844943097,413954288007559433755329699713866804710749652268151059918115348815925474842910720000)*t**7
           - QQ(4882321164854282623427463828745855894130208215961904469205260756604820743234704900167747753,12936071500236232304854053116058337647210926633379720622441104650497671088840960000)*t**6
           - QQ(36398103304520066098365558157422127347455927422509913596393052633155821154626830576085097433,25872143000472464609708106232116675294421853266759441244882209300995342177681920000)*t**5
           - QQ(168108082231614049052707339295479262031324376786405372698857619250210703675982492356828810819,58212321751063045371843239022262519412449169850208742800984970927239519899784320000)*t**4
           - QQ(5694176899498574510667890423110567593477487855183144378347226247962949388653159751849449037,1617008937529529038106756639507292205901365829172465077805138081312208886105120000)*t**3
           - QQ(154482622347268833757819824809033388503591365487934245386958884099214649755244381307907779,60637835157357338929003373981523457721301218593967440417692678049207833228942000)*t**2
           - QQ(2452813096069528207645703151222478123259511586701148682951852876484544822947007791153163,2425513406294293557160134959260938308852048743758697616707707121968313329157680)*t
           - QQ(34305265428126440542854669008203683099323146152358231964773310260498715579162112959703,202126117191191129763344579938411525737670728646558134725642260164026110763140),
         t**8 + QQ(693749860237914515552,67859264524169150569)*t**7
              + QQ(27761407182086143225024,610733380717522355121)*t**6
              + QQ(7785127652157884044288,67859264524169150569)*t**5
              + QQ(36567075214771261409792,203577793572507451707)*t**4
              + QQ(36336335165196147384320,203577793572507451707)*t**3
              + QQ(7452455676042754048000,67859264524169150569)*t**2
              + QQ(2593331082514399232000,67859264524169150569)*t
              + QQ(390399197427343360000,67859264524169150569)]

    G = [3725588592068034903797967297424801242396746870413359539263038139343329273586196480000*X -
         160420835591776763325581422211936558925462474417709511019228211783493866564923546661604487873*T**7 -
         1406108495478033395547109582678806497509499966197028487131115097902188374051595011248311352864*T**6 -
         5241326875850889518164640374668786338033653548841427557880599579174438246266263602956254030352*T**5 -
         10758917262823299139373269714910672770004760114329943852726887632013485035262879510837043892416*T**4 -
         13119383576444715672578819534846747735372132018341964647712009275306635391456880068261130581248*T**3 -
         9491412317016197146080450036267011389660653495578680036574753839055748080962214787557853941760*T**2 -
         3767520915562795326943800040277726397326609797172964377014046018280260848046603967211258368000*T -
         632314652371226552085897259159210286886724229880266931574701654721512325555116066073245696000,
         610733380717522355121*T**8 +
         6243748742141230639968*T**7 +
         27761407182086143225024*T**6 +
         70066148869420956398592*T**5 +
         109701225644313784229376*T**4 +
         109009005495588442152960*T**3 +
         67072101084384786432000*T**2 +
         23339979742629593088000*T +
         3513592776846090240000]

    assert [ f.clear_denoms()[1].set_ring(rZZ) for f in F ] == G

def test_PolyElement_cofactors():
    R, x, y = ring("x,y", ZZ)

    f, g = R(0), R(0)
    assert f.cofactors(g) == (0, 0, 0)

    f, g = R(2), R(0)
    assert f.cofactors(g) == (2, 1, 0)

    f, g = R(-2), R(0)
    assert f.cofactors(g) == (2, -1, 0)

    f, g = R(0), R(-2)
    assert f.cofactors(g) == (2, 0, -1)

    f, g = R(0), 2*x + 4
    assert f.cofactors(g) == (2*x + 4, 0, 1)

    f, g = 2*x + 4, R(0)
    assert f.cofactors(g) == (2*x + 4, 1, 0)

    f, g = R(2), R(2)
    assert f.cofactors(g) == (2, 1, 1)

    f, g = R(-2), R(2)
    assert f.cofactors(g) == (2, -1, 1)

    f, g = R(2), R(-2)
    assert f.cofactors(g) == (2, 1, -1)

    f, g = R(-2), R(-2)
    assert f.cofactors(g) == (2, -1, -1)

    f, g = x**2 + 2*x + 1, R(1)
    assert f.cofactors(g) == (1, x**2 + 2*x + 1, 1)

    f, g = x**2 + 2*x + 1, R(2)
    assert f.cofactors(g) == (1, x**2 + 2*x + 1, 2)

    f, g = 2*x**2 + 4*x + 2, R(2)
    assert f.cofactors(g) == (2, x**2 + 2*x + 1, 1)

    f, g = R(2), 2*x**2 + 4*x + 2
    assert f.cofactors(g) == (2, 1, x**2 + 2*x + 1)

    f, g = 2*x**2 + 4*x + 2, x + 1
    assert f.cofactors(g) == (x + 1, 2*x + 2, 1)

    f, g = x + 1, 2*x**2 + 4*x + 2
    assert f.cofactors(g) == (x + 1, 1, 2*x + 2)

    R, x, y, z, t = ring("x,y,z,t", ZZ)

    f, g = t**2 + 2*t + 1, 2*t + 2
    assert f.cofactors(g) == (t + 1, t + 1, 2)

    f, g = z**2*t**2 + 2*z**2*t + z**2 + z*t + z, t**2 + 2*t + 1
    h, cff, cfg = t + 1, z**2*t + z**2 + z, t + 1

    assert f.cofactors(g) == (h, cff, cfg)
    assert g.cofactors(f) == (h, cfg, cff)

    R, x, y = ring("x,y", QQ)

    f = QQ(1,2)*x**2 + x + QQ(1,2)
    g = QQ(1,2)*x + QQ(1,2)

    h = x + 1

    assert f.cofactors(g) == (h, g, QQ(1,2))
    assert g.cofactors(f) == (h, QQ(1,2), g)

    R, x, y = ring("x,y", RR)

    f = 2.1*x*y**2 - 2.1*x*y + 2.1*x
    g = 2.1*x**3
    h = 1.0*x

    assert f.cofactors(g) == (h, f/h, g/h)
    assert g.cofactors(f) == (h, g/h, f/h)

def test_PolyElement_gcd():
    R, x, y = ring("x,y", QQ)

    f = QQ(1,2)*x**2 + x + QQ(1,2)
    g = QQ(1,2)*x + QQ(1,2)

    assert f.gcd(g) == x + 1

def test_PolyElement_cancel():
    R, x, y = ring("x,y", ZZ)

    f = 2*x**3 + 4*x**2 + 2*x
    g = 3*x**2 + 3*x
    F = 2*x + 2
    G = 3

    assert f.cancel(g) == (F, G)

    assert (-f).cancel(g) == (-F, G)
    assert f.cancel(-g) == (-F, G)

    R, x, y = ring("x,y", QQ)

    f = QQ(1,2)*x**3 + x**2 + QQ(1,2)*x
    g = QQ(1,3)*x**2 + QQ(1,3)*x
    F = 3*x + 3
    G = 2

    assert f.cancel(g) == (F, G)

    assert (-f).cancel(g) == (-F, G)
    assert f.cancel(-g) == (-F, G)

def test_PolyElement_max_norm():
    R, x, y = ring("x,y", ZZ)

    assert R(0).max_norm() == 0
    assert R(1).max_norm() == 1

    assert (x**3 + 4*x**2 + 2*x + 3).max_norm() == 4

def test_PolyElement_evaluate():
    R, x = ring("x", ZZ)
    r = (x**3 + 4*x**2 + 2*x + 3).evaluate(x, 0)
    assert r == 3 and not isinstance(r, PolyElement)

    R, x, y, z = ring("x,y,z", ZZ)
    r = (x**3 + 4*x**2 + 2*x + 3).evaluate(x, 0)
    assert r == 3 and isinstance(r, PolyElement) and r.ring == R[1:]
    r = (x**3 + 4*x**2 + 2*x + 3).evaluate([(x, 0), (y, 0)])
    assert r == 3 and isinstance(r, PolyElement) and r.ring == R[2:]

def test_PolyElement_subs():
    R, x = ring("x", ZZ)
    r = (x**3 + 4*x**2 + 2*x + 3).subs(x, 0)
    assert r == 3 and isinstance(r, PolyElement) and r.ring == R

    R, x, y, z = ring("x,y,z", ZZ)
    r = (x**3 + 4*x**2 + 2*x + 3).subs(x, 0)
    assert r == 3 and isinstance(r, PolyElement) and r.ring == R
    r = (x**3 + 4*x**2 + 2*x + 3).subs([(x, 0), (y, 0)])
    assert r == 3 and isinstance(r, PolyElement) and r.ring == R

def test_PolyElement_compose():
    R, x = ring("x", ZZ)
    r = (x**3 + 4*x**2 + 2*x + 3).compose(x, 0)
    assert r == 3 and isinstance(r, PolyElement) and r.ring == R

    R, x, y, z = ring("x,y,z", ZZ)
    r = (x**3 + 4*x**2 + 2*x + 3).compose(x, 0)
    assert r == 3 and isinstance(r, PolyElement) and r.ring == R
    r = (x**3 + 4*x**2 + 2*x + 3).compose([(x, 0), (y, 0)])
    assert r == 3 and isinstance(r, PolyElement) and r.ring == R

    r = (x**3 + 4*x**2 + 2*x*y*z + 3).compose(x, y*z**2 - 1)
    q = (y*z**2 - 1)**3 + 4*(y*z**2 - 1)**2 + 2*(y*z**2 - 1)*y*z + 3
    assert r == q and isinstance(r, PolyElement) and r.ring == R

def test_PolyElement_is_():
    R, x,y = ring("x,y", QQ)

    assert (x - x).is_generator == False
    assert (x - x).is_ground == True
    assert (x - x).is_monomial == True
    assert (x - x).is_term == True

    assert (x - x + 1).is_generator == False
    assert (x - x + 1).is_ground == True
    assert (x - x + 1).is_monomial == True
    assert (x - x + 1).is_term == True

    assert x.is_generator == True
    assert x.is_ground == False
    assert x.is_monomial == True
    assert x.is_term == True

    assert (x*y).is_generator == False
    assert (x*y).is_ground == False
    assert (x*y).is_monomial == True
    assert (x*y).is_term == True

    assert (3*x).is_generator == False
    assert (3*x).is_ground == False
    assert (3*x).is_monomial == False
    assert (3*x).is_term == True

    assert (3*x + 1).is_generator == False
    assert (3*x + 1).is_ground == False
    assert (3*x + 1).is_monomial == False
    assert (3*x + 1).is_term == False
