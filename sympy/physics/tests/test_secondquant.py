from sympy.physics.secondquant import (
        Dagger, Bd, VarBosonicBasis, Bra, B, Ket, FixedBosonicBasis,
        matrix_rep, apply_operators, InnerProduct, commutator, KroneckerDelta,
        FockState, AnnihilateBoson, CreateBoson, BosonicOperator
        )

from sympy import (
    symbols, Symbol, sympify,
    sqrt, Rational, Sum, I, simplify,
    expand
)

def test_dagger():
    i, j, n, m = symbols('i j n m')
    assert Dagger(1) == 1
    assert Dagger(1.0) == 1.0
    assert Dagger(2*I) == -2*I
    assert Dagger(Rational(1,2)*I/3.0) == -Rational(1,2)*I/3.0
    assert Dagger(Ket([n])) == Bra([n])
    assert Dagger(B(0)) == Bd(0)
    assert Dagger(Bd(0)) == B(0)
    assert Dagger(B(n)) == Bd(n)
    assert Dagger(Bd(n)) == B(n)
    assert Dagger(B(0)+B(1)) == Bd(0) + Bd(1)
    assert Dagger(n*m) == Dagger(n)*Dagger(m) # n, m commute
    assert Dagger(B(n)*B(m)) == Bd(m)*Bd(n)
    assert Dagger(B(n)**10) == Dagger(B(n))**10

def test_operator():
    i, j = symbols('i j')
    o = BosonicOperator(i)
    assert o.state == i
    assert o.is_symbolic
    o = BosonicOperator(1)
    assert o.state == 1
    assert not o.is_symbolic

def test_create():
    i, j, n, m = symbols('i j n m')
    o = Bd(i)
    assert isinstance(o, CreateBoson)
    o = o.subs(i, j)
    assert o.atoms(Symbol) == set([j])
    o = Bd(0)
    assert o.apply_operator(Ket([n])) == sqrt(n+1)*Ket([n+1])
    o = Bd(n)
    assert o.apply_operator(Ket([n])) == o*Ket([n])

def test_annihilate():
    i, j, n, m = symbols('i j n m')
    o = B(i)
    assert isinstance(o, AnnihilateBoson)
    o = o.subs(i, j)
    assert o.atoms(Symbol) == set([j])
    o = B(0)
    assert o.apply_operator(Ket([n])) == sqrt(n)*Ket([n-1])
    o = B(n)
    assert o.apply_operator(Ket([n])) == o*Ket([n])

def test_basic_state():
    i, j, n, m = symbols('i j n m')
    s = FockState([0,1,2,3,4])
    assert len(s) == 5
    assert s.args[0] == tuple(range(5))
    assert s.up(0) == FockState([1,1,2,3,4])
    assert s.down(4) == FockState([0,1,2,3,3])
    for i in range(5):
        assert s.up(i).down(i) == s
    assert s.down(0) == 0
    for i in range(5):
        assert s[i] == i
    s = FockState([n,m])
    assert s.down(0) == FockState([n-1,m])
    assert s.up(0) == FockState([n+1,m])

def test_kronecker_delta():
    i, j, k = symbols('i j k')
    assert KroneckerDelta(i, i) == 1
    assert KroneckerDelta(i, i + 1) == 0
    assert KroneckerDelta(0, 0) == 1
    assert KroneckerDelta(0, 1) == 0
    # assert KroneckerDelta(i, i + k) == KroneckerDelta(1, k)
    assert KroneckerDelta(i + k, i + k) == 1
    assert KroneckerDelta(i + k, i + 1 + k) == 0
    assert KroneckerDelta(i, j).subs(dict(i=1, j=0)) == 0

# def Xtest_move1():
#     i, j = symbols('i j')
#     o = A(i)*C(j)
#     # This almost works, but has a minus sign wrong
#     assert move(o, 0, 1) == KroneckerDelta(i, j) + C(j)*A(i)
#
# def Xtest_move2():
#     i, j = symbols('i j')
#     o = C(j)*A(i)
#     # This almost works, but has a minus sign wrong
#     assert move(o, 0, 1) == -KroneckerDelta(i, j) + A(i)*C(j)

def test_basic_apply():
    n = symbols("n")
    e = B(0)*Ket([n])
    assert apply_operators(e) == sqrt(n)*Ket([n-1])
    e = Bd(0)*Ket([n])
    assert apply_operators(e) == sqrt(n+1)*Ket([n+1])

def test_commutation():
    n, m = symbols("n m")
    c = commutator(B(0), Bd(0))
    e = simplify(apply_operators(c*Ket([n])))
    assert e == Ket([n])
    c = commutator(B(0), B(1))
    e = simplify(apply_operators(c*Ket([n,m])))
    assert e == 0

def test_complex_apply():
    n, m = symbols("n m")
    o = Bd(0)*B(0)*Bd(1)*B(0)
    e = apply_operators(o*Ket([n,m]))
    answer = sqrt(n)*sqrt(m+1)*(-1+n)*Ket([-1+n,1+m])
    assert expand(e) == expand(answer)

def test_number_operator():
    n = symbols("n")
    o = Bd(0)*B(0)
    e = apply_operators(o*Ket([n]))
    assert e == n*Ket([n])

def test_inner_product():
    i, j, k, l = symbols('i j k l')
    s1 = Bra([0])
    s2 = Ket([1])
    assert InnerProduct(s1, Dagger(s1)) == 1
    assert InnerProduct(s1, s2) == 0
    s1 = Bra([i, j])
    s2 = Ket([k, l])
    r = InnerProduct(s1, s2)
    assert r == KroneckerDelta(i, k)*KroneckerDelta(j, l)

def test_symbolic_matrix_elements():
    n, m = symbols('n m')
    s1 = Bra([n])
    s2 = Ket([m])
    o = B(0)
    e = apply_operators(s1*o*s2)
    assert e == sqrt(m)*KroneckerDelta(n, m-1)

def test_matrix_elements():
    b = VarBosonicBasis(5)
    o = B(0)
    m = matrix_rep(o, b)
    for i in range(4):
        assert m[i, i+1] == sqrt(i+1)
    o = Bd(0)
    m = matrix_rep(o, b)
    for i in range(4):
        assert m[i+1, i] == sqrt(i+1)

def test_sho():
    n, m = symbols('n m')
    h_n = Bd(n)*B(n)*(n + Rational(1, 2))
    H = Sum(h_n, (n, 0, 5))
    o = H.doit()
    b = FixedBosonicBasis(2, 6)
    m = matrix_rep(o, b)
    # We need to double check these energy values to make sure that they
    # are correct and have the proper degeneracies!
    diag = [1, 2, 3, 3, 4, 5, 4, 5, 6, 7, 5, 6, 7, 8, 9, 6, 7, 8, 9, 10, 11]
    for i in range(len(diag)):
        assert diag[i] == m[i, i]
