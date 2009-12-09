from sympy.physics.secondquant import (
    Dagger, Bd, VarBosonicBasis, BBra, B, BKet, FixedBosonicBasis,
    matrix_rep, apply_operators, InnerProduct, Commutator, KroneckerDelta,
    FockState, AnnihilateBoson, CreateBoson, BosonicOperator,
    F, Fd, FKet, FBra, BosonState, CreateFermion, AnnihilateFermion,
    evaluate_deltas, AntiSymmetricTensor, contraction, NO, wicks,
    PermutationOperator, simplify_index_permutations, SymTuple,
        )

from sympy import (
    symbols, Symbol, sympify,
    sqrt, Rational, Sum, I, simplify,
    expand, Function
)


def test_SymTuple():
    t = (1,2,3,4)
    st =  SymTuple(t)
    assert set(t) == set(st)
    assert len(t) == len(st)
    assert set(t[:2]) == set(st[:2])
    assert isinstance(st[:], SymTuple)
    assert st == SymTuple((1,2,3,4))
    p,q,r,s = symbols('pqrs')
    t2=(p,q,r,s)
    st2 = SymTuple(t2)
    assert st2.atoms() == set(t2)

def test_PermutationOperator():
    p,q,r,s = symbols('pqrs')
    f,g,h,i = map(Function, 'fghi')
    P = PermutationOperator
    assert (P(p,q).get_permuted(f(p)*g(q)) ==
            -f(q)*g(p))
    expr = (f(p)*g(q)*h(r)*i(s)
        - f(q)*g(p)*h(r)*i(s)
        - f(p)*g(q)*h(s)*i(r)
        + f(q)*g(p)*h(s)*i(r))
    perms = [P(p,q),P(r,s)]
    assert (simplify_index_permutations(expr,perms) ==
        P(p,q)*P(r,s)*f(p)*g(q)*h(r)*i(s))



def test_dagger():
    i, j, n, m = symbols('i j n m')
    assert Dagger(1) == 1
    assert Dagger(1.0) == 1.0
    assert Dagger(2*I) == -2*I
    assert Dagger(Rational(1,2)*I/3.0) == -Rational(1,2)*I/3.0
    assert Dagger(BKet([n])) == BBra([n])
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
    assert o.apply_operator(BKet([n])) == sqrt(n+1)*BKet([n+1])
    o = Bd(n)
    assert o.apply_operator(BKet([n])) == o*BKet([n])

def test_annihilate():
    i, j, n, m = symbols('i j n m')
    o = B(i)
    assert isinstance(o, AnnihilateBoson)
    o = o.subs(i, j)
    assert o.atoms(Symbol) == set([j])
    o = B(0)
    assert o.apply_operator(BKet([n])) == sqrt(n)*BKet([n-1])
    o = B(n)
    assert o.apply_operator(BKet([n])) == o*BKet([n])

def test_basic_state():
    i, j, n, m = symbols('i j n m')
    s = BosonState([0,1,2,3,4])
    assert len(s) == 5
    assert s.args[0] == tuple(range(5))
    assert s.up(0) == BosonState([1,1,2,3,4])
    assert s.down(4) == BosonState([0,1,2,3,3])
    for i in range(5):
        assert s.up(i).down(i) == s
    assert s.down(0) == 0
    for i in range(5):
        assert s[i] == i
    s = BosonState([n,m])
    assert s.down(0) == BosonState([n-1,m])
    assert s.up(0) == BosonState([n+1,m])

def test_kronecker_delta():
    i, j, k = symbols('i j k')
    D = KroneckerDelta
    assert D(i, i) == 1
    assert D(i, i + 1) == 0
    assert D(0, 0) == 1
    assert D(0, 1) == 0
    # assert D(i, i + k) == D(0, k)
    assert D(i + k, i + k) == 1
    assert D(i + k, i + 1 + k) == 0
    assert D(i, j).subs(dict(i=1, j=0)) == 0
    assert D(i, j).subs(dict(i=3, j=3)) == 1


    i,j,k,l = symbols('ijkl',below_fermi=True,dummy=True)
    a,b,c,d = symbols('abcd',above_fermi=True, dummy=True)
    p,q,r,s = symbols('pqrs',dumy=True)

    assert D(i,a) == 0

    assert D(i,j).is_above_fermi == False
    assert D(a,b).is_above_fermi == True
    assert D(p,q).is_above_fermi == True
    assert D(i,q).is_above_fermi == False
    assert D(a,q).is_above_fermi == True

    assert D(i,j).is_below_fermi == True
    assert D(a,b).is_below_fermi == False
    assert D(p,q).is_below_fermi == True
    assert D(p,j).is_below_fermi == True
    assert D(q,b).is_below_fermi == False

    assert not D(i,q).indices_contain_equal_information
    assert not D(a,q).indices_contain_equal_information
    assert D(p,q).indices_contain_equal_information
    assert D(a,b).indices_contain_equal_information
    assert D(i,j).indices_contain_equal_information

    assert D(q,b).preferred_index == b
    assert D(q,b).killable_index == q
    assert D(q,i).preferred_index == i
    assert D(q,i).killable_index == q
    assert D(q,p).preferred_index == p
    assert D(q,p).killable_index == q


    EV = evaluate_deltas
    assert EV(D(a,q)*F(q)) == F(a)
    assert EV(D(i,q)*F(q)) == F(i)
    assert EV(D(a,q)*F(a)) == D(a,q)*F(a)
    assert EV(D(i,q)*F(i)) == D(i,q)*F(i)
    assert EV(D(a,b)*F(a)) == F(b)
    assert EV(D(a,b)*F(b)) == F(a)
    assert EV(D(i,j)*F(i)) == F(j)
    assert EV(D(i,j)*F(j)) == F(i)
    assert EV(D(p,q)*F(q)) == F(p)
    assert EV(D(p,q)*F(p)) == F(q)
    assert EV(D(p,j)*D(p,i)*F(i)) == F(j)
    assert EV(D(p,j)*D(p,i)*F(j)) == F(i)
    assert EV(D(p,q)*D(p,i))*F(i) == D(q,i)*F(i)


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
    e = B(0)*BKet([n])
    assert apply_operators(e) == sqrt(n)*BKet([n-1])
    e = Bd(0)*BKet([n])
    assert apply_operators(e) == sqrt(n+1)*BKet([n+1])


def test_complex_apply():
    n, m = symbols("n m")
    o = Bd(0)*B(0)*Bd(1)*B(0)
    e = apply_operators(o*BKet([n,m]))
    answer = sqrt(n)*sqrt(m+1)*(-1+n)*BKet([-1+n,1+m])
    assert expand(e) == expand(answer)

def test_number_operator():
    n = symbols("n")
    o = Bd(0)*B(0)
    e = apply_operators(o*BKet([n]))
    assert e == n*BKet([n])

def test_inner_product():
    i, j, k, l = symbols('i j k l')
    s1 = BBra([0])
    s2 = BKet([1])
    assert InnerProduct(s1, Dagger(s1)) == 1
    assert InnerProduct(s1, s2) == 0
    s1 = BBra([i, j])
    s2 = BKet([k, l])
    r = InnerProduct(s1, s2)
    assert r == KroneckerDelta(i, k)*KroneckerDelta(j, l)

def test_symbolic_matrix_elements():
    n, m = symbols('n m')
    s1 = BBra([n])
    s2 = BKet([m])
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
    o = H.doit(deep = False)
    b = FixedBosonicBasis(2, 6)
    m = matrix_rep(o, b)
    # We need to double check these energy values to make sure that they
    # are correct and have the proper degeneracies!
    diag = [1, 2, 3, 3, 4, 5, 4, 5, 6, 7, 5, 6, 7, 8, 9, 6, 7, 8, 9, 10, 11]
    for i in range(len(diag)):
        assert diag[i] == m[i, i]

def test_commutation():
    n, m = symbols("n m", above_fermi=True)
    c = Commutator(B(0), Bd(0))
    assert c == 1
    c = Commutator(Bd(0), B(0))
    assert c == -1
    c = Commutator(B(n), Bd(0))
    assert c == KroneckerDelta(n,0)
    c = Commutator(B(0), Bd(0))
    e = simplify(apply_operators(c*BKet([n])))
    assert e == BKet([n])
    c = Commutator(B(0), B(1))
    e = simplify(apply_operators(c*BKet([n,m])))
    assert e == 0


    c = Commutator(F(m), Fd(m))
    assert c == +1 - 2*NO(Fd(m)*F(m))
    c = Commutator(Fd(m), F(m))
    assert c == -1 + 2*NO(Fd(m)*F(m))

    C = Commutator
    X,Y,Z = symbols('XYZ',commutative=False)
    assert C(C(X,Y),Z) != 0
    assert C(C(X,Z),Y) != 0
    assert C(Y,C(X,Z)) != 0
    # assert (C(C(Y,Z),X).eval_nested() + C(C(Z,X),Y).eval_nested() + C(C(X,Y),Z).eval_nested()) == 0
    # assert (C(X,C(Y,Z)).eval_nested() + C(Y,C(Z,X)).eval_nested() + C(Z,C(X,Y)).eval_nested()) == 0

    i,j,k,l = symbols('ijkl',below_fermi=True)
    a,b,c,d = symbols('abcd',above_fermi=True)
    p,q,r,s = symbols('pqrs')
    D=KroneckerDelta

    assert C(Fd(a),F(i)) == -2*NO(F(i)*Fd(a))
    assert C(Fd(j),NO(Fd(a)*F(i))).doit() == -D(j,i)*Fd(a)
    assert C(Fd(a)*F(i),Fd(b)*F(j)).doit() == 0

def test_create_f():
    i, j, n, m = symbols('i j n m')
    o = Fd(i)
    assert isinstance(o, CreateFermion)
    o = o.subs(i, j)
    assert o.atoms(Symbol) == set([j])
    o = Fd(1)
    assert o.apply_operator(FKet([n])) == FKet([1,n])
    assert o.apply_operator(FKet([n])) ==-FKet([n,1])
    o = Fd(n)
    assert o.apply_operator(FKet([])) == FKet([n])

    vacuum = FKet([],fermi_level=4)
    assert vacuum == FKet([],fermi_level=4)

    i,j,k,l = symbols('ijkl',below_fermi=True)
    a,b,c,d = symbols('abcd',above_fermi=True)
    p,q,r,s = symbols('pqrs')

    assert Fd(i).apply_operator(FKet([i,j,k],4)) == FKet([j,k],4)
    assert Fd(a).apply_operator(FKet([i,b,k],4)) == FKet([a,i,b,k],4)



def test_annihilate_f():
    i, j, n, m = symbols('i j n m')
    o = F(i)
    assert isinstance(o, AnnihilateFermion)
    o = o.subs(i, j)
    assert o.atoms(Symbol) == set([j])
    o = F(1)
    assert o.apply_operator(FKet([1,n])) == FKet([n])
    assert o.apply_operator(FKet([n,1])) ==-FKet([n])
    o = F(n)
    assert o.apply_operator(FKet([n])) == FKet([])

    i,j,k,l = symbols('ijkl',below_fermi=True)
    a,b,c,d = symbols('abcd',above_fermi=True)
    p,q,r,s = symbols('pqrs')
    assert F(i).apply_operator(FKet([i,j,k],4)) == 0
    assert F(a).apply_operator(FKet([i,b,k],4)) == 0
    assert F(l).apply_operator(FKet([i,j,k],3)) == 0
    assert F(l).apply_operator(FKet([i,j,k],4)) == FKet([l,i,j,k],4)

def test_create_b():
    i, j, n, m = symbols('i j n m')
    o = Bd(i)
    assert isinstance(o, CreateBoson)
    o = o.subs(i, j)
    assert o.atoms(Symbol) == set([j])
    o = Bd(0)
    assert o.apply_operator(BKet([n])) == sqrt(n+1)*BKet([n+1])
    o = Bd(n)
    assert o.apply_operator(BKet([n])) == o*BKet([n])

def test_annihilate_b():
    i, j, n, m = symbols('i j n m')
    o = B(i)
    assert isinstance(o, AnnihilateBoson)
    o = o.subs(i, j)
    assert o.atoms(Symbol) == set([j])
    o = B(0)

def test_wicks():
    p,q,r,s = symbols('pqrs',above_fermi=True)

    # Testing for particles only

    str = F(p)*Fd(q)
    assert wicks(str) == NO(F(p)*Fd(q)) + KroneckerDelta(p,q)
    str = Fd(p)*F(q)
    assert wicks(str) == NO(Fd(p)*F(q))


    str = F(p)*Fd(q)*F(r)*Fd(s)
    nstr= wicks(str)
    fasit = NO(
    KroneckerDelta(p, q)*KroneckerDelta(r, s)
    + KroneckerDelta(p, q)*AnnihilateFermion(r)*CreateFermion(s)
    + KroneckerDelta(r, s)*AnnihilateFermion(p)*CreateFermion(q)
    - KroneckerDelta(p, s)*AnnihilateFermion(r)*CreateFermion(q)
    - AnnihilateFermion(p)*AnnihilateFermion(r)*CreateFermion(q)*CreateFermion(s))
    assert nstr == fasit

    assert (p*q*nstr).expand() == wicks(p*q*str)
    assert (nstr*p*q*2).expand() == wicks(str*p*q*2)


    # Testing CC equations particles and holes
    i,j,k,l = symbols('ijkl',below_fermi=True,dummy=True)
    a,b,c,d = symbols('abcd',above_fermi=True,dummy=True)
    p,q,r,s = symbols('pqrs',dummy=True)

    assert (wicks(F(a)*NO(F(i)*F(j))*Fd(b)) ==
            NO(F(a)*F(i)*F(j)*Fd(b)) +
            KroneckerDelta(a,b)*NO(F(i)*F(j)))
    assert (wicks(F(a)*NO(F(i)*F(j)*F(k))*Fd(b)) ==
            NO(F(a)*F(i)*F(j)*F(k)*Fd(b)) -
            KroneckerDelta(a,b)*NO(F(i)*F(j)*F(k)))


    expr = wicks(Fd(i)*NO(Fd(j)*F(k))*F(l))
    assert (expr ==
           -KroneckerDelta(i,k)*NO(Fd(j)*F(l)) -
            KroneckerDelta(j,l)*NO(Fd(i)*F(k)) -
            KroneckerDelta(i,k)*KroneckerDelta(j,l)+
            KroneckerDelta(i,l)*NO(Fd(j)*F(k)) +
            NO(Fd(i)*Fd(j)*F(k)*F(l)))
    expr = wicks(F(a)*NO(F(b)*Fd(c))*Fd(d))
    assert (expr ==
           -KroneckerDelta(a,c)*NO(F(b)*Fd(d)) -
            KroneckerDelta(b,d)*NO(F(a)*Fd(c)) -
            KroneckerDelta(a,c)*KroneckerDelta(b,d)+
            KroneckerDelta(a,d)*NO(F(b)*Fd(c)) +
            NO(F(a)*F(b)*Fd(c)*Fd(d)))


def test_NO():
    i,j,k,l = symbols('ijkl',below_fermi=True)
    a,b,c,d = symbols('abcd',above_fermi=True)
    p,q,r,s = symbols('pqrs', dummy=True)

    assert (NO(Fd(p)*F(q) + Fd(a)*F(b))==
       NO(Fd(p)*F(q)) + NO(Fd(a)*F(b)))
    assert (NO(Fd(i)*NO(F(j)*Fd(a))) ==
       NO(Fd(i)*F(j)*Fd(a)))
    assert NO(1) == 1
    assert NO(i) == i
    assert (NO(Fd(a)*Fd(b)*(F(c)+F(d))) ==
               NO(Fd(a)*Fd(b)*F(c)) +
               NO(Fd(a)*Fd(b)*F(d)))

    assert NO(Fd(a)*F(i))._remove_brackets()==Fd(a)*F(i)
    assert NO(Fd(a)*F(b))._remove_brackets()==Fd(a)*F(b)
    assert NO(F(j)*Fd(i))._remove_brackets()==F(j)*Fd(i)

    assert (NO(Fd(p)*F(q)).subs(Fd(p),Fd(a)+Fd(i)) ==
            NO(Fd(a)*F(q)) + NO(Fd(i)*F(q)))
    assert (NO(Fd(p)*F(q)).subs(F(q),F(a)+F(i)) ==
            NO(Fd(p)*F(a)) + NO(Fd(p)*F(i)))


    expr = NO(Fd(p)*F(q))._remove_brackets()
    assert wicks(expr) == NO(expr)

    assert NO(Fd(a)*F(b)) == - NO(F(b)*Fd(a))

    no =  NO(Fd(a)*F(i)*Fd(j)*F(b))
    assert no[0] == Fd(a)
    assert no[1] == F(i)
    assert no[2] == Fd(j)
    assert no[3] == F(b)
    l1 = [ ind for ind in no.iter_q_creators() ]
    assert l1 == [0,1]
    l2 = [ ind for ind in no.iter_q_annihilators() ]
    assert l2 == [3,2]
    assert no.get_subNO(1) == NO(Fd(a)*Fd(j)*F(b))


def test_contraction():
    i,j,k,l = symbols('ijkl',below_fermi=True)
    a,b,c,d = symbols('abcd',above_fermi=True)
    p,q,r,s = symbols('pqrs')
    assert contraction(Fd(i),F(j)) == KroneckerDelta(i,j)
    assert contraction(F(a),Fd(b)) == KroneckerDelta(a,b)
    assert contraction(F(a),Fd(i)) == 0
    assert contraction(Fd(a),F(i)) == 0
    assert contraction(F(i),Fd(a)) == 0
    assert contraction(Fd(i),F(a)) == 0
    assert contraction(Fd(i),F(p)) == KroneckerDelta(p,i)
    restr = evaluate_deltas(contraction(Fd(p),F(q)))
    assert restr.is_only_below_fermi
    restr = evaluate_deltas(contraction(F(p),Fd(q)))
    assert restr.is_only_above_fermi



def test_Tensors():
    i,j,k,l = symbols('ijkl',below_fermi=True,dummy=True)
    a,b,c,d = symbols('abcd',above_fermi=True,dummy=True)
    p,q,r,s = symbols('pqrs')

    AT= AntiSymmetricTensor
    assert AT('t',(a,b),(i,j)) == -AT('t',(b,a),(i,j))
    assert AT('t',(a,b),(i,j)) ==  AT('t',(b,a),(j,i))
    assert AT('t',(a,b),(i,j)) == -AT('t',(a,b),(j,i))
    assert AT('t',(a,a),(i,j)) == 0
    assert AT('t',(a,b),(i,i)) == 0
    assert AT('t',(a,b,c),(i,j)) == -AT('t',(b,a,c),(i,j))
    assert AT('t',(a,b,c),(i,j,k)) == AT('t',(b,a,c),(i,k,j))

    tabij = AT('t',(a,b),(i,j))
    assert a in tabij
    assert b in tabij
    assert i in tabij
    assert j in tabij
    assert tabij.subs(b,c) == AT('t',(a,c),(i,j))
    assert (2*tabij).subs(i,c) == 2*AT('t',(a,b),(c,j))


def test_fully_contracted():
    i,j,k,l = symbols('ijkl',below_fermi=True)
    a,b,c,d = symbols('abcd',above_fermi=True)
    p,q,r,s = symbols('pqrs', dummy=True)

    Fock = (AntiSymmetricTensor('f',(p,),(q,))*
            NO(Fd(p)*F(q)))
    V = (AntiSymmetricTensor('v',(p,q),(r,s))*
            NO(Fd(p)*Fd(q)*F(s)*F(r)))/4

    Fai=wicks(NO(Fd(i)*F(a))*Fock,
            keep_only_fully_contracted=True,
            simplify_kronecker_deltas=True)
    assert Fai == AntiSymmetricTensor('f',(a,),(i,))
    Vabij=wicks(NO(Fd(i)*Fd(j)*F(b)*F(a))*V,
            keep_only_fully_contracted=True,
            simplify_kronecker_deltas=True)
    assert Vabij==AntiSymmetricTensor('v',(a,b),(i,j))
