from sympy.logic.algorithms.euf_theory import _sym, _func, _apply, _Eq, EUFCongruenceClosure

def test_symbol_function_apply():
    x = _sym("x")
    y = _sym("y")
    f = _func("f", 1)
    g = _func("g", 2)
    fx = _apply(f, x)
    gyx = _apply(g, y, x)
    assert x.name == "x"
    assert f.name == "f" and f.arity == 1
    assert g.name == "g" and g.arity == 2
    assert fx.func == f and fx.args == (x,)
    assert gyx.func == g and gyx.args == (y, x)


def test_basic_equality():
    x, y, z = _sym("x"), _sym("y"), _sym("z")
    eqs = [_Eq(x, y), _Eq(y, z)]
    cc = EUFCongruenceClosure(eqs)
    assert cc.are_equal(x, y)
    assert cc.are_equal(y, z)
    assert cc.are_equal(x, z)
    w = _sym("w")
    assert not cc.are_equal(x, w)


def test_function_congruence():
    f = _func("f", 1)
    x, y = _sym("x"), _sym("y")
    fx = _apply(f, x)
    fy = _apply(f, y)
    eqs = [_Eq(x, y), _Eq(fx, _sym("a")), _Eq(fy, _sym("b"))]
    cc = EUFCongruenceClosure(eqs)
    assert cc.are_equal(fx, fy)
    assert cc.are_equal(_sym("a"), _sym("b"))


def test_unrelated_terms():
    f = _func("f", 1)
    g = _func("g", 1)
    x = _sym("x")
    fx = _apply(f, x)
    gx = _apply(g, x)
    eqs = [_Eq(fx, _sym("a")), _Eq(gx, _sym("b"))]
    cc = EUFCongruenceClosure(eqs)
    assert not cc.are_equal(fx, gx)
    assert not cc.are_equal(_sym("a"), _sym("b"))


def test_chain_of_equalities():
    letters = [_sym(chr(ord('a') + i)) for i in range(5)]
    eqs = [_Eq(letters[i], letters[i+1]) for i in range(len(letters)-1)]
    cc = EUFCongruenceClosure(eqs)
    for i in range(len(letters)):
        for j in range(len(letters)):
            assert cc.are_equal(letters[i], letters[j])


def test_function_same_args():
    f = _func("f", 1)
    x = _sym("x")
    fx1 = _apply(f, x)
    fx2 = _apply(f, x)
    cc = EUFCongruenceClosure([_Eq(fx1, fx2)])
    assert cc.are_equal(fx1, fx2)


def test_large_chain():
    letters = [_sym(chr(ord('a') + i)) for i in range(26)]
    eqs = [_Eq(letters[i], letters[i+1]) for i in range(25)]
    cc = EUFCongruenceClosure(eqs)
    for i in range(26):
        for j in range(26):
            assert cc.are_equal(letters[i], letters[j])


def test_function_congruence_with_merge():
    f = _func("f", 2)
    x, y, z = _sym("x"), _sym("y"), _sym("z")
    fxz = _apply(f, x, z)
    fyz = _apply(f, y, z)
    eqs = [_Eq(x, y), _Eq(fxz, _sym("a")), _Eq(fyz, _sym("b"))]
    cc = EUFCongruenceClosure(eqs)
    assert cc.are_equal(fxz, fyz)
    assert cc.are_equal(_sym("a"), _sym("b"))


def test_no_equalities():
    cc = EUFCongruenceClosure([])
    a, b = _sym("a"), _sym("b")
    assert not cc.are_equal(a, b)


def test_internal_structures_and_methods():
    f = _func("f", 1)
    x = _sym("x")
    fx = _apply(f, x)
    eqs = [_Eq(fx, _sym("a"))]
    cc = EUFCongruenceClosure(eqs)
    # Test _id_of and _term_of mapping
    assert cc._id_of[fx] == cc._id_of[fx]
    assert cc._term_of[cc._id_of[fx]] == fx
    # Test union-find
    i = cc._id_of[fx]
    assert cc._find(i) == cc._find(i)
    # Test lookup and use
    key = (f, tuple(cc._id_of[arg] for arg in fx.args))
    assert key in cc.lookup
    assert isinstance(cc.use[cc._id_of[x]], list)

