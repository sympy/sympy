from sympy import (S, Tuple, symbols, Interval, EmptySequence, oo, SeqPer\
                   , SeqFormula, SeqFunc, Lambda, sequence, SeqAdd, SeqMul)
from sympy.series.sequences import SeqExpr, SeqExprOp
from sympy.utilities.pytest import raises

x, y, z = symbols('x y z')
n, m = symbols('n m')


def test_EmptySequence():
    assert isinstance(S.EmptySequence, EmptySequence)

    assert S.EmptySequence.interval is S.EmptySet
    assert S.EmptySequence.length is S.Zero

    assert list(S.EmptySequence) == []


def test_SeqExpr():
    s = SeqExpr((1, 2, 3), (0, 10))

    assert isinstance(s, SeqExpr)
    assert s.gen == Tuple(1, 2, 3)
    assert s.interval == Interval(0, 10)
    assert s.start == 0
    assert s.stop == 10
    assert s.length == 11

    assert SeqExpr((1, 2, 3), (0, 10, 2)).length == 6
    assert SeqExpr((1, 2, 3), (0, oo)).length is oo

    assert SeqExpr((1, 2, 3), (oo, -oo)) is S.EmptySequence

    raises(ValueError, lambda: SeqExpr((1, 2, 3), (0, 1, 2, 3)))
    raises(ValueError, lambda: SeqExpr((1, 2, 3), (-oo, oo)))
    raises(ValueError, lambda: SeqExpr((1, 2, 3), (0, oo, oo)))


def test_SeqPer():
    s = SeqPer((1, 2, 3), (0, 5))

    assert isinstance(s, SeqPer)
    assert s.periodical == Tuple(1, 2, 3)
    assert s.period == 3
    assert s.coeff(3) == 1

    assert list(s) == [1, 2, 3, 1, 2, 3]
    assert s[:] == [1, 2, 3, 1, 2, 3]
    assert SeqPer((1, 2, 3), (0, 5, 2))[:] == [1, 3, 2]
    assert SeqPer((1, 2, 3), (-oo, 0))[0:6] == [1, 2, 3, 1, 2, 3]


def test_SeqFormula():
    s = SeqFormula((n**2, n), (0, 5))

    assert isinstance(s, SeqFormula)
    assert s.formula == n**2
    assert s.coeff(3) == 9

    assert list(s) == [i**2 for i in range(6)]
    assert s[:] == [i**2 for i in range(6)]
    assert SeqFormula((n**2, n), (0, 5, 2))[:] == [0, 4, 16]
    assert SeqFormula((n**2, n), (-oo, 0))[0:6] == [i**2 for i in range(6)]

    assert SeqFormula(n**2, (0, oo)) == SeqFormula((n**2, n), (0, oo))

    assert SeqFormula(n**2, (0, m)).subs(m, x) == SeqFormula(n**2, (0, x))
    assert SeqFormula((m*n**2, n), (0, oo)).subs(m, x) == \
        SeqFormula((x*n**2, n), (0, oo))


def test_SeqFunc():
    s = SeqFunc(Lambda(n, n**2), (0, 5))

    assert isinstance(s, SeqFunc)
    assert s.function == Lambda(n, n**2)
    assert s.coeff(3) == 9

    assert list(s) == [i**2 for i in range(6)]
    assert s[:] == [i**2 for i in range(6)]
    assert SeqFunc(Lambda(n, n**2), (0, 5, 2))[:] == [0, 4, 16]
    assert SeqFunc(Lambda(n, n**2), (-oo, 0))[0:6] == [i**2 for i in range(6)]


def test_sequence():
    form = SeqFormula((n**2, n), (0, 5))
    per = SeqPer((1, 2, 3), (0, 5))
    func = SeqFunc(Lambda(n, n**2), (0, 5))
    inter = SeqFormula((n**2, n))

    assert sequence(formula=(n**2, n), interval=(0, 5)) == form
    assert sequence(periodical=(1, 2, 3), interval=(0, 5)) == per
    assert sequence(func=Lambda(n, n**2), interval=(0, 5)) == func
    assert sequence(formula=(n**2, n)) == inter


def test_SeqExprOp():
    form = SeqFormula((n**2, n), (0, 10))
    per = SeqPer((1, 2, 3), (5, 10))
    func = SeqFunc(Lambda(m, m**2), (0, 10))

    s = SeqExprOp(form, per, func)
    assert s.gen == ((n**2, n), (1, 2, 3), Lambda(m, m**2))
    assert s.interval == Interval(5, 10)
    assert s.start == 5
    assert s.stop == 10
    assert s.length == 6
    assert s.variables == (n, m)


def test_SeqAdd():
    per = SeqPer((1, 2, 3))
    form = SeqFormula(n**2)
    func = SeqFunc(Lambda(n, n**2))

    per_bou = SeqPer((1, 2), (1, 5))
    form_bou = SeqFormula(n**2, (6, 10))
    func_bou = SeqFunc(Lambda(n, n**2), (0, 5))

    assert SeqAdd() == S.EmptySequence
    assert SeqAdd(S.EmptySequence) == S.EmptySequence
    assert SeqAdd(per) == per
    assert SeqAdd(per, S.EmptySequence) == per
    assert SeqAdd(per_bou, form_bou) == S.EmptySequence

    s = SeqAdd(per_bou, func_bou, evaluate=False)
    assert s.args == (func_bou, per_bou)
    assert s[:] == [2, 6, 10, 18, 26]
    assert list(s) == [2, 6, 10, 18, 26]

    assert isinstance(SeqAdd(per, per_bou, evaluate=False), SeqAdd)

    s1 =  SeqAdd(per, per_bou)
    assert isinstance(s1, SeqPer)
    assert s1 == SeqPer((2, 4, 4, 3, 3, 5), (1, 5))
    s2 = SeqAdd(func, func_bou)
    assert isinstance(s2, SeqFunc)
    assert s2 == SeqFunc(Lambda(n, 2*n**2), (0, 5))
    s3 = SeqAdd(form, form_bou)
    assert isinstance(s3, SeqFormula)
    assert s3 == SeqFormula(2*n**2, (6, 10))

    assert SeqAdd(form, form_bou, per) == \
        SeqAdd(per, SeqFormula(2*n**2, (6, 10)))
    assert SeqAdd(form, SeqAdd(form_bou, per)) == \
        SeqAdd(per, SeqFormula(2*n**2, (6, 10)))
    assert SeqAdd(per, SeqAdd(form, func), evaluate=False) == \
        SeqAdd(per, form, func)

    assert SeqAdd(SeqPer((1, 2)), SeqPer((1, 2))) == SeqPer((2, 4))


def test_SeqMul():
    per = SeqPer((1, 2, 3))
    form = SeqFormula(n**2)
    func = SeqFunc(Lambda(n, n**2))

    per_bou = SeqPer((1, 2), (1, 5))
    form_bou = SeqFormula(n**2, (6, 10))
    func_bou = SeqFunc(Lambda(n, n**2), (0, 5))

    assert SeqMul() == S.EmptySequence
    assert SeqMul(S.EmptySequence) == S.EmptySequence
    assert SeqMul(per) == per
    assert SeqMul(per, S.EmptySequence) == S.EmptySequence
    assert SeqMul(per_bou, form_bou) == S.EmptySequence

    s = SeqMul(per_bou, func_bou, evaluate=False)
    assert s.args == (func_bou, per_bou)
    assert s[:] == [1, 8, 9, 32, 25]
    assert list(s) == [1, 8, 9, 32, 25]

    assert isinstance(SeqMul(per, per_bou, evaluate=False), SeqMul)

    s1 =  SeqMul(per, per_bou)
    assert isinstance(s1, SeqPer)
    assert s1 == SeqPer((1, 4, 3, 2, 2, 6), (1, 5))
    s2 = SeqMul(func, func_bou)
    assert isinstance(s2, SeqFunc)
    assert s2 == SeqFunc(Lambda(n, n**4), (0, 5))
    s3 = SeqMul(form, form_bou)
    assert isinstance(s3, SeqFormula)
    assert s3 == SeqFormula(n**4, (6, 10))

    assert SeqMul(form, form_bou, per) == \
        SeqMul(per, SeqFormula(n**4, (6, 10)))
    assert SeqMul(form, SeqMul(form_bou, per)) == \
        SeqMul(per, SeqFormula(n**4, (6, 10)))
    assert SeqMul(per, SeqMul(form, func), evaluate=False) == \
        SeqMul(per, form, func)

    assert SeqMul(SeqPer((1, 2)), SeqPer((1, 2))) == SeqPer((1, 4))


def test_add():
    per = SeqPer((1, 2))
    form = SeqFormula(n**2)
    func = SeqFunc(Lambda(n, n**2))

    assert per.add(SeqPer((2, 3))) == SeqPer((3, 5))
    assert form.add(SeqFormula(n**3)) == SeqFormula(n**2 + n**3)
    assert func.add(SeqFunc(Lambda(n, n**3))) == \
                                        SeqFunc(Lambda(n, n**2 + n**3))

    assert per.add(form) == SeqAdd(per, form)
    assert per.add(form).add(func) == SeqAdd(per, form, func)


def test_sub():
    per = SeqPer((1, 2))
    form = SeqFormula(n**2)
    func = SeqFunc(Lambda(n, n**2))

    assert per.sub(SeqPer((2, 3))) == SeqPer((-1, -1))
    assert form.sub(SeqFormula(n**3)) == SeqFormula(n**2 - n**3)
    assert func.sub(SeqFunc(Lambda(n, n**3))) == \
                                        SeqFunc(Lambda(n, n**2 - n**3))

    assert per.sub(form) == SeqAdd(per, -form)
    assert per.sub(form).sub(func) == SeqAdd(per, -form, -func)


def test_mul__coeff_mul():
    assert SeqPer((1, 2)).coeff_mul(2) == SeqPer((2, 4))
    assert SeqFormula(n**2).coeff_mul(2) == SeqFormula(2*n**2)
    assert SeqFunc(Lambda(n, n**2)).coeff_mul(2) == SeqFunc(Lambda(n, 2*n**2))
    assert S.EmptySequence.coeff_mul(100) == S.EmptySequence

    assert SeqPer((1, 2)).mul(SeqPer((2, 3))) == SeqPer((2, 6))
    assert SeqFormula(n**2).mul(SeqFormula(n**3)) == SeqFormula(n**5)
    assert SeqFunc(Lambda(n, n**2)).mul(SeqFunc(Lambda(n, n**3))) == \
                                        SeqFunc(Lambda(n, n**5))

    assert S.EmptySequence.mul(SeqFormula(n**2)) == S.EmptySequence
    assert SeqFormula(n**2).mul(S.EmptySequence) == S.EmptySequence


def test_neg():
    assert -SeqPer((1, -2)) == SeqPer((-1, 2))
    assert -SeqFormula(n**2) == SeqFormula(-n**2)
    assert -SeqFunc(Lambda(n, n**2)) == SeqFunc(Lambda(n, -n**2))


def test_operations():
    per = SeqPer((1, 2))
    per2 = SeqPer((2, 4))
    form = SeqFormula(n**2)
    func = SeqFunc(Lambda(n, n**2))

    assert per + form + func == SeqAdd(per, form, func)
    assert per + form - func == SeqAdd(per, form, -func)
    assert per + form - S.EmptySequence == SeqAdd(per, form)
    assert per + per2 + form == SeqAdd(SeqPer((3, 6)), form)
    assert S.EmptySequence - per == -per
    assert form + form == SeqFormula(2*n**2)

    assert per * form * func == SeqMul(per, form, func)
    assert form * form == SeqFormula(n**4)
    assert form * -form == SeqFormula(-n**4)

    assert form * (per + func) == SeqMul(form, SeqAdd(per, func))
    assert form * (per + per) == SeqMul(form, per2)

    assert form.coeff_mul(n) == SeqFormula(n**3)
    assert per.coeff_mul(n) == SeqPer((n, 2*n))
    assert func.coeff_mul(n) == SeqFunc(Lambda(n, n**3))
