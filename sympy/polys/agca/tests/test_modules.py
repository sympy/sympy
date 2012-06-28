"""Test modules.py code."""

from sympy.polys.agca.modules import FreeModule, ModuleOrder, FreeModulePolyRing
from sympy.polys import CoercionFailed, QQ, lex, grlex, ilex, ZZ
from sympy.abc import x, y, z
from sympy.utilities.pytest import raises
from sympy import S

def test_FreeModuleElement():
    e = QQ[x].free_module(3).convert([1, x, x**2])
    f = [QQ[x].convert(1), QQ[x].convert(x), QQ[x].convert(x**2)]
    assert list(e) == f
    assert f[0] == e[0]
    assert f[1] == e[1]
    assert f[2] == e[2]
    raises(IndexError, lambda: e[3])

def test_FreeModule():
    M1 = FreeModule(QQ[x], 2)
    assert M1 == FreeModule(QQ[x], 2)
    assert M1 != FreeModule(QQ[y], 2)
    assert M1 != FreeModule(QQ[x], 3)
    M2 = FreeModule(QQ.poly_ring(x, order="ilex"), 2)

    assert [x, 1] in M1
    assert [x] not in M1
    assert [2, y] not in M1
    assert [1/(x + 1), 2] not in M1

    e = M1.convert([x, x**2 + 1])
    X = QQ[x].convert(x)
    assert e == [X, X**2 + 1]
    assert e == [x, x**2 + 1]
    assert 2*e == [2*x, 2*x**2 + 2]
    assert e*2 == [2*x, 2*x**2 + 2]
    assert e/2 == [x/2, (x**2 + 1)/2]
    assert x*e == [x**2, x**3 + x]
    assert e*x == [x**2, x**3 + x]
    assert X*e == [x**2, x**3 + x]
    assert e*X == [x**2, x**3 + x]

    assert [x, 1] in M2
    assert [x] not in M2
    assert [2, y] not in M2
    assert [1/(x + 1), 2] in M2

    e = M2.convert([x, x**2 + 1])
    X = QQ.poly_ring(x, order="ilex").convert(x)
    assert e == [X, X**2 + 1]
    assert e == [x, x**2 + 1]
    assert 2*e == [2*x, 2*x**2 + 2]
    assert e*2 == [2*x, 2*x**2 + 2]
    assert e/2 == [x/2, (x**2 + 1)/2]
    assert x*e == [x**2, x**3 + x]
    assert e*x == [x**2, x**3 + x]
    assert e/(1 + x) == [x/(1 + x), (x**2 + 1)/(1 + x)]
    assert X*e == [x**2, x**3 + x]
    assert e*X == [x**2, x**3 + x]

    M3 = FreeModule(QQ[x, y], 2)
    assert M3.convert(e) == M3.convert([x, x**2 + 1])

    raises(NotImplementedError, lambda: ZZ[x].free_module(2))
    raises(NotImplementedError, lambda: FreeModulePolyRing(ZZ, 2))

def test_ModuleOrder():
    o1 = ModuleOrder(lex, grlex)
    o2 = ModuleOrder(ilex, lex)

    assert o1 == ModuleOrder(lex, grlex)
    assert (o1 != ModuleOrder(lex, grlex)) is False
    assert o1 != o2

    assert o1((1, 2, 3)) == (1, (5, (2, 3)))
    assert o2((1, 2, 3)) == (-1, (2, 3))

def test_SubModulePolyRing_global():
    R = QQ[x, y]
    F = R.free_module(3)
    Fd = F.submodule([1, 0, 0], [1, 2, 0], [1, 2, 3])
    M = F.submodule([x**2 + y**2, 1, 0], [x, y, 1])

    assert F == Fd
    assert Fd == F
    assert F != M
    assert M != F
    assert Fd != M
    assert M != Fd
    assert Fd == F.submodule(*F.basis())

    assert Fd.is_full_module()
    assert not M.is_full_module()
    assert not Fd.is_zero()
    assert not M.is_zero()
    assert Fd.submodule().is_zero()

    assert M.contains([x**2 + y**2 + x, 1 + y, 1])
    assert not M.contains([x**2 + y**2 + x, 1 + y, 2])
    assert M.contains([y**2, 1 - x*y, -x])

    assert not F.submodule([1 + x, 0, 0]) == F.submodule([1, 0, 0])
    assert F.submodule([1, 0, 0], [0, 1, 0]).union(F.submodule([0, 0, 1])) == F

    m = F.convert([x**2 + y**2, 1, 0])
    n = M.convert(m)
    assert m.module is F
    assert n.module is M

    raises(ValueError, lambda: M.submodule([1, 0, 0]))

def test_SubModulePolyRing_local():
    R = QQ.poly_ring(x, y, order=ilex)
    F = R.free_module(3)
    Fd = F.submodule([1+x, 0, 0], [1+y, 2+2*y, 0], [1, 2, 3])
    M = F.submodule([x**2 + y**2, 1, 0], [x, y, 1])

    assert F == Fd
    assert Fd == F
    assert F != M
    assert M != F
    assert Fd != M
    assert M != Fd
    assert Fd == F.submodule(*F.basis())

    assert Fd.is_full_module()
    assert not M.is_full_module()
    assert not Fd.is_zero()
    assert not M.is_zero()
    assert Fd.submodule().is_zero()

    assert M.contains([x**2 + y**2 + x, 1 + y, 1])
    assert not M.contains([x**2 + y**2 + x, 1 + y, 2])
    assert M.contains([y**2, 1 - x*y, -x])

    assert F.submodule([1 + x, 0, 0]) == F.submodule([1, 0, 0])
    assert F.submodule([1, 0, 0], [0, 1, 0]).union(F.submodule([0, 0, 1 + x*y])) == F

    raises(ValueError, lambda: M.submodule([1, 0, 0]))

def test_SubModulePolyRing_nontriv_global():
    R = QQ[x, y, z]
    F = R.free_module(1)
    def contains(I, f):
        return F.submodule(*[[g] for g in I]).contains([f])

    assert contains([x, y], x)
    assert contains([x, y], x + y)
    assert not contains([x, y], 1)
    assert not contains([x, y], z)
    assert contains([x**2 + y, x**2 + x], x - y)
    assert not contains([x+y+z, x*y+x*z+y*z, x*y*z], x**2)
    assert contains([x+y+z, x*y+x*z+y*z, x*y*z], x**3)
    assert contains([x+y+z, x*y+x*z+y*z, x*y*z], x**4)
    assert not contains([x+y+z, x*y+x*z+y*z, x*y*z], x*y**2)
    assert contains([x+y+z, x*y+x*z+y*z, x*y*z], x**4 + y**3 + 2*z*y*x)
    assert contains([x+y+z, x*y+x*z+y*z, x*y*z], x*y*z)
    assert contains([x, 1+x+y, 5-7*y], 1)
    assert contains([x**3+y**3, y**3+z**3, z**3+x**3, x**2*y + x**2*z + y**2*z],
                    x**3)
    assert not contains([x**3+y**3, y**3+z**3, z**3+x**3, x**2*y + x**2*z + y**2*z],
                        x**2 + y**2)

    # compare local order
    assert not contains([x*(1+x+y), y*(1+z)], x)
    assert not contains([x*(1+x+y), y*(1+z)], x + y)

def test_SubModulePolyRing_nontriv_local():
    R = QQ.poly_ring(x, y, z, order=ilex)
    F = R.free_module(1)
    def contains(I, f):
        return F.submodule(*[[g] for g in I]).contains([f])

    assert contains([x, y], x)
    assert contains([x, y], x + y)
    assert not contains([x, y], 1)
    assert not contains([x, y], z)
    assert contains([x**2 + y, x**2 + x], x - y)
    assert not contains([x+y+z, x*y+x*z+y*z, x*y*z], x**2)
    assert contains([x*(1+x+y), y*(1+z)], x)
    assert contains([x*(1+x+y), y*(1+z)], x + y)

def test_syzygy():
    R = QQ[x, y, z]
    M = R.free_module(1).submodule([x*y], [y*z], [x*z])
    S = R.free_module(3).submodule([0, x, -y], [z, -x, 0])
    assert M.syzygy_module() == S

    F = R.free_module(3)
    assert F.submodule(*F.basis()).syzygy_module() == F.submodule()

def test_in_terms_of_generators():
    R = QQ.poly_ring(x, order="ilex")
    M = R.free_module(2).submodule([2*x, 0], [1, 2])
    assert M.in_terms_of_generators([x, x]) == [R.convert(S(1)/4), R.convert(x/2)]
    raises(ValueError, lambda: M.in_terms_of_generators([1, 0]))
