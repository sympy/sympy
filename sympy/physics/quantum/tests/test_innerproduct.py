from sympy import I, Integer, Symbol

from sympy.physics.quantum.innerproduct import InnerProduct
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.state import Bra, Ket, StateBase, TimeDepKet
from sympy.physics.quantum.spin import JzKet, JzBra

def test_innerproduct():
    k = Ket('k')
    b = Bra('b')
    ip = InnerProduct(b,k)
    assert isinstance(ip, InnerProduct)
    assert ip.bra == b
    assert ip.ket == k
    assert b*k == InnerProduct(b,k)
    assert k*(b*k)*b == k*InnerProduct(b,k)*b
    assert InnerProduct(b,k).subs(b,Dagger(k)) == Dagger(k)*k

def test_innerproduct_dagger():
    k = Ket('k')
    b = Bra('b')
    ip = b*k
    assert Dagger(ip) == Dagger(k)*Dagger(b)


class FooState(StateBase):
    pass


class FooKet(Ket, FooState):

    @property
    def dual_class(self):
        return FooBra

    def _eval_innerproduct_FooBra(self, bra):
        return Integer(1)

    def _eval_innerproduct_BarBra(self, bra):
        return I


class FooBra(Bra, FooState):
    @property
    def dual_class(self):
        return FooKet


class BarState(StateBase):
    pass


class BarKet(Ket, BarState):
    @property
    def dual_class(self):
        return BarBra


class BarBra(Bra, BarState):
    @property
    def dual_class(self):
        return BarKet


def test_doit():
    f = FooKet('foo')
    b = BarBra('bar')
    assert InnerProduct(b,f).doit() == I
    assert InnerProduct(Dagger(f),Dagger(b)).doit() == -I
    assert InnerProduct(Dagger(f),f).doit() == Integer(1)


def test_printing():
    from sympy import pretty, latex, srepr
    psi = Ket('psi')
    ip = Dagger(psi)*psi
    assert pretty(ip, use_unicode=True) == u'\u27e8\u03c8\u2758\u03c8\u27e9'
    assert latex(ip) == r"\left\langle \psi \middle| \psi \right\rangle"
    psi = TimeDepKet('psi','t')
    ip = Dagger(psi)*psi
    assert pretty(ip, use_unicode=True) == u'\u27e8\u03c8;t\u2758\u03c8;t\u27e9'
    assert latex(ip) == r"\left\langle \psi ; t \middle| \psi ; t \right\rangle"
    jket = JzKet('1/2','1/2')
    ip = Dagger(jket)*jket
    assert pretty(ip, use_unicode=True) == u'\u27e8z:1/2,1/2\u2758z:1/2,1/2\u27e9'
    assert latex(ip) == r"\left\langle z:1/2,1/2 \middle| z:1/2,1/2 \right\rangle"
    assert eval(srepr(ip)) == ip
