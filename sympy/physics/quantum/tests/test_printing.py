# -*- encoding: utf-8 -*-

from sympy.printing.latex import latex

from sympy.printing import sstr, srepr
def sT(expr, string):
    assert srepr(expr) == string
    assert eval(string) == expr

from sympy.printing.pretty import pretty as xpretty
def pretty(expr, order=None):
    """ASCII pretty-printing"""
    return xpretty(expr, order=order, use_unicode=False)

def upretty(expr, order=None):
    """Unicode pretty-printing"""
    return xpretty(expr, order=order, use_unicode=True)

from sympy import S, symbols, Symbol, Mul, Pow
from sympy import Rational, Integer

from sympy.physics.quantum.anticommutator import AntiCommutator
from sympy.physics.quantum.commutator import Commutator
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.qubit import Qubit
from sympy.physics.quantum.gate import IdentityGate
from sympy.physics.quantum.state import Bra, Ket, TimeDepKet, TimeDepBra
from sympy.physics.quantum.spin import JzKet, JzBra, JzOp, JxOp, JyKet, Jz, Jx
from sympy.physics.quantum.operator import Operator, OuterProduct
from sympy.physics.quantum.innerproduct import InnerProduct
from sympy.physics.quantum.tensorproduct import TensorProduct

A, B = symbols('A,B', commutative=False)

def test_anticommutator():
    ac = AntiCommutator(A,B)
    assert pretty(ac) == '{A,B}'
    assert upretty(ac) == u'{A,B}'
    assert latex(ac) == r'\left\{A,B\right\}'
    sT(ac, 'AntiCommutator(A,B)')


def test_commutator():
    A, B = symbols('A,B', commutative=False)
    c = Commutator(A,B)
    assert pretty(c) == '[A,B]'
    assert upretty(c) == u'[A,B]'
    assert latex(c) == r'\left[A,B\right]'
    sT(c, 'Commutator(A,B)')

def test_gate():
    q = Qubit('10101')
    g1 = IdentityGate(2)
    assert pretty(g1) == '1 \n 2'
    assert upretty(g1) == u'1 \n 2'
    #TODO issue 2510 assert pretty(g1*q) == '1 *|10101>\n 2        '
    assert upretty(g1*q) == u'1 ⋅❘10101⟩\n 2        '
    assert latex(g1) == r'1_{2}'
    assert latex(g1*q) == r'1_{2} {\left| 10101 \right\rangle }'
    sT(g1, 'IdentityGate(Integer(2))')
    sT(q, 'Qubit(Integer(1),Integer(0),Integer(1),Integer(0),Integer(1))')
    sT(g1*q, 'Mul(IdentityGate(Integer(2)), Qubit(Integer(1),Integer(0),Integer(1),Integer(0),Integer(1)))')

def test_inner_product():
    psi = Ket('psi')
    ip = Dagger(psi)*psi
    #TODO issue 2510 assert pretty(ip) == '<psi|psi>'
    assert upretty(ip) == u'⟨ψ❘ψ⟩'
    assert latex(ip) == r'\left\langle \psi \middle| \psi \right\rangle'
    sT(ip, "InnerProduct(Bra(Symbol('psi')),Ket(Symbol('psi')))")
    psi = TimeDepKet('psi','t')
    ip = Dagger(psi)*psi
    #TODO issue 2510 assert pretty(ip) == '<psi;t|psi;t>'
    assert upretty(ip) == u'⟨ψ;t❘ψ;t⟩'
    assert latex(ip) == r'\left\langle \psi ; t \middle| \psi ; t \right\rangle'
    sT(ip, "InnerProduct(TimeDepBra(Symbol('psi'),Symbol('t')),TimeDepKet(Symbol('psi'),Symbol('t')))")
    jket = JzKet(S(1)/2,S(1)/2)
    ip = Dagger(jket)*jket
    #TODO issue 2510 assert pretty(ip) == '<z:1/2,1/2|z:1/2,1/2>'
    assert upretty(ip) == u'⟨z:1/2,1/2❘z:1/2,1/2⟩'
    assert latex(ip) == r"\left\langle z:\frac{1}{2},\frac{1}{2} \middle| z:\frac{1}{2},\frac{1}{2} \right\rangle"
    sT(ip, 'InnerProduct(JzBra(Rational(1, 2),Rational(1, 2)),JzKet(Rational(1, 2),Rational(1, 2)))')

def test_operator():
    k = Ket('k')
    b = Bra('b')
    op = OuterProduct(k, b)
    A = Operator('A')
    D = Operator('D', Symbol('t'), S(1)/2)

    assert pretty(A) == 'A'
    assert upretty(A) == u'A'
    assert latex(A) == 'A'
    sT(A, "Operator(Symbol('A'))")
    #Dagger
    #TODO issue 2510 assert pretty(Dagger(A)) == ' +\nA '
    assert upretty(Dagger(A)) == u' †\nA '
    assert latex(Dagger(A)) == r'A^{\dag}'
    sT(Dagger(A), "Dagger(Operator(Symbol('A')))")
    #Outer product
    #TODO issue 2510 assert pretty(op) == u'|k><b|'
    assert upretty(op) == u'❘k⟩⟨b❘'
    assert latex(op) == r'{\left| k \right\rangle }{\left\langle b \right| }'
    sT(op, "OuterProduct(Ket(Symbol('k')),Bra(Symbol('b')))")
    #Inverse
    assert pretty(A.inv()) == u'1\n-\nA'
    assert upretty(A.inv()) == u'1\n─\nA'
    assert latex(A.inv()) == r'\frac{1}{A}'
    sT(A.inv(), "Pow(Operator(Symbol('A')), Integer(-1))")
    #two args
    assert pretty(D) == 'Operator(D,t,1/2)'
    assert upretty(D) == u'Operator(D,t,1/2)'
    assert latex(D) == r'Operator(D,t,\frac{1}{2})'
    sT(D,"Operator(Symbol('D'),Symbol('t'),Rational(1, 2))")

def test_qubit():
    q = Qubit('0101')
    #TODO issue 2510 assert pretty(q) == u'|0101>'
    assert upretty(q) == u'❘0101⟩'
    assert latex(q) == r"{\left| 0101 \right\rangle }"
    sT(q,'Qubit(Integer(0),Integer(1),Integer(0),Integer(1))')

def test_spin():
    jket = JzKet(S(1)/2,S(1)/2)
    #TODO issue 2510 assert pretty(jket) == u'|z:1/2,1/2>'
    assert upretty(jket) == u'❘z:1/2,1/2⟩'
    assert latex(jket) == r"{\left| z:\frac{1}{2},\frac{1}{2} \right\rangle }"
    sT(jket, 'JzKet(Rational(1, 2),Rational(1, 2))')

def test_ket():
    psi = Ket('psi')
    #TODO issue 2510 assert pretty(psi) == '|psi>'
    assert upretty(psi) == u'❘ψ⟩'
    #TODO issue 2510 assert pretty(Dagger(psi)) == u'<psi|'
    assert upretty(Dagger(psi)) == u'⟨ψ❘'
    assert latex(psi) == r'{\left| \psi \right\rangle }'
    assert latex(Dagger(psi)) == r'{\left\langle \psi \right| }'
    sT(psi, "Ket(Symbol('psi'))")
    sT(Dagger(psi), "Bra(Symbol('psi'))")

def test_tensor_product():
    pr = TensorProduct(Jz,Jx)*TensorProduct(JzKet(1,0),JyKet(1,1))
    #TODO issue 2510 assert upretty(pr) == u'J @ J *❘z:1,0⟩@ ❘y:1,1⟩\n z   x                 '
    assert upretty(pr) == u'J ⨂ J ⋅❘z:1,0⟩⨂ ❘y:1,1⟩\n z   x                 '
    assert latex(pr) == r"{J_z}\otimes {J_x} {{\left| z:1,0 \right\rangle }}\otimes {{\left| y:1,1 \right\rangle }}"
    sT(pr, "Mul(TensorProduct(JzOp(Symbol('J')), JxOp(Symbol('J'))), TensorProduct(JzKet(Integer(1),Integer(0)), JyKet(Integer(1),Integer(1))))")
