from sympy import *
from sympy.printing.mathml import mathml, MathMLPrinter
from xml.dom.minidom import parseString

x = Symbol('x')
mp = MathMLPrinter()

def test_mathml_core():
    mml_1 = mp._print(1+x)
    assert mml_1.nodeName == 'apply'
    nodes = mml_1.childNodes
    assert len(nodes) == 3
    assert nodes[0].nodeName == 'plus'
    assert nodes[0].hasChildNodes() == False
    assert nodes[0].nodeValue is None
    assert nodes[1].nodeName == 'cn'
    assert nodes[1].childNodes[0].nodeValue == '1'
    assert nodes[2].childNodes[0].nodeValue == 'x'

    mml_2 = mp._print(x**2)
    assert mml_2.nodeName == 'apply'
    nodes = mml_2.childNodes
    assert nodes[1].childNodes[0].nodeValue == 'x'
    assert nodes[2].childNodes[0].nodeValue == '2'

    mml_3 = mp._print(2*x)
    assert mml_3.nodeName == 'apply'
    nodes = mml_3.childNodes
    assert nodes[0].nodeName == 'times'
    assert nodes[1].childNodes[0].nodeValue == '2'
    assert nodes[2].childNodes[0].nodeValue == 'x'

def test_mathml_functions():
    mml_1 = mp._print(sin(x))
    assert mml_1.nodeName == 'apply'
    assert mml_1.childNodes[0].nodeName == 'sin'
    assert mml_1.childNodes[1].nodeName == 'ci'

    mml_2 = mp._print(diff(sin(x), x, evaluate=False))
    assert mml_2.nodeName == 'apply'
    assert mml_2.childNodes[0].nodeName == 'diff'
    assert mml_2.childNodes[1].nodeName == 'bvar'
    assert mml_2.childNodes[1].childNodes[0].nodeName == 'ci'  # below bvar there's <ci>x/ci>

def _test_mathml_limits():
    # XXX No unevaluated limits
    mml_1 = mp._print(limit(sin(x)/x, x, 0, evaluate=False))
    assert mml_1.childNodes[0].nodeName == 'limit'
    assert mml_1.childNodes[1].nodeName == 'bvar'
    assert mml_1.childNodes[1].childNodes[0].nodeName == 'ci'

def test_mathml_integrals():
    pass #TODO

def test_mathml_matrices():
    pass #TODO

def test_c2p():
    """This tests some optional routines that depend on libxslt1 (which is optional)"""
    try:
        from sympy.modules.mathml import c2p
        #assert c2p(f.mathml) == result
    except ImportError:
        pass
