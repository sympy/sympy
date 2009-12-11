from sympy import diff, Integral, Limit, sin, Symbol, Integer, Rational, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh, E, I, oo, pi, GoldenRatio, EulerGamma
from sympy.printing.mathml import mathml, MathMLPrinter
from xml.dom.minidom import parseString

x = Symbol('x')
mp = MathMLPrinter()

def test_printmethod():
    pass #TODO

def test_mathml_core():
    mml_1 = mp._print(1+x)
    assert mml_1.nodeName == 'apply'
    nodes = mml_1.childNodes
    assert len(nodes) == 3
    assert nodes[0].nodeName == 'plus'
    assert nodes[0].hasChildNodes() == False
    assert nodes[0].nodeValue is None
    assert nodes[1].nodeName in ['cn', 'ci']
    if nodes[1].nodeName == 'cn':
        assert nodes[1].childNodes[0].nodeValue == '1'
        assert nodes[2].childNodes[0].nodeValue == 'x'
    else:
        assert nodes[1].childNodes[0].nodeValue == 'x'
        assert nodes[2].childNodes[0].nodeValue == '1'

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

def test_mathml_limits():
    # XXX No unevaluated limits
    lim_fun = sin(x)/x
    mml_1 = mp._print(Limit(lim_fun, x, 0))
    assert mml_1.childNodes[0].nodeName == 'limit'
    assert mml_1.childNodes[1].nodeName == 'bvar'
    assert mml_1.childNodes[2].nodeName == 'lowlimit'
    assert mml_1.childNodes[3].toxml() == mp._print(lim_fun).toxml()

def test_mathml_integrals():
    integrand = x
    mml_1 = mp._print(Integral(integrand, (x, 0, 1)))
    assert mml_1.childNodes[0].nodeName == 'int'
    assert mml_1.childNodes[1].nodeName == 'bvar'
    assert mml_1.childNodes[2].nodeName == 'lowlimit'
    assert mml_1.childNodes[3].nodeName == 'uplimit'
    assert mml_1.childNodes[4].toxml() == mp._print(integrand).toxml()

def test_mathml_tuples():
    mml_1 = mp._print([2])
    assert mml_1.nodeName == 'list'
    assert mml_1.childNodes[0].nodeName == 'cn'
    assert len(mml_1.childNodes) == 1

    mml_2 = mp._print([2, Integer(1)])
    assert mml_2.nodeName == 'list'
    assert mml_2.childNodes[0].nodeName == 'cn'
    assert mml_2.childNodes[1].nodeName == 'cn'
    assert len(mml_2.childNodes) == 2

def test_mathml_matrices():
    pass #TODO

def test_mathml_add():
    mml = mp._print(x**5 - x**4 + x)
    assert mml.childNodes[0].nodeName == 'plus'
    assert mml.childNodes[1].childNodes[0].nodeName == 'minus'
    assert mml.childNodes[1].childNodes[1].nodeName == 'ci'
    assert mml.childNodes[2].childNodes[0].nodeName == 'power'

def test_mathml_Rational():
    mml_1 = mp._print(Rational(1,1))
    """should just return a number"""
    assert mml_1.nodeName == 'cn'

    mml_2 = mp._print(Rational(2,5))
    assert mml_2.childNodes[0].nodeName == 'divide'

def test_mathml_constants():
    mml = mp._print(I)
    assert mml.nodeName == 'imaginaryi'

    mml = mp._print(E)
    assert mml.nodeName == 'exponentiale'

    mml = mp._print(oo)
    assert mml.nodeName == 'infinity'

    mml = mp._print(pi)
    assert mml.nodeName == 'pi'

    assert mathml(GoldenRatio) == '<cn>\xcf\x86</cn>'

    mml = mathml(EulerGamma)
    assert mml == '<eulergamma/>'

def test_mathml_trig():
    mml = mp._print(sin(x))
    assert mml.childNodes[0].nodeName == 'sin'

    mml = mp._print(cos(x))
    assert mml.childNodes[0].nodeName == 'cos'

    mml = mp._print(tan(x))
    assert mml.childNodes[0].nodeName == 'tan'

    mml = mp._print(asin(x))
    assert mml.childNodes[0].nodeName == 'arcsin'

    mml = mp._print(acos(x))
    assert mml.childNodes[0].nodeName == 'arccos'

    mml = mp._print(atan(x))
    assert mml.childNodes[0].nodeName == 'arctan'

    mml = mp._print(sinh(x))
    assert mml.childNodes[0].nodeName == 'sinh'

    mml = mp._print(cosh(x))
    assert mml.childNodes[0].nodeName == 'cosh'

    mml = mp._print(tanh(x))
    assert mml.childNodes[0].nodeName == 'tanh'

    mml = mp._print(asinh(x))
    assert mml.childNodes[0].nodeName == 'arcsinh'

    mml = mp._print(atanh(x))
    assert mml.childNodes[0].nodeName == 'arctanh'

    mml = mp._print(acosh(x))
    assert mml.childNodes[0].nodeName == 'arccosh'

def test_c2p():
    """This tests some optional routines that depend on libxslt1 (which is optional)"""
    try:
        from sympy.modules.mathml import c2p
        #assert c2p(f.mathml) == result
    except ImportError:
        pass

def test_symbol():
    mml = mp._print(Symbol("x"))
    assert mml.nodeName == 'ci'
    assert mml.childNodes[0].nodeValue == 'x'
    del mml

    mml = mp._print(Symbol("x^2"))
    assert mml.nodeName == 'ci'
    assert mml.childNodes[0].nodeName == 'mml:msup'
    assert mml.childNodes[0].childNodes[0].nodeName == 'mml:mi'
    assert mml.childNodes[0].childNodes[0].childNodes[0].nodeValue == 'x'
    assert mml.childNodes[0].childNodes[1].nodeName == 'mml:mi'
    assert mml.childNodes[0].childNodes[1].childNodes[0].nodeValue == '2'
    del mml

    mml = mp._print(Symbol("x__2"))
    assert mml.nodeName == 'ci'
    assert mml.childNodes[0].nodeName == 'mml:msup'
    assert mml.childNodes[0].childNodes[0].nodeName == 'mml:mi'
    assert mml.childNodes[0].childNodes[0].childNodes[0].nodeValue == 'x'
    assert mml.childNodes[0].childNodes[1].nodeName == 'mml:mi'
    assert mml.childNodes[0].childNodes[1].childNodes[0].nodeValue == '2'
    del mml

    mml = mp._print(Symbol("x_2"))
    assert mml.nodeName == 'ci'
    assert mml.childNodes[0].nodeName == 'mml:msub'
    assert mml.childNodes[0].childNodes[0].nodeName == 'mml:mi'
    assert mml.childNodes[0].childNodes[0].childNodes[0].nodeValue == 'x'
    assert mml.childNodes[0].childNodes[1].nodeName == 'mml:mi'
    assert mml.childNodes[0].childNodes[1].childNodes[0].nodeValue == '2'
    del mml

    mml = mp._print(Symbol("x^3_2"))
    assert mml.nodeName == 'ci'
    assert mml.childNodes[0].nodeName == 'mml:msubsup'
    assert mml.childNodes[0].childNodes[0].nodeName == 'mml:mi'
    assert mml.childNodes[0].childNodes[0].childNodes[0].nodeValue == 'x'
    assert mml.childNodes[0].childNodes[1].nodeName == 'mml:mi'
    assert mml.childNodes[0].childNodes[1].childNodes[0].nodeValue == '2'
    assert mml.childNodes[0].childNodes[2].nodeName == 'mml:mi'
    assert mml.childNodes[0].childNodes[2].childNodes[0].nodeValue == '3'
    del mml

    mml = mp._print(Symbol("x__3_2"))
    assert mml.nodeName == 'ci'
    assert mml.childNodes[0].nodeName == 'mml:msubsup'
    assert mml.childNodes[0].childNodes[0].nodeName == 'mml:mi'
    assert mml.childNodes[0].childNodes[0].childNodes[0].nodeValue == 'x'
    assert mml.childNodes[0].childNodes[1].nodeName == 'mml:mi'
    assert mml.childNodes[0].childNodes[1].childNodes[0].nodeValue == '2'
    assert mml.childNodes[0].childNodes[2].nodeName == 'mml:mi'
    assert mml.childNodes[0].childNodes[2].childNodes[0].nodeValue == '3'
    del mml

    mml = mp._print(Symbol("x_2_a"))
    assert mml.nodeName == 'ci'
    assert mml.childNodes[0].nodeName == 'mml:msub'
    assert mml.childNodes[0].childNodes[0].nodeName == 'mml:mi'
    assert mml.childNodes[0].childNodes[0].childNodes[0].nodeValue == 'x'
    assert mml.childNodes[0].childNodes[1].nodeName == 'mml:mrow'
    assert mml.childNodes[0].childNodes[1].childNodes[0].nodeName == 'mml:mi'
    assert mml.childNodes[0].childNodes[1].childNodes[0].childNodes[0].nodeValue == '2'
    assert mml.childNodes[0].childNodes[1].childNodes[1].nodeName == 'mml:mo'
    assert mml.childNodes[0].childNodes[1].childNodes[1].childNodes[0].nodeValue == ' '
    assert mml.childNodes[0].childNodes[1].childNodes[2].nodeName == 'mml:mi'
    assert mml.childNodes[0].childNodes[1].childNodes[2].childNodes[0].nodeValue == 'a'
    del mml

    mml = mp._print(Symbol("x^2^a"))
    assert mml.nodeName == 'ci'
    assert mml.childNodes[0].nodeName == 'mml:msup'
    assert mml.childNodes[0].childNodes[0].nodeName == 'mml:mi'
    assert mml.childNodes[0].childNodes[0].childNodes[0].nodeValue == 'x'
    assert mml.childNodes[0].childNodes[1].nodeName == 'mml:mrow'
    assert mml.childNodes[0].childNodes[1].childNodes[0].nodeName == 'mml:mi'
    assert mml.childNodes[0].childNodes[1].childNodes[0].childNodes[0].nodeValue == '2'
    assert mml.childNodes[0].childNodes[1].childNodes[1].nodeName == 'mml:mo'
    assert mml.childNodes[0].childNodes[1].childNodes[1].childNodes[0].nodeValue == ' '
    assert mml.childNodes[0].childNodes[1].childNodes[2].nodeName == 'mml:mi'
    assert mml.childNodes[0].childNodes[1].childNodes[2].childNodes[0].nodeValue == 'a'
    del mml

    mml = mp._print(Symbol("x__2__a"))
    assert mml.nodeName == 'ci'
    assert mml.childNodes[0].nodeName == 'mml:msup'
    assert mml.childNodes[0].childNodes[0].nodeName == 'mml:mi'
    assert mml.childNodes[0].childNodes[0].childNodes[0].nodeValue == 'x'
    assert mml.childNodes[0].childNodes[1].nodeName == 'mml:mrow'
    assert mml.childNodes[0].childNodes[1].childNodes[0].nodeName == 'mml:mi'
    assert mml.childNodes[0].childNodes[1].childNodes[0].childNodes[0].nodeValue == '2'
    assert mml.childNodes[0].childNodes[1].childNodes[1].nodeName == 'mml:mo'
    assert mml.childNodes[0].childNodes[1].childNodes[1].childNodes[0].nodeValue == ' '
    assert mml.childNodes[0].childNodes[1].childNodes[2].nodeName == 'mml:mi'
    assert mml.childNodes[0].childNodes[1].childNodes[2].childNodes[0].nodeValue == 'a'
    del mml
