from sympy import sin, Function, symbols
from sympy.parsing.mathematica import mathematica, MathematicaParser
from sympy.core.sympify import sympify
from sympy.abc import n, w, x, y, z


def test_mathematica():
    d = {
        '- 6x': '-6*x',
        'Sin[x]^2': 'sin(x)**2',
        '2(x-1)': '2*(x-1)',
        '3y+8': '3*y+8',
        'ArcSin[2x+9(4-x)^2]/x': 'asin(2*x+9*(4-x)**2)/x',
        'x+y': 'x+y',
        '355/113': '355/113',
        '2.718281828': '2.718281828',
        'Sin[12]': 'sin(12)',
        'Exp[Log[4]]': 'exp(log(4))',
        '(x+1)(x+3)': '(x+1)*(x+3)',
        'Cos[ArcCos[3.6]]': 'cos(acos(3.6))',
        'Cos[x]==Sin[y]': 'cos(x)==sin(y)',
        '2*Sin[x+y]': '2*sin(x+y)',
        'Sin[x]+Cos[y]': 'sin(x)+cos(y)',
        'Sin[Cos[x]]': 'sin(cos(x))',
        '2*Sqrt[x+y]': '2*sqrt(x+y)',   # Test case from the issue 4259
        '+Sqrt[2]': 'sqrt(2)',
        '-Sqrt[2]': '-sqrt(2)',
        '-1/Sqrt[2]': '-1/sqrt(2)',
        '-(1/Sqrt[3])': '-(1/sqrt(3))',
        '1/(2*Sqrt[5])': '1/(2*sqrt(5))',
        'Mod[5,3]': 'Mod(5,3)',
        '-Mod[5,3]': '-Mod(5,3)',
        '(x+1)y': '(x+1)*y',
        'x(y+1)': 'x*(y+1)',
        'Sin[x]Cos[y]': 'sin(x)*cos(y)',
        'Sin[x]**2Cos[y]**2': 'sin(x)**2*cos(y)**2',
        'Cos[x]^2(1 - Cos[y]^2)': 'cos(x)**2*(1-cos(y)**2)',
        'x y': 'x*y',
        '2 x': '2*x',
        'x 8': 'x*8',
        '2 8': '2*8',
        '1 2 3': '1*2*3',
        ' -  2 *  Sqrt[  2 3 *   ( 1   +  5 ) ]  ': '-2*sqrt(2*3*(1+5))',
        'Log[2,4]': 'log(4,2)',
        'Log[Log[2,4],4]': 'log(4,log(4,2))',
        'Exp[Sqrt[2]^2Log[2, 8]]': 'exp(sqrt(2)**2*log(8,2))',
        'ArcSin[Cos[0]]': 'asin(cos(0))',
        'Log2[16]': 'log(16,2)',
        'Max[1,-2,3,-4]': 'Max(1,-2,3,-4)',
        'Min[1,-2,3]': 'Min(1,-2,3)',
        'Exp[I Pi/2]': 'exp(I*pi/2)',
        'ArcTan[x,y]': 'atan2(y,x)',
        'Pochhammer[x,y]': 'rf(x,y)',
        'ExpIntegralEi[x]': 'Ei(x)',
        'SinIntegral[x]': 'Si(x)',
        'CosIntegral[x]': 'Ci(x)',
        'AiryAi[x]': 'airyai(x)',
        'AiryAiPrime[5]': 'airyaiprime(5)',
        'AiryBi[x]' :'airybi(x)',
        'AiryBiPrime[7]' :'airybiprime(7)',
        'LogIntegral[4]':' li(4)',
        'PrimePi[7]': 'primepi(7)',
        'Prime[5]': 'prime(5)',
        'PrimeQ[5]': 'isprime(5)'
        }

    for e in d:
        assert mathematica(e) == sympify(d[e])


def test_parser_mathematica_tokenizer():
    parser = MathematicaParser()

    chain = lambda expr: parser._parse_tokenized_code(parser._tokenize_mathematica_code(expr))

    assert chain("a + b*c") == ["Plus", "a", ["Times", "b", "c"]]
    assert chain("a + b* c* d + 2 * e") == ["Plus", "a", ["Times", "b", "c", "d"], ["Times", "2", "e"]]
    assert chain("-a") == ["Times", -1, "a"]
    assert chain("a - b") == ["Plus", "a", ["Times", -1, "b"]]
    assert chain("a / b") == ["Times", "a", ["Power", "b", -1]]

    # Parentheses of various kinds, i.e. ( )  [ ]  [[ ]]  { }
    assert chain("(a + b) + c") == ["Plus", ["Plus", "a", "b"], "c"]
    assert chain(" a + (b + c) + d ") == ["Plus", "a", ["Plus", "b", "c"], "d"]
    assert chain("a * (b + c)") == ["Times", "a", ["Plus", "b", "c"]]
    assert chain("{a, b, 2, c}") == ["List", "a", "b", "2", "c"]
    assert chain("{a, {b, c}}") == ["List", "a", ["List", "b", "c"]]
    assert chain("a[b, c]") == ["a", "b", "c"]
    assert chain("a[[b, c]]") == ["Part", "a", "b", "c"]
    assert chain("a[[b, c[[d, {e,f}]]]]") == ["Part", "a", "b", ["Part", "c", "d", ["List", "e", "f"]]]
    assert chain("a[b[[c,d]]]") == ["a", ["Part", "b", "c", "d"]]

    # Flat operator:
    assert chain("a*b*c*d*e") == ["Times", "a", "b", "c", "d", "e"]
    assert chain("a +b + c+ d+e") == ["Plus", "a", "b", "c", "d", "e"]

    # Right priority operator:
    assert chain("a^b") == ["Power", "a", "b"]
    assert chain("a^b^c") == ["Power", "a", ["Power", "b", "c"]]
    assert chain("a^b^c^d") == ["Power", "a", ["Power", "b", ["Power", "c", "d"]]]

    # Left priority operator:
    assert chain("a/.b") == ["ReplaceAll", "a", "b"]
    assert chain("a/.b/.c/.d") == ["ReplaceAll", ["ReplaceAll", ["ReplaceAll", "a", "b"], "c"], "d"]


def test_parser_mathematica_exp_alt():
    parser = MathematicaParser()

    convert_chain2 = lambda expr: parser._convert_pylist_to_sympymform(parser._convert_fullform_to_pylist(expr))
    convert_chain3 = lambda expr: parser._convert_sympymform_to_sympy(convert_chain2(expr))

    Sin, Times, Plus, Power = symbols("Sin Times Plus Power", cls=Function)

    full_form1 = "Sin[Times[x, y]]"
    full_form2 = "Plus[Times[x, y], z]"
    full_form3 = "Sin[Times[x, Plus[y, z], Power[w, n]]]]"

    assert parser._convert_fullform_to_pylist(full_form1) == ["Sin", ["Times", "x", "y"]]
    assert parser._convert_fullform_to_pylist(full_form2) == ["Plus", ["Times", "x", "y"], "z"]
    assert parser._convert_fullform_to_pylist(full_form3) == ["Sin", ["Times", "x", ["Plus", "y", "z"], ["Power", "w", "n"]]]

    assert convert_chain2(full_form1) == Sin(Times(x, y))
    assert convert_chain2(full_form2) == Plus(Times(x, y), z)
    assert convert_chain2(full_form3) == Sin(Times(x, Plus(y, z), Power(w, n)))

    assert convert_chain3(full_form1) == sin(x*y)
    assert convert_chain3(full_form2) == x*y + z
    assert convert_chain3(full_form3) == sin(x*(y + z)*w**n)
