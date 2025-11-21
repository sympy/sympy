from sympy.parsing.mathematica import parse_mathematica
from sympy.parsing.mathematica import parse_mathematica
from sympy.core.sympify import sympify
from sympy import sin, cos
from sympy.abc import x

#@profile
def test_mathematica():
    test_result = [0, 0, 0]
    d = {
        '- 6x': '-6*x',
        'Sin[x]^2': 'sin(x)**2',
        '2(x-1)': '2*(x-1)',
        '3y+8': '3*y+8',
        'ArcSin[2x+9(4-x)^2]/x': 'asin(2*x+9*(4-x)**2)/x',
        'x+y': 'x+y',
        '355/113': '355/113',
        '2.718281828': '2.718281828',
        'Cos(1/2 * π)': 'Cos(π/2)',
        'Sin[12]': 'sin(12)',
        'Exp[Log[4]]': 'exp(log(4))',
        '(x+1)(x+3)': '(x+1)*(x+3)',
        'Cos[ArcCos[3.6]]': 'cos(acos(3.6))',
        'Cos[x]==Sin[y]': 'Eq(cos(x), sin(y))',
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
        'Sin[x]^2Cos[y]^2': 'sin(x)**2*cos(y)**2',
        'Cos[x]^2(1 - Cos[y]^2)': 'cos(x)**2*(1-cos(y)**2)',
        'x y': 'x*y',
        'x  y': 'x*y',
        '2 x': '2*x',
        'x 8': 'x*8',
        '2 8': '2*8',
        '4.x': '4.*x',
        '4. 3': '4.*3',
        '4. 3.': '4.*3.',
        '1 2 3': '1*2*3',
        ' -  2 *  Sqrt[  2 3 *   ( 1   +  5 ) ]  ': '-2*sqrt(2*3*(1+5))',
        'Log[2,4]': 'log(4,2)',
        'Log[Log[2,4],4]': 'log(4,log(4,2))',
        'Exp[Sqrt[2]^2Log[2, 8]]': 'exp(sqrt(2)**2*log(8,2))',
        '-2^-10*x':'x/1024',                         # Test case from the issue 24150
        'a^-b*c':'c/a**b',                           # Test case from the issue 24150
        '-2^-3':'-1/8',                              # Test case from the issue 24150
        'd-Log[2*4^-3]':'d + log(32)',               # Test case from the issue 24150
        'Tan[Pi]^Sin[3*Pi/2]':'zoo',                 # Test case from the issue 24150
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
        'AiryBi[x]': 'airybi(x)',
        'AiryBiPrime[7]': 'airybiprime(7)',
        'LogIntegral[4]': ' li(4)',
        'PrimePi[7]': 'primepi(7)',
        'Prime[5]': 'prime(5)',
        'PrimeQ[5]': 'isprime(5)'
        }

    for e in d:
        test_result[0] = test_result[0] + 1
        if(parse_mathematica(e) == sympify(d[e])): test_result[1] = test_result[1] + 1
        else: test_result[2] = test_result[2] + 1

    # The parsed form of this expression should not evaluate the Lambda object:
    test_result[0] = test_result[0] + 1
    if(parse_mathematica("Sin[#]^2 + Cos[#]^2 &[x]") == sin(x)**2 + cos(x)**2): test_result[1] = test_result[1] + 1
    else: test_result[2] = test_result[2] + 1

    return test_result


result = test_mathematica()
print("Run test: " + str(result[0]))
print("Pass test: " + str(result[1]))
print("Error test: " + str(result[2]))


#def parse(self, s):
#        s2 = self._from_mathematica_to_tokens(s)
#        s3 = self._from_tokens_to_fullformlist(s2)
#        s4 = self._from_fullformlist_to_sympy(s3)
#        return s4

## MiB ==> A mebibyte is equal to 2^20 or 1,048,576 bytes.
## zoo == complex infinity
## Testy typu: Stub je objekt, ktorý obsahuje preddefinované údaje a používa ich na odpovedanie na hovory počas testov. Používa sa, keď nemôžeme alebo nechceme zapojiť objekty, ktoré by odpovedali skutočnými údajmi alebo mali nežiaduce vedľajšie účinky.