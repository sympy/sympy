from __future__ import annotations
from sympy import sin, Function, symbols, Dummy, Lambda, cos, Symbol, factorial, S
from sympy.parsing.mathematica import (
    parse_mathematica, parse_mathematica_to_fullformlist, MathematicaParser)
from sympy.core.sympify import sympify
from sympy.abc import n, w, x, y, z
from sympy.testing.pytest import raises
from sympy.logic.boolalg import And, Or, Not


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
        'Cos(1/2 * π)': 'Cos*π/2',
        'Cos[1/2 * π]': 'cos(π/2)',
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
        'PrimeQ[5]': 'isprime(5)',
        'Rational[2,19]': 'Rational(2,19)',    # test case for issue 25716
        'Pi' : 'pi',  # test cases from issue 27868
        '3*Pi' : '3*pi',
        'Ωπ' : 'Ωπ',
        '3*Ωπ' : '3*Ωπ',
        '3 Ω π' : '3*Ω*π',
        'Pi*Ω' : 'pi*Ω',
        'Sqrt[2]*σ' : 'sqrt(2)*σ',
        'Log[e^2]' : 'log(e**2)',
        'Log[E^2]' : '2',
        'Log[ExponentialE^2]' : '2',
        '(3*數學)/數' : '(3*數學)/數',
        'I^2' : '-1',
        'ImaginaryI^2' : '-1',
        'ImaginaryJ^2' : '-1',
        '\\[Alpha]': 'α',
        'x\\[Beta]y': 'xβy',
        'x \\[Beta] y': 'x*β*y',
        'a + b\\[Gamma]\\[CapitalGamma]d': 'a + bγΓd',
        'a + b \\[Gamma] \\[CapitalGamma] d': 'a + b*γ*Γ*d',
        'a\\[LongEqual]b': 'Eq(a,b)',
        }

    for e in d:
        assert parse_mathematica(e) == sympify(d[e])

    # The parsed form of this expression should not evaluate the Lambda object:
    assert parse_mathematica("Sin[#]^2 + Cos[#]^2 &[x]") == sin(x)**2 + cos(x)**2

    d1, d2, d3 = symbols("d1:4", cls=Dummy)
    assert parse_mathematica("Sin[#] + Cos[#3] &").dummy_eq(Lambda((d1, d2, d3), sin(d1) + cos(d3)))
    assert parse_mathematica("Sin[#^2] &").dummy_eq(Lambda(d1, sin(d1**2)))
    assert parse_mathematica("Function[x, x^3]") == Lambda(x, x**3)
    assert parse_mathematica("Function[{x, y}, x^2 + y^2]") == Lambda((x, y), x**2 + y**2)


def test_parser_mathematica_tokenizer():
    # ``parse_mathematica_to_fullformlist`` is the public entry point for the
    # intermediate FullForm (as nested Python lists), i.e. the parse pipeline up
    # to (but excluding) the conversion into SymPy objects.
    chain = parse_mathematica_to_fullformlist

    # Basic patterns
    assert chain("x") == "x"
    assert chain("42") == "42"
    assert chain(".2") == ".2"
    # Unary plus is kept: Mathematica reads ``+x`` as ``Plus[x]``.
    assert chain("+x") == ["Plus", "x"]
    assert chain("-1") == "-1"
    assert chain("- 3") == "-3"
    assert chain("α") == "α"
    assert chain("α + β") == ["Plus", "α", "β"]
    assert chain("αβγ") == "αβγ"
    assert chain("α̇𝔟⃗𝒞̂") == "α̇𝔟⃗𝒞̂"
    assert chain("α β γ") == ["Times", "α", "β", "γ"]
    assert chain("μ1ν2") == "μ1ν2"
    assert chain("μ1 ν2") == ["Times", "μ1", "ν2"]
    assert chain("α + βγ") == ["Plus", "α", "βγ"]
    assert chain("α + β γ") == ["Plus", "α", ["Times", "β", "γ"]]
    assert chain("α̇ + 𝔟⃗ 𝒞̂") == ["Plus", "α̇", ["Times", "𝔟⃗", "𝒞̂"]]
    assert chain("+Sin[x]") == ["Plus", ["Sin", "x"]]
    assert chain("-Sin[x]") == ["Times", "-1", ["Sin", "x"]]
    assert chain("Cos(1/2 * π)") == ["Times", "Cos", ["Times", "1", ["Power", "2", "-1"], "π"]]
    assert chain("Cos[1/2 * π]") == ["Cos", ["Times", "1", ["Power", "2", "-1"], "π"]]
    assert chain("Cos[x]==Sin[y]") == ["Equal", ["Cos", "x"], ["Sin", "y"]]
    assert chain("Cos[x]!=Sin[y]") == ["Unequal", ["Cos", "x"], ["Sin", "y"]]
    assert chain("x(a+1)") == ["Times", "x", ["Plus", "a", "1"]]
    assert chain("(x)") == "x"
    assert chain("(+x)") == ["Plus", "x"]
    assert chain("-a") == ["Times", "-1", "a"]
    assert chain("(-x)") == ["Times", "-1", "x"]
    assert chain("(x + y)") == ["Plus", "x", "y"]
    assert chain("3 + 4") == ["Plus", "3", "4"]
    assert chain("a - 3") == ["Plus", "a", "-3"]
    assert chain("a - b") == ["Plus", "a", ["Times", "-1", "b"]]
    assert chain("7 * 8") == ["Times", "7", "8"]
    assert chain("a + b*c") == ["Plus", "a", ["Times", "b", "c"]]
    assert chain("a + b* c* d + 2 * e") == ["Plus", "a", ["Times", "b", "c", "d"], ["Times", "2", "e"]]
    assert chain("a / b") == ["Times", "a", ["Power", "b", "-1"]]

    # Missing asterisk (*) patterns:
    assert chain("x y") == ["Times", "x", "y"]
    assert chain("3 4") == ["Times", "3", "4"]
    assert chain("a[b] c") == ["Times", ["a", "b"], "c"]
    assert chain("(x) (y)") == ["Times", "x", "y"]
    assert chain("3 (a)") == ["Times", "3", "a"]
    assert chain("(a) b") == ["Times", "a", "b"]
    assert chain("4.2") == "4.2"
    assert chain("4 2") == ["Times", "4", "2"]
    assert chain("4  2") == ["Times", "4", "2"]
    assert chain("3 . 4") == ["Dot", "3", "4"]
    assert chain("4. 2") == ["Times", "4.", "2"]
    assert chain("x.y") == ["Dot", "x", "y"]
    assert chain("4.y") == ["Times", "4.", "y"]
    assert chain("4 .y") == ["Dot", "4", "y"]
    assert chain("x.4") == ["Times", "x", ".4"]
    assert chain("x0.3") == ["Times", "x0", ".3"]
    assert chain("x. 4") == ["Dot", "x", "4"]

    # Comments
    assert chain("a (* +b *) + c") == ["Plus", "a", "c"]
    assert chain("a (* + b *) + (**)c (* +d *) + e") == ["Plus", "a", "c", "e"]
    assert chain("""a + (*
    + b
    *) c + (* d
    *) e
    """) == ["Plus", "a", "c", "e"]

    # Operators couples + and -, * and / are mutually associative:
    # (i.e. expression gets flattened when mixing these operators)
    assert chain("a*b/c") == ["Times", "a", "b", ["Power", "c", "-1"]]
    assert chain("a/b*c") == ["Times", "a", ["Power", "b", "-1"], "c"]
    assert chain("a+b-c") == ["Plus", "a", "b", ["Times", "-1", "c"]]
    assert chain("a-b+c") == ["Plus", "a", ["Times", "-1", "b"], "c"]
    assert chain("-a + b -c ") == ["Plus", ["Times", "-1", "a"], "b", ["Times", "-1", "c"]]
    assert chain("a/b/c*d") == ["Times", "a", ["Power", "b", "-1"], ["Power", "c", "-1"], "d"]
    assert chain("a/b/c") == ["Times", "a", ["Power", "b", "-1"], ["Power", "c", "-1"]]
    assert chain("a-b-c") == ["Plus", "a", ["Times", "-1", "b"], ["Times", "-1", "c"]]
    assert chain("1/a") == ["Times", "1", ["Power", "a", "-1"]]
    assert chain("1/a/b") == ["Times", "1", ["Power", "a", "-1"], ["Power", "b", "-1"]]
    assert chain("-1/a*b") == ["Times", "-1", ["Power", "a", "-1"], "b"]

    # Enclosures of various kinds, i.e. ( )  [ ]  [[ ]]  { }
    assert chain("(a + b) + c") == ["Plus", ["Plus", "a", "b"], "c"]
    assert chain(" a + (b + c) + d ") == ["Plus", "a", ["Plus", "b", "c"], "d"]
    assert chain("a * (b + c)") == ["Times", "a", ["Plus", "b", "c"]]
    assert chain("a b (c d)") == ["Times", "a", "b", ["Times", "c", "d"]]
    assert chain("{a, b, 2, c}") == ["List", "a", "b", "2", "c"]
    assert chain("{a, {b, c}}") == ["List", "a", ["List", "b", "c"]]
    assert chain("{{a}}") == ["List", ["List", "a"]]
    assert chain("a[b, c]") == ["a", "b", "c"]
    assert chain("a[[b, c]]") == ["Part", "a", "b", "c"]
    assert chain("a[b[c]]") == ["a", ["b", "c"]]
    assert chain("a[[b, c[[d, {e,f}]]]]") == ["Part", "a", "b", ["Part", "c", "d", ["List", "e", "f"]]]
    assert chain("a[b[[c,d]]]") == ["a", ["Part", "b", "c", "d"]]
    assert chain("a[[b[c]]]") == ["Part", "a", ["b", "c"]]
    assert chain("a[[b[[c]]]]") == ["Part", "a", ["Part", "b", "c"]]
    assert chain("a[[b[c[[d]]]]]") == ["Part", "a", ["b", ["Part", "c", "d"]]]
    assert chain("a[b[[c[d]]]]") == ["a", ["Part", "b", ["c", "d"]]]
    assert chain("x[[a+1, b+2, c+3]]") == ["Part", "x", ["Plus", "a", "1"], ["Plus", "b", "2"], ["Plus", "c", "3"]]
    assert chain("x[a+1, b+2, c+3]") == ["x", ["Plus", "a", "1"], ["Plus", "b", "2"], ["Plus", "c", "3"]]
    assert chain("{a+1, b+2, c+3}") == ["List", ["Plus", "a", "1"], ["Plus", "b", "2"], ["Plus", "c", "3"]]

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

    # ``//`` is postfix application: ``x // f`` is ``f[x]`` (left associative).
    assert chain("a//b") == ["b", "a"]
    assert chain("a//b//c") == ["c", ["b", "a"]]
    assert chain("a//b//c//d") == ["d", ["c", ["b", "a"]]]

    # Not operator
    assert chain("!x") == ["Not", "x"]
    assert chain("!(a + b)") == ["Not", ["Plus", "a", "b"]]
    assert chain("!True") == ["Not", "True"]
    assert chain("!False") == ["Not", "False"]

    # Distinguish Not and factorial
    assert chain("x!") == ["Factorial", "x"]
    assert chain("(a + b)!") == ["Factorial", ["Plus", "a", "b"]]

    # Combination of Not with And/Or
    assert chain("!x && y") == ["And", ["Not", "x"], "y"]
    assert chain("x || !y") == ["Or", "x", ["Not", "y"]]
    assert chain("!x || !y") == ["Or", ["Not", "x"], ["Not", "y"]]

    # Enablement of implicit multiplication with factorial
    assert chain("x!y") == ["Times", ["Factorial", "x"], "y"]
    assert chain("x!y!") == ["Times", ["Factorial", "x"], ["Factorial", "y"]]

    # Compound expressions
    assert chain("a;b") == ["CompoundExpression", "a", "b"]
    assert chain("a;") == ["CompoundExpression", "a", "Null"]
    assert chain("a;b;") == ["CompoundExpression", "a", "b", "Null"]
    assert chain("a[b;c]") == ["a", ["CompoundExpression", "b", "c"]]
    assert chain("a[b,c;d,e]") == ["a", "b", ["CompoundExpression", "c", "d"], "e"]
    assert chain("a[b,c;,d]") == ["a", "b", ["CompoundExpression", "c", "Null"], "d"]

    # New lines
    assert chain("a\nb\n") == ["CompoundExpression", "a", "b"]
    assert chain("a\n\nb\n (c \nd)  \n") == ["CompoundExpression", "a", "b", ["Times", "c", "d"]]
    assert chain("\na; b\nc") == ["CompoundExpression", "a", "b", "c"]
    assert chain("a + \nb\n") == ["Plus", "a", "b"]
    assert chain("a\nb; c; d\n e; (f \n g); h + \n i") == ["CompoundExpression", "a", "b", "c", "d", "e", ["Times", "f", "g"], ["Plus", "h", "i"]]
    assert chain("\n{\na\nb; c; d\n e (f \n g); h + \n i\n\n}\n") == ["List", ["CompoundExpression", ["Times", "a", "b"], "c", ["Times", "d", "e", ["Times", "f", "g"]], ["Plus", "h", "i"]]]

    # Patterns
    assert chain("y_") == ["Pattern", "y", ["Blank"]]
    assert chain("y_.") == ["Optional", ["Pattern", "y", ["Blank"]]]
    assert chain("y__") == ["Pattern", "y", ["BlankSequence"]]
    assert chain("y___") == ["Pattern", "y", ["BlankNullSequence"]]
    assert chain("a[b_.,c_]") == ["a", ["Optional", ["Pattern", "b", ["Blank"]]], ["Pattern", "c", ["Blank"]]]
    assert chain("b_. c") == ["Times", ["Optional", ["Pattern", "b", ["Blank"]]], "c"]
    # Patterns in head position, i.e. ``Pattern[F, Blank[]][Pattern[x, Blank[]]]``
    assert chain("F_[x_]") == [["Pattern", "F", ["Blank"]], ["Pattern", "x", ["Blank"]]]
    assert chain("F_[x_, y_[z_]]") == [
        ["Pattern", "F", ["Blank"]],
        ["Pattern", "x", ["Blank"]],
        [["Pattern", "y", ["Blank"]], ["Pattern", "z", ["Blank"]]],
    ]

    # Slots for lambda functions.  A positional slot (``#3``) carries an
    # integer index; a named slot (``#name``) carries a string, which in the
    # FullForm list is ``["_Str", "name"]`` (matching Mathematica's Slot["name"]).
    assert chain("#") == ["Slot", "1"]
    assert chain("#3") == ["Slot", "3"]
    assert chain("#n") == ["Slot", ["_Str", "n"]]
    assert chain("#name") == ["Slot", ["_Str", "name"]]
    assert chain("##") == ["SlotSequence", "1"]
    assert chain("##3") == ["SlotSequence", "3"]
    # ``SlotSequence`` only takes an integer index, so ``##name`` is
    # ``SlotSequence[1] * name`` (implicit multiplication), not a named slot.
    assert chain("##a") == ["Times", ["SlotSequence", "1"], "a"]
    assert chain("##name") == ["Times", ["SlotSequence", "1"], "name"]

    # Lambda functions
    assert chain("x&") == ["Function", "x"]
    assert chain("#&") == ["Function", ["Slot", "1"]]
    assert chain("#+3&") == ["Function", ["Plus", ["Slot", "1"], "3"]]
    assert chain("#1 + #2&") == ["Function", ["Plus", ["Slot", "1"], ["Slot", "2"]]]
    assert chain("# + #&") == ["Function", ["Plus", ["Slot", "1"], ["Slot", "1"]]]
    assert chain("#&[x]") == [["Function", ["Slot", "1"]], "x"]
    assert chain("#1 + #2 & [x, y]") == [["Function", ["Plus", ["Slot", "1"], ["Slot", "2"]]], "x", "y"]
    assert chain("#1^2#2^3&") == ["Function", ["Times", ["Power", ["Slot", "1"], "2"], ["Power", ["Slot", "2"], "3"]]]

    # Assignment operators
    assert chain("a=b") == ["Set", "a", "b"]
    assert chain("a:=b") == ["SetDelayed", "a", "b"]
    assert chain("a+=b") == ["AddTo", "a", "b"]
    assert chain("a-=b") == ["SubtractFrom", "a", "b"]
    assert chain("a*=b") == ["TimesBy", "a", "b"]
    assert chain("a/=b") == ["DivideBy", "a", "b"]

    # Rules, replacement and conditions
    assert chain("a->b") == ["Rule", "a", "b"]
    assert chain("a:>b") == ["RuleDelayed", "a", "b"]
    assert chain("a/;b") == ["Condition", "a", "b"]

    # Comparison operators
    assert chain("a>b") == ["Greater", "a", "b"]
    assert chain("a>=b") == ["GreaterEqual", "a", "b"]
    assert chain("a<b") == ["Less", "a", "b"]
    assert chain("a<=b") == ["LessEqual", "a", "b"]
    assert chain("a===b") == ["SameQ", "a", "b"]
    assert chain("a=!=b") == ["UnsameQ", "a", "b"]

    # Function application operators.  ``@`` is prefix application and ``//``
    # is postfix; ``@@``/``@@@`` are Apply, ``/@``/``//@`` are Map/MapAll.
    assert chain("f@x") == ["f", "x"]
    assert chain("a@b@c") == ["a", ["b", "c"]]         # right associative
    assert chain("f@@x") == ["Apply", "f", "x"]
    assert chain("f@@@x") == ["Apply", "f", "x", ["List", "1"]]
    assert chain("f/@x") == ["Map", "f", "x"]
    assert chain("f//@x") == ["MapAll", "f", "x"]

    # Derivative: ``f'`` is ``Derivative[1][f]``; primes accumulate the order.
    assert chain("f'") == [["Derivative", "1"], "f"]
    assert chain("f''") == [["Derivative", "2"], "f"]
    assert chain("x'[t]") == [[["Derivative", "1"], "x"], "t"]
    # The prime binds to the whole pattern, not to the symbol inside it:
    # ``f_'[x_]`` is ``Derivative[1][Pattern[f, Blank[]]][Pattern[x, Blank[]]]``.
    assert chain("f_'[x_]") == [[["Derivative", "1"], ["Pattern", "f", ["Blank"]]],
                                ["Pattern", "x", ["Blank"]]]
    assert chain("f_''[x_]") == [[["Derivative", "2"], ["Pattern", "f", ["Blank"]]],
                                 ["Pattern", "x", ["Blank"]]]
    # The prime binds tighter than ``[``, so ``f'[x]`` is ``Derivative[1][f][x]``
    # and not ``Derivative[1][f[x]]``.  This has to keep holding when ``f'[x]``
    # is a subexpression rather than the whole input:
    d1fx = [[["Derivative", "1"], "f"], "x"]
    assert chain("f'[x]") == d1fx
    assert chain("f'[x]*g[x]") == ["Times", d1fx, ["g", "x"]]
    assert chain("f'[x]g[x]") == ["Times", d1fx, ["g", "x"]]
    assert chain("2 f'[x]") == ["Times", "2", d1fx]
    assert chain("f'[x]^2") == ["Power", d1fx, "2"]
    assert chain("f'[x][y]") == [d1fx, "y"]
    assert chain("f'[g'[x]]") == [[["Derivative", "1"], "f"],
                                  [[["Derivative", "1"], "g"], "x"]]
    # ... and equally for a run of several primes:
    d2fx = [[["Derivative", "2"], "f"], "x"]
    d3fx = [[["Derivative", "3"], "f"], "x"]
    assert chain("f''[x]") == d2fx
    assert chain("f'''[x]") == d3fx
    assert chain("f''[x]*g[x]") == ["Times", d2fx, ["g", "x"]]
    assert chain("2 f''[x]") == ["Times", "2", d2fx]
    assert chain("f''[x]^2") == ["Power", d2fx, "2"]
    assert chain("f''[x][y]") == [d2fx, "y"]
    assert chain("f'[x]+f''[x]+f'''[x]") == ["Plus", d1fx, d2fx, d3fx]
    assert chain("f''[x]/g'''[y]") == [
        "Times", d2fx, ["Power", [[["Derivative", "3"], "g"], "y"], "-1"]]
    # Conversely, a prime *after* a closed bracket takes the whole bracketed
    # expression as its operand: ``f[x]'`` is ``Derivative[1][f[x]]``.
    assert chain("f[x]'") == [["Derivative", "1"], ["f", "x"]]
    assert chain("f[x]'[y]") == [[["Derivative", "1"], ["f", "x"]], "y"]
    assert chain("f'[[1]]") == ["Part", [["Derivative", "1"], "f"], "1"]
    assert chain("(f+g)'[x]") == [[["Derivative", "1"], ["Plus", "f", "g"]], "x"]
    # Same deferral when the prime follows another still unapplied postfix op.
    assert chain("f!'") == [["Derivative", "1"], ["Factorial", "f"]]
    assert chain("f'!") == ["Factorial", [["Derivative", "1"], "f"]]

    # Postfix operators and Alternatives / PatternTest / Repeated
    assert chain("a--") == ["Decrement", "a"]
    assert chain("a!!") == ["Factorial2", "a"]
    assert chain("a|b|c") == ["Alternatives", "a", "b", "c"]
    assert chain("a?b") == ["PatternTest", "a", "b"]
    assert chain("a..") == ["Repeated", "a"]
    assert chain("a...") == ["RepeatedNull", "a"]

    # Span (``;;``) flattens into a single node with up to three arguments
    assert chain("a;;b") == ["Span", "a", "b"]
    assert chain("a;;b;;c") == ["Span", "a", "b", "c"]
    assert chain("a[[b;;c]]") == ["Part", "a", ["Span", "b", "c"]]

    # Typed blank via the infix ``_`` (e.g. ``x_Integer``)
    assert chain("x_Integer") == ["Pattern", "x", ["Blank", "Integer"]]

    # Strings inside Mathematica expressions:
    assert chain('"abc"') == ["_Str", "abc"]
    assert chain('"a\\"b"') == ["_Str", 'a"b']
    # This expression does not make sense mathematically, it's just testing the parser:
    assert chain('x + "abc" ^ 3') == ["Plus", "x", ["Power", ["_Str", "abc"], "3"]]
    assert chain('"a (* b *) c"') == ["_Str", "a (* b *) c"]
    assert chain('"a" (* b *) ') == ["_Str", "a"]
    assert chain('"a [ b] "') == ["_Str", "a [ b] "]
    raises(SyntaxError, lambda: chain('"'))
    raises(SyntaxError, lambda: chain('"\\"'))
    raises(SyntaxError, lambda: chain('"abc'))
    raises(SyntaxError, lambda: chain('"abc\\"def'))

    # Invalid expressions:
    raises(SyntaxError, lambda: chain("(,"))
    raises(SyntaxError, lambda: chain("()"))
    raises(SyntaxError, lambda: chain("a (* b"))
    raises(SyntaxError, lambda: parse_mathematica(r"\[Gamma[y"))
    raises(SyntaxError, lambda: parse_mathematica(r"\[Gamma\[Alpha]]"))
    raises(SyntaxError, lambda: parse_mathematica(r"]Gamma\["))
    raises(IndexError, lambda: parse_mathematica(r"]c"))
    raises(IndexError, lambda: parse_mathematica(r"]\[Gamma]"))


def test_parser_mathematica_exp_alt():
    parser = MathematicaParser()

    convert_chain2 = lambda expr: parser._from_fullformlist_to_fullformsympy(parser._from_fullform_to_fullformlist(expr))
    convert_chain3 = lambda expr: parser._from_fullformsympy_to_sympy(convert_chain2(expr))

    Sin, Times, Plus, Power = symbols("Sin Times Plus Power", cls=Function)

    full_form1 = "Sin[Times[x, y]]"
    full_form2 = "Plus[Times[x, y], z]"
    full_form3 = "Sin[Times[x, Plus[y, z], Power[w, n]]]]"
    full_form4 = "Rational[Rational[x, y], z]"

    assert parser._from_fullform_to_fullformlist(full_form1) == ["Sin", ["Times", "x", "y"]]
    assert parser._from_fullform_to_fullformlist(full_form2) == ["Plus", ["Times", "x", "y"], "z"]
    assert parser._from_fullform_to_fullformlist(full_form3) == ["Sin", ["Times", "x", ["Plus", "y", "z"], ["Power", "w", "n"]]]
    assert parser._from_fullform_to_fullformlist(full_form4) == ["Rational", ["Rational", "x", "y"], "z"]

    assert convert_chain2(full_form1) == Sin(Times(x, y))
    assert convert_chain2(full_form2) == Plus(Times(x, y), z)
    assert convert_chain2(full_form3) == Sin(Times(x, Plus(y, z), Power(w, n)))

    assert convert_chain3(full_form1) == sin(x*y)
    assert convert_chain3(full_form2) == x*y + z
    assert convert_chain3(full_form3) == sin(x*(y + z)*w**n)


def test_Mathematica_literal_regex():
    import sys
    import re
    from sympy.parsing.mathematica import MathematicaParser
    literal_regex = re.compile(MathematicaParser._literal)

    for c in map(chr, range(sys.maxunicode+1)):
        if c == "_" or (not c.isidentifier() and not f"x{c}".isidentifier()):
            assert not literal_regex.match(c)
            assert not literal_regex.fullmatch(f"x{c}")
        else:
            if c.isidentifier():
                assert literal_regex.fullmatch(c)
            if f"x{c}".isidentifier():
                assert literal_regex.fullmatch(f"x{c}")


def test_mathematica_not_operator():
    # Basic tests
    x = Symbol('x')
    assert parse_mathematica("!x") == Not(x)

    # And / Or combinations
    x1, x2 = Symbol('x1'), Symbol('x2')
    assert parse_mathematica("x1 && !x2") == And(x1, Not(x2))
    assert parse_mathematica("!x1 || !x2") == Or(Not(x1), Not(x2))

    # Constants
    assert parse_mathematica("!True") == S.false
    assert parse_mathematica("!False") == S.true

    # Factorial distinction
    assert parse_mathematica("x!") == factorial(x)

    # Factorial with implicit multiplication
    assert parse_mathematica("x!y") == y*factorial(x)
    assert parse_mathematica("x!y!") == factorial(x)*factorial(y)


def test_mathematica_star_operator():
    # \[Star], the literal ⋆ (U+22C6) and "Star[...]" are the same operator,
    # which is distinct from Times.  Assertions are made on the FullForm list so
    # they are not affected by any SymPy-level reordering.
    fullform = parse_mathematica_to_fullformlist
    assert fullform(r"a \[Star] b") == ["Star", "a", "b"]
    assert fullform("a \N{STAR OPERATOR} b") == ["Star", "a", "b"]
    assert fullform("Star[a, b]") == ["Star", "a", "b"]

    # Star is flat: a ⋆ b ⋆ c -> Star[a, b, c]
    assert fullform("a \N{STAR OPERATOR} b \N{STAR OPERATOR} c") == ["Star", "a", "b", "c"]

    # Precedence: Star binds looser than Times but tighter than Plus.
    assert fullform("a + b \N{STAR OPERATOR} c") == ["Plus", "a", ["Star", "b", "c"]]
    assert fullform("a*b \N{STAR OPERATOR} c*d") == \
        ["Star", ["Times", "a", "b"], ["Times", "c", "d"]]

    # Times is unaffected.
    assert fullform("a * b") == ["Times", "a", "b"]

    # End-to-end: with no built-in meaning Star becomes an undefined Function.
    a, b = symbols("a b")
    assert parse_mathematica(r"a \[Star] b") == Function("Star")(a, b)


def test_mathematica_operators_without_builtin_meaning():
    # Operators listed under "Operators without built-in meanings" in
    # https://reference.wolfram.com/language/tutorial/TextualInputAndOutput.html
    # Each must produce its own head instead of being dropped/parsed as Times.
    # Assertions are on the FullForm list (see parse_mathematica_to_fullformlist).
    fullform = parse_mathematica_to_fullformlist

    # Infix operators, both named-character and literal forms. CirclePlus,
    # CircleTimes, TildeTilde and LeftRightArrow are flat (x op y op z ->
    # op[x, y, z]).
    assert fullform(r"x \[CirclePlus] y") == ["CirclePlus", "x", "y"]
    assert fullform("x \N{CIRCLED PLUS} y \N{CIRCLED PLUS} z") == \
        ["CirclePlus", "x", "y", "z"]
    assert fullform(r"x \[CircleTimes] y") == ["CircleTimes", "x", "y"]
    assert fullform("x \N{CIRCLED TIMES} y \N{CIRCLED TIMES} z") == \
        ["CircleTimes", "x", "y", "z"]
    assert fullform(r"x \[TildeTilde] y") == ["TildeTilde", "x", "y"]
    assert fullform("x \N{ALMOST EQUAL TO} y \N{ALMOST EQUAL TO} z") == \
        ["TildeTilde", "x", "y", "z"]
    assert fullform(r"x \[LeftRightArrow] y") == ["LeftRightArrow", "x", "y"]
    assert fullform("x \N{LEFT RIGHT ARROW} y \N{LEFT RIGHT ARROW} z") == \
        ["LeftRightArrow", "x", "y", "z"]

    # Therefore is right-associative: x ∴ y ∴ z groups as x ∴ (y ∴ z).
    # https://reference.wolfram.com/language/ref/Therefore.html
    assert fullform(r"x \[Therefore] y") == ["Therefore", "x", "y"]
    assert fullform("x \N{THEREFORE} y \N{THEREFORE} z") == \
        ["Therefore", "x", ["Therefore", "y", "z"]]

    # Prefix operators.
    assert fullform(r"\[Del] x") == ["Del", "x"]
    assert fullform("\N{NABLA}x") == ["Del", "x"]
    assert fullform("\N{NABLA}(x + y)") == ["Del", ["Plus", "x", "y"]]
    assert fullform(r"\[Square] x") == ["Square", "x"]
    assert fullform("\N{WHITE SQUARE}x") == ["Square", "x"]

    # Matchfix AngleBracket, with one or more comma-separated arguments.
    assert fullform(r"\[LeftAngleBracket] x \[RightAngleBracket]") == \
        ["AngleBracket", "x"]
    assert fullform(r"\[LeftAngleBracket] a, b, c \[RightAngleBracket]") == \
        ["AngleBracket", "a", "b", "c"]

    # Precedence matches Wolfram: CirclePlus binds looser than Times but
    # tighter than Plus; CircleTimes binds tighter than Times; Del (prefix)
    # binds looser than Power.
    assert fullform("a + b \N{CIRCLED PLUS} c") == ["Plus", "a", ["CirclePlus", "b", "c"]]
    assert fullform("a \N{CIRCLED PLUS} b*c") == ["CirclePlus", "a", ["Times", "b", "c"]]
    assert fullform("a*b \N{CIRCLED TIMES} c*d") == \
        ["Times", "a", ["CircleTimes", "b", "c"], "d"]
    assert fullform("\N{NABLA}x^2") == ["Del", ["Power", "x", "2"]]

    # End-to-end: with no built-in meaning each head becomes an undefined
    # Function once converted to SymPy.
    x, y = symbols("x y")
    assert parse_mathematica(r"x \[CirclePlus] y") == Function("CirclePlus")(x, y)
    assert parse_mathematica(r"\[Del] x") == Function("Del")(x)
    assert parse_mathematica(r"\[LeftAngleBracket] x, y \[RightAngleBracket]") == \
        Function("AngleBracket")(x, y)


def test_parse_mathematica_to_fullformlist():
    # The function returns the intermediate FullForm as nested Python lists,
    # with heads and atoms represented as strings.
    assert parse_mathematica_to_fullformlist("Sin[x]^2 Tan[y]") == \
        ['Times', ['Power', ['Sin', 'x'], '2'], ['Tan', 'y']]
    assert parse_mathematica_to_fullformlist("x*(a + b)") == \
        ['Times', 'x', ['Plus', 'a', 'b']]
    assert parse_mathematica_to_fullformlist("F[7, 5, 3]") == ['F', '7', '5', '3']

    # Named characters are resolved before tokenizing, unlike the raw tokenizer
    # stages, so both the \[Name] form and the literal character are accepted.
    assert parse_mathematica_to_fullformlist(r"a \[CirclePlus] b") == \
        ['CirclePlus', 'a', 'b']
    assert parse_mathematica_to_fullformlist("a \N{STAR OPERATOR} b \N{STAR OPERATOR} c") == \
        ['Star', 'a', 'b', 'c']

    # The method form on the parser class is equivalent.
    assert MathematicaParser().parse_fullformlist("a + b") == ['Plus', 'a', 'b']


def test_mathematica_derivative_precedence():
    # ``'`` is positional: it takes as operand whatever precedes it, bounded by
    # operators binding tighter than itself.  All of these were checked against
    # Wolfram Mathematica's FullForm.
    chain = parse_mathematica_to_fullformlist
    d1f = [["Derivative", "1"], "f"]
    d2f = [["Derivative", "2"], "f"]

    # tighter than ``^``, ``*``, ``+``, ``.`` and ``;;`` -- the prime takes only
    # the token on its left
    assert chain("a^b'") == ["Power", "a", [["Derivative", "1"], "b"]]
    assert chain("a*b'") == ["Times", "a", [["Derivative", "1"], "b"]]
    assert chain("a+b'") == ["Plus", "a", [["Derivative", "1"], "b"]]
    assert chain("a.b'") == ["Dot", "a", [["Derivative", "1"], "b"]]
    # ... but looser than the ``@`` family, so there it takes the whole left side
    assert chain("a@b'") == [["Derivative", "1"], ["a", "b"]]
    assert chain("a@@b'") == [["Derivative", "1"], ["Apply", "a", "b"]]
    assert chain("a/@b'") == [["Derivative", "1"], ["Map", "a", "b"]]
    assert chain("a@b@c'") == [["Derivative", "1"], ["a", ["b", "c"]]]
    assert chain("c@f'[x]") == [[["Derivative", "1"], ["c", "f"]], "x"]
    # a prime *before* the application still binds to the head alone
    assert chain("f'@x") == [d1f, "x"]
    assert chain("f'^2") == ["Power", d1f, "2"]

    # ``!`` and ``[`` bind tighter than ``'`` and are applied first
    assert chain("f'!") == ["Factorial", d1f]
    assert chain("f!'") == [["Derivative", "1"], ["Factorial", "f"]]
    assert chain("f'[x]!") == ["Factorial", [d1f, "x"]]

    # ``?`` (PatternTest) binds tighter than ``'`` and than ``[``, but its left
    # operand is still the completed expression
    assert chain("f[x]?c") == ["PatternTest", ["f", "x"], "c"]
    assert chain("f'[x]?c") == ["PatternTest", [d1f, "x"], "c"]
    assert chain("f'[x][y]?c") == ["PatternTest", [[d1f, "x"], "y"], "c"]
    assert chain("a?b[c]") == [["PatternTest", "a", "b"], "c"]

    # stacked prefix/postfix operators need more than one pass
    assert chain("- -a") == ["Times", "-1", ["Times", "-1", "a"]]
    assert chain(r"\[Del]\[Del]a") == ["Del", ["Del", "a"]]
    assert chain("a&'") == [["Derivative", "1"], ["Function", "a"]]

    assert chain("f''[x]?c") == ["PatternTest", [d2f, "x"], "c"]


def test_mathematica_unary_plus():
    # ``+x`` is ``Plus[x]``, but a unary plus heading a term of a ``+`` chain is
    # absorbed: ``a + +b`` and ``+a + b`` are both plain ``Plus[a, b]``.
    chain = parse_mathematica_to_fullformlist
    assert chain("+a") == ["Plus", "a"]
    assert chain("+f[x]") == ["Plus", ["f", "x"]]
    assert chain("a*+b") == ["Times", "a", ["Plus", "b"]]
    assert chain("a + +b") == ["Plus", "a", "b"]
    assert chain("+a + b") == ["Plus", "a", "b"]
    assert chain("+a - b") == ["Plus", "a", ["Times", "-1", "b"]]
    # after an infix ``-`` the term is negated, so the unary plus is kept
    assert chain("a-+b") == ["Plus", "a", ["Times", "-1", ["Plus", "b"]]]
    assert chain("-+a") == ["Times", "-1", ["Plus", "a"]]


def test_mathematica_blank_needs_a_symbol():
    # ``x_``/``x_h`` build a Pattern only after a symbol.  After anything else
    # the blank stands on its own and merely multiplies: Mathematica reads
    # ``f[b]_c`` as ``Times[f[b], Blank[c]]``.
    chain = parse_mathematica_to_fullformlist
    assert chain("a_b") == ["Pattern", "a", ["Blank", "b"]]
    assert chain("x_Integer") == ["Pattern", "x", ["Blank", "Integer"]]
    assert chain("a_") == ["Pattern", "a", ["Blank"]]
    assert chain("f[b]_c") == ["Times", ["f", "b"], ["Blank", "c"]]
    assert chain("f[b]_") == ["Times", ["f", "b"], ["Blank"]]
    assert chain("f[b]__") == ["Times", ["f", "b"], ["BlankSequence"]]
    assert chain("f[b]___") == ["Times", ["f", "b"], ["BlankNullSequence"]]
    assert chain("f[b]_.") == ["Times", ["f", "b"], ["Optional", ["Blank"]]]
    assert chain("a'_") == ["Times", [["Derivative", "1"], "a"], ["Blank"]]
    assert chain("a____") == \
        ["Times", ["Pattern", "a", ["BlankNullSequence"]], ["Blank"]]
    assert chain("_?NumberQ") == ["PatternTest", ["Blank"], "NumberQ"]
    assert chain("f[x_?NumberQ]") == \
        ["f", ["PatternTest", ["Pattern", "x", ["Blank"]], "NumberQ"]]


def test_parse_mathematica_derivative():
    # Mathematica curries derivative heads (``f'[x]`` is ``Derivative[1][f][x]``);
    # SymPy has no curried heads, so the idiom is translated as a whole.
    f, g = Function("f"), Function("g")
    assert parse_mathematica("f'[x]") == f(x).diff(x)
    assert parse_mathematica("f''[x]") == f(x).diff(x, 2)
    assert parse_mathematica("f'[x]*g[x]") == f(x).diff(x)*g(x)
    assert parse_mathematica("2 f''[x] + g'[y]") == 2*f(x).diff(x, 2) + g(y).diff(y)
    assert parse_mathematica("Derivative[1][f][x]") == f(x).diff(x)
    assert parse_mathematica("Derivative[1,1][f][x,y]") == f(x, y).diff(x, y)
    assert parse_mathematica("Derivative[2,1][f][x,y]") == f(x, y).diff(x, x, y)
    assert parse_mathematica("Sin'[x]") == sin(x).diff(x)

    # a derivative taken at something that is not a plain symbol needs a Subs,
    # since the differentiation is with respect to the argument slot
    e = parse_mathematica("f'[2 x]")
    assert e.func.__name__ == "Subs"
    assert e.doit() == f(2*x).diff(x)/2

    # applied to nothing, the derivative operator itself is the value
    e = parse_mathematica("f'")
    assert isinstance(e, Lambda)
    assert e(x) == f(x).diff(x)
