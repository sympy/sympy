"""Transform a string with Python-like source code into SymPy expression. """

from sympy_tokenize import \
    generate_tokens, untokenize, TokenError, NUMBER, STRING, NAME, OP

from keyword import iskeyword
from StringIO import StringIO
import re

from sympy.core.basic import Basic, C

_re_repeated = re.compile(r"^(\d*)\.(\d*)\[(\d+)\]$")

def _add_factorial_tokens(name, result):
    if result == [] or result[-1][1] == '(':
        raise TokenError()

    beginning = [(NAME, name), (OP, '(')]
    end = [(OP, ')')]

    diff = 0
    length = len(result)

    for index, token in enumerate(result[::-1]):
        toknum, tokval = token
        i = length-index-1

        if tokval == ')':
            diff += 1
        elif tokval == '(':
            diff -= 1

        if diff == 0:
            if i-1 >= 0 and result[i-1][0] == NAME:
                return result[:i-1] + beginning + result[i-1:] + end
            else:
                return result[:i] + beginning + result[i:] + end

    return result

def _transform(s, local_dict, global_dict, rationalize, convert_xor):
    g = generate_tokens(StringIO(s).readline)

    result = []
    prevtoken = ''

    for toknum, tokval, _, _, _ in g:
        if toknum == NUMBER:
            number = tokval
            postfix = []

            if number.endswith('j') or number.endswith('J'):
                number = number[:-1]
                postfix = [(OP, '*'), (NAME, 'I')]

            if '.' in number or (('e' in number or 'E' in number) and
                    not (number.startswith('0x') or number.startswith('0X'))):
                match = _re_repeated.match(number)

                if match is not None:
                    # Clear repeating decimals, e.g. 3.4[31] -> (3 + 4/10 + 31/990)
                    pre, post, repetend = match.groups()

                    zeros = '0'*len(post)
                    post, repetends = [w.lstrip('0') for w in [post, repetend]] # or else interpreted as octal

                    a = pre or '0'
                    b, c = post or '0', '1' + zeros
                    d, e = repetends, ('9'*len(repetend)) + zeros

                    seq = [
                        (OP, '('),
                            (NAME, 'Integer'), (OP, '('), (NUMBER, a), (OP, ')'),
                        (OP, '+'),
                            (NAME, 'Rational'), (OP, '('), (NUMBER, b), (OP, ','), (NUMBER, c), (OP, ')'),
                        (OP, '+'),
                            (NAME, 'Rational'), (OP, '('), (NUMBER, d), (OP, ','), (NUMBER, e), (OP, ')'),
                        (OP, ')'),
                    ]
                elif rationalize:
                    seq = [(NAME, 'Rational'), (OP, '('), (STRING, repr(str(number))), (OP, ')')]
                else:
                    seq = [(NAME, 'Float'), (OP, '('), (NUMBER, repr(str(number))), (OP, ')')]
            else:
                seq = [(NAME, 'Integer'), (OP, '('), (NUMBER, number), (OP, ')')]

            result.extend(seq + postfix)
        elif toknum == NAME:
            name = tokval

            if name in ['True', 'False', 'None'] or iskeyword(name) or name in local_dict:
                result.append((NAME, name))
                continue
            elif name in global_dict:
                obj = global_dict[name]

                if isinstance(obj, (Basic, type)) or callable(obj):
                    result.append((NAME, name))
                    continue

            result.extend([
                (NAME, 'Symbol'),
                (OP, '('),
                (NAME, repr(str(name))),
                (OP, ')'),
            ])
        elif toknum == OP:
            op = tokval

            if op == '^' and convert_xor:
                result.append((OP, '**'))
            elif op == '!!':
                if prevtoken == '!' or prevtoken == '!!':
                    raise TokenError
                result = _add_factorial_tokens('factorial2', result)
            elif op == '!':
                if prevtoken == '!' or prevtoken == '!!':
                    raise TokenError
                result = _add_factorial_tokens('factorial', result)
            else:
                result.append((OP, op))
        else:
            result.append((toknum, tokval))

        prevtoken = tokval

    return untokenize(result)

def parse_expr(s, local_dict=None, rationalize=False, convert_xor=False):
    """
    Converts the string ``s`` to a SymPy expression, in ``local_dict``

    Examples
    ========

    >>> from sympy.parsing.sympy_parser import parse_expr

    >>> parse_expr("1/2")
    1/2
    >>> type(_)
    <class 'sympy.core.numbers.Half'>

    """
    if local_dict is None:
        local_dict = {}

    global_dict = {}
    exec 'from sympy import *' in global_dict

    # keep autosimplification from joining Integer or
    # minus sign into a Mul; this modification doesn't
    # prevent the 2-arg Mul from becoming and Add.
    hit = False
    if '(' in s:
        kern = '_kern'
        while kern in s:
            kern += "_"
        s = re.sub(r'(\d *\*|-) *\(', r'\1%s*(' % kern, s)
        hit = kern in s

    code = _transform(s.strip(), local_dict, global_dict, rationalize, convert_xor)
    expr = eval(code, global_dict, local_dict) # take local objects in preference

    if not hit:
        return expr
    try:
        return expr.xreplace({C.Symbol(kern): 1})
    except (TypeError, AttributeError):
        return expr
