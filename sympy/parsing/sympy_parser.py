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
        i = length - index - 1

        if tokval == ')':
            diff += 1
        elif tokval == '(':
            diff -= 1

        if diff == 0:
            if i - 1 >= 0 and result[i - 1][0] == NAME:
                return result[:i - 1] + beginning + result[i - 1:] + end
            else:
                return result[:i] + beginning + result[i:] + end

    return result

def _flatten(result):
    result2 = []
    for tok in result:
        if type(tok[0]) == tuple:
            result2.append(tok[0])
            result2.extend(tok[1])
        elif type(tok) == list:
            result2.extend(tok)
        else:
            result2.append(tok)
    return result2

def _implicit_multiplication_application(result):
    # Step 1: parenthesis grouping
    # Group sequences of '()' operators
    result2 = []
    stacks = []
    stacklevel = 0
    for token in result:
        if token[0] == OP:
            if token[1] == '(':
                stacks.append([])
                stacklevel += 1
            elif token[1] == ')':
                stacks[-1].append(token)
                stack = stacks.pop()
                if len(stacks) > 0:
                    stacks[-1].extend(stack)
                else:
                    temp = stack[1:-1]
                    temp = _implicit_multiplication_application(temp)
                    stack = [stack[0]] + temp + [stack[-1]]
                    result2.append(stack)
                stacklevel -= 1
                continue
        if stacklevel:
            stacks[-1].append(token)
        else:
            result2.append(token)
    # print 'STEP1', result2

    # Step 2: symbol/function application
    # Group NAME followed by a list
    result3 = []
    symbol = None
    for tok in result2:
        if type(tok) == tuple and tok[0] == NAME:
            symbol = tok
            result3.append(tok)
        elif type(tok) == list:
            if symbol:
                result3[-1] = (symbol, tok)
                symbol = None
            else:
                # Need some sort of deep-flatten
                # If we extend, the next step can't handle implicit paren
                # group multiplication
                # If we append, untokenize doesn't work
                result3.extend(tok)
        else:
            symbol = None
            result3.append(tok)
    # print 'STEP2', result3

    # Step 3: implicit multiplication
    # Look for two NAMEs in a row or NAME then paren group or paren group
    # then NAME
    result4 = []
    for tok, nextTok in zip(result3, result3[1:]):
        result4.append(tok)
        if type(tok[0]) == type(nextTok[0]) == tuple:
            if tok[0][0] == nextTok[0][0] == NAME:
                result4.append((OP, '*'))
        elif (type(tok[0]) == tuple and
              nextTok[0] == OP and nextTok[1] == '('):
            # Applied function followed by an open parenthesis
            result4.append((OP, '*'))
        elif (tok[0] == OP and tok[1] == ')' and
              type(nextTok[0]) == tuple):
            # Close parenthesis followed by an applied function
            result4.append((OP, '*'))
    result4.append(result3[-1])
    # print 'STEP3', result4

    # Step 4: implicit application
    # Look for a NAME that is not followed by an OP '(' and apply it to the
    # following token
    result5 = []
    appendParen = False
    for tok, nextTok in zip(result4, result4[1:]):
        result5.append(tok)
        if tok[0] == NAME and nextTok[0] != OP and nextTok[1] != '(':
            if tok[1] not in ('None', 'True', 'False'):
                # TODO: better way to handle this
                result5.append((OP, '('))
                appendParen = True
        elif appendParen:
            result5.append((OP, ')'))
            appendParen = False
    result5.append(result4[-1])
    # print 'STEP4', result5

    # Step 5: Flatten
    result6 = _flatten(result5)
    # print 'STEP5', result6

    return result6

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
                    post, repetends = [w.lstrip('0') for w in [post, repetend]]
                                                # or else interpreted as octal

                    a = pre or '0'
                    b, c = post or '0', '1' + zeros
                    d, e = repetends, ('9'*len(repetend)) + zeros

                    seq = [
                        (OP, '('),
                        (NAME,
                         'Integer'), (OP, '('), (NUMBER, a), (OP, ')'),
                        (OP, '+'),
                        (NAME, 'Rational'), (OP, '('), (
                            NUMBER, b), (OP, ','), (NUMBER, c), (OP, ')'),
                        (OP, '+'),
                        (NAME, 'Rational'), (OP, '('), (
                            NUMBER, d), (OP, ','), (NUMBER, e), (OP, ')'),
                        (OP, ')'),
                    ]
                elif rationalize:
                    seq = [(NAME, 'Rational'), (OP,
                            '('), (STRING, repr(str(number))), (OP, ')')]
                else:
                    seq = [(NAME, 'Float'), (OP, '('),
                           (NUMBER, repr(str(number))), (OP, ')')]
            else:
                seq = [(NAME, 'Integer'), (OP, '('), (
                    NUMBER, number), (OP, ')')]

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

            if name != '_kern':
                for var in name:
                    result.extend([
                        (NAME, 'Symbol'),
                        (OP, '('),
                        (NAME, repr(str(var))),
                        (OP, ')'),
                    ])
            else:
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

    result = _implicit_multiplication_application(result)

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
    # prevent the 2-arg Mul from becoming an Add, however.
    hit = False
    if '(' in s:
        kern = '_kern'
        while kern in s:
            kern += "_"
        s = re.sub(r'(\d *\*|-) *\(', r'\1%s*(' % kern, s)
        hit = kern in s

    code = _transform(
        s.strip(), local_dict, global_dict, rationalize, convert_xor)
    expr = eval(
        code, global_dict, local_dict)  # take local objects in preference

    if not hit:
        return expr
    rep = {C.Symbol(kern): 1}
    def _clear(expr):
        if hasattr(expr, 'xreplace'):
            return expr.xreplace(rep)
        elif isinstance(expr, (list, tuple, set)):
            return type(expr)([_clear(e) for e in expr])
        if hasattr(expr, 'subs'):
            return expr.subs(rep)
        return expr
    return _clear(expr)
