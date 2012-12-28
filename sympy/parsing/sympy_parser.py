"""Transform a string with Python-like source code into SymPy expression. """

from sympy_tokenize import \
    generate_tokens, untokenize, TokenError, \
    NUMBER, STRING, NAME, OP, ENDMARKER

from keyword import iskeyword
from StringIO import StringIO
import re
import unicodedata

from sympy.core.basic import Basic, C

_re_repeated = re.compile(r"^(\d*)\.(\d*)\[(\d+)\]$")
UNSPLITTABLE_TOKEN_NAMES = ['_kern']


def _token_splittable(token):
    """
    Predicate for whether a token name can be split into multiple tokens.

    A token is splittable if it does not contain an underscore charater and
    it is not the name of a Greek letter. This is used to implicitly convert
    expressions like 'xyz' into 'x*y*z'.
    """
    if '_' in token:
        return False
    elif token in UNSPLITTABLE_TOKEN_NAMES:
        return False
    else:
        try:
            return not unicodedata.lookup('GREEK SMALL LETTER ' + token)
        except KeyError:
            pass
    return True


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


class AppliedFunction(object):
    """
    A group of tokens representing a function and its arguments.

    `exponent` is for handling the shorthand sin^2, ln^2, etc.
    """
    def __init__(self, function, args, exponent=None):
        if exponent is None:
            exponent = []
        self.function = function
        self.args = args
        self.exponent = exponent
        self.items = ['function', 'args', 'exponent']

    def expand(self):
        """Return a list of tokens representing the function"""
        result = []
        result.append(self.function)
        result.extend(self.args)
        return result

    def __getitem__(self, index):
        return getattr(self, self.items[index])

    def __repr__(self):
        return "AppliedFunction(%s, %s, %s)" % (self.function, self.args,
                                                self.exponent)


class ParenthesisGroup(list):
    """List of tokens representing an expression in parentheses."""
    pass


def _flatten(result):
    result2 = []
    for tok in result:
        if isinstance(tok, AppliedFunction):
            result2.extend(tok.expand())
        else:
            result2.append(tok)
    return result2


def _group_parentheses(tokens, local_dict, global_dict):
    """Group tokens between parentheses with ParenthesisGroup.

    Also processes those tokens recursively.

    """
    result = []
    stacks = []
    stacklevel = 0
    for token in tokens:
        if token[0] == OP:
            if token[1] == '(':
                stacks.append(ParenthesisGroup([]))
                stacklevel += 1
            elif token[1] == ')':
                stacks[-1].append(token)
                stack = stacks.pop()

                if len(stacks) > 0:
                    # We don't recurse here since the upper-level stack
                    # would reprocess these tokens
                    stacks[-1].extend(stack)
                else:
                    # Recurse here to handle nested parentheses
                    # Strip off the outer parentheses to avoid an infinite loop
                    inner = stack[1:-1]
                    inner = implicit_multiplication_application(inner,
                                                                local_dict,
                                                                global_dict)
                    parenGroup = [stack[0]] + inner + [stack[-1]]
                    result.append(ParenthesisGroup(parenGroup))
                stacklevel -= 1
                continue
        if stacklevel:
            stacks[-1].append(token)
        else:
            result.append(token)
    return result


def _apply_functions(tokens, local_dict, global_dict):
    """Convert a NAME token + ParenthesisGroup into an AppliedFunction.

    Note that ParenthesisGroups, if not applied to any function, are
    converted back into lists of tokens.

    """
    result = []
    symbol = None
    for tok in tokens:
        if tok[0] == NAME:
            symbol = tok
            result.append(tok)
        elif isinstance(tok, ParenthesisGroup):
            if symbol:
                result[-1] = AppliedFunction(symbol, tok)
                symbol = None
            else:
                result.extend(tok)
        else:
            symbol = None
            result.append(tok)
    return result


def _split_symbols(tokens, local_dict, global_dict):
    result = []
    for tok in tokens:
        if isinstance(tok, AppliedFunction):
            if tok.function[1] == 'Symbol' and len(tok.args) == 3:
                # Middle token, get value, strip quotes
                symbol = tok.args[1][1][1:-1]
                if _token_splittable(symbol):
                    for char in symbol:
                        result.append(AppliedFunction(
                            tok.function,
                            [(OP, '('), (NAME, str(repr(char))), (OP, ')')]
                        ))
                    continue
        result.append(tok)
    return result


def _implicit_multiplication(tokens, local_dict, global_dict):
    """Implicitly adds '*' tokens.

    Cases:

    - Two AppliedFunctions next to each other ("sin(x)cos(x)")

    - AppliedFunction next to an open parenthesis ("sin x (cos x + 1)")

    - A close parenthesis next to an AppliedFunction ("(x+2)sin x")\

    - A closeparenthesis next to an open parenthesis ("(x+2)(x+3)")

    - An AppliedFunction next to an implicitly applied function ("sin(x)cos
      x")

    """
    result = []
    for tok, nextTok in zip(tokens, tokens[1:]):
        result.append(tok)
        if (isinstance(tok, AppliedFunction) and
                isinstance(nextTok, AppliedFunction)):
            result.append((OP, '*'))
        elif (isinstance(tok, AppliedFunction) and
              nextTok[0] == OP and nextTok[1] == '('):
            # Applied function followed by an open parenthesis
            result.append((OP, '*'))
        elif (tok[0] == OP and tok[1] == ')' and
              isinstance(nextTok, AppliedFunction)):
            # Close parenthesis followed by an applied function
            result.append((OP, '*'))
        elif (tok[0] == OP and tok[1] == ')' and
              nextTok[0] == NAME):
            # Close parenthesis followed by an implicitly applied function
            result.append((OP, '*'))
        elif (tok[0] == nextTok[0] == OP
              and tok[1] == ')' and nextTok[1] == '('):
            # Close parenthesis followed by an open parenthesis
            result.append((OP, '*'))
        elif (isinstance(tok, AppliedFunction) and nextTok[0] == NAME):
            # Applied function followed by implicitly applied function
            result.append((OP, '*'))
    result.append(tokens[-1])
    return result


def _implicit_application(tokens, local_dict, global_dict):
    """Adds parentheses as needed after functions.

    Also processes functions raised to powers."""
    result = []
    appendParen = 0  # number of closing parentheses to add
    skip = False  # number of tokens to delay before adding a ')' (to
                  # capture **, ^, etc.)
    exponent = None
    for tok, nextTok in zip(tokens, tokens[1:]):
        result.append(tok)
        if (tok[0] == NAME and
            nextTok[0] != OP and
                nextTok[0] != ENDMARKER):
            func = global_dict.get(tok[1])
            is_Function = getattr(func, 'is_Function', False)
            if (is_Function or
                (callable(func) and not hasattr(func, 'is_Function')) or
                    isinstance(nextTok, AppliedFunction)):
                result.append((OP, '('))
                appendParen += 1
        elif isinstance(tok, AppliedFunction) and not tok.args:
            # This occurs when we had a function raised to a power - it has
            # no arguments
            result[-1] = tok.function
            exponent = tok.exponent
            if nextTok[0] != OP and nextTok[1] != '(':
                result.append((OP, '('))
                appendParen += 1
        elif appendParen:
            if nextTok[0] == OP and nextTok[1] in ('^', '**'):
                skip = 1
                continue
            if skip:
                skip -= 1
                continue
            result.append((OP, ')'))
            appendParen -= 1
            if exponent and not appendParen:
                result.extend(exponent)
                exponent = None

    result.append(tokens[-1])

    if appendParen:
        result.extend([(OP, ')')] * appendParen)
    if exponent:
        result.extend(exponent)

    return result


def _function_exponents(tokens, local_dict, global_dict):
    """Preprocess functions raised to powers."""
    result = []
    need_exponent = False
    for tok, nextTok in zip(tokens, tokens[1:]):
        result.append(tok)
        if (tok[0] == NAME and nextTok[0] == OP
                and nextTok[1] in ('**', '^')):
            if getattr(global_dict.get(tok[1]), 'is_Function', False):
                result[-1] = AppliedFunction(tok, [])
                need_exponent = True
        elif need_exponent:
            del result[-1]
            result[-1].exponent.append(tok)
            if isinstance(tok, AppliedFunction):
                need_exponent = False
    result.append(tokens[-1])
    return result


def implicit_multiplication_application(result, local_dict, global_dict):
    """Allows a slightly relaxed syntax.

    - Parentheses for single-argument method calls are optional.

    - Multiplication is implicit.

    - Symbol names can be split (i.e. spaces are not needed between
      symbols).

    - Functions can be exponentiated.

    Example:

    >>> from sympy.parsing.sympy_parser import (parse_expr,
    ... standard_transformations, implicit_multiplication_application)
    >>> parse_expr("10sin**2 x**2 + 3xyz + tan theta",
    ... transformations=standard_transformations +
    ... (implicit_multiplication_application,))
    3*x*y*z + 10*sin(x**2)**2 + tan(theta)

    """
    # These are interdependent steps, so we don't expose them separately
    for step in (_group_parentheses,
                 _apply_functions,
                 _split_symbols,
                 _function_exponents,
                 _implicit_application,
                 _implicit_multiplication):
        result = step(result, local_dict, global_dict)

    result = _flatten(result)
    return result


def auto_symbol(tokens, local_dict, global_dict):
    """Inserts calls to ``Symbol`` for undefined variables."""
    result = []
    for toknum, tokval in tokens:
        if toknum == NAME:
            name = tokval

            if (name in ['True', 'False', 'None']
                or iskeyword(name)
                    or name in local_dict):
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
        else:
            result.append((toknum, tokval))

    return result


def factorial_notation(tokens, local_dict, global_dict):
    """Allows standard notation for factorial."""
    result = []
    prevtoken = ''
    for toknum, tokval in tokens:
        if toknum == OP:
            op = tokval

            if op == '!!':
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

    return result


def convert_xor(tokens, local_dict, global_dict):
    """Treats XOR, '^', as exponentiation, '**'."""
    result = []
    for toknum, tokval in tokens:
        if toknum == OP:
            if tokval == '^':
                result.append((OP, '**'))
            else:
                result.append((toknum, tokval))
        else:
            result.append((toknum, tokval))

    return result


def auto_number(tokens, local_dict, global_dict):
    """Converts numeric literals to use SymPy equivalents.

    Complex numbers use ``I``; integer literals use ``Integer``, float
    literals use ``Float``, and repeating decimals use ``Rational``.

    """
    result = []
    prevtoken = ''

    for toknum, tokval in tokens:
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
                else:
                    seq = [(NAME, 'Float'), (OP, '('),
                           (NUMBER, repr(str(number))), (OP, ')')]
            else:
                seq = [(NAME, 'Integer'), (OP, '('), (
                    NUMBER, number), (OP, ')')]

            result.extend(seq + postfix)
        else:
            result.append((toknum, tokval))

    return result


def rationalize(tokens, local_dict, global_dict):
    """Converts floats into ``Rational``s. Run AFTER ``auto_number``."""
    result = []
    passed_float = False
    for toknum, tokval in tokens:
        if toknum == NAME:
            if tokval == 'Float':
                passed_float = True
                tokval = 'Rational'
            result.append((toknum, tokval))
        elif passed_float == True and toknum == NUMBER:
            passed_float = False
            result.append((STRING, tokval))
        else:
            result.append((toknum, tokval))

    return result

standard_transformations = (auto_symbol, auto_number, factorial_notation)


def parse_expr(s, local_dict=None, transformations=standard_transformations):
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

    tokens = []
    input_code = StringIO(s.strip())
    for toknum, tokval, _, _, _ in generate_tokens(input_code.readline):
        tokens.append((toknum, tokval))

    for transform in transformations:
        tokens = transform(tokens, local_dict, global_dict)

    code = untokenize(tokens)
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
