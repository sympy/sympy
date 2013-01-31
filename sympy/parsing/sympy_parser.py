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

def _token_splittable(token):
    """
    Predicate for whether a token name can be split into multiple tokens.

    A token is splittable if it does not contain an underscore charater and
    it is not the name of a Greek letter. This is used to implicitly convert
    expressions like 'xyz' into 'x*y*z'.
    """
    if '_' in token:
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
    ... transformations=(standard_transformations +
    ... (implicit_multiplication_application,)))
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
    prevTok = (None, None)

    tokens.append((None, None))  # so zip traverses all tokens
    for tok, nextTok in zip(tokens, tokens[1:]):
        tokNum, tokVal = tok
        nextTokNum, nextTokVal = nextTok
        if tokNum == NAME:
            name = tokVal

            if (name in ['True', 'False', 'None']
                or iskeyword(name)
                or name in local_dict
                # Don't convert attribute access
                or (prevTok[0] == OP and prevTok[1] == '.')
                # Don't convert keyword arguments
                or (prevTok[0] == OP and prevTok[1] in ('(', ',')
                    and nextTokNum == OP and nextTokVal == '=')):
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
            result.append((tokNum, tokVal))

        prevTok = (tokNum, tokVal)

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
    """Treats XOR, ``^``, as exponentiation, ``**``."""
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
    """Converts floats into ``Rational``. Run AFTER ``auto_number``."""
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

#: Standard transformations for :func:`parse_expr`.
#: Inserts calls to :class:`Symbol`, :class:`Integer`, and other SymPy
#: datatypes and allows the use of standard factorial notation (e.g. ``x!``).
standard_transformations = (auto_symbol, auto_number, factorial_notation)


def stringify_expr(s, local_dict, global_dict, transformations):
    """
    Converts the string ``s`` to Python code, in ``local_dict``

    Generally, ``parse_expr`` should be used.
    """

    tokens = []
    input_code = StringIO(s.strip())
    for toknum, tokval, _, _, _ in generate_tokens(input_code.readline):
        tokens.append((toknum, tokval))

    for transform in transformations:
        tokens = transform(tokens, local_dict, global_dict)

    return untokenize(tokens)


def eval_expr(code, local_dict, global_dict):
    """
    Evaluate Python code generated by ``stringify_expr``.

    Generally, ``parse_expr`` should be used.
    """
    expr = eval(
        code, global_dict, local_dict)  # take local objects in preference

    return expr


def parse_expr(s, local_dict=None, transformations=standard_transformations,
               global_dict=None):
    """Converts the string ``s`` to a SymPy expression, in ``local_dict``

    Parameters
    ==========

    s : str
        The string to parse.

    local_dict : dict, optional
        A dictionary of local variables to use when parsing.

    global_dict : dict, optional
        A dictionary of global variables. By default, this is initialized
        with ``from sympy import *``; provide this parameter to override
        this behavior (for instance, to parse ``"Q & S"``).

    transformations : tuple, optional
        A tuple of transformation functions used to modify the tokens of the
        parsed expression before evaluation. The default transformations
        convert numeric literals into their SymPy equivalents, convert
        undefined variables into SymPy symbols, and allow the use of standard
        mathematical factorial notation (e.g. ``x!``).


    Examples
    ========

    >>> from sympy.parsing.sympy_parser import parse_expr
    >>> parse_expr("1/2")
    1/2
    >>> type(_)
    <class 'sympy.core.numbers.Half'>
    >>> from sympy.parsing.sympy_parser import standard_transformations,\\
    ... implicit_multiplication_application
    >>> transformations = (standard_transformations +
    ...     (implicit_multiplication_application,))
    >>> parse_expr("2x", transformations=transformations)
    2*x

    See Also
    ========

    stringify_expr, eval_expr, standard_transformations,
    implicit_multiplication_application

    """

    if local_dict is None:
        local_dict = {}

    if global_dict is None:
        global_dict = {}
        exec 'from sympy import *' in global_dict

    code = stringify_expr(s, local_dict, global_dict, transformations)
    return eval_expr(code, local_dict, global_dict)
