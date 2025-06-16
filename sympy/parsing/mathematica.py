from __future__ import annotations
import sys
import re
import typing
from itertools import product
from typing import Any, Callable

import sympy
from sympy import Mul, Add, Pow, Rational, log, exp, sqrt, cos, sin, tan, asin, acos, acot, asec, acsc, sinh, cosh, tanh, asinh, \
    acosh, atanh, acoth, asech, acsch, expand, im, flatten, polylog, cancel, expand_trig, sign, simplify, \
    UnevaluatedExpr, S, atan, atan2, Mod, Max, Min, rf, Ei, Si, Ci, airyai, airyaiprime, airybi, primepi, prime, \
    isprime, cot, sec, csc, csch, sech, coth, Function, E, I, pi, Tuple, GreaterThan, StrictGreaterThan, StrictLessThan, \
    LessThan, Equality, Or, And, Lambda, Integer, Dummy, symbols
from sympy.core.sympify import sympify, _sympify
from sympy.functions.special.bessel import airybiprime
from sympy.functions.special.error_functions import li
from sympy.utilities.exceptions import sympy_deprecation_warning


def mathematica(s, additional_translations=None):
    sympy_deprecation_warning(
        """The ``mathematica`` function for the Mathematica parser is now
deprecated. Use ``parse_mathematica`` instead.
The parameter ``additional_translation`` can be replaced by SymPy's
.replace( ) or .subs( ) methods on the output expression instead.""",
        deprecated_since_version="1.11",
        active_deprecations_target="mathematica-parser-new",
    )
    parser = MathematicaParser(additional_translations)
    return sympify(parser._parse_old(s))


def parse_mathematica(s):
    """
    Translate a string containing a Wolfram Mathematica expression to a SymPy
    expression.

    If the translator is unable to find a suitable SymPy expression, the
    ``FullForm`` of the Mathematica expression will be output, using SymPy
    ``Function`` objects as nodes of the syntax tree.

    Examples
    ========

    >>> from sympy.parsing.mathematica import parse_mathematica
    >>> parse_mathematica("Sin[x]^2 Tan[y]")
    sin(x)**2*tan(y)
    >>> e = parse_mathematica("F[7,5,3]")
    >>> e
    F(7, 5, 3)
    >>> from sympy import Function, Max, Min
    >>> e.replace(Function("F"), lambda *x: Max(*x)*Min(*x))
    21

    Both standard input form and Mathematica full form are supported:

    >>> parse_mathematica("x*(a + b)")
    x*(a + b)
    >>> parse_mathematica("Times[x, Plus[a, b]]")
    x*(a + b)

    To get a matrix from Wolfram's code:

    >>> m = parse_mathematica("{{a, b}, {c, d}}")
    >>> m
    ((a, b), (c, d))
    >>> from sympy import Matrix
    >>> Matrix(m)
    Matrix([
    [a, b],
    [c, d]])

    If the translation into equivalent SymPy expressions fails, an SymPy
    expression equivalent to Wolfram Mathematica's "FullForm" will be created:

    >>> parse_mathematica("x_.")
    Optional(Pattern(x, Blank()))
    >>> parse_mathematica("Plus @@ {x, y, z}")
    Apply(Plus, (x, y, z))
    >>> parse_mathematica("f[x_, 3] := x^3 /; x > 0")
    SetDelayed(f(Pattern(x, Blank()), 3), Condition(x**3, x > 0))
    """
    parser = MathematicaParser()
    return parser.parse(s)


def _parse_Function(*args):
    if len(args) == 1:
        arg = args[0]
        Slot = Function("Slot")
        slots = arg.atoms(Slot)
        numbers = [a.args[0] for a in slots]
        number_of_arguments = max(numbers)
        if isinstance(number_of_arguments, Integer):
            variables = symbols(f"dummy0:{number_of_arguments}", cls=Dummy)
            return Lambda(variables, arg.xreplace({Slot(i+1): v for i, v in enumerate(variables)}))
        return Lambda((), arg)
    elif len(args) == 2:
        variables = args[0]
        body = args[1]
        return Lambda(variables, body)
    else:
        raise SyntaxError("Function node expects 1 or 2 arguments")


def _deco(cls):
    cls._initialize_class()
    return cls


@_deco
class MathematicaParser:
    """
    An instance of this class converts a string of a Wolfram Mathematica
    expression to a SymPy expression.

    The main parser acts internally in three stages:

    1. tokenizer: tokenizes the Mathematica expression and adds the missing *
        operators. Handled by ``_from_mathematica_to_tokens(...)``
    2. full form list: sort the list of strings output by the tokenizer into a
        syntax tree of nested lists and strings, equivalent to Mathematica's
        ``FullForm`` expression output. This is handled by the function
        ``_from_tokens_to_fullformlist(...)``.
    3. SymPy expression: the syntax tree expressed as full form list is visited
        and the nodes with equivalent classes in SymPy are replaced. Unknown
        syntax tree nodes are cast to SymPy ``Function`` objects. This is
        handled by ``_from_fullformlist_to_sympy(...)``.

    """

    # left: Mathematica, right: SymPy
    CORRESPONDENCES = {
        'Sqrt[x]': 'sqrt(x)',
        'Rational[x,y]': 'Rational(x,y)',
        'Exp[x]': 'exp(x)',
        'Log[x]': 'log(x)',
        'Log[x,y]': 'log(y,x)',
        'Log2[x]': 'log(x,2)',
        'Log10[x]': 'log(x,10)',
        'Mod[x,y]': 'Mod(x,y)',
        'Max[*x]': 'Max(*x)',
        'Min[*x]': 'Min(*x)',
        'Pochhammer[x,y]':'rf(x,y)',
        'ArcTan[x,y]':'atan2(y,x)',
        'ExpIntegralEi[x]': 'Ei(x)',
        'SinIntegral[x]': 'Si(x)',
        'CosIntegral[x]': 'Ci(x)',
        'AiryAi[x]': 'airyai(x)',
        'AiryAiPrime[x]': 'airyaiprime(x)',
        'AiryBi[x]' :'airybi(x)',
        'AiryBiPrime[x]' :'airybiprime(x)',
        'LogIntegral[x]':' li(x)',
        'PrimePi[x]': 'primepi(x)',
        'Prime[x]': 'prime(x)',
        'PrimeQ[x]': 'isprime(x)'
    }

    # trigonometric, e.t.c.
    for arc, tri, h in product(('', 'Arc'), (
            'Sin', 'Cos', 'Tan', 'Cot', 'Sec', 'Csc'), ('', 'h')):
        fm = arc + tri + h + '[x]'
        if arc:  # arc func
            fs = 'a' + tri.lower() + h + '(x)'
        else:    # non-arc func
            fs = tri.lower() + h + '(x)'
        CORRESPONDENCES.update({fm: fs})

    REPLACEMENTS = {
        ' ': '',
        '^': '**',
        '{': '[',
        '}': ']',
    }

    RULES = {
        # a single whitespace to '*'
        'whitespace': (
            re.compile(r'''
                (?:(?<=[a-zA-Z\d])|(?<=\d\.))     # a letter or a number
                \s+                               # any number of whitespaces
                (?:(?=[a-zA-Z\d])|(?=\.\d))       # a letter or a number
                ''', re.VERBOSE),
            '*'),

        # add omitted '*' character
        'add*_1': (
            re.compile(r'''
                (?:(?<=[])\d])|(?<=\d\.))       # ], ) or a number
                                                # ''
                (?=[(a-zA-Z])                   # ( or a single letter
                ''', re.VERBOSE),
            '*'),

        # add omitted '*' character (variable letter preceding)
        'add*_2': (
            re.compile(r'''
                (?<=[a-zA-Z])       # a letter
                \(                  # ( as a character
                (?=.)               # any characters
                ''', re.VERBOSE),
            '*('),

        # convert 'Pi' to 'pi'
        'Pi': (
            re.compile(r'''
                (?:
                \A|(?<=[^a-zA-Z])
                )
                Pi                  # 'Pi' is 3.14159... in Mathematica
                (?=[^a-zA-Z])
                ''', re.VERBOSE),
            'pi'),
    }

    # Mathematica function name pattern
    FM_PATTERN = re.compile(r'''
                (?:
                \A|(?<=[^a-zA-Z])   # at the top or a non-letter
                )
                [A-Z][a-zA-Z\d]*    # Function
                (?=\[)              # [ as a character
                ''', re.VERBOSE)

    # list or matrix pattern (for future usage)
    ARG_MTRX_PATTERN = re.compile(r'''
                \{.*\}
                ''', re.VERBOSE)

    # regex string for function argument pattern
    ARGS_PATTERN_TEMPLATE = r'''
                (?:
                \A|(?<=[^a-zA-Z])
                )
                {arguments}         # model argument like x, y,...
                (?=[^a-zA-Z])
                '''

    # will contain transformed CORRESPONDENCES dictionary
    TRANSLATIONS: dict[tuple[str, int], dict[str, Any]] = {}

    # cache for a raw users' translation dictionary
    cache_original: dict[tuple[str, int], dict[str, Any]] = {}

    # cache for a compiled users' translation dictionary
    cache_compiled: dict[tuple[str, int], dict[str, Any]] = {}

    @classmethod
    def _initialize_class(cls):
        # get a transformed CORRESPONDENCES dictionary
        d = cls._compile_dictionary(cls.CORRESPONDENCES)
        cls.TRANSLATIONS.update(d)

    def __init__(self, additional_translations=None):
        self.translations = {}

        # update with TRANSLATIONS (class constant)
        self.translations.update(self.TRANSLATIONS)

        if additional_translations is None:
            additional_translations = {}

        # check the latest added translations
        if self.__class__.cache_original != additional_translations:
            if not isinstance(additional_translations, dict):
                raise ValueError('The argument must be dict type')

            # get a transformed additional_translations dictionary
            d = self._compile_dictionary(additional_translations)

            # update cache
            self.__class__.cache_original = additional_translations
            self.__class__.cache_compiled = d

        # merge user's own translations
        self.translations.update(self.__class__.cache_compiled)

    @classmethod
    def _compile_dictionary(cls, dic):
        # for return
        d = {}

        for fm, fs in dic.items():
            # check function form
            cls._check_input(fm)
            cls._check_input(fs)

            # uncover '*' hiding behind a whitespace
            fm = cls._apply_rules(fm, 'whitespace')
            fs = cls._apply_rules(fs, 'whitespace')

            # remove whitespace(s)
            fm = cls._replace(fm, ' ')
            fs = cls._replace(fs, ' ')

            # search Mathematica function name
            m = cls.FM_PATTERN.search(fm)

            # if no-hit
            if m is None:
                err = "'{f}' function form is invalid.".format(f=fm)
                raise ValueError(err)

            # get Mathematica function name like 'Log'
            fm_name = m.group()

            # get arguments of Mathematica function
            args, end = cls._get_args(m)

            # function side check. (e.g.) '2*Func[x]' is invalid.
            if m.start() != 0 or end != len(fm):
                err = "'{f}' function form is invalid.".format(f=fm)
                raise ValueError(err)

            # check the last argument's 1st character
            if args[-1][0] == '*':
                key_arg = '*'
            else:
                key_arg = len(args)

            key = (fm_name, key_arg)

            # convert '*x' to '\\*x' for regex
            re_args = [x if x[0] != '*' else '\\' + x for x in args]

            # for regex. Example: (?:(x|y|z))
            xyz = '(?:(' + '|'.join(re_args) + '))'

            # string for regex compile
            patStr = cls.ARGS_PATTERN_TEMPLATE.format(arguments=xyz)

            pat = re.compile(patStr, re.VERBOSE)

            # update dictionary
            d[key] = {}
            d[key]['fs'] = fs  # SymPy function template
            d[key]['args'] = args  # args are ['x', 'y'] for example
            d[key]['pat'] = pat

        return d

    def _convert_function(self, s):
        '''Parse Mathematica function to SymPy one'''

        # compiled regex object
        pat = self.FM_PATTERN

        scanned = ''                # converted string
        cur = 0                     # position cursor
        while True:
            m = pat.search(s)

            if m is None:
                # append the rest of string
                scanned += s
                break

            # get Mathematica function name
            fm = m.group()

            # get arguments, and the end position of fm function
            args, end = self._get_args(m)

            # the start position of fm function
            bgn = m.start()

            # convert Mathematica function to SymPy one
            s = self._convert_one_function(s, fm, args, bgn, end)

            # update cursor
            cur = bgn

            # append converted part
            scanned += s[:cur]

            # shrink s
            s = s[cur:]

        return scanned

    def _convert_one_function(self, s, fm, args, bgn, end):
        # no variable-length argument
        if (fm, len(args)) in self.translations:
            key = (fm, len(args))

            # x, y,... model arguments
            x_args = self.translations[key]['args']

            # make CORRESPONDENCES between model arguments and actual ones
            d = dict(zip(x_args, args))

        # with variable-length argument
        elif (fm, '*') in self.translations:
            key = (fm, '*')

            # x, y,..*args (model arguments)
            x_args = self.translations[key]['args']

            # make CORRESPONDENCES between model arguments and actual ones
            d = {}
            for i, x in enumerate(x_args):
                if x[0] == '*':
                    d[x] = ','.join(args[i:])
                    break
                d[x] = args[i]

        # out of self.translations
        else:
            err = "'{f}' is out of the whitelist.".format(f=fm)
            raise ValueError(err)

        # template string of converted function
        template = self.translations[key]['fs']

        # regex pattern for x_args
        pat = self.translations[key]['pat']

        scanned = ''
        cur = 0
        while True:
            m = pat.search(template)

            if m is None:
                scanned += template
                break

            # get model argument
            x = m.group()

            # get a start position of the model argument
            xbgn = m.start()

            # add the corresponding actual argument
            scanned += template[:xbgn] + d[x]

            # update cursor to the end of the model argument
            cur = m.end()

            # shrink template
            template = template[cur:]

        # update to swapped string
        s = s[:bgn] + scanned + s[end:]

        return s

    @classmethod
    def _get_args(cls, m):
        '''Get arguments of a Mathematica function'''

        s = m.string                # whole string
        anc = m.end() + 1           # pointing the first letter of arguments
        square, curly = [], []      # stack for brackets
        args = []

        # current cursor
        cur = anc
        for i, c in enumerate(s[anc:], anc):
            # extract one argument
            if c == ',' and (not square) and (not curly):
                args.append(s[cur:i])       # add an argument
                cur = i + 1                 # move cursor

            # handle list or matrix (for future usage)
            if c == '{':
                curly.append(c)
            elif c == '}':
                curly.pop()

            # seek corresponding ']' with skipping irrevant ones
            if c == '[':
                square.append(c)
            elif c == ']':
                if square:
                    square.pop()
                else:   # empty stack
                    args.append(s[cur:i])
                    break

        # the next position to ']' bracket (the function end)
        func_end = i + 1

        return args, func_end

    @classmethod
    def _replace(cls, s, bef):
        aft = cls.REPLACEMENTS[bef]
        s = s.replace(bef, aft)
        return s

    @classmethod
    def _apply_rules(cls, s, bef):
        pat, aft = cls.RULES[bef]
        return pat.sub(aft, s)

    @classmethod
    def _check_input(cls, s):
        for bracket in (('[', ']'), ('{', '}'), ('(', ')')):
            if s.count(bracket[0]) != s.count(bracket[1]):
                err = "'{f}' function form is invalid.".format(f=s)
                raise ValueError(err)

        if '{' in s:
            err = "Currently list is not supported."
            raise ValueError(err)

    def _parse_old(self, s):
        # input check
        self._check_input(s)

        # uncover '*' hiding behind a whitespace
        s = self._apply_rules(s, 'whitespace')

        # remove whitespace(s)
        s = self._replace(s, ' ')

        # add omitted '*' character
        s = self._apply_rules(s, 'add*_1')
        s = self._apply_rules(s, 'add*_2')

        # translate function
        s = self._convert_function(s)

        # '^' to '**'
        s = self._replace(s, '^')

        # 'Pi' to 'pi'
        s = self._apply_rules(s, 'Pi')

        # '{', '}' to '[', ']', respectively
#        s = cls._replace(s, '{')   # currently list is not taken into account
#        s = cls._replace(s, '}')

        return s

    def parse(self, s):
        s2 = named_characters_to_unicode(s)
        s3 = self._from_mathematica_to_tokens(s2)
        s4 = self._from_tokens_to_fullformlist(s3)
        s5 = self._from_fullformlist_to_sympy(s4)
        return s5

    INFIX = "Infix"
    PREFIX = "Prefix"
    POSTFIX = "Postfix"
    FLAT = "Flat"
    RIGHT = "Right"
    LEFT = "Left"

    _mathematica_op_precedence: list[tuple[str, str | None, dict[str, str | Callable]]] = [
        (POSTFIX, None, {";": lambda x: x + ["Null"] if isinstance(x, list) and x and x[0] == "CompoundExpression" else ["CompoundExpression", x, "Null"]}),
        (INFIX, FLAT, {";": "CompoundExpression"}),
        (INFIX, RIGHT, {"=": "Set", ":=": "SetDelayed", "+=": "AddTo", "-=": "SubtractFrom", "*=": "TimesBy", "/=": "DivideBy"}),
        (INFIX, LEFT, {"//": lambda x, y: [x, y]}),
        (POSTFIX, None, {"&": "Function"}),
        (INFIX, LEFT, {"/.": "ReplaceAll"}),
        (INFIX, RIGHT, {"->": "Rule", ":>": "RuleDelayed"}),
        (INFIX, LEFT, {"/;": "Condition"}),
        (INFIX, FLAT, {"|": "Alternatives"}),
        (POSTFIX, None, {"..": "Repeated", "...": "RepeatedNull"}),
        (INFIX, FLAT, {"||": "Or"}),
        (INFIX, FLAT, {"&&": "And"}),
        (PREFIX, None, {"!": "Not"}),
        (INFIX, FLAT, {"===": "SameQ", "=!=": "UnsameQ"}),
        (INFIX, FLAT, {"==": "Equal", "!=": "Unequal", "<=": "LessEqual", "<": "Less", ">=": "GreaterEqual", ">": "Greater"}),
        (INFIX, None, {";;": "Span"}),
        (INFIX, FLAT, {"+": "Plus", "-": "Plus"}),
        (INFIX, FLAT, {"*": "Times", "/": "Times"}),
        (INFIX, FLAT, {".": "Dot"}),
        (PREFIX, None, {"-": lambda x: MathematicaParser._get_neg(x),
                        "+": lambda x: x}),
        (INFIX, RIGHT, {"^": "Power"}),
        (INFIX, RIGHT, {"@@": "Apply", "/@": "Map", "//@": "MapAll", "@@@": lambda x, y: ["Apply", x, y, ["List", "1"]]}),
        (POSTFIX, None, {"'": "Derivative", "!": "Factorial", "!!": "Factorial2", "--": "Decrement"}),
        (INFIX, None, {"[": lambda x, y: [x, *y], "[[": lambda x, y: ["Part", x, *y]}),
        (PREFIX, None, {"{": lambda x: ["List", *x], "(": lambda x: x[0]}),
        (INFIX, None, {"?": "PatternTest"}),
        (POSTFIX, None, {
            "_": lambda x: ["Pattern", x, ["Blank"]],
            "_.": lambda x: ["Optional", ["Pattern", x, ["Blank"]]],
            "__": lambda x: ["Pattern", x, ["BlankSequence"]],
            "___": lambda x: ["Pattern", x, ["BlankNullSequence"]],
        }),
        (INFIX, None, {"_": lambda x, y: ["Pattern", x, ["Blank", y]]}),
        (PREFIX, None, {"#": "Slot", "##": "SlotSequence"}),
    ]

    _missing_arguments_default = {
        "#": lambda: ["Slot", "1"],
        "##": lambda: ["SlotSequence", "1"],
    }

    # This regex matches any valid python identifier — excluding
    # underscores, which Mathematica uses to denote patterns, and
    # therefore can't be part of a variable name.  The regex has the
    # form "[a][b]*", where `a` is the set of characters that can
    # start an identifier, and `b` is the set of characters that can
    # continue an identifier, which may also include numbers and
    # unicode combining characters.
    _literal = (
        "["
        + "".join(c for c in map(chr, range(sys.maxunicode+1)) if c!="_" and c.isidentifier())
        + "]["
        + "".join(c for c in map(chr, range(sys.maxunicode+1)) if c!="_" and ("x"+c).isidentifier())
        + "]*"
    )

    _number = r"(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)"

    _enclosure_open = ["(", "[", "[[", "{"]
    _enclosure_close = [")", "]", "]]", "}"]

    @classmethod
    def _get_neg(cls, x):
        return f"-{x}" if isinstance(x, str) and re.match(MathematicaParser._number, x) else ["Times", "-1", x]

    @classmethod
    def _get_inv(cls, x):
        return ["Power", x, "-1"]

    _regex_tokenizer = None

    def _get_tokenizer(self):
        if self._regex_tokenizer is not None:
            # Check if the regular expression has already been compiled:
            return self._regex_tokenizer
        tokens = [self._literal, self._number]
        tokens_escape = self._enclosure_open[:] + self._enclosure_close[:]
        for typ, strat, symdict in self._mathematica_op_precedence:
            for k in symdict:
                tokens_escape.append(k)
        tokens_escape.sort(key=lambda x: -len(x))
        tokens.extend(map(re.escape, tokens_escape))
        tokens.append(",")
        tokens.append("\n")
        tokenizer = re.compile("(" + "|".join(tokens) + ")")
        self._regex_tokenizer = tokenizer
        return self._regex_tokenizer

    def _from_mathematica_to_tokens(self, code: str):
        tokenizer = self._get_tokenizer()

        # Find strings:
        code_splits: list[str | list] = []
        while True:
            string_start = code.find("\"")
            if string_start == -1:
                if len(code) > 0:
                    code_splits.append(code)
                break
            match_end = re.search(r'(?<!\\)"', code[string_start+1:])
            if match_end is None:
                raise SyntaxError('mismatch in string "  " expression')
            string_end = string_start + match_end.start() + 1
            if string_start > 0:
                code_splits.append(code[:string_start])
            code_splits.append(["_Str", code[string_start+1:string_end].replace('\\"', '"')])
            code = code[string_end+1:]

        # Remove comments:
        for i, code_split in enumerate(code_splits):
            if isinstance(code_split, list):
                continue
            while True:
                pos_comment_start = code_split.find("(*")
                if pos_comment_start == -1:
                    break
                pos_comment_end = code_split.find("*)")
                if pos_comment_end == -1 or pos_comment_end < pos_comment_start:
                    raise SyntaxError("mismatch in comment (*  *) code")
                code_split = code_split[:pos_comment_start] + code_split[pos_comment_end+2:]
            code_splits[i] = code_split

        # Tokenize the input strings with a regular expression:
        def token_split(code):
            if isinstance(code, str):
                m = tokenizer.findall(code)
                if m or code.isascii():
                    return m
            return [code]

        token_lists = [token_split(code) for code in code_splits]
        tokens = [j for i in token_lists for j in i]
        # Remove newlines at the beginning
        while tokens and tokens[0] == "\n":
            tokens.pop(0)
        # Remove newlines at the end
        while tokens and tokens[-1] == "\n":
            tokens.pop(-1)

        return tokens

    def _is_op(self, token: str | list) -> bool:
        if isinstance(token, list):
            return False
        if re.match(self._literal, token):
            return False
        if re.match("-?" + self._number, token):
            return False
        return True

    def _is_valid_star1(self, token: str | list) -> bool:
        if token in (")", "}"):
            return True
        return not self._is_op(token)

    def _is_valid_star2(self, token: str | list) -> bool:
        if token in ("(", "{"):
            return True
        return not self._is_op(token)

    def _from_tokens_to_fullformlist(self, tokens: list):
        stack: list[list] = [[]]
        open_seq = []
        pointer: int = 0
        while pointer < len(tokens):
            token = tokens[pointer]
            if token in self._enclosure_open:
                stack[-1].append(token)
                open_seq.append(token)
                stack.append([])
            elif token == ",":
                if len(stack[-1]) == 0 and stack[-2][-1] == open_seq[-1]:
                    raise SyntaxError("%s cannot be followed by comma ," % open_seq[-1])
                stack[-1] = self._parse_after_braces(stack[-1])
                stack.append([])
            elif token in self._enclosure_close:
                ind = self._enclosure_close.index(token)
                if self._enclosure_open[ind] != open_seq[-1]:
                    unmatched_enclosure = SyntaxError("unmatched enclosure")
                    if token == "]]" and open_seq[-1] == "[":
                        if open_seq[-2] == "[":
                            # These two lines would be logically correct, but are
                            # unnecessary:
                            # token = "]"
                            # tokens[pointer] = "]"
                            tokens.insert(pointer+1, "]")
                        elif open_seq[-2] == "[[":
                            if tokens[pointer+1] == "]":
                                tokens[pointer+1] = "]]"
                            elif tokens[pointer+1] == "]]":
                                tokens[pointer+1] = "]]"
                                tokens.insert(pointer+2, "]")
                            else:
                                raise unmatched_enclosure
                    else:
                        raise unmatched_enclosure
                if len(stack[-1]) == 0 and stack[-2][-1] == "(":
                    raise SyntaxError("( ) not valid syntax")
                last_stack = self._parse_after_braces(stack[-1], True)
                stack[-1] = last_stack
                new_stack_element = []
                while stack[-1][-1] != open_seq[-1]:
                    new_stack_element.append(stack.pop())
                new_stack_element.reverse()
                if open_seq[-1] == "(" and len(new_stack_element) != 1:
                    raise SyntaxError("( must be followed by one expression, %i detected" % len(new_stack_element))
                stack[-1].append(new_stack_element)
                open_seq.pop(-1)
            else:
                stack[-1].append(token)
            pointer += 1
        if len(stack) != 1:
            raise RuntimeError("Stack should have only one element")
        return self._parse_after_braces(stack[0])

    def _util_remove_newlines(self, lines: list, tokens: list, inside_enclosure: bool):
        pointer = 0
        size = len(tokens)
        while pointer < size:
            token = tokens[pointer]
            if token == "\n":
                if inside_enclosure:
                    # Ignore newlines inside enclosures
                    tokens.pop(pointer)
                    size -= 1
                    continue
                if pointer == 0:
                    tokens.pop(0)
                    size -= 1
                    continue
                if pointer > 1:
                    try:
                        prev_expr = self._parse_after_braces(tokens[:pointer], inside_enclosure)
                    except SyntaxError:
                        tokens.pop(pointer)
                        size -= 1
                        continue
                else:
                    prev_expr = tokens[0]
                if len(prev_expr) > 0 and prev_expr[0] == "CompoundExpression":
                    lines.extend(prev_expr[1:])
                else:
                    lines.append(prev_expr)
                for i in range(pointer):
                    tokens.pop(0)
                size -= pointer
                pointer = 0
                continue
            pointer += 1

    def _util_add_missing_asterisks(self, tokens: list):
        size: int = len(tokens)
        pointer: int = 0
        while pointer < size:
            if (pointer > 0 and
                    self._is_valid_star1(tokens[pointer - 1]) and
                    self._is_valid_star2(tokens[pointer])):
                # This is a trick to add missing * operators in the expression,
                # `"*" in op_dict` makes sure the precedence level is the same as "*",
                # while `not self._is_op( ... )` makes sure this and the previous
                # expression are not operators.
                if tokens[pointer] == "(":
                    # ( has already been processed by now, replace:
                    tokens[pointer] = "*"
                    tokens[pointer + 1] = tokens[pointer + 1][0]
                else:
                    tokens.insert(pointer, "*")
                    pointer += 1
                    size += 1
            pointer += 1

    def _parse_after_braces(self, tokens: list, inside_enclosure: bool = False):
        op_dict: dict
        changed: bool = False
        lines: list = []

        self._util_remove_newlines(lines, tokens, inside_enclosure)

        for op_type, grouping_strat, op_dict in reversed(self._mathematica_op_precedence):
            if "*" in op_dict:
                self._util_add_missing_asterisks(tokens)
            size: int = len(tokens)
            pointer: int = 0
            while pointer < size:
                token = tokens[pointer]
                if isinstance(token, str) and token in op_dict:
                    op_name: str | Callable = op_dict[token]
                    node: list
                    first_index: int
                    if isinstance(op_name, str):
                        node = [op_name]
                        first_index = 1
                    else:
                        node = []
                        first_index = 0
                    if token in ("+", "-") and op_type == self.PREFIX and pointer > 0 and not self._is_op(tokens[pointer - 1]):
                        # Make sure that PREFIX + - don't match expressions like a + b or a - b,
                        # the INFIX + - are supposed to match that expression:
                        pointer += 1
                        continue
                    if op_type == self.INFIX:
                        if pointer == 0 or pointer == size - 1 or self._is_op(tokens[pointer - 1]) or self._is_op(tokens[pointer + 1]):
                            pointer += 1
                            continue
                    changed = True
                    tokens[pointer] = node
                    if op_type == self.INFIX:
                        arg1 = tokens.pop(pointer-1)
                        arg2 = tokens.pop(pointer)
                        if token == "/":
                            arg2 = self._get_inv(arg2)
                        elif token == "-":
                            arg2 = self._get_neg(arg2)
                        pointer -= 1
                        size -= 2
                        node.append(arg1)
                        node_p = node
                        if grouping_strat == self.FLAT:
                            while pointer + 2 < size and self._check_op_compatible(tokens[pointer+1], token):
                                node_p.append(arg2)
                                other_op = tokens.pop(pointer+1)
                                arg2 = tokens.pop(pointer+1)
                                if other_op == "/":
                                    arg2 = self._get_inv(arg2)
                                elif other_op == "-":
                                    arg2 = self._get_neg(arg2)
                                size -= 2
                            node_p.append(arg2)
                        elif grouping_strat == self.RIGHT:
                            while pointer + 2 < size and tokens[pointer+1] == token:
                                node_p.append([op_name, arg2])
                                node_p = node_p[-1]
                                tokens.pop(pointer+1)
                                arg2 = tokens.pop(pointer+1)
                                size -= 2
                            node_p.append(arg2)
                        elif grouping_strat == self.LEFT:
                            while pointer + 1 < size and tokens[pointer+1] == token:
                                if isinstance(op_name, str):
                                    node_p[first_index] = [op_name, node_p[first_index], arg2]
                                else:
                                    node_p[first_index] = op_name(node_p[first_index], arg2)
                                tokens.pop(pointer+1)
                                arg2 = tokens.pop(pointer+1)
                                size -= 2
                            node_p.append(arg2)
                        else:
                            node.append(arg2)
                    elif op_type == self.PREFIX:
                        if grouping_strat is not None:
                            raise TypeError("'Prefix' op_type should not have a grouping strat")
                        if pointer == size - 1 or self._is_op(tokens[pointer + 1]):
                            tokens[pointer] = self._missing_arguments_default[token]()
                        else:
                            node.append(tokens.pop(pointer+1))
                            size -= 1
                    elif op_type == self.POSTFIX:
                        if grouping_strat is not None:
                            raise TypeError("'Prefix' op_type should not have a grouping strat")
                        if pointer == 0 or self._is_op(tokens[pointer - 1]):
                            tokens[pointer] = self._missing_arguments_default[token]()
                        else:
                            node.append(tokens.pop(pointer-1))
                            pointer -= 1
                            size -= 1
                    if isinstance(op_name, Callable):  # type: ignore
                        op_call: Callable = typing.cast(Callable, op_name)
                        new_node = op_call(*node)
                        node.clear()
                        if isinstance(new_node, list):
                            node.extend(new_node)
                        else:
                            tokens[pointer] = new_node
                pointer += 1
        if len(tokens) > 1 or (len(lines) == 0 and len(tokens) == 0):
            if changed:
                # Trick to deal with cases in which an operator with lower
                # precedence should be transformed before an operator of higher
                # precedence. Such as in the case of `#&[x]` (that is
                # equivalent to `Lambda(d_, d_)(x)` in SymPy). In this case the
                # operator `&` has lower precedence than `[`, but needs to be
                # evaluated first because otherwise `# (&[x])` is not a valid
                # expression:
                return self._parse_after_braces(tokens, inside_enclosure)
            raise SyntaxError("unable to create a single AST for the expression")
        if len(lines) > 0:
            if tokens[0] and tokens[0][0] == "CompoundExpression":
                tokens = tokens[0][1:]
            compound_expression = ["CompoundExpression", *lines, *tokens]
            return compound_expression
        return tokens[0]

    def _check_op_compatible(self, op1: str, op2: str):
        if op1 == op2:
            return True
        muldiv = {"*", "/"}
        addsub = {"+", "-"}
        if op1 in muldiv and op2 in muldiv:
            return True
        if op1 in addsub and op2 in addsub:
            return True
        return False

    def _from_fullform_to_fullformlist(self, wmexpr: str):
        """
        Parses FullForm[Downvalues[]] generated by Mathematica
        """
        out: list = []
        stack = [out]
        generator = re.finditer(r'[\[\],]', wmexpr)
        last_pos = 0
        for match in generator:
            if match is None:
                break
            position = match.start()
            last_expr = wmexpr[last_pos:position].replace(',', '').replace(']', '').replace('[', '').strip()

            if match.group() == ',':
                if last_expr != '':
                    stack[-1].append(last_expr)
            elif match.group() == ']':
                if last_expr != '':
                    stack[-1].append(last_expr)
                stack.pop()
            elif match.group() == '[':
                stack[-1].append([last_expr])
                stack.append(stack[-1][-1])
            last_pos = match.end()
        return out[0]

    def _from_fullformlist_to_fullformsympy(self, pylist: list):
        from sympy import Function, Symbol

        def converter(expr):
            if isinstance(expr, list):
                if len(expr) > 0:
                    head = expr[0]
                    args = [converter(arg) for arg in expr[1:]]
                    return Function(head)(*args)
                else:
                    raise ValueError("Empty list of expressions")
            elif isinstance(expr, str):
                return Symbol(expr)
            else:
                return _sympify(expr)

        return converter(pylist)

    _node_conversions = {
        "Times": Mul,
        "Plus": Add,
        "Power": Pow,
        "Rational": Rational,
        "Log": lambda *a: log(*reversed(a)),
        "Log2": lambda x: log(x, 2),
        "Log10": lambda x: log(x, 10),
        "Exp": exp,
        "Sqrt": sqrt,

        "Sin": sin,
        "Cos": cos,
        "Tan": tan,
        "Cot": cot,
        "Sec": sec,
        "Csc": csc,

        "ArcSin": asin,
        "ArcCos": acos,
        "ArcTan": lambda *a: atan2(*reversed(a)) if len(a) == 2 else atan(*a),
        "ArcCot": acot,
        "ArcSec": asec,
        "ArcCsc": acsc,

        "Sinh": sinh,
        "Cosh": cosh,
        "Tanh": tanh,
        "Coth": coth,
        "Sech": sech,
        "Csch": csch,

        "ArcSinh": asinh,
        "ArcCosh": acosh,
        "ArcTanh": atanh,
        "ArcCoth": acoth,
        "ArcSech": asech,
        "ArcCsch": acsch,

        "Expand": expand,
        "Im": im,
        "Re": sympy.re,
        "Flatten": flatten,
        "Polylog": polylog,
        "Cancel": cancel,
        # Gamma=gamma,
        "TrigExpand": expand_trig,
        "Sign": sign,
        "Simplify": simplify,
        "Defer": UnevaluatedExpr,
        "Identity": S,
        # Sum=Sum_doit,
        # Module=With,
        # Block=With,
        "Null": lambda *a: S.Zero,
        "Mod": Mod,
        "Max": Max,
        "Min": Min,
        "Pochhammer": rf,
        "ExpIntegralEi": Ei,
        "SinIntegral": Si,
        "CosIntegral": Ci,
        "AiryAi": airyai,
        "AiryAiPrime": airyaiprime,
        "AiryBi": airybi,
        "AiryBiPrime": airybiprime,
        "LogIntegral": li,
        "PrimePi": primepi,
        "Prime": prime,
        "PrimeQ": isprime,

        "List": Tuple,
        "Greater": StrictGreaterThan,
        "GreaterEqual": GreaterThan,
        "Less": StrictLessThan,
        "LessEqual": LessThan,
        "Equal": Equality,
        "Or": Or,
        "And": And,

        "Function": _parse_Function,
    }

    _atom_conversions = {
        "I": I,
        "Pi": pi,
        "ExponentialE": E,
        "ImaginaryI": I,
        "ImaginaryJ": I,
    }

    def _from_fullformlist_to_sympy(self, full_form_list):

        def recurse(expr):
            if isinstance(expr, list):
                if isinstance(expr[0], list):
                    head = recurse(expr[0])
                else:
                    head = self._node_conversions.get(expr[0], Function(expr[0]))
                return head(*[recurse(arg) for arg in expr[1:]])
            else:
                return self._atom_conversions.get(expr, sympify(expr))

        return recurse(full_form_list)

    def _from_fullformsympy_to_sympy(self, mform):

        expr = mform
        for mma_form, sympy_node in self._node_conversions.items():
            expr = expr.replace(Function(mma_form), sympy_node)
        return expr


_named_characters = {
    "AAcute":"á",
    "ABar":"ā",
    "ACup":"ă",
    "ADoubleDot":"ä",
    "AE":"æ",
    "AGrave":"à",
    "AHat":"â",
    "Aleph":"ℵ",
    "AliasDelimiter":"AliasDelimiter",
    "AliasIndicator":"AliasIndicator",
    "AlignmentMarker":"AlignmentMarker",
    "Alpha":"α",
    "AltKey":"AltKey",
    "And":"∧",
    "Angle":"∠",
    "Angstrom":"Å",
    "AquariusSign":"♒",
    "AriesSign":"♈",
    "ARing":"å",
    "AscendingEllipsis":"⋰",
    "ATilde":"ã",
    "AutoLeftMatch":"AutoLeftMatch",
    "AutoOperand":"AutoOperand",
    "AutoPlaceholder":"AutoPlaceholder",
    "AutoRightMatch":"AutoRightMatch",
    "AutoSpace":"AutoSpace",
    "Backslash":"∖",
    "BeamedEighthNote":"♫",
    "BeamedSixteenthNote":"♬",
    "Because":"∵",
    "Bet":"ℶ",
    "Beta":"β",
    "BlackBishop":"♝",
    "BlackKing":"♚",
    "BlackKnight":"♞",
    "BlackPawn":"♟",
    "BlackQueen":"♛",
    "BlackRook":"♜",
    "Breve":"˘",
    "Bullet":"•",
    "CAcute":"ć",
    "CancerSign":"♋",
    "Cap":"⌢",
    "CapitalAAcute":"Á",
    "CapitalABar":"Ā",
    "CapitalACup":"Ă",
    "CapitalADoubleDot":"Ä",
    "CapitalAE":"Æ",
    "CapitalAGrave":"À",
    "CapitalAHat":"Â",
    "CapitalAlpha":"Α",
    "CapitalARing":"Å",
    "CapitalATilde":"Ã",
    "CapitalBeta":"Β",
    "CapitalCAcute":"Ć",
    "CapitalCCedilla":"Ç",
    "CapitalCHacek":"Č",
    "CapitalChi":"Χ",
    "CapitalDelta":"Δ",
    "CapitalDHacek":"Ď",
    "CapitalDifferentialD":"CapitalDifferentialD",
    "CapitalDigamma":"Ϝ",
    "CapitalEAcute":"É",
    "CapitalEBar":"Ē",
    "CapitalECup":"Ĕ",
    "CapitalEDoubleDot":"Ë",
    "CapitalEGrave":"È",
    "CapitalEHacek":"Ě",
    "CapitalEHat":"Ê",
    "CapitalEpsilon":"Ε",
    "CapitalEta":"Η",
    "CapitalEth":"Ð",
    "CapitalGamma":"Γ",
    "CapitalIAcute":"Í",
    "CapitalICup":"Ĭ",
    "CapitalIDoubleDot":"Ï",
    "CapitalIGrave":"Ì",
    "CapitalIHat":"Î",
    "CapitalIota":"Ι",
    "CapitalKappa":"Κ",
    "CapitalKoppa":"Ϟ",
    "CapitalLambda":"Λ",
    "CapitalLSlash":"Ł",
    "CapitalMu":"Μ",
    "CapitalNHacek":"Ň",
    "CapitalNTilde":"Ñ",
    "CapitalNu":"Ν",
    "CapitalOAcute":"Ó",
    "CapitalODoubleAcute":"Ő",
    "CapitalODoubleDot":"Ö",
    "CapitalOE":"Œ",
    "CapitalOGrave":"Ò",
    "CapitalOHat":"Ô",
    "CapitalOmega":"Ω",
    "CapitalOmicron":"Ο",
    "CapitalOSlash":"Ø",
    "CapitalOTilde":"Õ",
    "CapitalPhi":"Φ",
    "CapitalPi":"Π",
    "CapitalPsi":"Ψ",
    "CapitalRHacek":"Ř",
    "CapitalRho":"Ρ",
    "CapitalSampi":"Ϡ",
    "CapitalSHacek":"Š",
    "CapitalSigma":"Σ",
    "CapitalStigma":"Ϛ",
    "CapitalTau":"Τ",
    "CapitalTHacek":"Ť",
    "CapitalTheta":"Θ",
    "CapitalThorn":"Þ",
    "CapitalUAcute":"Ú",
    "CapitalUDoubleAcute":"Ű",
    "CapitalUDoubleDot":"Ü",
    "CapitalUGrave":"Ù",
    "CapitalUHat":"Û",
    "CapitalUpsilon":"Υ",
    "CapitalURing":"Ů",
    "CapitalXi":"Ξ",
    "CapitalYAcute":"Ý",
    "CapitalZeta":"Ζ",
    "CapitalZHacek":"Ž",
    "CapricornSign":"♑",
    "CCedilla":"ç",
    "Cedilla":"¸",
    "CenterDot":"·",
    "CenterEllipsis":"⋯",
    "Cent":"¢",
    "CHacek":"č",
    "CheckedBox":"☒",
    "Checkmark":"✓",
    "Chi":"χ",
    "CircleDot":"⊙",
    "CircleMinus":"⊖",
    "CirclePlus":"⊕",
    "CircleTimes":"⊗",
    "ClockwiseContourIntegral":"∲",
    "CloseCurlyDoubleQuote":"”",
    "CloseCurlyQuote":"’",
    "CloverLeaf":"⌘",
    "ClubSuit":"♣",
    "Colon":"∶",
    "CommandKey":"CommandKey",
    "Conditioned":"Conditioned",
    "Congruent":"≡",
    "Conjugate":"Conjugate",
    "ConjugateTranspose":"ConjugateTranspose",
    "ConstantC":"ConstantC",
    "Continuation":"Continuation",
    "ContourIntegral":"∮",
    "ControlKey":"ControlKey",
    "Coproduct":"∐",
    "Copyright":"©",
    "CounterClockwiseContourIntegral":"∳",
    "Cross":"Cross",
    "CupCap":"≍",
    "Cup":"⌣",
    "CurlyCapitalUpsilon":"ϒ",
    "CurlyEpsilon":"ε",
    "CurlyKappa":"ϰ",
    "CurlyPhi":"φ",
    "CurlyPi":"ϖ",
    "CurlyRho":"ϱ",
    "CurlyTheta":"ϑ",
    "Currency":"¤",
    "Dagger":"†",
    "Dalet":"ℸ",
    "Dash":"–",
    "Degree":"°",
    "DeleteKey":"DeleteKey",
    "Del":"∇",
    "Delta":"δ",
    "DescendingEllipsis":"⋱",
    "DHacek":"ď",
    "Diameter":"⌀",
    "Diamond":"⋄",
    "DiamondSuit":"♢",
    "DifferenceDelta":"DifferenceDelta",
    "DifferentialD":"DifferentialD",
    "Digamma":"ϝ",
    "DirectedEdge":"DirectedEdge",
    "DiscreteRatio":"DiscreteRatio",
    "DiscreteShift":"DiscreteShift",
    "DiscretionaryHyphen":"­",
    "DiscretionaryLineSeparator":" ",
    "DiscretionaryPageBreakAbove":"DiscretionaryPageBreakAbove",
    "DiscretionaryPageBreakBelow":"DiscretionaryPageBreakBelow",
    "DiscretionaryParagraphSeparator":" ",
    "Distributed":"Distributed",
    "Divides":"∣",
    "Divide":"÷",
    "DotEqual":"≐",
    "DotlessI":"ı",
    "DotlessJ":"DotlessJ",
    "DottedSquare":"DottedSquare",
    "DoubleContourIntegral":"∯",
    "DoubleDagger":"‡",
    "DoubledGamma":"DoubledGamma",
    "DoubleDot":"¨",
    "DoubleDownArrow":"⇓",
    "DoubledPi":"DoubledPi",
    "DoubleLeftArrow":"⇐",
    "DoubleLeftRightArrow":"⇔",
    "DoubleLeftTee":"⫤",
    "DoubleLongLeftArrow":"⟸",
    "DoubleLongLeftRightArrow":"⟺",
    "DoubleLongRightArrow":"⟹",
    "DoublePrime":"″",
    "DoubleRightArrow":"⇒",
    "DoubleRightTee":"⊨",
    "DoubleStruckA":"𝕒",
    "DoubleStruckB":"𝕓",
    "DoubleStruckC":"𝕔",
    "DoubleStruckCapitalA":"𝔸",
    "DoubleStruckCapitalB":"𝔹",
    "DoubleStruckCapitalC":"ℂ",
    "DoubleStruckCapitalD":"𝔻",
    "DoubleStruckCapitalE":"𝔼",
    "DoubleStruckCapitalF":"𝔽",
    "DoubleStruckCapitalG":"𝔾",
    "DoubleStruckCapitalH":"ℍ",
    "DoubleStruckCapitalI":"𝕀",
    "DoubleStruckCapitalJ":"𝕁",
    "DoubleStruckCapitalK":"𝕂",
    "DoubleStruckCapitalL":"𝕃",
    "DoubleStruckCapitalM":"𝕄",
    "DoubleStruckCapitalN":"ℕ",
    "DoubleStruckCapitalO":"𝕆",
    "DoubleStruckCapitalP":"ℙ",
    "DoubleStruckCapitalQ":"ℚ",
    "DoubleStruckCapitalR":"ℝ",
    "DoubleStruckCapitalS":"𝕊",
    "DoubleStruckCapitalT":"𝕋",
    "DoubleStruckCapitalU":"𝕌",
    "DoubleStruckCapitalV":"𝕍",
    "DoubleStruckCapitalW":"𝕎",
    "DoubleStruckCapitalX":"𝕏",
    "DoubleStruckCapitalY":"𝕐",
    "DoubleStruckCapitalZ":"ℤ",
    "DoubleStruckD":"𝕕",
    "DoubleStruckE":"𝕖",
    "DoubleStruckEight":"DoubleStruckEight",
    "DoubleStruckF":"𝕗",
    "DoubleStruckFive":"DoubleStruckFive",
    "DoubleStruckFour":"DoubleStruckFour",
    "DoubleStruckG":"𝕘",
    "DoubleStruckH":"𝕙",
    "DoubleStruckI":"𝕚",
    "DoubleStruckJ":"𝕛",
    "DoubleStruckK":"𝕜",
    "DoubleStruckL":"𝕝",
    "DoubleStruckM":"𝕞",
    "DoubleStruckN":"𝕟",
    "DoubleStruckNine":"DoubleStruckNine",
    "DoubleStruckO":"𝕠",
    "DoubleStruckOne":"DoubleStruckOne",
    "DoubleStruckP":"𝕡",
    "DoubleStruckQ":"𝕢",
    "DoubleStruckR":"𝕣",
    "DoubleStruckS":"𝕤",
    "DoubleStruckSeven":"DoubleStruckSeven",
    "DoubleStruckSix":"DoubleStruckSix",
    "DoubleStruckT":"𝕥",
    "DoubleStruckThree":"DoubleStruckThree",
    "DoubleStruckTwo":"DoubleStruckTwo",
    "DoubleStruckU":"𝕦",
    "DoubleStruckV":"𝕧",
    "DoubleStruckW":"𝕨",
    "DoubleStruckX":"𝕩",
    "DoubleStruckY":"𝕪",
    "DoubleStruckZ":"𝕫",
    "DoubleStruckZero":"DoubleStruckZero",
    "DoubleUpArrow":"⇑",
    "DoubleUpDownArrow":"⇕",
    "DoubleVerticalBar":"∥",
    "DownArrowBar":"⤓",
    "DownArrow":"↓",
    "DownArrowUpArrow":"⇵",
    "DownBreve":"DownBreve",
    "DownExclamation":"¡",
    "DownLeftRightVector":"⥐",
    "DownLeftTeeVector":"⥞",
    "DownLeftVector":"↽",
    "DownLeftVectorBar":"⥖",
    "DownPointer":"▾",
    "DownQuestion":"¿",
    "DownRightTeeVector":"⥟",
    "DownRightVector":"⇁",
    "DownRightVectorBar":"⥗",
    "DownTeeArrow":"↧",
    "DownTee":"⊤",
    "EAcute":"é",
    "Earth":"Earth",
    "EBar":"ē",
    "ECup":"ĕ",
    "EDoubleDot":"ë",
    "EGrave":"è",
    "EHacek":"ě",
    "EHat":"ê",
    "EighthNote":"♪",
    "Element":"∈",
    "Ellipsis":"…",
    "EmptyCircle":"○",
    "EmptyDiamond":"◇",
    "EmptyDownTriangle":"▽",
    "EmptyRectangle":"▯",
    "EmptySet":"∅",
    "EmptySmallCircle":"◦",
    "EmptySmallSquare":"◻",
    "EmptySquare":"□",
    "EmptyUpTriangle":"△",
    "EmptyVerySmallSquare":"▫",
    "EnterKey":"EnterKey",
    "EntityEnd":"EntityEnd",
    "EntityStart":"EntityStart",
    "Epsilon":"ϵ",
    "Equal":"Equal",
    "EqualTilde":"≂",
    "Equilibrium":"⇌",
    "Equivalent":"⧦",
    "ErrorIndicator":"ErrorIndicator",
    "EscapeKey":"EscapeKey",
    "Eta":"η",
    "Eth":"ð",
    "Euro":"€",
    "Exists":"∃",
    "ExponentialE":"ExponentialE",
    "FiLigature":"ﬁ",
    "FilledCircle":"●",
    "FilledDiamond":"◆",
    "FilledDownTriangle":"▼",
    "FilledLeftTriangle":"◀",
    "FilledRectangle":"▮",
    "FilledRightTriangle":"▶",
    "FilledSmallCircle":"FilledSmallCircle",
    "FilledSmallSquare":"■",
    "FilledSquare":"■",
    "FilledUpTriangle":"▲",
    "FilledVerySmallSquare":"▪",
    "FinalSigma":"ς",
    "FirstPage":"FirstPage",
    "FivePointedStar":"★",
    "Flat":"♭",
    "FlLigature":"ﬂ",
    "Florin":"ƒ",
    "ForAll":"∀",
    "FormalA":"ạ",
    "FormalAlpha":"α̣",
    "FormalB":"ḅ",
    "FormalBeta":"β̣",
    "FormalC":"c̣",
    "FormalCapitalA":"Ạ",
    "FormalCapitalAlpha":"Α̣",
    "FormalCapitalB":"Ḅ",
    "FormalCapitalBeta":"Β̣",
    "FormalCapitalC":"C̣",
    "FormalCapitalChi":"Χ̣",
    "FormalCapitalD":"Ḍ",
    "FormalCapitalDelta":"Δ̣",
    "FormalCapitalDigamma":"Ϝ̣",
    "FormalCapitalE":"Ẹ",
    "FormalCapitalEpsilon":"Ε̣",
    "FormalCapitalEta":"Η̣",
    "FormalCapitalF":"F̣",
    "FormalCapitalG":"G̣",
    "FormalCapitalGamma":"Γ̣",
    "FormalCapitalH":"Ḥ",
    "FormalCapitalI":"Ị",
    "FormalCapitalIota":"Ι̣",
    "FormalCapitalJ":"J̣",
    "FormalCapitalK":"Ḳ",
    "FormalCapitalKappa":"Κ̣",
    "FormalCapitalKoppa":"Ϟ̣",
    "FormalCapitalL":"Ḷ",
    "FormalCapitalLambda":"Λ̣",
    "FormalCapitalM":"Ṃ",
    "FormalCapitalMu":"Μ̣",
    "FormalCapitalN":"Ṇ",
    "FormalCapitalNu":"Ν̣",
    "FormalCapitalO":"Ọ",
    "FormalCapitalOmega":"Ω̣",
    "FormalCapitalOmicron":"Ο̣",
    "FormalCapitalP":"P̣",
    "FormalCapitalPhi":"Φ̣",
    "FormalCapitalPi":"Π̣",
    "FormalCapitalPsi":"Ψ̣",
    "FormalCapitalQ":"Q̣",
    "FormalCapitalR":"Ṛ",
    "FormalCapitalRho":"Ρ̣",
    "FormalCapitalS":"Ṣ",
    "FormalCapitalSampi":"Ϡ̣",
    "FormalCapitalSigma":"Σ̣",
    "FormalCapitalStigma":"Ϛ̣",
    "FormalCapitalT":"Ṭ",
    "FormalCapitalTau":"Τ̣",
    "FormalCapitalTheta":"Θ̣",
    "FormalCapitalU":"Ụ",
    "FormalCapitalUpsilon":"Υ̣",
    "FormalCapitalV":"Ṿ",
    "FormalCapitalW":"Ẉ",
    "FormalCapitalX":"X̣",
    "FormalCapitalXi":"Ξ̣",
    "FormalCapitalY":"Ỵ",
    "FormalCapitalZ":"Ẓ",
    "FormalCapitalZeta":"Ζ̣",
    "FormalChi":"χ̣",
    "FormalCurlyCapitalUpsilon":"ϒ̣",
    "FormalCurlyEpsilon":"ε̣",
    "FormalCurlyKappa":"ϰ̣",
    "FormalCurlyPhi":"φ̣",
    "FormalCurlyPi":"ϖ̣",
    "FormalCurlyRho":"ϱ̣",
    "FormalCurlyTheta":"ϑ̣",
    "FormalD":"ḍ",
    "FormalDelta":"δ̣",
    "FormalDigamma":"ϝ̣",
    "FormalE":"ẹ",
    "FormalEpsilon":"ϵ̣",
    "FormalEta":"η̣",
    "FormalF":"f̣",
    "FormalFinalSigma":"ς̣",
    "FormalG":"g̣",
    "FormalGamma":"γ̣",
    "FormalH":"ḥ",
    "FormalI":"ị",
    "FormalIota":"ι̣",
    "FormalJ":"j̣",
    "FormalK":"ḳ",
    "FormalKappa":"κ̣",
    "FormalKoppa":"ϟ̣",
    "FormalL":"ḷ",
    "FormalLambda":"λ̣",
    "FormalM":"ṃ",
    "FormalMu":"μ̣",
    "FormalN":"ṇ",
    "FormalNu":"ν̣",
    "FormalO":"ọ",
    "FormalOmega":"ω̣",
    "FormalOmicron":"ο̣",
    "FormalP":"p̣",
    "FormalPhi":"ϕ̣",
    "FormalPi":"π̣",
    "FormalPsi":"ψ̣",
    "FormalQ":"q̣",
    "FormalR":"ṛ",
    "FormalRho":"ρ̣",
    "FormalS":"ṣ",
    "FormalSampi":"ϡ̣",
    "FormalScriptA":"𝒶̣",
    "FormalScriptB":"𝒷̣",
    "FormalScriptC":"𝒸̣",
    "FormalScriptCapitalA":"𝒜̣",
    "FormalScriptCapitalB":"ℬ̣",
    "FormalScriptCapitalC":"𝒞̣",
    "FormalScriptCapitalD":"𝒟̣",
    "FormalScriptCapitalE":"ℰ̣",
    "FormalScriptCapitalF":"ℱ̣",
    "FormalScriptCapitalG":"𝒢̣",
    "FormalScriptCapitalH":"ℋ̣",
    "FormalScriptCapitalI":"ℐ̣",
    "FormalScriptCapitalJ":"𝒥̣",
    "FormalScriptCapitalK":"𝒦̣",
    "FormalScriptCapitalL":"ℒ̣",
    "FormalScriptCapitalM":"ℳ̣",
    "FormalScriptCapitalN":"𝒩̣",
    "FormalScriptCapitalO":"𝒪̣",
    "FormalScriptCapitalP":"𝒫̣",
    "FormalScriptCapitalQ":"𝒬̣",
    "FormalScriptCapitalR":"ℛ̣",
    "FormalScriptCapitalS":"𝒮̣",
    "FormalScriptCapitalT":"𝒯̣",
    "FormalScriptCapitalU":"𝒰̣",
    "FormalScriptCapitalV":"𝒱̣",
    "FormalScriptCapitalW":"𝒲̣",
    "FormalScriptCapitalX":"𝒳̣",
    "FormalScriptCapitalY":"𝒴̣",
    "FormalScriptCapitalZ":"𝒵̣",
    "FormalScriptD":"𝒹̣",
    "FormalScriptE":"ℯ̣",
    "FormalScriptF":"𝒻̣",
    "FormalScriptG":"ℊ̣",
    "FormalScriptH":"𝒽̣",
    "FormalScriptI":"𝒾̣",
    "FormalScriptJ":"𝒿̣",
    "FormalScriptK":"𝓀̣",
    "FormalScriptL":"𝓁̣",
    "FormalScriptM":"𝓂̣",
    "FormalScriptN":"𝓃̣",
    "FormalScriptO":"ℴ̣",
    "FormalScriptP":"𝓅̣",
    "FormalScriptQ":"𝓆̣",
    "FormalScriptR":"𝓇̣",
    "FormalScriptS":"𝓈̣",
    "FormalScriptT":"𝓉̣",
    "FormalScriptU":"𝓊̣",
    "FormalScriptV":"𝓋̣",
    "FormalScriptW":"𝓌̣",
    "FormalScriptX":"𝓍̣",
    "FormalScriptY":"𝓎̣",
    "FormalScriptZ":"𝓏̣",
    "FormalSigma":"σ̣",
    "FormalStigma":"ϛ̣",
    "FormalT":"ṭ",
    "FormalTau":"τ̣",
    "FormalTheta":"θ̣",
    "FormalU":"ụ",
    "FormalUpsilon":"υ̣",
    "FormalV":"ṿ",
    "FormalW":"ẉ",
    "FormalX":"x̣",
    "FormalXi":"ξ̣",
    "FormalY":"ỵ",
    "FormalZ":"ẓ",
    "FormalZeta":"ζ̣",
    "FreakedSmiley":"FreakedSmiley",
    "Function":"Function",
    "Gamma":"γ",
    "GeminiSign":"♊",
    "Gimel":"ℷ",
    "GothicA":"𝔞",
    "GothicB":"𝔟",
    "GothicC":"𝔠",
    "GothicCapitalA":"𝔄",
    "GothicCapitalB":"𝔅",
    "GothicCapitalC":"ℭ",
    "GothicCapitalD":"𝔇",
    "GothicCapitalE":"𝔈",
    "GothicCapitalF":"𝔉",
    "GothicCapitalG":"𝔊",
    "GothicCapitalH":"ℌ",
    "GothicCapitalI":"ℑ",
    "GothicCapitalJ":"𝔍",
    "GothicCapitalK":"𝔎",
    "GothicCapitalL":"𝔏",
    "GothicCapitalM":"𝔐",
    "GothicCapitalN":"𝔑",
    "GothicCapitalO":"𝔒",
    "GothicCapitalP":"𝔓",
    "GothicCapitalQ":"𝔔",
    "GothicCapitalR":"ℜ",
    "GothicCapitalS":"𝔖",
    "GothicCapitalT":"𝔗",
    "GothicCapitalU":"𝔘",
    "GothicCapitalV":"𝔙",
    "GothicCapitalW":"𝔚",
    "GothicCapitalX":"𝔛",
    "GothicCapitalY":"𝔜",
    "GothicCapitalZ":"ℨ",
    "GothicD":"𝔡",
    "GothicE":"𝔢",
    "GothicEight":"GothicEight",
    "GothicF":"𝔣",
    "GothicFive":"GothicFive",
    "GothicFour":"GothicFour",
    "GothicG":"𝔤",
    "GothicH":"𝔥",
    "GothicI":"𝔦",
    "GothicJ":"𝔧",
    "GothicK":"𝔨",
    "GothicL":"𝔩",
    "GothicM":"𝔪",
    "GothicN":"𝔫",
    "GothicNine":"GothicNine",
    "GothicO":"𝔬",
    "GothicOne":"GothicOne",
    "GothicP":"𝔭",
    "GothicQ":"𝔮",
    "GothicR":"𝔯",
    "GothicS":"𝔰",
    "GothicSeven":"GothicSeven",
    "GothicSix":"GothicSix",
    "GothicT":"𝔱",
    "GothicThree":"GothicThree",
    "GothicTwo":"GothicTwo",
    "GothicU":"𝔲",
    "GothicV":"𝔳",
    "GothicW":"𝔴",
    "GothicX":"𝔵",
    "GothicY":"𝔶",
    "GothicZ":"𝔷",
    "GothicZero":"GothicZero",
    "GrayCircle":"●",
    "GraySquare":"■",
    "GreaterEqualLess":"⋛",
    "GreaterEqual":"≥",
    "GreaterFullEqual":"≧",
    "GreaterGreater":"≫",
    "GreaterLess":"≷",
    "GreaterSlantEqual":"⩾",
    "GreaterTilde":"≳",
    "Hacek":"ˇ",
    "HappySmiley":"☺",
    "HBar":"ℏ",
    "HeartSuit":"♡",
    "HermitianConjugate":"HermitianConjugate",
    "HorizontalLine":"─",
    "HumpDownHump":"≎",
    "HumpEqual":"≏",
    "Hyphen":"‐",
    "IAcute":"í",
    "ICup":"ĭ",
    "IDoubleDot":"ï",
    "IGrave":"ì",
    "IHat":"î",
    "ImaginaryI":"ImaginaryI",
    "ImaginaryJ":"ImaginaryJ",
    "ImplicitPlus":"ImplicitPlus",
    "Implies":"Implies",
    "IndentingNewLine":"IndentingNewLine",
    "Infinity":"∞",
    "Integral":"∫",
    "Intersection":"⋂",
    "InvisibleApplication":"\u2061",
    "InvisibleComma":"\u2063",
    "InvisiblePostfixScriptBase":"InvisiblePostfixScriptBase",
    "InvisiblePrefixScriptBase":"InvisiblePrefixScriptBase",
    "InvisibleSpace":"\u200B",
    "InvisibleTimes":"\u2062",
    "Iota":"ι",
    "Jupiter":"♃",
    "Kappa":"κ",
    "KernelIcon":"KernelIcon",
    "Koppa":"ϟ",
    "Lambda":"λ",
    "LastPage":"LastPage",
    "LeftAngleBracket":"〈",
    "LeftArrowBar":"⇤",
    "LeftArrow":"←",
    "LeftArrowRightArrow":"⇆",
    "LeftAssociation":"LeftAssociation",
    "LeftBracketingBar":"LeftBracketingBar",
    "LeftCeiling":"⌈",
    "LeftDoubleBracket":"〚",
    "LeftDoubleBracketingBar":"LeftDoubleBracketingBar",
    "LeftDownTeeVector":"⥡",
    "LeftDownVectorBar":"⥙",
    "LeftDownVector":"⇃",
    "LeftFloor":"⌊",
    "LeftGuillemet":"«",
    "LeftModified":"LeftModified",
    "LeftPointer":"◂",
    "LeftRightArrow":"↔",
    "LeftRightVector":"⥎",
    "LeftSkeleton":"LeftSkeleton",
    "LeftTee":"⊣",
    "LeftTeeArrow":"↤",
    "LeftTeeVector":"⥚",
    "LeftTriangle":"⊲",
    "LeftTriangleBar":"⧏",
    "LeftTriangleEqual":"⊴",
    "LeftUpDownVector":"⥑",
    "LeftUpTeeVector":"⥠",
    "LeftUpVector":"↿",
    "LeftUpVectorBar":"⥘",
    "LeftVector":"↼",
    "LeftVectorBar":"⥒",
    "LeoSign":"♌",
    "LessEqual":"≤",
    "LessEqualGreater":"⋚",
    "LessFullEqual":"≦",
    "LessGreater":"≶",
    "LessLess":"≪",
    "LessSlantEqual":"⩽",
    "LessTilde":"≲",
    "LetterSpace":"LetterSpace",
    "LibraSign":"♎",
    "LightBulb":"LightBulb",
    "LineSeparator":" ",
    "LongDash":"—",
    "LongEqual":"LongEqual",
    "LongLeftArrow":"⟵",
    "LongLeftRightArrow":"⟷",
    "LongRightArrow":"⟶",
    "LowerLeftArrow":"↙",
    "LowerRightArrow":"↘",
    "LSlash":"ł",
    "Mars":"♂",
    "MathematicaIcon":"MathematicaIcon",
    "MeasuredAngle":"∡",
    "MediumSpace":"\u2005",
    "Mercury":"☿",
    "Mho":"℧",
    "Micro":"µ",
    "MinusPlus":"∓",
    "Mu":"μ",
    "Nand":"⊼",
    "Natural":"♮",
    "NegativeMediumSpace":"NegativeMediumSpace",
    "NegativeThickSpace":"NegativeThickSpace",
    "NegativeThinSpace":"NegativeThinSpace",
    "NegativeVeryThinSpace":"NegativeVeryThinSpace",
    "Neptune":"♆",
    "NestedGreaterGreater":"⪢",
    "NestedLessLess":"⪡",
    "NeutralSmiley":"NeutralSmiley",
    "NewLine":"\n",
    "NHacek":"ň",
    "NoBreak":"NoBreak",
    "NonBreakingSpace":"\xa0",
    "Nor":"⊽",
    "NotCongruent":"≢",
    "NotCupCap":"≭",
    "NotDoubleVerticalBar":"∦",
    "NotElement":"∉",
    "NotEqual":"≠",
    "NotEqualTilde":"NotEqualTilde",
    "NotExists":"∄",
    "NotGreater":"≯",
    "NotGreaterEqual":"≱",
    "NotGreaterFullEqual":"≩",
    "NotGreaterGreater":"NotGreaterGreater",
    "NotGreaterLess":"≹",
    "NotGreaterSlantEqual":"NotGreaterSlantEqual",
    "NotGreaterTilde":"≵",
    "NotHumpDownHump":"NotHumpDownHump",
    "NotHumpEqual":"NotHumpEqual",
    "NotLeftTriangle":"⋪",
    "NotLeftTriangleBar":"NotLeftTriangleBar",
    "NotLeftTriangleEqual":"⋬",
    "NotLessEqual":"≰",
    "NotLessFullEqual":"≨",
    "NotLessGreater":"≸",
    "NotLess":"≮",
    "NotLessLess":"NotLessLess",
    "NotLessSlantEqual":"NotLessSlantEqual",
    "NotLessTilde":"≴",
    "Not":"¬",
    "NotNestedGreaterGreater":"NotNestedGreaterGreater",
    "NotNestedLessLess":"NotNestedLessLess",
    "NotPrecedes":"⊀",
    "NotPrecedesEqual":"NotPrecedesEqual",
    "NotPrecedesSlantEqual":"⋠",
    "NotPrecedesTilde":"⋨",
    "NotReverseElement":"∌",
    "NotRightTriangle":"⋫",
    "NotRightTriangleBar":"NotRightTriangleBar",
    "NotRightTriangleEqual":"⋭",
    "NotSquareSubset":"NotSquareSubset",
    "NotSquareSubsetEqual":"⋢",
    "NotSquareSuperset":"NotSquareSuperset",
    "NotSquareSupersetEqual":"⋣",
    "NotSubset":"⊄",
    "NotSubsetEqual":"⊈",
    "NotSucceeds":"⊁",
    "NotSucceedsEqual":"NotSucceedsEqual",
    "NotSucceedsSlantEqual":"⋡",
    "NotSucceedsTilde":"⋩",
    "NotSuperset":"⊅",
    "NotSupersetEqual":"⊉",
    "NotTilde":"≁",
    "NotTildeEqual":"≄",
    "NotTildeFullEqual":"≇",
    "NotTildeTilde":"≉",
    "NotVerticalBar":"NotVerticalBar",
    "NTilde":"ñ",
    "Nu":"ν",
    "NumberSign":"NumberSign",
    "OAcute":"ó",
    "ODoubleAcute":"ő",
    "ODoubleDot":"ö",
    "OE":"œ",
    "OGrave":"ò",
    "OHat":"ô",
    "Omega":"ω",
    "Omicron":"ο",
    "OpenCurlyDoubleQuote":"“",
    "OpenCurlyQuote":"‘",
    "OptionKey":"OptionKey",
    "Or":"∨",
    "OSlash":"ø",
    "OTilde":"õ",
    "OverBrace":"︷",
    "OverBracket":"⎴",
    "OverParenthesis":"︵",
    "Paragraph":"¶",
    "ParagraphSeparator":" ",
    "PartialD":"∂",
    "PermutationProduct":"PermutationProduct",
    "Perpendicular":"⟂",
    "Phi":"ϕ",
    "Pi":"π",
    "Piecewise":"Piecewise",
    "PiscesSign":"♓",
    "Placeholder":"Placeholder",
    "PlusMinus":"±",
    "Pluto":"♇",
    "Precedes":"≺",
    "PrecedesEqual":"⪯",
    "PrecedesSlantEqual":"≼",
    "PrecedesTilde":"≾",
    "Prime":"′",
    "Product":"∏",
    "Proportion":"∷",
    "Proportional":"∝",
    "Psi":"ψ",
    "QuarterNote":"♩",
    "RawAmpersand":"&",
    "RawAt":"@",
    "RawBackquote":"`",
    "RawBackslash":"\\",
    "RawColon":":",
    "RawComma":",",
    "RawDash":"-",
    "RawDollar":"$",
    "RawDot":".",
    "RawDoubleQuote":"\"",
    "RawEqual":"=",
    "RawExclamation":"!",
    "RawGreater":">",
    "RawLeftBrace":"{",
    "RawLeftBracket":"[",
    "RawLeftParenthesis":"(",
    "RawLess":"<",
    "RawNumberSign":"#",
    "RawPercent":"%",
    "RawPlus":"+",
    "RawQuestion":"?",
    "RawQuote":"'",
    "RawReturn":"RawReturn",
    "RawRightBrace":"}",
    "RawRightBracket":"]",
    "RawRightParenthesis":")",
    "RawSemicolon":";",
    "RawSlash":"/",
    "RawSpace":"RawSpace",
    "RawStar":"*",
    "RawTab":"RawTab",
    "RawTilde":"~",
    "RawUnderscore":"_",
    "RawVerticalBar":"|",
    "RawWedge":"^",
    "RegisteredTrademark":"®",
    "ReturnIndicator":"↵",
    "ReturnKey":"ReturnKey",
    "ReverseDoublePrime":"‶",
    "ReverseElement":"∋",
    "ReverseEquilibrium":"⇋",
    "ReversePrime":"‵",
    "ReverseUpEquilibrium":"⥯",
    "RHacek":"ř",
    "Rho":"ρ",
    "RightAngle":"∟",
    "RightAngleBracket":"〉",
    "RightArrow":"→",
    "RightArrowBar":"⇥",
    "RightArrowLeftArrow":"⇄",
    "RightAssociation":"RightAssociation",
    "RightBracketingBar":"RightBracketingBar",
    "RightCeiling":"⌉",
    "RightDoubleBracket":"〛",
    "RightDoubleBracketingBar":"RightDoubleBracketingBar",
    "RightDownTeeVector":"⥝",
    "RightDownVector":"⇂",
    "RightDownVectorBar":"⥕",
    "RightFloor":"⌋",
    "RightGuillemet":"»",
    "RightModified":"RightModified",
    "RightPointer":"▸",
    "RightSkeleton":"RightSkeleton",
    "RightTee":"⊢",
    "RightTeeArrow":"↦",
    "RightTeeVector":"⥛",
    "RightTriangle":"⊳",
    "RightTriangleBar":"⧐",
    "RightTriangleEqual":"⊵",
    "RightUpDownVector":"⥏",
    "RightUpTeeVector":"⥜",
    "RightUpVector":"↾",
    "RightUpVectorBar":"⥔",
    "RightVector":"⇀",
    "RightVectorBar":"⥓",
    "RoundImplies":"⥰",
    "RoundSpaceIndicator":"RoundSpaceIndicator",
    "Rule":"Rule",
    "RuleDelayed":"RuleDelayed",
    "SadSmiley":"☹",
    "SagittariusSign":"♐",
    "Sampi":"ϡ",
    "Saturn":"♄",
    "ScorpioSign":"♏",
    "ScriptA":"𝒶",
    "ScriptB":"𝒷",
    "ScriptC":"𝒸",
    "ScriptCapitalA":"𝒜",
    "ScriptCapitalB":"ℬ",
    "ScriptCapitalC":"𝒞",
    "ScriptCapitalD":"𝒟",
    "ScriptCapitalE":"ℰ",
    "ScriptCapitalF":"ℱ",
    "ScriptCapitalG":"𝒢",
    "ScriptCapitalH":"ℋ",
    "ScriptCapitalI":"ℐ",
    "ScriptCapitalJ":"𝒥",
    "ScriptCapitalK":"𝒦",
    "ScriptCapitalL":"ℒ",
    "ScriptCapitalM":"ℳ",
    "ScriptCapitalN":"𝒩",
    "ScriptCapitalO":"𝒪",
    "ScriptCapitalP":"𝒫",
    "ScriptCapitalQ":"𝒬",
    "ScriptCapitalR":"ℛ",
    "ScriptCapitalS":"𝒮",
    "ScriptCapitalT":"𝒯",
    "ScriptCapitalU":"𝒰",
    "ScriptCapitalV":"𝒱",
    "ScriptCapitalW":"𝒲",
    "ScriptCapitalX":"𝒳",
    "ScriptCapitalY":"𝒴",
    "ScriptCapitalZ":"𝒵",
    "ScriptD":"𝒹",
    "ScriptDotlessI":"ScriptDotlessI",
    "ScriptDotlessJ":"ScriptDotlessJ",
    "ScriptE":"ℯ",
    "ScriptEight":"ScriptEight",
    "ScriptF":"𝒻",
    "ScriptFive":"ScriptFive",
    "ScriptFour":"ScriptFour",
    "ScriptG":"ℊ",
    "ScriptH":"𝒽",
    "ScriptI":"𝒾",
    "ScriptJ":"𝒿",
    "ScriptK":"𝓀",
    "ScriptL":"ℓ",
    "ScriptM":"𝓂",
    "ScriptN":"𝓃",
    "ScriptNine":"ScriptNine",
    "ScriptO":"ℴ",
    "ScriptOne":"ScriptOne",
    "ScriptP":"𝓅",
    "ScriptQ":"𝓆",
    "ScriptR":"𝓇",
    "ScriptS":"𝓈",
    "ScriptSeven":"ScriptSeven",
    "ScriptSix":"ScriptSix",
    "ScriptT":"𝓉",
    "ScriptThree":"ScriptThree",
    "ScriptTwo":"ScriptTwo",
    "ScriptU":"𝓊",
    "ScriptV":"𝓋",
    "ScriptW":"𝓌",
    "ScriptX":"𝓍",
    "ScriptY":"𝓎",
    "ScriptZ":"𝓏",
    "ScriptZero":"ScriptZero",
    "Section":"§",
    "SelectionPlaceholder":"SelectionPlaceholder",
    "SHacek":"š",
    "Sharp":"♯",
    "ShortDownArrow":"ShortDownArrow",
    "ShortLeftArrow":"ShortLeftArrow",
    "ShortRightArrow":"ShortRightArrow",
    "ShortUpArrow":"ShortUpArrow",
    "Sigma":"σ",
    "SixPointedStar":"✶",
    "SkeletonIndicator":"⁃",
    "SmallCircle":"∘",
    "SpaceIndicator":"␣",
    "SpaceKey":"SpaceKey",
    "SpadeSuit":"♠",
    "SpanFromAbove":"SpanFromAbove",
    "SpanFromBoth":"SpanFromBoth",
    "SpanFromLeft":"SpanFromLeft",
    "SphericalAngle":"∢",
    "Sqrt":"√",
    "Square":"Square",
    "SquareIntersection":"⊓",
    "SquareSubset":"⊏",
    "SquareSubsetEqual":"⊑",
    "SquareSuperset":"⊐",
    "SquareSupersetEqual":"⊒",
    "SquareUnion":"⊔",
    "Star":"⋆",
    "Sterling":"£",
    "Stigma":"ϛ",
    "Subset":"⊂",
    "SubsetEqual":"⊆",
    "Succeeds":"≻",
    "SucceedsEqual":"⪰",
    "SucceedsSlantEqual":"≽",
    "SucceedsTilde":"≿",
    "SuchThat":"∍",
    "Sum":"∑",
    "Superset":"⊃",
    "SupersetEqual":"⊇",
    "SystemEnterKey":"SystemEnterKey",
    "SystemsModelDelay":"SystemsModelDelay",
    "SZ":"ß",
    "TabKey":"TabKey",
    "Tau":"τ",
    "TaurusSign":"♉",
    "TensorProduct":"TensorProduct",
    "TensorWedge":"TensorWedge",
    "THacek":"ť",
    "Therefore":"∴",
    "Theta":"θ",
    "ThickSpace":"\u2004",
    "ThinSpace":"\u2009",
    "Thorn":"þ",
    "Tilde":"∼",
    "TildeEqual":"≃",
    "TildeFullEqual":"≅",
    "TildeTilde":"≈",
    "Times":"×",
    "Trademark":"™",
    "Transpose":"Transpose",
    "TripleDot":"TripleDot",
    "UAcute":"ú",
    "UDoubleAcute":"ű",
    "UDoubleDot":"ü",
    "UGrave":"ù",
    "UHat":"û",
    "UnderBrace":"︸",
    "UnderBracket":"⎵",
    "UnderParenthesis":"︶",
    "UndirectedEdge":"UndirectedEdge",
    "Union":"⋃",
    "UnionPlus":"⊎",
    "UpArrow":"↑",
    "UpArrowBar":"⤒",
    "UpArrowDownArrow":"⇅",
    "UpDownArrow":"↕",
    "UpEquilibrium":"⥮",
    "UpperLeftArrow":"↖",
    "UpperRightArrow":"↗",
    "UpPointer":"▴",
    "Upsilon":"υ",
    "UpTee":"⊥",
    "UpTeeArrow":"↥",
    "Uranus":"⛢",
    "URing":"ů",
    "Vee":"⋁",
    "Venus":"♀",
    "VerticalBar":"VerticalBar",
    "VerticalEllipsis":"⋮",
    "VerticalLine":"│",
    "VerticalSeparator":"VerticalSeparator",
    "VerticalTilde":"≀",
    "VeryThinSpace":"\u200a",
    "VirgoSign":"♍",
    "WarningSign":"WarningSign",
    "WatchIcon":"⌚",
    "Wedge":"⋀",
    "WeierstrassP":"℘",
    "WhiteBishop":"♗",
    "WhiteKing":"♔",
    "WhiteKnight":"♘",
    "WhitePawn":"♙",
    "WhiteQueen":"♕",
    "WhiteRook":"♖",
    "Wolf":"Wolf",
    "WolframLanguageLogo":"WolframLanguageLogo",
    "WolframLanguageLogoCircle":"WolframLanguageLogoCircle",
    "Xi":"ξ",
    "Xnor":"Xnor",
    "Xor":"⊻",
    "YAcute":"ý",
    "YDoubleDot":"ÿ",
    "Yen":"¥",
    "Zeta":"ζ",
    "ZHacek":"ž",
}

def named_characters_to_unicode(s: str) -> str:
    """
    Convert Mathematica's named characters to SymPy equivalents.

    The list of named characters is available at

        https://reference.wolfram.com/language/guide/ListingOfNamedCharacters.html
    """
    # Mathematica's named characters always start with `\[`, end with
    # `]`, and have only characters in [a-zA-Z] in between.
    if r"\[" in s:  # Don't bother if there's no `\[`
        pattern = r"\\\[([a-zA-Z]+)\]"
        def replace(match):
            name = match.group(1)
            return _named_characters.get(name, match.group(0))
        s = re.sub(pattern, replace, s)
    return s
