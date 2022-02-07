import typing
from typing import Any, Dict as tDict, Tuple as tTuple, List, Optional, Union as tUnion, Callable

from itertools import product
import re

import sympy

from sympy import Mul, Add, Pow, log, exp, sqrt, cos, sin, tan, asin, acos, acot, asec, acsc, sinh, cosh, tanh, asinh, \
    acosh, atanh, acoth, asech, acsch, expand, im, flatten, polylog, cancel, expand_trig, sign, simplify, \
    UnevaluatedExpr, S
from sympy.core.sympify import sympify, _sympify


def mathematica(s, additional_translations=None):
    """
    Users can add their own translation dictionary.
    variable-length argument needs '*' character.

    Examples
    ========

    >>> from sympy.parsing.mathematica import mathematica
    >>> mathematica('Log3[9]', {'Log3[x]':'log(x,3)'})
    2
    >>> mathematica('F[7,5,3]', {'F[*x]':'Max(*x)*Min(*x)'})
    21

    """

    parser = MathematicaParser(additional_translations)
    return sympify(parser.parse(s))


def _deco(cls):
    cls._initialize_class()
    return cls


@_deco
class MathematicaParser:
    '''An instance of this class converts a string of a basic Mathematica
    expression to SymPy style. Output is string type.'''

    # left: Mathematica, right: SymPy
    CORRESPONDENCES = {
        'Sqrt[x]': 'sqrt(x)',
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
    TRANSLATIONS = {}  # type: tDict[tTuple[str, int], tDict[str, Any]]

    # cache for a raw users' translation dictionary
    cache_original = {}  # type: tDict[tTuple[str, int], tDict[str, Any]]

    # cache for a compiled users' translation dictionary
    cache_compiled = {}  # type: tDict[tTuple[str, int], tDict[str, Any]]

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
            d = {k: v for k, v in zip(x_args, args)}

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
        square, curly = [], []      # stack for brakets
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

    def parse(self, s):
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

    INFIX = "Infix"
    PREFIX = "Prefix"
    POSTFIX = "Postfix"
    FLAT = "Flat"
    RIGHT = "Right"
    LEFT = "Left"

    _mathematica_op_precedence: List[tTuple[str, Optional[str], tDict[str, tUnion[str, Callable]]]] = [
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

    _literal = r"[A-Za-z][A-Za-z0-9]*"
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
        tokenizer = re.compile("(" + "|".join(tokens) + ")")
        self._regex_tokenizer = tokenizer
        return self._regex_tokenizer

    def _tokenize_mathematica_code(self, code: str):
        tokenizer = self._get_tokenizer()

        # uncover '*' hiding behind a whitespace
        code = self._apply_rules(code, 'whitespace')

        # Remove white spaces:
        nol1 = r"(?:[A-Za-z0-9]|(?<=\d)\.)"
        nol2 = r"(?:[A-Za-z0-9]|\.(?=\d))"
        code = re.sub(rf"(?<={nol1})\s+(?={nol2})", "*", code)
        code = re.sub(rf"(?<=[])])\s+(?={nol2})", "*", code)
        code = re.sub(rf"(?<={nol1})\s+(?=[(])", "*", code)
        code = re.sub(r"(?<=\))\s+(?=\()", "*", code)

        # Handle patterns like "x.3" ==> "x*.3" and "x1.2" ==> "x1*.2"
        code = re.sub(r"(?<=[A-Za-z])(\d*)\.(?=\d)", r"\1*.", code)

        # add omitted '*' character
        code = self._apply_rules(code, 'add*_1')
        code = self._apply_rules(code, 'add*_2')

        tokens = tokenizer.findall(code)
        return tokens

    def _is_op(self, token: tUnion[str, list]) -> bool:
        if isinstance(token, list):
            return False
        if re.match(self._literal, token):
            return False
        if re.match("-?" + self._number, token):
            return False
        return True

    def _parse_tokenized_code(self, tokens: list):
        stack: List[list] = [[]]
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
                last_stack = self._parse_after_braces(stack[-1])
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
        assert len(stack) == 1
        return self._parse_after_braces(stack[0])

    def _parse_after_braces(self, tokens: list):
        op_dict: dict
        token_len_start: int = len(tokens)
        for op_type, flattening_strat, op_dict in reversed(self._mathematica_op_precedence):
            size: int = len(tokens)
            pointer: int = 0
            while pointer < size:
                token = tokens[pointer]
                if isinstance(token, str) and token in op_dict:
                    op_name: tUnion[str, Callable] = op_dict[token]
                    node: list
                    first_index: int
                    if isinstance(op_name, str):
                        node = [op_name]
                        first_index = 1
                    else:
                        node = []
                        first_index = 0
                    if op_type == self.PREFIX and pointer > 0 and not self._is_op(tokens[pointer - 1]):
                        pointer += 1
                        continue
                    if op_type == self.POSTFIX and pointer < size - 1 and not self._is_op(tokens[pointer + 1]):
                        pointer += 1
                        continue
                    if op_type == self.INFIX:
                        if pointer == 0 or pointer == size - 1 or self._is_op(tokens[pointer - 1]) or self._is_op(tokens[pointer + 1]):
                            pointer += 1
                            continue
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
                        if flattening_strat == self.FLAT:
                            while pointer + 1 < size and self._check_op_compatible(tokens[pointer+1], token):
                                node_p.append(arg2)
                                other_op = tokens.pop(pointer+1)
                                arg2 = tokens.pop(pointer+1)
                                if other_op == "/":
                                    arg2 = self._get_inv(arg2)
                                elif other_op == "-":
                                    arg2 = self._get_neg(arg2)
                                size -= 2
                            node_p.append(arg2)
                        elif flattening_strat == self.RIGHT:
                            while pointer + 1 < size and tokens[pointer+1] == token:
                                node_p.append([op_name, arg2])
                                node_p = node_p[-1]
                                tokens.pop(pointer+1)
                                arg2 = tokens.pop(pointer+1)
                                size -= 2
                            node_p.append(arg2)
                        elif flattening_strat == self.LEFT:
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
                        assert flattening_strat is None
                        if pointer == size - 1 or self._is_op(tokens[pointer + 1]):
                            tokens[pointer] = self._missing_arguments_default[token]()
                        else:
                            node.append(tokens.pop(pointer+1))
                            size -= 1
                    elif op_type == self.POSTFIX:
                        assert flattening_strat is None
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
        if len(tokens) != 1:
            if len(tokens) < token_len_start:
                # Trick to deal with cases in which an operator with lower
                # precedence should be transformed before an operator of higher
                # precedence. Such as in the case of `#&[x]` (that is
                # equivalent to `Lambda(d_, d_)(x)` in SymPy). In this case the
                # operator `&` has lower precedence than `[`, but needs to be
                # evaluated first because otherwise `# (&[x])` is not a valid
                # expression:
                return self._parse_after_braces(tokens)
            raise SyntaxError("unable to create a single AST for the expression")
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

    def _convert_fullform_to_pylist(self, wmexpr: str):
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

    def _convert_pylist_to_sympymform(self, pylist: list):
        from sympy import Function, Symbol

        def converter(expr):
            if isinstance(expr, list):
                if len(expr) > 0:
                    head = expr[0]
                    args = [converter(arg) for arg in expr[1:]]
                    return Function(head)(*args)
                else:
                    raise ValueError("error")
            elif isinstance(expr, str):
                return Symbol(expr)
            else:
                return _sympify(expr)

        return converter(pylist)

    _node_conversions = dict(
        Times=(Mul, False),
        Plus=(Add, False),
        Power=(Pow, False),
        Log=(log, False),
        Exp=(exp, False),
        Sqrt=(sqrt, False),
        Cos=(cos, False),
        Sin=(sin, False),
        Tan=(tan, False),
        Cot=(tan, True),
        Sec=(cos, True),
        Csc=(sin, True),
        ArcSin=(asin, False),
        ArcCos=(acos, False),
        # ArcTan=(atan, False),
        ArcCot=(acot, False),
        ArcSec=(asec, False),
        ArcCsc=(acsc, False),
        Sinh=(sinh, False),
        Cosh=(cosh, False),
        Tanh=(tanh, False),
        Coth=(tanh, True),
        Sech=(cosh, True),
        sech=(cosh, True),
        Csch=(sinh, True),
        csch=(sinh, True),
        ArcSinh=(asinh, False),
        ArcCosh=(acosh, False),
        ArcTanh=(atanh, False),
        ArcCoth=(acoth, False),
        ArcSech=(asech, False),
        ArcCsch=(acsch, False),
        Expand=(expand, False),
        Im=(im, False),
        Re=(sympy.re, False),
        Flatten=(flatten, False),
        Polylog=(polylog, False),
        Cancel=(cancel, False),
        # Gamma=(gamma, False),
        TrigExpand=(expand_trig, False),
        Sign=(sign, False),
        Simplify=(simplify, False),
        Defer=(UnevaluatedExpr, False),
        Identity=(S, False),
        # Sum=(Sum_doit, False),
        # Module=(With, False),
        # Block=(With, False),
        Null=(lambda *a: S.Zero, False),
    )

    def _convert_sympymform_to_sympy(self, mform):
        from sympy import Function

        expr = mform
        for mma_form, (sympy_node, invert) in self._node_conversions.items():
            if invert:
                expr = expr.replace(Function(mma_form), lambda *args: S.One/sympy_node(*args))
            else:
                expr = expr.replace(Function(mma_form), sympy_node)
        return expr
