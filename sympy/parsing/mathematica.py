from __future__ import annotations
import re
import typing
from itertools import product
from typing import Any, Callable

import sympy
from sympy import Mul, Add, Pow, Rational, log, exp, sqrt, cos, sin, tan, asin, acos, acot, asec, acsc, sinh, cosh, tanh, asinh, \
    acosh, atanh, acoth, asech, acsch, expand, im, flatten, polylog, cancel, expand_trig, sign, simplify, \
    UnevaluatedExpr, S, atan, atan2, Mod, Max, Min, rf, Ei, Si, Ci, airyai, airyaiprime, airybi, primepi, prime, \
    isprime, cot, sec, csc, csch, sech, coth, Function, E, I, pi, Tuple, GreaterThan, StrictGreaterThan, StrictLessThan, \
    LessThan, Equality, Or, And, Lambda, Integer, Dummy, symbols, Not, factorial, \
    Derivative, Subs, Symbol
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

    Derivatives written with a prime are translated as such:

    >>> parse_mathematica("f'[x] + f''[x]")
    Derivative(f(x), x) + Derivative(f(x), (x, 2))
    >>> parse_mathematica("Sin'[x]")
    cos(x)

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


def parse_mathematica_to_fullformlist(s):
    """
    Translate a string containing a Wolfram Mathematica expression to its
    ``FullForm`` represented as nested Python lists.

    This is the intermediate representation produced by
    :func:`parse_mathematica` before it is turned into a SymPy expression:
    each node is a list whose first element is the head (as a string) followed
    by its arguments, and atoms are left as strings. It is useful when the SymPy
    conversion is undesired, for example to inspect the parsed syntax tree or to
    perform a custom conversion.

    Examples
    ========

    >>> from sympy.parsing.mathematica import parse_mathematica_to_fullformlist
    >>> parse_mathematica_to_fullformlist("Sin[x]^2 Tan[y]")
    ['Times', ['Power', ['Sin', 'x'], '2'], ['Tan', 'y']]
    >>> parse_mathematica_to_fullformlist("x*(a + b)")
    ['Times', 'x', ['Plus', 'a', 'b']]
    >>> parse_mathematica_to_fullformlist("F[7, 5, 3]")
    ['F', '7', '5', '3']
    """
    parser = MathematicaParser()
    return parser.parse_fullformlist(s)


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


def _literal_character_ranges(initial):
    """List all characters that may be used in identifiers

    Mathematica identifiers cannot contain underscores.  This function
    returns a regex range matching all characters that may be used in
    a Python identifier, excluding underscores.  The `initial` argument
    is a bool indicating whether this range describes valid initial
    characters (True) or valid subsequent characters (False).

    Mathematica additionally allows ``$`` anywhere in a symbol name
    (``$Version``) and ``` ` ``` as the context separator (``Global`sym``);
    neither is a Python identifier character, so both are added explicitly.
    """
    import sys
    def valid(character):
        if character == "_":
            return False
        elif character == "$":
            return True
        elif character == "`":
            return not initial
        elif initial:
            return character.isidentifier()
        else:
            return ("x" + character).isidentifier()

    # Note that 0 corresponds to the NUL character, which is not an
    # identifier, so we can use it as a signal that we haven't found a
    # valid character yet.
    characters = ""
    range_begin_character = ""
    range_begin_index = 0
    range_end_character = ""
    range_end_index = 0

    for codepoint in range(sys.maxunicode + 1):
        character = chr(codepoint)
        if valid(character):
            if range_begin_index == 0:
                range_begin_character = character
                range_begin_index = codepoint
            range_end_character = character
            range_end_index = codepoint
        elif range_begin_index != 0:
            # We have reached the end of a range (length 1, 2, or n),
            # so we add it to our list of characters.
            if range_begin_index == range_end_index:
                characters += range_begin_character
            elif range_begin_index + 1 == range_end_index:
                characters += range_begin_character + range_end_character
            else:
                characters += f"{range_begin_character}-{range_end_character}"
            # Reset the range
            range_begin_character = ""
            range_begin_index = 0

    # In case the loop above ended on a valid character we add it / close the range here
    if range_begin_index != 0:
        if range_begin_index == range_end_index:
            characters += range_begin_character
        elif range_begin_index + 1 == range_end_index:
            characters += range_begin_character + range_end_character
        else:
            characters += f"{range_begin_character}-{range_end_character}"

    return characters


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
        s4 = self.parse_fullformlist(s)
        s5 = self._from_fullformlist_to_sympy(s4)
        return s5

    def parse_fullformlist(self, s):
        s2 = named_characters_to_unicode(s)
        s3 = self._from_mathematica_to_tokens(s2)
        s4 = self._from_tokens_to_fullformlist(s3)
        return s4

    INFIX = "Infix"
    PREFIX = "Prefix"
    POSTFIX = "Postfix"
    FLAT = "Flat"
    RIGHT = "Right"
    LEFT = "Left"

    _mathematica_op_precedence: list[tuple[str, str | None, dict[str, str | Callable]]] = [
        (POSTFIX, None, {";": lambda x: x + ["Null"] if isinstance(x, list) and x and x[0] == "CompoundExpression" else ["CompoundExpression", x, "Null"]}),
        (INFIX, FLAT, {";": "CompoundExpression"}),
        # ``=`` and ``:=`` bind looser than the compound assignments, so
        # ``a = b += c`` is ``Set[a, AddTo[b, c]]`` but ``a += b = c`` is
        # ``Set[AddTo[a, b], c]``.
        (INFIX, RIGHT, {"=": "Set", ":=": "SetDelayed"}),
        (INFIX, RIGHT, {"+=": "AddTo", "-=": "SubtractFrom", "*=": "TimesBy", "/=": "DivideBy"}),
        (INFIX, RIGHT, {"\N{THEREFORE}": "Therefore"}),
        # ``x // y`` is postfix function application, equivalent to ``y[x]``
        # (e.g. ``expr // Simplify`` means ``Simplify[expr]``).
        (INFIX, LEFT, {"//": lambda x, y: [y, x]}),
        (POSTFIX, None, {"&": "Function"}),
        (INFIX, LEFT, {"/.": "ReplaceAll"}),
        # ``\[Rule]`` (U+F522) and ``\[RuleDelayed]`` (U+F51F) are the character
        # forms of ``->`` and ``:>``; accepting them here means both the
        # ``\[Rule]`` spelling and the character itself parse.
        (INFIX, RIGHT, {"->": "Rule", ":>": "RuleDelayed",
                        "\uf522": "Rule", "\uf51f": "RuleDelayed"}),
        (INFIX, FLAT, {"\N{LEFT RIGHT ARROW}": "LeftRightArrow"}),
        (INFIX, LEFT, {"/;": "Condition"}),
        # ``p : v`` names a pattern when ``p`` is a symbol (``s:{__}`` is
        # ``Pattern[s, List[BlankSequence[]]]``) and supplies a default value
        # otherwise (``x_:1`` is ``Optional[Pattern[x, Blank[]], 1]``).  Being
        # left associative, ``a:b:c`` is ``Optional[Pattern[a, b], c]``.
        (INFIX, LEFT, {":": lambda x, y: (["Pattern", x, y]
                                          if MathematicaParser._is_symbol(x)
                                          else ["Optional", x, y])}),
        (INFIX, FLAT, {"|": "Alternatives"}),
        (POSTFIX, None, {"..": "Repeated", "...": "RepeatedNull"}),
        (INFIX, FLAT, {"||": "Or"}),
        (INFIX, FLAT, {"&&": "And"}),
        (PREFIX, None, {"!": "Not"}),
        (INFIX, FLAT, {"===": "SameQ", "=!=": "UnsameQ"}),
        (INFIX, FLAT, {"==": "Equal", "!=": "Unequal", "<=": "LessEqual", "<": "Less", ">=": "GreaterEqual", ">": "Greater"}),
        (INFIX, FLAT, {"\N{ALMOST EQUAL TO}": "TildeTilde"}),
        (INFIX, FLAT, {";;": "Span"}),
        (INFIX, FLAT, {"+": "Plus", "-": "Plus"}),
        (INFIX, FLAT, {"\N{CIRCLED PLUS}": "CirclePlus"}),
        (INFIX, FLAT, {"\N{STAR OPERATOR}": "Star"}),
        (INFIX, FLAT, {"*": "Times", "/": "Times"}),
        (INFIX, FLAT, {"\N{CIRCLED TIMES}": "CircleTimes"}),
        (INFIX, FLAT, {".": "Dot"}),
        (PREFIX, None, {"-": lambda x: MathematicaParser._get_neg(x),
                        # Unary ``+`` is not a no-op: ``+x`` is ``Plus[x]``.
                        "+": lambda x: ["Plus", x]}),
        (PREFIX, None, {"\N{NABLA}": "Del", "\N{WHITE SQUARE}": "Square"}),
        (INFIX, RIGHT, {"^": "Power"}),
        # ``f'`` is ``Derivative[1][f]``; each extra prime raises the order, so
        # ``f''`` is ``Derivative[2][f]`` (not ``Derivative[Derivative[f]]``).
        # The prime binds tighter than ``^`` but looser than the ``@`` family and
        # than ``[``: ``a^b'`` is ``Power[a, Derivative[1][b]]`` while ``a@b'`` is
        # ``Derivative[1][a[b]]``.
        (POSTFIX, None, {"'": lambda x: MathematicaParser._get_derivative(x)}),
        (INFIX, RIGHT, {"@@": "Apply", "/@": "Map", "//@": "MapAll", "@@@": lambda x, y: ["Apply", x, y, ["List", "1"]]}),
        # ``f @ x`` is prefix function application, equivalent to ``f[x]``.  It
        # binds tighter than ``^``, ``'`` and ``@@`` but looser than ``[``.
        (INFIX, RIGHT, {"@": lambda x, y: [x, y]}),
        (POSTFIX, None, {"!": "Factorial", "!!": "Factorial2", "--": "Decrement"}),
        (INFIX, None, {"[": lambda x, y: [x, *y], "[[": lambda x, y: ["Part", x, *y]}),
        (PREFIX, None, {"{": lambda x: ["List", *x], "(": lambda x: x[0],
                        "\N{LEFT-POINTING ANGLE BRACKET}": lambda x: ["AngleBracket", *x]}),
        (INFIX, None, {"?": "PatternTest"}),
        (POSTFIX, None, {
            "_": lambda x: ["Pattern", x, ["Blank"]],
            "_.": lambda x: ["Optional", ["Pattern", x, ["Blank"]]],
            "__": lambda x: ["Pattern", x, ["BlankSequence"]],
            "___": lambda x: ["Pattern", x, ["BlankNullSequence"]],
        }),
        # ``_h``/``__h``/``___h`` restrict the blank to expressions with head
        # ``h``.  These are handled directly in ``_parse_after_braces``, which
        # also covers the anonymous forms the INFIX machinery cannot express
        # (``_h`` has no left operand at all); the entries are kept here to
        # declare the tokens and their precedence level.
        (INFIX, None, {"_": lambda x, y: ["Pattern", x, ["Blank", y]],
                       "__": lambda x, y: ["Pattern", x, ["BlankSequence", y]],
                       "___": lambda x, y: ["Pattern", x, ["BlankNullSequence", y]]}),
        (PREFIX, None, {"#": "Slot", "##": "SlotSequence"}),
        # ``s::name`` is ``MessageName[s, "name"]`` -- the name is a string, and
        # further ``::`` sections simply add more of them.  It binds tighter than
        # anything else, ``a::b[c]`` being ``MessageName[a, "b"][c]``.
        (INFIX, FLAT, {"::": lambda *a: ["MessageName", a[0],
                                         *[["_Str", n] if isinstance(n, str) else n
                                           for n in a[1:]]]}),
    ]

    # Contexts Mathematica has on ``$ContextPath`` in a default session, and so
    # strips from a symbol name.  Any other context is part of the name.
    _default_contexts = ("System`", "Global`")

    # Tokens that bind to a neighbour only when written without a space between
    # them; see ``token_split`` in ``_from_mathematica_to_tokens``.
    _whitespace_sensitive = ("_", "_.", "__", "___", "#", "##")

    # Heads built by each blank token, used for ``_h``, ``x_h`` and friends.
    _blank_heads = {"_": "Blank", "__": "BlankSequence", "___": "BlankNullSequence"}

    # How each bracketing prefix turns its group into a node.  Mirrors the PREFIX
    # entries above, for ``_complete_left_operand``.
    _enclosure_prefixes = {
        "(": lambda group: group[0],
        "{": lambda group: ["List", *group],
        "\N{LEFT-POINTING ANGLE BRACKET}": lambda group: ["AngleBracket", *group],
    }

    _missing_arguments_default = {
        "#": lambda: ["Slot", "1"],
        "##": lambda: ["SlotSequence", "1"],
        # A ``_`` construct only builds a Pattern when it follows a symbol.
        # Everywhere else it stands on its own as a bare blank -- see
        # ``_is_symbol``.
        "_": lambda: ["Blank"],
        "_.": lambda: ["Optional", ["Blank"]],
        "__": lambda: ["BlankSequence"],
        "___": lambda: ["BlankNullSequence"],
    }

    _blanks = ("_", "_.", "__", "___")

    # This regex matches any valid python identifier -- excluding
    # underscores, which Mathematica uses to denote patterns, and
    # therefore can't be part of a variable name.  The regex has the
    # form "[a][b]*", where `a` is the set of characters that can
    # start an identifier, and `b` is the set of characters that can
    # continue an identifier, which may also include numbers and
    # unicode combining characters.
    _literal = (
        "["
        + _literal_character_ranges(True)
        + "]["
        + _literal_character_ranges(False)
        + "]*"
    )

    _number = r"(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)"

    _enclosure_open = ["(", "[", "[[", "{", "\N{LEFT-POINTING ANGLE BRACKET}"]
    _enclosure_close = [")", "]", "]]", "}", "\N{RIGHT-POINTING ANGLE BRACKET}"]

    @classmethod
    def _get_neg(cls, x):
        return f"-{x}" if isinstance(x, str) and re.match(MathematicaParser._number, x) else ["Times", "-1", x]

    @classmethod
    def _get_inv(cls, x):
        return ["Power", x, "-1"]

    @classmethod
    def _get_derivative(cls, x):
        # ``x'`` is ``Derivative[1][x]``.  Stacking primes raises the order:
        # ``x''`` is ``Derivative[2][x]``, not ``Derivative[Derivative[x]]``.
        if (isinstance(x, list) and len(x) == 2 and isinstance(x[0], list)
                and len(x[0]) == 2 and x[0][0] == "Derivative"):
            return [["Derivative", str(int(x[0][1]) + 1)], x[1]]
        return [["Derivative", "1"], x]

    _infix_tighter_than_derivative: tuple | None = None

    @classmethod
    def _tighter_than_derivative(cls) -> tuple:
        """Infix operators that bind tighter than ``'``, brackets excluded.

        A prime takes as its operand everything to its left up to the first
        operator looser than itself, so in ``c @ f'[x]`` it captures the whole
        ``c @ f`` and the result is ``Derivative[1][c[f]][x]``.  Brackets are left
        out because those are completed rather than waited for.
        """
        if cls._infix_tighter_than_derivative is None:
            levels = cls._mathematica_op_precedence
            after = False
            found = []
            for op_type, _strat, op_dict in levels:
                if after and op_type == cls.INFIX:
                    found.extend(t for t in op_dict if t not in ("[", "[["))
                if "'" in op_dict:
                    after = True
            cls._infix_tighter_than_derivative = tuple(found)
        return cls._infix_tighter_than_derivative

    def _complete_left_operand(self, tokens: list, pointer: int) -> int:
        """Finish building the operand sitting just left of ``tokens[pointer]``.

        ``_parse_after_braces`` visits the precedence levels from the tightest to
        the loosest, but function application (``f[x]``, ``f[[i]]``) and ``'``
        (Derivative) are *positional*: they attach to whatever happens to precede
        them.  An operator that binds tighter than they do can therefore be
        reached while the operand on its left is only half assembled -- in
        ``f'[x]?c`` the ``?`` binds tighter than both, yet its left operand is the
        whole ``Derivative[1][f][x]``, which at that point is still the three
        loose tokens ``f``, ``'``, ``[``.

        Collapse those pending pieces, innermost first, and return how many tokens
        were removed so the caller can fix up its cursor.  Levels looser than
        ``'`` and ``[`` never see anything left to do here, so this is a no-op for
        them.
        """
        removed = 0
        while pointer > 0:
            if tokens[pointer - 1] == "'":
                start = pointer - 1
                while start > 0 and tokens[start - 1] == "'":
                    start -= 1
                # the primed expression may itself be a pending application
                inner = self._complete_left_operand(tokens, start)
                if inner:
                    removed += inner
                    pointer -= inner
                    continue
                if start == 0 or self._is_op(tokens[start - 1]):
                    return removed
                if start > 1 and tokens[start - 2] in self._tighter_than_derivative():
                    # The prime reaches further left than this single token; let
                    # that tighter operator be applied first.
                    return removed
                head = tokens[start - 1]
                for _ in range(pointer - start):
                    head = self._get_derivative(head)
                tokens[start - 1:pointer] = [head]
                removed += pointer - start
                pointer = start
            elif pointer > 1 and tokens[pointer - 2] in ("[", "[["):
                # the head of the application may need completing first
                inner = self._complete_left_operand(tokens, pointer - 2)
                if inner:
                    removed += inner
                    pointer -= inner
                    continue
                if pointer < 3 or self._is_op(tokens[pointer - 3]):
                    return removed
                head, args = tokens[pointer - 3], tokens[pointer - 1]
                if tokens[pointer - 2] == "[":
                    node = [head, *args]
                else:
                    node = ["Part", head, *args]
                tokens[pointer - 3:pointer] = [node]
                removed += 2
                pointer -= 2
            elif (pointer > 1 and isinstance(tokens[pointer - 2], str)
                    and tokens[pointer - 2] in self._enclosure_prefixes):
                # A bracketing prefix (``(``, ``{``, ``\[LeftAngleBracket]``) is
                # applied at a looser level than ``?`` or ``_``, so those would
                # otherwise take the raw group as their operand and ``(a+b)?c``
                # would lose the ``Plus``.
                build = self._enclosure_prefixes[tokens[pointer - 2]]
                tokens[pointer - 2:pointer] = [build(tokens[pointer - 1])]
                removed += 1
                pointer -= 1
            else:
                return removed
        return removed

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
                matches = list(tokenizer.finditer(code))
                if matches or code.isascii():
                    out: list = []
                    for i, match in enumerate(matches):
                        tok = match.group()
                        if i > 0 and matches[i - 1].end() != match.start():
                            # Whitespace is dropped from here on, but for blanks
                            # and slots it carries meaning: ``x_h`` is a pattern
                            # while ``x _ h`` is a product of three things, and
                            # ``#a`` is ``Slot["a"]`` while ``# a`` is
                            # ``Slot[1] a``.  Make the implied ``*`` explicit
                            # while the spacing is still visible.
                            prev = matches[i - 1].group()
                            if ((tok in self._whitespace_sensitive and not self._is_op(prev))
                                    or (prev in self._whitespace_sensitive and not self._is_op(tok))):
                                out.append("*")
                        out.append(tok)
                    return out
            return [code]

        token_lists = [token_split(code) for code in code_splits]
        tokens = [j for i in token_lists for j in i]
        # Remove newlines at the beginning
        while tokens and tokens[0] == "\n":
            tokens.pop(0)
        # Remove newlines at the end
        while tokens and tokens[-1] == "\n":
            tokens.pop(-1)

        # Resolve context-qualified symbols the way Mathematica does.  Only the
        # contexts on ``$ContextPath`` are dropped, which by default means
        # ``System``` and ``Global```: ``System`Plus`` is ``Plus``, but
        # ``Foo`Private`x`` keeps its context and stays ``Foo`Private`x``.
        for i, token in enumerate(tokens):
            if (isinstance(token, str) and token.startswith(self._default_contexts)
                    and self._is_symbol(token)):
                name = token.split("`", 1)[1]
                if name:
                    tokens[i] = name

        return tokens

    def _is_op(self, token: str | list) -> bool:
        if isinstance(token, list):
            return False
        if re.match(self._literal, token):
            return False
        if re.match("-?" + self._number, token):
            return False
        return True

    @classmethod
    def _is_symbol(cls, token: str | list) -> bool:
        # ``x_``, ``x__`` and ``x_h`` name the pattern they build only when ``x``
        # is a symbol.  Mathematica reads ``f[x]_`` as ``Times[f[x], Blank[]]``
        # and ``f[x]_c`` as ``Times[f[x], Blank[c]]``, not as patterns.  The same
        # rule decides between ``Pattern`` and ``Optional`` for ``p : v``.
        return isinstance(token, str) and re.fullmatch(cls._literal, token) is not None

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
                    if token == "+" and op_type == self.PREFIX:
                        # A lone unary plus is kept (``+x`` is ``Plus[x]``), but one
                        # heading a term of a ``+`` chain is absorbed: both
                        # ``a + +b`` and ``+a + b`` are plain ``Plus[a, b]``.  After
                        # an infix ``-`` it is kept, the term being negated:
                        # ``a - +b`` is ``Plus[a, Times[-1, Plus[b]]]``.
                        after_infix_plus = (pointer > 1 and tokens[pointer - 1] == "+"
                                            and not self._is_op(tokens[pointer - 2]))
                        heads_chain = (pointer == 0 and pointer + 2 < size
                                       and tokens[pointer + 2] in ("+", "-"))
                        if after_infix_plus or heads_chain:
                            tokens.pop(pointer)
                            size -= 1
                            changed = True
                            continue
                    # ``_h`` restricts a blank to head ``h``.  A symbol on its left
                    # names the resulting pattern (``x_h`` is
                    # ``Pattern[x, Blank[h]]``); with anything else on the left --
                    # another expression, an operator, or nothing at all -- the
                    # blank is anonymous and merely multiplies, so ``_h`` is
                    # ``Blank[h]`` and ``f[x]_h`` is ``Times[f[x], Blank[h]]``.
                    # This is handled here rather than through the INFIX machinery
                    # because that always consumes an operand on the left.
                    if op_type == self.INFIX and token in self._blank_heads:
                        if pointer == size - 1 or self._is_op(tokens[pointer + 1]):
                            pointer += 1
                            continue
                        blank = [self._blank_heads[token], tokens[pointer + 1]]
                        if pointer > 0 and self._is_symbol(tokens[pointer - 1]):
                            tokens[pointer - 1:pointer + 2] = [
                                ["Pattern", tokens[pointer - 1], blank]]
                            size -= 2
                        else:
                            tokens[pointer:pointer + 2] = [blank]
                            size -= 1
                            pointer += 1
                        changed = True
                        continue
                    if op_type in (self.INFIX, self.POSTFIX):
                        # Both take an operand on their left, which may still be
                        # a pending application or a run of ``'``.
                        absorbed = self._complete_left_operand(tokens, pointer)
                        if absorbed:
                            pointer -= absorbed
                            size -= absorbed
                            changed = True
                    if op_type == self.INFIX:
                        if pointer == 0 or pointer == size - 1 or self._is_op(tokens[pointer - 1]) or self._is_op(tokens[pointer + 1]):
                            pointer += 1
                            continue
                    # Special case: "!" without preceding operand is PREFIX Not, not POSTFIX Factorial
                    if token == "!" and op_type == self.POSTFIX:
                        if pointer == 0 or self._is_op(tokens[pointer - 1]):
                            pointer += 1
                            continue
                    # An operator whose operand is still an unapplied operator --
                    # the ``'`` in ``a&'``, the outer ``-`` in ``- -a`` -- has to
                    # wait for that one to be built.  Unless a default operand is
                    # defined for it (``#`` alone is ``Slot[1]``), leave it to the
                    # re-parsing loop at the end of this method.
                    if token not in self._missing_arguments_default:
                        if op_type == self.POSTFIX and (pointer == 0 or self._is_op(tokens[pointer - 1])):
                            pointer += 1
                            continue
                        if op_type == self.PREFIX and (pointer == size - 1 or self._is_op(tokens[pointer + 1])):
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
                            # Collect the whole right-associative chain of
                            # operands, then nest from the right.  The outermost
                            # application is left for the ``op_name(*node)`` call
                            # below (needed when ``op_name`` is a callable, e.g.
                            # ``@`` or ``@@@``); for a string head the node is
                            # built here directly.  The chain continues across any
                            # operator of this level, not just repetitions of the
                            # same one, since they share a precedence: ``a = b := c``
                            # is ``Set[a, SetDelayed[b, c]]``.  Each link keeps its
                            # own head, so they are collected alongside the operands.
                            operands = [arg1, arg2]
                            chain = [token]
                            while (pointer + 2 < size
                                   and isinstance(tokens[pointer + 1], str)
                                   and tokens[pointer + 1] in op_dict):
                                chain.append(tokens.pop(pointer + 1))
                                operands.append(tokens.pop(pointer + 1))
                                size -= 2
                            rest = operands[-1]
                            for i in range(len(chain) - 1, 0, -1):
                                link = op_dict[chain[i]]
                                rest = ([link, operands[i], rest] if isinstance(link, str)
                                        else link(operands[i], rest))
                            node.clear()
                            if isinstance(op_name, str):
                                node.extend([op_name, operands[0], rest])
                            else:
                                node.extend([operands[0], rest])
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
                        if token in ("#", "##"):
                            # Slot (``#``) and SlotSequence (``##``) have special
                            # argument rules that mirror the Wolfram Language:
                            #   #        -> Slot[1]
                            #   #3       -> Slot[3]           (positional: integer)
                            #   #name    -> Slot["name"]      (named: a string, i.e.
                            #                                  ["_Str", "name"] here)
                            #   ##       -> SlotSequence[1]
                            #   ##3      -> SlotSequence[3]
                            #   ##name   -> SlotSequence[1]*name  (SlotSequence only
                            #                                      takes integer indices,
                            #                                      so ``name`` is a
                            #                                      separate factor)
                            nxt = tokens[pointer + 1] if pointer + 1 < size else None
                            is_number = isinstance(nxt, str) and re.fullmatch(self._number, nxt) is not None
                            if nxt is None or self._is_op(nxt) or (token == "##" and not is_number):
                                tokens[pointer] = self._missing_arguments_default[token]()
                            elif is_number:
                                node.append(tokens.pop(pointer + 1))
                                size -= 1
                            else:  # ``#`` followed by a named slot
                                node.append(["_Str", tokens.pop(pointer + 1)])
                                size -= 1
                        elif pointer == size - 1 or self._is_op(tokens[pointer + 1]):
                            tokens[pointer] = self._missing_arguments_default[token]()
                        else:
                            node.append(tokens.pop(pointer+1))
                            size -= 1
                    elif op_type == self.POSTFIX:
                        if grouping_strat is not None:
                            raise TypeError("'Prefix' op_type should not have a grouping strat")
                        if (pointer == 0 or self._is_op(tokens[pointer - 1])
                                or (token in self._blanks
                                    and not self._is_symbol(tokens[pointer - 1]))):
                            # Stands on its own -- no operand is consumed, so the
                            # ``op_name`` call below must be skipped too.
                            tokens[pointer] = self._missing_arguments_default[token]()
                            pointer += 1
                            continue
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
        "Not": Not,
        "Function": _parse_Function,
        "Factorial": factorial,
    }

    _atom_conversions = {
        "I": I,
        "Pi": pi,
        "ExponentialE": E,
        "ImaginaryI": I,
        "ImaginaryJ": I,
    }

    def _derivative_to_sympy(self, expr: list, recurse: Callable):
        """Translate Mathematica's ``Derivative[n][f]`` idiom, or return ``None``.

        Mathematica curries the head of a derivative -- ``f'[x]`` is
        ``Derivative[1][f][x]``, the operator ``Derivative[1]`` applied to ``f``
        and the result applied to ``x``.  SymPy has no curried heads (an applied
        function is not callable again), so the whole idiom has to be recognised
        in one go rather than node by node.
        """
        applied = (isinstance(expr[0], list) and len(expr[0]) == 2
                   and isinstance(expr[0][0], list) and expr[0][0][0] == "Derivative")
        bare = isinstance(expr[0], list) and expr[0][0] == "Derivative" and len(expr) == 2
        if not (applied or bare):
            return None
        spec, fname = expr[0] if applied else (expr[0], expr[1])
        orders = spec[1:] if bare else expr[0][0][1:]
        if not all(isinstance(o, str) and o.isdigit() for o in orders):
            return None
        orders = [int(o) for o in orders]
        func = (self._node_conversions.get(fname, Function(fname))
                if isinstance(fname, str) else recurse(fname))
        if bare:
            # Nothing is applied to it, so the derivative operator itself is the
            # value: ``f'`` becomes ``Lambda(_x, Derivative(f(_x), _x))``.
            slots = [Dummy() for _ in orders]
            return Lambda(tuple(slots),
                          Derivative(func(*slots), *zip(slots, orders)).doit(deep=False))
        args = [recurse(a) for a in expr[1:]]
        if len(orders) != len(args):
            return None
        if len(set(args)) == len(args) and all(isinstance(a, Symbol) for a in args):
            # ``doit`` so that a known function is actually differentiated, as in
            # Mathematica: ``Sin'[x]`` is ``Cos[x]``.  It is a no-op on an
            # undefined function, and ``deep=False`` keeps the arguments intact.
            return Derivative(func(*args), *zip(args, orders)).doit(deep=False)
        # ``f'`` differentiates with respect to the argument *slot*, so anything
        # that is not a plain symbol has to be substituted in afterwards:
        # ``f'[2 x]`` is ``Subs(Derivative(f(_u), _u), _u, 2*x)``.
        slots = [Dummy() for _ in args]
        return Subs(Derivative(func(*slots), *zip(slots, orders)),
                    slots, args).doit(deep=False)

    def _from_fullformlist_to_sympy(self, full_form_list):

        def recurse(expr):
            if isinstance(expr, list):
                derivative = self._derivative_to_sympy(expr, recurse)
                if derivative is not None:
                    return derivative
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


def named_characters_to_unicode(s: str) -> str:
    """
    Convert Mathematica's named characters to SymPy equivalents.

    The list of named characters is available at

        https://reference.wolfram.com/language/guide/ListingOfNamedCharacters.html
    """
    from .mathematica_named_characters import mathematica_named_characters
    # Mathematica's named characters always start with `\[`, end with
    # `]`, and have only characters in [a-zA-Z] in between.
    if r"\[" in s:  # Don't bother if there's no `\[`
        pattern = r"\\\[([a-zA-Z]+)\]"
        def replace(match):
            name = match.group(1)
            if name not in mathematica_named_characters:
                raise ValueError(f"Unknown Mathematica named character: {name}")
            return mathematica_named_characters[name]
        s = re.sub(pattern, replace, s)
    if r"\[" in s:
        raise SyntaxError(f"Unmatched '\\[' in '{s}'")
    return s
