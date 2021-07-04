from typing import Any, Dict, Tuple
from itertools import product
import re
from sympy import sympify


def maple(s, additional_translations=None):
    parser = MapleParser(additional_translations)
    return sympify(parser.parse(s))

def _deco(cls):
    cls._initialize_class()
    return cls

@_deco
class MapleParser:
    '''
    An instance of this class converts a string of a basic Maple
    expression to SymPy style. Output is string type.

    Examples
    ========

    >>> from sympy.parsing.maple import maple
    >>> maple("diff(y(x) ,x) + x")
    x + Derivative(y(x), x)

    >>> maple("x**2*diff(y(x), x, x) + sqrt(x)")
    sqrt(x) + x**2*Derivative(y(x), (x, 2))

    >>> maple("(y(x)*sqrt(y(x)^2+x^2)+(y(x)^2-x^2) \
        *sin(alpha)-2*x*y(x)*cos(alpha))*diff(y(x), \
        x)+x*sqrt(y(x)^2+x^2)+2*x*y(x)*sin(alpha)+ \
        (y(x)^2-x^2)*cos(alpha)")
    x*sqrt(x**2 + y(x)**2) + 2*x*y(x)*sin(alpha) + (-x**2 + y(x)**2)*cos(alpha) + (-2*x*y(x)*cos(alpha) + (-x**2 + y(x)**2)*sin(alpha) + sqrt(x**2 + y(x)**2)*y(x))*Derivative(y(x), x)
    '''

    # left: Maple, right: SymPy
    CORRESPONDENCES = {
        'log(x,y)': 'log(y,x)',
        'ln(x)': 'log(x)',
        'log10(x)': 'log(x,10)',
        'mod(x,y)': 'Mod(x,y)',
        'max(*x)': 'Max(*x)',
        'min(*x)': 'Min(*x)',
        'dsolve(*x)': 'dsolve(*x)',
        'Int(y, x)': 'Integral(y, x)',
        'pochhammer(x,y)':'rf(x,y)',
        'arctan(x,y)':'atan2(y,x)',
        'BesselJ(x,y)':'besselj(x,y)',
        'BesselY(x,y)':'bessely(x,y)'
    }

    functions = {k.split('(')[0]: len(v)-len(k) for k, v in CORRESPONDENCES.items()}

    # trigonometric, e.t.c.
    for arc, tri, h in product(('', 'arc'), (
            'sin', 'cos', 'tan', 'cot', 'sec', 'csc'), ('', 'h')):
        fm = arc + tri + h + '(x)'
        if arc:  # arc func
            fs = 'a' + tri.lower() + h + '(x)'
        else:    # non-arc func
            fs = tri.lower() + h + '(x)'
        CORRESPONDENCES.update({fm: fs})

    REPLACEMENTS = {
        ' ': '',
        '^': '**',
    }

    RULES = {
        # a single whitespace to '*'
        'whitespace': (
            re.compile(r'''
                (?<=[a-zA-Z\d])     # a letter or a number
                \                   # a whitespace
                (?=[a-zA-Z\d])      # a letter or a number
                ''', re.VERBOSE),
            '*'),
    }

    # Maple function name pattern
    FM_PATTERN = re.compile(r'''
                (?:
                \A|(?<=[^a-zA-Z])   # at the top or a non-letter
                )
                [a-zA-Z\d]*         # Function
                (?=\()              # ( as a character
                ''', re.VERBOSE)

    # list or matrix pattern (for future usage)
    ARG_MTRX_PATTERN = re.compile(r'''
                \[.*\]
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
    TRANSLATIONS = {}  # type: Dict[Tuple[str, int], Dict[str, Any]]

    # cache for a raw users' translation dictionary
    cache_original = {}  # type: Dict[Tuple[str, int], Dict[str, Any]]

    # cache for a compiled users' translation dictionary
    cache_compiled = {}  # type: Dict[Tuple[str, int], Dict[str, Any]]

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

            # search Maple function name
            m = cls.FM_PATTERN.search(fm)

            # if no-hit
            if m is None:
                err = "'{f}' function form is invalid.".format(f=fm)
                raise ValueError(err)

            # get Maple function name like 'Log'
            fm_name = m.group()

            # get arguments of Maple function
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
        '''Parse Maple function to SymPy one'''

        for fn in self.functions:
            pat = re.compile(fn, re.VERBOSE)
            matches = pat.finditer(s)

            for i, m in enumerate(matches):
                if m is None:
                    continue

                # get arguments, and the end position of fm function
                args, end = self._get_args(m)

                # the start position of fm function
                bgn = m.start()
                bgn += i*self.functions[fn]
                end += i*self.functions[fn]
                # convert Maple function to SymPy one
                s = self._convert_one_function(s, fn, args, bgn, end)
        return s

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
        '''Get arguments of a Maple function'''

        s = m.string                # whole string
        anc = m.end() + 1           # pointing the first letter of arguments
        square, curly = [], []      # stack for brakets
        args = []

        # current cursor
        cur = anc
        for i, c in enumerate(s[anc:], anc):
            # extract one argument
            if c == ',' and (not square) and (not curly):
                args.append(s[cur:i].strip())       # add an argument
                cur = i + 1                 # move cursor

            # handle list or matrix (for future usage)
            if c == '[':
                curly.append(c)
            elif c == ']':
                curly.pop()

            # seek corresponding ')' with skipping irrevant ones
            if c == '(':
                square.append(c)
            elif c == ')':
                if square:
                    square.pop()
                else:   # empty stack
                    args.append(s[cur:i].strip())
                    break

        # the next position to ')' bracket (the function end)
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

        # translate function
        s = self._convert_function(s)

        # '^' to '**'
        s = self._replace(s, '^')

        return s
