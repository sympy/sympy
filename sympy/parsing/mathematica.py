from __future__ import print_function, division

from collections import defaultdict
from itertools import product
import re
from sympy import sympify


def mathematica(s):
    return sympify(parse(s))


def _deco(cls):
    cls._initialize_class()
    return cls


@_deco
class parse(object):
    '''parse class converts basic Mathematica expression to Python.'''

    # left: Mathematica, right: SymPy
    correspondence = {
        'Sqrt[x]': 'sqrt(x)',
        'Exp[x]': 'exp(x)',
        'Log[x]': 'log(x)',
        'Log[x,y]': 'log(y,x)',
        'Log2[x]': 'log(x,2)',
        'Log10[x]': 'log(x,10)',
        'Mod[x,y]': 'Mod(x,y)',
        'Max[*x]': 'Max(*x)',
        'Min[*x]': 'Min(*x)',
    }

    # trigonometric, e.t.c.
    for arc, tri, h in product(('', 'Arc'), (
            'Sin', 'Cos', 'Tan', 'Cot', 'Sec', 'Csc'), ('', 'h')):
        fm = arc + tri + h + '[x]'
        if arc:  # arc func
            fs = 'a' + tri.lower() + h + '(x)'
        else:    # non-arc func
            fs = tri.lower() + h + '(x)'
        correspondence.update({fm: fs})

    replaces = {
        ' ': '',
        '^': '**',
        '{': '[',
        '}': ']',
    }

    rules = {
        # a single whitespace to '*'
        'whitespace': (
            re.compile(r'''
                (?<=[a-zA-Z\d])     # a letter or a number
                \                   # a whitespace
                (?=[a-zA-Z\d])      # a letter or a number
                ''', re.VERBOSE),
            '*'),

        # add omitted '*' character
        'add*_1': (
            re.compile(r'''
                (?<=[])\d])         # ], ) or a number
                                    # ''
                (?=[(a-zA-Z])       # ( or a single letter
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

    # mathematica function name pattern
    fm_pattern = re.compile(r'''
                (?:
                \A|(?<=[^a-zA-Z])   # at the top or a non-letter
                )
                [A-Z][a-zA-Z\d]*    # Function
                (?=\[)              # [ as a character
                ''', re.VERBOSE)

    # list or matrix pattern (for future usage)
    arg_mtrx_pattern = re.compile(r'''
                \{.*\}
                ''', re.VERBOSE)

    # regex string for function argument pattern
    arg_pattern_template = r'''
                (?:
                \A|(?<=[^a-zA-Z])
                )
                {arguments}         # model argument like x, y,...
                (?=[^a-zA-Z])
                '''

    def __new__(cls, s):
        return cls.eval(s)

    @classmethod
    def _initialize_class(cls):
        '''Executed by a class decorator (_deco)'''

        cls.T = defaultdict(dict)
        cls._update_translation_dictionary(cls.correspondence)

    @classmethod
    def _update_translation_dictionary(cls, dic):
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

            # search function name
            m = cls.fm_pattern.search(fm)

            # if no-hit
            if m is None:
                err = "'{f}' function form is invalid.".format(f=fm)
                raise ValueError(err)

            # get mathematica function name like 'Log'
            fm_name = m.group()

            # get arguments of mathematica function
            args, end = cls._get_args(m)

            # funciton side check. (ex) '2*Func[x]' is invalid.
            if m.start() != 0 or end != len(fm):
                err = "'{f}' function form is invalid.".format(f=fm)
                raise ValueError(err)

            # check the last argument's 1st character
            if args[-1][0] == '*':
                key_arg = '*'
            else:
                key_arg = len(args)

            key = (fm_name, key_arg)

            # SymPy function template
            cls.T[key]['fs'] = fs

            # args are ['x', 'y'] for example
            cls.T[key]['args'] = args

            # convert '*x' to '\\*x' for regex
            re_args = [x if x[0] != '*' else '\\' + x for x in args]

            # for regex. Example: (?:(x|y|z))
            xyz = '(?:(' + '|'.join(re_args) + '))'

            # string for regex compile
            patStr = cls.arg_pattern_template.format(arguments=xyz)

            # compile pattern of arguments
            pat = re.compile(patStr, re.VERBOSE)

            cls.T[key]['pat'] = pat

    @classmethod
    def add_translation(cls, user_dic):
        '''Users can add their own translation dictionary
        # Example 1
        In [1]: user_dic
        Out[1]: {'Log3[x]': 'log(x, 3)'}

        In [2]: parse.add_translation(user_dic)

        In [3]: mathematica('Log3[81]')
        Out[3]: 4

        # Example 2
        In [1]: user_dic
        Out[1]: {'Func[a]': 'sqrt(a)*exp(1/a)**2'}

        In [2]: parse.add_translation(user_dic)

        In [3]: mathematica('Func[3]')
        Out[3]: sqrt(3)*exp(2/3)

        # Example 3
        In [1]: user_dic
        Out[1]: {'Func[x, y, *z]': 'Max(*z)*(x**2 + y**2)'}

        In [2]: parse.add_translation(user_dic)

        In [3]: mathematica('Func[2,3,1,2,3,4,5,6,7,8,9,10]')
        Out[3]: 130

        variable-length argument needs '*' character '''

        if not isinstance(user_dic, dict):
            raise ValueError('argument must be dict type')
        cls._update_translation_dictionary(user_dic)

    @classmethod
    def _convert_function(cls, s):
        '''Parse Mathematica function to sympy one'''

        # compiled regex object
        pat = cls.fm_pattern

        scanned = ''                # converted string
        cur = 0                     # position cursor
        while True:
            m = pat.search(s)

            if m is None:
                # append the rest of string
                scanned += s
                break

            # get mathematica function name
            fm = m.group()

            # get arguments, and the end position of fm function
            args, end = cls._get_args(m)

            # the start position of fm function
            bgn = m.start()

            # convert mathematica function to sympy one
            s = cls._convert_one_function(s, fm, args, bgn, end)

            # update cursor
            cur = bgn

            # append converted part
            scanned += s[:cur]

            # shrink s
            s = s[cur:]

        return scanned

    @classmethod
    def _convert_one_function(cls, s, fm, args, bgn, end):
        # no variable-length argument
        if (fm, len(args)) in cls.T:
            key = (fm, len(args))

            # x, y,... model arguments
            x_args = cls.T[key]['args']

            # make correspondence between model arguments and actual ones
            d = {k: v for k, v in zip(x_args, args)}

        # with variable-length argument
        elif (fm, '*') in cls.T:
            key = (fm, '*')

            # x, y,..*args (model arguments)
            x_args = cls.T[key]['args']

            # make correspondence between model arguments and actual ones
            d = {}
            for i, x in enumerate(x_args):
                if x[0] == '*':
                    d[x] = ','.join(args[i:])
                    break
                d[x] = args[i]

        # out of cls.T
        else:
            err = "'{f}' is out of the whitelist.".format(f=fm)
            raise ValueError(err)

        # template string of converted function
        template = cls.T[key]['fs']

        # regex pattern for x_args
        pat = cls.T[key]['pat']

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

            # add the correspoinding actual argument
            scanned += template[:xbgn] + d[x]

            # update cursor to the end of the model arugment
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

        # the end postion of function + 1 (the next position to ']' bracket)
        func_end = i + 1

        return args, func_end

    @classmethod
    def _replace(cls, s, bef):
        aft = cls.replaces[bef]
        s = s.replace(bef, aft)
        return s

    @classmethod
    def _apply_rules(cls, s, bef):
        pat, aft = cls.rules[bef]
        return pat.sub(aft, s)

    @classmethod
    def _check_input(cls, s):
        for bracket in (('[', ']'), ('{', '}'), ('(', ')')):
            if s.count(bracket[0]) != s.count(bracket[1]):
                err = "'{f}' function form is invalid.".format(f=s)
                raise ValueError(err)

        if '{' in s:
                err = "Currently list is not suported.".format(f=s)
                raise ValueError(err)

    @classmethod
    def eval(cls, s):
        # input check
        cls._check_input(s)

        # uncover '*' hiding behind a whitespace
        s = cls._apply_rules(s, 'whitespace')

        # remove whitespace(s)
        s = cls._replace(s, ' ')

        # add omitted '*' character
        s = cls._apply_rules(s, 'add*_1')
        s = cls._apply_rules(s, 'add*_2')

        # translate function
        s = cls._convert_function(s)

        # '^' to '**'
        s = cls._replace(s, '^')

        # 'Pi' to 'pi'
        s = cls._apply_rules(s, 'Pi')

        # '{', '}' to '[', ']', respectively
#        s = cls._replace(s, '{')   # currently list is not taken into account
#        s = cls._replace(s, '}')

        return s
