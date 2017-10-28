from __future__ import print_function, division

from itertools import product
import re
from sympy import sympify


def mathematica(s):
    return sympify(parse(s))


class parse(object):
    '''parse class converts basic Mathematica expression to Python.'''

    # translate
    F = {
        # (func in Mathematica, number of arguments): func in sympy
        ('Mod', 2): 'Mod',
        ('Sqrt', 1): 'sqrt',
        ('Exp', 1): 'exp',
        ('Log', 1): 'log',
    }

    # trigonometric, e.t.c.
    for arc, tri, h in product(('', 'Arc'), (
            'Sin', 'Cos', 'Tan', 'Cot', 'Sec', 'Csc'), ('', 'h')):

        fm = arc + tri + h

        if arc:  # arc func
            fs = 'a' + tri.lower() + h

        else:    # non-arc func
            fs = tri.lower() + h

        F.update({(fm, 1): fs})

    F_swap = {
        # when you say log(2,4) in Mathematica, base is 2.
        ('Log', 2): 'log',
    }

    replaces = {
        ' ': '',
        '^': '**',
        '[': '(',
        ']': ')',
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

        # Add omitted '*' character
        'add*_1': (
            re.compile(r'''
                (?<=[])\d])         # ], ) or a number
                                    # ''
                (?=[(a-zA-Z])       # ( or a single letter
                ''', re.VERBOSE),
            '*'),

        # Add omitted '*' character (variable letter preceding)
        'add*_2': (
            re.compile(r'''
                (?<=[a-zA-Z])       # a letter
                \(                  # ( as a character
                (?=.)               # any characters
                ''', re.VERBOSE),
            '*('),
    }

    func_anchor = re.compile(r'''
                (?:
                \A|(?<=[^a-zA-Z])   # at the top or a non-letter
                )
                [A-Z][a-zA-Z]*      # Function
                (?=\[)              # [ as a character
                ''', re.VERBOSE)

    def __new__(cls, s):
        return cls.eval(s)

    @classmethod
    def _convert_function(cls, s):
        '''Parse Mathematica function to sympy one'''

        pat = cls.func_anchor       # compiled regex object

        scanned = ''                 # converted string
        cur = 0                     # position cursor
        while True:
            m = pat.search(s)

            if m is None:
                # append the rest of string
                scanned += s
                break

            # get function name
            fm = m.group()

            # get arguments of function
            args, ancs = cls._get_args(m)

            # convert Mathematica function to sympy one
            s = cls._convert(pat, s, fm, args, ancs)

            # update cursor
            cur = m.end()

            # append converted part
            scanned += s[:cur]

            # shrink s
            s = s[cur:]

        return scanned

    @staticmethod
    def _get_args(m):
        '''Get arguments of Mathematica function'''

        s = m.string                # whole string
        anc = m.end() + 1           # pointing the first letter of arguments
        square, curly = [], []      # stack for brakets
        args, arg_ancs = [], []

        for i, c in enumerate(s[anc:], anc):
            # extract one argument
            if c == ',' and (not square) and (not curly):
                args.append(s[anc:i])       # add an argument
                arg_ancs.append((anc, i))   # add anchor position
                anc = i + 1                 # move anchor

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
                    args.append(s[anc:i])
                    arg_ancs.append((anc, i))
                    break

        return args, arg_ancs

    @classmethod
    def _convert(cls, pat, s, fm, args, ancs):
        # make key
        key = (fm, len(args))

        # convert function in F dictionary
        if key in cls.F:
            fs = cls.F[key]

        # convert function in F_swap dictionary
        elif key in cls.F_swap:
            fs = cls.F_swap[key]

            # swap argument.
            s = cls._swap_args(s, args, ancs)

        # just downcase function name
        else:
            fs = fm.lower()

        # replace Mathematica function with sympy function just once
        s = pat.sub(fs, s, count=1)

        return s

    @staticmethod
    def _swap_args(s, args, ancs):
        arg1, arg2 = args
        (bgn1, _), (_, end2) = ancs
        s = s[:bgn1] + arg2 + ',' + arg1 + s[end2:]
        return s

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
    def eval(cls, s):
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

        # '[', ']' to '(', ')', respectively
        s = cls._replace(s, '[')
        s = cls._replace(s, ']')

        return s
