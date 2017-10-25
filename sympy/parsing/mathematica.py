from __future__ import print_function, division

import re
from sympy import sympify


def mathematica(s):
    return sympify(parse(s))


class _Parse(object):
    '''An instance of this class converts basic Mathematica expression to
    Python.'''

    replaces1 = (
        (' ', ''),
        ('^', '**'),
    )

    replaces2 = (
        ('[', '('),
        (']', ')'),
    )

    rules1 = (
        # 'Arc' to 'a'
        (r'''
                (               # group1
                .*              # any characters
                [^a-zA-Z]       # a single character except alphabet
                )
                Arc             # detect 'Arc'
                (               # group2
                [a-zA-Z]        # a single alphabet
                .*              # any characters
                )
                ''',
         lambda m: m.group(1) + 'a' + m.group(2)),

        # 'Arc' (at the top) to 'a'
        (r'''
                \AArc           # detect 'Arc' at the top of the string
                (               # group1
                [a-zA-Z]        # a single alphabet
                .*              # any characters
                )
                ''',
         lambda m: 'a' + m.group(1)),

        # Add omitted '*' character
        (r'''
                (               # group1
                .*              # any characters
                [])\d]          # ], ) or a single number
                )
                (               # group2
                [(a-zA-Z]       # ( or a single alphabet
                .*              # any characters
                )
                ''',
         lambda m: m.group(1) + '*' + m.group(2)),

        # Add omitted '*' character (variable letter preceding)
        (r'''
                (               # group1
                .*              # any characters
                [a-zA-Z]        # a single alphabet
                )
                \(              # ( as a character
                (.*)            # group2, any characters
                ''',
         lambda m: m.group(1) + '*(' + m.group(2)),
    )

    rules2 = (
        # 'mod' to 'Mod'
        (r'''
                (               # group1
                .*              # any characters
                [^a-zA-Z]       # a single character except alphabet
                )
                mod             # detect 'mod'
                (.*)            # group2, any characters
                ''',
         lambda m: m.group(1) + 'Mod' + m.group(2)),

        # 'mod' (at the top) to 'Mod'
        (r'''
                \Amod           # detect 'mod' at the top of the string
                (.*)            # group1, any characters
                ''',
         lambda m: 'Mod' + m.group(1)),
    )

    @staticmethod
    def _replace(s, replaces):
        for bef, aft in replaces:
            s = s.replace(bef, aft)
        return s

    @staticmethod
    def _apply_rules(s, rules):
        for rule, action in rules:
            pat = re.compile(rule, re.VERBOSE)  # VERBOSE: for readable code
            while True:
                m = re.search(pat, s)
                if m:
                    s = action(m)
                else:
                    break
        return s

    def __call__(self, s):
        # '^' to '**' and remove Whitespace(s)
        s = self._replace(s, self.replaces1)

        # 'Arc' to 'a' and add omitted '*' character
        s = self._apply_rules(s, self.rules1)

        # convert to lower letters
        s = s.lower()

        # 'mod' to 'Mod'
        s = self._apply_rules(s, self.rules2)

        # '[', ']' to '(', ')', respectively
        s = self._replace(s, self.replaces2)

        return s


# instantiate _Parse
parse = _Parse()
