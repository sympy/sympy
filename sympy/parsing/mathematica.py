from __future__ import print_function, division

import re
from sympy import sympify


def mathematica(s):
    return sympify(parse(s))


def parse(s):
    s = s.strip()

    # Begin rules
    rules = (
        # Arithmetic operation between a constant and a function
        (r"""
            \A          # Matches only at the start of the string
            (           # group1 
            [+-]?       # + or - sign with 0 or 1 repetition 
            \d+         # any number with 1 or more repetition 
            [*/+^-]     # either of arithmetic operands 
            )
            (           # group1
            \w+         # alphanumeric character and the underscore
            \[          # [ as a character
            [^\]]+      # character except ] with 1 or more repetittion
            [^[]*       # character except [ with 0 or more repetittion
            \]          # ] as a character
            ) 
            \Z          # Matches only at the end of the string
            """,
        lambda m: m.group(1) + parse(m.group(2))),

        # Arithmetic operation between two functions
        (r"""
            \A
            (           # group1
            \w+         
            \[          
            [^\]]+      
            [^\[]*      
            \]          
            ) 
            (           # group2
            [*/+-^]     
            ) 
            (           # group3
            \w+
            \[
            [^\]]+
            [^\[]*
            \]
            )
            \Z
            """,
        lambda m: parse(m.group(1)) + m.group(2) + parse(m.group(3))),

        # Function call
        (r"""
            \A
            ([+-]?)     # group1, + or - sign with 0 or 1 repetition 
            (\w+)       # group2 
            \[          
            (           # group3
            [^\]]+      
            [^[]*       
            ) 
            \]          
            \Z
            """,  
        lambda m: m.group(1) + translateFunction(
            m.group(2)) + "(" + parse(m.group(3)) + ")"),

        # Parenthesized implied multiplication
        (r"""
            \(          # ( as a character
            (.+)        # group1, any character with 1 or more repetition
            \)          # ) as a character
            \(
            (.+)        # group2
            \)
            """,  
        lambda m: "(" + parse(m.group(1)) + ")*(" + parse(m.group(2)) + ")"),

        # Parenthesized expression
        (r"""
            \A
            ([+-]?)     # group1, + or - sign with 0 or 1 repetition
            \(      
            (.+)        # group2
            \)      
            \Z
            """,  
        lambda m: m.group(1) + "(" + parse(m.group(2)) + ")"),

        # Implied multiplication - a(b)
        (r"""
            \A
            (           # group1
            .*          
            [\w.]       # either of alphanumeric character, the underscore or dot
            )
            \(          
            (.+)        # group2
            \)      
            \Z
            """, 
        lambda m: parse(m.group(1)) + "*(" + parse(m.group(2)) + ")"),

        # Implied multiplication - (a)b
        (r"""
            \A
            \(      
            (.+)        # group1
            \)      
            (           # group2
            [\w.]
            .*
            )
            \Z
            """,  
        lambda m: "(" + parse(m.group(1)) + ")*" + parse(m.group(2))),

        # Implied multiplication - 2a
        (r"""
            \A
            (           # group1
            -?          # - sign with 0 or 1 repetition
            \ *         # space with 0 or more repetition
            [\d.]+     # either of number or dot with 1 or more repetition
            )
            (           # group2
            [a-zA-Z]    
            .*
            )
            \Z
            """,  
        lambda m: parse(m.group(1)) + "*" + parse(m.group(2))),

        # Infix operator
        (r"""
            \A
            (           # group1
            [^=]+       # any characters except = with 1 or more repetition
            )
            (           # group2
            [*/+^=-]    # operands and = 
            =?          # = with 0 or 1 repetition
            )
            (.+)        # group3
            \Z
            """,  
        lambda m: parse(m.group(1)) + translateOperator(
            m.group(2)) + parse(m.group(3))),
        )
    # End rules

    for rule, action in rules:
        rule = re.compile(rule, re.VERBOSE)
        m = re.match(rule, s)
        if m:
            return action(m)

    return s


def translateFunction(s):
    if s.startswith("Arc"):
        return "a" + s[3:]
    return s.lower()


def translateOperator(s):
    dictionary = {'^': '**'}
    if s in dictionary:
        return dictionary[s]
    return s
