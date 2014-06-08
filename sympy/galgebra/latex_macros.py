#!/usr/bin/python

import pyparsing # make sure you have this installed
latex_chrs = (pyparsing.Word(pyparsing.alphanums) | '+' | '-' | "\\" | '*')
braces = pyparsing.nestedExpr( '{', '}', content=latex_chrs)

def parse_braces(s):
    res =  braces.parseString(s)
    print res.asList()
    return


macros = '\\newcommand{\\f}[2]{{#1}\\left ( {#2} \\right )}\n' \
         '\\newcommand{\\bfrac{}{}}[2]{{\displaystyle \\frac{#1}{#2}}}'
latex = '\\f{\\sin}{\\theta}\n' \
        '\\bfrac{1}{2}'

parse_braces(latex)



