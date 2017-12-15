import os

from sympy.external import import_module
from sympy import (
    Abs,
    Add,
    asin,
    cos,
    cosh,
    cot,
    csc,
    Derivative,
    E,
    Eq,
    factorial,
    GreaterThan,
    Integral,
    LessThan,
    Limit,
    log,
    Mul,
    oo,
    Pow,
    Product,
    root,
    sec,
    sech,
    sin,
    sinh,
    sqrt,
    StrictGreaterThan,
    StrictLessThan,
    Sum,
    Symbol,
    tan,
    tanh,
)


_PARSER = None


def parse_latex(latex_str, debug=False):
    global _PARSER
    if _PARSER is None:
        _PARSER = LaTeXParser()
    return _PARSER.parse(latex_str, debug=debug)


def _Minus(p1, p3, **kw):
    return Add(p1, Mul(-1, p3, **kw), **kw)


def _Div(p1, p3, **kw):
    return Mul(p1, Pow(p3, -1, **kw), **kw)


class LaTeXParser(object):
    def __init__(self, **kw):
        lex = import_module('ply.lex', __import__kwargs={'fromlist': ['_']})
        yacc = import_module('ply.yacc', __import__kwargs={'fromlist': ['_']})

        if None in [lex, yacc]:
            return

        self.debug = kw.get('debug', 0)
        self.names = {}
        try:
            modname = os.path.split(os.path.splitext(__file__)[0])[
                1] + "_" + self.__class__.__name__
        except:
            modname = "parser" + "_" + self.__class__.__name__
        self.debugfile = modname + ".dbg"
        self.tabmodule = modname + "_" + "parsetab"
        lex.lex(module=self, debug=self.debug)
        yacc.yacc(module=self,
                  debug=self.debug,
                  debugfile=self.debugfile,
                  tabmodule=self.tabmodule)

        self._yacc = yacc

    def parse(self, latex_str, debug=False):
        return self._yacc.parse(latex_str, debug=debug)

    class ParseError(Exception):
        def __init__(self, msg, offset, p):
            self.msg = msg
            self.offset = offset
            self.p = p

        def __repr__(self):
            return "ParseError(%r, %r): %s" % (self.msg, self.offset, self.p)

        def __str__(self):
            return "%s at position %s: %s" % (self.msg, self.offset + 1, self.p)

    binary_ops = {
        '-': _Minus,
        '*': Mul,
        '/': _Div,
        '\\cdot': Mul,
        '\\div': _Div,
        '\\geq': GreaterThan,
        '\\leq': LessThan,
        '\\times': Mul,
        '^': Pow,
        '+': Add,
        '<': StrictLessThan,
        '<=': LessThan,
        '=': Eq,
        '>': StrictGreaterThan,
        '>=': GreaterThan,
    }

    functions = {
        '\\arcsin': asin,
        '\\cos': cos,
        '\\cosh': cosh,
        '\\cot': cot,
        '\\csc': csc,
        '\\ln': lambda x, **kw: log(x, E, **kw),
        '\\log': lambda x, **kw: log(x, 10, **kw),
        '\\sec': sec,
        '\\sech': sech,
        '\\sin': sin,
        '\\sinh': sinh,
        '\\sqrt': lambda x, **kw: sqrt(x),
        '\\tan': tan,
        '\\tanh': tanh,
    }

    constants = {
        '\\infty': oo,
    }

    tokens = (
        '_ARCSIN',
        '_CDOT',
        '_COS',
        '_COSH',
        '_COT',
        '_CSC',
        '_DIV',
        '_FRAC',
        '_GEQ',
        '_INFTY',
        '_INT',
        '_LEQ',
        '_LIM',
        '_LOG', '_LN',
        '_LONG_RIGHT_ARROW',
        '_RIGHT_ARROW',
        '_SEC',
        '_SIN',
        '_SINH',
        '_SQRT',
        '_TAN',
        '_TANH',
        '_TIMES',
        '_TO',
        'ADD',
        'BANG',
        'CARET',
        'CMD',
        'DIFFERENTIAL',
        'DIVIDE',
        'EQUAL',
        'FLOAT',
        'GREATER_EQUAL',
        'GREATER',
        'INTEGER',
        'LEFT_BRACE',
        'LEFT_BRACKET',
        'LEFT_PAREN',
        'LESS_EQUAL',
        'LESS',
        'LETTER',
        'MULTIPLY',
        'PIPE',
        'RIGHT_BRACE',
        'RIGHT_BRACKET',
        'RIGHT_PAREN',
        'SUBTRACT',
        'UNDERSCORE',
    )

    # Separating
    t_CARET = r'\^'
    t_LEFT_BRACE = r'\{'
    t_LEFT_BRACKET = r'\['
    t_LEFT_PAREN = r'\('
    t_PIPE = r'\|'
    t_RIGHT_BRACE = r'\}'
    t_RIGHT_BRACKET = r'\]'
    t_RIGHT_PAREN = r'\)'
    t_UNDERSCORE = r'_'

    # Operating
    t__CDOT = r'\\cdot'
    t__DIV = r'\\div'
    t__TIMES = r'\\times'
    t_ADD = r'\+'
    t_DIVIDE = r'/'
    t_EQUAL = r'='
    t_GREATER = r'>'
    t_GREATER_EQUAL = r'>='
    t_LESS = r'<'
    t_LESS_EQUAL = r'<='
    t_MULTIPLY = r'\*'
    t_SUBTRACT = r'\-'

    # Calling
    t__ARCSIN = r'\\arcsin'
    t__COS = r'\\cos'
    t__COSH = r'\\cosh'
    t__COT = r'\\cot'
    t__CSC = r'\\csc'
    t__INT = r'\\int'
    t__LN = r'\\ln'
    t__LOG = r'\\log'
    t__SEC = r'\\sec'
    t__SIN = r'\\sin'
    t__SINH = r'\\sinh'
    t__SQRT = r'\\sqrt'
    t__TAN = r'\\tan'
    t__TANH = r'\\tanh'

    # Fractioning
    t__FRAC = r'\\frac'
    t_DIFFERENTIAL = r'd'

    # Limiting
    t__LIM = r'\\lim'
    t__LONG_RIGHT_ARROW = r'\\[lL]ongrightarrow'
    t__RIGHT_ARROW = r'\\[rR]ightarrow'
    t__TO = r'\\to'

    # Ignoring
    t_ignore = ' \t'

    # Naming
    t_LETTER = r'[a-zA-Z]'

    # Unknown variable
    t_CMD = (r'\\(?!('
             + '|'.join([k[1:] for k in functions.keys()]) +
             r'|cdot|times|div'
             r'|int|sqrt'
             r'|frac'
             r'|lim|to|[rR]ightarrow|([lL]ong)rightarrow'
             r'|infty'
             r'))([a-zA-Z\d][a-zA-Z\d]*)')

    t__INFTY = r'\\infty'

    # Postfix
    t_BANG = r'\!'

    # Numbers
    def t_FLOAT(self, t):
        r'\d+[\.,]\d+'
        t.value = float(t.value)
        return t

    def t_INTEGER(self, t):
        r'\d+'
        t.value = int(t.value)
        return t

    # Leftovers
    def t_newline(self, t):
        r'\n+'
        t.lexer.lineno += t.value.count("\n")

    # Erroring
    def t_error(self, t):
        print("T_ERROR", t.lexer.lineno, t.lexpos, t)
        raise self.ParseError("unknown character", t.lexpos, t)

    # Precedence rules
    precedence = (
        ('left', 'DIFFERENTIAL'),
        ('left', '_LIM', '_FRAC', 'PIPE'),
        ('left', '_SIN', '_COS', ),
        ('left', 'ADD', 'SUBTRACT'),
        ('left', 'MULTIPLY', '_CDOT', '_TIMES', 'DIVIDE', '_DIV'),
        ('left', 'CARET', 'UNDERSCORE'),
        ('right', 'NEGATIVE'),
    )

    def p_statement_expr(self, p):
        '''statement : expression'''
        p[0] = p[1]

    def p_expression_integral(self, p):
        '''statement : _INT expression diff_expression
                     | _INT diff_expression'''
        print("INTEGRAL", p)
        if len(p) == 4:
            p[0] = Integral(p[3], p[2])
        elif len(p) == 3:
            p[0] = Integral(1, p[2])

    def p_expression_diff_expr(self, p):
        ''' diff_expression : DIFFERENTIAL expression'''
        p[0] = p[2]

    def p_expression_pipe(self, p):
        ''' expression : PIPE expression PIPE'''
        print("PIPE", p[2])
        p[0] = Abs(p[2])

    def p_expression_limit(self, p):
        ''' expression : _LIM UNDERSCORE LEFT_BRACE expression_limit_to RIGHT_BRACE expression'''
        print("LIMIT")
        var, approaching, direction = list(p[4])
        print("\tVAR", var)
        print("\tAPPROACHING", approaching)
        print("\tDIRECTION", direction)
        p[0] = Limit(p[6], var, approaching, dir=direction)

    def p_expression_limit_to(self, p):
        ''' expression_limit_to : name_expression _TO expression
                                | name_expression _RIGHT_ARROW expression
                                | name_expression _LONG_RIGHT_ARROW expression
                                | name_expression _TO expression CARET LEFT_BRACE direction RIGHT_BRACE
                                | name_expression _RIGHT_ARROW expression CARET LEFT_BRACE direction RIGHT_BRACE
                                | name_expression _LONG_RIGHT_ARROW expression CARET LEFT_BRACE direction RIGHT_BRACE '''
        p[0] = [p[1], p[3], p[6] if len(p) == 8 else '+']
        print("LIMIT_TO", list(p[0]))

    def p_expression_direction(self, p):
        ''' direction : ADD
                      | SUBTRACT
        '''
        print("DIRECTION", str(p[1]))
        p[0] = p[1]

    def p_expression_frac(self, p):
        ''' expression : _FRAC LEFT_BRACE DIFFERENTIAL RIGHT_BRACE LEFT_BRACE diff_expression RIGHT_BRACE expression
                       | _FRAC LEFT_BRACE diff_expression RIGHT_BRACE LEFT_BRACE diff_expression RIGHT_BRACE
                       | _FRAC LEFT_BRACE expression RIGHT_BRACE LEFT_BRACE expression RIGHT_BRACE '''
        print("FRAC", list(p))
        if len(p) == 8:
            p[0] = _Div(p[3], p[6], evaluate=False)
        elif len(p) == 10:
            p[0] = Derivative(p[9], p[7])
        elif len(p) == 11:
            p[0] = Derivative(p[4], 8)

    def p_expression_func_subp(self, p):
        ''' expression : _SIN CARET expression expression
                       | _SIN UNDERSCORE expression expression
                       | _LOG UNDERSCORE expression expression
                       | _LOG CARET expression expression
        '''
        print("FUNC_SUBP")
        if p[1] == '\\sin':
            print("\t", p[4])
            p[0] = asin(p[4], evaluate=False)
        if p[1] == '\\log':
            print("\t", p[1], p[4], p[3])
            p[0] = log(p[4], p[3], evaluate=False)

    def p_expression_func(self, p):
        '''expression : _SIN expression
                      | _TAN expression
                      | _COS expression
                      | _CSC expression
                      | _SEC expression
                      | _COT expression
                      | _COSH expression
                      | _SINH expression
                      | _TANH expression
                      | _ARCSIN expression
                      | _LOG expression
                      | _LN expression
                      | _SQRT expression'''
        p[0] = self.functions[p[1]](p[2], evaluate=False)

    def p_expression_binop(self, p):
        '''expression : expression ADD expression
                      | expression SUBTRACT expression
                      | expression MULTIPLY expression
                      | expression _CDOT expression
                      | expression _TIMES expression
                      | expression DIVIDE expression
                      | expression _DIV expression
                      | expression CARET expression
                      | expression EQUAL expression
                      | expression GREATER expression
                      | expression GREATER_EQUAL expression
                      | expression _GEQ expression
                      | expression LESS expression
                      | expression LESS_EQUAL expression
                      | expression _LEQ expression'''
        p[0] = self.binary_ops[p[2]](p[1], p[3], evaluate=False)

    def p_expression_bang(self, p):
        '''expression : expression BANG'''
        p[0] = factorial(p[1], evaluate=False)

    def p_expression_group_mul(self, p):
        '''expression : expression expression'''
        p[0] = Mul(p[1], p[2], evaluate=False)

    def p_expression_negative(self, p):
        '''expression : SUBTRACT expression %prec NEGATIVE'''
        p[0] = Mul(-1, p[2], evaluate=False)

    def p_expression_group(self, p):
        '''expression : LEFT_PAREN expression RIGHT_PAREN
                      | LEFT_BRACE expression RIGHT_BRACE
                      | LEFT_BRACKET expression RIGHT_BRACKET'''
        p[0] = p[2]

    def p_expression_number(self, p):
        '''expression_number : FLOAT
                             | INTEGER
                             | constant'''
        p[0] = p[1]

    def p_expression_constant(self, p):
        '''constant : _INFTY'''
        p[0] = self.constants[p[1]]

    def p_expression_number_expr(self, p):
        '''expression : expression_number'''
        p[0] = p[1]

    def p_expression_name(self, p):
        '''name_expression : LETTER
                           | CMD'''
        if p[1][0] == '\\':
            p[0] = Symbol(p[1][1:])
        else:
            p[0] = Symbol(p[1])

    def p_expression_name_expression(self, p):
        '''expression : name_expression'''
        p[0] = p[1]

    def p_error(self, p):
        raise self.ParseError(("Syntax error at '%s'" % p.value), 0, p)
