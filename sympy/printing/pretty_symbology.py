from unicodedata import lookup as U
import re


# prefix conventions when constructing tables
# L   - LATIN     i
# G   - GREEK     beta
# D   - DIGIT     0
# S   - SYMBOL    +


__all__ = ['greek','sub','sup','vobj','hobj','pretty_symbol']


_use_unicode = False

def pretty_use_unicode(flag = None):
    """Set whether pretty-printer should use unicode by default"""
    global _use_unicode
    if flag is None:
        return _use_unicode

    use_unicode_prev = _use_unicode
    _use_unicode = flag
    return use_unicode_prev


def xstr(*args):
    """call str or unicode depending on current mode"""
    if _use_unicode:
        return unicode(*args)
    else:
        return str(*args)

# GREEK
g   = lambda l: U('GREEK SMALL LETTER %s' % l.upper())
G   = lambda l: U('GREEK CAPITAL LETTER %s' % l.upper())

# XXX lambda <-> lamda
greek_letters = [
    'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta', 'theta',
    'iota', 'kappa', 'lamda', 'mu', 'nu', 'xi', 'omicron', 'pi', 'rho',
    'sigma', 'tau', 'upsilon', 'phi', 'chi', 'psi', 'omega' ]

# {}  greek letter -> (g,G)
greek = dict([(l, (g(l), G(l))) for l in greek_letters])


digit_2txt = {
    '0' :   'ZERO',
    '1' :   'ONE',
    '2' :   'TWO',
    '3' :   'THREE',
    '4' :   'FOUR',
    '5' :   'FIVE',
    '6' :   'SIX',
    '7' :   'SEVEN',
    '8' :   'EIGHT',
    '9' :   'NINE',
}

symb_2txt = {
    '+' :   'PLUS SIGN',
    '-' :   'MINUS',
    '=' :   'EQUALS SIGN',
    '(' :   'LEFT PARENTHESIS',
    ')' :   'RIGHT PARENTHESIS',
    '[' :   'LEFT SQUARE BRACKET',
    ']' :   'RIGHT SQUARE BRACKET',
    '{' :   'LEFT CURLY BRACKET',
    '}' :   'RIGHT CURLY BRACKET',

    # non-std
    'sum':  'SUMMATION',
    'int':  'INTEGRAL',
}

# SUBSCRIPT & SUPERSCRIPT
LSUB = lambda letter: U('LATIN SUBSCRIPT SMALL LETTER %s' % letter.upper())
GSUB = lambda letter: U('GREEK SUBSCRIPT SMALL LETTER %s' % letter.upper())
DSUB = lambda digit:  U('SUBSCRIPT %s' % digit_2txt[digit])
SSUB = lambda symb:   U('SUBSCRIPT %s' % symb_2txt[symb])

LSUP = lambda letter: U('SUPERSCRIPT LATIN SMALL LETTER %s' % letter.upper())
DSUP = lambda digit:  U('SUPERSCRIPT %s' % digit_2txt[digit])
SSUP = lambda symb:   U('SUPERSCRIPT %s' % symb_2txt[symb])

sub = {}    # symb -> subscript symbol
sup = {}    # symb -> superscript symbol

# latin subscripts
for l in 'aeioruvx':
    sub[l] = LSUB(l)

for l in 'in':
    sup[l] = LSUP(l)

for g in ['beta', 'gamma', 'rho', 'phi', 'chi']:
    sub[g] = GSUB(g)

for d in [str(i) for i in range(10)]:
    sub[d] = DSUB(d)
    sup[d] = DSUP(d)

for s in '+-=()':
    sub[s] = SSUB(s)
    sup[s] = SSUP(s)


# VERTICAL OBJECTS
HUP = lambda symb: U('%s UPPER HOOK'    % symb_2txt[symb])
CUP = lambda symb: U('%s UPPER CORNER'  % symb_2txt[symb])
MID = lambda symb: U('%s MIDDLE PIECE'  % symb_2txt[symb])
EXT = lambda symb: U('%s EXTENSION'     % symb_2txt[symb])
HLO = lambda symb: U('%s LOWER HOOK'    % symb_2txt[symb])
CLO = lambda symb: U('%s LOWER CORNER'  % symb_2txt[symb])
TOP = lambda symb: U('%s TOP'           % symb_2txt[symb])
BOT = lambda symb: U('%s BOTTOM'        % symb_2txt[symb])

# {} '('  ->  (extension, start, end, middle) 1-character
_xobj_unicode = {

    # vertical symbols
    #          ext       top       bot        mid           c1
    '(' :   (( EXT('('), HUP('('), HLO('(') ),              '('),
    ')' :   (( EXT(')'), HUP(')'), HLO(')') ),              ')'),
    '[' :   (( EXT('['), CUP('['), CLO('[') ),              '['),
    ']' :   (( EXT(']'), CUP(']'), CLO(']') ),              ']'),
    '{' :   (( EXT('('), HUP('{'), HLO('{'),  MID('{')  ),  '{'),   # XXX EXT is wrong
    '}' :   (( EXT(')'), HUP('}'), HLO('}'),  MID('}')  ),  '}'),   # XXX EXT is wrong
    '|' :   U('BOX DRAWINGS LIGHT VERTICAL'),

    'int':  (( EXT('int'), U('TOP HALF INTEGRAL'), U('BOTTOM HALF INTEGRAL') ), U('INTEGRAL')),
   #'sum':  ( U('N-ARY SUMMATION'), TOP('sum'), None, None, BOT('sum')     ),


    # horizontal objects
    #'-' :  '-',
    '-' :   U('BOX DRAWINGS LIGHT HORIZONTAL'),
    '_' :   U('HORIZONTAL SCAN LINE-9'),        # XXX symbol ok?

    # diagonal objects '\' & '/' ?
    '/' :   U('BOX DRAWINGS LIGHT DIAGONAL UPPER RIGHT TO LOWER LEFT'),
    '\\':   U('BOX DRAWINGS LIGHT DIAGONAL UPPER LEFT TO LOWER RIGHT'),
}

_xobj_ascii = {
    # vertical symbols
    #          ext  top   bot   mid         c1
    '(' :   (( '|', '/',  '\\'  ),          '('),
    ')' :   (( '|', '\\', '/'   ),          ')'),
    '[' :   (( '|', '-',  '-'   ),          '['),
    ']' :   (( '|', '-',  '-'   ),          ']'),
    '{' :   (( '|', '/',  '\\', '<' ),      '{'),
    '}' :   (( '|', '\\', '/',  '>' ),      '}'),
    '|' :   '|',

    'int':  ( ' | ', '  /', '/  ' ),

    # horizontal objects
    '-' :   '-',
    '_' :   '_',

    # diagonal objects '\' & '/' ?
    '/' :   '/',
    '\\':   '\\',
}


def xobj(symb, length):
    """Construct spatial object of given length.

    return: [] of equal-length strings
    """

    assert length > 0

    if _use_unicode:
        _xobj = _xobj_unicode
    else:
        _xobj = _xobj_ascii

    vinfo = _xobj[symb]

    c1 = top = bot = mid = None

    if not isinstance(vinfo, tuple):        # 1 entry
        ext = vinfo
    else:
        if isinstance(vinfo[0], tuple):     # (vlong), c1
            vlong = vinfo[0]
            c1    = vinfo[1]
        else:                               # (vlong), c1
            vlong = vinfo

        ext = vlong[0]

        try:
            top = vlong[1]
            bot = vlong[2]
            mid = vlong[3]
        except IndexError:
            pass

    if c1  is None:  c1  = ext
    if top is None:  top = ext
    if bot is None:  bot = ext
    if mid is not None:
        if (length % 2) == 0:
            raise ValueError('xobj: expect length = 2*k+1')
    else:
        mid = ext

    if length == 1:
        return c1


    res = []
    next= (length-2)//2
    nmid= (length-2) - next*2

    res += [top]
    res += [ext]*next
    res += [mid]*nmid
    res += [ext]*next
    res += [bot]

    return res


def vobj(symb, height):
    return '\n'.join( xobj(symb, height) )

def hobj(symb, width):
    return ''.join( xobj(symb, width) )

# RADICAL
# n -> symbol
root = {
    2   :   U('SQUARE ROOT'),   # U('RADICAL SYMBOL BOTTOM')
    3   :   U('CUBE ROOT'),
    4   :   U('FOURTH ROOT'),
}


# RATIONAL
VF  = lambda txt:   U('VULGAR FRACTION %s' % txt)

# (p,q) -> symbol
frac = {
    (1,2)   :   VF('ONE HALF'),
    (1,3)   :   VF('ONE THIRD'),
    (2,3)   :   VF('TWO THIRDS'),
    (1,4)   :   VF('ONE QUARTER'),
    (3,4)   :   VF('THREE QUARTERS'),
    (1,5)   :   VF('ONE FIFTH'),
    (2,5)   :   VF('TWO FIFTHS'),
    (3,5)   :   VF('THREE FIFTHS'),
    (4,5)   :   VF('FOUR FIFTHS'),
    (1,6)   :   VF('ONE SIXTH'),
    (5,6)   :   VF('FIVE SIXTHS'),
    (1,8)   :   VF('ONE EIGHTH'),
    (3,8)   :   VF('THREE EIGHTHS'),
    (5,8)   :   VF('FIVE EIGHTHS'),
    (7,8)   :   VF('SEVEN EIGHTHS'),
}


# RELATIONAL
relations = {
    '=='    : ( '=',    '='),
    '<'     : ( '<',    '<'),
    '<='    : ('<=',    U('LESS-THAN OR EQUAL TO')),
    '>='    : ('>=',    U('GREATER-THAN OR EQUAL TO')),
    '!='    : ('!=',    U('NOT EQUAL TO')),
}


def xrel(rel_name):
    op = relations[rel_name]

    if _use_unicode:
        return op[1]
    else:
        return op[0]


# SYMBOLS

atoms_table = {
    # class         how-to-display
    'Exp1'              :   U('SCRIPT SMALL E'),
    'Pi'                :   U('GREEK SMALL LETTER PI'),
    'Infinity'          :   U('INFINITY'),
    'NegativeInfinity'  :   '-'+U('INFINITY'),
    'ImaginaryUnit'     :   U('GREEK SMALL LETTER IOTA'),
    #'ImaginaryUnit'     :   U('MATHEMATICAL ITALIC SMALL I'),
    #'ImaginaryUnit'     :   U('DOUBLE-STRUCK ITALIC SMALL I'),
}

def pretty_atom(atom_name, default=None):
    """return pretty representation of an atom"""
    if _use_unicode:
        return atoms_table[atom_name]
    else:
        if default is not None:
            return default

        raise KeyError('only unicode')  # send it default printer

def pretty_symbol(symb_name):
    """return pretty representation of a symbol"""
    # let's split symb_name into symbol + index
    # UC: beta1
    # UC: f_beta

    if not _use_unicode:
        return symb_name

        #          name       ^        sup     _    sub or digit
    m = re.match('([^^_\d]+)(\^)?((?(2)[^_]+))(_)?((?(4).+|\d*))$', symb_name)
    if m is None:
        #print 'DON\'T MATCH'
        return symb_name

    #print m.groups()
    name, ssup, isup, ssub, isub = m.groups('')

    # let's prettify name
    gG = greek.get(name.lower())
    if gG is not None:
        name = name.islower() and gG[0] or gG[1]

    # let's pretty sup/sub
    psup = sup.get(isup)
    psub = sub.get(isub)

    res = name
    if psup is not None:
        res += psup
    else:
        res += '%s%s' % (ssup, isup)

    if psub is not None:
        res += psub
    else:
        res += '%s%s' % (ssub, isub)

    return res
