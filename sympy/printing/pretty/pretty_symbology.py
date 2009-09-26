"""Symbolic primitives + unicode/ascii abstraction  for pretty.py"""

import sys
warnings = ''

# first, setup unicodedate environment
try:
    import unicodedata

    # Python2.4 unicodedata misses some symbols, like subscript 'i', etc,
    # and we still want SymPy to be fully functional under Python2.4
    if sys.hexversion < 0x02050000:
        unicodedata_missing = {
            'GREEK SUBSCRIPT SMALL LETTER BETA' : u'\u1d66',
            'GREEK SUBSCRIPT SMALL LETTER GAMMA': u'\u1d67',
            'GREEK SUBSCRIPT SMALL LETTER RHO'  : u'\u1d68',
            'GREEK SUBSCRIPT SMALL LETTER PHI'  : u'\u1d69',
            'GREEK SUBSCRIPT SMALL LETTER CHI'  : u'\u1d6a',

            'LATIN SUBSCRIPT SMALL LETTER A'    : u'\u2090',
            'LATIN SUBSCRIPT SMALL LETTER E'    : u'\u2091',
            'LATIN SUBSCRIPT SMALL LETTER I'    : u'\u1d62',
            'LATIN SUBSCRIPT SMALL LETTER O'    : u'\u2092',
            'LATIN SUBSCRIPT SMALL LETTER R'    : u'\u1d63',
            'LATIN SUBSCRIPT SMALL LETTER U'    : u'\u1d64',
            'LATIN SUBSCRIPT SMALL LETTER V'    : u'\u1d65',
            'LATIN SUBSCRIPT SMALL LETTER X'    : u'\u2093',
        }
    else:
        unicodedata_missing = {}

    def U(name):
        """unicode character by name or None if not found"""
        try:
            u = unicodedata.lookup(name)
        except KeyError:
            u = unicodedata_missing.get(name)

            if u is None:
                global warnings
                warnings += 'W: no \'%s\' in unocodedata\n' % name

        return u

except ImportError:
    warnings += 'W: no unicodedata available\n'
    U = lambda name: None

from sympy.printing.conventions import split_super_sub
import re


# prefix conventions when constructing tables
# L   - LATIN     i
# G   - GREEK     beta
# D   - DIGIT     0
# S   - SYMBOL    +


__all__ = ['greek','sub','sup','xsym','vobj','hobj','pretty_symbol']


_use_unicode = False

def pretty_use_unicode(flag = None):
    """Set whether pretty-printer should use unicode by default"""
    global _use_unicode
    global warnings
    if flag is None:
        return _use_unicode

    if flag and warnings:
        # print warnings (if any) on first unicode usage
        print "I: pprint -- we are going to use unicode, but there are following problems:"
        print warnings
        warnings = ''

    use_unicode_prev = _use_unicode
    _use_unicode = flag
    return use_unicode_prev

def pretty_try_use_unicode():
    """See if unicode output is available and leverage it if possible"""

    try:
        symbols = []

        # see, if we can represent greek alphabet
        for g,G in greek.itervalues():
            symbols.append(g)
            symbols.append(G)

        # and atoms
        symbols += atoms_table.values()

        for s in symbols:
            if s is None:
                return  # common symbols not present!

            encoding = getattr(sys.stdout, 'encoding', None)

            # this happens when e.g. stdout is redirected through a pipe, or is
            # e.g. a cStringIO.StringO
            if encoding is None:
                return  # sys.stdout has no encoding

            # try to encode
            s.encode(encoding)

    except UnicodeEncodeError:
        pass
    else:
        pretty_use_unicode(True)


def xstr(*args):
    """call str or unicode depending on current mode"""
    if _use_unicode:
        return unicode(*args)
    else:
        return str(*args)

# COMPATIBILITY TWEAKS
def fixup_tables():
    # python2.4 unicodedata lacks some definitions

    for d in sub, sup:
        for k in d.keys():
            if d[k] is None:
                del d[k]


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
    '{}' :  'CURLY BRACKET',
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
    '(' :   (( EXT('('),  HUP('('), HLO('(') ),              '('),
    ')' :   (( EXT(')'),  HUP(')'), HLO(')') ),              ')'),
    '[' :   (( EXT('['),  CUP('['), CLO('[') ),              '['),
    ']' :   (( EXT(']'),  CUP(']'), CLO(']') ),              ']'),
    '{' :   (( EXT('{}'), HUP('{'), HLO('{'),  MID('{')  ),  '{'),
    '}' :   (( EXT('{}'), HUP('}'), HLO('}'),  MID('}')  ),  '}'),
    '|' :   U('BOX DRAWINGS LIGHT VERTICAL'),

    'lfloor' : (( EXT('['), EXT('['), CLO('[') ), U('LEFT FLOOR')),
    'rfloor' : (( EXT(']'), EXT(']'), CLO(']') ), U('RIGHT FLOOR')),
    'lceil'  : (( EXT('['), CUP('['), EXT('[') ), U('LEFT CEILING')),
    'rceil'  : (( EXT(']'), CUP(']'), EXT(']') ), U('RIGHT CEILING')),

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

# XXX this looks ugly
#   '[' :   (( '|', '-',  '-'   ),          '['),
#   ']' :   (( '|', '-',  '-'   ),          ']'),
# XXX not so ugly :(
    '[' :   (( '[', '[',  '['   ),          '['),
    ']' :   (( ']', ']',  ']'   ),          ']'),

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

    # TODO robustify when no unicodedat available
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
            # even height, but we have to print it somehow anyway...
            # XXX is it ok?
            length += 1

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
    """Construct vertical object of a given height

       see: xobj
    """
    return '\n'.join( xobj(symb, height) )

def hobj(symb, width):
    """Construct horizontal object of a given width

       see: xobj
    """
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


# atom symbols
_xsym = {
    '=='    : ( '=',    '='),
    '<'     : ( '<',    '<'),
    '<='    : ('<=',    U('LESS-THAN OR EQUAL TO')),
    '>='    : ('>=',    U('GREATER-THAN OR EQUAL TO')),
    '!='    : ('!=',    U('NOT EQUAL TO')),
    '*'     : ('*',     U('DOT OPERATOR')),
}


def xsym(sym):
    """get symbology for a 'character'"""
    op = _xsym[sym]

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
    'NegativeInfinity'  :   U('INFINITY') and ('-'+U('INFINITY')),  # XXX what to do here
    #'ImaginaryUnit'     :   U('GREEK SMALL LETTER IOTA'),
    #'ImaginaryUnit'     :   U('MATHEMATICAL ITALIC SMALL I'),
    'ImaginaryUnit'     :   U('DOUBLE-STRUCK ITALIC SMALL I'),
    'EmptySet'          :   U('EMPTY SET'),
    'Union'             :   U('UNION')
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

    name, sups, subs = split_super_sub(symb_name)

    # let's prettify name
    gG = greek.get(name.lower())
    if gG is not None:
        if name.islower():
            greek_name = greek.get(name.lower())[0]
        else:
            greek_name = greek.get(name.lower())[1]
        # some letters may not be available
        if greek_name is not None:
            name = greek_name

    # Let's prettify sups/subs. If it fails at one of them, pretty sups/subs are
    # not used at all.
    def pretty_list(l, mapping):
        result = []
        for s in l:
            pretty = mapping.get(s)
            if pretty is None:
                try: # match by separate characters
                    pretty = ''.join([mapping[c] for c in s])
                except KeyError:
                    return None
            result.append(pretty)
        return result

    pretty_sups = pretty_list(sups, sup)
    if pretty_sups is not None:
        pretty_subs = pretty_list(subs, sub)
    else:
        pretty_subs = None

    # glue the results into one string
    if pretty_subs is None: # nice formatting of sups/subs did not work
        if len(sups) > 0:
            sups_result = '^' + '^'.join(sups)
        else:
            sups_result = ''
        if len(subs) > 0:
            subs_result = '_' + '_'.join(subs)
        else:
            subs_result = ''
    else:
        sups_result = ' '.join(pretty_sups)
        subs_result = ' '.join(pretty_subs)


    return ''.join([name, sups_result, subs_result])


# final fixup
fixup_tables()
