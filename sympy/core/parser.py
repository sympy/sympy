"""
Defines parser for the following syntax rules:


<expr> = <lambda-test>
<lambda-test> = [ lambda <identifier-list> : ] <xor-test>
<xor-test> = [ <xor-test> <xor-op> ] <or-test>
<or-test> = [ <or-test> <or-op> ] <and-test>
<and_test> = [ <and-test> <and-op> ] <not-test>
<not-test> = [ <not-op> ] <relational>
<relational> = [ <arith> <rel-op> ] <arith>
<arith> = [ <arith> <add-op> ] <term>
                 | <factor>
<factor> = [ <add-op> ] <term>
<term> = [ <term> <mult-op> ] <power>
<power> = <primary> [ <power-op> <power> ]

<primary> = <atom>
            | <attr-ref>
            | <slicing>
            | <call>

<atom> = <identifier> | <literal> | <parenth>
<identifier> = <letter> [ <alphanumeric_character> ]...
<literal> = <int-literal>
            | <float-literal>
            | <logical-literal>
<parenth> = ( <expr-list> )

<attr-ref> = <primary> . <identifier>

<slicing> = <primary> [ <subscript-list> ]
<subscript> = <expr> | <slice>
<slice> = [ <expr> ] : [ <expr> ] [ : <expr> ]

<call> = <primary> ( [ <argument-list> ] )
<argument> = [ <identifier> = ] <expr>

"""

__all__ = ['Expr','Lambda_Test', 'XOr_Test', 'Or_Test','And_Test','Not_Test','Relational','Arith',
           'Factor','Term','Power','Primary','Identifier','Int_Literal',
           'Float_Literal','Logical_Literal','Parenth','Attr_Ref',
           'Slicing','Slice','Call','Argument','Subscript',
           'Subscript_List','Argument_List','Expr_List']

import re

from splitline import string_replace_map

from basic import Basic

from pattern_tools import Pattern

###############################################################################
############################### PATTERNS ######################################
###############################################################################

name = Pattern('<name>', r'[A-Z]\w*',flags=re.I)
digit_string = Pattern('<digit-string>',r'\d+')
significand = digit_string + '.' + ~digit_string | '.' + digit_string
exponent_letter = Pattern('<exponent-letter>',r'[E]',flags=re.I)
sign = Pattern('<sign>',r'[+-]')
signed_digit_string = ~sign + digit_string
exponent = signed_digit_string
float_string = significand + ~(exponent_letter + exponent) | \
               digit_string + exponent_letter + exponent

power_op = Pattern('<power-op>',r'(?<![*])[*]{2}(?![*])')
mult_op = Pattern('<mult-op>',r'(?<![*])[*](?![*])|(?<![/])[/](?![/])')
add_op = Pattern('<add-op>',r'[+-]')

rel_op = Pattern('<rel-op>',r'[=]{2}|[!][=]|[<][>]|[<][=]|[<]|[>][=]|[>]|\bin\b|\bnot\s+in\b')
not_op = Pattern('<not-op>',r'\bnot\b|[~]')
and_op = Pattern('<and-op>',r'\band\b|[&]')
or_op = Pattern('<or-op>',r'\bor\b|[|]')
xor_op = Pattern('<xor-op>',r'\bxor\b|[\^]')

lambda_name = Pattern('<lambda>',r'\blambda\b')
abs_lambda = abs(lambda_name)

int_literal_constant_named = signed_digit_string.named('value')
abs_int_literal_constant_named = abs(int_literal_constant_named)

float_literal_constant = significand + ~(exponent_letter + exponent) |\
                         digit_string + exponent_letter + exponent
signed_float_literal_constant = ~sign + float_literal_constant
float_literal_constant_named = signed_float_literal_constant.named('value')
abs_float_literal_constant_named = abs(float_literal_constant_named)

logical_literal_constant_named = Pattern('<value>',r'\bTrue\b|\bFalse\b', flags=re.I).named()
abs_logical_literal_constant_named = abs(logical_literal_constant_named)

abs_name = abs(name)


###############################################################################
############################### BASE CLASSES ##################################
###############################################################################


class NoMatchError(Exception):
    pass


class Base(object):
    """ Base class for symbolic syntax rules.

    All Base classes have the following attributes:
      .string - original string to construct a class instance
    """
    subclasses = {}
    
    def __new__(cls, string):
        match = cls.__dict__.get('match', None)
        errmsg = '%s: %r' % (cls.__name__, string)
        if match is not None:
            try:
                result = cls.match(string)
            except NoMatchError, msg:
                result = None
        else:
            result = None
        if isinstance(result, tuple):
            obj = object.__new__(cls)
            obj.string = string
            if hasattr(cls, 'init'): obj.init(*result)
            return obj
        elif isinstance(result, Base):
            return result
        elif result is None:
            for subcls in Base.subclasses.get(cls.__name__,[]):
                try:
                    obj = subcls(string)
                except NoMatchError, msg:
                    obj = None
                if obj is not None:
                    return obj
        raise NoMatchError,errmsg

    def init(self, *items):
        self.items = items
        return

    def torepr(self):
        return '%s(%s)' % (self.__class__.__name__, ', '.join(map(repr,self.items)))

    def compare(self, other):
        return cmp(self.items,other.items)

    def __str__(self): return self.tostr()

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, str(self))

    def __cmp__(self, other):
        if self is other: return 0
        if not isinstance(other, self.__class__): return cmp(self.__class__, other.__class__)
        return self.compare(other)

    def tosymbolic(self, commutative=True):
        raise NotImplementedError,`self`


class SequenceBase(Base):
    """
    <sequence-base> = <obj>, <obj> [ , <obj> ]...
    """

    def match(separator, subcls, string):
        line, repmap = string_replace_map(string)
        if isinstance(separator, str):
            splitted = line.split(separator)
        else:
            splitted = separator[1].split(line)
            separator = separator[0]
        if len(splitted)<=1: return
        lst = []
        for p in splitted:
            lst.append(subcls(repmap(p.strip())))
        return separator, tuple(lst)
    match = staticmethod(match)

    def init(self, separator, items):
        self.separator = separator
        self.items = items
        return

    def tostr(self):
        s = self.separator
        if s==',': s = s + ' '
        elif s==' ': pass
        else: s = ' ' + s + ' '
        return s.join(map(str, self.items))

    def torepr(self):
        return '%s(%r, %r)' % (self.__class__.__name__, self.separator, self.items)

    def compare(self, other):
        return cmp((self.separator,self.items),(other.separator,self.items))

    def tosymbolic(self, commutative=True):
        return [item.tosymbolic(commutative) for item in self.items]


class UnaryOpBase(Base):
    """
    <unary-op-base> = <unary-op> <rhs>
    """

    def tostr(self):
        return '%s %s' % tuple(self.items)

    def match(op_pattern, rhs_cls, string):
        m = op_pattern.match(string)
        if not m: return
        rhs = string[m.end():].lstrip()
        if not rhs: return
        op = string[:m.end()].rstrip()
        return op, rhs_cls(rhs)
    match = staticmethod(match)

    def tosymbolic(self, commutative=True):
        op = self.items[0]
        mth = {'+':'__pos__','-':'__neg__',
               '~':'__invert__','not':'__invert__',
               }.get(op, None)
        if mth is not None:
            rhs = self.items[1].tosymbolic(commutative)
            return getattr(rhs,mth)()
        return Base.tosymbolic(self, commutative)


class BinaryOpBase(Base):
    """
    <binary-op-base> = <lhs> <op> <rhs>
    <op> is searched from right by default.
    """

    def match(lhs_cls, op_pattern, rhs_cls, string, right=True):
        line, repmap = string_replace_map(string)
        if isinstance(op_pattern, str):
            if right:
                t = line.rsplit(op_pattern,1)
            else:
                t = line.split(op_pattern,1)
            if len(t)!=2: return
            lhs, rhs = t[0].rstrip(), t[1].lstrip()
            op = op_pattern
        else:
            if right:
                t = op_pattern.rsplit(line)
            else:
                t = op_pattern.lsplit(line)
            if t is None or len(t)!=3: return
            lhs, op, rhs = t
            lhs = lhs.rstrip()
            rhs = rhs.lstrip()
            op = op
        if not lhs: return
        if not rhs: return
        lhs_obj = lhs_cls(repmap(lhs))
        rhs_obj = rhs_cls(repmap(rhs))
        return lhs_obj, op, rhs_obj
    match = staticmethod(match)

    def tostr(self):
        return '%s %s %s' % tuple(self.items)

    def tosymbolic(self, commutative=True):
        op = self.items[1]
        mth = {'+':'__add__','-':'__sub__',
               '*':'__mul__','/':'__div__',
               '**':'__pow__','>':'__gt__',
               '<':'__lt__','>=':'__ge__',
               '<=':'__le__','==':'__eq__',
               '<>':'__ne__','!=':'__ne__',
               'and':'__and__','&':'__and__',
               'or':'__or__','|':'__or__',
               'xor':'__xor__','^':'__xor__',
               }.get(op, None)
        if not commutative:
            if mth=='__mul__': mth = 'ncmul'
            if mth=='__div__': mth = 'ncdiv'

        if mth is not None:
            lhs = self.items[0].tosymbolic(commutative)
            rhs = self.items[2].tosymbolic(commutative)
            return getattr(lhs,mth)(rhs)
        return Base.tosymbolic(self. commutative)


class KeywordValueBase(Base):
    """
    <keyword-value-base> = [ <lhs> = ] <rhs>
    """

    def match(lhs_cls, rhs_cls, string, require_lhs = True, upper_lhs = False):
        if require_lhs and '=' not in string: return
        if isinstance(lhs_cls, (list, tuple)):
            for s in lhs_cls:
                try:
                    obj = KeywordValueBase.match(s, rhs_cls, string,
                                                 require_lhs=require_lhs, upper_lhs=upper_lhs)
                except NoMatchError:
                    obj = None
                if obj is not None: return obj
            return obj
        lhs,rhs = string.split('=',1)
        lhs = lhs.rstrip()
        rhs = rhs.lstrip()
        if not rhs: return
        if not lhs:
            if require_lhs: return
            return None, rhs_cls(rhs)
        if isinstance(lhs_cls, str):
            if upper_lhs:
                lhs = lhs.upper()
            if lhs_cls!=lhs: return
            return lhs, rhs_cls(rhs)
        return lhs_cls(lhs),rhs_cls(rhs)
    match = staticmethod(match)

    def tostr(self):
        if self.items[0] is None: return str(self.items[1])
        return '%s = %s' % tuple(self.items)


class BracketBase(Base):
    """
    <bracket-base> = <left-bracket-base> <something> <right-bracket>
    """

    def match(brackets, cls, string, require_cls=True):
        i = len(brackets)/2
        left = brackets[:i]
        right = brackets[-i:]
        if string.startswith(left) and string.endswith(right):
            line = string[i:-i].strip()
            if not line:
                if require_cls:
                    return
                return left,None,right
            return left,cls(line),right
        return
    match = staticmethod(match)

    def tostr(self):
        if self.items[1] is None:
            return '%s%s' % (self.items[0], self.items[2])
        return '%s%s%s' % tuple(self.items)


class CallBase(Base):
    """
    <call-base> = <lhs> ( [ <rhs> ] )
    """

    def match(lhs_cls, rhs_cls, string, upper_lhs = False, require_rhs=False, parens='()'):
        if not string.endswith(parens[1]): return
        line, repmap = string_replace_map(string)
        i = line.rfind(parens[0])
        if i==-1: return
        lhs = line[:i].rstrip()
        if not lhs: return
        rhs = line[i+1:-1].strip()
        lhs = repmap(lhs)
        if upper_lhs:
            lhs = lhs.upper()
        rhs = repmap(rhs)
        if isinstance(lhs_cls, str):
            if lhs_cls!=lhs: return
        else:
            lhs = lhs_cls(lhs)
        if rhs:
            if isinstance(rhs_cls, str):
                if rhs_cls!=rhs: return
            else:
                rhs = rhs_cls(rhs)
            return lhs, rhs
        elif require_rhs:
            return
        return lhs, None
    match = staticmethod(match)

    def tostr(self):
        if self.items[1] is None: return '%s()' % (self.items[0])
        return '%s(%s)' % (self.items[0], self.items[1])


class StringBase(Base):
    """
    <string-base> = <xyz>
    """

    def match(pattern, string):
        if isinstance(pattern, (list,tuple)):
            for p in pattern:
                obj = StringBase.match(p, string)
                if obj is not None: return obj
            return
        if isinstance(pattern, str):
            if len(pattern)==len(string) and pattern==string: return string,
            return
        if pattern.match(string): return string,
        return
    match = staticmethod(match)

    def init(self, string):
        self.string = string
        return

    def tostr(self): return str(self.string)

    def torepr(self): return '%s(%r)' % (self.__class__.__name__, self.string)

    def compare(self, other):
        return cmp(self.string,other.string)


class NumberBase(Base):
    """
    <number-base> = <number>
    """

    def match(number_pattern, string):
        m = number_pattern.match(string)
        if m is None: return
        return m.group('value').upper(),
    match = staticmethod(match)

    def tostr(self):
        return str(self.items[0])

    def compare(self, other):
        return cmp(self.items[0], other.items[0])

###############################################################################
############################## USEFUL CLASSES #################################
###############################################################################


class Expr(Base):
    """
    <expr> = <lambda-test>
    """
    subclass_names = ['Lambda_Test']

    def match(string): return
    match = staticmethod(match)

class Lambda_Test(Base):
    """
    <lambda-test> = [ lambda <identifier-list> : ] <or-test>
    """
    
    subclass_names = ['Or_Test']
    use_names = ['Identifier_List']

    def match(string):
        if not string.startswith('lambda'): return
        i = string.find(':')
        if i==-1: return
        return Identifier_List(string[6:i].strip()), Or_Test(string[i+1:].lstrip())
    match = staticmethod(match)

    def tostr(self):
        return 'lambda %s : %s' % tuple(self.items)

    def tosymbolic(self, commutative=True):
        ids = self.items[0].tosymbolic(commutative)
        if not isinstance(ids, list):
            ids = [ids]
        expr = self.items[1].tosymbolic(commutative)
        return Basic.Lambda(expr, *ids)


        
class Or_Test(BinaryOpBase):
    """
    <or-test> = [ <or-test> <or-op> ] <XOr_Test>
    <or-op>  = or | \|
    """
    
    subclass_names = ['XOr_Test']
    use_names = ['Or_Test']

    def match(string):
        return BinaryOpBase.match(Or_Test, or_op.named(), XOr_Test, string)
    match = staticmethod(match)

class XOr_Test(BinaryOpBase):
    """
    <xor-test> = [ <xor-test> <xor-op> ] <And_Test>
    <xor-op>  = xor | ^
    """
    
    subclass_names = ['And_Test']
    use_names = ['XOr_Test']

    def match(string):
        return BinaryOpBase.match(XOr_Test, xor_op.named(), And_Test, string)
    match = staticmethod(match)

class And_Test(BinaryOpBase):
    """
    <And_Test> = [ <And_Test> <and-op> ] <not-test>
    <and-op> = and | &
    """
    subclass_names = ['Not_Test']
    use_names = ['And_Test','Not_Test']

    def match(string):
        return BinaryOpBase.match(And_Test,and_op.named(),Not_Test,string)
    match = staticmethod(match)


class Not_Test(UnaryOpBase):
    """
    <not-test> = [ <not-op> ] <Relational>
    <not-op> = not | ~
    """
    subclass_names = ['Relational']
    use_names = []

    def match(string):
        return UnaryOpBase.match(not_op.named(),Relational,string)
    match = staticmethod(match)


class Relational(BinaryOpBase):
    """
    <Relational> = [ <arith> <rel-op> ] <arith>
    <rel-op> = == | <> | != | < | <= | > | >= | in | not in
    """
    subclass_names = ['Arith']
    use_names = []

    def match(string):
        return BinaryOpBase.match(Arith,rel_op.named(),Arith,string)
    match = staticmethod(match)


class Arith(BinaryOpBase):
    """
    <arith> = [ <arith> <add-op> ] <term>
              | <factor>
    <add-op>   = +
                 | -
    """
    subclass_names = ['Factor']
    use_names = ['Arith']

    def match(string):
        return BinaryOpBase.match(Arith,add_op.named(),Term,string)
    match = staticmethod(match)


class Factor(UnaryOpBase):
    """
    <factor> = [ <add-op> ] <term>
    <add-op>   = +
                 | -
    """
    subclass_names = ['Term']
    use_names = []

    def match(string):
        return UnaryOpBase.match(add_op.named(),Term,string)
    match = staticmethod(match)


class Term(BinaryOpBase):
    """
    <term> = [ <term> <mult-op> ] <power>
    <mult-op>  = *
                 | /
    """
    subclass_names = ['Power']
    use_names = ['Term','Power']

    def match(string):
        return BinaryOpBase.match(Term,mult_op.named(),Power,string)
    match = staticmethod(match)


class Power(BinaryOpBase): # R704
    """
    <power> = <primary> [ <power-op> <power> ]
    <power-op> = **
    """
    subclass_names = ['Primary']
    use_names = ['Power']

    def match(string):
        returnv = BinaryOpBase.match(Primary,power_op.named(),Power,string,right=False)
        return returnv
    match = staticmethod(match)


class Primary(Base): # R701
    """
    <primary> = <atom>
                | <attrref>
                | <slicing>
                | <call>
    <atom> = <identifier> | <literal> | <parenth>
    <literal> = <int-literal>
            | <float-literal>
            | <logical-literal>
    """
    subclass_names = [
                      'Attr_Ref',
                      'Slicing',
                      'Call',
                      'Int_Literal',
                      'Float_Literal',
                      'Logical_Literal',
                      'Parenth',
                      'Identifier',
                      ]


class Parenth(BracketBase):
    """
    <parenth> = ( <expr-list> )
    """
    subclass_names = []
    use_names = ['Expr_List']

    def match(string):
        return BracketBase.match('()', Expr_List, string)
    match = staticmethod(match)

    def tosymbolic(self, commutative=True):
        return self.items[1].tosymbolic(commutative)


class Identifier(StringBase):
    """
    <identifier> = <letter> [ <alphanumeric_character> ]...
    """
    subclass_names = []

    def match(string):
        return StringBase.match(abs_name, string)
    match = staticmethod(match)

    def tosymbolic(self, commutative=True):
        obj = Basic.singleton.get(self.string)
        if obj is not None:
            return obj()
        return Basic.Symbol(self.string)


class Int_Literal(NumberBase):
    """
    <int-literal> = [ <sign> ] <digit-string>
    """
    subclass_names = []

    def match(string):
        return NumberBase.match(abs_int_literal_constant_named, string)
    match = staticmethod(match)

    def tosymbolic(self, commutative=True):
        return Basic.Integer(int(self.items[0]))


class Float_Literal(NumberBase):
    """
    """
    subclass_names = []

    def match(string):
        return NumberBase.match(abs_float_literal_constant_named, string)
    match = staticmethod(match)

    def tosymbolic(self, commutative=True):
        return Basic.Real(self.items[0])

class Logical_Literal(NumberBase):
    """
    <logical-literal-constant> = True | False
    """
    subclass_names = []

    def match(string):
        m = abs_logical_literal_constant_named.match(string)
        if m is None: return
        return m.group('value'),
    match = staticmethod(match)

    def tosymbolic(self, commutative=True):
        return Basic.Symbol(self.items[0].upper())


class Attr_Ref(Base):
    """
    <attr-ref> = <primary> . <identifier>
    """
    subclass_names = []
    use_names = ['Primary','Identifier']

    def match(string):
        if '.' not in string: return
        i = string.rindex('.')
        attrname = Identifier(string[i+1:].lstrip())
        return Primary(string[:i].rstrip()), attrname
    match = staticmethod(match)

    def tostr(self):
        return '%s.%s' % (self.items[0], self.items[1])

    def tosymbolic(self, commutative=True):
        lhs = self.items[0].tosymbolic(commutative)
        rhs = self.items[1]
        return getattr(lhs, rhs.string)


class Slicing(CallBase):
    """
    <slicing> = <primary> <left-square-bracket> <subscript-list> <right-square-bracket>
    """
    subclass_names = []
    use_names = ['Primary','Subscript_List']

    def match(string):
        return CallBase.match(Primary, Subscript_List, string, require_rhs=True, parens = '[]')
    match = staticmethod(match)

    def tostr(self):
        return '%s[%s]' % (self.items[0], self.items[1])

    def tosymbolic(self, commutative=True):
        lhs = self.items[0].tosymbolic(commutative)
        rhs = self.items[1].tosymbolic(commutative)
        return lhs[rhs]


class Subscript(Base):
    """
    <subscript> = <expr> | <slice>
    """
    subclass_names = ['Slice','Expr']


class Slice(Base):
    """
    <slice> = [ <expr> ] : [ <expr> ] [ : <expr> ]
    """
    subclass_names = []
    use_names = ['Expr']

    def match(string):
        line, repmap = string_replace_map(string)
        t = line.split(':')
        if len(t)<=1 or len(t)>3: return
        lhs_obj,rhs_obj, stride_obj = None, None, None
        if len(t)==2:
            lhs,rhs = t[0].rstrip(),t[1].lstrip()
        else:
            lhs,rhs,stride = t[0].rstrip(),t[1].strip(),t[2].lstrip()
            if stride:
                stride_obj = Expr(repmap(stride))
        if lhs:
            lhs_obj = Expr(repmap(lhs))
        if rhs:
            rhs_obj = Expr(repmap(rhs))
        return lhs_obj, rhs_obj, stride_obj
    match = staticmethod(match)

    def tostr(self):
        s = ''
        if self.items[0] is not None:
            s += str(self.items[0]) + ' :'
        else:
            s += ':'
        if self.items[1] is not None:
            s += ' ' + str(self.items[1])
        if self.items[2] is not None:
            s += ' : ' + str(self.items[2])
        return s

    def tosymbolic(self, commutative=True):
        s = []
        for i in self.items:
            if i is None:
                s.append(i)
            else:
                s.append(i.tosymbolic(commutative))
        return slice(*s)


class Call(CallBase):
    """
    <call> = <primary> ( [ <argument-list> ] )
    """
    subclass_names = []
    use_names = ['Primary', 'Argument_List']

    def match(string):
        return CallBase.match(Primary, Argument_List, string)
    match = staticmethod(match)

    def tosymbolic(self, commutative=True):
        lhs = self.items[0].tosymbolic(commutative = False)
        rhs = self.items[1]
        if rhs is None: return lhs()
        rhs = rhs.tosymbolic(commutative)
        if isinstance(rhs, list):
            return lhs(*rhs)
        return lhs(rhs)

class Argument(KeywordValueBase):
    """
    <argument> = [ <identifier> = ] <expr>
    """
    subclass_names = ['Expr']
    use_names = ['Identifier']

    def match(string):
        return KeywordValueBase.match(Identifier, Expr, string)
    match = staticmethod(match)

    def tosymbolic(self, commutative=True):
        s = self.items[0].tosymbolic(commutative)
        e = self.items[1].tosymbolic(commutative)
        if isinstance(e,(list, tuple)):
            return Basic.Range(s, *e)
        return Basic.Keyword(s, e)

###############################################################################
################ GENERATE Scalar_, _List, _Name CLASSES #######################
###############################################################################

ClassType = type(Base)
_names = dir()
for clsname in _names:
    cls = eval(clsname)
    if not (isinstance(cls, ClassType) and issubclass(cls, Base) and not cls.__name__.endswith('Base')): continue
    names = getattr(cls, 'subclass_names', []) + getattr(cls, 'use_names', [])
    for n in names:
        if n in _names: continue
        if n.endswith('_List'):
            _names.append(n)
            n = n[:-5]
            #print 'Generating %s_List' % (n)
            exec '''\
class %s_List(SequenceBase):
    subclass_names = [\'%s\']
    use_names = []
    def match(string): return SequenceBase.match(r\',\', %s, string)
    match = staticmethod(match)
''' % (n, n, n)

Base_classes = {}
for clsname in dir():
    cls = eval(clsname)
    if isinstance(cls, ClassType) and issubclass(cls, Base) and not cls.__name__.endswith('Base'):
        Base_classes[cls.__name__] = cls

###############################################################################
##################### OPTIMIZE subclass_names tree ############################
###############################################################################

if 1: # Optimize subclass tree:

    def _rpl_list(clsname):
        if not Base_classes.has_key(clsname):
            print 'Not implemented:',clsname
            return [] # remove this code when all classes are implemented
        cls = Base_classes[clsname]
        if cls.__dict__.has_key('match'): return [clsname]
        l = []
        for n in getattr(cls,'subclass_names',[]):
            l1 = _rpl_list(n)
            for n1 in l1:
                if n1 not in l:
                    l.append(n1)
        return l

    for cls in Base_classes.values():
        if not hasattr(cls, 'subclass_names'): continue
        opt_subclass_names = []
        for n in cls.subclass_names:
            for n1 in _rpl_list(n):
                if n1 not in opt_subclass_names:  opt_subclass_names.append(n1)
        if not opt_subclass_names==cls.subclass_names:
            #print cls.__name__,':',', '.join(cls.subclass_names),'->',', '.join(opt_subclass_names)
            cls.subclass_names[:] = opt_subclass_names
        #else:
        #    print cls.__name__,':',opt_subclass_names


# Initialize Base.subclasses dictionary:
for clsname, cls in Base_classes.items():
    subclass_names = getattr(cls, 'subclass_names', None)
    if subclass_names is None:
        print '%s class is missing subclass_names list' % (clsname)
        continue
    try:
        l = Base.subclasses[clsname]
    except KeyError:
        Base.subclasses[clsname] = l = []
    for n in subclass_names:
        if Base_classes.has_key(n):
            l.append(Base_classes[n])
        else:
            print '%s not implemented needed by %s' % (n,clsname)

#EOF
