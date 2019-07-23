"""
The classes used here are for the internal use of assumptions system
only and should not be used anywhere else as these don't posses the
signatures common to SymPy objects. For general use of logic constructs
please refer to sympy.logic classes And, Or, Not, etc.
"""
from itertools import combinations, product
from sympy import S
from sympy.logic.boolalg import BooleanFunction, Or, And, Not


# Early encoding for literals
literals_store = dict()


class Literal(object):
    """
    The smallest element of a CNF object
    """
    literals_count = 0

    def __new__(cls, lit, is_not=False):
        if type(lit).__name__ == 'Not':
            lit = lit.args[0]
            is_not = True
        k = literals_store.get(lit, None)
        if k is not None:
            code = k
        elif lit in (S.false, False):
            code = 0
        else:
            cls.literals_count += 1
            code = cls.literals_count
        obj = super(Literal, cls).__new__(cls)
        obj.code = code
        obj.lit = lit
        literals_store[lit] = obj.code
        obj.is_Not = is_not
        return obj

    @property
    def arg(self):
        sign = -1 if self.is_Not else 1
        return sign*self.code

    def __invert__(self):
        is_not = not self.is_Not
        return Literal(self.lit, is_not=is_not)

    def __str__(self):
        return '%s(%s, %s)' % (type(self).__name__, self.lit, self.is_Not)

    __repr__ = __str__

    def __eq__(self, other):
        return self.arg == other.arg

    def __hash__(self):
        h = hash((type(self).__name__, + self.arg))
        return h

    def __lt__(self, other):
        return self.arg < other.arg


def literal_sort_key(lit):
    _type = type(lit).__name__
    if _type == 'Literal':
        return lit.arg

    if _type == 'AND':
        return -2, lit.args

    if _type == 'OR':
        return -1, lit.args


class OR(object):
    """
    A low-level implementation for Or
    """
    def __init__(self, *args):
        self._args = args

    @property
    def args(self):
        return sorted(self._args, key=literal_sort_key)

    def __invert__(self):
        return AND(*[~arg for arg in self._args])

    def __hash__(self):
        return hash((type(self).__name__,) + tuple(self.args))

    def __eq__(self, other):
        return self.args == other.args

    def __str__(self):
        s = '(' + ' | '.join([str(arg) for arg in self.args]) + ')'
        return s

    __repr__ = __str__


class AND(object):
    """
    A low-level implementation for And
    """
    def __init__(self, *args):
        self._args = args

    def __invert__(self):
        return OR(*[~arg for arg in self._args])

    @property
    def args(self):
        return sorted(self._args, key=literal_sort_key)

    def __hash__(self):
        return hash((type(self).__name__,) + tuple(self.args))

    def __eq__(self, other):
        return self.args == other.args

    def __str__(self):
        s = '('+' & '.join([str(arg) for arg in self.args])+')'
        return s

    __repr__ = __str__


def to_NNF(expr):
    """
    Generates the Negation Normal Form of any boolean expression in terms
    of AND, OR, and Literal objects.
    """
    if not isinstance(expr, BooleanFunction):
        return Literal(expr)

    klass = type(expr).__name__
    if klass == 'Not':
        arg = expr.args[0]
        tmp = to_NNF(arg)  # Strategy: negate the NNF of expr
        return ~tmp

    if klass == 'Or':
        return OR(*[to_NNF(x) for x in Or.make_args(expr)])

    if klass == 'And':
        return AND(*[to_NNF(x) for x in And.make_args(expr)])

    if klass == 'Nand':
        tmp = AND(*[to_NNF(x) for x in expr.args])
        return ~tmp

    if klass == 'Nor':
        tmp = OR(*[to_NNF(x) for x in expr.args])
        return ~tmp

    if klass == 'Xor':
        cnfs = []
        for i in range(0, len(expr.args) + 1, 2):
            for neg in combinations(expr.args, i):
                clause = [~to_NNF(s) if s in neg else to_NNF(s)
                          for s in expr.args]
                cnfs.append(OR(*clause))
        return AND(*cnfs)

    if klass == 'Xnor':
        cnfs = []
        for i in range(0, len(expr.args) + 1, 2):
            for neg in combinations(expr.args, i):
                clause = [~to_NNF(s) if s in neg else to_NNF(s)
                          for s in expr.args]
                cnfs.append(OR(*clause))
        return ~AND(*cnfs)

    if klass == 'Implies':
        L, R = map(to_NNF, expr.args)
        return OR(~L, R)

    if klass == 'Equivalent':
        cnfs = []
        for a, b in zip(expr.args, expr.args[1:]):
            a = to_NNF(a)
            b = to_NNF(b)
            cnfs.append(OR(~a, b))
        a = to_NNF(expr.args[-1])
        b = to_NNF(expr.args[0])
        cnfs.append(OR(~a, b))
        return AND(*cnfs)

    if klass == 'ITE':
        L, M, R = map(CNF.to_CNF, expr.args)
        return AND(OR(~L, M), OR(L, R))


def distribute_AND_over_OR(expr):
    """
    Distributes AND over OR in the NNF expression.
    Returns the result( Conjunctive Normal Form of expression)
    as a CNF object.
    """
    if not isinstance(expr, (AND, OR)):
        tmp = set()
        tmp.add(frozenset((expr,)))
        return CNF(tmp)

    klass = type(expr).__name__
    if klass == 'OR':
        return CNF.all_or(*[distribute_AND_over_OR(arg)
                            for arg in expr._args])

    if klass == 'AND':
        return CNF.all_and(*[distribute_AND_over_OR(arg)
                             for arg in expr._args])


class CNF(object):
    """
    Class to represent CNF of a Boolean expression.
    Consists of set of clauses, which themselves are stored as
    frozenset of Literal objects.
    """
    def __init__(self, clauses=None):
        if not clauses:
            clauses = set()
        self.clauses = clauses

    def add(self, prop):
        clauses = CNF.to_CNF(prop).clauses
        self.add_clauses(clauses)

    def __str__(self):
        s = ' & '.join(
            ['(' + ' | '.join([str(lit) for lit in clause]) +')'
            for clause in self.clauses]
        )
        return s

    def extend(self, props):
        for p in props:
            self.add(p)
        return self

    def copy(self):
        return CNF(set(self.clauses))

    def add_clauses(self, clauses):
        self.clauses |= clauses

    @classmethod
    def from_prop(cls, prop):
        res = CNF()
        res.add(prop)
        return res

    def __iand__(self, other):
        self.add_clauses(other.clauses)
        return self

    def all_predicates(self):
        predicates = set()
        for c in self.clauses:
            predicates |= {arg.lit for arg in c}
        return predicates

    def _or(self, cnf):
        clauses = set()
        for a, b in product(self.clauses, cnf.clauses):
            tmp = set(a)
            for t in b:
                tmp.add(t)
            clauses.add(frozenset(tmp))
        return CNF(clauses)

    def _and(self, cnf):
        clauses = self.clauses.union(cnf.clauses)
        return CNF(clauses)

    def _not(self):
        clss = list(self.clauses)
        ll = set()
        for x in clss[-1]:
            ll.add(frozenset((~x,)))
        ll = CNF(ll)

        for rest in clss[:-1]:
            p = set()
            for x in rest:
                p.add(frozenset((~x,)))
            ll = ll._or(CNF(p))
        return ll

    @classmethod
    def all_or(cls, *cnfs):
        b = cnfs[0].copy()
        for rest in cnfs[1:]:
            b = b._or(rest)
        return b

    @classmethod
    def all_and(cls, *cnfs):
        b = cnfs[0].copy()
        for rest in cnfs[1:]:
            b = b._and(rest)
        return b

    @classmethod
    def to_CNF(cls, expr):
        expr = to_NNF(expr)
        expr = distribute_AND_over_OR(expr)
        return expr

    @classmethod
    def CNF_to_cnf(cls, cnf):
        """
        Converts CNF object to SymPy's boolean expression
        retaining the form of expression.
        """
        def remove_literal(arg):
            return Not(arg.lit) if arg.is_Not else arg.lit

        return And(*(Or(*(remove_literal(arg) for arg in clause)) for clause in cnf.clauses))


class EncodedCNF(object):
    """
    Class for encoding the CNF expression.
    """
    def __init__(self, data=None, encoding=None):
        if not data and not encoding:
            data = list()
            encoding = dict()
        self.data = data
        self.encoding = encoding
        self._symbols = list(encoding.keys())

    def from_cnf(self, cnf):
        self._symbols = list(cnf.all_predicates())
        n = len(self._symbols)
        self.encoding = dict(list(zip(self._symbols, list(range(1, n + 1)))))
        self.data = [self.encode(clause) for clause in cnf.clauses]

    @property
    def symbols(self):
        return self._symbols

    @property
    def variables(self):
        return range(1, len(self._symbols) + 1)

    def copy(self):
        new_data = [set(clause) for clause in self.data]
        return EncodedCNF(new_data, dict(self.encoding))

    def add_prop(self, prop):
        cnf = CNF.from_prop(prop)
        self.add_from_cnf(cnf)

    def add_from_cnf(self, cnf):
        clauses = [self.encode(clause) for clause in cnf.clauses]
        self.data += clauses

    def encode_arg(self, arg):
        literal = arg.lit
        value = self.encoding.get(literal, None)
        if value is None:
            n = len(self._symbols)
            self._symbols.append(literal)
            value = self.encoding[literal] = n + 1
        if arg.is_Not:
            return -value
        else:
            return value

    def encode(self, clause):
        return {self.encode_arg(arg) if arg.code else arg.code for arg in clause}
