from __future__ import absolute_import

class Domain(object):
    def __init__(self):
        self.types = {}

    def register_type(
            self,
            type,
            head=None,
            child_getter=None,
            child_length=None, **_3to2kwargs
    ):
        if 'commutative' in _3to2kwargs: commutative = _3to2kwargs['commutative']; del _3to2kwargs['commutative']
        else: commutative = False
        if 'associative' in _3to2kwargs: associative = _3to2kwargs['associative']; del _3to2kwargs['associative']
        else: associative = False
        if 'arity' in _3to2kwargs: arity = _3to2kwargs['arity']; del _3to2kwargs['arity']
        else: arity = None
        self.types[type] = TypeInfo(head or type, child_getter, child_length, arity, associative, commutative)

    def is_symbol(self, expr):
        return type(expr) in self.types and self.types[type(expr)].child_getter is None

    def is_operation(self, expr):
        return type(expr) in self.types and self.types[type(expr)].child_getter is not None

    def is_associative_operation(self, expr):
        return type(expr) in self.types and self.types[type(expr)].is_associative(expr)

    def is_commutative_operation(self, expr):
        return type(expr) in self.types and self.types[type(expr)].is_associative(expr)


class TypeInfo(object):
    def __init__(self, head, child_getter, child_length, arity=None, associative=False, commutative=False):
        self.arity = arity
        self.head = head
        self.child_getter = child_getter
        self.child_length = child_length
        self.associative = associative
        self.commutative = commutative

    def get_head(self, expr):
        return self.head(expr) if callable(self.head) else self.head

    def get_children(self, expr):
        return self.child_getter(expr) if callable(self.child_getter) else None

    def get_child_count(self, expr):
        return self.child_length(expr) if callable(self.child_length) else None

    def get_arity(self, expr):
        return self.arity(expr) if callable(self.arity) else self.arity

    def is_associative(self, expr):
        return self.associative(expr) if callable(self.associative) else self.associative

    def is_commutative(self, expr):
        return self.commutative(expr) if callable(self.commutative) else self.commutative
