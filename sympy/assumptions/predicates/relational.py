from sympy.assumptions import Predicate
from sympy.core.relational import Relational, Eq, Ne, Ge, Gt, Le, Lt
from sympy.multipledispatch import Dispatcher


class BinaryRelation(Predicate):
    """
    Base class for binary relation predicates
    """
    handler = None

    def __call__(self, *args):
        if len(args) != 2:
            raise ValueError("Binary relation takes two arguments, but got %s." % len(args))
        return super().__call__(*args)

    def eval(self, args, assumptions=True):
        rel = self.as_Relational()
        result = rel(*args).simplify()
        if isinstance(result, Relational):
            return None
        return result


class EqualityPredicate(BinaryRelation):
    """
    Binary predicate for ``==``

    """
    # TODO: Add examples
    name = 'eq'
    handler = Dispatcher("EqualityHandler", doc="Handler for key 'eq'.")

    def as_Relational(self):
        return Eq


class UnequalityPredicate(BinaryRelation):
    """
    Binary predicate for ``!=``.

    """
    # TODO: Add examples
    name = 'ne'
    handler = Dispatcher("UnequalityHandler", doc="Handler for key 'ne'.")

    def as_Relational(self):
        return Ne


class GreaterThanPredicate(BinaryRelation):
    """
    Binary predicate for ``>=``.

    """
    # TODO: Add examples
    name = 'ge'
    handler = Dispatcher("GreaterThanHandler", doc="Handler for key 'ge'.")

    def as_Relational(self):
        return Ge


class LessThanPredicate(BinaryRelation):
    """
    Binary predicate for ``<=``.

    """
    # TODO: Add examples
    name = 'le'
    handler = Dispatcher("LessThanHandler", doc="Handler for key 'le'.")

    def as_Relational(self):
        return Le


class StrictGreaterThanPredicate(BinaryRelation):
    """
    Binary predicate for ``>``.

    """
    # TODO: Add examples
    name = 'gt'
    handler = Dispatcher("StrictGreaterThanHandler", doc="Handler for key 'gt'.")

    def as_Relational(self):
        return Gt


class StrictLessThanPredicate(BinaryRelation):
    """
    Binary predicate for ``<``.

    """
    # TODO: Add examples
    name = 'lt'
    handler = Dispatcher("StrictLessThanHandler", doc="Handler for key 'lt'.")

    def as_Relational(self):
        return Lt
