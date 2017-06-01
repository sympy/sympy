from __future__ import absolute_import
from .expressions import Operation, Wildcard, AssociativeOperation, CommutativeOperation, SymbolWildcard, Pattern

__all__ = [
    u'is_constant', u'is_syntactic', u'get_head', u'match_head', u'preorder_iter', u'preorder_iter_with_position',
    u'is_anonymous', u'contains_variables_from_set', u'register_operation_factory', u'create_operation_expression'
]


def is_constant(expression):
    u"""Check if the given expression is constant, i.e. it does not contain Wildcards."""
    if isinstance(expression, Operation):
        return all(is_constant(o) for o in expression)
    return not isinstance(expression, Wildcard)


def is_syntactic(expression):
    u"""
    Check if the given expression is syntactic, i.e. it does not contain sequence wildcards or
    associative/commutative operations.
    """
    if isinstance(expression, (AssociativeOperation, CommutativeOperation)):
        return False
    if isinstance(expression, Operation):
        return all(is_syntactic(o) for o in expression)
    if isinstance(expression, Wildcard):
        return expression.fixed_size
    return True


def get_head(expression):
    u"""Returns the given expression's head."""
    if isinstance(expression, Wildcard):
        if isinstance(expression, SymbolWildcard):
            return expression.symbol_type
        return None
    return type(expression)


def match_head(subject, pattern):
    u"""Checks if the head of subject matches the pattern's head."""
    if isinstance(pattern, Pattern):
        pattern = pattern.expression
    pattern_head = get_head(pattern)
    if pattern_head is None:
        return True
    subject_head = get_head(subject)
    assert subject_head is not None
    return issubclass(subject_head, pattern_head)


def preorder_iter(expression):
    u"""Iterate over the expression in preorder."""
    yield expression
    if isinstance(expression, Operation):
        for operand in expression:
            for i in preorder_iter(operand):
                yield i


def preorder_iter_with_position(expression):
    u"""Iterate over the expression in preorder.

    Also yields the position of each subexpression.
    """
    yield expression, ()
    if isinstance(expression, Operation):
        for i, operand in enumerate(expression):
            for child, pos in preorder_iter_with_position(operand):
                yield child, (i, ) + pos


def is_anonymous(expression):
    u"""Returns True iff the expression does not contain any variables."""
    if hasattr(expression, u'variable_name') and expression.variable_name:
        return False
    if isinstance(expression, Operation):
        return all(is_anonymous(o) for o in expression)
    return True


def contains_variables_from_set(expression, variables):
    u"""Returns True iff the expression contains any of the variables from the given set."""
    if hasattr(expression, u'variable_name') and expression.variable_name in variables:
        return True
    if isinstance(expression, Operation):
        return any(contains_variables_from_set(o, variables) for o in expression)
    return False


def simple_operation_factory(op, args, variable_name):
    return type(op)(args)


_operation_factories = {
    list: simple_operation_factory,
    tuple: simple_operation_factory,
    set: simple_operation_factory,
    frozenset: simple_operation_factory,
    # TODO: Add support for dicts
}


def register_operation_factory(operation, factory):
    _operation_factories[operation] = factory


def create_operation_expression(old_operation, new_operands, variable_name=True):
    operation = type(old_operation)
    for parent in operation.__mro__:
        if parent in _operation_factories:
            return _operation_factories[parent](old_operation, new_operands, variable_name)
    if variable_name is True:
        variable_name = getattr(old_operation, u'variable_name', None)
    if variable_name is False:
        return operation(*new_operands)
    return operation(*new_operands, variable_name=variable_name)
