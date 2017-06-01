# -*- coding: utf-8 -*-
u"""This module contains various functions for working with expressions.

- With `substitute()` you can replace occurrences of variables with an expression or sequence of expressions.
- With `replace()` you can replace a subexpression at a specific position with a different expression or
  sequence of expressions.
- With `replace_many()` works the same as `replace()`, but you can replace multiple positions at once.
- With `replace_all()` you can apply a set of replacement rules repeatedly to an expression.
- With `is_match()` you can check whether a pattern matches a subject expression.
"""

from __future__ import absolute_import
import itertools
import math
from typing import Callable, List, NamedTuple, Sequence, Tuple, Union, Iterable

from multiset import Multiset

from .expressions.expressions import (
    Expression, Operation, Pattern, Wildcard, SymbolWildcard, AssociativeOperation, CommutativeOperation
)
from .expressions.substitution import Substitution
from .expressions.functions import preorder_iter_with_position, create_operation_expression
#from .matching.one_to_one import match

__all__ = [u'substitute', u'replace', u'replace_all', u'replace_many', u'is_match', u'ReplacementRule']

Replacement = Union[Expression, List[Expression]]


def substitute(expression, substitution):
    u"""Replaces variables in the given *expression* using the given *substitution*.

    >>> print(substitute(f(x_), {'x': a}))
    f(a)

    If nothing was substituted, the original expression is returned:

    >>> expression = f(x_)
    >>> result = substitute(expression, {'y': a})
    >>> print(result)
    f(x_)
    >>> expression is result
    True

    Note that this function returns a list of expressions iff the expression is a variable and its substitution
    is a list of expressions. In other cases were a substitution is a list of expressions, the expressions will
    be integrated as operands in the surrounding operation:

    >>> print(substitute(f(x_, c), {'x': [a, b]}))
    f(a, b, c)

    If you substitute with a `Multiset` of values, they will be sorted:

    >>> replacement = Multiset([b, a, b])
    >>> print(substitute(f(x_, c), {'x': replacement}))
    f(a, b, b, c)

    Parameters:
        expression:
            An expression in which variables are substituted.
        substitution:
            A substitution dictionary. The key is the name of the variable,
            the value either an expression or a list of expression to use as a replacement for
            the variable.

    Returns:
        The expression resulting from applying the substitution.
    """
    if isinstance(expression, Pattern):
        expression = expression.expression
    return _substitute(expression, substitution)[0]


def _substitute(expression, substitution):
    if getattr(expression, u'variable_name', False) and expression.variable_name in substitution:
        return substitution[expression.variable_name], True
    elif isinstance(expression, Operation):
        any_replaced = False
        new_operands = []
        for operand in expression:
            result, replaced = _substitute(operand, substitution)
            if replaced:
                any_replaced = True
            if isinstance(result, Expression):
                new_operands.append(result)
            elif isinstance(result, Multiset):
                new_operands.extend(sorted(result))
            else:
                new_operands.extend(result)
        if any_replaced:
            return create_operation_expression(expression, new_operands), True

    return expression, False


def replace(expression, position, replacement):
    u"""Replaces the subexpression of `expression` at the given `position` with the given `replacement`.

    The original `expression` itself is not modified, but a modified copy is returned. If the replacement
    is a list of expressions, it will be expanded into the list of operands of the respective operation:

    >>> print(replace(f(a), (0, ), [b, c]))
    f(b, c)

    Parameters:
        expression:
            An :class:`Expression` where a (sub)expression is to be replaced.
        position:
            A tuple of indices, e.g. the empty tuple refers to the `expression` itself,
            `(0, )` refers to the first child (operand) of the `expression`, `(0, 0)` to the first
            child of the first child etc.
        replacement:
            Either an :class:`Expression` or a list of :class:`Expression`\s to be
            inserted into the `expression` instead of the original expression at that `position`.

    Returns:
        The resulting expression from the replacement.

    Raises:
        IndexError: If the position is invalid or out of range.
    """

    if len(position) == 0:
        return replacement
    if not isinstance(expression, Operation):
        raise IndexError(u"Invalid position {!r} for expression {!s}".format(position, expression))
    if position[0] >= len(expression):
        raise IndexError(u"Position {!r} out of range for expression {!s}".format(position, expression))
    pos = position[0]
    operands = list(expression)
    subexpr = replace(operands[pos], position[1:], replacement)
    if isinstance(subexpr, Sequence):
        new_operands = tuple(operands[:pos]) + tuple(subexpr) + tuple(operands[pos + 1:])
        return create_operation_expression(expression, new_operands)
    operands[pos] = subexpr
    return create_operation_expression(expression, operands)


def replace_many(expression, replacements):
    u"""Replaces the subexpressions of *expression* at the given positions with the given replacements.

    The original *expression* itself is not modified, but a modified copy is returned. If the replacement
    is a sequence of expressions, it will be expanded into the list of operands of the respective operation.

    This function works the same as `replace`, but allows multiple positions to be replaced at the same time.
    However, compared to just replacing each position individually with `replace`, this does work when positions are
    modified due to replacing a position with a sequence:

    >>> expr = f(a, b)
    >>> expected_result = replace_many(expr, [((0, ), [c, c]), ((1, ), a)])
    >>> print(expected_result)
    f(c, c, a)

    However, using `replace` for one position at a time gives the wrong result:

    >>> step1 = replace(expr, (0, ), [c, c])
    >>> print(step1)
    f(c, c, b)
    >>> step2 = replace(step1, (1, ), a)
    >>> print(step2)
    f(c, a, b)

    Parameters:
        expression:
            An :class:`Expression` where a (sub)expression is to be replaced.
        replacements:
            A collection of tuples consisting of a position in the expression and a replacement for that position.
            With just a single replacement pair, this is equivalent to using `replace`:

            >>> replace(a, (), b) == replace_many(a, [((), b)])
            True

    Returns:
        The resulting expression from the replacements.

    Raises:
        IndexError: If a position is invalid or out of range or if you try to replace a subterm of a term you are
        already replacing.
    """

    if len(replacements) == 0:
        return expression
    replacements = sorted(replacements)
    if len(replacements[0][0]) == 0:
        if len(replacements) > 1:
            raise IndexError(
                u"Cannot replace child positions for expression {}, got {!r}".format(expression, replacements[1:])
            )
        return replacements[0][1]
    if len(replacements) == 1:
        return replace(expression, replacements[0][0], replacements[0][1])
    if not isinstance(expression, Operation):
        raise IndexError(u"Invalid replacements {!r} for expression {!s}".format(replacements, expression))
    operands = list(expression)
    new_operands = []
    last_index = 0
    for index, group in itertools.groupby(replacements, lambda r: r[0][0]):
        new_operands.extend(operands[last_index:index])
        replacements = [(pos[1:], r) for pos, r in group]
        if len(replacements) == 1:
            replacement = replace(operands[index], replacements[0][0], replacements[0][1])
        else:
            replacement = replace_many(operands[index], replacements)
        if isinstance(replacement, Expression):
            new_operands.append(replacement)
        else:
            new_operands.extend(replacement)
        last_index = index + 1
    new_operands.extend(operands[last_index:len(operands)])
    return create_operation_expression(expression, new_operands)


ReplacementRule = NamedTuple(u'ReplacementRule', [(u'pattern', Pattern), (u'replacement', Callable[..., Expression])])


def replace_all(expression, rules, max_count=float('inf')):
    u"""Replace all occurrences of the patterns according to the replacement rules.

    A replacement rule consists of a *pattern*, that is matched against any subexpression
    of the expression. If a match is found, the *replacement* callback of the rule is called with
    the variables from the match substitution. Whatever the callback returns is used as a replacement for the
    matched subexpression. This can either be a single expression or a sequence of expressions, which is then
    integrated into the surrounding operation in place of the subexpression.

    Note that the pattern can therefore not be a single sequence variable/wildcard, because only single expressions
    will be matched.

    Args:
        expression:
            The expression to which the replacement rules are applied.
        rules:
            A collection of replacement rules that are applied to the expression.
        max_count:
            If given, at most *max_count* applications of the rules are performed. Otherwise, the rules
            are applied until there is no more match. If the set of replacement rules is not confluent,
            the replacement might not terminate without a *max_count* set.

    Returns:
        The resulting expression after the application of the replacement rules. This can also be a sequence of
        expressions, if the root expression is replaced with a sequence of expressions by a rule.
    """
    # TODO Fix the use of head, does not work for head = None
    rules = [ReplacementRule(pattern, replacement) for pattern, replacement in rules]
    expression = expression
    replaced = True
    replace_count = 0
    while replaced and replace_count < max_count:
        replaced = False
        for subexpr, pos in preorder_iter_with_position(expression):
            for pattern, replacement in rules:
                try:
                    subst = next(match(subexpr, pattern))
                    result = replacement(**subst)
                    expression = replace(expression, pos, result)
                    replaced = True
                    break
                except StopIteration:
                    pass
            if replaced:
                break
        replace_count += 1

    return expression


def is_match(subject, pattern):
    u"""
    Check whether the given *subject* matches given *pattern*.

    Args:
        subject:
            The subject.
        pattern:
            The pattern.

    Returns:
        True iff the subject matches the pattern.
    """
    return any(True for _ in match(subject, pattern))
