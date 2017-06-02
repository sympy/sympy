u"""Contains several pattern constraint classes.

A pattern constraint is used to further filter which subjects a pattern matches.

The most common use would be the :class:`CustomConstraint`, which wraps a lambda or function to act as a constraint:

>>> a_symbol_constraint = CustomConstraint(lambda x: x.name.startswith('a'))
>>> pattern = Pattern(x_, a_symbol_constraint)
>>> is_match(Symbol('a1'), pattern)
True
>>> is_match(Symbol('b1'), pattern)
False

There is also the :class:`EqualVariablesConstraint` which will try to unify the substitutions of the variables and only
match if it succeeds:

>>> equal_constraint = EqualVariablesConstraint('x', 'y')
>>> pattern = Pattern(f(x_, y_), equal_constraint)
>>> is_match(f(a, a), pattern)
True
>>> is_match(f(a, b), pattern)
False

You can also create a subclass of the :class:`Constraint` class to create your own custom constraint type.
"""

from __future__ import absolute_import

import funcsigs

from collections import OrderedDict

from . import substitution
from ..utils import get_short_lambda_source

__all__ = [u'Constraint', u'EqualVariablesConstraint', u'CustomConstraint']


class Constraint(object):  # pylint: disable=too-few-public-methods
    u"""Base for pattern constraints.

    A constraint is essentially a callback, that receives the match :class:`Substitution` and returns a :class:`bool`
    indicating whether the match is valid.

    You have to override all the abstract methods if you wish to create your own subclass.
    """

    def __call__(self, match):  # pylint: disable=missing-raises-doc
        u"""Return True, iff the constraint is fulfilled by the substitution.

        Override this in your subclass to define the actual constraint behavior.

        Args:
            match:
                The (current) match substitution. Note that the matching is done from left to right, so not all
                variables may have a value yet. You need to override `variables` so that the constraint gets
                called once all the variables it depends on have a value assigned to them.

        Returns:
            True, iff the constraint is fulfilled by the substitution.
        """
        raise NotImplementedError

    def __eq__(self, other):
        u"""Constraints need to be equatable."""
        raise NotImplementedError

    def __hash__(self):
        u"""Constraints need to be hashable."""
        raise NotImplementedError

    @property
    def variables(self):
        u"""The names of the variables the constraint depends upon.

        Used by matchers to decide when a constraint can be evaluated (which is when all
        the dependency variables have been assigned a value). If the set is empty, the constraint will
        only be evaluated once the whole match is complete.
        """
        return frozenset()

    def with_renamed_vars(self, renaming):  # pylint: disable=missing-raises-doc
        u"""Return a *copy* of the constraint with renamed variables.
        This is called when the variables in the expression are renamed and hence the ones in the constraint have to be
        renamed as well. A later invocation of :meth:`__call__` will have the new variable names.
        You will have to implement this if your constraint needs to use the variables of the match substitution.
        Note that this can be called multiple times and you might have to account for that.
        Also, this should not modify the original constraint but rather return a copy.
        Args:
            renaming:
                A dictionary mapping old names to new names.
        Returns:
            A copy of the constraint with renamed variables.
        """
        raise NotImplementedError


class EqualVariablesConstraint(Constraint):  # pylint: disable=too-few-public-methods
    u"""A constraint that ensure multiple variables are equal.

    The constraint tries to unify the substitutions for the variables and is fulfilled iff that succeeds.
    """

    def __init__(self, *variables):
        u"""
        Args:
            *variables: The names of the variables to check for equality.
        """
        self._variables = frozenset(variables)

    @property
    def variables(self):
        return self._variables

    def __call__(self, match):
        subst = substitution.Substitution()
        for name in self._variables:
            try:
                subst.try_add_variable(u'_', match[name])
            except ValueError:
                return False
        return True

    def __str__(self):
        return u'({!s})'.format(u' == '.join(sorted(self._variables)))

    def __repr__(self):
        return u'EqualVariablesConstraint({!s})'.format(u' == '.join(sorted(self._variables)))

    def __eq__(self, other):
        return isinstance(other, EqualVariablesConstraint) and self._variables == other._variables

    def __hash__(self):
        return hash(self._variables)

    def with_renamed_vars(self, renaming):
        return EqualVariablesConstraint(*(renaming.get(v, v) for v in self.variables))


class CustomConstraint(Constraint):  # pylint: disable=too-few-public-methods
    u"""Wrapper for lambdas of functions as constraints.

    The parameter names have to be the same as the the variable names in the expression:

    >>> constraint = CustomConstraint(lambda x, y: x.name < y.name)
    >>> pattern = Pattern(f(x_, y_), constraint)
    >>> is_match(f(a, b), pattern)
    True
    >>> is_match(f(b, a), pattern)
    False

    The ordering of the parameters is not important. You only need to have the parameters needed for the constraint,
    not all variables occurring in the pattern.

    Note, that the matching happens from left left to right, so not all variables may have been assigned a value when
    constraint is called. For constraints over multiple variables you should attach the constraint to the last
    variable occurring in the pattern or a surrounding operation.
    """

    def __init__(self, constraint):
        u"""
        Args:
            constraint:
                The constraint callback.

        Raises:
            ValueError:
                If the callback has positional-only or variable parameters (*args and **kwargs).
        """
        self.constraint = constraint
        signature = funcsigs.signature(constraint)

        self._variables = OrderedDict()

        for param in signature.parameters.values():
            if param.kind == funcsigs.Parameter.POSITIONAL_OR_KEYWORD or param.kind == funcsigs.Parameter.KEYWORD_ONLY:
                self._variables[param.name] = param.name
            elif param.kind == funcsigs.Parameter.VAR_KEYWORD:
                raise ValueError(u"Constraint cannot have variable keyword arguments ({})".format(param.name))
            else:
                raise ValueError(
                    u"Constraint cannot have positional-only or variable positional arguments ({})".format(param.name)
                )

    @property
    def variables(self):
        return frozenset(self._variables.values())

    def __call__(self, match):
        args = dict((name, match[var_name]) for name, var_name in self._variables.items())

        return self.constraint(**args)

    def _get_name(self):
        try:
            return get_short_lambda_source(self.constraint) or self.constraint.__name__
        except Exception:
            return u'UNKNOWN'

    def __str__(self):
        return u'({!s})'.format(self._get_name())

    def __repr__(self):
        return u'CustomConstraint({!s})'.format(self._get_name())

    def __eq__(self, other):
        return (
            isinstance(other, CustomConstraint) and self.constraint == other.constraint and
            self._variables == other._variables
        )

    def __hash__(self):
        return hash(self.constraint)

    def with_renamed_vars(self, renaming):
        cc = CustomConstraint(self.constraint)
        for param_name, old_name in list(cc._variables.items()):
            cc._variables[param_name] = renaming.get(old_name, old_name)
        return cc
