from __future__ import absolute_import
from typing import Dict, Tuple, Union, cast

from multiset import Multiset

from . import expressions

__all__ = [u'Substitution']

VariableReplacement = Union[Tuple[u'expressions.Expression', ...], Multiset, u'expressions.Expression']


class Substitution(Dict[str, VariableReplacement]):
    """Special :class:`dict` for substitutions with nicer formatting.

    The key is a variable's name and the value the replacement for it.
    """

    def try_add_variable(self, variable_name, replacement):
        """Try to add the variable with its replacement to the substitution.

        This considers an existing replacement and will only succeed if the new replacement
        can be merged with the old replacement. Merging can occur if either the two replacements
        are equivalent. Replacements can also be merged if the old replacement for the variable_name was
        unordered (i.e. a :class:`~.Multiset`) and the new one is an equivalent ordered version of it:

        Args:
            variable:
                The name of the variable to add.
            replacement:
                The replacement for the variable.

        Raises:
            ValueError:
                if the variable cannot be merged because it conflicts with the existing
                substitution for the variable_name.
        """
        if variable_name not in self:
            self[variable_name] = replacement.copy() if isinstance(replacement, Multiset) else replacement
        else:
            existing_value = self[variable_name]

            if isinstance(existing_value, tuple):
                if isinstance(replacement, Multiset):
                    if Multiset(existing_value) != replacement:
                        raise ValueError
                elif replacement != existing_value:
                    raise ValueError
            elif isinstance(existing_value, Multiset):
                if isinstance(replacement, expressions.Expression):
                    raise ValueError
                compare_value = Multiset(replacement)
                if existing_value == compare_value:
                    if not isinstance(replacement, Multiset):
                        self[variable_name] = replacement
                else:
                    raise ValueError
            elif replacement != existing_value:
                raise ValueError

    def union_with_variable(self, variable, replacement):
        """Try to create a new substitution with the given variable added.

        See :meth:`try_add_variable` for a version of this method that modifies the substitution
        in place.

        Args:
            variable_name:
                The name of the variable to add.
            replacement:
                The substitution for the variable.

        Returns:
            The new substitution with the variable_name added or merged.

        Raises:
            ValueError:
                if the variable cannot be merged because it conflicts with the existing
                substitution for the variable.
        """
        new_subst = Substitution(self)
        new_subst.try_add_variable(variable, replacement)
        return new_subst

    def extract_substitution(self, subject, pattern):
        """Extract the variable substitution for the given pattern and subject.

        This assumes that subject and pattern already match when being considered as linear.
        Also, they both must be :term:`syntactic`, as sequence variables cannot be handled here.
        All that this method does is checking whether all the substitutions for the variables can be unified.
        So, in case it returns ``False``, the substitution is invalid for the match.

        ..warning::

            This method mutates the substitution and will even do so in case the extraction fails.

            Create a copy before using this method if you need to preserve the original substitution.

        Args:
            subject:
                A :term:`syntactic` subject that matches the pattern.
            pattern:
                A :term:`syntactic` pattern that matches the subject.

        Returns:
            ``True`` iff the substitution could be extracted successfully.
        """
        if getattr(pattern, u'variable_name', False):
            try:
                self.try_add_variable(pattern.variable_name, subject)
            except ValueError:
                return False
            return True
        elif isinstance(pattern, expressions.Operation):
            assert isinstance(subject, type(pattern))
            assert len(subject) == len(pattern)
            op_expression = cast(expressions.Operation, subject)
            for subj, patt in zip(op_expression, pattern):
                if not self.extract_substitution(subj, patt):
                    return False
        return True

    def union(self, *others):
        """Try to merge the substitutions.

        If a variable occurs in multiple substitutions, try to merge the replacements.
        See :meth:`union_with_variable` to see how replacements are merged.

        Does not modify any of the original substitutions.

        Args:
            others:
                The other substitutions to merge with this one.

        Returns:
            The new substitution with the other substitutions merged.

        Raises:
            ValueError:
                if a variable occurs in multiple substitutions but cannot be merged because the
                substitutions conflict.
        """
        new_subst = Substitution(self)
        for other in others:
            for variable_name, replacement in other.items():
                new_subst.try_add_variable(variable_name, replacement)
        return new_subst

    def rename(self, renaming):
        """Return a copy of the substitution with renamed variables.

        Args:
            renaming:
                A dictionary mapping old variable names to new ones.

        Returns:
            A copy of the substitution where variable names have been replaced according to the given renaming
            dictionary. Names that are not contained in the dictionary are left unchanged.
        """
        return Substitution((renaming.get(name, name), value) for name, value in self.items())

    @staticmethod
    def _match_value_repr_str(value
                             ):  # pragma: no cover
        if isinstance(value, (list, tuple)):
            return u'({!s})'.format(u', '.join(str(x) for x in value))
        if isinstance(value, Multiset):
            return u'{{{!s}}}'.format(u', '.join(str(x) for x in sorted(value)))
        return str(value)

    def __str__(self):
        return u'{{{}}}'.format(
            u', '.join(u'{!s} ↦ {!s}'.format(k, self._match_value_repr_str(v)) for k, v in sorted(self.items()))
        )

    def __repr__(self):
        return u'{{{}}}'.format(u', '.join(u'{!r}: {!r}'.format(k, v) for k, v in sorted(self.items())))

    def __copy__(self):
        return type(self)(self)
