"""Bosonic and fermionic quantum operators.

"""

from warnings import warn

from sympy.core.compatibility import u
from sympy import Add, Mul, Pow, Integer, exp, sqrt, conjugate
from sympy.physics.quantum import Operator, Commutator, AntiCommutator, Dagger
from sympy.physics.quantum import HilbertSpace, FockSpace, Ket, Bra
from sympy.functions.special.tensor_functions import KroneckerDelta


__all__ = [
    'BosonOperator',
    'BosonFockKet',
    'BosonFockBra',
    'BosonCoherentKet',
    'BosonCoherentBra',
    'FermionOperator',
    'FermionFockKet',
    'FermionFockBra',
    'normal_order',
    'normal_ordered_form'
]


class BosonOperator(Operator):
    """A bosonic operator that satisfies [a, Dagger(a)] == 1.

    Parameters
    ==========

    name : str
        A string that labels the bosonic mode.

    annihilation : bool
        A bool that indicates if the bosonic operator is an annihilation (True,
        default value) or creation operator (False)

    Examples
    ========

    >>> from sympy.physics.quantum import Dagger, Commutator, BosonOperator
    >>> a = BosonOperator("a")
    >>> Commutator(a, Dagger(a)).doit()
    1
    """

    @property
    def name(self):
        return self.args[0]

    @property
    def is_annihilation(self):
        return bool(self.args[1])

    @classmethod
    def default_args(self):
        return ("a", True)

    def __new__(cls, *args, **hints):
        if not len(args) in [1, 2]:
            raise ValueError('1 or 2 parameters expected, got %s' % args)

        if len(args) == 1:
            args = (args[0], Integer(1))

        if len(args) == 2:
            args = (args[0], Integer(args[1]))

        return Operator.__new__(cls, *args)

    def _eval_commutator(self, other, **hints):

        if isinstance(other, BosonOperator):
            if self.name == other.name:
                # [a^\dagger, a] = -1
                if not self.is_annihilation and other.is_annihilation:
                    return Integer(-1)

        return None

    def _eval_adjoint(self):
        return BosonOperator(str(self.name), not self.is_annihilation)

    def _print_contents_latex(self, printer, *args):
        if self.is_annihilation:
            return r'{%s}' % str(self.name)
        else:
            return r'{{%s}^\dag}' % str(self.name)

    def _print_contents(self, printer, *args):
        if self.is_annihilation:
            return r'%s' % str(self.name)
        else:
            return r'Dagger(%s)' % str(self.name)

    def _print_contents_pretty(self, printer, *args):
        from sympy.printing.pretty.stringpict import prettyForm
        pform = printer._print(self.args[0], *args)
        if self.is_annihilation:
            return pform**prettyForm(u('\u2020'))
        else:
            return pform


class FermionOperator(Operator):
    """A fermionic operator that satisfies {c, Dagger(c)} == 1.

    Parameters
    ==========

    name : str
        A string that labels the fermionic mode.

    annihilation : bool
        A bool that indicates if the fermionic operator is an annihilation
        (True, default value) or creation operator (False)

    Examples
    ========

    >>> from sympy.physics.quantum import Dagger, AntiCommutator
    >>> from sympy.physics.quantum import FermionOperator
    >>> c = FermionOperator("c")
    >>> AntiCommutator(c, Dagger(c)).doit()
    1
    """
    @property
    def name(self):
        return self.args[0]

    @property
    def is_annihilation(self):
        return bool(self.args[1])

    @classmethod
    def default_args(self):
        return ("c", True)

    def __new__(cls, *args, **hints):
        if not len(args) in [1, 2]:
            raise ValueError('1 or 2 parameters expected, got %s' % args)

        if len(args) == 1:
            args = (args[0], Integer(1))

        if len(args) == 2:
            args = (args[0], Integer(args[1]))

        return Operator.__new__(cls, *args)

    def _eval_anticommutator(self, other, **hints):

        if isinstance(other, FermionOperator):
            if self.name == other.name:
                # {a^\dagger, a} = 1
                if not self.is_annihilation and other.is_annihilation:
                    return Integer(1)

        return None

    def _eval_adjoint(self):
        return FermionOperator(str(self.name), not self.is_annihilation)

    def _print_contents_latex(self, printer, *args):
        if self.is_annihilation:
            return r'{%s}' % str(self.name)
        else:
            return r'{{%s}^\dag}' % str(self.name)

    def _print_contents(self, printer, *args):
        if self.is_annihilation:
            return r'%s' % str(self.name)
        else:
            return r'Dagger(%s)' % str(self.name)

    def _print_contents_pretty(self, printer, *args):
        from sympy.printing.pretty.stringpict import prettyForm
        pform = printer._print(self.args[0], *args)
        if self.is_annihilation:
            return pform**prettyForm(u('\u2020'))
        else:
            return pform


def _expand_powers(factors):
    """
    Helper function for normal_ordered_form and normal_order: Expand a
    power expression to a multiplication expression so that that the
    expression can be handled by the normal ordering functions.
    """

    new_factors = []
    for factor in factors.args:
        if isinstance(factor, Pow):
            for n in range(factor.args[1]):
                new_factors.append(factor.args[0])
        else:
            new_factors.append(factor)

    return new_factors


def _normal_ordered_form_factor(product, independent=False, recursive_limit=10,
                                _recursive_depth=0):
    """
    Helper function for normal_ordered_form_factor: Write multiplication
    expression with bosonic or fermionic operators on normally ordered form,
    using the bosonic and fermionic commutation relations. The resulting
    operator expression is equivalent to the argument, but will in general be
    a sum of operator products instead of a simple product.
    """

    factors = _expand_powers(product)

    new_factors = []
    n = 0
    while n < len(factors) - 1:

        if isinstance(factors[n], BosonOperator) and factors[n].is_annihilation:
            # boson
            if not isinstance(factors[n + 1], BosonOperator):
                new_factors.append(factors[n])
            else:
                if factors[n + 1].is_annihilation:
                    new_factors.append(factors[n])
                else:
                    if factors[n].args[0] != factors[n + 1].args[0]:
                        if independent:
                            c = 0
                        else:
                            c = Commutator(factors[n], factors[n + 1])
                        new_factors.append(factors[n + 1] * factors[n] + c)
                    else:
                        c = Commutator(factors[n], factors[n + 1])
                        new_factors.append(
                            factors[n + 1] * factors[n] + c.doit())
                    n += 1

        elif (isinstance(factors[n], FermionOperator) and
              factors[n].is_annihilation):
            # fermion
            if not isinstance(factors[n + 1], FermionOperator):
                new_factors.append(factors[n])
            else:
                if factors[n + 1].is_annihilation:
                    new_factors.append(factors[n])
                else:
                    if factors[n].args[0] != factors[n + 1].args[0]:
                        if independent:
                            c = 0
                        else:
                            c = AntiCommutator(factors[n], factors[n + 1])
                        new_factors.append(-factors[n + 1] * factors[n] + c)
                    else:
                        c = AntiCommutator(factors[n], factors[n + 1])
                        new_factors.append(
                            -factors[n + 1] * factors[n] + c.doit())
                    n += 1

        else:
            new_factors.append(factors[n])

        n += 1

    if n == len(factors) - 1:
        new_factors.append(factors[-1])

    if new_factors == factors:
        return product
    else:
        expr = Mul(*new_factors).expand()
        return normal_ordered_form(expr,
                                   recursive_limit=recursive_limit,
                                   _recursive_depth=_recursive_depth + 1,
                                   independent=independent)


def _normal_ordered_form_terms(expr, independent=False, recursive_limit=10,
                               _recursive_depth=0):
    """
    Helper function for normal_ordered_form: loop through each term in an
    addition expression and call _normal_ordered_form_factor to perform the
    factor to an normally ordered expression.
    """

    new_terms = []
    for term in expr.args:
        if isinstance(term, Mul):
            new_term = _normal_ordered_form_factor(
                term, recursive_limit=recursive_limit,
                _recursive_depth=_recursive_depth, independent=independent)
            new_terms.append(new_term)
        else:
            new_terms.append(term)

    return Add(*new_terms)


def normal_ordered_form(expr, independent=False, recursive_limit=10,
                        _recursive_depth=0):
    """Write an expression with bosonic or fermionic operators on normal
    ordered form, where each term is normally ordered. Note that this
    normal ordered form is equivalent to the original expression.

    Parameters
    ==========

    expr : expression
        The expression write on normal ordered form.

    recursive_limit : int (default 10)
        The number of allowed recursive applications of the function.

    Examples
    ========

    >>> from sympy.physics.quantum import normal_ordered_form
    >>> from sympy.physics.quantum import BosonOperator, Dagger
    >>> a = BosonOperator("a")
    >>> normal_ordered_form(a * Dagger(a))
    1 + Dagger(a)*a
    """

    if _recursive_depth > recursive_limit:
        warn.warning("Warning: too many recursions, aborting")
        return expr

    if isinstance(expr, Add):
        return _normal_ordered_form_terms(expr,
                                          recursive_limit=recursive_limit,
                                          _recursive_depth=_recursive_depth,
                                          independent=independent)
    elif isinstance(expr, Mul):
        return _normal_ordered_form_factor(expr,
                                           recursive_limit=recursive_limit,
                                           _recursive_depth=_recursive_depth,
                                           independent=independent)
    else:
        return expr


def _normal_order_factor(product, recursive_limit=10, _recursive_depth=0):
    """
    Helper function for normal_order: Normal order a multiplication expression
    with bosonic or fermionic operators. In general the resulting operator
    expression will not be equivalent to original product.
    """

    factors = _expand_powers(product)

    n = 0
    new_factors = []
    while n < len(factors) - 1:

        if isinstance(factors[n], BosonOperator) and factors[n].is_annihilation:
            # boson
            if not isinstance(factors[n + 1], BosonOperator):
                new_factors.append(factors[n])
            else:
                if factors[n + 1].is_annihilation:
                    new_factors.append(factors[n])
                else:
                    if factors[n].args[0] != factors[n + 1].args[0]:
                        new_factors.append(factors[n + 1] * factors[n])
                    else:
                        new_factors.append(factors[n + 1] * factors[n])
                    n += 1

        elif (isinstance(factors[n], FermionOperator) and
              factors[n].is_annihilation):
            # fermion
            if not isinstance(factors[n + 1], FermionOperator):
                new_factors.append(factors[n])
            else:
                if factors[n + 1].is_annihilation:
                    new_factors.append(factors[n])
                else:
                    if factors[n].args[0] != factors[n + 1].args[0]:
                        new_factors.append(-factors[n + 1] * factors[n])
                    else:
                        new_factors.append(-factors[n + 1] * factors[n])
                    n += 1

        else:
            new_factors.append(factors[n])

        n += 1

    if n == len(factors) - 1:
        new_factors.append(factors[-1])

    if new_factors == factors:
        return product
    else:
        expr = Mul(*new_factors).expand()
        return normal_order(expr,
                            recursive_limit=recursive_limit,
                            _recursive_depth=_recursive_depth + 1)


def _normal_order_terms(expr, recursive_limit=10, _recursive_depth=0):
    """
    Helper function for normal_order: look through each term in an addition
    expression and call _normal_order_factor to perform the normal ordering
    on the factors.
    """

    new_terms = []
    for term in expr.args:
        if isinstance(term, Mul):
            new_term = _normal_order_factor(term,
                                            recursive_limit=recursive_limit,
                                            _recursive_depth=_recursive_depth)
            new_terms.append(new_term)
        else:
            new_terms.append(term)

    return Add(*new_terms)


def normal_order(expr, recursive_limit=10, _recursive_depth=0):
    """Normal order an expression with bosonic or fermionic operators. Note
    that this normal order is not equivalent to the original expression, but
    the creation and annihilation operators in each term in expr is reordered
    so that the expression becomes normal ordered.

    Parameters
    ==========

    expr : expression
        The expression to normal order.

    recursive_limit : int (default 10)
        The number of allowed recursive applications of the function.

    Examples
    ========

    >>> from sympy.physics.quantum import normal_order, BosonOperator, Dagger
    >>> a = BosonOperator("a")
    >>> normal_order(a * Dagger(a))
    Dagger(a)*a
    """
    if _recursive_depth > recursive_limit:
        warn.warning("Warning: too many recursions, aborting")
        return expr

    if isinstance(expr, Add):
        return _normal_order_terms(expr, recursive_limit=recursive_limit,
                                   _recursive_depth=_recursive_depth)
    elif isinstance(expr, Mul):
        return _normal_order_factor(expr, recursive_limit=recursive_limit,
                                    _recursive_depth=_recursive_depth)
    else:
        return expr


class BosonFockKet(Ket):
    """Fock state ket for a bosonic mode.

    Parameters
    ==========

    n : Number
        The Fock state number.

    """

    def __new__(cls, n):
        return Ket.__new__(cls, n)

    @property
    def n(self):
        return self.label[0]

    @classmethod
    def dual_class(self):
        return BosonFockBra

    @classmethod
    def _eval_hilbert_space(cls, label):
        return FockSpace()

    def _eval_innerproduct_BosonFockBra(self, bra, **hints):
        return KroneckerDelta(self.n, bra.n)

    def _apply_operator_BosonOperator(self, op, **options):
        if op.is_annihilation:
            return sqrt(self.n) * BosonFockKet(self.n - 1)
        else:
            return sqrt(self.n + 1) * BosonFockKet(self.n + 1)


class BosonFockBra(Bra):
    """Fock state bra for a bosonic mode.

    Parameters
    ==========

    n : Number
        The Fock state number.

    """

    def __new__(cls, n):
        return Bra.__new__(cls, n)

    @property
    def n(self):
        return self.label[0]

    @classmethod
    def dual_class(self):
        return BosonFockKet

    @classmethod
    def _eval_hilbert_space(cls, label):
        return FockSpace()


class BosonCoherentKet(Ket):
    """Coherent state ket for a bosonic mode.

    Parameters
    ==========

    alpha : Number, Symbol
        The complex amplitude of the coherent state.

    """

    def __new__(cls, alpha):
        return Ket.__new__(cls, alpha)

    @property
    def alpha(self):
        return self.label[0]

    @classmethod
    def dual_class(self):
        return BosonCoherentBra

    @classmethod
    def _eval_hilbert_space(cls, label):
        return HilbertSpace()

    def _eval_innerproduct_BosonCoherentBra(self, bra, **hints):
        if self.alpha == bra.alpha:
            return Integer(1)
        else:
            return exp(-(abs(self.alpha)**2 + abs(bra.alpha)**2 - 2 * conjugate(bra.alpha) * self.alpha)/2)

    def _apply_operator_BosonOperator(self, op, **options):
        if op.is_annihilation:
            return self.alpha * self
        else:
            return None


class BosonCoherentBra(Bra):
    """Coherent state bra for a bosonic mode.

    Parameters
    ==========

    alpha : Number, Symbol
        The complex amplitude of the coherent state.

    """

    def __new__(cls, alpha):
        return Bra.__new__(cls, alpha)

    @property
    def alpha(self):
        return self.label[0]

    @classmethod
    def dual_class(self):
        return BosonCoherentKet

    def _apply_operator_BosonOperator(self, op, **options):
        if not op.is_annihilation:
            return self.alpha * self
        else:
            return None


class FermionFockKet(Ket):
    """Fock state ket for a fermionic mode.

    Parameters
    ==========

    n : Number
        The Fock state number.

    """

    def __new__(cls, n):
        if not n in [0, 1]:
            raise ValueError("n must be 0 or 1")
        return Ket.__new__(cls, n)

    @property
    def n(self):
        return self.label[0]

    @classmethod
    def dual_class(self):
        return FermionFockBra

    @classmethod
    def _eval_hilbert_space(cls, label):
        return HilbertSpace()

    def _eval_innerproduct_FermionFockBra(self, bra, **hints):
        return KroneckerDelta(self.n, bra.n)

    def _apply_operator_FermionOperator(self, op, **options):
        if op.is_annihilation:
            if self.n == 1:
                return FermionFockKet(0)
            else:
                return Integer(0)
        else:
            if self.n == 0:
                return FermionFockKet(1)
            else:
                return Integer(0)


class FermionFockBra(Bra):
    """Fock state bra for a fermionic mode.

    Parameters
    ==========

    n : Number
        The Fock state number.

    """

    def __new__(cls, n):
        if not n in [0, 1]:
            raise ValueError("n must be 0 or 1")
        return Bra.__new__(cls, n)

    @property
    def n(self):
        return self.label[0]

    @classmethod
    def dual_class(self):
        return FermionFockKet
