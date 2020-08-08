from sympy.core import Tuple
from sympy.core.sympify import _sympify
from sympy.sets import Union
from sympy.map import (
    Map, AppliedMap,
    VectorAdditionOperator, Addition,
    ScalarMultiplicationOperator, VectorMultiplicationOperator, Multiplication,
    ExponentElement,
)
from ..module import VectorSpace

__all__ = [
    'FunctionAdditionOperator', 'FunctionAddition',
    'FunctionScalarMultiplicationOperator',
    'FunctionMultiplication',
]

class FunctionAdditionOperator(VectorAdditionOperator):
    r"""
    Class for addition operator between functions whose codomains are vector space [1].

    Explanation
    ===========

    Let $V$ be a vector space defined over a field $F$, and all elements of a set $S$
    are functions whose codomains are $V$. Then, these functions form a vector space
    $V^{*}$ over $F$. In $V^{*}$, elements of $S$ are vectors and element of $F$ are
    scalars.
    To be called, function addition operator needs keyword arguments *sv_mul*. 
    This is to perform the collection of vectors, i.e. $f + f = \left(1+1 \right) f$.

    Parameters
    ==========

    domain : ProductSet of FunctionSet

    codomain : FunctionSet
        This FunctionSet's ``codomain`` attribute must return a vector space.

    identity : Map
        Identity element

    Examples
    ========

    >>> from sympy import (
    ... S, Set, Map, VectorAdditionOperator, AbelianGroup,
    ... ScalarMultiplicationOperator, VectorSpace,
    ... FunctionSet, ConstantMap, FunctionAdditionOperator,
    ... FunctionScalarMultiplicationOperator
    ... )

    Construct underlying structures

    >>> F = S.RealsField
    >>> X = Set('X')
    >>> f, g = Map('f', domain=X, codomain=F), Map('g', domain=X, codomain=F)

    >>> A = Set('A')
    >>> e = A.element('e')
    >>> vadd = VectorAdditionOperator(A**2, A, e)
    >>> G = AbelianGroup('G', (A,), (vadd,))
    >>> smul = ScalarMultiplicationOperator(F*G, G)
    >>> V = VectorSpace('V', (F, G), (smul,))

    >>> fs = FunctionSet(domain=X, codomain=V)
    >>> zerofunc = ConstantMap(F.add_op.identity, domain=X)

    >>> fadd = FunctionAdditionOperator(fs**2, fs, zerofunc)
    >>> fG = AbelianGroup('fG', (fs,), (fadd,))
    >>> fsmul = FunctionScalarMultiplicationOperator(F*fG, fG)

    Function addition

    >>> fadd(f, g, sv_mul=fsmul)
    f + g : X -> R
    >>> fadd(f, f, sv_mul=fsmul, evaluate=True)
    2*f : X -> R

    Constructing a functions pace makes applying the operator easier

    >>> from sympy import VectorSpace
    >>> FS = VectorSpace('FS', (F, fG), (fsmul,))

    >>> FS.add(f, g)
    f + g : X -> R
    >>> FS.add(f, f, evaluate=True)
    2*f : X -> R

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Function_space

    """

    @property
    def ss_add(self):
        return self.codomain.codomain.ring.add_op

    @property
    def ss_mul(self):
        return self.codomain.codomain.ring.mul_op

    def __call__(self, *args, sv_mul, evaluate=False, **kwargs):
        return FunctionAddition(
            self, args, (sv_mul,), evaluate=evaluate, **kwargs
        )

    def apply(self, *args, aux, **kwargs):
        sv_mul, = aux
        ss_add, ss_mul = self.ss_add, self.ss_mul
        evaluate = kwargs.get('evaluate', False)

        # sympify the arguments
        args = [_sympify(a) for a in args]

        if evaluate:
            args = self.addition_process(args, sv_mul, ss_add, ss_mul)

            if not args:
                return self.identity
            elif len(args) == 1:
                return args[0]

        # return Addition class with processed arguments
        args = Tuple(*[_sympify(a) for a in args])
        aux = Tuple(*[_sympify(a) for a in aux])
        result = super(AppliedMap, FunctionAddition).__new__(
            FunctionAddition, self, args, aux,
        )
        return result

class FunctionAddition(Addition, Map):
    """
    Added function

    Parameters
    ==========

    mapping : Map

    args : tuple of arguments
        Arguments applied to *mapping*

    aux : tuple of Maps
        Auxillary operators

    evaluate : bool, optional
        If True, returns evaluated application of *map* to *args*

    Examples
    ========

    >>> from sympy import (
    ... S, Set, Map, VectorAdditionOperator, AbelianGroup,
    ... ScalarMultiplicationOperator, VectorSpace,
    ... FunctionSet, ConstantMap, FunctionAdditionOperator,
    ... FunctionScalarMultiplicationOperator, VectorSpace
    ... )

    Construct underlying structures

    >>> F = S.RealsField
    >>> X = Set('X')
    >>> f, g = Map('f', domain=X, codomain=F), Map('g', domain=X, codomain=F)

    >>> A = Set('A')
    >>> e = A.element('e')
    >>> vadd = VectorAdditionOperator(A**2, A, e)
    >>> G = AbelianGroup('G', (A,), (vadd,))
    >>> smul = ScalarMultiplicationOperator(F*G, G)
    >>> V = VectorSpace('V', (F, G), (smul,))

    >>> fs = FunctionSet(domain=X, codomain=V)
    >>> zerofunc = ConstantMap(F.add_op.identity, domain=X)

    >>> fadd = FunctionAdditionOperator(fs**2, fs, zerofunc)
    >>> fG = AbelianGroup('fG', (fs,), (fadd,))
    >>> fsmul = FunctionScalarMultiplicationOperator(F*fG, fG)

    >>> FS = VectorSpace('FS', (F, fG), (fsmul,))
    >>> addf = FS.add(f, g)

    >>> addf
    f + g : X -> R

    >>> x = X.element('x')
    >>> addf(x)
    (f + g)(x)

    >>> addf(x, evaluate=True)
    f(x) + g(x)

    """

    @property
    def domain(self):
        return Union(*[f.domain for f in self.arguments])

    @property
    def codomain(self):
        return Union(*[f.codomain for f in self.arguments])

    def eval(self, *args):
        add = self.map.codomain.codomain.add
        terms = [a(*args) for a in self.arguments]
        return add(*terms)

class FunctionScalarMultiplicationOperator(ScalarMultiplicationOperator):
    r"""
    Class for multiplication operator between scalar function
    whose codomain is vector space [1].

    Explanation
    ===========

    Let $V$ be a vector space defined over a field $F$, and all elements of a set $S$
    are functions whose codomains are $V$. Then, these functions form a vector space
    $V^{*}$ over $F$. In $V^{*}$, elements of $S$ are vectors and element of $F$ are
    scalars.

    Parameters
    ==========

    domain : ProductSet of Ring and AbelianGroup

    codomain : AbelianGroup

    Examples
    ========

    >>> from sympy import (
    ... S, Set, Map, VectorAdditionOperator, AbelianGroup,
    ... ScalarMultiplicationOperator, VectorSpace,
    ... FunctionSet, ConstantMap, FunctionAdditionOperator,
    ... FunctionScalarMultiplicationOperator
    ... )

    Construct underlying structures

    >>> F = S.RealsField
    >>> X = Set('X')
    >>> f = Map('f', domain=X, codomain=F)

    >>> A = Set('A')
    >>> e = A.element('e')
    >>> vadd = VectorAdditionOperator(A**2, A, e)
    >>> G = AbelianGroup('G', (A,), (vadd,))
    >>> smul = ScalarMultiplicationOperator(F*G, G)
    >>> V = VectorSpace('V', (F, G), (smul,))

    >>> fs = FunctionSet(domain=X, codomain=V)
    >>> zerofunc = ConstantMap(F.add_op.identity, domain=X)

    >>> fadd = FunctionAdditionOperator(fs**2, fs, zerofunc)
    >>> fG = AbelianGroup('fG', (fs,), (fadd,))
    >>> fsmul = FunctionScalarMultiplicationOperator(F*fG, fG)

    Multiplication of function with scalar

    >>> fsmul(2, f)
    2*f : X -> R

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Function_space

    """
    def __call__(self, a, b, evaluate=False, **kwargs):
        kwargs.update(evaluate=evaluate)
        return FunctionMultiplication(self, (a, b), (), **kwargs)

    def apply(self, a, b, aux, **kwargs):
        evaluate = kwargs.get('evaluate', False)

        # sympify the arguments
        args = [_sympify(i) for i in (a, b)]

        # order a and b so that a is scalar and b is vector
        a, = [i for i in args if self.scalar_ring.contains(i) == True]
        b, = [i for i in args if self.vector_group.contains(i) == True]

        if evaluate:
            args = self.multiplication_process(args)

            if len(args) == 1:
                return args[0]

        # return Multiplication class with processed arguments
        args = Tuple(*[_sympify(a) for a in args])
        aux = Tuple(*[_sympify(a) for a in aux])
        result = super(AppliedMap, FunctionMultiplication).__new__(
            FunctionMultiplication, self, args, aux
        )
        return result

class FunctionMultiplication(Multiplication, Map):
    """
    Multiplied function

    Parameters
    ==========

    mapping : Map

    args : tuple of arguments
        Arguments applied to *mapping*

    aux : tuple of Maps
        Auxillary operators

    evaluate : bool, optional
        If True, returns evaluated application of *map* to *args*

    Examples
    ========

    >>> from sympy import (
    ... S, Set, Map, VectorAdditionOperator, AbelianGroup,
    ... ScalarMultiplicationOperator, VectorSpace,
    ... FunctionSet, ConstantMap, FunctionAdditionOperator,
    ... FunctionScalarMultiplicationOperator
    ... )

    Construct underlying structures

    >>> F = S.RealsField
    >>> X = Set('X')
    >>> f = Map('f', domain=X, codomain=F)

    >>> A = Set('A')
    >>> e = A.element('e')
    >>> vadd = VectorAdditionOperator(A**2, A, e)
    >>> G = AbelianGroup('G', (A,), (vadd,))
    >>> smul = ScalarMultiplicationOperator(F*G, G)
    >>> V = VectorSpace('V', (F, G), (smul,))

    >>> fs = FunctionSet(domain=X, codomain=V)
    >>> zerofunc = ConstantMap(F.add_op.identity, domain=X)

    >>> fadd = FunctionAdditionOperator(fs**2, fs, zerofunc)
    >>> fG = AbelianGroup('fG', (fs,), (fadd,))
    >>> fsmul = FunctionScalarMultiplicationOperator(F*fG, fG)
    >>> mulf = fsmul(2, f)

    >>> mulf
    2*f : X -> R

    >>> x = X.element('x')
    >>> mulf(x)
    (2*f)(x)

    >>> mulf(x, evaluate=True)
    2*f(x)

    """
    @property
    def domain(self):
        vectors = [a for a in self.arguments if self.map.codomain.contains(a) == True]
        return Union(*[f.domain for f in vectors])

    @property
    def codomain(self):
        vectors = [a for a in self.arguments if self.map.codomain.contains(a) == True]
        return Union(*[f.codomain for f in vectors])

    def eval(self, *args):
        mul = self.map.scalar_ring.mul
        terms = []
        for a in self.arguments:
            if self.map.codomain.contains(a) == True:
                terms.append(a(*args))
            else:
                terms.append(a)
        return mul(*terms)
