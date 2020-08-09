from sympy.assumptions import ask, Q
from sympy.core import Tuple
from sympy.core.sympify import _sympify
from sympy.sets import Union
from sympy.map import (
    Map, AppliedMap,
    VectorAdditionOperator, Addition,
    ScalarMultiplicationOperator, VectorMultiplicationOperator, Multiplication,
    ExponentElement, ConstantMap, ExponentElement, InverseElement,
)

__all__ = [
    'FunctionAdditionOperator', 'FunctionAddition',
    'FunctionScalarMultiplicationOperator',
    'FunctionMultiplication',
    'FunctionVectorMultiplicationOperator',
    'FunctionExponent', 'ReciprocalFunction',
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
    ... S, Map, FunctionSet, ConstantMap,
    ... FunctionAdditionOperator, AbelianGroup,
    ... FunctionScalarMultiplicationOperator
    ... )

    Construct underlying structures.
    Although function space needs vector space as codomain, scalar field
    can be given if no vector is involved.

    >>> F = S.RealsField
    >>> X = S.Reals
    >>> f, g = Map('f', domain=X, codomain=F), Map('g', domain=X, codomain=F)

    >>> fs = FunctionSet(domain=X, codomain=F)
    >>> zerofunc = ConstantMap(F.add_op.identity, domain=X)
    >>> fadd = FunctionAdditionOperator(fs**2, fs, zerofunc)
    >>> fG = AbelianGroup('fG', (fs,), (fadd,))

    >>> fsmul = FunctionScalarMultiplicationOperator(F*fG, fG)

    Function addition

    >>> fadd(f, g, sv_mul=fsmul)
    f + g : Reals -> R
    >>> fadd(f, f, sv_mul=fsmul, evaluate=True)
    2*f : Reals -> R

    Constructing a function space makes applying the operator easier

    >>> from sympy import VectorSpace
    >>> FS = VectorSpace('FS', (F, fG), (fsmul,))

    >>> FS.add(f, g)
    f + g : Reals -> R
    >>> FS.add(f, f, evaluate=True)
    2*f : Reals -> R

    Using infix operator constructs the suitable structure and returns evaluated result.

    >>> f + f
    2*f : Reals -> R

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
    ... FunctionScalarMultiplicationOperator
    ... )

    Construct underlying structures

    >>> F = S.RealsField
    >>> X = Set('X')

    In this example, functions return vector. Therefore, full vector space
    must be constructed first.

    >>> A = Set('A')
    >>> e = A.element('e')
    >>> vadd = VectorAdditionOperator(A**2, A, e)
    >>> G = AbelianGroup('G', (A,), (vadd,))
    >>> smul = ScalarMultiplicationOperator(F*G, G)
    >>> V = VectorSpace('V', (F, G), (smul,)) # vector space which is function's codomain

    >>> f, g = Map('f', domain=X, codomain=G), Map('g', domain=X, codomain=G)

    >>> fs = FunctionSet(domain=X, codomain=V)
    >>> zerofunc = ConstantMap(F.add_op.identity, domain=X)

    >>> fadd = FunctionAdditionOperator(fs**2, fs, zerofunc)
    >>> fG = AbelianGroup('fG', (fs,), (fadd,))
    >>> fsmul = FunctionScalarMultiplicationOperator(F*fG, fG)

    >>> FS = VectorSpace('FS', (F, fG), (fsmul,))
    >>> addf = FS.add(f, g)

    Added function

    >>> addf
    f + g : X -> G

    Argument can be applied to added function

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

    @property
    def codomain_structure(self):
        return self.map.codomain.codomain

    def _contained(self, other):
        return Map._contained(self, other)

    def eval(self, *args):
        add = self.map.codomain.codomain.add
        terms = [a(*args) for a in self.arguments]
        return add(*terms)

class FunctionScalarMultiplicationOperator(ScalarMultiplicationOperator):
    r"""
    Class for multiplication operator between scalar and function
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
    ... S, Map, FunctionSet, ConstantMap,
    ... FunctionAdditionOperator, AbelianGroup,
    ... FunctionScalarMultiplicationOperator
    ... )

    Construct underlying structures.
    Although function space needs vector space as codomain, scalar field
    can be given if no vector is involved.

    >>> F = S.RealsField
    >>> X = S.Reals
    >>> f = Map('f', domain=X, codomain=F)

    >>> fs = FunctionSet(domain=X, codomain=F)
    >>> zerofunc = ConstantMap(F.add_op.identity, domain=X)
    >>> fadd = FunctionAdditionOperator(fs**2, fs, zerofunc)
    >>> fG = AbelianGroup('fG', (fs,), (fadd,))

    >>> fsmul = FunctionScalarMultiplicationOperator(F*fG, fG)

    Multiplication of function with scalar

    >>> fsmul(2, f)
    2*f : Reals -> R

    Using infix operator constructs the suitable structure and returns evaluated result.

    >>> 2*f
    2*f : Reals -> R

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Function_space

    """
    @property
    def codomain_structure(self):
        return self.codomain.domain.codomain

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

    Construct underlying structures.
    In this example, functions return vector. Therefore, full vector space
    must be constructed first.

    >>> F = S.RealsField
    >>> X = Set('X')

    >>> A = Set('A')
    >>> e = A.element('e')
    >>> vadd = VectorAdditionOperator(A**2, A, e)
    >>> G = AbelianGroup('G', (A,), (vadd,))
    >>> smul = ScalarMultiplicationOperator(F*G, G)
    >>> V = VectorSpace('V', (F, G), (smul,)) # vector space which is function's codomain

    >>> f, g = Map('f', domain=X, codomain=G), Map('g', domain=X, codomain=G)

    >>> fs = FunctionSet(domain=X, codomain=V)
    >>> zerofunc = ConstantMap(F.add_op.identity, domain=X)

    >>> fadd = FunctionAdditionOperator(fs**2, fs, zerofunc)
    >>> fG = AbelianGroup('fG', (fs,), (fadd,))
    >>> fsmul = FunctionScalarMultiplicationOperator(F*fG, fG)

    >>> FS = VectorSpace('FS', (F, fG), (fsmul,))
    >>> mulf = FS.mul(2, f)

    Multiplied function

    >>> mulf
    2*f : X -> G

    Argument can be applied to multiplied function

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

    @property
    def codomain_structure(self):
        return self.map.codomain_structure

    def _contained(self, other):
        return Map._contained(self, other)

    def eval(self, *args):
        mul = self.codomain_structure.mul
        terms = []
        for a in self.arguments:
            if self.map.codomain.contains(a) == True:
                terms.append(a(*args))
            else:
                terms.append(a)
        return mul(*terms)

class FunctionVectorMultiplicationOperator(VectorMultiplicationOperator):
    r"""
    Class for multiplication operator between two functions
    whose codomains are algebraic algebra [1].

    Parameters
    ==========

    domain : ProductSet of two VectorSpace

    codomain : VectorSpace

    Examples
    ========

    >>> from sympy import (
    ... S, Map, FunctionSet, ConstantMap,
    ... FunctionAdditionOperator, AbelianGroup,
    ... FunctionScalarMultiplicationOperator,
    ... VectorSpace, FunctionVectorMultiplicationOperator,
    ... )

    Construct underlying structures.
    Although function space needs vector space as codomain, scalar field
    can be given if no vector is involved.

    >>> F = S.RealsField
    >>> X = S.Reals
    >>> f, g = Map('f', domain=X, codomain=F), Map('g', domain=X, codomain=F)

    >>> fs = FunctionSet(domain=X, codomain=F)
    >>> zerofunc = ConstantMap(F.add_op.identity, domain=X)
    >>> fadd = FunctionAdditionOperator(fs**2, fs, zerofunc)
    >>> fG = AbelianGroup('fG', (fs,), (fadd,))

    >>> fsmul = FunctionScalarMultiplicationOperator(F*fG, fG)

    >>> FS = VectorSpace('FS', (F, fG), (fsmul,))
    >>> ff_mul = FunctionVectorMultiplicationOperator(FS*FS, fG)

    >>> ff_mul(f, g)
    f*g : Reals -> R

    >>> x = X.element('x')
    >>> ff_mul(f, g)(x)
    (f*g)(x)
    >>> ff_mul(f, g)(x, evaluate=True)
    f(x)*g(x)

    Using infix operator constructs the suitable structure and returns evaluated result.

    >>> f*g
    f*g : Reals -> R

    >>> f*f
    f**2 : Reals -> R

    """

    @property
    def codomain_structure(self):
        return self.codomain.domain.codomain

    @property
    def identity(self):
        iden_scalar = self.codomain_structure.ring.mul_op.identity
        domain = self.codomain.domain.domain
        return ConstantMap(iden_scalar, domain)

    def __call__(self, *args, evaluate=False, **kwargs):
        kwargs.update(evaluate=evaluate)
        return FunctionMultiplication(self, args, (), **kwargs)

    def apply(self, *args, aux, **kwargs):
        evaluate = kwargs.get('evaluate', False)

        # sympify the arguments
        args = [_sympify(i) for i in args]

        if evaluate:
            args = self.multiplication_process(args)

            if len(args) == 0:
                return self.identity
            if len(args) == 1:
                return args[0]

        # return Multiplication class with processed arguments
        args = Tuple(*[_sympify(a) for a in args])
        aux = Tuple(*[_sympify(a) for a in aux])
        result = super(AppliedMap, FunctionMultiplication).__new__(
            FunctionMultiplication, self, args, aux
        )
        return result

    def is_scalarfunction(self, f):
        """
        Return ``True`` if *f* is a function which returns scalar.
        """
        if f.codomain.is_subset(self.codomain_structure.ring):
            return True
        return False

    def collect_scalarfunc(self, seq):

        scalarfuncs = []
        vectorfuncs = []
        for f in seq:
            if self.is_scalarfunction(f):
                scalarfuncs.append(f)
            else:
                vectorfuncs.append(f)
        combined_sf = self.assoc_comm_process(scalarfuncs)
        return combined_sf + vectorfuncs

    def multiplication_process(self, seq):
        seq = self.remove_identity(seq)
        seq = self.collect_scalarfunc(seq)
        return seq

    def _expop_apply(self, f, n, **kwargs):
        result = super()._expop_apply(f, n, **kwargs)
        if result is None:
            result = super(Map, FunctionExponent).__new__(
                FunctionExponent, self.exponent_operator(), Tuple(f, n)
            )
        return result

    def _invop_apply(self, f, **kwargs):
        if self.is_scalarfunction(f):
            return self.exponent_operator()(f, -1, **kwargs)
        result = super()._invop_apply(f, **kwargs)
        if result is None:
            result = super(Map, ReciprocalFunction).__new__(
                ReciprocalFunction, self.inverse_operator(), Tuple(f)
            )
        return result

class FunctionExponent(ExponentElement, Map):
    """
    Multiplicative exponentiation of function.

    """
    @property
    def domain(self):
        return self.base.domain

    @property
    def codomain(self):
        return self.base.codomain

    @property
    def codomain_structure(self):
        return self.map.base_op.codomain_structure

    def _contained(self, other):
        return Map._contained(self, other)

    def eval(self, *args):
        pow = self.codomain_structure.pow
        return pow(self.base(*args), self.exp)

class ReciprocalFunction(InverseElement, Map):
    """
    Multiplicative inverse of function, where the function returns
    vector and the vector-vector product is not associative.

    """
    @property
    def domain(self):
        return self.base.domain

    @property
    def codomain(self):
        return self.base.codomain

    @property
    def codomain_structure(self):
        return self.map.base_op.codomain_structure

    def _contained(self, other):
        return Map._contained(self, other)

    def eval(self, *args):
        pow = self.codomain_structure.pow
        return pow(self.base(*args), -1)
