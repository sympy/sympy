from sympy.core import Pow, S, Basic
from sympy.core.add import Add, _addsort
from sympy.core.mul import Mul, _mulsort
from sympy.core.cache import cacheit
from sympy.core.sympify import _sympify
from sympy.sets import Intersection, Union
from .map import Map, ConstantMap

__all__ = [
    'MapAdd', 'MapMul', 'MapPow',
]

class MapAdd(Map, Add):
    r"""
    Vector addition between functions.

    Explanation
    ===========

    Let $V$ be a vector space, and $V_F$ is a set of functions whose codomain is $V$.
    Then, $V_F$ forms a vector space and addition between functions can be defined: 
    $f + g : x \mapsto f(x)+g(x)$ where $f \in V_F$, $g \in V_F$ and 
    $x \in \textrm{dom}(f) \cap \textrm{dom}(g)$.

    .. note::
       Currently, ``MapAdd`` expects every functions to return scalar.

    Examples
    ========

    >>> from sympy import Sin, Cos, S, pi
    >>> sin, cos = Sin(S.Reals), Cos(S.Reals)

    >>> sin + cos
    sin + cos : Reals -> Reals

    >>> (sin + cos)(pi)
    (sin + cos)(pi)

    >>> (sin + cos)(pi, evaluate=True)
    -1

    """

    identity = None # identity depends on the domain

    @cacheit
    def __new__(cls, *args, evaluate=False, **options):

        if not args:
            raise ValueError("No argument given to %s" % cls)

        domain = Intersection(*[a.domain for a in args])
        codomain = Union(*[a.codomain for a in args])

        if evaluate:
            args, _, _ = cls.flatten(list(args))

        if len(args) == 0:
            return ConstantMap(0, domain=domain)
        if len(args) == 1:
            return args[0]

        obj = cls._from_args(args)
        obj._domain = domain
        obj._codomain = codomain
        return obj

    @classmethod
    def _from_args(cls, args, **kwargs):
        if not args:
            raise ValueError("No argument given to %s" % cls)
        elif len(args) == 1:
            return args[0]

        obj = Basic.__new__(cls, *args)
        obj.is_commutative = None
        return obj

    def _new_rawargs(cls, *args, **kwargs):
        if not args:
            raise ValueError("No argument given to %s" % cls)
        elif len(args) == 1:
            return args[0]

        obj = Basic.__new__(cls, *args)
        obj.is_commutative = None
        return obj

    @property
    def domain(self):
        return self._domain

    @property
    def codomain(self):
        return self._codomain

    @classmethod
    def flatten(cls, seq):
        terms = {}
        coeff = S.Zero # will be converted to ConstantMap

        for o in seq:

            if isinstance(o, cls):
                seq.extend(o.args)
                continue

            elif isinstance(o, ConstantMap):
                coeff += o.output
                continue

            elif isinstance(o, MapMul):
                c, s = o.as_coeff_Mul()

            else:
                c, s = S.One, o

            if s in terms:
                terms[s] += c
            else:
                terms[s] = c

        newseq = []
        for s, c in terms.items():
            if c.is_zero:
                continue
            elif c is S.One:
                newseq.append(s)
            else:
                if isinstance(s, MapMul):
                    cs = s._new_rawargs(*((c,) + s.args))
                    newseq.append(cs)
                elif isinstance(s, cls):
                    newseq.append(MapMul(c, s, evaluate=False))
                else:
                    newseq.append(MapMul(c, s, evaluate=True))

        _addsort(newseq)

        if coeff is not S.Zero:
            newseq.insert(0, ConstantMap(coeff))

        return newseq, [], None

    def eval(self, *args, **kwargs):
        kwargs.update(evaluate=True)
        terms = [a(*args, **kwargs) for a in self.args]
        return Add(*terms, evaluate=True)

    def _eval_derivative_n_times(self, index, count):
        args = [a.diff((index, count), evaluate=True) for a in self.args]
        return self.func(*args)

class MapMul(Map, Mul):
    r"""
    Scalar multiplication between scalar and function, and vector multiplication
    between function and function.

    Explanation
    ===========

    Let $V$ be an algebra [1] over a field $F$ and $V_F$ is a set of functions whose
    codomain is $V$. Then, $V_F$ forms an algebra over $F$. Therefore, scalar multiplication
    of function can be defined: $c \cdot f : x \mapsto c \cdot f(x)$ where $f \in V_F$,
    $c \in F$ and $x \in \textrm{dom}(f)$. Also, bilinear product between functions can be
    defined: $f \cdot g : x \mapsto f(x) \cdot g(x)$ where $f \in V_F$, $g \in V_F$ and 
    $x \in \textrm{dom}(f) \cap \textrm{dom}(g)$.


    .. note::
       Currently, ``MapMul`` expects every functions to return scalar.
       This means every arguments of ``MapMul`` are expected to be commutable.

    Examples
    ========

    >>> from sympy import Sin, Cos, S, pi
    >>> sin, cos = Sin(S.Reals), Cos(S.Reals)

    >>> sin*cos
    sin*cos : Reals -> Reals

    >>> (sin*cos)(pi/4)
    (sin*cos)(pi/4)

    >>> (sin * cos)(pi/4, evaluate=True)
    1/2

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Algebra_over_a_field

    """

    identity = None # identity depends on the domain

    @cacheit
    def __new__(cls, *args, evaluate=False, **options):

        if not args:
            raise ValueError("No argument given to %s" % cls)

        args = list(map(_sympify, args))
        domain = Intersection(*[a.domain for a in args if hasattr(a, 'domain')])
        codomain = Union(*[a.codomain for a in args if hasattr(a, 'codomain')])

        if evaluate:
            args, _, _ = cls.flatten(args)

        if len(args) == 0:
            return ConstantMap(1, domain=domain)
        if len(args) == 1:
            return args[0]

        obj = cls._from_args(args)
        obj._domain = domain
        obj._codomain = codomain
        return obj

    @classmethod
    def _from_args(cls, args, **kwargs):
        if not args:
            raise ValueError("No argument given to %s" % cls)
        elif len(args) == 1:
            return args[0]

        obj = Basic.__new__(cls, *args)
        obj.is_commutative = None
        return obj

    def _new_rawargs(cls, *args, **kwargs):
        if not args:
            raise ValueError("No argument given to %s" % cls)
        elif len(args) == 1:
            return args[0]

        obj = Basic.__new__(cls, *args)
        obj.is_commutative = None
        return obj

    @property
    def domain(self):
        return self._domain

    @property
    def codomain(self):
        return self._codomain

    @classmethod
    def flatten(cls, seq):
        coeff = S.One
        terms = {}
        newseq = []

        for o in seq:

            if not isinstance(o, Map) or isinstance(o, ConstantMap):

                if isinstance(o, ConstantMap):
                    o = o.output

                coeff *= o
                if coeff == 0:
                    return [ConstantMap(0)], [], None
                continue

            elif isinstance(o, cls):
                seq.extend(o.args)
                continue

            b, e = o.as_base_exp()
            if b in terms:
                terms[b] += e
            else:
                terms[b] = e

        for b, e in terms.items():
            if e.is_zero:
                continue
            elif e is S.One:
                newseq.append(b)
            else:
                newseq.append(MapPow(b, e, evaluate=True))

        _mulsort(newseq)

        if not newseq:
            newseq = [ConstantMap(coeff)]
        elif coeff is not S.One:
            newseq.insert(0, coeff)

        return newseq, [], None

    def as_coeff_Mul(self, *args, **kwargs):
        coeff, args = self.args[0], self.args[1:]
        if not isinstance(coeff, Map):
            return coeff, self.func(*args)
        return S.One, self

    def eval(self, *args, **kwargs):
        kwargs.update(evaluate=True)
        terms = []
        for a in self.args:
            if isinstance(a, Map):
                term = a(*args, **kwargs)
            else:
                term = a
            terms.append(term)
        return Mul(*terms, evaluate=True)

    def fdiff(self, index):
        terms = []
        for i, ith in enumerate(self.args):
            if not isinstance(ith, Map):
                continue
            deriv_ith = ith.diff(index, evaluate=True)
            args = [*self.args[:i]] + [deriv_ith] + [*self.args[i+1:]]
            term = self.func(*args)
            terms.append(term)
        return MapAdd(*terms)

class MapPow(Map, Pow):
    r"""
    Multiplicational power of function.

    Explanation
    ===========

    If multiplication between functions is defined, multiplicational power of
    function is naturally defined as $f^{n} : x \mapsto {f(x)}^{n}$.

    Examples
    ========

    >>> from sympy import Sin, Cos, S, pi
    >>> sin, cos = Sin(S.Reals), Cos(S.Reals)

    >>> sin**2
    sin**2 : Reals -> Reals

    >>> (sin**2)(pi/4)
    (sin**2)(pi/4)

    >>> (sin**2)(pi/4, evaluate=True)
    1/2

    """

    @cacheit
    def __new__(cls, b, e, evaluate=False):
        e = _sympify(e)
        if evaluate:
            obj = cls._new_eval(b, e)
            if obj is not None:
                return obj
        return super(Pow, MapPow).__new__(cls, b, e)

    @classmethod
    def _new_eval(cls, b, e):
        if e is S.Zero:
            return ConstantMap(S.One)
        elif e is S.One:
            return b
        obj = b._eval_power(e)
        if obj is not None:
            return obj

    @property
    def domain(self):
        return self.base.domain

    @property
    def codomain(self):
        return self.base.codomain

    def as_base_exp(self):
        return self.args

    def eval(self, *args, **kwargs):
        kwargs.update(evaluate=True)
        terms = []
        for a in self.args:
            if isinstance(a, Map):
                term = a(*args, **kwargs)
            else:
                term = a
            terms.append(term)
        return Pow(*terms, evaluate=True)

    def fdiff(self, index):
        dbase = self.base.diff(index, evaluate=True)
        return self * dbase * self.exp/self.base
