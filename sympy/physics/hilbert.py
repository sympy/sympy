"""
Hilbert spaces for quantum mechanics [1].

References
==========

[1] http://en.wikipedia.org/wiki/Hilbert_space
"""

from sympy import Expr, Interval, oo, sympify

class HilbertSpace(Expr):
    """
    An abstract Hilbert space for quantum mechanics.

    Examples
    ========

    """

    description = 'General abstract Hilbert space.'

    def __new__(cls):
        obj = Expr.__new__(cls, **{'commutative': False})
        return obj

    @property
    def dimension(self):
        """Return the Hilbert dimension of the space."""
        raise NotImplementedError('This Hilbert space has no dimension.')

    def __add__(self, other):
        return DirectSumHilbertSpace(self, other)

    def __radd__(self, other):
        return DirectSumHilbertSpace(other, self)

    def __mul__(self, other):
        return TensorProductHilbertSpace(self, other)

    def __rmul__(self, other):
        return TensorProductHilbertSpace(other, self)

    def __pow__(self, other, mod=None):
        if mod is not None:
            raise ValueError('The third argument to __pow__ is not supported for Hilbert spaces.')
        return DirectPowerHilbertSpace(self, other)

    def __contains__(self, other):
        """Is the operator or state in this Hilbert space."""
        hs = other.hilbert_space
        if hs == self:
            return True
        else:
            return False

class l2(HilbertSpace):
    """
    l2 Hilbert space of any dimension.

    l2 (little-ell-two) is a Hilbert space of complex valued vectors. The
    number of components of the vectors in the space is the dimension of
    the Hilbert space. l2 can have a dimension of infinity, but it is a
    countable infinity (vector components labeled by the natural numbers).

    A classic example of an l2 space is spin-1/2, which is l2(2). For spin-s,
    the space is l2(2*s+1). Quantum computing with N qubits is done with
    the direct product space l2(2)**N.

    Examples
    ========

    Let's say you want to create an l2(2) Hilbert space (2 dimensional). All
    you need to do simply type these commands in:

        >>> from sympy import Expr, Interval, oo, sympify
        >>> from sympy.physics.hilbert import l2
        >>> hs = l2(2)
        >>> hs
        l2(2)

    And creating l2 Hilbert spaces of different dimensions (e.g. infinite or
    symbolic) is quite simple as well:

        >>> from sympy import oo
        >>> from sympy.abc import x
        >>> hs = l2(oo)
        >>> hs
        l2(oo)
        >>> hs = l2(x)
        >>> hs
        l2(x)

    If you want to call an l2 Hilbert space's dimension, simply follow these
    commands:

        >>> hs = l2(42)
        >>> hs.dimension
        42
    """

    def __new__(cls, dimension):
        dimension = sympify(dimension)
        if not (dimension.is_Integer or dimension is oo or dimension.is_Symbol):
            raise TypeError('l2 dimension must be an integer, oo or a Symbol: %r' % dimension)
        obj = Expr.__new__(cls, dimension, **{'commutative': False})
        return obj

    @property
    def dimension(self):
        return self.args[0]

    @property
    def description(self):
        return 'Hilbert space of length %s complex valued vectors.' % str(self.dimension)

    def _sympyrepr(self, printer, *args):
        return "l2(%s)" % printer._print(self.dimension, *args)

    def _sympystr(self, printer, *args):
        return "l2(%s)" % printer._print(self.dimension, *args)

class L2(HilbertSpace):
    """
    The Hilbert space of square integrable functions on an interval.

    L2 (big-ell-two) has a dimension of infinity. L2 is different than
    l2(oo) because the elements of L2 have an uncountable number of components,
    that is, they are functions.

    Examples
    ========

    If you want to create an L2 space over certain intervals (integers are used
    here, but you can use oo or symbols), use commands such as these:

        >>> from sympy import Expr, Interval, oo, sympify
        >>> from sympy.physics.hilbert import L2
        >>> HS = L2(Interval(-2, 42))
        >>> HS
        L2([-2, 42])

    If you want to call an L2 Hilbert space's dimension or interval, simply follow
    these commands:

        >>> HS.interval
        [-2, 42]
        >>> HS.dimension
        oo
    """

    def __new__(cls, interval):
        if not isinstance(interval, Interval):
            raise TypeError('L2 interval must be an Interval instance: %r' % interval)
        obj = Expr.__new__(cls, interval, **{'commutative': False})
        return obj

    @property
    def dimension(self):
        return oo

    @property
    def interval(self):
        return self.args[0]

    @property
    def description(self):
        return 'Hilbert space of square integrable functions on the interval %s.' % str(self.interval)

    def _sympyrepr(self, printer, *args):
        return "L2(%s)" % printer._print(self.interval, *args)

    def _sympystr(self, printer, *args):
        return "L2(%s)" % printer._print(self.interval, *args)

class FockSpace(HilbertSpace):
    """
    The Hilbert space for second quantization and field theory.

    Technically, this Hilbert space is a symmetrized/anti-symmetrized infinite
    direct sum of direct products of single particle Hilbert spaces [1]. This
    is a mess, so we have a class to represent it directly.

    Examples
    ========

    Creating a Fock space is quite simple:

        >>> from sympy.physics.hilbert import FockSpace
        >>> fs = FockSpace()
        >>> fs
        FS()

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Fock_space
    """

    def __new__(cls):
        obj = Expr.__new__(cls, **{'commutative': False})
        return obj

    @property
    def dimension(self):
        return oo

    @property
    def description(self):
        return 'Hilbert space of Fock space.'

    def _sympyrepr(self, printer, *args):
        return "FockSpace()"

    def _sympystr(self, printer, *args):
        return "FS()"

class TensorProductHilbertSpace(HilbertSpace):
    """
    A direct or tensor product of Hilbert spaces [1].

    Examples
    ========

    The tensor product between two Hilbert spaces is represented by the
    operator "*" Only the same type of Hilbert space with the same
    dimension and/or interval will be combined into a direct power,
    otherwise the tensor product behaves as follows:

        >>> from sympy import Expr, Interval, oo, sympify
        >>> from sympy.physics.hilbert import l2, L2, FockSpace, TensorProductHilbertSpace, DirectPowerHilbertSpace
        >>> hs = l2(2)
        >>> HS = L2(Interval(-42, 42))
        >>> fs = FockSpace()
        >>> tensor_product = hs*HS*fs
        >>> tensor_product
        l2(2)*L2([-42, 42])*FS()

    Here is how the properties work:

        >>> tensor_product.dimension
        oo
        >>> tensor_product.spaces
        (l2(2), L2([-42, 42]), FS())

    Identical Hilbert spaces will be combined into one direct power:

        >>> ms = hs*hs*hs
        >>> ms
        l2(2)**(3)

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Hilbert_space#Tensor_products
    """

    def __new__(cls, *args):
        r = cls.eval(args)
        if isinstance(r, Expr):
            return r
        obj = Expr.__new__(cls, *args, **{'commutative': False})
        return obj

    @classmethod
    def eval(cls, args):
        """Evaluates the direct product."""
        new_args = []
        recall = False
        #flatten arguments
        for arg in args:
            if isinstance(arg, TensorProductHilbertSpace):
                new_args.extend(arg.args)
                recall = True
            elif isinstance(arg, (HilbertSpace, DirectPowerHilbertSpace)):
                new_args.append(arg)
            else:
                raise TypeError('Hilbert spaces can only be multiplied by other Hilbert spaces: %r' % arg)
        #combine like arguments into direct powers
        comb_args = []
        prev_arg = None
        for new_arg in new_args:
            if prev_arg != None:
                if isinstance(new_arg, DirectPowerHilbertSpace) and isinstance(prev_arg, DirectPowerHilbertSpace) and new_arg.base == prev_arg.base:
                    prev_arg = new_arg.base**(new_arg.exp+prev_arg.exp)
                elif isinstance(new_arg, DirectPowerHilbertSpace) and new_arg.base == prev_arg:
                    prev_arg = prev_arg**(new_arg.exp+1)
                elif isinstance(prev_arg, DirectPowerHilbertSpace) and new_arg == prev_arg.base:
                    prev_arg = new_arg**(prev_arg.exp+1)
                elif new_arg == prev_arg:
                    prev_arg = new_arg**2
                else:
                    comb_args.append(prev_arg)
                    prev_arg = new_arg
            elif prev_arg == None:
                prev_arg = new_arg
        comb_args.append(prev_arg)
        if recall:
            return TensorProductHilbertSpace(*comb_args)
        elif len(comb_args) == 1:
            return DirectPowerHilbertSpace(comb_args[0].base, comb_args[0].exp)
        else:
            return None

    @property
    def dimension(self):
        arg_list = [arg.dimension for arg in self.args]
        if oo in arg_list:
            return oo
        else:
            return reduce(lambda x,y: x*y, arg_list)

    @property
    def description(self):
        return "A direct product Hilbert space."

    @property
    def spaces(self):
        """A tuple of the Hilbert spaces in this tensor product."""
        return self.args

    def _spaces_printer(self, printer, *args):
        spaces_strs = []
        for arg in self.args:
            s = printer._print(arg, *args)
            if isinstance(arg, DirectSumHilbertSpace):
                s = '(%s)' % s
            spaces_strs.append(s)
        return spaces_strs

    def _sympyrepr(self, printer, *args):
        spaces_reprs = self._spaces_printer(printer, *args)
        return "TensorProductHilbertSpace(%s)" % ','.join(spaces_reprs)

    def _sympystr(self, printer, *args):
        spaces_strs = self._spaces_printer(printer, *args)
        return '*'.join(spaces_strs)

class DirectSumHilbertSpace(HilbertSpace):
    """
    A direct sum of Hilbert spaces [1].

    Examples
    ========

    This class uses the "+" operator to represent direct sums between
    different Hilbert spaces. Here is a basic example:

        >>> from sympy import Expr, Interval, oo, sympify
        >>> from sympy.physics.hilbert import l2, L2, DirectSumHilbertSpace
        >>> from sympy.abc import x
        >>> hs = l2(x)
        >>> HS = L2(Interval(0, 1))
        >>> direct_sum = hs+HS
        >>> direct_sum
        l2(x)+L2([0, 1])

    Here is how the properties work:

        >>> direct_sum.dimension
        oo
        >>> direct_sum.spaces
        (l2(x), L2([0, 1]))

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Hilbert_space#Direct_sums
    """
    def __new__(cls, *args):
        r = cls.eval(args)
        if isinstance(r, Expr):
            return r
        obj = Expr.__new__(cls, *args, **{'commutative': False})
        return obj

    @classmethod
    def eval(cls, args):
        """Evaluates the direct product."""
        new_args = []
        recall = False
        #flatten arguments
        for arg in args:
            if isinstance(arg, DirectSumHilbertSpace):
                new_args.extend(arg.args)
                recall = True
            elif isinstance(arg, HilbertSpace):
                new_args.append(arg)
            else:
                raise TypeError('Hilbert spaces can only be summed with other Hilbert spaces: %r' % arg)
        if recall:
            return DirectSumHilbertSpace(*new_args)
        else:
            return None

    @property
    def dimension(self):
        arg_list = [arg.dimension for arg in self.args]
        if oo in arg_list:
            return oo
        else:
            return reduce(lambda x,y: x+y, arg_list)

    @property
    def description(self):
        return "A direct sum Hilbert space."

    @property
    def spaces(self):
        """A tuple of the Hilbert spaces in this direct sum."""
        return self.args

    def _sympyrepr(self, printer, *args):
        spaces_reprs = [printer._print(arg, *args) for arg in self.args]
        return "DirectSumHilbertSpace(%s)" % ','.join(spaces_reprs)

    def _sympystr(self, printer, *args):
        spaces_strs = [printer._print(arg, *args) for arg in self.args]
        return '+'.join(spaces_strs)

class DirectPowerHilbertSpace(HilbertSpace):
    """
    An exponentiated (iterated tensor/direct product) Hilbert space [1].

    Examples
    ========

    Direct powers (repeated tensor products) are represented by the
    operator "**" Identical Hilbert spaces that are multiplied together
    will be automatically combined into a single direct power. Here
    are some examples:

        >>> from sympy import Expr, Interval, oo, sympify
        >>> from sympy.physics.hilbert import l2, L2, DirectPowerHilbertSpace, TensorProductHilbertSpace
        >>> from sympy.abc import x
        >>> p1 = l2(3)**2
        >>> p1
        l2(3)**(2)

        >>> hs = l2(2)
        >>> direct_power = hs*hs*(hs**2)*hs**x
        >>> direct_power
        l2(2)**(4 + x)

        >>> HS = L2(Interval(-21, 21))
        >>> dir_pow = HS**(42+x)
        >>> dir_pow
        L2([-21, 21])**(42 + x)

    You can check certain properties of direct powers such as their
    dimensions, bases, exponents, etc:

        >>> p1.base
        l2(3)
        >>> p1.exp
        2
        >>> p1.dimension
        9

        >>> direct_power.dimension
        2**(4 + x)

        >>> dir_pow.base
        L2([-21, 21])
        >>> dir_pow.exp
        42 + x
        >>> dir_pow.dimension
        oo

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Hilbert_space#Tensor_products
    """

    def __new__(cls, *args):
        r = cls.eval(args)
        if isinstance(r, Expr):
            return r
        return Expr.__new__(cls, *r, **{'commutative': False})

    @classmethod
    def eval(cls, args):
        new_args = args[0], sympify(args[1])
        exp = new_args[1]
        #simplify hs**1 -> hs
        if exp == 1:
            return args[0]
        #simplify hs**0 -> 1
        if exp == 0:
            return sympify(1)
        #check (and allow) for hs**(x+42+y...) case
        if len(exp.atoms()) == 1:
            if not (exp.is_Integer and exp >= 0 or exp.is_Symbol):
                raise ValueError('Hilbert spaces can only be raised to positive integer powers or symbols: %r' % exp)
        for power in exp.atoms():
            if not (power.is_Integer or power.is_Symbol):
                raise ValueError('Hilbert spaces can only be raised to positive integer powers or symbols: %r' % exp)
        return new_args

    @property
    def base(self):
        return self.args[0]

    @property
    def exp(self):
        return self.args[1]

    @property
    def dimension(self):
        if self.base.dimension == oo:
            return oo
        else:
            return self.base.dimension**self.exp

    @property
    def description(self):
        return "An exponentiated Hilbert space."

    def _sympyrepr(self, printer, *args):
        return "DirectPowerHilbertSpace(%s,%s)" % (printer._print(self.base, *args), printer._print(self.exp, *args))

    def _sympystr(self, printer, *args):
        return "%s**(%s)" % (printer._print(self.base, *args), printer._print(self.exp, *args))
