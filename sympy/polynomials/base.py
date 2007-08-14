"""Module providing the class Polynomial and low-level functions"""

from sympy.core import *
# Use (faster) Singleton comparisons.
from sympy.core.basic import S
# Need numbers.gcd, for content.
from sympy.core import numbers
# To determine coefficient type:
from sympy.core.numbers import NumberSymbol, ImaginaryUnit
from sympy.utilities import *


# This is the list of possible rings the coefficients could lie in,
# ordered by inclusion.
coeff_rings = ['int', 'rat', 'real', 'cplx', 'sym']


# Global default, probably the most efficient for division?
default_order = 'grevlex' 


# Local exception type.
class PolynomialException(Exception):
    pass


class Polynomial(Basic):
    """Unified representation of all kinds of polynomials.

    Usage:
    ======
        Most of the time, the Polynomial instances are created of a
        SymPy expression, which is a polynomial.

        Optionally, the user can give a 'var' argument containing the
        list of variables (in order), so that he can view 2*x + x*y as
        a polynomial in x only, for example. If not given, the
        occuring symbols are extracted automatically and sorted
        alphabetically. 

        The (optional) argument 'order' defines the monomial order,
        which is of importance to division algorithms and Groebner
        bases and defaults to 'grevlex', that is graded reverse
        lexicographic ordering. Other options are:
        'lex' - lexicographic order
        'grlex' - graded lexicographic order
        '1-el' - first elimination order

        Alternatively, a Polynomial can be instantiated by giving its
        coefficients and exponents in nested tuples, alone or in
        addidition to the sympy_expr. Here, no consistency checks are
        done as this mode is intended for internal use mostly.

        The built Polynomial instances support arithmetic, like
        addition and multiplication, but fall back to the underlying
        SymPy expression, when a non-Polynomial is encountered.

        The SymPy expression of a Polynomial f can be accessed through
        f.sympy_expr. The coefficients and exponents are held in
        f.coeffs and f.var and f.order hold the respective arguments.

    Notes:
    ======
        Computes the coefficients with the exponents in sparse
        representation for faster arithmetic and leading terms etc.
        Tries to be compatible with other SymPy expressions, for
        example, by forwarding most attributes like assumptions to the
        underlying SymPy expression. 

    Examples:
    =========
        >>> x, y = symbols('xy')                  
        >>> f = Polynomial(x + 1)
        >>> f.sympy_expr
        1 + x
        >>> f.coeffs
        ((1, 1), (1, 0))
        >>> f.var
        [x]
        >>> f.order
        'grevlex'
        >>> print f
        1 + x
        >>> f
        Polynomial(1 + x, ((1, 1), (1, 0)), [x], 'grevlex')
        >>> g = Polynomial(y**2 - x*y)       
        >>> s = f + g                          
        >>> s.var == [x, y]                  
        True
        >>> bool(s == y**2 - x*y + x + 1)
        True
        >>> h = Polynomial(g.sympy_expr, var=y)
        >>> g.coeffs
        ((-1, 1, 1), (1, 0, 2))
        >>> h.coeffs
        ((1, 2), (-x, 1))

    Also see L{sympy2coefficients}, L{coefficients2sympy}.

    """
    
    def __new__(cls, sympy_expr=None, coeffs=None, var=None, order=None,
                **assumptions):
        obj = Basic.__new__(cls)

        if sympy_expr is None and coeffs is None:
            raise PolynomialException("No polynomial data given!")
        elif sympy_expr is not None and coeffs is None:
            # Polynomial is contructed by the sympy expression.
            sympy_expr = sympify(sympy_expr).expand()
            if var is None:
                # Automatically extract the variables from the expression.
                var = list(sympy_expr.atoms(type=Symbol))
                var.sort()
            if isinstance(var, Symbol):
                var = [var]
            obj.var = var
            if var and not sympy_expr.is_polynomial(*var):
                raise PolynomialException("%s is not a polynomial!"
                                          % sympy_expr)
            obj.sympy_expr = sympy_expr
            if order is None:
                order = default_order
            obj.order = order
            obj.coeffs = sympy2coefficients(sympy_expr, var, order)
        elif sympy_expr is None and coeffs is not None:
            # This polynomial is constructed by its coeffs.
            # No sympify is used, all Terms are assumed to be
            # instances of Basic.
            obj.coeffs = coeffs
            if var is None:
                raise PolynomialException("Ambiguous data, "
                                          + "please specify the var.")
            if isinstance(var, Symbol):
                var = [var]
            obj.var = var
            if order is None:
                order = default_order
            obj.order = order
            obj.sympy_expr = coefficients2sympy(coeffs, var)
        else:
            # Both the sympy expression and the coeffs are given.
            # No sanity check done, which would be slower than just
            # giving one of them.
            obj.sympy_expr = sympy_expr
            obj.coeffs = coeffs
            if var is None:
                var = list(sympy_expr.atoms(type=Symbol))
                var.sort()
            if isinstance(var, Symbol):
                var = [var]
            obj.var = var
            if order is None:
                order = default_order
            obj.order = order
        return obj


    def __getattribute__(self, name):
        """Redirect most attributes to the underlying SymPy expression."""

        # Check if the attribute belongs to the Polynomial itself:
        if name not in ("__class__",
                        "__dict__",
                        "__new__",
                        "__getattribute__",
                        "__str__",
                        "__repr__",
                        "__eq__",
                        "__ne__",
                        "__pos__",
                        "__neg__",
                        "__add__",
                        "__radd__",
                        "__sub__",
                        "__rsub__",
                        "__mul__",
                        "__rmul__",
                        "__call__",
                        "as_integer",
                        "as_monic",
                        "as_primitive",
                        "content",
                        "diff",
                        "leading_coeff",
                        "leading_term",
                        "nth_coeff",
                        "coeffs",
                        "order",
                        "sympy_expr",
                        "var"):
            try:
                # This fails when the Polynomial is instantiated.
                se = object.__getattribute__(self, "__dict__")["sympy_expr"]
                # This uses the SymPy expressions' attributes
                return object.__getattribute__(se, name)
            except KeyError:
                # The .sympy_expr doesn't yet exist.
                pass

        # This uses the Polynomial's attributes
        return object.__getattribute__(self, name)


    def __str__(self):
        """Return only the SymPy expression to be human-readable."""
        return str(self.sympy_expr)
    

    def __repr__(self):
        """Returns a string that could be used to reconstruct this object."""
        return "Polynomial(%s, %s, %s, %s)" % (repr(self.sympy_expr),
                  repr(self.coeffs), repr(self.var), repr(self.order))


    def __eq__(self, other):
        """Equality is restricted to equality of SymPy expression.

        This overwrites the == operator. Other attributes such as
        variables or monomial order are not compared.

        """

        if isinstance(other, Polynomial):
            return self.sympy_expr == other.sympy_expr
        else:
            return self.sympy_expr == other


    def __ne__(self, other):
        """Also see L{__eq__}."""
        if isinstance(other, Polynomial):
            return self.sympy_expr != other.sympy_expr
        else:
            return self.sympy_expr != other


    def __pos__(self):
        """Just returns the Polynomial."""
        return self


    def __neg__(self):
        """Returns the Polynomial multiplied by -1."""
        return Polynomial(sympy_expr=-self.sympy_expr,
                          coeffs=tuple([(-term[0],) + term[1:]
                                        for term in self.coeffs]),
                          var=self.var,
                          order=self.order)


    def __add__(self, other):
        """Overwrites the + operator.

        Implements an addition algorithm for instances of Polynomial
        with matching variables, using the coefficients and exponents,
        but falls back to the SymPy expressions otherwise. It even
        returns a non-Polynomial object when encountering one.

        Also see L{Polynomial}, L{__mul__}.

        """

        # Fall back to sympy expression, when other is no Polynomial.
        if not isinstance(other, Polynomial):
            return self.sympy_expr + other

        # When coeffs are not compatible, use the sum of the
        # SymPy expressions.
        if self.var != other.var \
           or (len(self.var) > 1 and self.order != other.order):
            return Polynomial(self.sympy_expr + other.sympy_expr)

        # Check if one is 0, then return the other.
        if self.sympy_expr is S.Zero:
            return other
        if other.sympy_expr is S.Zero:
            return self
        
        # Now we are going to do the addition using the coeffs.
        # Merge the terms of self and other:
        result_list = []
        s, o = self.coeffs, other.coeffs
        i, j = 0, 0
        while i < len(s) and j < len(o):
            if (s[i][1:] == o[j][1:]):
                c = s[i][0] + o[j][0]
                if c is not S.Zero:
                    result_list.append((c,) + s[i][1:])
                i += 1
                j += 1
            elif term_cmp(s[i], o[j], self.order) > 0:
                result_list.append(s[i])
                i += 1
            else:
                result_list.append(o[j])
                j += 1
        # Append the rest.
        result_list += s[i:]
        result_list += o[j:]
        # Check if something was appended to the (empty) result.
        if len(result_list) == 0:
            return Polynomial(S.Zero, var=self.var, order=self.order)

        return Polynomial(coeffs=tuple(result_list), var=self.var,
                          order=self.order)
        

    def __radd__(self, other):
        """Also see L{__add__}."""
        return self.__add__(other)


    def __sub__(self, other):
        """Also see L{__add__},  L{__neg__}."""
        return self.__add__(-other)


    def __rsub__(self, other):
        """Also see L{__add__}, L{__neg__}"""
        return (-self).__add__(other)


    def __mul__(self, other):
        """Overwrites the * operator.

        Implements a multiplication algorithm for instances of Polynomial
        with matching variables, using the coefficients and exponents,
        but falls back to the SymPy expressions otherwise. It even
        returns a non-Polynomial object when encountering one.

        Also see L{Polynomial}, L{__add__}.

        """

        # Fall back to SymPy expression, if other is no Polynomial.
        if not isinstance(other, Polynomial):
            return self.sympy_expr * other

        # When coeffs are not compatible, use the product of the
        # SymPy expressions.
        if self.var != other.var \
               or (len(self.var) > 1 and self.order != other.order):
            return Polynomial(self.sympy_expr * other.sympy_expr)

        # Check if one is 0, then return 0.x
        if self.sympy_expr is S.Zero:
            return self
        if other.sympy_expr is S.Zero:
            return other

        # Check if one is 1, then return other.
        if self.sympy_expr is S.One:
            return other
        if other.sympy_expr is S.One:
            return self
        
        # Now we are going to do the multiplication using the coefficients.
        result_dict = {}
        for self_term in self.coeffs:
            for other_term in other.coeffs:
                key = tuple([i+j for i,j in zip(self_term[1:], other_term[1:])])
                if result_dict.has_key(key):
                    result_dict[key] += (self_term[0]*other_term[0]).expand()
                else:
                    result_dict[key] = (self_term[0]*other_term[0]).expand()
        result_list = [(result_dict[key],) + key for key in result_dict
                       if result_dict[key] is not S.Zero]
        result_list.sort(cmp=lambda x,y:term_cmp(x,y,self.order), reverse=True)

        return Polynomial(coeffs=tuple(result_list),
                          var=self.var,
                          order=self.order)
    

    def __rmul__(self, other):
        """Also see L{__mul__}."""
        return self.__mul__(other)


    def __pow__(self, other):
        """Overwrites the ** operator."""
        if isinstance(other, Polynomial):
            other = other.sympy_expr
        if other is S.Zero:
            return Polynomial(S.One, var=self.var, order=self.order)
        elif other is S.One:
            return self
        elif isinstance(other, Integer) and other.is_positive:
            # TODO: Implement efficient power algorithm using coeffs?
            return Polynomial(self.sympy_expr**other,
                              var=self.var, order=self.order)
        else:
            # Probably not a polynomial, fall back to sympy expression.
            return self.sympy_expr**other


    def __call__(self, *point):
        """Evaluate the polynomial function at a specific point.

        Usage:
        ======
            Give an arbitrary argument for each variable and get the
            SymPy expression with the variables substituted.

        Examples:
        =========
            >>> x, y = symbols('xy')
            >>> f = Polynomial(2*x - y)
            >>> f(1, 7)
            -5
            >>> f(3*x, x)
            5*x
            
        """
        
        if len(point) != len(self.var):
            raise PolynomialException('No proper input for evaluation.')
        result = self.sympy_expr
        for v, x in zip(self.var, point):
            result = result.subs(v, x)
        return result


    def as_integer(self):
        """Return the polynomial with integer coefficients.

        Usage:
        ======
            Starting from an instance of Polynomial with only rational
            coefficients, this function multiplies it with the common
            denominator. The result is a tuple consisting of the
            factor applied to the coefficients and a new instance of
            Polynomial in the integers.

        Example:
        ========
            >>> x, y = symbols('xy')
            >>> f = Polynomial(x/6 + y/4 + x*y/3)
            >>> denominator, f = f.as_integer()
            >>> print denominator
            12
            >>> print f
            2*x + 3*y + 4*x*y
            >>> denominator, f = f.as_integer()
            >>> print denominator
            1
            >>> print f
            2*x + 3*y + 4*x*y

        Also see L{as_monic}, L{as_primitive}.

        """
        
        denom = S.One
        for term in self.coeffs:
            if not isinstance(term[0], Rational):
                print "%s is no rational coefficient!" % term[0]
                return S.Zero, self
            else:
                # Compute the least common multiple of the denominators:
                denom = term[0].q*denom/numbers.gcd(int(denom), int(term[0].q))
        if denom is S.One:
            return S.One, self
        else:
            return (denom,
                    Polynomial(coeffs=tuple([((term[0]*denom),) + term[1:]
                                             for term in self.coeffs]),
                               var=self.var, order=self.order))


    def as_monic(self):
        """Return the polynomial with leading coefficient 1.

        Usage:
        ======
            Starting with any instance of Polynomial, this returns the
            former leading coefficient and a new Polynomial which is
            monic.

        Examples:
        =========
            >>> x, y = symbols('xy')
            >>> f = Polynomial(x/2 + y/4 + x*y/3)
            >>> leadcoeff, f = f.as_monic()
            >>> print leadcoeff
            1/3
            >>> print f
            (3/2)*x + (3/4)*y + x*y
            >>> leadcoeff, f = f.as_monic()
            >>> print leadcoeff
            1
            >>> print f
            (3/2)*x + (3/4)*y + x*y

        Also see L{as_integer}, L{as_primitive}, L{leading_coeff}.

        """

        lc = self.leading_coeff()
        return lc, Polynomial(coeffs=tuple([((term[0]/lc).expand(),) + term[1:]
                                            for term in self.coeffs]),
                              var=self.var, order=self.order)


    def as_primitive(self):
        """Return the content and a primitive Pxolynomial.

        Usage:
        ======
            Starting with any instance of Polynomial, this returns the
            content, that is, the greatest common divisor of the
            (integer) coefficients, and a new Polynomial which is
            primitive, that is, of content 1. Only works for integer
            coefficients.

        Examples:
        =========
            >>> x, y = symbols('xy')
            >>> f = Polynomial(6*x + 20*y + 4*x*y)
            >>> content, f = f.as_primitive()
            >>> print content
            2
            >>> print f
            3*x + 10*y + 2*x*y
            >>> content, f = f.as_primitive()
            >>> print content
            1
            >>> print f
            3*x + 10*y + 2*x*y

        Also see L{as_integer}, L{as_monic}, L{content}.

        """

        c = self.content()
        if c is S.Zero:
            return S.Zero, self
        elif c is S.One:
            return S.One, self
        else:
            return c, Polynomial(coeffs=tuple([(term[0]/c,) + term[1:]
                                               for term in self.coeffs]),
                                 var=self.var, order=self.order)


    def content(self):
        """Return the content of a Polynomial.

        Usage:
        ======
            Returns the content, that is, the greatest common divisor of the
            (integer) coefficients.

        Examples:
        =========
            >>> x, y = symbols('xy')
            >>> f = Polynomial(6*x + 20*y + 4*x*y)
            >>> f.content()
            2

        Also see L{as_primitive}, L{leading_coeff}.

        """

        result = 0
        for term in self.coeffs:
            if not isinstance(term[0], Integer):
                print "%s is no integer coefficient!" % term[0]
                return S.Zero
            result = numbers.gcd(result, abs(int(term[0])))
        return Integer(result)


    def diff(self, variable):
        """Derivative of a Polynomial.

        Usage:
        ======
            Returns a new instance of Polynomial which is the partial
            derivative by the given variable.

        Examples:
        =========
            >>> x, y, z = symbols('xyz')
            >>> f = Polynomial(6*x + 20*y + 4*x*y)
            >>> fx = f.diff(x)
            >>> print fx
            6 + 4*y
            >>> fz = f.diff(z)
            >>> print fz
            0

        """

        if not variable in self.var:
            return Polynomial(S.Zero, var=self.var, order=self.order)
        for i, v in enumerate(self.var):
            if v is variable:
                i += 1 # Variables begin at index 1 in coeffs.
                break
        result_list = [((term[0]*term[i]).expand(),) + term[1:i]
                       + ((term[i]-1),) + term[i+1:]
                       for term in self.coeffs if term[i].is_positive]
        if len(result_list) == 0:
            return Polynomial(sympy_expr=S.Zero, var=self.var,
                              order=self.order)
        else:
            result_list.sort(cmp=lambda x,y:term_cmp(x,y,self.order),
                             reverse=True)
            return Polynomial(coeffs=tuple(result_list), var=self.var,
                              order=self.order)


    def leading_coeff(self):
        """Return the leading coefficient of a Polynomial.

        Usage:
        ======
            This gives the coefficient, that is, non-symbolic part, of
            the leading term, according to the monomial order, or
            simply highest degree, in the univariate case.

        Examples:
        =========
            >>> x, y = symbols('xy')
            >>> f = Polynomial(6*x + 20*y + 4*x*y)
            >>> f.leading_coeff()
            4

        Also see L{as_monic}, L{leading_term}, L{nth_coeff}.

        """

        return self.coeffs[0][0]


    def leading_term(self):
        """Return the leading term of a Polynomial.

        Usage:
        ======
            The leading term, according to the monomial order, or
            simply highest degree, in the univariate case.

        Examples:
        =========
            >>> x, y = symbols('xy')
            >>> f = Polynomial(6*x + 20*y + 4*x*y)
            >>> print f.leading_term()
            4*x*y

        Also see L{leading_coeff}.

        """

        return Polynomial(coeffs=(self.coeffs[0],), var=self.var,
                          order=self.order)
        

    def nth_coeff(self, *exponent):
        """Return a specific coefficient of a Polynomial.

        Usage:
        ======
            This gives the coefficient, that is, non-symbolic part, of
            the term with matching exponents, or 0, if it doesn't appear.

        Examples:
        =========
            >>> x, y = symbols('xy')
            >>> f = Polynomial(6*x + 20*y + 4*x*y)
            >>> f.nth_coeff(1, 0)
            6
            >>> f.nth_coeff(1, 1)
            4
            >>> f.nth_coeff(0, 0)
            0

        Also see L{leading_coeff}.

        """

        for term in self.coeffs:
            if term[1:] == exponent:
                return term[0]
        else: # No term with matching exponent found.
            return S.Zero


def sympy2coefficients(sympy_expr, var, order):
    """Return the tuple of coefficients and exponents.

    Usage:
    ======
        This functions computes the tuples of coefficients and
        exponents from a given SymPy expression. This expression is
        assumed to be expanded already. The arguments 'var' and
        'order' define the occuring variables and the monomial order
        to use, respectively.

        Normally, the user would never call this function himself, it
        is rather used within the Polynomial class.

    Examples:
    =========
        >>> x, y = symbols('xy')
        >>> sympy2coefficients(2*x + 3, [x], 'lex')
        ((2, 1), (3, 0))
        >>> sympy2coefficients(x**2*y + 4*y, [x, y], 'lex')
        ((1, 2, 1), (4, 0, 1))
        >>> sympy2coefficients(x**2*y + 4*y, [y], 'lex')
        ((4 + x**2, 1),)

    Also see L{Polynomial}, L{sympy2coefficients}.

    """

    result_dict = {}
    if isinstance(sympy_expr, Add):
        terms = sympy_expr[:]
    else:
        terms = [sympy_expr]
    for term in terms:
        if isinstance(term, Mul):
            factors = term[:]
        else:
            factors = [term]
        c = S.One
        exponents = [S.Zero]*len(var)
        for factor in factors:
            # Check if any of the variables occur.
            if filter(lambda x:x in var, factor.atoms(type=Symbol)):
                if isinstance(factor, Pow) \
                   and isinstance(factor.base, Symbol) \
                   and factor.exp.is_integer \
                   and factor.exp.is_positive:
                    exponents[var.index(factor.base)] += factor.exp
                elif isinstance(factor, Symbol):
                    exponents[var.index(factor)] += 1
                else:
                    raise PolynomialException("%s is not a polynomial!"
                                              % sympy_expr)
            else: # The factor is relativly constant.
                c *= factor
        exponents = tuple(exponents)
        if result_dict.has_key(exponents):
            result_dict[exponents] += c
        else:
            result_dict[exponents] = c

    coefficient_list = [(result_dict[key],) + key
                        for key in result_dict.keys()
                        if result_dict[key] is not S.Zero]
    coefficient_list.sort(cmp=lambda x,y:term_cmp(x, y, order), reverse=True)
    if len(coefficient_list) == 0:
        coefficient_list = [tuple([S.Zero]*(len(var) + 1))]
    return tuple(coefficient_list)


def coefficients2sympy(coeffs, var):
    """Return the SymPy expression of given coefficients and exponents.

    Usage:
    ======
        This functions computes the original SymPy expression from the
        tuples of coefficients and exponents. The argument 'var'
        defines the occuring variables.

        Normally, the user would never call this function himself, it
        is rather used within the Polynomial class.

    Examples:
    =========
        >>> x, y = symbols('xy')
        >>> coefficients2sympy(((2, 1), (3, 0)), [x])
        3 + 2*x
        >>> coefficients2sympy(((1, 2, 1), (4, 0, 1)), [x, y])
        4*y + y*x**2
        >>> coefficients2sympy(((1, 2, 1), (4, 0, 1)), [y, x])
        4*x + x*y**2

    Also see L{Polynomial}, L{sympy2coefficients}.

    """

    if len(coeffs) == 0:
        raise PolynomialException('Bad coefficient list.')
    elif len(coeffs[0]) != len(var) + 1:
        raise PolynomialException('Wrong number of var given.')

    result = S.Zero
    for term in coeffs:
        c = term[0]
        for i, v in enumerate(var):
            c *= v**term[i+1]
        result += c
    return result


# Simple helper functions common to several algorithms.

def reverse(t):
    """Return a tuple with reversed order"""
    return tuple([t[i] for i in range(len(t)-1, 0, -1)])


def term_cmp(a, b, order):
    """Compares tuples occuring in the Polynomial's coeffs."""
    if order == 'lex':
        return cmp(a[1:], b[1:])
    elif order == 'grlex':
        return cmp((sum(a[1:]),) + a[1:], (sum(b[1:]),) + b[1:])
    elif order == 'grevlex':
        return cmp((sum(a[1:]),) + reverse(map(lambda l:-l, a[1:])),
                   (sum(b[1:]),) + reverse(map(lambda l:-l, b[1:])))
    elif order == '1-el':
        return cmp((a[1], sum(a[2:])) + reverse(map(lambda l:-l,a[2:])),
                   (b[1], sum(b[2:])) + reverse(map(lambda l:-l,b[2:])))
    else:
        raise PolynomialException(str(order) + 'is not an implemented order.')


def term_mult(a, b):
    """Multiplication of a tuple representing a term.

    a and b are assumed to be tuples of some Polynomial's coeffs of
    same length.
    
    """

    return ((a[0]*b[0]).expand(),) \
           + tuple(map(lambda (x,y): x+y, zip(a[1:], b[1:])))

def term_div(a, b):
    """Division of a tuple representing a term.

    a and b are assumed to be tuples of some Polynomial's coeffs of
    same length.
    
    """

    return ((a[0]/b[0]).expand(),) \
           + tuple(map(lambda (x,y): x-y, zip(a[1:], b[1:])))


def term_is_mult(a, b):
    """Return True if a is a multiple of b.

    a and b are assumed to be tuples some Polynomial's coeffs same
    length.

    """

    return all([x.is_nonnegative for x in term_div(a, b)[1:]])


def term_lcm(a, b):
    """Least common multiple of tuples representing terms.

    a and b are assumed to be tuples some Polynomial's coeffs same
    length. The coefficient is set to 1.

    """

    return (S.One,) + tuple([max(aa, bb) for aa, bb in zip(a[1:], b[1:])])


def merge_var(*a):
    """Return a sorted list of the symbols in the arguments"""
    result = []
    for var in a:
        for sym in var:
            if not sym in result:
                result.append(sym)
    result.sort()
    return result


def coeff_ring(atom):
    """Determine the number type of some atom, or some list of atoms."""
    if isinstance(atom, (Number, NumberSymbol, ImaginaryUnit)) \
        or (isinstance(atom, (Add, Mul))
            and all(map(lambda a:isinstance(a, (Number, NumberSymbol,
                                                ImaginaryUnit)), atom[:]))) \
        or (isinstance(atom, Pow) and isinstance(atom.base, (Number,
                                                             NumberSymbol,
                                                             ImaginaryUnit))
            and isinstance(atom.exp, (Number, NumberSymbol, ImaginaryUnit))):
        if atom.is_integer:
            return 'int'
        elif atom.is_rational:
            return 'rat'
        elif atom.is_real:
            return 'real'
        else:
            return 'cplx'
    elif isinstance(atom, list):
        # Get the coefficient ring of each atom and look for the worst case.
        result = 'int'
        for a in atom:
            cr = coeff_ring(a)
            assert type(cr) == str
            if coeff_rings.index(cr) > coeff_rings.index(result):
                result = cr
        return result
    else:
        return 'sym'


def get_numbers(atom):
    """Extracts the non-symbolic part of an expression."""
    # TODO: Merge with coeff_ring, without recursion!
    result = []
    if isinstance(atom, (Number, NumberSymbol, ImaginaryUnit)) \
        or (isinstance(atom, (Add, Mul))
            and all(map(lambda a:isinstance(a, (Number,
                                                NumberSymbol,
                                                ImaginaryUnit)),
                        atom[:]))) \
        or (isinstance(atom, Pow) and isinstance(atom.base, (Number,
                                                             NumberSymbol,
                                                             ImaginaryUnit))
            and isinstance(atom.exp, (Number, NumberSymbol, ImaginaryUnit))):
        return [atom]
    elif isinstance(atom, (Add, Mul)):
        for a in atom:
            result.append(get_numbers(a))
    return result


def integer_divisors(n):
    """Returns a list of all positive integer divisors of n."""
    n = abs(n)
    r = []
    for i in range(1, n/2+1):
        if n % i == 0:
            r.append(i)
    r.append(n)
    return r
