"""Module providing the class Polynomial and low-level functions"""

class PolynomialException(Exception):
    pass

coeff_rings = ['int', 'rat', 'real', 'cplx', 'sym']
default_order = 'grevlex' # Global default, most efficient for division?

from sympy.polynomials.common import *

class Polynomial(Basic):
    """Polynomial representation in coefficient list form.

    Offers more efficient arithmetic and a unified way to handle
    polynomials and access their coefficient list as well as keep
    compatibility with other Basic objects. Input is automatically
    expanded on instantiation

    Examples:
    >>> x = Symbol('x')                  
    >>> y = Symbol('y')                  
    >>> f = Polynomial(x + 1)            
    >>> g = Polynomial(y**2 - x*y)       
    >>> s = f+g                          
    >>> s.var == [x, y]                  
    True
    >>> bool(s == y**2 - x*y + x + 1)
    True

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
                var = list(sympy_expr.atoms(type=Symbol))
                var.sort()
            if isinstance(var, Symbol):
                var = [var]
            obj.var = var
            if not (isinstance(sympy_expr, Basic)
                    and sympy_expr.is_polynomial(*var)):
                raise PolynomialException("%s is not a polynomial!" % sympy_expr)
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
                        "__neq__",
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
                pass

        # This uses the Polynomial's attributes
        return object.__getattribute__(self, name)

    def __str__(self):
        return str(self.sympy_expr)
    
    def __repr__(self):
        return "Polynomial(%s, %s, %s, %s)" % (repr(self.sympy_expr),
                  repr(self.coeffs), repr(self.var), repr(self.order))

    def __eq__(self, other):
        if isinstance(other, Polynomial):
            return self.sympy_expr == other.sympy_expr
        else:
            return self.sympy_expr == other

    def __neq__(self, other):
        if isinstance(other, Polynomial):
            return self.sympy_expr != other.sympy_expr
        else:
            return self.sympy_expr != other

    def __pos__(self):
        return self

    def __neg__(self):
        return Polynomial(sympy_expr=-self.sympy_expr,
                          coeffs=tuple([(-term[0],) + term[1:]
                                        for term in self.coeffs]),
                          var=self.var,
                          order=self.order)

    def __add__(self, other):
        # Fall back to sympy expression, when other is no Polynomial.
        if not isinstance(other, Polynomial):
            return self.sympy_expr + other

        # When coeffs are not compatible, use the sum of the
        # SymPy expressions.
        if self.var != other.var \
           or (len(self.var) > 1 and self.order != other.order):
            return Polynomial(self.sympy_expr + other.sympy_expr)

        # Check if one is 0.
        if self.coeffs[0][0] is S.Zero:
            return other
        if other.coeffs[0][0] is S.Zero:
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
        return self.__add__(other)

    def __sub__(self, other):
        return self.__add__(-other)

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, other):
        # Fall back to SymPy expression, if other is no Polynomial
        if not isinstance(other, Polynomial):
            return self.sympy_expr * other

        # When coeffs are not compatible, use the product of the
        # SymPy expressions.
        if self.var != other.var \
               or (len(self.var) > 1 and self.order != other.order):
            return Polynomial(self.sympy_expr * other.sympy_expr)

        # Check if one is 0.
        if self.coeffs[0][0] is S.Zero:
            return self
        if other.coeffs[0][0] is S.Zero:
            return other
        
        # Now we are going to do the multiplication using the coeffs.
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

        return Polynomial(sympy_expr=None,
                          coeffs=tuple(result_list),
                          var=self.var,
                          order=self.order)
    
    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, other):
        if isinstance(other, Polynomial):
            other = other.sympy_expr
        if other is S.Zero:
            return Polynomial(sympy_expr=S.One, coeffs=None,
                              var=self.var,
                              order=self.order)
        elif other is S.One:
            return self
        elif isinstance(other, Integer) and other.is_positive:
            # TODO: Implement efficient power algorithm using coeffs?
            return Polynomial(sympy_expr=self.sympy_expr**other,
                              var=self.var, order=self.order)
        else:
            # Probably not a polynomial, fall back to sympy expression.
            return self.sympy_expr**other

    def __call__(self, *point):
        # TODO: Implement Horner's method?
        if len(point) != len(self.var):
            raise PolynomialException('No proper input for evaluation.')
        result = self.sympy_expr
        for v, x in zip(self.var, point):
            result = result.subs(v, x)
        return result

    def as_integer(self):
        denom = S.One
        for term in self.coeffs:
            if not isinstance(term[0], Rational):
                print "%s is no rational coefficient!" % term[0]
                return S.Zero, self
            else: # Compute the lcm:
                denom = term[0].q*denom/numbers.gcd(int(denom), int(term[0].q))
        if denom is S.One:
            return S.One, self
        else:
            return denom, Polynomial(coeffs=tuple([((term[0]*denom),) + term[1:]
                                                   for term in self.coeffs]),
                                     var=self.var, order=self.order)

    def as_monic(self):
        lc = self.leading_coeff()
        return lc, Polynomial(coeffs=tuple([((term[0]/lc).expand(),) + term[1:]
                                            for term in self.coeffs]),
                              var=self.var, order=self.order)

    def as_primitive(self):
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
        result = 0
        for term in self.coeffs:
            if not isinstance(term[0], Integer):
                print "%s is no integer coefficient!" % term[0]
                return S.Zero
            result = numbers.gcd(result, abs(int(term[0])))
        return Integer(result)

    def diff(self, variable):
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
        return self.coeffs[0][0]

    def leading_term(self):
        return Polynomial(coeffs=(self.coeffs[0],), var=self.var,
                          order=self.order)
        
    def nth_coeff(self, *exponent):
        for term in self.coeffs:
            if term[1:] == exponent:
                return term[0]
        else: # No term with matching exponent found.
            return S.Zero

def sympy2coefficients(sympy_expr, var, order):
    """Return the tuple of coefficients and exponents.

    Currently, lexicographic ('lex'), graded lex ('grlex'), graded
    reverse lex ('grevlex') and 1-elimination ('1-el') orders are implemented.
    The list of var determines the order of the var.
    The input is assumed to be an expanded polynomial.

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
    """Makes a sympy expression out of a coefficient list."""
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
