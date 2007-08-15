"""Univariate polynomials in sparse representation"""


class SparsePolynomial(object):
    """Dictionary-based polynomial representation. (abstract class)"""

    coeff_type = int
    zero = 0

    def __init__(self, coeffs={}):
        self.coeffs = coeffs
        if coeffs:
            self.degree = max(coeffs.iterkeys())
        else:
            self.degree = -1

    def __repr__(self):
        return "%s(%s)" % (self.__class__, repr(self.coeffs))

    def __getitem__(self, item):
        """Returns the nth coefficient."""
        return self.coeffs.get(item, self.zero)

    def __eq__(self, other):
        return self.coeffs == other.coeffs

    def __ne__(self, other):
        return self.coeffs != other.coeffs

    def __nonzero__(self):
        return bool(self.coeffs)

    def __pos__(self):
        return self

    def __neg__(self):
        result_dict = {}
        for key, value in self.coeffs.iteritems():
            result_dict[key] = -value
        return self.__class__(result_dict)

    def scale(self, coefficient, exponent=0):
        result_dict = {}
        if coefficient:
            for e, c in self.coeffs.iteritems():
                result_dict[e + exponent] = c * coefficient
        return self.__class__(result_dict)

    def __add__(self, other):
        result_dict = {}
        result_dict.update(self.coeffs)
        for exponent, coefficient in other.coeffs.iteritems():
            if exponent in result_dict:
                value = result_dict[exponent] + coefficient
                if value:
                    result_dict[exponent] = value
                else:
                    del result_dict[exponent]
            else:
                result_dict[exponent] = coefficient
        return self.__class__(result_dict)

    def __sub__(self, other):
        result_dict = {}
        result_dict.update(self.coeffs)
        for exponent, coefficient in other.coeffs.iteritems():
            if exponent in result_dict:
                value = result_dict[exponent] - coefficient
                if value:
                    result_dict[exponent] = value
                else:
                    del result_dict[exponent]
            else:
                result_dict[exponent] = - coefficient
        return self.__class__(result_dict)

    def __mul__(self, other):
        result = self.__class__()
        for e, c in other.coeffs.iteritems():
            result += self.scale(c, e)
        return result

    def __pow__(self, exponent):
        """Repeated Squaring."""
        assert isinstance(exponent, (int, long)) and exponent >= 0
        if exponent == 0:
            return self.__class__({0: self.coeff_type(1)})
        elif exponent == 1:
            return self
        binary_repr = []
        while exponent:
            if exponent % 2:
                binary_repr.insert(0, 1)
                exponent = (exponent - 1) / 2
            else:
                binary_repr.insert(0, 0)
                exponent /= 2
        result = self
        for k in binary_repr[1:]:
            result *= result
            if k:
                result *= self
        return result

    def diff(self):
        result_dict = {}
        for e, c in self.coeffs.iteritems():
            if e > 0:
                coeff = c*self.coeff_type(e)
                if coeff:
                    result_dict[e-1] = coeff
        return self.__class__(result_dict)

    def evaluate(self, point):
        result = self.zero
        for e, c in self.coeffs.iteritems():
            result += c*point**e
        return result
