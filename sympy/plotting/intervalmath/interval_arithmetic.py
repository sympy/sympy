"""
Interval Arithmetic for plotting.
This does not implement interval arithmetic accurately and
hence cannot be used for other purposes. If you want to use interval
arithmetic, use mpmath's interval arithmetic.
"""


from sympy.external import import_module
np = import_module('numpy')

# This module implements interval arithmetic using numpy and
#python floating points. The rounding up and down is not handled
#and hence this is not an accurate implementation of interval
#arithmetic.


#Q: Why use numpy? Why not simply use mpmath's interval arithmetic?
#A: mpmath's interval arithmetic simulates a floating point unit
#and hence is slow, while numpy evaluations are orders of magnitude
#faster.

#Q: Why create a seperate class for intervals? Why not use sympy's
#Interval Sets?
#A: The functionalities that will be required for plotting is quite
#different from what Interval Sets implement.

#Q: Why is rounding up and down according to IEEE754 not handled?
#A: It is not possible to do it in both numpy and python. An external
#library has to used, which defeats the whole purpose ie. speed. Also
#rounding is handled for very few functions in those libraries.

#Q Will my plots be affected?
#A It will not affect most of the plots. The interval arithmetic
#module based suffers the same problems as that of floating point
#arithmetic. Plotting based on mpmath will also be implemented for
#plots that would require high precision.

class interval(object):
    """ Represents an interval containing floating points as start and
    end of the interval."""
    #The is_valid variable tracks whether the interval obtained as the result
    #result of the function is in the domain and is continuous.
    #
    #True: Represents the interval result of a function is continuous and
    #in the domain of the function.
    #
    #False: The interval argument of the function was not in the domain of
    #the function, hence the is_valid of the result interval is False
    #
    #None: The function was not continuous over the interval or
    #      The function's argument interval is partly in the domain of the
    #        function

    #The comparision of two intervals returns a tuple of two 3-valued logic
    #values.

    #The first value determines the comparision as follows:
    #True: If the comparision is True throughout the intervals.
    #False: If the comparision is False throughout the intervals.
    #None: If the comparision is True for some part of the intervals.

    #The second value is determined as follows:
    #True: If both the intervals in comparision are valid.
    #False: If atleast one of the intervals is False, else
    #None

    def __init__(self, *args, **kwargs):
        self.is_valid = kwargs.pop('is_valid', True)
        if len(args) == 1:
            if isinstance(args[0], interval):
                self.start, self.end = args[0].start, args[0].end
            else:
                self.start = float(args[0])
                self.end = float(args[0])
        elif len(args) == 2:
            if args[0] < args[1]:
                self.start = float(args[0])
                self.end = float(args[1])
            else:
                self.start = float(args[1])
                self.end = float(args[0])

        else:
            raise ValueError("interval takes a maximum of two float values as arguments")

    @property
    def mid(self):
        return (self.start + self.end) / 2.0

    @property
    def width(self):
        return self.end - self.start

    def __repr__(self):
        return "interval(%f, %f)" % (self.start, self.end)

    def __str__(self):
        return "[%f, %f]" % (self.start, self.end)

    def __lt__(self, other):
        if isinstance(other, (int, float)):
            if self.end < other:
                return (True, self.is_valid)
            elif self.start > other:
                return (False, self.is_valid)
            else:
                return (None, self.is_valid)

        if isinstance(other, interval):
            if self.is_valid is False or other.is_valid is False:
                valid = False
            elif self.is_valid is None or other.is_valid is None:
                valid = None
            else:
                valid = True
            if self.end < other. start:
                return (True, valid)
            if self.start > other.end:
                return (False, valid)
            return (None, valid)
        else:
            return NotImplemented

    def __gt__(self, other):
        if isinstance(other, (int, float)):
            other = interval(other)
            return other.__lt__(self)
        elif isinstance(other, interval):
            return other.__lt__(self)
        else:
            return NotImplemented

    def __eq__(self, other):
        if isinstance(other, (int, float)):
            if self.start == other and self.end == other:
                return (True, self.is_valid)
            if self.__lt__(other)[0] is not None:
                return (False, self.is_valid)
            else:
                return (None, self.is_valid)
        if isinstance(other, interval):
            if self.is_valid is False or other.is_valid is False:
                valid = False
            elif self.is_valid is None or other.is_valid is None:
                valid = None
            else:
                valid = True
            if self.start == other.start and self.end == other.end:
                return (True, valid)
            if self.__lt__(other)[0] is not None:
                return (False, valid)
            return (None, valid)
        else:
            return NotImplemented

    def __ne__(self, other):
        if isinstance(other, (int, float)):
            if self.start == other and self.end == other:
                return (False, self.is_valid)
            if self.__lt__(other)[0] is not None:
                return (True, self.is_valid)
            else:
                return (None, self.is_valid)
        if isinstance(other, interval):
            if self.is_valid is False or other.is_valid is False:
                valid = False
            elif self.is_valid is None or other.is_valid is None:
                valid = None
            else:
                valid = True
            if self.start == other.start and self.end == other.end:
                return (False, valid)
            if not self.__lt__(other)[0] == None:
                return (True, valid)
            return (None, valid)
        else:
            return NotImplemented

    def __le__(self, other):
        if isinstance(other, (int, float)):
            if self.end <= other:
                return (True, self.is_valid)
            if self.start > other:
                return (False, self.is_valid)
            else:
                return (None, self.is_valid)

        if isinstance(other, interval):
            if self.is_valid is False or other.is_valid is False:
                valid = False
            elif self.is_valid is None or other.is_valid is None:
                valid = None
            else:
                valid = True
            if self.end <= other.start:
                return (True, valid)
            if self.start > other.end:
                return (False, valid)
            return (None, valid)
        else:
            return NotImplemented

    def __ge__(self, other):
        if isinstance(other, (int, float)):
            other = interval(other)
            return other.__le__(self)
        elif isinstance(other, interval):
            return other.__le__(self)

    def __add__(self, other):
        if isinstance(other, (int, float)):
            if self.is_valid:
                return interval(self.start + other, self.end + other)
            else:
                start = self.start + other
                end = self.end + other
                return interval(start, end, is_valid=self.is_valid)

        elif isinstance(other, interval):
            start = self.start + other.start
            end = self.end + other.end
            if self.is_valid and other.is_valid:
                return interval(start, end)
            elif self.is_valid is False or other.is_valid is False:
                return interval(start, end, is_valid=False)
            else:
                return interval(start, end, is_valid=None)
        else:
            raise NotImplemented

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, (int, float)):
            if self.is_valid:
                return interval(self.start - other, self.end - other)
            else:
                start = self.start - other
                end = self.end - other
                return interval(start, end, is_valid=self.is_valid)
        elif isinstance(other, interval):
            start = self.start - other.end
            end = self.end - other.start
            if self.is_valid and other.is_valid:
                return interval(self.start - other.end, self.end - other.start)
            elif self.is_valid is False or other.is_valid is False:
                return interval(start, end, is_valid=False)
            else:
                return interval(start, end, is_valid=None)
        else:
            raise NotImplemented

    def __rsub__(self, other):
        if isinstance(other, (int, float)):
            if self.is_valid:
                return interval(other - self.end, other - self.start)
            else:
                start = other - self.end
                end = other - self.start
                return interval(start, end, is_valid=False)
        elif isinstance(other, interval):
            start = other.start - self.end
            end = other.end - start
            if self.is_valid and other.is_valid:
                return interval(start, end)
            elif self.is_valid is False or other.is_valid is False:
                return interval(start, end, is_valid=False)
            else:
                return interval(start, end, is_valid=None)

        else:
            raise NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)

    def __neg__(self):
        if self.is_valid:
            return interval(-self.end, -self.start)
        else:
            return interval(-self.end, -self.start, is_valid=False)

    def __mul__(self, other):
        if isinstance(other, (interval, int, float)):
            other = interval(other)
            valid = True
            if self.is_valid is False or other.is_valid is False:
                valid = False
            elif self.is_valid is None or other.is_valid is None:
                valid = None

            if self in interval(0):
                #handle 0 * inf cases
                if not np.isfinite(other.start) or not np.isfinite(other.end):
                    return interval(-np.inf, np.inf, is_valid=valid)

            if other in interval(0):
                #handle 0 * inf cases
                if not np.isfinite(self.start) or not np.isfinite(self.end):
                    return interval(-np.inf, np.inf, is_valid=valid)

            if self.start >= 0:
                #positive * positive
                if other.start >= 0:
                    start = self.start * other.start
                    end = self.end * other.end
                    if np.isnan(start):
                        start = 0
                    if np.isnan(end):
                        end = np.inf
                    return interval(start, end, is_valid=valid)

                #positive * negative
                elif other.end <= 0:
                    start = self.end * other.start
                    end = self.start * other.end
                    if np.isnan(start):
                        start = -np.inf
                    if np.isnan(end):
                        end = 0
                    return interval(start, end, is_valid=valid)

                #positive * both signs
                else:
                    start = self.end * other.start
                    end = self.end * other.end
                    if np.isnan(start):
                        start = -np.inf
                    if np.isnan(end):
                        end = np.inf
                    return interval(start, end, is_valid=valid)
            elif self.end <= 0:
                # negative * positive
                if other.start >= 0:
                    start = self.start * other.end
                    end = self.end * other.start
                    if np.isnan(start):
                        start = -np.inf
                    if np.isnan(end):
                        end = 0
                    return interval(start, end, is_valid=valid)
                #negative * negative
                elif other.end <= 0:
                    start = self.end * other.end
                    end = self.start * other.start
                    if np.isnan(start):
                        start = 0
                    if np.isnan(end):
                        end = np.inf
                    return interval(start, end, is_valid=valid)
                #negative * both sign
                else:
                    start = self.start * other.end
                    end = self.start * other.start
                    if np.isnan(start):
                        start = -np.inf
                    if np.isnan(end):
                        end = np.inf
                    return interval(start, end, is_valid=valid)
            else:
                #both signs * positive
                if other.start >= 0:
                    start = self.start * other.end
                    end = self.end * other.end
                    if np.isnan(start):
                        start = -np.inf
                    if np.isnan(end):
                        end = np.inf
                    return interval(start, end, is_valid=valid)
                #both signs * negative
                elif other.end <= 0:
                    start = self.end * other.start
                    end = self.start * other.start
                    if np.isnan(start):
                        start = -np.inf
                    if np.isnan(end):
                        end = np.inf
                    return interval(start, end, is_valid=valid)
                #both signs * both signs
                else:
                    inters = []
                    inters.append(self.start * other.start)
                    inters.append(self.end * other.start)
                    inters.append(self.start * other.end)
                    inters.append(self.end * other.end)
                    if any(np.isnan(inter) for inter in inters):
                        start = -np.inf
                        end = np.inf
                    else:
                        start = max(inters)
                        end = min(inters)
                    return interval(start, end, is_valid=valid)
        else:
            return NotImplemented

    def __contains__(self, other):
        other = interval(other)
        return self.start <= other.start and other.end <= self.end

    def __rdiv__(self, other):
        if isinstance(other, (int, float)):
            other = interval(other)
            return other.__div__(self)
        elif isinstance(other, interval):
            return other.__div__(self)
        else:
            NotImplemented

    def __truediv__(self, other):
        return self.__div__(other)

    def __rtruediv__(self, other):
        return self.__rdiv__(other)

    def __div__(self, other):
        #Both None and False are handled
        if not self.is_valid:
            #Don'other divide as the value is not valid
            return interval(-np.inf, np.inf, is_valid=self.is_valid)
        if isinstance(other, (int, float)):
            if not self.is_valid:
                #Don'other divide as the value is not valid
                return interval(-np.inf, np.inf, is_valid=self.is_valid)
            else:
                if other == 0:
                    #Divide by zero encountered. valid nowhere
                    return interval(-np.inf, np.inf, is_valid=False)
                else:
                    return interval(self.start / other, self.end / other)

        elif isinstance(other, interval):
            if other.is_valid is False or self.is_valid is False:
                return interval(-np.inf, np.inf, is_valid=False)
            elif other.is_valid is None or self.is_valid is None:
                return interval(-np.inf, np.inf, is_valid=None)
            else:
                if self in interval(0):
                    if interval(0) in other:
                        return interval(-np.inf, np.inf, is_valid=None)
                    return interval(0)
               #denominator contains both signs, ie being divided by zero
               #return the whole real line with is_valid = None
                if other.start <= 0 and other.end >= 0:
                    return interval(-np.inf, np.inf, is_valid=None)

                #denominator negative
                if other.end < 0:
                    self = -self
                    other = -other

                #denominator positive
                if self.start >= 0:
                    start = self.start / other.end
                    end = self.end / other.start
                    return interval(start, end)
                elif self.end <= 0:
                    start = self.start / other.start
                    end = self.end / other.end
                    return interval(start, end)
                else:
                    #both signs
                    start = self.start / other. start
                    end = self.end / other.start
                    return interval(start, end)

    def __pow__(self, other):
        #Implements only power to an integer.
        if not self.is_valid:
            return self
        if isinstance(other, interval):
            return NotImplemented
        elif isinstance(other, (float, int)):
            if other == int(other):
                if other < 0:
                    return interval(1, 1) / self.__pow__(-other)
                else:
                    if other & 1:
                        return interval(self.start ** other, self.end ** other)
                    else:
                        #both non - positive
                        if self.end <= 0:
                            return interval(self.end ** other, self.start ** other)
                        elif self.start >= 0:
                            return interval(self.start ** other, self.end ** other)
                        else:
                            return interval(0, max(self.start ** other, self.end ** other))
            else:
                return NotImplemented
