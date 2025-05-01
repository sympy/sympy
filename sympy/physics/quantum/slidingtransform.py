"""A sliding window tranform to process Mul-like expressions."""

from itertools import tee

from sympy.core.expr import Expr
from sympy.core.mul import Mul
from sympy.multipledispatch.dispatcher import Dispatcher

from sympy.utilities.misc import debug

__all__ = [
    'SlidingTransform',
    'DipatchingSlidingTransform'
]


def _split_on_condition(seq, condition):
    l1, l2 = tee((condition(item), item) for item in seq)
    return tuple(i for p, i in l1 if p), tuple(i for p, i in l2 if not p)


def _split_cnc(seq):
    c, nc = _split_on_condition(seq, lambda x: x.is_commutative)
    return c, nc


class SlidingTransform(object):
 
    def __init__(self, unary=None, binary=None, reverse=False):
        self._unary = unary
        self._binary = binary
        self._reverse = reverse

    @property
    def reverse(self):
        return self._reverse

    @property
    def unary(self):
        return self._unary

    @property
    def binary(self):
        return self._binary

    def __call__(self, expr, **options):
        if not isinstance(expr, Mul):
            raise TypeError('The SlidingTransform only works on Mul instances.')

        input = list(expr.args)
        output = []
        # if self.reverse:
        #     debug('ST1: ', input, output)
        if self._unary is not None:
            while len(input) > 0:
                next = input.pop(0)
                transformed = self._unary(next, **options)
                if transformed is None:
                    output.append(next)
                else:
                    output.extend(transformed)
            input = output
            output = []
            # debug('qapply_Mul unary: ', input)

        # if self.reverse:
        #     debug('ST2: ', input, output)
        c_parts = []
        if self._binary is not None:
            if not self.reverse:
                # Continue as long as we have at least 2 elements
                while len(input) > 1:
                    # Get first two elements
                    lhs = input.pop(0)
                    rhs = input[0]  # Look at second element without popping yet

                    transformed = self._binary(lhs, rhs, **options)

                    if transformed is None:
                        # If transform returns None, append first element
                        output.append(lhs)
                    else:
                        c_t, nc_t = _split_cnc(transformed)
                        c_parts.extend(c_t)
                        # This item was transformed, pop and discard
                        input.pop(0)
                        if nc_t:
                            # The last item goes back to be transformed again
                            input.insert(0, nc_t[-1])
                            # All other items go directly into the result
                            output.extend(nc_t[:-1])

                # Append any remaining element
                if input:
                    output.append(input[0])

            else: # reverse = True case
                while len(input) > 1:
                    # Get first two elements
                    rhs = input.pop(-1)
                    lhs = input[-1]

                    # Make sure that lhs and rhs are commutative
                    if rhs.is_commutative:
                        c_parts.append(rhs)
                        continue
                    if lhs.is_commutative:
                        c_parts.append(lhs)
                        input.pop(-1)
                        input.append(rhs)
                        continue

                    transformed = self._binary(lhs, rhs)

                    if transformed is None:
                        output.insert(0, rhs)
                    else:
                        c_t, nc_t = _split_cnc(transformed)
                        c_parts.extend(c_t)
                        input.pop(-1)
                        if nc_t:
                            input.append(nc_t[0])
                            output[0:0] = nc_t[1:]

                if input:
                    output.insert(0, input[-1])
            # if self.reverse:
            #     debug('ST3: ', c_parts, input, output)

            input = output
            output = []

        if self.reverse:
            return Mul(*c_parts)*Mul(*input)
        else:
            return Mul(*c_parts)*Mul._from_args(input, is_commutative=False)
