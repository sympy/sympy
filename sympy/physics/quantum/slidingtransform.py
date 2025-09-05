"""A sliding window tranform to process Mul-like expressions."""

from itertools import tee

from sympy.core.mul import Mul
from sympy.core.singleton import S

from sympy.utilities.misc import debug

__all__ = [
    'SlidingTransform',
]


def _split_on_condition(seq, condition):
    l1, l2 = tee((condition(item), item) for item in seq)
    return tuple(i for p, i in l1 if p), tuple(i for p, i in l2 if not p)


def _split_cnc(seq):
    c, nc = _split_on_condition(seq, lambda x: x.is_commutative)
    return c, nc


class SlidingTransform(object):

    def __init__(self, unary=None, binary=None, reverse=False, from_args=True):
        self._unary = unary
        self._binary = binary
        self._reverse = reverse
        self._from_args = from_args

    @property
    def reverse(self):
        return self._reverse

    @property
    def unary(self):
        return self._unary

    @property
    def binary(self):
        return self._binary

    @property
    def from_args(self):
        return self._from_args

    def __call__(self, expr, **options):
        debug(f"SlidingTransform(): {expr} {self.reverse} {self.from_args}")
        if not isinstance(expr, Mul):
            raise TypeError('The SlidingTransform only works on Mul instances.')

        input = list(expr.args)
        output = []
        c_parts = []
        if self._unary is not None:
            while len(input) > 0:
                next = input.pop(0)
                if next.is_commutative:
                    c_parts.append(next)
                    continue
                transformed = self._unary(next, **options)
                if transformed == (S.Zero,):
                    return S.Zero
                elif transformed is None:
                    output.append(next)
                else:
                    c, nc = _split_cnc(transformed)
                    output.extend(nc)
                    c_parts.extend(c)
            input = output
            output = []

        if self._binary is not None:
            if not self.reverse:
                # Continue as long as we have at least 2 elements
                while len(input) > 1:
                    # Get first two elements
                    lhs = input.pop(0)
                    rhs = input[0]  # Look at second element without popping yet
                    # Make sure that lhs and rhs are noncommutative
                    if lhs.is_commutative:
                        c_parts.append(lhs)
                        continue
                    if rhs.is_commutative:
                        c_parts.append(rhs)
                        input.pop(0)
                        input.append(lhs)
                        continue

                    transformed = self._binary(lhs, rhs, **options)

                    if transformed is None:
                        # If transform returns None, append first element
                        output.append(lhs)
                    else:
                        c_t, nc_t = _split_cnc(transformed)
                        if S.Zero in c_parts:
                            # Return immediately if we get a zero
                            return S.Zero
                        c_parts.extend(c_t)
                        # This rhs was transformed, pop and discard
                        input.pop(0)
                        if nc_t:
                            # The last item goes back to be transformed again
                            input.insert(0, nc_t[-1])
                            # All other items go directly into the result
                            output.extend(nc_t[:-1])

                # Append any remaining element
                if input:
                    last = input[0]
                    if last.is_commutative:
                        c_parts.append(last)
                    else:
                        output.append(input[0])

            else: # reverse = True case
                while len(input) > 1:
                    # Get first two elements
                    rhs = input.pop(-1)
                    lhs = input[-1]
                    # Make sure that lhs and rhs are noncommutative
                    if rhs.is_commutative:
                        c_parts.append(rhs)
                        continue
                    if lhs.is_commutative:
                        c_parts.append(lhs)
                        input.pop(-1)
                        input.append(rhs)
                        continue

                    transformed = self._binary(lhs, rhs, **options)

                    if transformed is None:
                        output.insert(0, rhs)
                    else:
                        c_t, nc_t = _split_cnc(transformed)
                        if S.Zero in c_parts:
                            # Return immediately if we get a zero
                            return S.Zero
                        c_parts.extend(c_t)
                        # The lhs item was transformed, pop and discard
                        input.pop(-1)
                        if nc_t:
                            # The first item goes back to be transformed again
                            input.append(nc_t[0])
                            # All other items go directly into the result
                            output[0:0] = nc_t[1:]

                if input:
                    last = input[-1]
                    if last.is_commutative:
                        c_parts.append(last)
                    else:
                        output.insert(0, input[-1])

            input = output
            output = []

        nc_parts = tuple(input)
        c_parts = tuple(c_parts)
                        
        # In the logic below, we go out of our ways to avoid triggering the post-processor logic
        # unless it is absolutely needed by using Mul._from_args and always validating when we
        # have an actual Mul

        # Handle the commutative part
        if len(c_parts) > 1:
            m = Mul(*c_parts)
            if isinstance(m, Mul):
                c_result_args = m.args
            else:
                c_result_args = (m,)
        else:
            c_result_args = c_parts

        # Handle the non-commutative part
        if len(nc_parts) > 1:
            if self.from_args:
                nc_result_args = nc_parts
            else:
                # This runs post-processors and may mutate the Mul into something else
                m = Mul(*nc_parts)
                if isinstance(m, Mul):
                    nc_result_args = m.args
                else:
                    nc_result_args = (m,)
        else:
            nc_result_args = nc_parts
           
        # Combine the commutative and non-commutative parts
        is_commutative = False if nc_result_args else True
        result = Mul._from_args(c_result_args + nc_result_args, is_commutative=is_commutative)
        return result
