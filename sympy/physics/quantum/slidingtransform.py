"""A sliding window tranform to process Mul-like expressions."""

from sympy.core.expr import Expr
from sympy.core.mul import Mul
from sympy.multipledispatch.dispatcher import Dispatcher


class SlidingTransform(object):
 
    def __init__(self, unary=None, binary=None, reverse=False):
        self._unary = unary
        self._binary = binary
        self._reverse = reverse

    @property
    def reverse(self):
        return self._reverse

    def __call__(self, expr):
        if not isinstance(expr, Mul):
            raise TypeError('The SlidingTransform only works on Mul instances.')

        input = list(expr.args)
        output = []

        if self._unary is not None:
            while len(input) > 0:
                next = input.pop(0)
                transformed = self._unary(next)
                if transformed is None:
                    output.append(next)
                else:
                    output.extend(transformed)
            input = output
            output = []


        if self._binary is not None:
            if not self.reverse:
                # Continue as long as we have at least 2 elements
                while len(input) > 1:
                    # Get first two elements
                    lhs = input.pop(0)
                    rhs = input[0]  # Look at second element without popping yet

                    transformed = self._binary(lhs, rhs)

                    if transformed is None:
                        # If transform returns None, append first element
                        output.append(lhs)
                    else:
                        # This item was transformed, pop and discard
                        input.pop(0)
                        # The last item goes back to be transformed again
                        input.insert(0, transformed[-1])
                        # All other items go directly into the result
                        output.extend(transformed[:-1])

                # Append any remaining element
                if input:
                    output.append(input[0])

            else:
                while len(input) > 1:
                    # Get first two elements
                    rhs = input.pop(-1)
                    lhs = input[-1]

                    transformed = self._binary(lhs, rhs)

                    if transformed is None:
                        output.insert(0, rhs)
                    else:
                        input.pop(-1)
                        input.append(transformed[0])
                        output[0:0] = transformed[1:]

                if input:
                    output.insert(0, input[-1])

            input = output
            output = []

        return Mul._from_args(input, is_commutative=False)


class DipatchingSlidingTransform(SlidingTransform):

    def __init__(self, unary=False, binary=True):
        if unary:
            unary = Dispatcher('_sliding_transform_unary')
            unary.register((Expr), lambda x: None)
        else:
            unary = None

        if binary:
            binary = Dispatcher('_sliding_transform_binary')
            binary.register((Expr, Expr), lambda x, y: None)
        else:
            binary = None

        super().__init__(unary, binary)

    @property
    def unary(self):
        return self._unary

    @property
    def binary(self):
        return self._binary
