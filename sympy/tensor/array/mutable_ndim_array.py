import collections

from sympy import Integer
from sympy.matrices import Matrix
from sympy.tensor.array.ndim_array import NDimArray


class MutableNDimArray(NDimArray):

    @classmethod
    def _scan_iterable_shape(cls, iterable):
        def f(pointer):
            if not isinstance(pointer, collections.Iterable):
                return [pointer], ()

            result = []
            elems, shapes = zip(*[f(i) for i in pointer])
            if len(set(shapes)) != 1:
                raise ValueError("could not determine shape unambiguously")
            for i in elems:
                result.extend(i)
            return result, (len(shapes),)+shapes[0]

        return f(iterable)

    @classmethod
    def _handle_ndarray_creation_inputs(cls, *args, **kwargs):

        # Construct N-dim array from an iterable (numpy arrays included):
        if len(args) == 1 and isinstance(args[0], collections.Iterable):
            iterable, shape = cls._scan_iterable_shape(args[0])

        # Construct N-dim array from a Matrix:
        elif len(args) == 1 and isinstance(args[0], Matrix):
            shape = args[0].shape
            iterable = args[0]

        # Construct N-dim array from another N-dim array:
        elif len(args) == 1 and isinstance(args[0], NDimArray):
            shape = args[0].shape
            iterable = args[0]

        # Construct NDimArray(dim1, dim2, dim3, ... , optional_iterable)
        elif len(args) > 1:
            shape = []
            for arg in args:
                if not isinstance(arg, (int, Integer)):
                    iterable = arg
                    break
                shape.append(arg)
            assert len(args) == len(shape) + 1

        elif len(args) == 0:
            shape = ()
            iterable = ()
        else:
            raise TypeError("Data type not understood")

        return tuple(shape), iterable
