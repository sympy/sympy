from sympy import tensorproduct, MutableDenseNDimArray
from sympy.tensor.tensor import (TensExpr, TensMul, TensorIndex)


class PartialDerivative(TensExpr):
    """
    Partial derivative for tensor expressions.

    Examples
    ========

    >>> from sympy.tensor.tensor import TensorIndexType, TensorHead
    >>> from sympy.tensor.toperators import PartialDerivative
    >>> from sympy import symbols
    >>> L = TensorIndexType("L")
    >>> A = TensorHead("A", [L])
    >>> i, j = symbols("i j")

    >>> expr = PartialDerivative(A(i), A(j))
    >>> expr
    PartialDerivative(A(i), A(j))

    The ``PartialDerivative`` object behaves like a tensorial expression:

    >>> expr.get_indices()
    [i, -j]

    Indices can be contracted:

    >>> expr = PartialDerivative(A(i), A(i))
    >>> expr
    PartialDerivative(A(L_0), A(L_0))
    >>> expr.get_indices()
    [L_0, -L_0]
    """

    def __new__(cls, expr, *variables):

        # Flatten:
        if isinstance(expr, PartialDerivative):
            variables = expr.variables + variables
            expr = expr.expr

        # TODO: check that all variables have rank 1.

        args, indices, free, dum = cls._contract_indices_for_derivative(
            expr, variables)

        obj = TensExpr.__new__(cls, *args)

        obj._indices = indices
        obj._free = free
        obj._dum = dum
        return obj

    @classmethod
    def _contract_indices_for_derivative(cls, expr, variables):
        variables_opposite_valence = []
        for i in variables:
            i_free_indices = i.get_free_indices()
            variables_opposite_valence.append(i.xreplace({k: -k for k in i_free_indices}))

        args, indices, free, dum = TensMul._tensMul_contract_indices(
            [expr] + variables_opposite_valence, replace_indices=True)

        for i in range(1, len(args)):
            i_indices = args[i].get_free_indices()
            args[i] = args[i].xreplace({k: -k for k in i_indices})

        return args, indices, free, dum

    def doit(self):
        args, indices, free, dum = self._contract_indices_for_derivative(self.expr, self.variables)

        obj = self.func(*args)
        obj._indices = indices
        obj._free = free
        obj._dum = dum
        return obj

    def get_indices(self):
        return self._indices

    def get_free_indices(self):
        free = sorted(self._free, key=lambda x: x[1])
        return [i[0] for i in free]

    @property
    def expr(self):
        return self.args[0]

    @property
    def variables(self):
        return self.args[1:]

    def _extract_data(self, replacement_dict):
        from .array import derive_by_array, tensorcontraction
        indices, array = self.expr._extract_data(replacement_dict)
        for variable in self.variables:
            var_indices, var_array = variable._extract_data(replacement_dict)
            var_indices = [-i for i in var_indices]
            coeff_array, var_array = zip(*[i.as_coeff_Mul() for i in var_array])
            array = derive_by_array(array, var_array)
            array = array.as_mutable()  # type: MutableDenseNDimArray
            varindex = var_indices[0]  # type: TensorIndex
            # Remove coefficients of base vector:
            coeff_index = [0] + [slice(None) for i in range(len(indices))]
            for i, coeff in enumerate(coeff_array):
                coeff_index[0] = i
                array[tuple(coeff_index)] /= coeff
            if -varindex in indices:
                pos = indices.index(-varindex)
                array = tensorcontraction(array, (0, pos+1))
                indices.pop(pos)
            else:
                indices.append(varindex)
        return indices, array
