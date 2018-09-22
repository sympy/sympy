from sympy.tensor.tensor import (TensExpr, TensMul)
from sympy import S, diag


class PartialDerivative(TensExpr):

    def __new__(cls, expr, *variables):

        # Flatten:
        if isinstance(expr, PartialDerivative):
            variables = expr.variables + variables
            expr = expr.expr

        # TODO: check that all variables have rank 1.

        args, indices, free, dum = TensMul._tensMul_contract_indices([expr] +
            list(variables), replace_indices=True)

        obj = TensExpr.__new__(cls, *args)

        obj._indices = indices
        obj._free = free
        obj._dum = dum
        return obj

    def doit(self):
        args, indices, free, dum = TensMul._tensMul_contract_indices(self.args)

        obj = self.func(*args)
        obj._indices = indices
        obj._free = free
        obj._dum = dum
        return obj

    def get_indices(self):
        return self._indices

    @property
    def expr(self):
        return self.args[0]

    @property
    def variables(self):
        return self.args[1:]

    def _extract_data(self, replacement_dict):
        from .array import derive_by_array, tensorcontraction, tensorproduct
        indices, array = self.expr._extract_data(replacement_dict)
        #import pdb; pdb.set_trace()
        for variable in self.variables:
            var_indices, var_array = variable._extract_data(replacement_dict)
            coeff_array, var_array = zip(*[i.as_coeff_Mul() for i in var_array])
            array = derive_by_array(array, var_array)
            varindex = var_indices[0]
            # TODO: find a more efficient way:
            array = tensorcontraction(tensorproduct(diag(*[S.One/i for i in coeff_array]), array), (1, 2))
            if -varindex in indices:
                pos = indices.index(-varindex)
                array = tensorcontraction(array, (0, pos+1))
                indices.pop(pos)
            else:
                indices.append(varindex)
        return indices, array
