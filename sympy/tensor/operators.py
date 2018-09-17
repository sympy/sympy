from sympy.tensor.tensor import (TensExpr, TensMul)


class PartialDerivative(TensExpr):

    def __new__(cls, expr, *variables):

        # Flatten:
        if isinstance(expr, PartialDerivative):
            variables = expr.variables + variables
            expr = expr.expr

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
