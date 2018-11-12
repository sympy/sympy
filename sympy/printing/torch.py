from sympy.printing.pycode import AbstractPythonCodePrinter


class TorchPrinter(AbstractPythonCodePrinter):
    printmethod = "_torchcode"

    def _print_MatMul(self, expr):
        return self._expand_fold_binary_op("torch.mm", expr.args)

    def _print_MatPow(self, expr):
        return self._expand_fold_binary_op("torch.mm", [expr.base]*expr.exp)

    def _print_MatrixBase(self, expr):
        data = "["+", ".join(["["+", ".join([self._print(j) for j in i])+"]" for i in expr.tolist()])+"]"
        return "%s(%s)" % (
            self._module_format("torch.FloatTensor"),
            str(data)
        )

    def _print_CodegenArrayTensorProduct(self, expr):
        array_list = [j for i, arg in enumerate(expr.args) for j in
                (self._print(arg), "[%i, %i]" % (2*i, 2*i+1))]
        letters = self._get_letter_generator_for_einsum()
        contraction_string = ",".join(["".join([next(letters) for j in range(i)]) for i in expr.subranks])
        return '%s("%s", [%s])' % (
                self._module_format('torch.einsum'),
                contraction_string,
                ", ".join([self._print(arg) for arg in expr.args])
        )

    def _print_CodegenArrayContraction(self, expr):
        from sympy.codegen.array_utils import CodegenArrayTensorProduct
        base = expr.expr
        contraction_indices = expr.contraction_indices
        contraction_string, letters_free, letters_dum = self._get_einsum_string(base.subranks, contraction_indices)

        if len(contraction_indices) == 0:
            return self._print(base)
        if isinstance(base, CodegenArrayTensorProduct):
            elems = ["%s" % (self._print(arg)) for arg in base.args]
            return "%s(\"%s\", [%s])" % (
                self._module_format("torch.einsum"),
                contraction_string,
                ", ".join(elems)
            )
        raise NotImplementedError()

    def _print_CodegenArrayDiagonal(self, expr):
        from sympy.codegen.array_utils import CodegenArrayTensorProduct
        diagonal_indices = list(expr.diagonal_indices)
        if len(diagonal_indices) > 1:
            # TODO: this should be handled in sympy.codegen.array_utils,
            # possibly by creating the possibility of unfolding the
            # CodegenArrayDiagonal object into nested ones. Same reasoning for
            # the array contraction.
            raise NotImplementedError
        if len(diagonal_indices[0]) != 2:
            raise NotImplementedError
        if isinstance(expr.expr, CodegenArrayTensorProduct):
            subranks = expr.expr.subranks
            elems = expr.expr.args
        else:
            subranks = expr.subranks
            elems = [expr.expr]
        diagonal_string, letters_free, letters_dum = self._get_einsum_string(subranks, diagonal_indices)
        elems = [self._print(i) for i in elems]
        return '%s("%s", [%s])' % (
            self._module_format("torch.einsum"),
            "{0}->{1}{2}".format(diagonal_string, "".join(letters_free), "".join(letters_dum)),
            ", ".join(elems)
        )

    def _print_CodegenArrayPermuteDims(self, expr):
        return "%s.permute(%s)" % (
            self._print(expr.expr),
            ", ".join([self._print(i) for i in expr.permutation.args[0]]),
        )

    def _print_CodegenArrayElementwiseAdd(self, expr):
        return self._expand_fold_binary_op('torch.add', expr.args)


def torch_code(expr):
    printer = TorchPrinter()
    return printer.doprint(expr)
