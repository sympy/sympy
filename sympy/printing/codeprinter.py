from __future__ import print_function, division

from sympy.core import C, Add, Mul, Pow, S
from sympy.core.compatibility import default_sort_key
from sympy.core.mul import _keep_coeff
from sympy.printing.str import StrPrinter
from sympy.printing.precedence import precedence


class AssignmentError(Exception):
    """
    Raised if an assignment variable for a loop is missing.
    """
    pass


class CodePrinter(StrPrinter):
    """
    The base class for code-printing subclasses.
    """

    _operators = {
        'and': '&&',
        'or': '||',
        'not': '!',
    }

    def _doprint_a_piece(self, expr, assign_to=None):
        # Here we print an expression that may contain Indexed objects, they
        # correspond to arrays in the generated code.  The low-level implementation
        # involves looping over array elements and possibly storing results in temporary
        # variables or accumulate it in the assign_to object.

        lhs_printed = self._print(assign_to)
        lines = []

        # Setup loops over non-dummy indices  --  all terms need these
        if self._settings.get('contract', True):
            indices = self.get_expression_indices(expr, assign_to)
        else:
            indices = []
        openloop, closeloop = self._get_loop_opening_ending(indices)

        # Setup loops over dummy indices  --  each term needs separate treatment
        from sympy.tensor import get_contraction_structure
        if self._settings.get('contract', True):
            d = get_contraction_structure(expr)
        else:
            d = {None: (expr,)}

        # terms with no summations first
        if None in d:
            text = CodePrinter.doprint(self, Add(*d[None]))
        else:
            # If all terms have summations we must initialize array to Zero
            text = CodePrinter.doprint(self, 0)

        # skip redundant assignments
        if text != lhs_printed:
            lines.extend(openloop)
            if assign_to is not None:
                text = self._get_statement("%s = %s" % (lhs_printed, text))
            lines.append(text)
            lines.extend(closeloop)

        for dummies in d:
            # then terms with summations
            if isinstance(dummies, tuple):
                indices = self._sort_optimized(dummies, expr)
                openloop_d, closeloop_d = self._get_loop_opening_ending(
                    indices)

                for term in d[dummies]:
                    if term in d and not ([list(f.keys()) for f in d[term]]
                            == [[None] for f in d[term]]):
                        # If one factor in the term has it's own internal
                        # contractions, those must be computed first.
                        # (temporary variables?)
                        raise NotImplementedError(
                            "FIXME: no support for contractions in factor yet")
                    else:

                        # We need the lhs expression as an accumulator for
                        # the loops, i.e
                        #
                        # for (int d=0; d < dim; d++){
                        #    lhs[] = lhs[] + term[][d]
                        # }           ^.................. the accumulator
                        #
                        # We check if the expression already contains the
                        # lhs, and raise an exception if it does, as that
                        # syntax is currently undefined.  FIXME: What would be
                        # a good interpretation?
                        if assign_to is None:
                            raise AssignmentError(
                                "need assignment variable for loops")
                        if term.has(assign_to):
                            raise ValueError("FIXME: lhs present in rhs,\
                                this is undefined in CCodePrinter")

                        lines.extend(openloop)
                        lines.extend(openloop_d)
                        text = "%s = %s" % (lhs_printed, CodePrinter.doprint(
                            self, assign_to + term))
                        lines.append(self._get_statement(text))
                        lines.extend(closeloop_d)
                        lines.extend(closeloop)

        return lines

    def get_expression_indices(self, expr, assign_to):
        from sympy.tensor import get_indices, get_contraction_structure
        rinds, junk = get_indices(expr)
        linds, junk = get_indices(assign_to)

        # support broadcast of scalar
        if linds and not rinds:
            rinds = linds
        if rinds != linds:
            raise ValueError("lhs indices must match non-dummy"
                    " rhs indices in %s" % expr)

        return self._sort_optimized(rinds, assign_to)

    def _sort_optimized(self, indices, expr):

        if not indices:
            return []

        # determine optimized loop order by giving a score to each index
        # the index with the highest score are put in the innermost loop.
        score_table = {}
        for i in indices:
            score_table[i] = 0

        arrays = expr.atoms(C.Indexed)
        for arr in arrays:
            for p, ind in enumerate(arr.indices):
                try:
                    score_table[ind] += self._rate_index_position(p)
                except KeyError:
                    pass

        return sorted(indices, key=lambda x: score_table[x])

    def _print_NumberSymbol(self, expr):
        # A Number symbol that is not implemented here or with _printmethod
        # is registered and evaluated
        self._number_symbols.add((expr,
            self._print(expr.evalf(self._settings["precision"]))))
        return str(expr)

    def _print_Dummy(self, expr):
        # dummies must be printed as unique symbols
        return "%s_%i" % (expr.name, expr.dummy_index)  # Dummy

    _print_Catalan = _print_NumberSymbol
    _print_EulerGamma = _print_NumberSymbol
    _print_GoldenRatio = _print_NumberSymbol

    def _print_And(self, expr):
        PREC = precedence(expr)
        return (" %s " % self._operators['and']).join(self.parenthesize(a, PREC)
                for a in sorted(expr.args, key=default_sort_key))

    def _print_Or(self, expr):
        PREC = precedence(expr)
        return (" %s " % self._operators['or']).join(self.parenthesize(a, PREC)
                for a in sorted(expr.args, key=default_sort_key))

    def _print_Xor(self, expr):
        if self._operators.get('xor') is None:
            return self._print_not_supported(expr)
        PREC = precedence(expr)
        return (" %s " % self._operators['xor']).join(self.parenthesize(a, PREC)
                for a in expr.args)

    def _print_Equivalent(self, expr):
        if self._operators.get('equivalent') is None:
            return self._print_not_supported(expr)
        PREC = precedence(expr)
        return (" %s " % self._operators['equivalent']).join(self.parenthesize(a, PREC)
                for a in expr.args)

    def _print_Not(self, expr):
        PREC = precedence(expr)
        return self._operators['not'] + self.parenthesize(expr.args[0], PREC)

    def _print_Mul(self, expr):

        prec = precedence(expr)

        c, e = expr.as_coeff_Mul()
        if c < 0:
            expr = _keep_coeff(-c, e)
            sign = "-"
        else:
            sign = ""

        a = []  # items in the numerator
        b = []  # items that are in the denominator (if any)

        if self.order not in ('old', 'none'):
            args = expr.as_ordered_factors()
        else:
            # use make_args in case expr was something like -x -> x
            args = Mul.make_args(expr)

        # Gather args for numerator/denominator
        for item in args:
            if item.is_commutative and item.is_Pow and item.exp.is_Rational and item.exp.is_negative:
                if item.exp != -1:
                    b.append(Pow(item.base, -item.exp, evaluate=False))
                else:
                    b.append(Pow(item.base, -item.exp))
            else:
                a.append(item)

        a = a or [S.One]

        a_str = [self.parenthesize(x, prec) for x in a]
        b_str = [self.parenthesize(x, prec) for x in b]

        if len(b) == 0:
            return sign + '*'.join(a_str)
        elif len(b) == 1:
            if len(a) == 1 and not (a[0].is_Atom or a[0].is_Add):
                return sign + "%s/" % a_str[0] + '*'.join(b_str)
            else:
                return sign + '*'.join(a_str) + "/%s" % b_str[0]
        else:
            return sign + '*'.join(a_str) + "/(%s)" % '*'.join(b_str)

    def _print_not_supported(self, expr):
        self._not_supported.add(expr)
        return self.emptyPrinter(expr)

    # The following can not be simply translated into C or Fortran
    _print_Basic = _print_not_supported
    _print_ComplexInfinity = _print_not_supported
    _print_Derivative = _print_not_supported
    _print_dict = _print_not_supported
    _print_ExprCondPair = _print_not_supported
    _print_GeometryEntity = _print_not_supported
    _print_Infinity = _print_not_supported
    _print_Integral = _print_not_supported
    _print_Interval = _print_not_supported
    _print_Limit = _print_not_supported
    _print_list = _print_not_supported
    _print_Matrix = _print_not_supported
    _print_DeferredVector = _print_not_supported
    _print_NaN = _print_not_supported
    _print_NegativeInfinity = _print_not_supported
    _print_Normal = _print_not_supported
    _print_Order = _print_not_supported
    _print_PDF = _print_not_supported
    _print_RootOf = _print_not_supported
    _print_RootsOf = _print_not_supported
    _print_RootSum = _print_not_supported
    _print_Sample = _print_not_supported
    _print_SparseMatrix = _print_not_supported
    _print_tuple = _print_not_supported
    _print_Uniform = _print_not_supported
    _print_Unit = _print_not_supported
    _print_Wild = _print_not_supported
    _print_WildFunction = _print_not_supported
