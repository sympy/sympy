from sympy.core import S, Add
from sympy.printing.str import StrPrinter
from sympy.tensor import Idx, Indexed, get_indices, get_contraction_structure

class CodePrinter(StrPrinter):

    def _doprint_a_piece(self, expr, assign_to=None):
        # Here we print an expression that may contain Indexed objects, they
        # correspond to arrays in the generated code.  The low-level implementation
        # involves looping over array elements and possibly storing results in temporary
        # variables or accumulate it in the assign_to object.

        lhs_printed = self._print(assign_to)
        lines = []

        # Setup loops over non-dummy indices  --  all terms need these
        indices = self.get_expression_indices(expr, assign_to)
        openloop, closeloop, junk = self._get_loop_opening_ending_ints(indices)

        # Setup loops over dummy indices  --  each term needs separate treatment
        d = get_contraction_structure(expr)

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
                openloop_d, closeloop_d, junk = self._get_loop_opening_ending_ints(indices)

                for term in d[dummies]:
                    if term in d and not ([f.keys() for f in d[term]]
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
                        if term.has(assign_to):
                            raise(ValueError("FIXME: lhs present in rhs,\
                                this is undefined in CCodePrinter"))

                        lines.extend(openloop)
                        lines.extend(openloop_d)
                        text = "%s = %s" % (lhs_printed, CodePrinter.doprint(self, assign_to + term))
                        lines.append(self._get_statement(text))
                        lines.extend(closeloop_d)
                        lines.extend(closeloop)

        return lines

    def get_expression_indices(self, expr, assign_to):
        rinds, junk = get_indices(expr)
        linds, junk = get_indices(assign_to)

        # support broadcast of scalar
        if linds and not rinds:
            rinds = linds
        if rinds != linds:
            raise ValueError("lhs indices must match non-dummy rhs indices")

        return self._sort_optimized(rinds, assign_to)

    def _sort_optimized(self, indices, expr):

        if not indices:
            return []

        # determine optimized loop order by giving a score to each index
        # the index with the highest score are put in the innermost loop.
        score_table = {}
        for i in indices:
            score_table[i] = 0

        arrays = expr.atoms(Indexed)
        for arr in arrays:
            for p, ind in enumerate(arr.indices):
                try:
                    score_table[ind] += self._rate_index_position(p)
                except KeyError:
                    pass

        return sorted(indices, key=lambda x: score_table[x])

