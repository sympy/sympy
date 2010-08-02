"""
C code printer

The CCodePrinter converts single sympy expressions into single C expressions,
using the functions defined in math.h where possible.

A complete code generator, which uses ccode extensively, can be found in
sympy.utilities.codegen. The codegen module can be used to generate complete
source code files that are compilable without further modifications.
"""

from str import StrPrinter
from sympy.printing.precedence import precedence
from sympy.core import S, Basic, Add, Symbol, NumberSymbol
from sympy.functions import Piecewise, piecewise_fold
from sympy.tensor import Idx, Indexed
from sympy.tensor.index_methods import get_indices, get_contraction_structure


# dictionary mapping sympy function to (argument_conditions, C_function).
# Used in CCodePrinter._print_Function(self)
known_functions = {
        "ceiling": [(lambda x: True, "ceil")],
        "abs": [(lambda x: not x.is_integer, "fabs")],
        }

class CCodePrinter(StrPrinter):
    """A printer to convert python expressions to strings of c code"""
    printmethod = "_ccode"

    _default_settings = {
        'order': None,
        'full_prec': 'auto',
        'precision': 15,
        'user_functions': {},
        'human': True,
    }

    def __init__(self, settings={}):
        """Register function mappings supplied by user"""
        StrPrinter.__init__(self, settings)
        self.known_functions = dict(known_functions)
        userfuncs = settings.get('user_functions', {})
        for k,v in userfuncs.items():
            if not isinstance(v, tuple):
                userfuncs[k] = (lambda *x: True, v)
        self.known_functions.update(userfuncs)

    def doprint(self, expr, assign_to=None):

        if isinstance(assign_to, basestring):
            assign_to = Symbol(assign_to)
        elif not isinstance(assign_to, (Basic, type(None))):
            raise TypeError("CCodePrinter cannot assign to object of type %s"%
                    type(result_variable))

        # keep a set of expressions that are not strictly translatable to C
        # and number constants that must be declared and initialized
        not_c = self._not_c = set()
        self._number_symbols = set()

        # We treat top level Piecewise here to get if tests outside loops
        lines = []
        if isinstance(expr, Piecewise):
            for i, (e, c) in enumerate(expr.args):
                if i == 0:
                    lines.append("if (%s) {" % self._print(c))
                elif i == len(expr.args)-1 and c == True:
                    lines.append("else {")
                else:
                    lines.append("else if (%s) {" % self._print(c))
                code0 = self._doprint_a_piece(e, assign_to)
                lines.extend(code0)
                lines.append("}")
        else:
            code0 = self._doprint_a_piece(expr, assign_to)
            lines.extend(code0)

        # format the output
        if self._settings["human"]:
            frontlines = []
            if len(not_c) > 0:
                frontlines.append("// Not C:")
                for expr in sorted(not_c, key=str):
                    frontlines.append("// %s" % expr)
            for name, value in sorted(self._number_symbols, key=str):
                frontlines.append("double const %s = %s;" % (name, value))
            lines = frontlines + lines
            lines = "\n".join(lines)
            result = self.indent_code(lines)
        else:
            lines = self.indent_code("\n".join(lines))
            result = self._number_symbols, not_c, lines
        del self._not_c
        del self._number_symbols
        return result

    def _doprint_a_piece(self, expr, assign_to=None):
        # Here we print an expression that may contain Indexed objects, they
        # correspond to arrays in the generated code.  The low-level implementation
        # involves looping over array elements and possibly storing results in temporary
        # variables or accumulate it in the assign_to object.

        rc, rnc = get_indices(expr)
        lc, lnc = get_indices(assign_to)

        # support broadcast of scalar
        if lc + lnc and not rc + rnc:
            rc = lc
            rnc = lnc

        if rc + rnc != lc + lnc:
            raise ValueError("lhs indices must match rhs indices")

        # Setup loops over non-dummy indices  --  all terms need these
        openloop, closeloop, junk = self._get_loop_opening_ending_ints(rc + rnc)

        lhs_printed = self._print(assign_to)
        lines = []

        # Setup loops over dummy indices  --  each term needs separate treatment
        d = get_contraction_structure(expr)

        # terms with no summations first
        if None in d:
            text = StrPrinter.doprint(self, Add(*d[None]))
        else:
            # If all terms have summations we must initialize array to Zero
            text = StrPrinter.doprint(self, S.Zero)

        lines.extend(openloop)
        if assign_to is not None:
            text = "%s = %s;" % (lhs_printed, text)
        lines.append(text)
        lines.extend(closeloop)

        for dummies in d:
            # then terms with summations
            if isinstance(dummies, tuple):
                openloop_d, closeloop_d, junk = self._get_loop_opening_ending_ints(dummies)

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
                        text = "%s = %s;" % (lhs_printed, StrPrinter.doprint(self, assign_to + term))
                        lines.append(text)
                        lines.extend(closeloop_d)
                        lines.extend(closeloop)

        return lines

    def _get_loop_opening_ending_ints(self, indices):
        """Returns a tuple (open_lines, close_lines) containing lists of codelines
        """
        # FIXME: sort indices in an optimized way
        open_lines = []
        close_lines = []
        local_ints = []

        loopstart = "for (int %(var)s=%(start)s; %(var)s<%(end)s; %(var)s++){"
        for i in indices:
            # C arrays start at 0 and end at dimension-1
            open_lines.append(loopstart % {
                'var': i.label,
                'start': i.lower,
                'end': i.upper + 1})
            close_lines.append("}")
            local_ints.append(i)
        return open_lines, close_lines, local_ints

    def _print_Pow(self, expr):
        PREC = precedence(expr)
        if expr.exp is S.NegativeOne:
            return '1.0/%s'%(self.parenthesize(expr.base, PREC))
        elif expr.exp == 0.5:
            return 'sqrt(%s)' % self._print(expr.base)
        else:
            return 'pow(%s, %s)'%(self._print(expr.base),
                                 self._print(expr.exp))

    def _print_Rational(self, expr):
        p, q = int(expr.p), int(expr.q)
        return '%d.0/%d.0' % (p, q)

    def _print_Indexed(self, expr):
        # calculate index for 1d array
        dims = expr.dimensions
        inds = [ i.label for i in expr.indices ]
        elem = S.Zero
        offset = S.One
        for i in reversed(range(expr.rank)):
            elem += offset*inds[i]
            offset *= dims[i]
        return "%s[%s]" % (self._print(expr.stem.label), self._print(elem))

    def _print_Exp1(self, expr):
        return "M_E"

    def _print_Pi(self, expr):
        return 'M_PI'

    def _print_Infinity(self, expr):
        return 'HUGE_VAL'

    def _print_NegativeInfinity(self, expr):
        return '-HUGE_VAL'

    def _print_Piecewise(self, expr):
        # This method is called only for inline if constructs
        # Top level piecewise is handled in doprint()
        ecpairs = ["(%s) {\n%s\n}\n" % (self._print(c), self._print(e)) \
                       for e, c in expr.args[:-1]]
        last_line = ""
        if expr.args[-1].cond == True:
            last_line = "else {\n%s\n}" % self._print(expr.args[-1].expr)
        else:
            ecpairs.append("(%s) {\n%s\n" % \
                           (self._print(expr.args[-1].cond),
                            self._print(expr.args[-1].expr)))
        code = "if %s" + last_line
        return code % "else if ".join(ecpairs)

    def _print_And(self, expr):
        PREC = precedence(expr)
        return '&&'.join(self.parenthesize(a, PREC) for a in expr.args)

    def _print_Or(self, expr):
        PREC = precedence(expr)
        return '||'.join(self.parenthesize(a, PREC) for a in expr.args)

    def _print_Not(self, expr):
        PREC = precedence(expr)
        return '!'+self.parenthesize(expr.args[0], PREC)

    def _print_Function(self, expr):
        if expr.func.__name__ in self.known_functions:
            cond_cfunc = self.known_functions[expr.func.__name__]
            for cond, cfunc in cond_cfunc:
                if cond(*expr.args):
                    return "%s(%s)" % (cfunc, self.stringify(expr.args, ", "))
        return StrPrinter._print_Function(self, expr)

    def _print_NumberSymbol(self, expr):
        # A Number symbol that is not implemented here or with _printmethod
        # is registered and printed with str().
        self._number_symbols.add((str(expr), expr.evalf(self._settings["precision"])))
        return str(expr)

    _print_Catalan = _print_NumberSymbol
    _print_EulerGamma = _print_NumberSymbol
    _print_GoldenRatio = _print_NumberSymbol

    def _print_not_c(self, expr):
        self._not_c.add(expr)
        return self.emptyPrinter(expr)

    # The following can not be simply translated into C.
    _print_Basic = _print_not_c
    _print_ComplexInfinity = _print_not_c
    _print_Derivative = _print_not_c
    _print_dict = _print_not_c
    _print_Dummy = _print_not_c
    _print_ExprCondPair = _print_not_c
    _print_GeometryEntity = _print_not_c
    _print_Integral = _print_not_c
    _print_Interval = _print_not_c
    _print_Limit = _print_not_c
    _print_list = _print_not_c
    _print_Matrix = _print_not_c
    _print_DeferredVector = _print_not_c
    _print_NaN = _print_not_c
    _print_Normal = _print_not_c
    _print_Order = _print_not_c
    _print_PDF = _print_not_c
    _print_RootOf = _print_not_c
    _print_RootsOf = _print_not_c
    _print_RootSum = _print_not_c
    _print_Sample = _print_not_c
    _print_SMatrix = _print_not_c
    _print_tuple = _print_not_c
    _print_Uniform = _print_not_c
    _print_Unit = _print_not_c
    _print_Wild = _print_not_c
    _print_WildFunction = _print_not_c

    def indent_code(self, code):
        """Accepts a string of code or a list of code lines"""

        if isinstance(code, basestring):
           code_lines = self.indent_code(code.splitlines())
           return '\n'.join(code_lines)

        tab = "   "
        inc_token = ('{', '(', '{\n', '(\n')
        dec_token = ('}', ')')

        code = [ line.lstrip() for line in code ]

        from sympy.utilities.iterables import any  # 2.4 support
        increase = [ int(any(map(line.endswith, inc_token))) for line in code ]
        decrease = [ int(any(map(line.startswith, dec_token))) for line in code ]

        pretty = []
        level = 0
        for n, line in enumerate(code):
            if not line:
                pretty.append(line)
                continue
            level -= decrease[n]
            pretty.append("%s%s" % (tab*level, line))
            level += increase[n]
        return pretty


def ccode(expr, assign_to=None, **settings):
    r"""Converts an expr to a string of c code

        Arguments:
          expr  --  a sympy expression to be converted

        Optional arguments:
          precision  --  the precision for numbers such as pi [default=15]
          user_functions  --  A dictionary where keys are FunctionClass instances
                              and values are there string representations.
                              Alternatively, the dictionary value can be a list
                              of tuples i.e. [(argument_test, cfunction_string)].
                              See below for examples.
          human  --  If True, the result is a single string that may contain
                     some constant declarations for the number symbols. If
                     False, the same information is returned in a more
                     programmer-friendly data structure.

        >>> from sympy import ccode, symbols, Rational, sin
        >>> x, tau = symbols(["x", "tau"])
        >>> ccode((2*tau)**Rational(7,2))
        '8*sqrt(2)*pow(tau, 7.0/2.0)'
        >>> ccode(sin(x), assign_to="s")
        's = sin(x);'


    """
    return CCodePrinter(settings).doprint(expr, assign_to)

def print_ccode(expr, **settings):
    """Prints C representation of the given expression."""
    print ccode(expr, **settings)
