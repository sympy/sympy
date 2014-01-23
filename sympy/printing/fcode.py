"""
Fortran code printer

The FCodePrinter converts single sympy expressions into single Fortran
expressions, using the functions defined in the Fortran 77 standard where
possible. Some useful pointers to Fortran can be found on wikipedia:

http://en.wikipedia.org/wiki/Fortran

Most of the code below is based on the "Professional Programmer\'s Guide to
Fortran77" by Clive G. Page:

http://www.star.le.ac.uk/~cgp/prof77.html

Fortran is a case-insensitive language. This might cause trouble because
SymPy is case sensitive. The implementation below does not care and leaves
the responsibility for generating properly cased Fortran code to the user.
"""

from __future__ import print_function, division

import string

from sympy.core import S, C, Add, N
from sympy.core.compatibility import string_types
from sympy.printing.codeprinter import CodePrinter
from sympy.printing.precedence import precedence

class FCodePrinter(CodePrinter):
    """A printer to convert sympy expressions to strings of Fortran code"""
    printmethod = "_fcode"

    _default_settings = {
        'order': None,
        'full_prec': 'auto',
        'assign_to': None,
        'precision': 15,
        'user_functions': {},
        'human': True,
        'source_format': 'fixed',
        'contract': True,
    }

    _implicit_functions = set([
        "sin", "cos", "tan", "asin", "acos", "atan", "atan2", "sinh",
        "cosh", "tanh", "sqrt", "log", "exp", "erf", "Abs", "sign", "conjugate",
    ])

    _operators = {
        'and': '.and.',
        'or': '.or.',
        'xor': '.neqv.',
        'equivalent': '.eqv.',
        'not': '.not. ',
    }

    _relationals = {
        '!=': '/=',
    }

    def __init__(self, settings=None):
        CodePrinter.__init__(self, settings)
        self._init_leading_padding()
        assign_to = self._settings['assign_to']
        if isinstance(assign_to, string_types):
            self._settings['assign_to'] = C.Symbol(assign_to)
        elif not isinstance(assign_to, (C.Basic, type(None))):
            raise TypeError("FCodePrinter cannot assign to object of type %s" %
                    type(assign_to))

    def _rate_index_position(self, p):
        """function to calculate score based on position among indices

        This method is used to sort loops in an optimized order, see
        CodePrinter._sort_optimized()
        """
        return -p*5

    def _get_statement(self, codestring):
        return codestring

    def _init_leading_padding(self):
        # leading columns depend on fixed or free format
        if self._settings['source_format'] == 'fixed':
            self._lead_code = "      "
            self._lead_cont = "     @ "
            self._lead_comment = "C     "
        elif self._settings['source_format'] == 'free':
            self._lead_code = ""
            self._lead_cont = "      "
            self._lead_comment = "! "
        else:
            raise ValueError(
                "Unknown source format: %s" % self._settings[
                'source_format']
            )

    def _pad_leading_columns(self, lines):
        result = []
        for line in lines:
            if line.startswith('!'):
                result.append(self._lead_comment + line[1:].lstrip())
            else:
                result.append(self._lead_code + line)
        return result

    def _get_loop_opening_ending(self, indices):
        """Returns a tuple (open_lines, close_lines) containing lists of codelines
        """
        open_lines = []
        close_lines = []
        for i in indices:
            # fortran arrays start at 1 and end at dimension
            var, start, stop = map(self._print,
                    [i.label, i.lower + 1, i.upper + 1])
            open_lines.append("do %s = %s, %s" % (var, start, stop))
            close_lines.append("end do")
        return open_lines, close_lines

    def doprint(self, expr):
        """Returns Fortran code for expr (as a string)"""
        # find all number symbols
        self._number_symbols = set()

        # keep a set of expressions that are not strictly translatable to
        # Fortran.
        self._not_supported = set()

        lines = []
        from sympy.functions import Piecewise
        if isinstance(expr, Piecewise):
            # support for top-level Piecewise function
            for i, (e, c) in enumerate(expr.args):
                if i == 0:
                    lines.append("if (%s) then" % self._print(c))
                elif i == len(expr.args) - 1 and c == True:
                    lines.append("else")
                else:
                    lines.append("else if (%s) then" % self._print(c))
                lines.extend(
                    self._doprint_a_piece(e, self._settings['assign_to']))
            lines.append("end if")
        else:
            lines.extend(
                self._doprint_a_piece(expr, self._settings['assign_to']))

        # format the output
        if self._settings["human"]:
            frontlines = []
            if len(self._not_supported) > 0:
                frontlines.append("! Not Fortran:")
                for expr in sorted(self._not_supported, key=self._print):
                    frontlines.append("! %s" % repr(expr))
            for name, value in sorted(self._number_symbols, key=str):
                frontlines.append("parameter (%s = %s)" % (str(name), value))
            frontlines.extend(lines)
            lines = frontlines
            lines = self.indent_code(lines)
            lines = self._wrap_fortran(lines)
            result = "\n".join(lines)
        else:
            lines = self.indent_code(lines)
            lines = self._wrap_fortran(lines)
            result = self._number_symbols, self._not_supported, "\n".join(
                lines)

        del self._not_supported
        del self._number_symbols
        return result

    def _print_Add(self, expr):
        # purpose: print complex numbers nicely in Fortran.
        # collect the purely real and purely imaginary parts:
        pure_real = []
        pure_imaginary = []
        mixed = []
        for arg in expr.args:
            if arg.is_number and arg.is_real:
                pure_real.append(arg)
            elif arg.is_number and arg.is_imaginary:
                pure_imaginary.append(arg)
            else:
                mixed.append(arg)
        if len(pure_imaginary) > 0:
            if len(mixed) > 0:
                PREC = precedence(expr)
                term = Add(*mixed)
                t = self._print(term)
                if t.startswith('-'):
                    sign = "-"
                    t = t[1:]
                else:
                    sign = "+"
                if precedence(term) < PREC:
                    t = "(%s)" % t

                return "cmplx(%s,%s) %s %s" % (
                    self._print(Add(*pure_real)),
                    self._print(-S.ImaginaryUnit*Add(*pure_imaginary)),
                    sign, t,
                )
            else:
                return "cmplx(%s,%s)" % (
                    self._print(Add(*pure_real)),
                    self._print(-S.ImaginaryUnit*Add(*pure_imaginary)),
                )
        else:
            return CodePrinter._print_Add(self, expr)

    def _print_Function(self, expr):
        name = self._settings["user_functions"].get(expr.__class__)
        eargs = expr.args
        if name is None:
            from sympy.functions import conjugate
            if expr.func == conjugate:
                name = "conjg"
            else:
                name = expr.func.__name__
            if hasattr(expr, '_imp_') and isinstance(expr._imp_, C.Lambda):
                # inlined function.
                # the expression is printed with _print to avoid loops
                return self._print(expr._imp_(*eargs))
            if expr.func.__name__ not in self._implicit_functions:
                self._not_supported.add(expr)
            else:
                # convert all args to floats
                eargs = map(N, eargs)
        return "%s(%s)" % (name, self.stringify(eargs, ", "))

    _print_factorial = _print_Function

    def _print_ImaginaryUnit(self, expr):
        # purpose: print complex numbers nicely in Fortran.
        return "cmplx(0,1)"

    def _print_int(self, expr):
        return str(expr)

    def _print_Mul(self, expr):
        # purpose: print complex numbers nicely in Fortran.
        if expr.is_number and expr.is_imaginary:
            return "cmplx(0,%s)" % (
                self._print(-S.ImaginaryUnit*expr)
            )
        else:
            return CodePrinter._print_Mul(self, expr)

    _print_Exp1 = CodePrinter._print_NumberSymbol
    _print_Pi = CodePrinter._print_NumberSymbol

    def _print_Pow(self, expr):
        PREC = precedence(expr)
        if expr.exp == -1:
            return '1.0/%s' % (self.parenthesize(expr.base, PREC))
        elif expr.exp == 0.5:
            if expr.base.is_integer:
                # Fortan intrinsic sqrt() does not accept integer argument
                if expr.base.is_Number:
                    return 'sqrt(%s.0d0)' % self._print(expr.base)
                else:
                    return 'sqrt(dble(%s))' % self._print(expr.base)
            else:
                return 'sqrt(%s)' % self._print(expr.base)
        else:
            return CodePrinter._print_Pow(self, expr)

    def _print_Rational(self, expr):
        p, q = int(expr.p), int(expr.q)
        return "%d.0d0/%d.0d0" % (p, q)

    def _print_Float(self, expr):
        printed = CodePrinter._print_Float(self, expr)
        e = printed.find('e')
        if e > -1:
            return "%sd%s" % (printed[:e], printed[e + 1:])
        return "%sd0" % printed

    def _print_Indexed(self, expr):
        inds = [ self._print(i) for i in expr.indices ]
        return "%s(%s)" % (self._print(expr.base.label), ", ".join(inds))

    def _print_Idx(self, expr):
        return self._print(expr.label)

    def _wrap_fortran(self, lines):
        """Wrap long Fortran lines

           Argument:
             lines  --  a list of lines (without \\n character)

           A comment line is split at white space. Code lines are split with a more
           complex rule to give nice results.
        """
        # routine to find split point in a code line
        my_alnum = set("_+-." + string.digits + string.ascii_letters)
        my_white = set(" \t()")

        def split_pos_code(line, endpos):
            if len(line) <= endpos:
                return len(line)
            pos = endpos
            split = lambda pos: \
                (line[pos] in my_alnum and line[pos - 1] not in my_alnum) or \
                (line[pos] not in my_alnum and line[pos - 1] in my_alnum) or \
                (line[pos] in my_white and line[pos - 1] not in my_white) or \
                (line[pos] not in my_white and line[pos - 1] in my_white)
            while not split(pos):
                pos -= 1
                if pos == 0:
                    return endpos
            return pos
        # split line by line and add the splitted lines to result
        result = []
        if self._settings['source_format'] == 'free':
            trailing = ' &'
        else:
            trailing = ''
        for line in lines:
            if line.startswith(self._lead_comment):
                # comment line
                if len(line) > 72:
                    pos = line.rfind(" ", 6, 72)
                    if pos == -1:
                        pos = 72
                    hunk = line[:pos]
                    line = line[pos:].lstrip()
                    result.append(hunk)
                    while len(line) > 0:
                        pos = line.rfind(" ", 0, 66)
                        if pos == -1 or len(line) < 66:
                            pos = 66
                        hunk = line[:pos]
                        line = line[pos:].lstrip()
                        result.append("%s%s" % (self._lead_comment, hunk))
                else:
                    result.append(line)
            elif line.startswith(self._lead_code):
                # code line
                pos = split_pos_code(line, 72)
                hunk = line[:pos].rstrip()
                line = line[pos:].lstrip()
                if line:
                    hunk += trailing
                result.append(hunk)
                while len(line) > 0:
                    pos = split_pos_code(line, 65)
                    hunk = line[:pos].rstrip()
                    line = line[pos:].lstrip()
                    if line:
                        hunk += trailing
                    result.append("%s%s" % (self._lead_cont, hunk))
            else:
                result.append(line)
        return result

    def indent_code(self, code):
        """Accepts a string of code or a list of code lines"""
        if isinstance(code, string_types):
            code_lines = self.indent_code(code.splitlines(True))
            return ''.join(code_lines)

        free = self._settings['source_format'] == 'free'
        code = [ line.lstrip(' \t') for line in code ]

        inc_keyword = ('do ', 'if(', 'if ', 'do\n', 'else')
        dec_keyword = ('end do', 'enddo', 'end if', 'endif', 'else')

        increase = [ int(any(map(line.startswith, inc_keyword)))
                     for line in code ]
        decrease = [ int(any(map(line.startswith, dec_keyword)))
                     for line in code ]
        continuation = [ int(any(map(line.endswith, ['&', '&\n'])))
                         for line in code ]

        level = 0
        cont_padding = 0
        tabwidth = 3
        new_code = []
        for i, line in enumerate(code):
            if line == '' or line == '\n':
                new_code.append(line)
                continue
            level -= decrease[i]

            if free:
                padding = " "*(level*tabwidth + cont_padding)
            else:
                padding = " "*level*tabwidth

            line = "%s%s" % (padding, line)
            if not free:
                line = self._pad_leading_columns([line])[0]

            new_code.append(line)

            if continuation[i]:
                cont_padding = 2*tabwidth
            else:
                cont_padding = 0
            level += increase[i]

        if not free:
            return self._wrap_fortran(new_code)
        return new_code


def fcode(expr, **settings):
    """Converts an expr to a string of Fortran 77 code

       Parameters
       ==========

       expr : sympy.core.Expr
           a sympy expression to be converted
       assign_to : optional
           When given, the argument is used as the name of the
           variable to which the Fortran expression is assigned.
           (This is helpful in case of line-wrapping.)
       precision : optional
           the precision for numbers such as pi [default=15]
       user_functions : optional
           A dictionary where keys are FunctionClass instances and values
           are there string representations.
       human : optional
           If True, the result is a single string that may contain some
           parameter statements for the number symbols. If False, the same
           information is returned in a more programmer-friendly data
           structure.
       source_format : optional
           The source format can be either 'fixed' or 'free'.
           [default='fixed']
       contract: optional
           If True, `Indexed` instances are assumed to obey
           tensor contraction rules and the corresponding nested
           loops over indices are generated. Setting contract = False
           will not generate loops, instead the user is responsible
           to provide values for the indices in the code. [default=True]

       Examples
       ========

       >>> from sympy import fcode, symbols, Rational, pi, sin
       >>> x, tau = symbols('x,tau')
       >>> fcode((2*tau)**Rational(7,2))
       '      8*sqrt(2.0d0)*tau**(7.0d0/2.0d0)'
       >>> fcode(sin(x), assign_to="s")
       '      s = sin(x)'
       >>> print(fcode(pi))
             parameter (pi = 3.14159265358979d0)
             pi
       >>> from sympy import Eq, IndexedBase, Idx
       >>> len_y = 5
       >>> y = IndexedBase('y', shape=(len_y,))
       >>> t = IndexedBase('t', shape=(len_y,))
       >>> Dy = IndexedBase('Dy', shape=(len_y-1,))
       >>> i = Idx('i', len_y-1)
       >>> e=Eq(Dy[i], (y[i+1]-y[i])/(t[i+1]-t[i]))
       >>> fcode(e.rhs, assign_to=e.lhs, contract=False)
       '      Dy(i) = (y(i + 1) - y(i))*1.0/(t(i + 1) - t(i))'

    """
    # run the printer
    printer = FCodePrinter(settings)
    return printer.doprint(expr)


def print_fcode(expr, **settings):
    """Prints the Fortran representation of the given expression.

       See fcode for the meaning of the optional arguments.
    """
    print(fcode(expr, **settings))
