"""
Octave/Matlab code printer
"""

from __future__ import print_function, division
from sympy.core import C, Add, Mul, Pow, S, Rational
from sympy.core.compatibility import string_types
from sympy.core.mul import _keep_coeff
from sympy.printing.codeprinter import CodePrinter, Assignment
from sympy.printing.str import StrPrinter
from sympy.printing.precedence import precedence
# FIXME: causes errors on import
#from sympy.functions import sqrt

# dictionary mapping sympy function to (argument_conditions, C_function).
# Used in OctaveCodePrinter._print_Function(self)
# FIXME: almost certainly incomplete, perhaps better to default yet!
known_functions = {
    "sin": "sin",
    "cos": "cos",
    "tan": "tan",
    "asin": "asin",
    "acos": "acos",
    "atan": "atan",
    "atan2": "atan2",
    "sinh": "sinh",
    "cosh": "cosh",
    "tanh": "tanh",
    "asinh": "asinh",
    "acosh": "acosh",
    "atanh": "atanh",
    "log": "log",
    "exp": "exp",
    "erf": "erf",
    "gamma": "gamma",
    "sign": "sign",
    "floor": "floor",
    "Abs": "abs",       # different names after this line
    "ceiling": "ceil",
    "conjugate": "conj",
}
# FIXME: test and/or merge these in
# csc
# sec
# cot
# coth
# acot
# acoth
# factorial
# erfc
# erfinv
# erfcinv
# erfi
# heaviside|Heaviside
# dirac|DiracDelta


class OctaveCodePrinter(CodePrinter):
    """
    A printer to convert expressions to strings of Octave/Matlab code.
    """
    printmethod = "_octave"
    language = "Octave"

    _operators = {
        'and': '&&',
        'or': '||',
        'not': '~',
    }

    _default_settings = {
        'order': None,
        'full_prec': 'auto',
        'precision': 15,
        'user_functions': {},
        'human': True,
        'contract': True
    }
    # FIXME: contract is for expressing tensors as loops (if True), or
    # just assignment (if False).  Needs tests for tensors, borrow
    # some from C?

    def __init__(self, settings={}):
        super(OctaveCodePrinter, self).__init__(settings)
        self.known_functions = dict(known_functions)
        userfuncs = settings.get('user_functions', {})
        self.known_functions.update(userfuncs)

    def _rate_index_position(self, p):
        return p*5

    def _get_statement(self, codestring):
        return "%s;" % codestring

    def _get_comment(self, text):
        return "% {0}".format(text)

    def _declare_number_const(self, name, value):
        return "{0} = {1};".format(name, value)

    def _format_code(self, lines):
        return self.indent_code(lines)

    def _traverse_matrix_indices(self, mat):
        # Octave uses Fortran order (column-major)
        # FIXME: need tests for matrix traverse stuff
        rows, cols = mat.shape
        return ((i, j) for j in range(cols) for i in range(rows))

    def _get_loop_opening_ending(self, indices):
        open_lines = []
        close_lines = []
        for i in indices:
            # Octave arrays start at 1 and end at dimension
            var, start, stop = map(self._print,
                    [i.label, i.lower + 1, i.upper + 1])
            open_lines.append("for %s = %s:%s" % (var, start, stop))
            close_lines.append("end")
        return open_lines, close_lines


    def _print_Mul(self, expr):
        # print complex numbers nicely in Octave
        if (expr.is_number and expr.is_imaginary and
                expr.as_coeff_Mul()[0].is_integer):
            return "%si" % self._print(-S.ImaginaryUnit*expr)

        # cribbed from str.py
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
            if (item.is_commutative and item.is_Pow and item.exp.is_Rational
                    and item.exp.is_negative):
                if item.exp != -1:
                    b.append(Pow(item.base, -item.exp, evaluate=False))
                else:
                    b.append(Pow(item.base, -item.exp))
            elif item.is_Rational and item is not S.Infinity:
                if item.p != 1:
                    a.append(Rational(item.p))
                if item.q != 1:
                    b.append(Rational(item.q))
            else:
                a.append(item)

        a = a or [S.One]

        a_str = list(map(lambda x: self.parenthesize(x, prec), a))
        b_str = list(map(lambda x: self.parenthesize(x, prec), b))

        # from here it differs from str.py to deal with "*" and ".*"
        def multjoin(a, a_str):
            # here we probably are assuming the constants will come first
            r = a_str[0]
            #for (ai, ai_str) in zip(a, a_str)[1:]:
            #    if ai.is_number
            for i in range(1, len(a)):
                if a[i-1].is_number:
                    mulsym = '*'
                else:
                    mulsym = '.*'
                r = r + mulsym + a_str[i]
            return r

        if len(b) == 0:
            return sign + multjoin(a, a_str)
        elif len(b) == 1:
            if b[0].is_number:
                divsym = '/'
            else:
                divsym = './'
            return sign + multjoin(a, a_str) + divsym + b_str[0]
        else:
            if all([bi.is_number for bi in b]):
                divsym = '/'
            else:
                divsym = './'
            return (sign + multjoin(a, a_str) +
                    divsym + "(%s)" % multjoin(b, b_str))

    # FIXME: need tests that use HadamardProduct: user needs this to
    # mix matrix symbols and .* products?
    #def _print_MatMul(self, expr):
    #def _print_HadamardProduct(self, expr):


    def _print_Pow(self, expr):
        if all([x.is_number for x in expr.args]):
            sym = '^'
        else:
            sym = '.^'
        return super(OctaveCodePrinter, self)._print_Pow(expr, powsymbol=sym)

    def _print_MatPow(self, expr):
        return super(OctaveCodePrinter, self)._print_MatPow(expr, powsymbol='^')

    def _print_Pi(self, expr):
        return 'pi'

    def _print_ImaginaryUnit(self, expr):
        return "1i"

    def _print_Exp1(self, expr):
        return "exp(1)"

    def _print_GoldenRatio(self, expr):
        # FIXME: how to do better, e.g., for octave_code(2*GoldenRatio)?
        #return self._print((1+sqrt(S(5)))/2)
        return "(1+sqrt(5))/2"

    def _print_NumberSymbol(self, expr):
        # FIXME: perhaps this should be based on Assignment?  Like in
        # Piecewise.  This form is better for inline code, but the
        # CodePrinter implementation is nicer for longer programs.
        return "%.18g" % float(expr)

    def _print_Infinity(self, expr):
        return 'inf'

    def _print_NegativeInfinity(self, expr):
        return '-inf'

    def _print_NaN(self, expr):
        return 'NaN'

    def _print_list(self, expr):
        return '{' + ', '.join(self._print(a) for a in expr) + '}'
    _print_tuple = _print_list
    _print_Tuple = _print_list

    # FIXME: ccode (and fcode?) need better bool support
    def _print_BooleanTrue(self, expr):
        return "true"

    def _print_BooleanFalse(self, expr):
        return "false"

    def _print_bool(self, expr):
        return str(expr).lower()

    def _print_MatrixBase(self, A):
        # Handle zero dimensions:
        if (A.rows, A.cols) == (0, 0):
            return '[]'
        elif A.rows == 0 or A.cols == 0:
            return 'zeros(%s, %s)' % (A.rows, A.cols)
        elif (A.rows, A.cols) == (1, 1):
            # Octave does not distinguish between scalars and 1x1 matrices
            return self._print(A[0, 0])
        elif A.rows == 1:
            return "[%s]" % A.table(self, rowstart='', rowend='', colsep=' ')
        elif A.cols == 1:
            # FIXME: not ideal, makes each equispaced
            return "[%s]" % A.table(self, rowstart='', rowend='',
                                    rowsep='; ', colsep=' ')
        return "[%s]" % A.table(self, rowstart='', rowend='',
                                rowsep='; ...\n', colsep=' ')
    # FIXME: see my prosposed change for _print_NumberSymbol, same here
    _print_SparseMatrix = \
        _print_MutableSparseMatrix = \
        _print_ImmutableSparseMatrix = \
        _print_Matrix = \
        _print_DenseMatrix = \
        _print_MutableDenseMatrix = \
        _print_ImmutableMatrix = \
        _print_ImmutableDenseMatrix = \
        _print_MatrixBase

    def _print_MatrixElement(self, expr):
        return self._print(expr.parent) + '(%s, %s)'%(expr.i+1, expr.j+1)

    def _print_Indexed(self, expr):
        inds = [ self._print(i) for i in expr.indices ]
        return "%s(%s)" % (self._print(expr.base.label), ", ".join(inds))

    def _print_Idx(self, expr):
        return self._print(expr.label)

    def _print_Identity(self, expr):
        return "eye(%s)" % self._print(expr.shape[0])

    def _print_Piecewise(self, expr):
        if expr.args[-1].cond != True:
            # We need the last conditional to be a True, otherwise the resulting
            # function may not return a result.
            raise ValueError("All Piecewise expressions must contain an "
                             "(expr, True) statement to be used as a default "
                             "condition. Without one, the generated "
                             "expression may not evaluate to anything under "
                             "some condition.")
        lines = []
        # FIXME: for Octave, user might want to force the inline mode (its
        # better for vector operations for example).  How to expose this?
        # Choosing just based on Assignment seems too coarse.
        if expr.has(Assignment):
            for i, (e, c) in enumerate(expr.args):
                if i == 0:
                    lines.append("if (%s)" % self._print(c))
                elif i == len(expr.args) - 1 and c == True:
                    lines.append("else")
                else:
                    lines.append("elseif (%s)" % self._print(c))
                code0 = self._print(e)
                lines.append(code0)
                if i == len(expr.args) - 1:
                    lines.append("end")
            return "\n".join(lines)
        else:
            # This Piecewise was used in an expression, so do inline
            # where each cond, expr pair is like a nested Horner form:
            #   (condition) .* (expr) + (not cond) .* (<others>)
            # FIXME: ccode says some things won't work inline, true here?
            ecpairs = ["({0}).*({1}) + (~({0})).*(".format
                       (self._print(c), self._print(e))
                       for e, c in expr.args[:-1]]
            elast = "%s" % self._print(expr.args[-1].expr)
            pw = " ...\n".join(ecpairs) + elast + ")"*len(ecpairs)
            # Note: current need these outer brackets for 2*pw.  Would be
            # nicer to teach parenthesize() to do this for us when needed!
            return "(" + pw + ")"

    def indent_code(self, code):
        """Accepts a string of code or a list of code lines"""

        # code mostly copied from ccode
        if isinstance(code, string_types):
            code_lines = self.indent_code(code.splitlines(True))
            return ''.join(code_lines)

        tab = "  "
        inc_token = ('function ', 'if ', 'elseif ', 'else', 'for')
        dec_token = ('end')

        # pre-strip left-space from the code
        code = [ line.lstrip(' \t') for line in code ]

        increase = [ int(any(map(line.startswith, inc_token)))
                     for line in code ]
        decrease = [ int(any(map(line.startswith, dec_token)))
                     for line in code ]

        pretty = []
        level = 0
        for n, line in enumerate(code):
            if line == '' or line == '\n':
                pretty.append(line)
                continue
            level -= decrease[n]
            pretty.append("%s%s" % (tab*level, line))
            level += increase[n]
        return pretty


def octave_code(expr, assign_to=None, **settings):
    r"""Converts an expr to a string of Octave/Matlab code

    Parameters
    ==========

    expr : Expr
        A sympy expression to be converted.
    assign_to : optional
        When given, the argument is used as the name of the variable to which
        the expression is assigned.  Can be a string, ``Symbol``,
        ``MatrixSymbol``, or ``Indexed`` type.  This can be helpful for
        expressions that generate multi-line statements.
    precision : integer, optional
        The precision for numbers such as pi [default=15].
        [FIXME: is this used?  remove?]
    user_functions : dict, optional
        A dictionary where keys are ``FunctionClass`` instances and values are
        their string representations.  Alternatively, the dictionary value can
        be a list of tuples i.e. [(argument_test, cfunction_string)].  See
        below for examples.
    human : bool, optional
        If True, the result is a single string that may contain some constant
        declarations for the number symbols.  If False, the same information is
        returned in a tuple of (symbols_to_declare, not_supported_functions,
        code_text).  [default=True].
    contract: bool, optional
        If True, ``Indexed`` instances are assumed to obey tensor contraction
        rules and the corresponding nested loops over indices are generated.
        Setting contract=False will not generate loops, instead the user is
        responsible to provide values for the indices in the code.
        [default=True].
        [FIXME: these are untested]

    Examples
    ========

    >>> from sympy import octave_code, symbols, sin
    >>> x = symbols('x')
    >>> octave_code(sin(x).series(x).removeO())
    'x.^5/120 - x.^3/6 + x'

    >>> from sympy import Rational, ceiling, Abs
    >>> x, y, tau = symbols("x, y, tau")
    >>> octave_code((2*tau)**Rational(7, 2))
    '8*sqrt(2)*tau.^(7/2)'

    Note that element-wise (Hadamard) operations are used by default between
    symbols.  This is because its very common in Octave to write "vectorized"
    code.  It is harmless if the values are scalars.

    >>> octave_code(sin(pi*x*y), assign_to="s")
    's = sin(pi*x.*y);'

    If you need a matrix product "*" or matrix power "^", you can specify the
    symbol as a ``MatrixSymbol``.

    >>> from sympy import Symbol, MatrixSymbol
    >>> n = Symbol('n', integer=True, positive=True)
    >>> A = MatrixSymbol('A', n, n)
    >>> octave_code(3*pi*A**3)
    '(3*pi)*A^3

    Unfortunately, there is currently there is no easy way to specify scalar
    symbols (other than 1x1 Matrix), so sometimes the code might have some
    minor cosmetic issues.  For example, here presumably x and y are scalars
    and a human being might write "(x^2*y)*A^3":

    >>> octave_code(x**2*y*A**3)
    '(x.^2.*y)*A^3

    Matrices are supported.  They can be assigned to a string using
    ``assign_to`` or to a ``MatrixSymbol``.  The latter must have the same
    dimensions.
    [FIXME: currently can also be assigned to a ``Symbol``.]

    >>> from sympy import Matrix, MatrixSymbol
    >>> mat = Matrix([[x**2, sin(x)], [x*y, ceiling(x)]])
    >>> print(octave_code(mat, assign_to='A'))
    A = [x.^2  sin(x); ...
    x.*y ceil(x)];

    Contrast this with:

    >>> A = MatrixSymbol('A', 2, 2)
    >>> print(octave_code(mat, A))
    A(1, 1) = x.^2;
    A(2, 1) = x.*y;
    A(1, 2) = sin(x);
    A(2, 2) = ceil(x);

    ``Piecewise`` expressions can be dealt with in two ways using either
    conditionals or logical masking.  Currently, if an ``assign_to`` variable
    is provided then an ``if`` statement is created, otherwise logical masking
    is used.  A future version might offer a more customizable choice [FIXME].
    Note that if the ``Piecewise`` lacks a default term, represented by
    ``(expr, True)`` then an error will be thrown.  This is to prevent
    generating an expression that may not evaluate to anything.

    >>> from sympy import Piecewise
    >>> pw = Piecewise((x + 1, x > 0), (x, True))
    >>> print(octave_code(pw, assign_to=tau))
    if (x > 0)
      tau = x + 1;
    else
      tau = x;
    end
    >>> octave_code(pw)
    '((x > 0).*(x + 1) + (~(x > 0)).*(x))'

    Note that any expression that can be generated normally can also exist
    inside a Matrix:

    >>> mat = Matrix([[x**2, pw, sin(x)]])
    >>> octave_code(mat, assign_to='A')
    'A = [x.^2 ((x > 0).*(x + 1) + (~(x > 0)).*(x)) sin(x)];'
    >>> A = MatrixSymbol('A', 1, 3)
    >>> print(octave_code(pw, assign_to=A))
    A(1, 1) = x.^2;
    if (x > 0)
      A(1, 2) = x + 1;
    else
      A(1, 2) = x;
    end
    A(1, 3) = sin(x);

    Custom printing can be defined for certain types by passing a dictionary of
    "type" : "function" to the ``user_functions`` kwarg.  Alternatively, the
    dictionary value can be a list of tuples i.e., [(argument_test,
    cfunction_string)].  This can be used to call a custom Octave function.

    >>> f = Function('f')
    >>> g = Function('g')
    >>> custom_functions = {
    ...   "f": "existing_octave_fcn",
    ...   "g": [(lambda x: x.is_Matrix, "my_mat_fcn"),
    ...         (lambda x: not x.is_Matrix, "my_fcn")]
    ... }
    >>> mat = Matrix([[1, x]])
    >>> octave_code(f(x) + g(x) + g(mat), user_functions=custom_functions)
    'existing_octave_fcn(x) + my_fcn(x) + my_mat_fcn([1 x])'

    [FIXME: test loops with ``Indexed`` types and add here, see ``ccode``]
    """
    return OctaveCodePrinter(settings).doprint(expr, assign_to)


def print_octave_code(expr, **settings):
    """Prints the Octave/Matlab representation of the given expression.

       See `octave_code` for the meaning of the optional arguments.
    """
    print(octave_code(expr, **settings))
