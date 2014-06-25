"""
module for generating C, C++, Fortran77, Fortran90 and python routines that
evaluate sympy expressions. This module is work in progress. Only the
milestones with a '+' character in the list below have been completed.

--- How is sympy.utilities.codegen different from sympy.printing.ccode? ---

We considered the idea to extend the printing routines for sympy functions in
such a way that it prints complete compilable code, but this leads to a few
unsurmountable issues that can only be tackled with dedicated code generator:

- For C, one needs both a code and a header file, while the printing routines
  generate just one string. This code generator can be extended to support
  .pyf files for f2py.

- SymPy functions are not concerned with programming-technical issues, such
  as input, output and input-output arguments. Other examples are contiguous
  or non-contiguous arrays, including headers of other libraries such as gsl
  or others.

- It is highly interesting to evaluate several sympy functions in one C
  routine, eventually sharing common intermediate results with the help
  of the cse routine. This is more than just printing.

- From the programming perspective, expressions with constants should be
  evaluated in the code generator as much as possible. This is different
  for printing.

--- Basic assumptions ---

* A generic Routine data structure describes the routine that must be
  translated into C/Fortran/... code. This data structure covers all
  features present in one or more of the supported languages.

* Descendants from the CodeGen class transform multiple Routine instances
  into compilable code. Each derived class translates into a specific
  language.

* In many cases, one wants a simple workflow. The friendly functions in the
  last part are a simple api on top of the Routine/CodeGen stuff. They are
  easier to use, but are less powerful.

--- Milestones ---

+ First working version with scalar input arguments, generating C code,
  tests
+ Friendly functions that are easier to use than the rigorous
  Routine/CodeGen workflow.
+ Integer and Real numbers as input and output
+ Output arguments
+ InputOutput arguments
+ Sort input/output arguments properly
+ Contiguous array arguments (numpy matrices)
+ Also generate .pyf code for f2py (in autowrap module)
+ Isolate constants and evaluate them beforehand in double precision
+ Fortran 90

- Common Subexpression Elimination
- User defined comments in the generated code
- Optional extra include lines for libraries/objects that can eval special
  functions
- Test other C compilers and libraries: gcc, tcc, libtcc, gcc+gsl, ...
- Contiguous array arguments (sympy matrices)
- Non-contiguous array arguments (sympy matrices)
- ccode must raise an error when it encounters something that can not be
  translated into c. ccode(integrate(sin(x)/x, x)) does not make sense.
- Complex numbers as input and output
- A default complex datatype
- Include extra information in the header: date, user, hostname, sha1
  hash, ...
- Fortran 77
- C++
- Python
- ...

"""

from __future__ import print_function, division

import os

from sympy import __version__ as sympy_version
from sympy.core import Symbol, S, Expr, Tuple, Equality, Function
from sympy.core.compatibility import is_sequence, StringIO, string_types
from sympy.printing.codeprinter import AssignmentError
from sympy.printing.ccode import ccode, CCodePrinter
from sympy.printing.fcode import fcode, FCodePrinter
from sympy.tensor import Idx, Indexed, IndexedBase


__all__ = [
    # description of routines
    "Routine", "DataType", "default_datatypes", "get_default_datatype",
    "Argument", "InputArgument", "Result",
    # routines -> code
    "CodeGen", "CCodeGen", "FCodeGen",
    # friendly functions
    "codegen",
]


#
# Description of routines
#


class Routine(object):
    """Generic description of an evaluation routine for a set of sympy expressions.

       A CodeGen class can translate instances of this class into C/Fortran/...
       code. The routine specification covers all the features present in these
       languages. The CodeGen part must raise an exception when certain features
       are not present in the target language. For example, multiple return
       values are possible in Python, but not in C or Fortran. Another example:
       Fortran and Python support complex numbers, while C does not.
    """
    def __init__(self, name, expr, argument_sequence=None):
        """Initialize a Routine instance.

        ``name``
            A string with the name of this routine in the generated code
        ``expr``
            The sympy expression that the Routine instance will represent.  If
            given a list or tuple of expressions, the routine will be
            considered to have multiple return values.
        ``argument_sequence``
            Optional list/tuple containing arguments for the routine in a
            preferred order.  If omitted, arguments will be ordered
            alphabetically, but with all input aguments first, and then output
            or in-out arguments.

        A decision about whether to use output arguments or return values,
        is made depending on the mathematical expressions.  For an expression
        of type Equality, the left hand side is made into an OutputArgument
        (or an InOutArgument if appropriate).  Else, the calculated
        expression is the return values of the routine.

        A tuple of exressions can be used to create a routine with both
        return value(s) and output argument(s).

        """
        arg_list = []

        if is_sequence(expr):
            if not expr:
                raise ValueError("No expression given")
            expressions = Tuple(*expr)
        else:
            expressions = Tuple(expr)

        # local variables
        local_vars = set([i.label for i in expressions.atoms(Idx)])

        # symbols that should be arguments
        symbols = expressions.free_symbols - local_vars

        # Decide whether to use output argument or return value
        return_val = []
        output_args = []
        for expr in expressions:
            if isinstance(expr, Equality):
                out_arg = expr.lhs
                expr = expr.rhs
                if isinstance(out_arg, Indexed):
                    dims = tuple([ (S.Zero, dim - 1) for dim in out_arg.shape])
                    symbol = out_arg.base.label
                elif isinstance(out_arg, Symbol):
                    dims = []
                    symbol = out_arg
                else:
                    raise CodeGenError(
                        "Only Indexed or Symbol can define output arguments")

                if expr.has(symbol):
                    output_args.append(
                        InOutArgument(symbol, out_arg, expr, dimensions=dims))
                else:
                    output_args.append(OutputArgument(
                        symbol, out_arg, expr, dimensions=dims))

                # avoid duplicate arguments
                symbols.remove(symbol)
            else:
                return_val.append(Result(expr))

        # setup input argument list
        array_symbols = {}
        for array in expressions.atoms(Indexed):
            array_symbols[array.base.label] = array

        for symbol in sorted(symbols, key=str):
            if symbol in array_symbols:
                dims = []
                array = array_symbols[symbol]
                for dim in array.shape:
                    dims.append((S.Zero, dim - 1))
                metadata = {'dimensions': dims}
            else:
                metadata = {}

            arg_list.append(InputArgument(symbol, **metadata))

        output_args.sort(key=lambda x: str(x.name))
        arg_list.extend(output_args)

        if argument_sequence is not None:
            # if the user has supplied IndexedBase instances, we'll accept that
            new_sequence = []
            for arg in argument_sequence:
                if isinstance(arg, IndexedBase):
                    new_sequence.append(arg.label)
                else:
                    new_sequence.append(arg)
            argument_sequence = new_sequence

            missing = [x for x in arg_list if x.name not in argument_sequence]
            if missing:
                raise CodeGenArgumentListError("Argument list didn't specify: %s" %
                        ", ".join([str(m.name) for m in missing]), missing)

            # create redundant arguments to produce the requested sequence
            name_arg_dict = dict([(x.name, x) for x in arg_list])
            new_args = []
            for symbol in argument_sequence:
                try:
                    new_args.append(name_arg_dict[symbol])
                except KeyError:
                    new_args.append(InputArgument(symbol))
            arg_list = new_args

        self.name = name
        self.arguments = arg_list
        self.results = return_val
        self.local_vars = local_vars

    @property
    def variables(self):
        """Returns a set containing all variables possibly used in this routine.

        For routines with unnamed return values, the dummies that may or may
        not be used will be included in the set.
        """
        v = set(self.local_vars)
        for arg in self.arguments:
            v.add(arg.name)
        for res in self.results:
            v.add(res.result_var)
        return v

    @property
    def result_variables(self):
        """Returns a list of OutputArgument, InOutArgument and Result.

        If return values are present, they are at the end ot the list.
        """
        args = [arg for arg in self.arguments if isinstance(
            arg, (OutputArgument, InOutArgument))]
        args.extend(self.results)
        return args


class DataType(object):
    """Holds strings for a certain datatype in different programming languages."""
    def __init__(self, cname, fname, pyname):
        self.cname = cname
        self.fname = fname
        self.pyname = pyname


default_datatypes = {
    "int": DataType("int", "INTEGER*4", "int"),
    "float": DataType("double", "REAL*8", "float")
}


def get_default_datatype(expr):
    """Derives a decent data type based on the assumptions on the expression."""
    if expr.is_integer:
        return default_datatypes["int"]
    else:
        return default_datatypes["float"]


class Variable(object):
    """Represents a typed variable."""

    def __init__(self, name, datatype=None, dimensions=None, precision=None):
        """Initializes a Variable instance

           name  --  must be of class Symbol
           datatype  --  When not given, the data type will be guessed based
                         on the assumptions on the symbol argument.
           dimension  --  If present, the argument is interpreted as an array.
                          Dimensions must be a sequence containing tuples, i.e.
                          (lower, upper) bounds for each index of the array
           precision  --  FIXME
        """
        if not isinstance(name, Symbol):
            raise TypeError("The first argument must be a sympy symbol.")
        if datatype is None:
            datatype = get_default_datatype(name)
        elif not isinstance(datatype, DataType):
            raise TypeError("The (optional) `datatype' argument must be an instance of the DataType class.")
        if dimensions and not isinstance(dimensions, (tuple, list)):
            raise TypeError(
                "The dimension argument must be a sequence of tuples")

        self._name = name
        self._datatype = {
            'C': datatype.cname,
            'FORTRAN': datatype.fname,
            'PYTHON': datatype.pyname
        }
        self.dimensions = dimensions
        self.precision = precision

    @property
    def name(self):
        return self._name

    def get_datatype(self, language):
        """Returns the datatype string for the requested langage.

            >>> from sympy import Symbol
            >>> from sympy.utilities.codegen import Variable
            >>> x = Variable(Symbol('x'))
            >>> x.get_datatype('c')
            'double'
            >>> x.get_datatype('fortran')
            'REAL*8'
        """
        try:
            return self._datatype[language.upper()]
        except KeyError:
            raise CodeGenError("Has datatypes for languages: %s" %
                    ", ".join(self._datatype))


class Argument(Variable):
    """An abstract Argument data structure: a name and a data type.

       This structure is refined in the descendants below.
    """

    def __init__(self, name, datatype=None, dimensions=None, precision=None):
        """ See docstring of Variable.__init__
        """

        Variable.__init__(self, name, datatype, dimensions, precision)


class InputArgument(Argument):
    pass


class ResultBase(object):
    """Base class for all ``outgoing'' information from a routine

       Objects of this class stores a sympy expression, and a sympy object
       representing a result variable that will be used in the generated code
       only if necessary.
   """
    def __init__(self, expr, result_var):
        self.expr = expr
        self.result_var = result_var


class OutputArgument(Argument, ResultBase):
    """OutputArgument are always initialized in the routine
    """
    def __init__(self, name, result_var, expr, datatype=None, dimensions=None, precision=None):
        """ See docstring of Variable.__init__
        """
        Argument.__init__(self, name, datatype, dimensions, precision)
        ResultBase.__init__(self, expr, result_var)


class InOutArgument(Argument, ResultBase):
    """InOutArgument are never initialized in the routine
    """

    def __init__(self, name, result_var, expr, datatype=None, dimensions=None, precision=None):
        """ See docstring of Variable.__init__
        """
        Argument.__init__(self, name, datatype, dimensions, precision)
        ResultBase.__init__(self, expr, result_var)


class Result(ResultBase):
    """An expression for a scalar return value.

       The name result is used to avoid conflicts with the reserved word
       'return' in the python language. It is also shorter than ReturnValue.

    """

    def __init__(self, expr, datatype=None, precision=None):
        """Initialize a (scalar) return value.

           The second argument is optional. When not given, the data type will
           be guessed based on the assumptions on the expression argument.
        """
        if not isinstance(expr, Expr):
            raise TypeError("The first argument must be a sympy expression.")

        temp_var = Variable(Symbol('result_%s' % hash(expr)),
                datatype=datatype, dimensions=None, precision=precision)
        ResultBase.__init__(self, expr, temp_var.name)
        self._temp_variable = temp_var

    def get_datatype(self, language):
        return self._temp_variable.get_datatype(language)


#
# Transformation of routine objects into code
#

class CodeGen(object):
    """Abstract class for the code generators."""

    def __init__(self, project="project"):
        """Initialize a code generator.

           Derived classes will offer more options that affect the generated
           code.
        """
        self.project = project

    def write(self, routines, prefix, to_files=False, header=True, empty=True):
        """Writes all the source code files for the given routines.

            The generate source is returned as a list of (filename, contents)
            tuples, or is written to files (see options). Each filename consists
            of the given prefix, appended with an appropriate extension.

            ``routines``
                A list of Routine instances to be written
            ``prefix``
                The prefix for the output files
            ``to_files``
                When True, the output is effectively written to files.
                [DEFAULT=False] Otherwise, a list of (filename, contents)
                tuples is returned.
            ``header``
                When True, a header comment is included on top of each source
                file. [DEFAULT=True]
            ``empty``
                When True, empty lines are included to structure the source
                files. [DEFAULT=True]

        """
        if to_files:
            for dump_fn in self.dump_fns:
                filename = "%s.%s" % (prefix, dump_fn.extension)
                with open(filename, "w") as f:
                    dump_fn(self, routines, f, prefix, header, empty)
        else:
            result = []
            for dump_fn in self.dump_fns:
                filename = "%s.%s" % (prefix, dump_fn.extension)
                contents = StringIO()
                dump_fn(self, routines, contents, prefix, header, empty)
                result.append((filename, contents.getvalue()))
            return result

    def dump_code(self, routines, f, prefix, header=True, empty=True):
        """Write the code file by calling language specific methods in correct order

        The generated file contains all the definitions of the routines in
        low-level code and refers to the header file if appropriate.

        :Arguments:

        routines
            A list of Routine instances
        f
            A file-like object to write the file to
        prefix
            The filename prefix, used to refer to the proper header file. Only
            the basename of the prefix is used.

        :Optional arguments:

        header
            When True, a header comment is included on top of each source file.
            [DEFAULT=True]
        empty
            When True, empty lines are included to structure the source files.
            [DEFAULT=True]
        """

        code_lines = self._preprocessor_statements(prefix)

        for routine in routines:
            if empty:
                code_lines.append("\n")
            code_lines.extend(self._get_routine_opening(routine))
            code_lines.extend(self._declare_arguments(routine))
            code_lines.extend(self._declare_locals(routine))
            if empty:
                code_lines.append("\n")
            code_lines.extend(self._call_printer(routine))
            if empty:
                code_lines.append("\n")
            code_lines.extend(self._get_routine_ending(routine))

        code_lines = self._indent_code(''.join(code_lines))

        if header:
            code_lines = ''.join(self._get_header() + [code_lines])

        if code_lines:
            f.write(code_lines)


class CodeGenError(Exception):
    pass


class CodeGenArgumentListError(Exception):
    @property
    def missing_args(self):
        return self.args[1]


header_comment = """Code generated with sympy %(version)s

See http://www.sympy.org/ for more information.

This file is part of '%(project)s'
"""


class CCodeGen(CodeGen):
    """
    Generator for C code

    The .write() method inherited from CodeGen will output a code file and an
    inteface file, <prefix>.c and <prefix>.h respectively.
    """

    code_extension = "c"
    interface_extension = "h"

    def _get_header(self):
        """Writes a common header for the generated files."""
        code_lines = []
        code_lines.append("/" + "*"*78 + '\n')
        tmp = header_comment % {"version": sympy_version,
            "project": self.project}
        for line in tmp.splitlines():
            code_lines.append(" *%s*\n" % line.center(76))
        code_lines.append(" " + "*"*78 + "/\n")
        return code_lines

    def get_prototype(self, routine):
        """Returns a string for the function prototype for the given routine.

           If the routine has multiple result objects, an CodeGenError is
           raised.

           See: http://en.wikipedia.org/wiki/Function_prototype
        """
        if len(routine.results) > 1:
            raise CodeGenError("C only supports a single or no return value.")
        elif len(routine.results) == 1:
            ctype = routine.results[0].get_datatype('C')
        else:
            ctype = "void"

        type_args = []
        for arg in routine.arguments:
            name = ccode(arg.name)
            if arg.dimensions:
                type_args.append((arg.get_datatype('C'), "*%s" % name))
            elif isinstance(arg, ResultBase):
                type_args.append((arg.get_datatype('C'), "&%s" % name))
            else:
                type_args.append((arg.get_datatype('C'), name))
        arguments = ", ".join([ "%s %s" % t for t in type_args])
        return "%s %s(%s)" % (ctype, routine.name, arguments)

    def _preprocessor_statements(self, prefix):
        code_lines = []
        code_lines.append("#include \"%s.h\"\n" % os.path.basename(prefix))
        code_lines.append("#include <math.h>\n")
        return code_lines

    def _get_routine_opening(self, routine):
        prototype = self.get_prototype(routine)
        return ["%s {\n" % prototype]

    def _declare_arguments(self, routine):
        # arguments are declared in prototype
        return []

    def _declare_locals(self, routine):
        # loop variables are declared in loop statement
        return []

    def _call_printer(self, routine):
        code_lines = []
        for result in routine.result_variables:
            if isinstance(result, Result):
                assign_to = None
            elif isinstance(result, (OutputArgument, InOutArgument)):
                assign_to = result.result_var

            try:
                constants, not_c, c_expr = ccode(
                    result.expr, assign_to=assign_to, human=False)
            except AssignmentError:
                assign_to = result.result_var
                code_lines.append(
                    "%s %s;\n" % (result.get_datatype('c'), str(assign_to)))
                constants, not_c, c_expr = ccode(
                    result.expr, assign_to=assign_to, human=False)

            for name, value in sorted(constants, key=str):
                code_lines.append("double const %s = %s;\n" % (name, value))
            if assign_to:
                code_lines.append("%s\n" % c_expr)
            else:
                code_lines.append("   return %s;\n" % c_expr)
        return code_lines

    def _indent_code(self, codelines):
        p = CCodePrinter()
        return p.indent_code(codelines)

    def _get_routine_ending(self, routine):
        return ["}\n"]

    def dump_c(self, routines, f, prefix, header=True, empty=True):
        self.dump_code(routines, f, prefix, header, empty)
    dump_c.extension = code_extension
    dump_c.__doc__ = CodeGen.dump_code.__doc__

    def dump_h(self, routines, f, prefix, header=True, empty=True):
        """Writes the C header file.

           This file contains all the function declarations.

           :Arguments:

           routines
                A list of Routine instances
           f
                A file-like object to write the file to
           prefix
                The filename prefix, used to construct the include guards.

           :Optional arguments:

           header
                When True, a header comment is included on top of each source
                file. [DEFAULT=True]
           empty
                When True, empty lines are included to structure the source
                files. [DEFAULT=True]
        """
        if header:
            print(''.join(self._get_header()), file=f)
        guard_name = "%s__%s__H" % (self.project.replace(
            " ", "_").upper(), prefix.replace("/", "_").upper())
        # include guards
        if empty:
            print(file=f)
        print("#ifndef %s" % guard_name, file=f)
        print("#define %s" % guard_name, file=f)
        if empty:
            print(file=f)
        # declaration of the function prototypes
        for routine in routines:
            prototype = self.get_prototype(routine)
            print("%s;" % prototype, file=f)
        # end if include guards
        if empty:
            print(file=f)
        print("#endif", file=f)
        if empty:
            print(file=f)
    dump_h.extension = interface_extension

    # This list of dump functions is used by CodeGen.write to know which dump
    # functions it has to call.
    dump_fns = [dump_c, dump_h]


class FCodeGen(CodeGen):
    """
    Generator for Fortran 95 code

    The .write() method inherited from CodeGen will output a code file and an
    inteface file, <prefix>.f90 and <prefix>.h respectively.
    """

    code_extension = "f90"
    interface_extension = "h"

    def __init__(self, project='project'):
        CodeGen.__init__(self, project)

    def _get_symbol(self, s):
        """returns the symbol as fcode print it"""
        return fcode(s).strip()

    def _get_header(self):
        """Writes a common header for the generated files."""
        code_lines = []
        code_lines.append("!" + "*"*78 + '\n')
        tmp = header_comment % {"version": sympy_version,
            "project": self.project}
        for line in tmp.splitlines():
            code_lines.append("!*%s*\n" % line.center(76))
        code_lines.append("!" + "*"*78 + '\n')
        return code_lines

    def _preprocessor_statements(self, prefix):
        return []

    def _get_routine_opening(self, routine):
        """
        Returns the opening statements of the fortran routine
        """
        code_list = []
        if len(routine.results) > 1:
            raise CodeGenError(
                "Fortran only supports a single or no return value.")
        elif len(routine.results) == 1:
            result = routine.results[0]
            code_list.append(result.get_datatype('fortran'))
            code_list.append("function")
        else:
            code_list.append("subroutine")

        args = ", ".join("%s" % self._get_symbol(arg.name)
                for arg in routine.arguments)

        # name of the routine + arguments
        code_list.append("%s(%s)\n" % (routine.name, args))
        code_list = [ " ".join(code_list) ]

        code_list.append('implicit none\n')
        return code_list

    def _declare_arguments(self, routine):
        # argument type declarations
        code_list = []
        array_list = []
        scalar_list = []
        for arg in routine.arguments:

            if isinstance(arg, InputArgument):
                typeinfo = "%s, intent(in)" % arg.get_datatype('fortran')
            elif isinstance(arg, InOutArgument):
                typeinfo = "%s, intent(inout)" % arg.get_datatype('fortran')
            elif isinstance(arg, OutputArgument):
                typeinfo = "%s, intent(out)" % arg.get_datatype('fortran')
            else:
                raise CodeGenError("Unkown Argument type: %s" % type(arg))

            fprint = self._get_symbol

            if arg.dimensions:
                # fortran arrays start at 1
                dimstr = ", ".join(["%s:%s" % (
                    fprint(dim[0] + 1), fprint(dim[1] + 1))
                    for dim in arg.dimensions])
                typeinfo += ", dimension(%s)" % dimstr
                array_list.append("%s :: %s\n" % (typeinfo, fprint(arg.name)))
            else:
                scalar_list.append("%s :: %s\n" % (typeinfo, fprint(arg.name)))

        # scalars first, because they can be used in array declarations
        code_list.extend(scalar_list)
        code_list.extend(array_list)

        return code_list

    def _declare_locals(self, routine):
        code_list = []
        for var in sorted(routine.local_vars, key=str):
            typeinfo = get_default_datatype(var)
            code_list.append("%s :: %s\n" % (
                typeinfo.fname, self._get_symbol(var)))
        return code_list

    def _get_routine_ending(self, routine):
        """
        Returns the closing statements of the fortran routine
        """
        if len(routine.results) == 1:
            return ["end function\n"]
        else:
            return ["end subroutine\n"]

    def get_interface(self, routine):
        """Returns a string for the function interface for the given routine and
           a single result object, which can be None.

           If the routine has multiple result objects, a CodeGenError is
           raised.

           See: http://en.wikipedia.org/wiki/Function_prototype

        """
        prototype = [ "interface\n" ]
        prototype.extend(self._get_routine_opening(routine))
        prototype.extend(self._declare_arguments(routine))
        prototype.extend(self._get_routine_ending(routine))
        prototype.append("end interface\n")

        return "".join(prototype)

    def _call_printer(self, routine):
        declarations = []
        code_lines = []
        for result in routine.result_variables:
            if isinstance(result, Result):
                assign_to = routine.name
            elif isinstance(result, (OutputArgument, InOutArgument)):
                assign_to = result.result_var

            constants, not_fortran, f_expr = fcode(result.expr,
                assign_to=assign_to, source_format='free', human=False)

            for obj, v in sorted(constants, key=str):
                t = get_default_datatype(obj)
                declarations.append(
                    "%s, parameter :: %s = %s\n" % (t.fname, obj, v))
            for obj in sorted(not_fortran, key=str):
                t = get_default_datatype(obj)
                if isinstance(obj, Function):
                    name = obj.func
                else:
                    name = obj
                declarations.append("%s :: %s\n" % (t.fname, name))

            code_lines.append("%s\n" % f_expr)
        return declarations + code_lines

    def _indent_code(self, codelines):
        p = FCodePrinter({'source_format': 'free', 'human': False})
        return p.indent_code(codelines)

    def dump_f95(self, routines, f, prefix, header=True, empty=True):
        # check that symbols are unique with ignorecase
        for r in routines:
            lowercase = set([str(x).lower() for x in r.variables])
            orig_case = set([str(x) for x in r.variables])
            if len(lowercase) < len(orig_case):
                raise CodeGenError("Fortran ignores case. Got symbols: %s" %
                        (", ".join([str(var) for var in r.variables])))
        self.dump_code(routines, f, prefix, header, empty)
    dump_f95.extension = code_extension
    dump_f95.__doc__ = CodeGen.dump_code.__doc__

    def dump_h(self, routines, f, prefix, header=True, empty=True):
        """Writes the interface to a header file.

           This file contains all the function declarations.

           :Arguments:

           routines
                A list of Routine instances
           f
                A file-like object to write the file to
           prefix
                The filename prefix

           :Optional arguments:

           header
                When True, a header comment is included on top of each source
                file. [DEFAULT=True]
           empty
                When True, empty lines are included to structure the source
                files. [DEFAULT=True]
        """
        if header:
            print(''.join(self._get_header()), file=f)
        if empty:
            print(file=f)
        # declaration of the function prototypes
        for routine in routines:
            prototype = self.get_interface(routine)
            f.write(prototype)
        if empty:
            print(file=f)
    dump_h.extension = interface_extension

    # This list of dump functions is used by CodeGen.write to know which dump
    # functions it has to call.
    dump_fns = [dump_f95, dump_h]


def get_code_generator(language, project):
    CodeGenClass = {"C": CCodeGen, "F95": FCodeGen}.get(language.upper())
    if CodeGenClass is None:
        raise ValueError("Language '%s' is not supported." % language)
    return CodeGenClass(project)


#
# Friendly functions
#


def codegen(
    name_expr, language, prefix, project="project", to_files=False, header=True, empty=True,
        argument_sequence=None):
    """Write source code for the given expressions in the given language.

    :Mandatory Arguments:

    ``name_expr``
        A single (name, expression) tuple or a list of (name, expression)
        tuples. Each tuple corresponds to a routine.  If the expression is an
        equality (an instance of class Equality) the left hand side is
        considered an output argument.
    ``language``
            A string that indicates the source code language. This is case
            insensitive. For the moment, only 'C' and 'F95' is supported.
    ``prefix``
            A prefix for the names of the files that contain the source code.
            Proper (language dependent) suffixes will be appended.

    :Optional Arguments:

    ``project``
        A project name, used for making unique preprocessor instructions.
        [DEFAULT="project"]
    ``to_files``
        When True, the code will be written to one or more files with the given
        prefix, otherwise strings with the names and contents of these files
        are returned. [DEFAULT=False]
    ``header``
        When True, a header is written on top of each source file.
        [DEFAULT=True]
    ``empty``
        When True, empty lines are used to structure the code.  [DEFAULT=True]
    ``argument_sequence``
        sequence of arguments for the routine in a preferred order.  A
        CodeGenError is raised if required arguments are missing.  Redundant
        arguments are used without warning.

        If omitted, arguments will be ordered alphabetically, but with all
        input aguments first, and then output or in-out arguments.

    >>> from sympy.utilities.codegen import codegen
    >>> from sympy.abc import x, y, z
    >>> [(c_name, c_code), (h_name, c_header)] = codegen(
    ...     ("f", x+y*z), "C", "test", header=False, empty=False)
    >>> print(c_name)
    test.c
    >>> print(c_code)
    #include "test.h"
    #include <math.h>
    double f(double x, double y, double z) {
      return x + y*z;
    }
    >>> print(h_name)
    test.h
    >>> print(c_header)
    #ifndef PROJECT__TEST__H
    #define PROJECT__TEST__H
    double f(double x, double y, double z);
    #endif

    """

    # Initialize the code generator.
    code_gen = get_code_generator(language, project)

    # Construct the routines based on the name_expression pairs.
    #  mainly the input arguments require some work
    routines = []
    if isinstance(name_expr[0], string_types):
        # single tuple is given, turn it into a singleton list with a tuple.
        name_expr = [name_expr]

    for name, expr in name_expr:
        routines.append(Routine(name, expr, argument_sequence))

    # Write the code.
    return code_gen.write(routines, prefix, to_files, header, empty)
