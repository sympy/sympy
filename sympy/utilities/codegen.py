"""
module for generating C, C++, Fortran77, Fortran90 and python routines that
evaluate sympy expressions. This module is work in progress. Only the
milestones with a '+' character in the list below have been completed.


--- How is sympy.utilities.codegen different from sympy.printing.ccode? ---

We considered the idea to extend the printing routines for sympy functions in
such a way that it prints complete compilable code, but this leads to a few
unsurmountable issues that can only be tackled with dedicated code generator:

- For C, one needs both a code and a header file, while the printing routines
generate just one string. This code generator can be extended to support .pyf
files for f2py.

- Sympy functions are not concerned with programming-technical issues, such as
input, output and input-output arguments. Other examples are contiguous or
non-contiguous arrays, including headers of other libraries such as gsl or others.

- It is highly interesting to evaluate several sympy functions in one C routine,
eventually sharing common intermediate results with the help of the cse routine.
This is more than just printing.

- From the programming perspective, expressions with constants should be
evaluated in the code generator as much as possible. This is different for
printing.


--- Basic assumptions ---

* A generic Routine data structure describes the routine that must be translated
  into C/Fortran/... code. This data structure covers all features present in
  one or more of the supported languages.

* Descendants from the CodeGen class transform multiple Routine instances into
  compilable code. Each derived class translates into a specific language.

* In many cases, one wants a simple workflow. The friendly functions in the last
  part are a simple api on top of the Routine/CodeGen stuff. They are easier to
  use, but are less powerful.


--- Milestones ---

+ First working version with scalar input arguments, generating C code, tests
+ Friendly functions that are easier to use than the rigorous Routine/CodeGen
  workflow.
+ Integer and Real numbers as input and output
- Optional extra include lines for libraries/objects that can eval special
  functions
- Test other C compilers and libraries: gcc, tcc, libtcc, gcc+gsl, ...
- Output arguments
- InputOutput arguments
- Sort input/output arguments properly
- Contiguous array arguments (sympy matrices)
- Non-contiguous array arguments (sympy matrices)
- ccode must raise an error when it encounters something that can not be
  translated into c. ccode(integrate(sin(x)/x, x)) does not make sense.
- Complex numbers as input and output
- Also generate .pyf code for f2py
- A default complex datatype
- Include extra information in the header: date, user, hostname, sha1 hash, ...
- Isolate constants and evaluate them beforehand in double precision
- Common Subexpression Elimination
- User defined comments in the generated code
- Fortran 77
- Fortran 90
- C++
- Python
- ...
"""


from sympy.core.symbol import Symbol
from sympy.core.expr import Expr
from sympy.core.containers import Tuple
from sympy.printing.ccode import ccode
from sympy.printing.fcode import fcode
from sympy.tensor import Idx, IndexedElement
from sympy.core.relational import Equality
from sympy.utilities import flatten

from StringIO import StringIO
import sympy
import os


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
    def __init__(self, name, expr):
        """Initialize a Routine instance.

           Arguments:
             name  --  A string with the name of this routine in the generated
                       code
             expr  --  The sympy expression that the Routine instance will represent.
                       If given a list or tuple of expressions, the routine
                       will be considered to have multiple return values.

            A decision about whether to use output arguments or return values,
            is made depending on the mathematical expressions.  For an expression
            of type Equality, the left hand side is made into an OutputArgument
            (or an InOutArgument if appropriate).  Else, the calculated
            expression is the return values of the routine.

            A tuple of exressions can be used to create a routine with both
            return value(s) and output argument(s).

        """
        arg_list = []

        if isinstance(expr, (list, tuple)):
            if not expr:
                raise ValueError("No expression given")
            expressions = Tuple(*expr)
        else:
            expressions = Tuple(expr)

        # local variables
        local_vars = set([i.label for i in expressions.atoms(Idx)])

        # symbols that should be arguments
        symbols = expressions.atoms(Symbol) - local_vars

        # Decide whether to use output argument or return value
        return_val = []
        output_args = []
        for expr in expressions:
            if isinstance(expr, Equality):
                out_arg = expr.lhs
                expr = expr.rhs
                if isinstance(out_arg, IndexedElement):
                    dims = out_arg.ranges
                    symbol = out_arg.stem.label
                elif isinstance(out_arg, Symbol):
                    dims = []
                    symbol = out_arg
                else:
                    raise CodeGenError("Only IndexedElement or Symbol can define output arguments")

                if expr.has(symbol):
                    output_args.append(InOutArgument(symbol, out_arg, expr, dimensions=dims))
                else:
                    output_args.append(OutputArgument(symbol, out_arg, expr, dimensions=dims))

                # avoid duplicate arguments
                symbols.remove(symbol)
            else:
                return_val.append(Result(expr))

        # setup input argument list
        array_symbols = {}
        for array in expressions.atoms(IndexedElement):
            array_symbols[array.stem.label] = array

        for symbol in sorted(symbols):
            if symbol in array_symbols:
                dims = []
                array = array_symbols[symbol]
                for i in array.indices:
                    if i.lower == None:
                        dims.append((S.Zero, S.Zero))
                    else:
                        dims.append((i.lower, i.upper))
                metadata = {'dimensions': dims}
            else:
                metadata = {}

            arg_list.append(InputArgument(symbol, **metadata))

        arg_list.extend(output_args)

        self.name = name
        self.arguments = arg_list
        self.results = return_val
        self.local_vars = local_vars

    @property
    def result_variables(self):
        """Returns a tuple of OutputArgument, InOutArgument and Result."""
        args = [arg for arg in self.arguments if isinstance(arg, (OutputArgument, InOutArgument))]
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


class Argument(object):
    """An abstract Argument data structure: a name and a data type.

       This structure is refined in the descendants below.
    """

    def __init__(self, name, datatype=None, dimensions=None, precision=None):
        """Initialize an input argument.

           name  --  must be of class Symbol
           datatype  --  When not given, the data type will be guessed based
                         on the assumptions on the symbol argument.
           dimension  --  If present, the argument is interpreted as an array.
                          Dimensions must be a sequence containing tuples
                          of first and last index in the range.
           precision  --  FIXME
        """

        if not isinstance(name, Symbol):
            raise TypeError("The first argument must be a sympy symbol.")
        if datatype is None:
            datatype = get_default_datatype(name)
        elif not isinstance(datatype, DataType):
            raise TypeError("The (optional) `datatype' argument must be an instance of the DataType class.")
        if dimensions and not isinstance(dimensions, (tuple, list)):
            raise TypeError("The dimension argument must be a sequence of tuples")
        self.name = name
        self.datatype = datatype
        self.dimensions = dimensions
        self.precision = precision

    def get_symbols(self):
        """Returns a set of all symbols related to this argument.

        Scalar arguments return themselves in a set, while array arguments return
        the array variable as well as all symbols the specifiy dimensions.
        """
        if self.dimensions:
            symbs = set(flatten(self.dimensions))
            symbs.add(self.name)
            return symbs
        else:
            return set([self.name])

class InputArgument(Argument):
    pass

class ResultBase(object):

    @property
    def needs_initialization(self):
        return self._need_initialization

    def _prepare_expr(self):
        """Depending on the expression, it may need to be prepared for loops

        Example:

        The math expression denoting a matrix vector product is

        y(i) = A(i,j)*x(j)          (imlicit summation over j)

        Obviously, the code that calculates the matrix vector product must change
        the right hand side a little in order to store intermediate results
        correctly ....................
                                      |
        y = 0                         |
        do i=1,m                      |
            do j = 1,n                V
                y(i) = A(i,j)*x(i) + y(i)
            end do j
        end do i

        This would not be necessary if the left hand side has all indices that are on
        the right hand side.  E.g. for a dyadic product:

        do i=1,m
            do j = 1,n
                A(i,j) = x(i)*y(j)
            end do j
        end do i

        """
        rhs_loops = set([ i.label for i in self.expr.atoms(Idx) ])
        lhs_loops = set([ i.label for i in self.result_var.atoms(Idx) ])
        if rhs_loops - lhs_loops:
            self.expr = self.expr + self.result_var

class OutputArgument(Argument, ResultBase):
    """OutputArgument are always initialized in the routine
    """
    _need_initialization = True
    def __init__(self, name, result_var, expr, datatype=None, dimensions=None, precision=None):
        Argument.__init__(self, name, datatype, dimensions, precision)
        self.expr = expr
        self.result_var = result_var
        self._prepare_expr()

class InOutArgument(Argument, ResultBase):
    """InOutArgument are never initialized in the routine
    """
    _need_initialization = False

    def __init__(self, name, result_var, expr, datatype=None, dimensions=None, precision=None):
        Argument.__init__(self, name, datatype, dimensions, precision)
        self.expr = expr
        self.result_var = result_var


class Result(ResultBase):
    """An expression for a scalar return value.

       The name result is used to avoid conflicts with the reserved word
       'return' in the python language. It is also shorter than ReturnValue.
    """
    _need_initialization = False

    def __init__(self, expr, datatype=None):
        """Initialize a (scalar) return value.

           The second argument is optional. When not given, the data type will
           be guessed based on the assumptions on the expression argument.
        """
        if not isinstance(expr, Expr):
            raise TypeError("The first argument must be a sympy expression.")
        if datatype is None:
            datatype = get_default_datatype(expr)
        elif not isinstance(datatype, DataType):
            raise TypeError("The (optional) second argument must be an instance of the DataType class.")
        self.expr = expr
        self.datatype = datatype
        self.result_var = Symbol('result_%s'%hash(expr))
        self._prepare_expr()

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

           Arguments:
             routines  --  A list of Routine instances to be written
             prefix  --  The prefix for the output files

           Optional arguments:
             to_files  --  When True, the output is effectively written to
                           files. [DEFAULT=False] Otherwise, a list of
                           (filename, contents) tuples is returned.
             header  --  When True, a header comment is included on top of each
                         source file. [DEFAULT=True]
             empty  --  When True, empty lines are included to structure the
                        source files. [DEFAULT=True]
        """
        if to_files:
            for dump_fn in self.dump_fns:
                filename = "%s.%s" % (prefix, dump_fn.extension)
                f = file(filename, "w")
                dump_fn(self, routines, f, prefix, header, empty)
                f.close()
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

           Arguments:
             routines  --  a list of Routine instances
             f  --  a file-like object to write the file to
             prefix  --  the filename prefix, used to refer to the proper header
                         file. Only the basename of the prefix is used.

           Optional arguments:
             header  --  When True, a header comment is included on top of each
                         source file. [DEFAULT=True]
             empty  --  When True, empty lines are included to structure the
                        source files. [DEFAULT=True]
        """

        code_lines = []
        if header:
            code_lines.extend(self._get_header())

        code_lines.extend(self._preprosessor_statements(prefix))
        for routine in routines:
            if empty: code_lines.append("\n")
            code_lines.extend(self._get_routine_opening(routine))
            code_lines.extend(self._declare_arguments(routine))
            code_lines.extend(self._declare_locals(routine))
            if empty: code_lines.append("\n")
            code_lines.extend(self._call_printer(routine))
            if empty: code_lines.append("\n")
            code_lines.extend(self._get_routine_ending(routine))

        print >> f, ''.join(code_lines),

class CodeGenError(Exception):
    pass


header_comment = """Code generated with sympy %(version)s

See http://www.sympy.org/ for more information.

This file is part of '%(project)s'
"""


class CCodeGen(CodeGen):

    def _get_header(self):
        """Writes a common header for the generated files."""
        code_lines = []
        code_lines.append("/" + "*"*78 + '\n')
        tmp = header_comment % {"version": sympy.__version__, "project": self.project}
        for line in tmp.splitlines():
            code_lines.append(" *%s*\n" % line.center(76))
        code_lines.append(" " + "*"*78 + "/\n")
        return code_lines

    def _get_prototype(self, routine):
        """Returns a string for the function prototype for the given routine.

           If the routine has multiple result objects, an CodeGenError is
           raised.

           See: http://en.wikipedia.org/wiki/Function_prototype
        """
        if len(routine.results) > 1:
            raise CodeGenError("C only supports a single or no return value.")
        elif len(routine.results) == 1:
            ctype = routine.results[0].datatype.cname
        else:
            ctype = "void"
        arguments = ", ".join([ "%s %s" % (arg.datatype.cname, arg.name)
                for arg in routine.arguments ])
        return "%s %s(%s)" % (ctype, routine.name, arguments)

    def _preprosessor_statements(self, prefix):
        code_lines = []
        code_lines.append("#include \"%s.h\"\n" % os.path.basename(prefix))
        code_lines.append("#include <math.h>\n")
        return code_lines

    def _get_routine_opening(self, routine):
        prototype = self._get_prototype(routine)
        return ["%s {\n" % prototype]

    def _declare_arguments(self, routine):
        # arguments are declared in prototype
        return []

    def _declare_locals(self, routine):
        code_list = []
        for var in sorted(routine.local_vars, key=str):
            typeinfo = get_default_datatype(var)
            code_list.append("%s %s\n" % (typeinfo.cname, var))
        return code_list

    def _call_printer(self, routine):
        for result in routine.result_variables:
            if isinstance(result, Result):
                return ["   return %s;\n" % ccode(result.expr)]
            elif isinstance(result, (OutputArgument, InOutArgument)):
                raise NotImplementedError

    def _get_routine_ending(self, routine):
        return ["}\n"]

    def dump_c(self, routines, f, prefix, header=True, empty=True):
        self.dump_code(routines, f, prefix, header, empty)
    dump_c.extension = "c"

    def dump_h(self, routines, f, prefix, header=True, empty=True):
        """Writes the C header file.

           This file contains all the function declarations.

           Arguments:
             routines  --  a list of Routine instances
             f  --  a file-like object to write the file to
             prefix  --  the filename prefix, used to construct the include
                         guards.

           Optional arguments:
             header  --  When True, a header comment is included on top of each
                         source file. [DEFAULT=True]
             empty  --  When True, empty lines are included to structure the
                        source files. [DEFAULT=True]
        """
        if header:
            print >> f, ''.join(self._get_header())
        guard_name = "%s__%s__H" % (self.project.replace(" ", "_").upper(), prefix.replace("/", "_").upper())
        # include guards
        if empty: print >> f
        print >> f, "#ifndef %s" % guard_name
        print >> f, "#define %s" % guard_name
        if empty: print >> f
        # declaration of the function prototypes
        for routine in routines:
            prototype = self._get_prototype(routine)
            print >> f, "%s;" % prototype
        # end if include guards
        if empty: print >> f
        print >> f, "#endif"
        if empty: print >> f
    dump_h.extension = "h"

    # This list of dump functions is used by CodeGen.write to know which dump
    # functions it has to call.
    dump_fns = [dump_c, dump_h]

class FCodeGen(CodeGen):
    """
    Generator for Fortran 95 code
    """

    def __init__(self, project='project'):
        CodeGen.__init__(self, project)

    def _get_header(self):
        """Writes a common header for the generated files."""
        code_lines = []
        code_lines.append("!" + "*"*78 + '\n')
        tmp = header_comment % {"version": sympy.__version__, "project": self.project}
        for line in tmp.splitlines():
            code_lines.append("!*%s*\n" % line.center(76))
        code_lines.append("!" + "*"*78 + '\n')
        return code_lines

    def _preprosessor_statements(self, prefix):
        return []

    def _get_routine_opening(self, routine):
        """
        Returns the opening statements of the fortran routine
        """
        code_list = []
        if len(routine.results) > 1:
            raise CodeGenError("Fortran only supports a single or no return value.")
        elif len(routine.results) == 1:
            result = routine.results[0]
            code_list.append(result.datatype.fname)
            code_list.append("function")
        else:
            code_list.append("subroutine")

        # name of the routine + arguments
        code_list.append("%s(%s)\n" % (routine.name,
            ", ".join("%s" % arg.name for arg in routine.arguments)))
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
                typeinfo = "%s, intent(in)" % arg.datatype.fname
            elif isinstance(arg, InOutArgument):
                typeinfo = "%s, intent(inout)" % arg.datatype.fname
            elif isinstance(arg, OutputArgument):
                typeinfo = "%s, intent(out)" % arg.datatype.fname
            else:
                raise CodeGenError("Unkown Argument type: %s"%type(arg))

            if arg.dimensions:
                # fortran arrays start at 1
                dimstr = ", ".join(["%s:%s"%(dim[0]+1, dim[1]+1)
                    for dim in arg.dimensions])
                typeinfo += ", dimension(%s)" % dimstr
                array_list.append("%s :: %s\n" % (typeinfo, arg.name))
            else:
                scalar_list.append("%s :: %s\n" % (typeinfo, arg.name))

        # scalars first, because they can be used in array declarations
        code_list.extend(scalar_list)
        code_list.extend(array_list)

        return code_list

    def _declare_locals(self, routine):
        code_list = []
        for var in sorted(routine.local_vars, key=str):
            typeinfo = get_default_datatype(var)
            code_list.append("%s :: %s\n" % (typeinfo.fname, var))

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

    def _init_resultvars(self, routine):
        """Returns codelines that intialize the result variables if applicable."""
        code_lines = []
        for arg in routine.arguments:
            if isinstance(arg, OutputArgument):
                if arg.datatype.fname == 'REAL*8':
                    code_lines.append("%s = 0.d0\n" % arg.name)
                elif arg.datatype.fname == 'INTEGER*4':
                    code_lines.append("%s = 0\n" % arg.name)
                else:
                    raise NotImplementedError
        if routine.results:
            code_lines.append("%s = 0.d0\n" % routine.name)

        return code_lines

    def _call_printer(self, routine):
        code_lines = []
        for result in routine.result_variables:
            if isinstance(result, Result):
                assign_to = routine.name
            elif isinstance(result, (OutputArgument, InOutArgument)):
                assign_to = result.result_var

            constants, not_fortran, f_expr = fcode(result.expr,
                assign_to=assign_to, source_format='free', human=False)
            code_lines.extend(constants)
            if result.needs_initialization:
                code_lines.extend(self._init_resultvars(routine))
            code_lines.append("%s\n" % f_expr)
        return code_lines

    def dump_f95(self, routines, f, prefix, header=True, empty=True):
        self.dump_code(routines, f, prefix, header, empty)
    dump_f95.extension = "f90"

    def dump_h(self, routines, f, prefix, header=True, empty=True):
        """Writes the interface  header file.

           This file contains all the function declarations.

           Arguments:
             routines  --  a list of Routine instances
             f  --  a file-like object to write the file to
             prefix  --  the filename prefix, used to construct the include
                         guards.

           Optional arguments:
             header  --  When True, a header comment is included on top of each
                         source file. [DEFAULT=True]
             empty  --  When True, empty lines are included to structure the
                        source files. [DEFAULT=True]
        """
        if header:
            print >> f, ''.join(self._get_header())
        if empty: print >> f
        # declaration of the function prototypes
        for routine in routines:
            prototype  = self.get_interface(routine)
            print >> f, prototype,
        if empty: print >> f
    dump_h.extension = "h"

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


def codegen(name_expr, language, prefix, project="project", to_files=False, header=True, empty=True):
    """Write source code for the given expressions in the given language.

       Mandatory Arguments:
         name_expr  --  A single (name, expression) tuple or a list of
                        (name, expression) tuples. Each tuple corresponds to a
                        routine.  If the expression is an equality (an
                        insance of class Equality) the left hand side is considered
                        an output argument.
         language  --  A string that indicates the source code language. This
                       is case insensitive. For the moment, only 'C' is
                       supported.
         prefix  --  A prefix for the names of the files that contain the source
                     code. Proper (language dependent) suffixes will be
                     appended.

       Optional Arguments:
         project  --  A project name, used for making unique preprocessor
                      instructions. [DEFAULT="project"]
         to_files  --  When True, the code will be written to one or more files
                       with the given prefix, otherwise strings with the names
                       and contents of these files are returned. [DEFAULT=False]
         header  --  When True, a header is written on top of each source file.
                     [DEFAULT=True]
         empty  --  When True, empty lines are used to structure the code.
                    [DEFAULT=True]

       >>> from sympy import symbols
       >>> from sympy.utilities.codegen import codegen
       >>> from sympy.abc import x, y, z
       >>> [(c_name, c_code), (h_name, c_header)] = \\
       ...     codegen(("f", x+y*z), "C", "test", header=False, empty=False)
       >>> print c_name
       test.c
       >>> print c_code,
       #include "test.h"
       #include <math.h>
       double f(double x, double y, double z) {
         return x + y*z;
       }
       >>> print h_name
       test.h
       >>> print c_header,
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
    if isinstance(name_expr[0], basestring):
        # single tuple is given, turn it into a singleton list with a tuple.
        name_expr = [name_expr]

    for name, expr in name_expr:
        routines.append(Routine(name, expr))

    # Write the code.
    return code_gen.write(routines, prefix, to_files, header, empty)
