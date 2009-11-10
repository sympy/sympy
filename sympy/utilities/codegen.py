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
- Isolate constants and evaluate them beforehand in double precission
- Common Subexpression Elimination
- User defined comments in the generated code
- Fortran 77
- Fortran 90
- C++
- Python
- ...
"""


from sympy.core.symbol import Symbol
from sympy.core.basic import Basic
from sympy.utilities.iterables import postorder_traversal
from sympy.printing.ccode import ccode

from StringIO import StringIO
import sympy, os


__all__ = [
    # description of routines
    "Routine", "DataType", "default_datatypes", "get_default_datatype",
    "Argument", "InputArgument", "Result",
    # routines -> code
    "CodeGen", "CCodeGen",
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
    def __init__(self, name, arguments, results):
        """Initialize a Routine instance.

           Arguments:
             name  --  A string with the name of this routine in the generated
                       code
             arguments  --  A list of all input arguments. (TODO: add output
                            and input/output arguments). See InputArgument class.
             results  -- The return values. See Result class.
        """
        if len(results) == 0:
            # This condition will be updated when we support output arguments.
            raise ValueError("At least one result is required.")
        self.name = name
        self.arguments = arguments
        self.results = results


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

    def __init__(self, name, datatype):
        self.name = name
        self.datatype = datatype


class InputArgument(Argument):
    """A (scalar) input argument."""
    def __init__(self, symbol, datatype=None):
        """Initialize a (scalar) input argument.

           The second argument is optional. When not given, the data type will
           be guessed based on the assumptions on the symbol argument.
        """
        if not isinstance(symbol, Symbol):
            raise TypeError("The first argument must be a sympy symbol.")
        if datatype is None:
            datatype = get_default_datatype(symbol)
        elif not isinstance(datatype, DataType):
            raise TypeError("The (optional) second argument must be an instance of the DataType class.")
        self.symbol = symbol
        Argument.__init__(self, self.symbol.name, datatype)


class Result(object):
    """An expression for a return value.

       The name result is used to avoid conflicts with the reserved word
       'return' in the python language. It is also shorter than ReturnValue.
    """
    def __init__(self, expr, datatype=None):
        """Initialize a (scalar) return value.

           The second argument is optional. When not given, the data type will
           be guessed based on the assumptions on the expression argument.
        """
        if not isinstance(expr, Basic):
            raise TypeError("The first argument must be a sympy expression.")
        if datatype is None:
            datatype = get_default_datatype(expr)
        elif not isinstance(datatype, DataType):
            raise TypeError("The (optional) second argument must be an instance of the DataType class.")
        self.expr = expr
        self.datatype = datatype


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


class CodeGenError(Exception):
    pass


header_comment = """Code generated with sympy %(version)s

See http://www.sympy.org/ for more information.

This file is part of '%(project)s'
"""


class CCodeGen(CodeGen):
    def _dump_header(self, f):
        """Writes a common header for the .c and the .h file."""
        print >> f, "/****************************************************************************** "
        tmp = header_comment % {"version": sympy.__version__, "project": self.project}
        for line in tmp.splitlines():
            print >> f, " *%s* " % line.center(76)
        print >> f, " ******************************************************************************/"

    def get_prototype_result(self, routine):
        """Returns a string for the function prototype for the given routine and
           a single result object, which can be None.

           If the routine has multiple result objects, an CodeGenError is
           raised.

           See: http://en.wikipedia.org/wiki/Function_prototype
        """
        prototype = []
        if len(routine.results) > 1:
            raise CodeGenError("C only supports a single or no return value.")
        elif len(routine.results) == 1:
            result = routine.results[0]
            prototype.append(result.datatype.cname)
        else:
            result = None
            prototype.append("void")
        # name of the routine + arguments + curly opening brackets
        prototype.append("%s(%s)" % (
            routine.name,
            ", ".join("%s %s" % (arg.datatype.cname, arg.name) for arg in routine.arguments)
        ))
        return " ".join(prototype), result

    def dump_c(self, routines, f, prefix, header=True, empty=True):
        """Write the C code file.

           This file contains all the definitions of the routines in c code and
           refers to the header file.

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
        if header:
            self._dump_header(f)
        if empty: print >> f
        print >> f, "#include \"%s.h\"" % os.path.basename(prefix)
        print >> f, "#include <math.h>"
        if empty: print >> f
        for routine in routines:
            # function definitions.
            prototype, result = self.get_prototype_result(routine)
            print >> f, "%s {" % prototype
            # return value
            if result is not None:
                print >> f, "  return %s;" % ccode(result.expr)
            # curly closing brackets
            print >> f, "}"
            if empty: print >> f
        if empty: print >> f
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
            self._dump_header(f)
        guard_name = "%s__%s__H" % (self.project.replace(" ", "_").upper(), prefix.replace("/", "_").upper())
        # include guards
        if empty: print >> f
        print >> f, "#ifndef %s" % guard_name
        print >> f, "#define %s" % guard_name
        if empty: print >> f
        # declaration of the function prototypes
        for routine in routines:
            prototype, result = self.get_prototype_result(routine)
            print >> f, "%s;" % prototype
        # end if include guards
        if empty: print >> f
        print >> f, "#endif"
        if empty: print >> f
    dump_h.extension = "h"

    # This list of dump functions is used by CodeGen.write to know which dump
    # functions it has to call.
    dump_fns = [dump_c, dump_h]


#
# Friendly functions
#


def codegen(name_expr, language, prefix, project="project", to_files=False, header=True, empty=True):
    """Write source code for the given expressions in the given language.

       Mandatory Arguments:
         name_expr  --  A single (name, expression) tuple or a list of
                        (name, expression) tuples. Each tuple corresponds to a
                        routine
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
    CodeGenClass = {"C": CCodeGen}.get(language.upper())
    if CodeGenClass is None:
        raise ValueError("Language '%s' is not supported." % language)
    code_gen = CodeGenClass(project)

    # Construct the routines based on the name_expression pairs.
    #  mainly the input arguments require some work
    routines = []
    if isinstance(name_expr[0], basestring):
        # single tuple is given, turn it into a singleton list with a tuple.
        name_expr = [name_expr]
    for name, expr in name_expr:
        symbols = set([])
        for sub in postorder_traversal(expr):
            if isinstance(sub, Symbol):
                symbols.add(sub)
        routines.append(Routine(name, [InputArgument(symbol) for symbol in sorted(symbols)], [Result(expr)]))

    # Write the code.
    return code_gen.write(routines, prefix, to_files, header, empty)

