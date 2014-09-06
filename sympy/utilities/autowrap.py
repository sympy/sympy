"""Module for compiling codegen output, and wrap the binary for use in
python.

.. note:: To use the autowrap module it must first be imported

   >>> from sympy.utilities.autowrap import autowrap

This module provides a common interface for different external backends, such
as f2py, fwrap, Cython, SWIG(?) etc. (Currently only f2py and Cython are
implemented) The goal is to provide access to compiled binaries of acceptable
performance with a one-button user interface, i.e.

    >>> from sympy.abc import x,y
    >>> expr = ((x - y)**(25)).expand()
    >>> binary_callable = autowrap(expr)
    >>> binary_callable(1, 2)
    -1.0

The callable returned from autowrap() is a binary python function, not a
SymPy object.  If it is desired to use the compiled function in symbolic
expressions, it is better to use binary_function() which returns a SymPy
Function object.  The binary callable is attached as the _imp_ attribute and
invoked when a numerical evaluation is requested with evalf(), or with
lambdify().

    >>> from sympy.utilities.autowrap import binary_function
    >>> f = binary_function('f', expr)
    >>> 2*f(x, y) + y
    y + 2*f(x, y)
    >>> (2*f(x, y) + y).evalf(2, subs={x: 1, y:2})
    0.e-110

The idea is that a SymPy user will primarily be interested in working with
mathematical expressions, and should not have to learn details about wrapping
tools in order to evaluate expressions numerically, even if they are
computationally expensive.

When is this useful?

    1) For computations on large arrays, Python iterations may be too slow,
       and depending on the mathematical expression, it may be difficult to
       exploit the advanced index operations provided by NumPy.

    2) For *really* long expressions that will be called repeatedly, the
       compiled binary should be significantly faster than SymPy's .evalf()

    3) If you are generating code with the codegen utility in order to use
       it in another project, the automatic python wrappers let you test the
       binaries immediately from within SymPy.

    4) To create customized ufuncs for use with numpy arrays.
       See *ufuncify*.

When is this module NOT the best approach?

    1) If you are really concerned about speed or memory optimizations,
       you will probably get better results by working directly with the
       wrapper tools and the low level code.  However, the files generated
       by this utility may provide a useful starting point and reference
       code. Temporary files will be left intact if you supply the keyword
       tempdir="path/to/files/".

    2) If the array computation can be handled easily by numpy, and you
       don't need the binaries for another project.

"""

from __future__ import print_function, division

_doctest_depends_on = { 'exe': ('f2py', 'gfortran'), 'modules': ('numpy',)}

import sys
import os
import shutil
import tempfile
from subprocess import STDOUT, CalledProcessError

from sympy.core.cache import cacheit
from sympy.core.compatibility import check_output
from sympy.utilities.codegen import (
    get_code_generator, Routine, OutputArgument, InOutArgument, InputArgument,
    CodeGenArgumentListError, Result, ResultBase
)
from sympy.utilities.lambdify import implemented_function
from sympy.utilities.decorator import doctest_depends_on
from sympy import C


class CodeWrapError(Exception):
    pass


class CodeWrapper:
    """Base Class for code wrappers"""
    _filename = "wrapped_code"
    _module_basename = "wrapper_module"
    _module_counter = 0

    @property
    def filename(self):
        return "%s_%s" % (self._filename, CodeWrapper._module_counter)

    @property
    def module_name(self):
        return "%s_%s" % (self._module_basename, CodeWrapper._module_counter)

    def __init__(self, generator, filepath=None, flags=[], verbose=False):
        """
        generator -- the code generator to use
        """
        self.generator = generator
        self.filepath = filepath
        self.flags = flags
        self.quiet = not verbose

    @property
    def include_header(self):
        return bool(self.filepath)

    @property
    def include_empty(self):
        return bool(self.filepath)

    def _generate_code(self, main_routine, routines):
        routines.append(main_routine)
        self.generator.write(
            routines, self.filename, True, self.include_header,
            self.include_empty)

    def wrap_code(self, routine, helpers=[]):

        workdir = self.filepath or tempfile.mkdtemp("_sympy_compile")
        if not os.access(workdir, os.F_OK):
            os.mkdir(workdir)
        oldwork = os.getcwd()
        os.chdir(workdir)
        try:
            sys.path.append(workdir)
            self._generate_code(routine, helpers)
            self._prepare_files(routine)
            self._process_files(routine)
            mod = __import__(self.module_name)
        finally:
            sys.path.remove(workdir)
            CodeWrapper._module_counter += 1
            os.chdir(oldwork)
            if not self.filepath:
                shutil.rmtree(workdir)

        return self._get_wrapped_function(mod)

    def _process_files(self, routine):
        command = self.command
        command.extend(self.flags)
        try:
            retoutput = check_output(command, stderr=STDOUT)
        except CalledProcessError as e:
            raise CodeWrapError(
                "Error while executing command: %s. Command output is:\n%s" % (
                    " ".join(command), e.output.decode()))
        if not self.quiet:
            print(retoutput)


class DummyWrapper(CodeWrapper):
    """Class used for testing independent of backends """

    template = """# dummy module for testing of SymPy
def %(name)s():
    return "%(expr)s"
%(name)s.args = "%(args)s"
%(name)s.returns = "%(retvals)s"
"""

    def _prepare_files(self, routine):
        return

    def _generate_code(self, routine, helpers):
        with open('%s.py' % self.module_name, 'w') as f:
            printed = ", ".join(
                [str(res.expr) for res in routine.result_variables])
            # convert OutputArguments to return value like f2py
            args = filter(lambda x: not isinstance(
                x, OutputArgument), routine.arguments)
            retvals = []
            for val in routine.result_variables:
                if isinstance(val, Result):
                    retvals.append('nameless')
                else:
                    retvals.append(val.result_var)

            print(DummyWrapper.template % {
                'name': routine.name,
                'expr': printed,
                'args': ", ".join([str(a.name) for a in args]),
                'retvals': ", ".join([str(val) for val in retvals])
            }, end="", file=f)

    def _process_files(self, routine):
        return

    @classmethod
    def _get_wrapped_function(cls, mod):
        return mod.autofunc


class CythonCodeWrapper(CodeWrapper):
    """Wrapper that uses Cython"""

    setup_template = (
            "from distutils.core import setup\n"
            "from distutils.extension import Extension\n"
            "from Cython.Distutils import build_ext\n"
            "\n"
            "setup(\n"
            "    cmdclass = {{'build_ext': build_ext}},\n"
            "    ext_modules = [Extension({ext_args}, extra_compile_args=['-std=c99'])]\n"
            "        )")

    pyx_imports = (
            "import numpy as np\n"
            "cimport numpy as np\n\n")

    pyx_header = (
            "cdef extern from '{header_file}.h':\n"
            "    {prototype}\n\n")

    pyx_func = (
            "def {name}_c({arg_string}):\n"
            "\n"
            "{declarations}"
            "{body}")

    _need_numpy = False

    @property
    def command(self):
        command = [sys.executable, "setup.py", "build_ext", "--inplace"]
        return command

    def _prepare_files(self, routine):
        pyxfilename = self.module_name + '.pyx'
        codefilename = "%s.%s" % (self.filename, self.generator.code_extension)

        # pyx
        with open(pyxfilename, 'w') as f:
            self.dump_pyx([routine], f, self.filename)

        # setup.py
        ext_args = [repr(self.module_name), repr([pyxfilename, codefilename])]
        with open('setup.py', 'w') as f:
            f.write(self.setup_template.format(ext_args=", ".join(ext_args)))

    @classmethod
    def _get_wrapped_function(cls, mod):
        return mod.autofunc_c

    def dump_pyx(self, routines, f, prefix):
        """Write a Cython file with python wrappers

        This file contains all the definitions of the routines in c code and
        refers to the header file.

        Arguments
        ---------
        routines
            List of Routine instances
        f
            File-like object to write the file to
        prefix
            The filename prefix, used to refer to the proper header file.
            Only the basename of the prefix is used.
        """
        headers = []
        functions = []
        for routine in routines:
            prototype = self.generator.get_prototype(routine)

            # C Function Header Import
            headers.append(self.pyx_header.format(header_file=prefix,
                    prototype=prototype))

            # Partition the C function arguments into categories
            py_rets, py_args, py_loc, py_inf = self._partition_args(routine.arguments)

            # Function prototype
            name = routine.name
            arg_string = ", ".join(self._prototype_arg(arg) for arg in py_args)

            # Local Declarations
            local_decs = []
            for arg, val in py_inf.items():
                proto = self._prototype_arg(arg)
                mat, ind = val
                local_decs.append("    cdef {0} = {1}.shape[{2}]".format(proto, mat, ind))
            local_decs.extend(["    cdef {0}".format(self._declare_arg(a)) for a in py_loc])
            declarations = "\n".join(local_decs)
            if declarations:
                declarations = declarations + "\n"

            # Function Body
            args_c = ", ".join([self._call_arg(a) for a in routine.arguments])
            rets = ", ".join([str(r.name) for r in py_rets])
            if routine.results:
                body = '    return %s(%s)' % (routine.name, args_c)
                if rets:
                    body = body + ', ' + rets
            else:
                body = '    %s(%s)\n' % (routine.name, args_c)
                body = body + '    return ' + rets

            functions.append(self.pyx_func.format(name=name, arg_string=arg_string,
                    declarations=declarations, body=body))

        # Write text to file
        if self._need_numpy:
            # Only import numpy if required
            f.write(self.pyx_imports)
        f.write('\n'.join(headers))
        f.write('\n'.join(functions))

    def _partition_args(self, args):
        """Group function arguments into categories."""
        py_args = []
        py_returns = []
        py_locals = []
        py_inferred = {}
        for arg in args:
            if isinstance(arg, OutputArgument):
                py_returns.append(arg)
                py_locals.append(arg)
            elif isinstance(arg, InOutArgument):
                py_returns.append(arg)
                py_args.append(arg)
            else:
                py_args.append(arg)
        # Find arguments that are array dimensions. These can be inferred
        # locally in the Cython code.
            if isinstance(arg, (InputArgument, InOutArgument)) and arg.dimensions:
                dims = [d[1] + 1 for d in arg.dimensions]
                sym_dims = [(i, d) for (i, d) in enumerate(dims) if isinstance(d, C.Symbol)]
                for (i, d) in sym_dims:
                    py_inferred[d] = (arg.name, i)
        for arg in args:
            if arg.name in py_inferred:
                py_inferred[arg] = py_inferred.pop(arg.name)
        # Filter inferred arguments from py_args
        py_args = [a for a in py_args if not a in py_inferred]
        return py_returns, py_args, py_locals, py_inferred

    def _prototype_arg(self, arg):
        mat_dec = "np.ndarray[{mtype}, ndim={ndim}] {name}"
        np_types = {'double': 'np.double_t',
                    'int': 'np.int_t'}
        t = arg.get_datatype('c')
        if arg.dimensions:
            self._need_numpy = True
            ndim = len(arg.dimensions)
            mtype = np_types[t]
            return mat_dec.format(mtype=mtype, ndim=ndim, name=arg.name)
        else:
            return "%s %s" % (t, str(arg.name))

    def _declare_arg(self, arg):
        proto = self._prototype_arg(arg)
        if arg.dimensions:
            shape = '(' + ','.join(str(i[1] + 1) for i in arg.dimensions) + ')'
            return proto + " = np.empty({shape})".format(shape=shape)
        else:
            return proto + " = 0"

    def _call_arg(self, arg):
        if arg.dimensions:
            t = arg.get_datatype('c')
            return "<{0}*> {1}.data".format(t, arg.name)
        elif isinstance(arg, ResultBase):
            return "&{0}".format(arg.name)
        else:
            return str(arg.name)


class F2PyCodeWrapper(CodeWrapper):
    """Wrapper that uses f2py"""

    @property
    def command(self):
        filename = self.filename + '.' + self.generator.code_extension
        args = ['-c', '-m', self.module_name, filename]
        command = [sys.executable, "-c", "import numpy.f2py as f2py2e;f2py2e.main()"]+args
        return command

    def _prepare_files(self, routine):
        pass

    @classmethod
    def _get_wrapped_function(cls, mod):
        return mod.autofunc


def _get_code_wrapper_class(backend):
    wrappers = { 'F2PY': F2PyCodeWrapper, 'CYTHON': CythonCodeWrapper,
        'DUMMY': DummyWrapper}
    return wrappers[backend.upper()]

@cacheit
@doctest_depends_on(exe=('f2py', 'gfortran'), modules=('numpy',))
def autowrap(
    expr, language='F95', backend='f2py', tempdir=None, args=None, flags=(),
        verbose=False, helpers=()):
    """Generates python callable binaries based on the math expression.

    expr
        The SymPy expression that should be wrapped as a binary routine

    :Optional arguments:

    language
        The programming language to use, currently 'C' or 'F95'
    backend
        The wrapper backend to use, currently f2py or Cython
    tempdir
        Path to directory for temporary files.  If this argument is supplied,
        the generated code and the wrapper input files are left intact in the
        specified path.
    args
        Sequence of the formal parameters of the generated code, if ommited the
        function signature is determined by the code generator.
    flags
        Tuple of additional option flags that will be passed to the backend
    verbose
        If True, autowrap will not mute the command line backends.  This can be
        helpful for debugging.
    helpers
        Used to define auxillary expressions needed for the main expr.  If the
        main expression need to do call a specialized function it should be put
        in the ``helpers`` tuple.  Autowrap will then make sure that the
        compiled main expression can link to the helper routine.  Items should
        also be tuples with (<funtion_name>, <sympy_expression>, <arguments>).
        It is mandatory to supply an argument sequence to helper routines.

    >>> from sympy.abc import x, y, z
    >>> from sympy.utilities.autowrap import autowrap
    >>> expr = ((x - y + z)**(13)).expand()
    >>> binary_func = autowrap(expr)
    >>> binary_func(1, 4, 2)
    -1.0

    """

    code_generator = get_code_generator(language, "autowrap")
    CodeWrapperClass = _get_code_wrapper_class(backend)
    code_wrapper = CodeWrapperClass(code_generator, tempdir, flags, verbose)
    try:
        routine = Routine('autofunc', expr, args)
    except CodeGenArgumentListError as e:
        # if all missing arguments are for pure output, we simply attach them
        # at the end and try again, because the wrappers will silently convert
        # them to return values anyway.
        new_args = []
        for missing in e.missing_args:
            if not isinstance(missing, OutputArgument):
                raise
            new_args.append(missing.name)
        routine = Routine('autofunc', expr, args + new_args)

    helps = []
    for name, expr, args in helpers:
        helps.append(Routine(name, expr, args))

    return code_wrapper.wrap_code(routine, helpers=helps)


@doctest_depends_on(exe=('f2py', 'gfortran'), modules=('numpy',))
def binary_function(symfunc, expr, **kwargs):
    """Returns a sympy function with expr as binary implementation

    This is a convenience function that automates the steps needed to
    autowrap the SymPy expression and attaching it to a Function object
    with implemented_function().

    >>> from sympy.abc import x, y
    >>> from sympy.utilities.autowrap import binary_function
    >>> expr = ((x - y)**(25)).expand()
    >>> f = binary_function('f', expr)
    >>> type(f)
    <class 'sympy.core.function.UndefinedFunction'>
    >>> 2*f(x, y)
    2*f(x, y)
    >>> f(x, y).evalf(2, subs={x: 1, y: 2})
    -1.0
    """
    binary = autowrap(expr, **kwargs)
    return implemented_function(symfunc, binary)

@doctest_depends_on(exe=('f2py', 'gfortran'), modules=('numpy',))
def ufuncify(args, expr, **kwargs):
    """
    Generates a binary ufunc-like lambda function for numpy arrays

    ``args``
        Either a Symbol or a tuple of symbols. Specifies the argument sequence
        for the ufunc-like function.

    ``expr``
        A SymPy expression that defines the element wise operation

    ``kwargs``
        Optional keyword arguments are forwarded to autowrap().

    The returned function can only act on one array at a time, as only the
    first argument accept arrays as input.

    .. Note:: a *proper* numpy ufunc is required to support broadcasting, type
       casting and more.  The function returned here, may not qualify for
       numpy's definition of a ufunc.  That why we use the term ufunc-like.

    References
    ==========
    [1] http://docs.scipy.org/doc/numpy/reference/ufuncs.html

    Examples
    ========

    >>> from sympy.utilities.autowrap import ufuncify
    >>> from sympy.abc import x, y
    >>> import numpy as np
    >>> f = ufuncify([x, y], y + x**2)
    >>> f([1, 2, 3], 2)
    [ 3.  6.  11.]
    >>> a = f(np.arange(5), 3)
    >>> isinstance(a, np.ndarray)
    True
    >>> print a
    [ 3. 4. 7. 12. 19.]

    """
    y = C.IndexedBase(C.Dummy('y'))
    x = C.IndexedBase(C.Dummy('x'))
    m = C.Dummy('m', integer=True)
    i = C.Dummy('i', integer=True)
    i = C.Idx(i, m)
    l = C.Lambda(args, expr)
    f = implemented_function('f', l)

    if isinstance(args, C.Symbol):
        args = [args]
    else:
        args = list(args)

    # ensure correct order of arguments
    kwargs['args'] = [y, x] + args[1:] + [m]

    # first argument accepts an array
    args[0] = x[i]
    return autowrap(C.Equality(y[i], f(*args)), **kwargs)
