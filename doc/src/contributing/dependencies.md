# Dependencies

This page lists the hard and optional dependencies of SymPy.

There are several packages that, when installed, can enable certain additional
SymPy functionality. Most users and contributors will not need to install any
of the packages mentioned below (except for the hard dependencies), unless
they intend to use or contribute to the parts of SymPy that can use those
packages.

Every dependency listed below can be installed with conda via
[conda-forge](https://conda-forge.org/), and most can also be installed with
`pip`.

This page does not list packages which themselves depend on SymPy, only those
packages that SymPy depends on. An incomplete list of packages that depend on
SymPy can be found on the [main SymPy
webpage](https://www.sympy.org/en/index.html), and a more complete list can be
found on
[GitHub](https://github.com/sympy/sympy/network/dependents?dependent_type=PACKAGE)
or [libraries.io](https://libraries.io/pypi/sympy/dependents).

(hard-dependencies)=
## Hard Dependencies

SymPy only has one hard dependency, which is required for it to work: mpmath.

(dependencies-mpmath)=
- **mpmath**: [mpmath](https://mpmath.org/) is a pure Python package for
  arbitrary precision arithmetic. It is used under the hood whenever SymPy
  calculates the floating-point value of a function, e.g., when using
  [evalf](sympy.core.evalf.EvalfMixin.evalf).

  SymPy cannot function without mpmath and will fail to import if it is not
  installed. If you get an error like

  ```pytb
  ImportError: SymPy now depends on mpmath as an external library. See
  https://docs.sympy.org/latest/install.html#mpmath for more information.
  ```

  this means that you did not install mpmath correctly. [This
  page](mpmath-install) explains how to install it.

  Most methods of installing SymPy, such as the ones outlined in the
  [installation](installation) guide, will install mpmath automatically. You
  typically only need to install mpmath manually if you did not actually
  install SymPy, e.g., if you are developing directly on SymPy in the git
  repository.

(optional-dependencies)=
## Optional Dependencies

These dependencies are not required to use SymPy. The vast majority of SymPy
functions do not require them, however, a few functions such as plotting and
automatic wrapping of code generated functions require additional dependencies
to function.

Additionally, as a contributor, when running the SymPy tests, some tests will
be skipped if a dependency they require is not installed. The [GitHub Actions
CI](https://github.com/sympy/sympy/actions) which is run on every SymPy pull
request will automatically install these dependencies in the
"optional-dependencies" build, but you may wish to install them locally if you
are working on a part of SymPy that uses them.

### Recommended Optional Dependencies

These dependencies are not required for SymPy to function, but it is
recommended that all users install them if they can, as they will improve the
general performance of SymPy.

- **gmpy2**: [gmpy2](https://gmpy2.readthedocs.io/en/latest/) is a
  Python wrapper for the [GMP multiple-precision
  library](https://gmplib.org/). It provides large integers that are faster
  than the built-in Python `int`. When gmpy2 is installed, it is used
  automatically by certain core functions that operate on integers, such as
  the [polys](polys-docs). See {ref}`polys-domainsref` for more details. SymPy
  uses `gmpy2` automatically when it is installed. No further action is
  required to enable it.

  The polys themselves are used by many parts of SymPy, such as the
  integration algorithms, simplification algorithms like `collect()` and
  `factor()`, the matrices, and some parts of the core. Thus, installing
  `gmpy2` can speed up many parts of SymPy. It is not a required dependency of
  SymPy because it makes use of a non-Python library (GMP), which is also
  non-BSD licensed. However, we recommended all users who are able to to
  install `gmpy2` to get a better SymPy experience.

### Interactive Use

SymPy is designed to be used both interactively and as a library. When used
interactively, SymPy is able to interface with IPython and Jupyter notebooks.

- **IPython**: The {func}`~.init_session` function and `isympy` command will
  automatically start IPython if it is installed. In addition to the usual
  benefits of using [IPython](https://ipython.org/), this enables interactive
  plotting with matplotlib. Also some flags such as `auto_symbols` and
  `auto_int_to_Integer` will only work in IPython.

  The `IPython` package is required to run some of the tests in sympy/interactive.

- **Jupyter Notebook and Qt Console**: SymPy expressions automatically print
  using MathJax in the [Jupyter Notebook](https://jupyter.org/) and with LaTeX
  [Qt Console](https://qtconsole.readthedocs.io/en/stable/) (if
  [LaTeX](dependencies-latex) is installed).

### Printing

The {func}`~.preview` function automatically converts SymPy expressions into
images rendered with LaTeX. `preview()` can either save the image to a file or
show it with a viewer.

(dependencies-latex)=
- **LaTeX**: A $\mathrm{\LaTeX}$ distributions such as [TeXLive](https://tug.org/texlive/) or
[MiKTeX](https://miktex.org/) is required for {func}`~.preview` to function.

### Parsing

Several functions in the {mod}`sympy.parsing` submodule require external
dependencies to function. Note that not all parsers require external modules
at this time. The Python ({func}`~.parse_expr`), Mathematica
({func}`~.parse_mathematica`), and Maxima ({func}`~.parse_maxima`) parsers do not
require any external dependencies.

- **antlr-python-runtime**: [ANTLR](https://www.antlr.org/) can be used for the
  {func}`LaTeX parser <sympy.parsing.latex.parse_latex>`, and is used in the
  [Autolev](autolev_parser) parsers. They both require the ANTLR Python
  runtime to be installed. The package for this is called
  `antlr-python-runtime` with conda and `antlr4-python3-runtime` with pip.
  Also be aware that the version of the ANTLR Python runtime must match the
  version that was used to compile the LaTeX and Autolev parsers (4.10).

- **lark**: [Lark](https://lark-parser.readthedocs.io/en/stable/) can be used
  as an alternative backend for the {func}`LaTeX parser <sympy.parsing.latex.parse_latex>`.

- **Clang Python Bindings**: The C parser (`sympy.parsing.c.parse_c`) requires
  the Clang Python bindings. The package for this is called `python-clang`
  with conda and `clang` with pip.

- **lfortran**: The Fortran parser (in `sympy.parsing.fortran`) requires
  [LFortran](https://lfortran.org/).

### Logic

The {func}`~.satisfiable` function includes a pure Python implementation of
the DPLL satisfiability algorithm. But it can optionally use faster C SAT
solvers if they are installed. Note that `satisfiable()` is also used by
{func}`~.ask`.

- **pycosat**: [Pycosat](https://pypi.org/project/pycosat/) is used
  automatically if it is installed. The use of pycosat can be forced by using
  `satisfiable(algorithm='pycosat')`.

- **pysat**: [Pysat](https://pysathq.github.io/) is a library which wraps many
  SAT solvers. It can also be used as a backend to `satisfiable()`. Presently,
  only [Minisat](http://minisat.se/MiniSat.html) is implemented, using
  `satisfiable(algorithm=minisat22')`.

### Plotting

The {mod}`sympy.plotting.plot` module makes heavy use of external plotting
libraries to render plots. The primarily plotting module that is supported is
Matplotlib.

- **matplotlib**: Most plotting functionality requires the
  [Matplotlib](https://matplotlib.org/) plotting library. Without Matplotlib
  installed, most plotting functions will either fail or give rudimentary
  [text plots](textplot).

- **pyglet**: SymPy has a submodule {mod}`sympy.plotting.pygletplot` that can
  be used to interface with the [pyglet](https://pyglet.org/) module to do 2D
  and 3D plotting.

(dependencies-lambdify)=
### lambdify

{func}`~.lambdify` is a function that converts SymPy expressions into functions
that can be evaluated numerically using various libraries as backends.
`lambdify` is the primary vehicle by which users interface between SymPy and
these libraries. It is the standard way to convert a symbolic SymPy expression
into an evaluable numeric function.

In principle, `lambdify` can interface with any external library if the user
passes in an appropriate namespace dictionary as the third argument, but by
default, `lambdify` is aware of several popular numeric Python libraries.
These libraries are enabled as backends in `lambdify` with built-in
translations to convert SymPy expressions into the appropriate functions for
those libraries.

- **NumPy**: By default, if it is installed, `lambdify` creates functions
  using [NumPy](https://numpy.org/) (if NumPy is not installed, `lambdify`
  produces functions using the standard library
  [math](https://docs.python.org/3/library/math.html) module, although this
  behavior is primarily provided for backwards compatibility).

- **SciPy**: If [SciPy](https://scipy.org/) is installed, `lambdify` will use
  it automatically. SciPy is needed to lambdify certain [special
  functions](https://docs.scipy.org/doc/scipy/reference/special.html) that are
  not included in NumPy.

- **CuPy**: [CuPy](https://cupy.dev/) is a library that provides a NumPy
  compatible interface for CUDA GPUs. `lambdify` can produce CuPy compatible
  functions using `lambdify(modules='cupy')`.

- **Jax**: [JAX](https://github.com/google/jax) is a library that uses XLA
  to compile and run NumPy programs on GPUs and TPUs. `lambdify` can produce
  JAX compatibly functions using `lambdify(modules='jax')`.

- **TensorFlow**: [TensorFlow](https://www.tensorflow.org/) is a popular
  machine learning library. `lambdify` can produce TensorFlow compatible
  functions using `lambdify(modules='tensorflow')`.

- **NumExpr**: [NumExpr](https://github.com/pydata/numexpr) is a fast
  numerical expression evaluator for NumPy. `lambdify` can produce NumExpr
  compatible functions using `lambdify(modules='numexpr')`.

- **mpmath**: `lambdify` can also produce mpmath compatible functions. Note
  that mpmath is already a [required dependency](dependencies-mpmath) of
  SymPy. This functionality is useful for converting a SymPy expression to a
  function for use with pure mpmath.

### Code Generation

SymPy can [generate code](codegen_prose) for a large number of languages by
converting SymPy expressions into valid code for those languages. It also has
functionality for some languages to automatically compile and run the code.

Note that the dependencies below are **not** a list of supported languages
that SymPy can generate code for. Rather it is a list of packages that SymPy
can interface with in some way. For most languages that SymPy supports code
generation, it simply generates a string representing the code for that
language, so no dependency on that language is required to use the code
generation functionality. A dependency is typically only required for features
that automatically take the generated code and compile it to a function that
can be used within Python. Note that {func}`~.lambdify` is a special case of
this, but its dependencies are listed [above](dependencies-lambdify).

#### Autowrap

- **NumPy**: [NumPy](https://numpy.org/) and, optionally, its subpackage
  [f2py](https://numpy.org/doc/stable/f2py/), can be used to generate Python
  functions using the {func}`~.autowrap` or {func}`~.ufuncify`
  functions.

- **Cython**: [Cython](https://cython.org/) can be used as a backend for
  {func}`~.autowrap` or {func}`~.ufuncify`. Cython is also used in some of the
  `sympy.codegen` tests to compile some examples.

(dependencies-compilers)=
- **Compilers**: {func}`~.autowrap`, {func}`~.ufuncify`, and
  related functions rely on a compiler to compile the generated code to a
  function. Most standard C, C++, and Fortran compilers are supported,
  including [Clang/LLVM](https://clang.llvm.org/),
  [GCC](https://gcc.gnu.org/), and
  [ifort](https://en.wikipedia.org/wiki/Intel_Fortran_Compiler).

#### Code Printers

Most code printers generate Python strings, and therefore do not require the
given library or language compiler as a dependency. However, a few code
printers generate Python functions instead of strings:

- **Aesara**: The {mod}`sympy.printing.aesaracode` module contains functions
  to convert SymPy expressions into a functions using the
  [Aeseara](https://aesara.readthedocs.io/en/latest) (previously Theano)
  library. The Aesara code generation functions return Aesara graph objects.

- **llvmlite**: The `sympy.printing.llvmjitcode` module supports generating LLVM Jit
  from a SymPy expression. The functions make use of
  [llvmlite](https://llvmlite.readthedocs.io/en/latest/), a Python wrapper
  around [LLVM](https://llvm.org/). The `llvm_callable()` function
  generates callable functions.

- **TensorFlow**: The `sympy.printing.tensorflow` module supports generating
  functions using the [TensorFlow](https://www.tensorflow.org/), a popular
  machine learning library. Unlike the above two examples, `tensorflow_code()`
  function **does** generate Python strings. However, `tensorflow` is imported
  if available in order to automatically detect the TensorFlow version. If it
  is not installed, the `tensorflow_code()` function assumes the latest
  supported version of TensorFlow.

#### Testing-Only Dependencies

- **Wurlitzer**: [Wurlitzer](https://github.com/minrk/wurlitzer) is a Python
  package that allows capturing output from C extensions. It is used by some
  of the tests in the `sympy.codegen` submodule. It is only used by the test
  suite. It is not used by any end-user functionality. If it is not installed,
  some tests will be skipped.

- **Cython**: [Cython](https://cython.org/) is also used in
  some of the `sympy.codegen` tests to compile some examples.


- **Compilers**: The various [compilers](dependencies-compilers) mentioned
  above are used in some of the codegen and autowrap tests if they are
  installed.

### Statistics

The {func}`sympy.stats.sample` function uses an external library to produce
samples from the given distribution. At least one of the following libraries
is required to use the sampling functionality of `sympy.stats`.

- **SciPy**: `sample(library='scipy')` is the default. This uses [scipy.stats](https://docs.scipy.org/doc/scipy/reference/stats.html).

- **NumPy**: `sample(library='numpy')` uses the [NumPy
  random module](https://numpy.org/doc/stable/reference/random/index.html).

- **pymc**: `sample(library='pymc')` uses
  [PyMC](https://www.pymc.io/) to do sampling.

### Optional SymEngine Backend

- **python-symengine**: [SymEngine](https://symengine.org/) is a fast symbolic
  manipulation library, written in C++. The SymEngine Python bindings may be
  used as an optional backend for SymPy core. To do this, first install the
  SymEngine Python bindings (with `pip install symengine` or `conda install -c
  conda-forge python-symengine`) and run SymPy with the `USE_SYMENGINE=1`
  environment variable.

  Presently, the SymEngine backend is only used by the
  [sympy.physics.mechanics](physics_mechanics) and
  [sympy.liealgebras](lie-algebras) modules, although you can also interface
  with SymPy's SymEngine backend directly by importing things from
  `sympy.core.backend`:

  ```
  >>> from sympy.core.backend import Symbol
  >>> # This will create a SymEngine Symbol object if the USE_SYMENGINE
  >>> # environment variable is configured. Otherwise it will be an ordinary
  >>> # SymPy Symbol object.
  >>> x = Symbol('x')
  ```

  SymEngine backend support is still experimental, so certain SymPy functions
  may not work correctly when it is enabled.

### Sage

[Sage](https://www.sagemath.org/) is an open source mathematics software that
incorporates a large number of open source mathematics libraries. SymPy is one
of the libraries used by Sage.

Most of the code that interfaces between SymPy and Sage is in Sage itself, but
a few `_sage_` methods in SymPy that do some very basic setting up of the
Sage/SymPy wrappers. These methods should typically only be called by Sage
itself.

## Development Dependencies

Typical development on SymPy does not require any additional dependencies
beyond Python and mpmath.

### Getting the Source Code

- **git**: The [SymPy source code](https://github.com/sympy/sympy) uses the
  [git](https://git-scm.com/) version control system. See the [installation
  guide](installation-git) and the [contributor guide](devsetup)
  for instructions on how to get the development version of SymPy from git.

### Running the Tests

The base SymPy tests do not require any additional dependencies, however most
of the above dependencies may be required for some tests to run. Tests that
depend on optional dependencies should be skipped when they are not installed,
either by using the `sympy.testing.pytest.skip()` function or by setting `skip
= True` to skip the entire test file. Optional modules in tests and SymPy
library code should be imported with `import_module()`.

- **pytest**: [Pytest](https://docs.pytest.org/en/latest/) is not a required dependency
  for the SymPy test suite. SymPy has its own test runner, which can be
  accessed via the `bin/test` script in the SymPy source directory or the
  {func}`~.test` function.

  However, if you prefer to use pytest, you can use it to run the tests
  instead of the SymPy test runner. Tests in SymPy should use the wrappers in
  {mod}`sympy.testing.pytest` instead of using pytest functions directly.

- **Cloudpickle**: The [cloudpickle](https://github.com/cloudpipe/cloudpickle)
  package can be used to more effectively pickle SymPy objects than the
  built-in Python [pickle](https://docs.python.org/3/library/pickle.html).
  Some tests in `sympy.utilities.tests.test_pickling.py` depend on cloudpickle
  to run. It is not otherwise required for any SymPy function.

- **hypothesis**: [Hypothesis](https://github.com/HypothesisWorks/hypothesis/tree/master)
  is a required dependency for the SymPy test suit.

### Building the Documentation

Building the documentation requires several additional dependencies. [This
page](build-the-documentation) outlines these dependencies and how to install
them. It is only necessary to install these dependencies if you are
contributing documentation to SymPy and want to check that the HTML or PDF
documentation renders correctly. If you only want to view the documentation
for the development version of SymPy, development builds of the docs are
hosted online at https://docs.sympy.org/dev/index.html.

### Running the Benchmarks

The benchmarks for SymPy are hosted at
https://github.com/sympy/sympy_benchmarks. The
[README](https://github.com/sympy/sympy_benchmarks#readme) in that repository
explains how to run the benchmarks.

Note that the benchmarks are also run automatically on the [GitHub Actions
CI](https://github.com/sympy/sympy/actions), so it is generally not necessary
to run them yourself as a contributor unless you want to reproduce the
benchmarks results on your computer or add a new benchmark to the suite.

- **asv**: [Airspeed Velocity](https://asv.readthedocs.io/en/stable/) is the
  package used for running the benchmarks. Note that the package name that you
  install is called `asv`.
