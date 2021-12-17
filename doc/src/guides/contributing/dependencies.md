# Dependencies

There are several packages that, when installed, can enable certain additional
SymPy functionality. Most users and contributors will not need to install any
of the packages mentioned below (except for the hard dependencies), unless
they intend to use the part(s) of SymPy that can use those packages.

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

## Hard Dependencies

SymPy only has one hard dependency, which is required for it to work: mpmath.

### mpmath

[mpmath](https://mpmath.org/) is a pure Python package for arbitrary precision
arithmetic. It is used under the hood whenever SymPy calculates the
floating-point value of a function, e.g., when using
[evalf](sympy.core.evalf.EvalfMixin.evalf).

SymPy cannot function with mpmath and will fail to import if it is not
installed. If you get an error like

```
ImportError: SymPy now depends on mpmath as an external library. See
https://docs.sympy.org/latest/install.html#mpmath for more information.
```

this means that you did not install mpmath correctly. [This
page](mpmath-install) explains how to install it.

Most methods of installing SymPy, such as the ones outlined in the
[installation](installation) guide, will install mpmath automatically. You
typically only need to install mpmath manually if you did not actually install
SymPy, e.g., if you are developing directly on SymPy in the git repository.

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

#### gmpy2

[gmpy2](https://gmpy2.readthedocs.io/en/latest/intro.html) is a Python wrapper
for the GMP multiple-precision library. It provides large integers that are
faster than the built-in Python `int`. When gmpy2 is installed, it is used
automatically by certain core functions that operate on integers, such as the
[polys](polys-docs). See {ref}`polys-domainsref` for more details.

### Interactive Use

#### IPython

#### Jupyter Notebook

### LaTeX Parsing

#### antlr-python-runtime

#### python-clang

### Logic

#### pycosat

#### pysat

### Plotting

#### matplotlib

#### pyglet

### lambdify

#### numpy

#### scipy

#### cupy

#### tensorflow

#### numexpr

### Code Generation

#### NumPy

#### Cython

#### Aesara

#### pymc3

#### llvmlite

#### Wurlitzer

[Wurlitzer](https://github.com/minrk/wurlitzer) is a Python package that
allows capturing output from C extensions. It is used by some of the tests in
the `sympy.codegen` submodule. It is only used by the test suite. It is not
used by any end-user functionality. If it is not installed, some tests will be
skipped.

### Optional SymEngine Backend

#### python-symengine

[SymEngine](https://symengine.org/) is a fast symbolic manipulation library,
written in C++. The SymEngine Python bindings may be used as an optional
backend for SymPy core. To do this, first install the SymEngine Python
bindings (with `pip install symengine` or `conda install -c conda-forge
python-symengine`) and run SymPy with the `USE_SYMENGINE=1` environment variable.

Presently, the SymEngine backend is only used by the
[sympy.physics.mechanics](classical_mechanics) and
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

### Experimental Rubi Integrator

#### MatchPy

[MatchPy](https://matchpy.readthedocs.io/en/latest/) is a library for doing
pattern matching. It is used in the experimental sympy.integrals.rubi module,
but presently, it is not used anywhere else in SymPy. SymPy and MatchPy are
able to interface with each other.



### Other

#### cloudpickle

The [cloudpickle](https://github.com/cloudpipe/cloudpickle) package can be
used to more effectively pickle SymPy objects than the built-in Python [pickle](https://docs.python.org/3/library/pickle.html).
Some tests in `sympy.utilities.tests.test_pickling.py` depend on cloudpickle
to run. It is not otherwise required for any SymPy function.

#### Sage


## Development Dependencies

Typical development on SymPy does not require any additional dependencies
beyond Python and mpmath.

### git

The [SymPy source code](https://github.com/sympy/sympy) uses the
[git](https://git-scm.com/) version control system. See
the [installation guide](installation-git) and [development
workflow](https://github.com/sympy/sympy/wiki/Development-workflow#set-up-git)
for instructions on how to get the development version of SymPy from git.

### Building the Documentation

Building the documentation requires several additional dependencies.
[This page](build-docs) outlines these dependencies and how to install
them. It is only necessary to install these dependencies if you are
contributing documentation to SymPy and want to check that the HTML or PDF
documentation renders correctly. Documentation for the development version of
SymPy is hosted online at https://docs.sympy.org/dev/index.html.

### Running the Benchmarks

The benchmarks for SymPy are hosted at
https://github.com/sympy/sympy_benchmarks. The README in that repository
explains how to run the benchmarks.

Note that the benchmarks are also run automatically on the [GitHub Actions
CI](https://github.com/sympy/sympy/actions), so it is not generally necessary
to run them yourself as a contributor unless you want to reproduce the
benchmarks results on your computer or add a new benchmark to the suite.

#### asv

[Airspeed Velocity](https://asv.readthedocs.io/en/stable/) is the package used
for running the benchmarks. Note that the package name that you install is
`asv`.

### Tests only dependencies

#### Cloudpickle
