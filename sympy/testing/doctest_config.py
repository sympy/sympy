"""The single source of truth for files that need special doctest treatment.

The list of files that are excluded from the ordinary doctest run used to be
maintained independently in three places: ``bin/doctest``,
``sympy/testing/runtests.py`` and ``bin/test_optional_dependencies.py``. Those
lists drifted apart, and worse, a list entry can silently stop matching
anything when the file it refers to is moved or deleted. The test runner then
simply has nothing to filter out, so CI stays green while a document is no
longer tested at all (or, in the other direction, a document that needs an
optional dependency is suddenly tested without it).

All such files are now declared here. A file is either

* in :data:`SKIPPED_FILES`, because it is never part of the ordinary doctest
  run -- the value explains why; or
* in :data:`OPTIONAL_FILES`, because its doctests need optional dependencies.
  Such a file is skipped by the ordinary doctest run when the dependencies are
  missing and is run by ``bin/test_optional_dependencies.py`` when they are
  installed.

Paths are relative to the root of the repository and always use ``/`` as the
separator, also on Windows. Use ``convert_to_native_paths`` from
``sympy.testing.runtests`` (or :func:`native_path` here) before comparing them
with real paths.
"""
from __future__ import annotations

import os

from sympy.external import import_module


# Files that are never part of the ordinary doctest run, mapped to the reason
# why they are not.
SKIPPED_FILES = {
    "doc/src/modules/plotting.rst":
        "generates live plots",
    "doc/src/explanation/modules/physics/mechanics/autolev_parser.rst":
        "needs the missing double_pendulum.al fixture",
    "sympy/conftest.py":
        "depends on pytest",
    "sympy/core/compatibility.py":
        "backwards compatibility shim, importing it triggers a deprecation warning",
    "sympy/core/trace.py":
        "backwards compatibility shim, importing it triggers a deprecation warning",
    "sympy/galgebra.py":
        "no longer part of SymPy",
    "sympy/parsing/autolev/_antlr/autolevlexer.py":
        "generated code",
    "sympy/parsing/autolev/_antlr/autolevlistener.py":
        "generated code",
    "sympy/parsing/autolev/_antlr/autolevparser.py":
        "generated code",
    "sympy/parsing/autolev/test-examples/pydy-example-repo/chaos_pendulum.py":
        "generated example",
    "sympy/parsing/autolev/test-examples/pydy-example-repo/double_pendulum.py":
        "generated example",
    "sympy/parsing/autolev/test-examples/pydy-example-repo/mass_spring_damper.py":
        "generated example",
    "sympy/parsing/autolev/test-examples/pydy-example-repo/non_min_pendulum.py":
        "generated example",
    "sympy/parsing/autolev/test-examples/ruletest1.py":
        "generated example",
    "sympy/parsing/autolev/test-examples/ruletest2.py":
        "generated example",
    "sympy/parsing/autolev/test-examples/ruletest3.py":
        "generated example",
    "sympy/parsing/autolev/test-examples/ruletest4.py":
        "generated example",
    "sympy/parsing/autolev/test-examples/ruletest5.py":
        "generated example",
    "sympy/parsing/autolev/test-examples/ruletest6.py":
        "generated example",
    "sympy/parsing/autolev/test-examples/ruletest7.py":
        "generated example",
    "sympy/parsing/autolev/test-examples/ruletest8.py":
        "generated example",
    "sympy/parsing/autolev/test-examples/ruletest9.py":
        "generated example",
    "sympy/parsing/autolev/test-examples/ruletest10.py":
        "generated example",
    "sympy/parsing/autolev/test-examples/ruletest11.py":
        "generated example",
    "sympy/parsing/autolev/test-examples/ruletest12.py":
        "generated example",
    "sympy/parsing/latex/_antlr/latexlexer.py":
        "generated code",
    "sympy/parsing/latex/_antlr/latexparser.py":
        "generated code",
    "sympy/plotting/pygletplot/__init__.py":
        "crashes on some systems",
    "sympy/plotting/pygletplot/plot.py":
        "crashes on some systems",
    "sympy/testing/randtest.py":
        "backwards compatibility shim, importing it triggers a deprecation warning",
    "sympy/this.py":
        "prints text",
    "sympy/utilities/autowrap.py":
        "disabled because of doctest failures",
    "sympy/utilities/pytest.py":
        "deprecated stub to be removed",
    "sympy/utilities/randtest.py":
        "deprecated stub to be removed",
    "sympy/utilities/runtests.py":
        "deprecated stub to be removed",
    "sympy/utilities/tmpfiles.py":
        "deprecated stub to be removed",
}

# Files whose doctests need optional dependencies. The file is skipped by the
# ordinary doctest run unless every listed dependency can be imported. The
# documentation files are then run by bin/test_optional_dependencies.py, which
# installs those dependencies.
OPTIONAL_FILES = {
    "doc/src/explanation/best-practices.md":
        ("numpy",),
    "doc/src/explanation/modules/physics/biomechanics/biomechanics.rst":
        ("numpy", "scipy", "matplotlib"),
    "doc/src/guides/solving/solve-numerically.md":
        ("numpy", "scipy"),
    "doc/src/guides/solving/solve-ode.md":
        ("numpy", "scipy"),
    "doc/src/modules/diffgeom.rst":
        ("numpy",),
    "doc/src/modules/numeric-computation.rst":
        ("numpy", "aesara", "cupy", "jax"),
    "doc/src/tutorials/physics/biomechanics/biomechanical-model-example.rst":
        ("numpy", "scipy", "matplotlib"),
    "sympy/parsing/autolev/__init__.py":
        ("antlr4",),
    "sympy/parsing/latex/_parse_latex_antlr.py":
        ("antlr4",),
    "sympy/parsing/sym_expr.py":
        ("lfortran",),
    "sympy/plotting/experimental_lambdify.py":
        ("numpy",),
    "sympy/plotting/plot_implicit.py":
        ("numpy",),
    "sympy/plotting/pygletplot":
        ("pyglet",),
    "sympy/printing/aesaracode.py":
        ("aesara",),
}

# Files that are skipped when running on CI even though their dependencies
# happen to be installed there.
ON_CI_SKIPPED_FILES = {
    "sympy/plotting/pygletplot":
        "crashes in CI environments",
}


def native_path(path):
    """Convert a path of this configuration to a normalized native path."""
    path = path.replace("/", os.sep).replace("\\", os.sep)
    return os.path.normcase(os.path.normpath(path))


def get_sympy_dir():
    """Return the root directory of the repository."""
    this_file = os.path.abspath(__file__)
    sympy_dir = os.path.join(os.path.dirname(this_file), "..", "..")
    return os.path.normcase(os.path.normpath(sympy_dir))


def _available(dependency):
    return import_module(dependency) is not None


def regular_doctest_blacklist(available=None):
    """Return the files to exclude from the ordinary doctest run.

    ``available`` optionally maps a dependency name to a boolean telling
    whether that dependency can be imported. Anything not in that mapping is
    looked up with ``sympy.external.import_module``.
    """
    available = dict(available or {})
    blacklist = list(SKIPPED_FILES)

    for path, dependencies in OPTIONAL_FILES.items():
        if not all(available.get(dep, _available(dep)) for dep in dependencies):
            blacklist.append(path)

    if os.getenv("CI", None):
        blacklist.extend(ON_CI_SKIPPED_FILES)

    return blacklist


def optional_doctest_files(available=None):
    """Return the documentation files to test with optional dependencies.

    These are the documentation files of :data:`OPTIONAL_FILES` whose
    dependencies are all installed. The ordinary doctest run skips them in that
    case, so ``bin/test_optional_dependencies.py`` has to test them instead.
    """
    available = dict(available or {})
    return [path for path, dependencies in OPTIONAL_FILES.items()
            if path.startswith("doc/")
            and all(available.get(dep, _available(dep)) for dep in dependencies)
            and path not in ON_CI_SKIPPED_FILES]


def validate_doctest_config(sympy_dir=None):
    """Return a list of problems found in this configuration.

    An empty list means that every configured path exists in the repository.
    A path that no longer exists is reported, because the test runner silently
    matches no file with it.
    """
    sympy_dir = get_sympy_dir() if sympy_dir is None else sympy_dir
    problems = []

    for path in SKIPPED_FILES:
        if not os.path.exists(os.path.join(sympy_dir, native_path(path))):
            problems.append(f"{path} does not exist")

    for path in list(OPTIONAL_FILES) + list(ON_CI_SKIPPED_FILES):
        if not os.path.exists(os.path.join(sympy_dir, native_path(path))):
            problems.append(f"{path} does not exist")

    for path in ON_CI_SKIPPED_FILES:
        if path not in OPTIONAL_FILES:
            problems.append(f"{path} is only skipped on CI but has no dependencies")

    for path, reason in SKIPPED_FILES.items():
        if not reason:
            problems.append(f"{path} has no reason for being skipped")
        if path in OPTIONAL_FILES:
            problems.append(f"{path} is both skipped and optional")

    return problems
