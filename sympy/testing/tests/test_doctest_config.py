from __future__ import annotations

import os

from sympy.testing.doctest_config import (
    OPTIONAL_FILES,
    SKIPPED_FILES,
    get_sympy_dir,
    native_path,
    optional_doctest_files,
    regular_doctest_blacklist,
    validate_doctest_config,
)


DEPENDENCIES = ["numpy", "scipy", "matplotlib", "aesara", "cupy", "jax",
                "antlr4", "lfortran", "pyglet"]

BIOMECHANICS_DOC = ("doc/src/tutorials/physics/biomechanics/"
                    "biomechanical-model-example.rst")
PYGLETPLOT = "sympy/plotting/pygletplot"


def test_doctest_config_paths_exist():
    # A path that no longer exists silently matches no file, so the document
    # it refers to is either not tested at all or is tested unexpectedly.
    assert validate_doctest_config(get_sympy_dir()) == []


def test_native_path():
    expected = os.path.normcase(os.path.join("doc", "src", "index.rst"))
    assert native_path("doc/src/index.rst") == expected
    assert native_path(r"doc\src\index.rst") == expected


def test_regular_doctest_blacklist_follows_dependencies():
    none_installed = dict.fromkeys(DEPENDENCIES, False)
    all_installed = dict.fromkeys(DEPENDENCIES, True)

    blacklist = regular_doctest_blacklist(available=none_installed)
    assert blacklist.count(BIOMECHANICS_DOC) == 1

    blacklist = regular_doctest_blacklist(available=all_installed)
    assert BIOMECHANICS_DOC not in blacklist

    # A file without dependencies is skipped either way.
    assert "sympy/this.py" in blacklist


def test_regular_doctest_blacklist_skips_pygletplot_on_ci():
    available = dict.fromkeys(DEPENDENCIES, True)
    ci = os.getenv("CI", None)
    os.environ["CI"] = "true"
    try:
        assert PYGLETPLOT in regular_doctest_blacklist(available=available)
    finally:
        if ci is None:
            del os.environ["CI"]
        else:
            os.environ["CI"] = ci


def test_optional_doctest_files_follows_dependencies():
    none_installed = dict.fromkeys(DEPENDENCIES, False)
    all_installed = dict.fromkeys(DEPENDENCIES, True)

    files = optional_doctest_files(available=all_installed)
    biomechanics = ("doc/src/explanation/modules/physics/biomechanics/"
                    "biomechanics.rst")
    assert biomechanics in files
    # Only documentation is listed here; the modules have their own runner.
    assert not any(f.startswith("sympy/") for f in files)

    assert optional_doctest_files(available=none_installed) == []


def test_skipped_and_optional_files_do_not_overlap():
    assert not set(SKIPPED_FILES) & set(OPTIONAL_FILES)
