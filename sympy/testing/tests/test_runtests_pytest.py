import pathlib

import pytest

from sympy.testing.runtests_pytest import (
    make_absolute_path,
    sympy_dir,
    update_args_with_paths,
    update_args_with_rootdir,
)


def test_update_args_with_rootdir():
    """`--rootdir` and directory three above this added as arguments."""
    args = update_args_with_rootdir([])
    assert args == ['--rootdir', str(pathlib.Path(__file__).parents[3])]


class TestMakeAbsolutePath:

    @staticmethod
    @pytest.mark.parametrize(
        'partial_path', ['sympy', 'sympy/core', 'sympy/nonexistant_directory'],
    )
    def test_valid_partial_path(partial_path):
        """Paths that start with `sympy` are valid."""
        _ = make_absolute_path(partial_path)

    @staticmethod
    @pytest.mark.parametrize(
        'partial_path', ['not_sympy', 'also/not/sympy'],
    )
    def test_invalid_partial_path_raises_value_error(partial_path):
        """A `ValueError` is raises on paths that don't start with `sympy`."""
        with pytest.raises(ValueError):
            _ = make_absolute_path(partial_path)


class TestUpdateArgsWithPaths:

    @staticmethod
    def test_no_paths():
        """If no paths are passed, only `sympy` and `doc/src` are appended.

        `sympy` and `doc/src` are the `testpaths` stated in `pytest.ini`. They
        need to be manually added as if any path-related arguments are passed
        to `pytest.main` then the settings in `pytest.ini` may be ignored.

        """
        paths = []
        args = update_args_with_paths(paths=paths, args=[])
        expected = [
            str(pathlib.Path(sympy_dir(), 'sympy')),
            str(pathlib.Path(sympy_dir(), 'doc/src')),
        ]
        assert args == expected

    @staticmethod
    @pytest.mark.parametrize(
        'path',
        ['sympy/core/tests/test_basic.py', '_basic']
    )
    def test_run_one_file(path):
        """Single files/paths, full or partial, are matched correctly."""
        args = update_args_with_paths(paths=[path], args=[])
        expected = [
            str(pathlib.Path(sympy_dir(), 'sympy/core/tests/test_basic.py')),
        ]
        assert args == expected

    @staticmethod
    def test_run_partial_path_from_root():
        """Partial paths from the root directly are matched correctly."""
        args = update_args_with_paths(paths=['sympy/functions'], args=[])
        expected = [str(pathlib.Path(sympy_dir(), 'sympy/functions'))]
        assert args == expected

    @staticmethod
    def test_run_multiple_paths_from_root():
        """Multiple paths, partial or full, are matched correctly."""
        paths = ['sympy/core/tests/test_basic.py', 'sympy/functions']
        args = update_args_with_paths(paths=paths, args=[])
        expected = [
            str(pathlib.Path(sympy_dir(), 'sympy/core/tests/test_basic.py')),
            str(pathlib.Path(sympy_dir(), 'sympy/functions')),
        ]
        assert args == expected

    @staticmethod
    @pytest.mark.parametrize(
        'paths, expected_paths',
        [
            (
                ['/core', '/util'],
                [
                    'sympy/core',
                    'sympy/logic/utilities',
                    'sympy/utilities',
                    'doc/src/modules/utilities',
                    'doc/src/reference/public/utilities',
                ]
            ),
        ]
    )
    @pytest.mark.skip(reason='no way of guaranteeing result when run in CI')
    def test_run_multiple_paths_from_non_root(paths, expected_paths):
        """Multiple partial paths are matched correctly."""
        args = update_args_with_paths(paths=paths, args=[])
        expected_paths = [
            str(pathlib.Path(sympy_dir(), path)) for path in expected_paths
        ]
        assert args == expected_paths
