from __future__ import print_function

import os
import shutil
import sys

# Add $INSTALLDIR/bin as a module path so we can import its files as modules.
my_dir = os.path.dirname(os.path.abspath(__file__))
bin_dir = os.path.join(my_dir, "..")
bin_dir = os.path.join(bin_dir, "..")
bin_dir = os.path.join(bin_dir, "..")
bin_dir = os.path.join(bin_dir, "bin")
bin_dir = os.path.normpath(bin_dir)
sys.path.insert(0, bin_dir)

# Import the main diagnose_imports program as a module, so we can test its
# individual functions.
import diagnose_imports

# Undo the path manipulation so that further test code does not inadvertently
# import something from $INSTALLDIR/bin
del sys.path[0]


### Tests to check that command-line options find their way to the subprocess


def pass_options(main_option_text):
    """Pass main_text through main process option parser, subprocess option
    generation, and subprocess option parser.
    Return the result of each stage so the caller can inspect it."""
    main_options = \
        diagnose_imports.parse_main_process_options(main_option_text, False)
    sub_text = diagnose_imports.make_subprocess_options(main_options).as_list()
    sub_options = diagnose_imports.SubprocessOptions(sub_text)
    return main_options, sub_text, sub_options


def test_imports_option():
    main_options, sub_text, sub_options = pass_options(['asdf', 'qwert'])
    assert main_options.imports == ['asdf', 'qwert']
    assert sub_text == ['Siasdf', 'Siqwert']
    assert sub_options.imports == ['asdf', 'qwert']


def test_exclude_option():
    def do_test(main_text):
        main_options, sub_text, sub_options = pass_options(main_text)
        assert main_options.excludes == ['qwert']
        assert sub_text == ['S-qwert']
        assert sub_options.excludes == ['qwert']
    do_test(['-X', 'qwert'])
    do_test(['--exclude', 'qwert'])


def test_include_option():
    def do_test(main_text):
        main_options, sub_text, sub_options = pass_options(main_text)
        assert main_options.includes == ['qwert']
        assert sub_text == ['S+qwert']
        assert sub_options.includes == ['qwert']
    do_test(['-I', 'qwert'])
    do_test(['--include', 'qwert'])


def test_filter_mode_option():
    def do_test(main_text, expected_sub_text, expected_value):
        main_options, sub_text, sub_options = pass_options(main_text)
        assert main_options.filter_mode == expected_value
        assert sub_text == expected_sub_text
        assert sub_options.filter_mode == expected_value
    do_test([], [], 'conservative')
    do_test(['--filter-mode=conservative'], [], 'conservative')
    do_test(['--filter-mode=aggressive'], ['Sfaggressive'], 'aggressive')
    do_test(['--filter-mode=local'], ['Sflocal'], 'local')


def test_duplicate_option():
    def do_test(main_text, expected_sub_text, expected_value):
        main_options, sub_text, sub_options = pass_options(main_text)
        assert main_options.duplicate == expected_value
        assert sub_text == expected_sub_text
        assert sub_options.duplicate == expected_value
    do_test(['--duplicate'], ['Rd'], True)
    do_test([], [], False)


def test_redefinition_option():
    def do_test(main_text, expected_sub_text, expected_value):
        main_options, sub_text, sub_options = pass_options(main_text)
        assert main_options.redefinition == expected_value
        assert sub_text == expected_sub_text
        assert sub_options.redefinition == expected_value
    do_test(['--redefinition'], ['Rr'], True)
    do_test([], [], False)


def test_indirect_option():
    def do_test(main_text, expected_sub_text, expected_value):
        main_options, sub_text, sub_options = pass_options(main_text)
        assert main_options.indirect == expected_value
        assert sub_text == expected_sub_text
        assert sub_options.indirect == expected_value
    do_test(['--indirect'], ['Ri'], True)
    do_test([], [], False)


def test_origins_option():
    def do_test(main_text, expected_sub_text, expected_value):
        main_options, sub_text, sub_options = pass_options(main_text)
        assert main_options.origins == expected_value
        assert sub_text == expected_sub_text
        assert sub_options.origins == expected_value
    do_test(['--origins'], ['Ro'], True)
    do_test([], [], False)


def test_trace_option():
    def do_test(main_text, expected_sub_text, expected_value):
        main_options, sub_text, sub_options = pass_options(main_text)
        assert main_options.trace == expected_value
        assert sub_text == expected_sub_text
        assert sub_options.trace == expected_value
    do_test(['--trace'], ['Rt'], True)
    do_test([], [], False)


def test_subcommand_option():
    def do_test(main_text, expected_sub_text, expected_value):
        main_options, sub_text, sub_options = pass_options(main_text)
        assert main_options.subcommand == expected_value
        assert sub_text == expected_sub_text
        assert sub_options.subcommand == expected_value
    do_test(['--subcommand'], ['Rs'], True)
    do_test([], [], False)


def test_sort_option():
    def do_test(main_text, expected_sub_text, expected_value):
        main_options, sub_text, sub_options = pass_options(main_text)
        assert main_options.sort == expected_value
        assert sub_text == expected_sub_text
        assert sub_options.sort == expected_value
    do_test([], [], 'no')
    do_test(['--sort=no'], [], 'no')
    do_test(['--sort=importer'], ['Osimporter'], 'importer')
    do_test(['--sort=imported'], ['Osimported'], 'imported')


def test_indent_option():
    def do_test(main_text, expected_sub_text, expected_value):
        main_options, sub_text, sub_options = pass_options(main_text)
        assert main_options.indent == expected_value
        assert sub_text == expected_sub_text
        assert sub_options.indent == expected_value
    do_test([], [], 2)
    do_test(['--indent=4'], ['Oi4'], 4)


def test_report_option():
    def do_test(main_text, expected_sub_text, expected_value):
        main_options, sub_text, sub_options = pass_options(main_text)
        assert main_options.report == expected_value
        assert sub_text == expected_sub_text
        assert sub_options.report == expected_value
    do_test([], [], 'global')
    do_test(['--report=global'], [], 'global')
    do_test(['--report=local'], ['Orlocal'], 'local')
    do_test(['--report=calls'], ['Orcalls'], 'calls')


def test_dir_option():
    def do_test(main_text, expected_sub_text, expected_value):
        main_options, sub_text, sub_options = pass_options(main_text)
        assert main_options.dirs == expected_value
        assert sub_text == expected_sub_text
        assert sub_options.dirs == expected_value
    do_test(
        [],
        [],
        [])
    do_test(
        ['--dir=foo', 'bar'],
        ['Odfoo=bar'],
        [('foo', 'bar')])
    do_test(
        ['--dir=foo', 'bar', '--dir=fum', 'baz'],
        ['Odfoo=bar', 'Odfum=baz'],
        [('foo', 'bar'), ('fum', 'baz')])


def test_headers_option():
    def do_test(main_text, expected_sub_text, expected_main, expected_sub):
        main_options, sub_text, sub_options = pass_options(main_text)
        assert main_options.headers == expected_main
        assert sub_text == expected_sub_text
        assert sub_options.headers == expected_sub
    do_test(
        ['--duplicate'],
        ['Rd'],
        None, False)
    do_test(
        ['--duplicate', '--trace'],
        ['Rd', 'Rt', 'Oh'],
        None, True)
    do_test(
        ['--duplicate', '--headers'],
        ['Rd', 'Oh'],
        True, True)
    do_test(
        ['--duplicate', '--trace', '--headers'],
        ['Rd', 'Rt', 'Oh'],
        True, True)
    do_test(
        ['--duplicate', '--no-headers'],
        ['Rd'],
        False, False)
    do_test(
        ['--duplicate', '--trace', '--no-headers'],
        ['Rd', 'Rt'],
        False, False)


### Tests to check that imports are recorded as expected

# We cannot import sympy itself for testing purposes: First, an unknown number
# of its modules and packages might already have been imported for previous
# tests and for the test harness itself.
# So we need to set up a directory of test data that we can import.
# As an added bonus, that directory has not been imported (it's inside a freshly
# created directory that no other code in SymPy or the tests has any business
# with), so we don't have to worry about polluting anybody's namespace, or about
# getting our own imports polluted by anybody else's activities.
# Also, we can import directly into the running process, so we can inspect the
# log directly instead of having to run a separate process and interpret its
# output.

test_data_path = os.path.normpath(os.path.join(my_dir, 'import_test_data'))
shutil.rmtree(test_data_path, ignore_errors=True)
os.mkdir(test_data_path)

test_data = """
### __init__.py
import deepimport1
import import_through_function

def f():
    deepimport1.f()
    import_through_function.f()

### deepimport1.py
import deepimport2

def f():
    deepimport2.f()

### deepimport2.py
def f():
    pass

### import_through_function.py

def f():
    pass

def importing_function():
    import imported_through_function
    imported_through_function.f()

importing_function()

### imported_through_function.py

def f():
    pass
"""
for filedata in test_data.split('\n### '):
    filedata = filedata.strip()
    if filedata:
        filename, contents = filedata.split('\n', 1)
        filename = os.path.normpath(os.path.join(test_data_path, filename))
        with open(filename, mode='w', buffering=0) as fileobj:
            fileobj.write(contents)
            fileobj.write('\n')

# WARNING: The Eclipse plugin for Python (PyDev) manipulates the call stack to
# insert its own hooks and wrappers, and also imports modules for its own
# purposes.
# This will inject PyDev's own modules into the import log while debugging
# run_import.
diagnose_imports.run_import(['sympy.utilities.tests.import_test_data'])

def test_log_exists():
    log = diagnose_imports.import_log
    print (log)
    assert log is not None
