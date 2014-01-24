#!/usr/bin/env python

"""
Import diagnostics. Run bin/diagnose_imports.py --help for details.
"""

### imports ####################################################################

from __future__ import print_function

# Doing an import a second time is different from doing it the first time.
# Since we want to diagnose import processing inside SymPy, we cannot import
# any SymPy modules; in fact this test is not a normal SymPy test module,
# it is a standalone program.

# We do some imports of the Python standard modules.
# This does obscure the import diagnostics for these modules.
# This is acceptable since the main purpose of this program is to diagnose
# the import of SymPy modules.

# Can't import sympy.core.compatibility.
# Duplicate of its logic what we need.
import sys
PY3 = sys.version_info[0] > 2 # from sympy.core.compatibility
if PY3:
    import builtins
else:
    import __builtin__ as builtins

import optparse # No argparse in Python 2.6
import collections
import inspect
import os.path

### Option processing ##########################################################

def parse_options():
    option_parser = optparse.OptionParser(
        usage=
            "%prog [options]\n"
            "Run diagnostics on the imports in a software package.\n"
            "Later options override earlier ones.",
        epilog=
            'If Python encounters the same module again during an import, '
            'it will do the import (and it will be recorded by this program), '
            'but Python will not repeat running the module\'s initialization '
            'code and import statements contained in it.'
            '-- '
            'This can not work well for "as" imports '
            '("from foo import bar as baz") '
            'because Python does not provide renaming information '
            'to import hooks. '
            'WORKAROUND: Avoid renaming imports such as the one above; '
            'instead, do "import foo" and use "foo.bar" instead of "baz".')

    option_group = optparse.OptionGroup(
        option_parser,
        'Scope options',
        'Define which imports to analyse:')
    option_group.add_option(
        '--import',
        action='store',
        default='sympy',
        dest='import_', # import is a keyword, can't use options.import
        metavar="MODULE",
        help=
            'What module to import and analyse')
    option_group.add_option(
        '--recursive',
        action='store',
        type='choice',
        choices=['yes', 'no'],
        default='yes',
        metavar="yes|no",
        help=
            'Whether to analyse imports done inside imports')
    option_group.add_option(
        '--exclude',
        action='append',
        default=['sympy.mpmath'],
        help=
            'A module to avoid analyzing; '
            'naming a module excludes submodules as well; '
            'by default, sympy.mpmath is excluded')
    option_group.add_option(
        '--include',
        action='append',
        default=[],
        help=
            'A module to analyze anyway, even if it is listed with --include; '
            'naming a module includes submodules as well')
    option_parser.add_option_group(option_group)

    option_group = optparse.OptionGroup(
        option_parser,
        'Analysis options',
        'Define what kind of analysis to run:')
    option_group.add_option(
        '--duplicate',
        action='store_true',
        default=False,
        help=
            'print names imported more than once with the same definition '
            '(such imports are useless)')
    option_group.add_option(
        '--redefinition',
        action='store_true',
        default=False,
        help=
            'print names imported with two different definitions '
            '(such imports can be misleading, '
            'particularly in SymPy '
            'with its large number of global symbols)')
    option_group.add_option(
        '--indirect',
        action='store_true',
        default=False,
        help=
            'print names imported from a module '
            'other than the one that the name was defined in '
            '(every indirect import makes it harder '
            'to find the origin of a name, or to deal with import cycles)')
    option_group.add_option(
        '--problems',
        action='store_true',
        default=False,
        help=
            'all analyses that are relevant to a SymPy release; '
            'currently, these are --duplicate, --redefinition, and --indirect')
    option_group.add_option(
        '--origins',
        action='store_true',
        default=False,
        help=
            'print the defining module for each name, sorted by name')
    option_group.add_option(
        '--trace',
        action='store_true',
        default=False,
        help=
            'print a raw trace of import activities')
    option_parser.add_option_group(option_group)

    option_group = optparse.OptionGroup(
        option_parser,
        'Output options')
    option_group.add_option(
        '--sort',
        choices=['no', 'importer', 'imported'],
        default='no',
        metavar="no|importer|imported",
        help=
            'Sort output lines; '
            'no: do not sort, '
            'i.e. print information in the order that it became available; '
            'importer: sort by importing module\'s name; '
            'imported: sorty by imported module\'s name')
    option_group.add_option(
        '--indent',
        type='int',
        default=2,
        help=
          'How many blanks to indent the output for a subelement; '
          'default %default')
    option_group.add_option(
        '--report',
        type='choice',
        choices=['global', 'local', 'calls'],
        default='calls',
        metavar="global|local|calls",
        help=
            'What lines of code to report; '
            'global: module level imports; '
            'local: all imports (module and function level); '
            'calls: all imports plus all function calls leading to them')
    option_parser.add_option_group(option_group)

    (options, args) = option_parser.parse_args()
    if args:
        option_parser.error('Unexpected arguments: %s' % ' '.join(args))
    if options.problems:
        options.duplicate = True
        options.redefinition = True
        options.indirect = True
    if (not options.duplicate
        and not options.redefinition
        and not options.indirect
        and not options.problems
        and not options.origins
        and not options.trace):
        option_parser.error("At least one analysis option must be given")
    return options

### Running the import #########################################################

def super_dir(file_name, levels=1):
    result = os.path.abspath(file_name)
    for _ in range(levels):
        result = os.path.join(result, '..')
    result = os.path.normpath(result)
    return result

def sub_dir(dir_name, additional_name):
    result = os.path.abspath(dir_name)
    result = os.path.join(result, additional_name)
    result = os.path.normpath(result)
    return result

sympy_dir = super_dir(__file__, 3)
sys.path.insert(0, sympy_dir)

import_log = None
"""The ImportLogRecord of the top-level import analyzed."""

# Information from a call stack frame that was active during an import statement.
# This will be mostly import statements with the occasional function call.
# The latter get introduced in imports inside functions.
ImportLogRecord = collections.namedtuple("ImportLogRecord", [
    'file', # File of statement
    'line_number', # Line number of statement
    'imported_module', # None if not an import statement
    'name_list', # Names in "from" clause
        # None if there was no such clause (or not an import statement)
        # List of all names in the module's dict if the name list contained '*'
    'children', # List of ImportLogRecords that have self as parent
    ])

import_log = None

_active_log_records = []
"""The currently active ImportLogRecords."""

_active_import_stack_frames = []
"""The call stack frames that point to import statements."""

def _import_wrapper(module, globals=globals(), locals=[], fromlist=None, level=-1):
    """This is a replacement for __import__.

    It records the current state in import_log and calls _builtin_import.
    It is used by recording_import."""

    # Find out if there are any stack frames between the next outer import and
    # the top of _activeLogRecords. If there are, we're inside an import inside
    # one or more function calls and want to record these stack frames with a
    # module and name_list of None.
    global _active_log_records
    global _active_import_stack_frames
    my_frame = sys._getframe()
    importer_frame = sys._getframe(1)
    frames = []
    if len(_active_import_stack_frames) > 0:
        frame = importer_frame.f_back
        while frame != _active_import_stack_frames[-1]:
            if frame.f_code.co_filename != my_frame.f_code.co_filename:
                # Record only if it's not code from this file
                frames.append(frame)
            frame = frame.f_back
        for frame in reversed(frames):
            _active_log_records.append(ImportLogRecord(
                file=frame.f_code.co_filename,
                line_number=frame.f_lineno,
                imported_module=None,
                name_list=None,
                children=[] # Incremental fill-in during nested imports
                ))
            if len(_active_log_records) > 1:
                _active_log_records[-2].children.append(_active_log_records[-1])
    (file, lineno, _, _, _) = inspect.getframeinfo(importer_frame)
    _active_log_records.append(ImportLogRecord(
        file=file,
        line_number=lineno,
        imported_module=module,
        name_list=[], # Filled below
        children=[], # Incremental fill-in during nested imports
        ))
    if len(_active_log_records) > 1:
        _active_log_records[-2].children.append(_active_log_records[-1])
    _active_import_stack_frames.append(importer_frame)
    # Run the import and record the outcome
    result = None
    try:
        result = _builtin_import(module, globals, locals, fromlist, level)
    finally:
        try:
            # Fill in name_list
            if fromlist != None:
                # fromlist == None means we have "import <module>"
                # instead of "from <module> import <name>, ..."
                # These don't add to the module's global namespace, and hence are
                # irrelevant for most (but maybe not all) reports.
                if '*' in fromlist:
                    # Starred import, equivalent to importing all names that the
                    # module defines. These are the keys of result.__dict__.
                    _active_log_records[-1].name_list.extend(result.__dict__.iterkeys())
                else:
                    _active_log_records[-1].name_list.extend(fromlist)
        finally:
            global import_log
            for _ in range(len(frames) + 1):
                import_log = _active_log_records.pop()
            _active_import_stack_frames.pop()
    # Return the outcome from the original __import__ in _builtin_import
    return result

def recording_import(module):
    global import_log
    import_log = None
    global _builtin_import
    _builtin_import = builtins.__import__
    builtins.__import__ = _import_wrapper
    global _active_log_records
    _active_log_records = []
    global _active_import_stack_frames
    _active_import_stack_frames = []
    try:
        __import__(module)
    finally:
        builtins.__import__ = _builtin_import


### Generating output ##########################################################

def printable_filename(file_name, paths):
    """file_name, with the values from paths replaced by the keys.

        >>> printable_filename(
        ...     '/home/john/workspace/sympy/install/sympy.py',
        ...     {'SYMPY': '/home/john/workspace/sympy/install'})
        $SYMPY_DIR/sympy.py
    """
    for (path_name, path) in paths:
        if path == os.path.commonprefix([path, file_name]):
            return os.path.join('$' + path_name, os.path.relpath(file_name, path))
    return file_name

def dump_log(log_record, depth, options):
    if log_record == None:
        return
    str = '%*s%s line %d' % (
          depth * options.indent, '',
          printable_filename(log_record.file, options.paths),
          log_record.line_number,
          )
    if log_record.imported_module != None:
        str += ': '
        if log_record.name_list != None:
            str += 'from %s import %s' % (
                log_record.imported_module,
                ', '.join(log_record.name_list))
        else:
            str += 'import %s' % log_record.imported_module
    print(str)
    if log_record.children != None:
        for child in log_record.children:
            dump_log(child, depth + 1, options)

def dump_origins(options):
    pass # FIXME implement

def report_problems(options):
    pass # FIXME implement

def report(options):
    global import_log
    if options.trace:
        dump_log(import_log, 0, options)
    if options.origins:
        dump_origins(options)
    if options.duplicate or options.redefinition or options.indirect:
        report_problems(options)

### Main program ###############################################################

if __name__ == "__main__":
    options = parse_options()
    sympy_install_dir = super_dir(sympy_dir)
    options.paths = [
        ('PYTHON', super_dir(sys.executable, 2)),
        ('SYMPY_BIN', sub_dir(sympy_install_dir, 'bin')),
        ('SYMPY', sympy_dir),
        ('SYMPY_INSTALL', sympy_install_dir),
        ]
    print(options.paths)
    recording_import(options.import_)
    report(options)

### Currently unused code ######################################################

def in_module(a, b):
    """Is a the same module as or a submodule of b?"""
    return a == b or a != None and b != None and a.startswith(b + '.')

def is_relevant(module):
    """Is module relevant for import checking?

    Only imports between relevant modules will be checked."""
    return in_module(module, 'sympy') and not in_module(module, 'sympy.mpmath')

def is_package_init_file(filename):
    return (filename.endswith('__init__.py')
        or filename.endswith('__init__.pyc')
        or filename.endswith('__init__.pyo'))

sorted_messages = []

def msg(msg, *args):
    global options, sorted_messages
    if options.by_process:
        print(msg % args)
    else:
        sorted_messages.append(msg % args)
