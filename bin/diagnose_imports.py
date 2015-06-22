#!/usr/bin/env python
"""Import diagnostics.

Traces and reports the direct and indirect imports that happen due to a given
module import.

See the output of diagnose_import --help for details."""

# Commands for testing:
# cd ~/Projekte/sympy-workspace/sympy-project/sympy
# export PYTHONPATH=.
# bin/diagnose_imports sympy --dir SYMPY /home/jo/Projekte/sympy-workspace/sympy-project/sympy/sympy --dir PYTON /usr/lib/python2.7 --exclude sympy.mpmath.libmp --origins 2>&1 | less -SN

# This works by installing a hook that gets called at every import.
# The hook records a full trace of all import activity, which is then walked by
# the recording pass to extract the requested information.
#
# A module's initialization code (__init__.py, usually) is run only once, during
# the first import; if we import a module before the hook is installed, we'll be
# unable to import that module and see imports during initialization.
# sys and __builtin__/builtins are preimported by Python, so importing them from
# application code, as an exception, is okay since it does not lose anything.
# (We're lucky that these two modules are exactly what we need to install the
# import hook.)

# The overall strategy is to do the import in a subprocess that starts with a
# fresh Pyhon runtime environment, with nothing imported yet.
# The main process does everything before the import, the subprocess does the
# import and anything that comes after it.
# The main process uses command-line parameters to instruct the subprocess.

import sys


UseDefault = object()
"""Parameter for Option to indicate that it should provide the default value
for the parameter. (None is used to indicated "do not use this parameter at
all", which is different."""


class Option(object):
    """Defines an option for the main program as well as its equivalent for the
    subprocess.

    Instance attributes:
        attr_name
            Attribute name to use for storing this option's value in a DTO.
        short_opt
            Short option character to use on the command line.
        long_opt
            Long option string to use on the command line.
        meta_var
            Name by which this option's value goes in help texts.
        help_text
            Help text to print for this option.
    """

    def __init__(self,
            attr_name,
            short_opt=UseDefault,
            long_opt=UseDefault,
            meta_var=UseDefault,
            action=UseDefault,
            type_=UseDefault,
            default=UseDefault,
            choices=UseDefault,
            nargs=UseDefault,
            help_text=None):
        self.attr_name = attr_name
        if short_opt is UseDefault:
            self.short_opt = '-' + attr_name[0:1]
        elif short_opt is None:
            self.short_opt = None
        else:
            self.short_opt = short_opt
        if long_opt is UseDefault:
            self.long_opt = '--' + attr_name
        elif long_opt is None:
            self.long_opt = None
        else:
            self.long_opt = long_opt
        if meta_var is UseDefault:
            self.meta_var = attr_name.upper()
        else:
            self.meta_var = meta_var
        self.action = action
        self.type_ = type_
        self.default = default
        self.choices = choices
        self.nargs = nargs
        assert help_text is not None
        self.help_text = help_text

    def add_to_option_group(self, option_group):
        """Enters self into an optparse option_group for interpreting the main
        program's command line options.

        The caller must configure the optparse object so that it will store the
        option value in an attribute with name self.attr_name).

        (All options of diagnose_imports happen to be in an option group.)"""
        if self.short_opt:
            if self.long_opt:
                args = (self.short_opt, self.long_opt)
            else:
                args = (self.short_opt, )
        else:
            args = (self.long_opt, )
        kwargs = dict()
        def add_kwarg(key, value):
            if value is not UseDefault:
                kwargs[key] = value
        add_kwarg('dest', self.attr_name)
        add_kwarg('metavar', self.meta_var)
        add_kwarg('action', self.action)
        add_kwarg('type', self.type_)
        add_kwarg('default', self.default)
        add_kwarg('choices', self.choices)
        add_kwarg('nargs', self.nargs)
        add_kwarg('help', self.help_text)
        option_group.add_option(*args, **kwargs)

    def append_to_subprocess_options(self, option_dto, subprocess_options):
        """Takes the value of the attribute named self.attr_name in
        option_dto, and appends it to the string list subprocess_options so that
        load_from_subprocess_options can retrieve that option value for its
        option DTO easily (i.e. without the help of optparse)."""
        raise NotImplementedError()

    def load_from_subprocess_options(self, subprocess_options, option_dto):
        """This is the inverse function of append_to_subprocess_options."""
        raise NotImplementedError()


class BooleanOption(Option):

    def __init__(self,
            attr_name,
            short_opt=UseDefault,
            long_opt=UseDefault,
            help_text=None):
        Option.__init__(self,
            attr_name,
            short_opt=short_opt,
            long_opt=long_opt,
            action='store_true',
            default=False,
            help_text=help_text)


class TriBooleanOption(Option):

    def __init__(self,
            attr_name,
            short_opt=UseDefault,
            long_opt=UseDefault,
            help_text=None):
        Option.__init__(self,
            attr_name,
            short_opt=short_opt,
            long_opt=long_opt,
            action='store_true',
            default=None,
            help_text=help_text)


class ListOption(Option):

    def __init__(self,
            attr_name,
            short_opt=UseDefault,
            long_opt=UseDefault,
            meta_var=UseDefault,
            help_text=None):
        Option.__init__(self,
            attr_name,
            short_opt=short_opt,
            long_opt=long_opt,
            meta_var=meta_var,
            action='append',
            default=[],
            help_text=help_text)


class ChoiceOption(Option):

    def __init__(self,
            attr_name,
            short_opt=UseDefault,
            long_opt=UseDefault,
            choices=None,
            help_text=None):
        assert choices
        assert help_text
        Option.__init__(self,
            attr_name=attr_name,
            short_opt=short_opt,
            long_opt=long_opt,
            meta_var='|'.join(choices),
            default=choices[0],
            choices=choices,
            help_text=help_text)


class IntOption(Option):

    def __init__(self,
            attr_name,
            short_opt=UseDefault,
            long_opt=UseDefault,
            meta_var=UseDefault,
            default=None,
            help_text=None):
        assert default is not None
        assert help_text
        Option.__init__(self,
            attr_name=attr_name,
            short_opt=short_opt,
            long_opt=long_opt,
            meta_var=meta_var,
            type_='int',
            help_text=help_text)
        self.default = default


class PairListOption(Option):

    def __init__(self,
            attr_name,
            short_opt=UseDefault,
            long_opt=UseDefault,
            meta_var=None,
            meta_var_2=None,
            help_text=None):
        assert meta_var is not None
        assert meta_var_2 is not None
        assert help_text is not None
        Option.__init__(self,
            attr_name=attr_name,
            short_opt=short_opt,
            long_opt=long_opt,
            meta_var = meta_var + ' ' + meta_var_2,
            action='append',
            default=[],
            nargs=2,
            help_text=help_text)
        self.meta_var_2 = meta_var_2


class ExcludeOption(ListOption):

    def __init__(self):
        ListOption.__init__(self,
            attr_name='excludes',
            short_opt='-X',
            long_opt='--exclude',
            meta_var='EXCLUDED',
            help_text=
                'do not analyze EXCLUDED, '
                'nor any of its submodules except those in --include '
                '(defaults to no excludes)')


class IncludeOption(ListOption):

    def __init__(self):
        ListOption.__init__(self,
            attr_name='includes',
            short_opt='-I',
            long_opt='--include',
            meta_var='INCLUDED',
            help_text=
                'analyze INCLUDED, '
                'and all its submodules not in --exclude '
                '(defaults to no includes)')

class FilterModeOption(ChoiceOption):

    def __init__(self):
        ChoiceOption.__init__(self,
            attr_name='filter_mode',
            short_opt=None,
            long_opt='--filter-mode',
            choices=['conservative', 'aggressive', 'local'],
            help_text=
                'what to do if include/exclude filtering has a filtered module '
                'import an unfiltered one; '
                'conservative: include the filtered module anyway; '
                'agressive: exclude the unfiltered module anyway; '
                'local: stick with the filter settings - the unfiltered module '
                'will be shown as being imported directly from the next '
                'unfiltered module up the import chain '
                '(defaults to %default)')

class DuplicateOption(BooleanOption):

    def __init__(self):
        BooleanOption.__init__(self,
            attr_name='duplicate',
            help_text=
                'in the error report, '
                'include names imported more than once with the same definition '
                '(such imports are useless)')


class RedefinitionOption(BooleanOption):

    def __init__(self):
        BooleanOption.__init__(self,
            attr_name='redefinition',
            help_text=
                'in the error report, '
                'include names imported with two different definitions '
                '(such a practice can become confusing '
                'as the number of importable names increases)')


class IndirectOption(BooleanOption):

    def __init__(self):
        BooleanOption.__init__(self,
            attr_name='indirect',
            help_text=
                'in the error report, '
                'include names imported from a module '
                'other than the one '
                'that the name was originally defined in '
                '(indirect imports make it harder '
                'to find the origin of a name, '
                'hampering the diagnosis and resolution of import cycles)')


class OriginsOption(BooleanOption):

    def __init__(self):
        BooleanOption.__init__(self,
            attr_name='origins',
            help_text=
                'print the origins report, '
                'which gives the the defining module for each name, '
                'sorted by name')


class TraceOption(BooleanOption):

    def __init__(self):
        BooleanOption.__init__(self,
            attr_name='trace',
            help_text=
                'print the trace report, '
                'which is a raw trace of import activities')


class SubcommandOption(BooleanOption):

    def __init__(self):
        BooleanOption.__init__(self,
            attr_name='subcommand',
            help_text=
                'print the subcommand report, which shows the internal options '
                'passed to the subprocess used to test-run the imports; '
                'this is for debugging purposes')


class SortOption(ChoiceOption):

    def __init__(self):
        ChoiceOption.__init__(self,
            attr_name='sort',
            short_opt=None,
            long_opt=UseDefault,
            choices=['no', 'importer', 'imported'],
            help_text=
                'sort order for the error report; '
                'no: do not sort; '
                'importer: sort by importing module\'s name; '
                'imported: sort by imported module\'s name '
                '(defaults to %default)')


class IndentOption(IntOption):

    def __init__(self):
        IntOption.__init__(self,
            attr_name='indent',
            short_opt=None,
            default=2,
            help_text=
                'How many blanks to indent the output for a subelement '
                '(defaults to %default)')


class ReportOption(ChoiceOption):

    def __init__(self):
        ChoiceOption.__init__(self,
            attr_name='report',
            short_opt=None,
            choices=['global', 'local', 'calls'],
            help_text=
                'what lines of code to report; '
                '"global": module level imports; '
                '"local": all imports (module and function level); '
                '"calls": all imports plus all function calls leading to them '
                '(defaults to %default)')


class DirOption(PairListOption):

    def __init__(self):
        PairListOption.__init__(self,
            attr_name='dirs',
            short_opt=None,
            long_opt='--dir',
            meta_var='NAME',
            meta_var_2='DIR',
            help_text=
                'display any directory prefix DIR as $NAME '
                '(can be repeated)')


class HeadersOption(TriBooleanOption):

    def __init__(self):
        TriBooleanOption.__init__(self,
            attr_name='headers',
            short_opt=None,
            help_text=
                'print a header before each report '
                '(default: no headers if one report, '
                'headers if there are multiple reports)')


class SubprocessOptions(object):
    """Defines the command-line options for the subprocess."""
    # The command line is a series of Xyv with:
    #   X = S for scope options, R for report options, O for output options
    #   y = letter to further select the actual option to set
    #   v = value to pass in (not for yes/no options)
    # *All* yes/no options are present on the command line.
    # Lists are transmitted as a series of Options, e.g. Si<module> gives a
    # module to import and analyse, S-<module> one module to exclude from the
    # reports.

    def as_list(self):
        result = []
        ##################################
        # ==> append_to_subprocess_options
        ##################################
        # Scope options
        result.extend(map(lambda module: 'Si%s' % module, self.imports))
        result.extend(map(lambda module: 'S-%s' % module, self.excludes))
        result.extend(map(lambda module: 'S+%s' % module, self.includes))
        if self.filter_mode != 'conservative':
            result.append('Sf%s' % self.filter_mode)
        # Report options
        if self.duplicate: result.append('Rd')
        if self.redefinition: result.append('Rr')
        if self.indirect: result.append('Ri')
        if self.origins: result.append('Ro')
        if self.trace: result.append('Rt')
        if self.subcommand: result.append('Rs')
        # Output options
        if self.sort != 'no': result.append('Os%s' % self.sort)
        if self.indent != 2: result.append('Oi%s' % self.indent)
        if self.report != 'global': result.append('Or%s' % self.report)
        result.extend(map(lambda dir_name: 'Od%s=%s' % dir_name, self.dirs))
        if self.headers: result.append('Oh')
        return result

    def __init__(self, args=[]):
        ##############
        # ==> __init__
        ##############
        self.args = args # just for __repr__
        # Scope options
        self.imports = []
        self.excludes = []
        self.includes = []
        self.filter_mode = 'conservative'
        # Report options
        self.duplicate = False
        self.redefinition = False
        self.indirect = False
        self.origins = False
        self.trace = False
        self.subcommand = False
        # Output options
        self.sort = 'no'
        self.indent = 2
        self.report = 'global'
        self.dirs = []
        self.headers = False
        for arg in args:
            ##############################################################
            # option == 'Xx' ==> subprocess_option_code = 'Xx' in __init__
            # part after : ==> load_from_subprocess_option
            ##############################################################
            option = arg[0:2]
            value = arg[2:]
            # Scope options
            if option == 'Si': self.imports.append(value)
            elif option == 'S-': self.excludes.append(value)
            elif option == 'S+': self.includes.append(value)
            elif option == 'Sf': self.filter_mode = value
            # Report options
            elif option == 'Rd': self.duplicate = True
            elif option == 'Rr': self.redefinition = True
            elif option == 'Ri': self.indirect = True
            elif option == 'Ro': self.origins = True
            elif option == 'Rt': self.trace = True
            elif option == 'Rs': self.subcommand = True
            # Output options
            elif option == 'Os': self.sort = value
            elif option == 'Oi': self.indent = int(value)
            elif option == 'Or': self.report = value
            elif option == 'Od': self.dirs.append(tuple(value.split('=', 1)))
            elif option == 'Oh': self.headers = True
    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.args)

def parse_main_process_options(args, exit_program=True):

    import optparse

    class OptionParser(optparse.OptionParser):

        def __init__(self, *args, **kwargs):
            optparse.OptionParser.__init__(self, *args, **kwargs)
            self.exit_program = False

        def error(self, msg):
            if self.exit_program:
                optparse.OptionParser.error(self, msg)
            else:
                raise AssertionError(msg)

    report_options = []
    option_parser = OptionParser(
        usage=
            "%prog [OPTION...] MODULE...\n"
            "Diagnose imports triggered by importing the given MODULEs.\n"
            "Later OPTIONs silently override earlier ones.\n"
            "\n"
            "There is no option to extend the directory search path; use the environment\n"
            "variable PYTHONPATH for that.\n"
            "\n"
            "You may see multiple imports of the same module, but module initialization\n"
            "happens only during the initial import. That means that any nested imports\n"
            "during initialization happen and are reported only during the initial import.\n"
            "\n"
            "This tool can not work well for \"as\" imports (\"from foo import bar as baz\")\n"
            "because Python does not provide renaming information to the import tracing\n"
            "machinery.\n"
            "The easiest workaround is to avoid renaming (i.e. \"import foo\", then use\n"
            "\"foo.bar\" instead of \"baz\".)\n"
            "\n"
            "If you suspect that this tool is not reporting all imports, run it with\n"
            "--trace and compare the output with that from \"python -v -c import MODULE\"."
            )
        # Well, there *are* workarounds for "as" imports, it's just that every
        # known one is a bit brittle.
        # Approach 1: Retrieve the caller's frame from the call stack, retrieve
        # the source code, parse its import statement. We'd have to deal with
        # missing sources, and we'd want much more elaborate command-line
        # options to let the user specify where to find any sources.
        # Approach 2: Retrieve the namespace, record the visible names, run the
        # import, retrieve the namespace again and compare with the prerecorded
        # namespace contents. Needs to be careful to not capture names from a
        # nested import. Functions that modify themselves would show up as
        # reinterpretations if called from module initialization code other than
        # the module that imports them; there might be more tricky corner cases.
    option_parser.exit_program = exit_program

    option_group = optparse.OptionGroup(
        option_parser,
        'Scope options',
        'Define which imports to analyze.')
    ExcludeOption().add_to_option_group(option_group)
    IncludeOption().add_to_option_group(option_group)
    FilterModeOption().add_to_option_group(option_group)
    option_parser.add_option_group(option_group)

    option_group = optparse.OptionGroup(
        option_parser,
        'Reports options',
        'Define what reports to print.')
    report_options.append('duplicate')
    DuplicateOption().add_to_option_group(option_group)
    report_options.append('redefinition')
    RedefinitionOption().add_to_option_group(option_group)
    report_options.append('indirect')
    IndirectOption().add_to_option_group(option_group)
    report_options.append('origins')
    OriginsOption().add_to_option_group(option_group)
    report_options.append('trace')
    TraceOption().add_to_option_group(option_group)
    report_options.append('subcommand')
    SubcommandOption().add_to_option_group(option_group)
    option_parser.add_option_group(option_group)

    option_group = optparse.OptionGroup(
        option_parser,
        'Output options')
    SortOption().add_to_option_group(option_group)
    IndentOption().add_to_option_group(option_group)
    ReportOption().add_to_option_group(option_group)
    DirOption().add_to_option_group(option_group)
    HeadersOption().add_to_option_group(option_group)
    option_group.add_option(
        '--no-headers',
        dest='headers',
        action='store_false',
        default=None,
        help=
            'print reports without headers '
            '(default: no headers if one report, '
            'headers if there are multiple reports)'
        )
    option_parser.add_option_group(option_group)

    options, imports = option_parser.parse_args(args)
    options.imports = imports
    options.report_option_count = \
        sum(1 for option in report_options if getattr(options, option))
    options.parser = option_parser
    return options


def validate_options(options):
    # Check that we have something to do
    if not options.imports:
        options.parser.error('No module to import given')
    # Check that no module is both --included and --excluded
    if set(options.includes) & set(options.excludes):
        options.parser.error(
            '"%s" both on --include and --exclude'
            % '", "'.join(set(options.include) & set(options.exclude)))
    # Check that at least one report option is present
    if options.report_option_count == 0:
        options.parser.error(
            "At least one report option must be given (try --help)")


def make_subprocess_options(options):
    subprocess_options = SubprocessOptions()
    # Scope options
    subprocess_options.imports = options.imports
    subprocess_options.excludes = options.excludes
    subprocess_options.includes = options.includes
    subprocess_options.filter_mode = options.filter_mode
    # Report options
    subprocess_options.duplicate = options.duplicate
    subprocess_options.redefinition = options.redefinition
    subprocess_options.indirect = options.indirect
    subprocess_options.origins = options.origins
    subprocess_options.trace = options.trace
    subprocess_options.subcommand = options.subcommand
    # Output options
    subprocess_options.sort = options.sort
    subprocess_options.indent = options.indent
    subprocess_options.report = options.report
    subprocess_options.dirs = options.dirs
    if options.headers is None:
        subprocess_options.headers = options.report_option_count > 1
    else:
        subprocess_options.headers = options.headers
    return subprocess_options


def prepare_diagnostics():
    """Main process: Inspect the command line using optparse, then start a
    subprocess using this same program with a command line that's easy to
    inspect without using optparse.

    The goal is to allow the subprocess to do the imports to diagnose without
    having to import other modules first."""

    options = parse_main_process_options(sys.argv[1:])
    validate_options(options)

    # Option parsing and validation will do a hard program abort if anything was
    # wrong with the command line.
    # So we can now assume everything is okay with the options object.
    import subprocess
    process = subprocess.Popen(
        [
            sys.executable, # the Python interpreter
            __file__, # this program
            '--!', # "run as subprocess" marker
        ] + make_subprocess_options(options).as_list(),
        bufsize = -1)
    process.wait()


class ImportLogRecord(object):
    """Import statements and their callers, as a call tree.

    Attributes:
    * file_name: Absolute name of source file (None if from command line)
    * line of import (index if from command line)
    * module: imported module name if import statement, None otherwise.
    * names: The list of imported names if import statement.
        Has all names of the module if the import list contained a star.
        An empty list if the statement was not an import statement.
    * children
    * seq: sequence number (chronologically ascending)

    Two kinds of stack frames may be recorded in an ImportLogRecord:
    * An import statement.
    * A function call that directly or indirectly calls an import statement.

    NOTE: The canonical implementation for this would use namedtuple, but
    that would require an import.
    Also, namedtuples like to be immutable, and we will want to update
    fields."""

    last_seq = 0
    '''Last seq attribute issued for an ImportLogRecord object'''

    def __init__(self,
        file_name, module, line,
        imported_module = None, imported_names = [],
        children = [],
        seq = None
    ):
        self.file_name = file_name
        self.module = module
        self.line = line
        self.imported_module = imported_module
        self.imported_names = []
        self.children = []
        if seq == None:
            self.__class__.last_seq = self.__class__.last_seq + 1
            self.seq = self.__class__.last_seq
        else:
            self.seq = seq

    def __repr__(self):
        return "%s(file_name=%r,module=%r,line=%r,\nimported_module=%r,\nimported_names=%r,\nchildren=%r,\nseq=%r)" % (
            self.__class__.__name__,
            self.file_name, self.module, self.line,
            self.imported_module,
            self.imported_names,
            self.children,
            self.seq)

# The stack frame that started the current import run
import_root_frame = None

# An ImportLogRecord for each module to diagnose.
# This intentionally is the same type as that of ImportLogRecord.children.
import_log = []

# For each name, a list of ImportLogRecords where the name was imported.
origins = {}

def run_import(imports):
    """Run the import, with a logging hook installed.
    Return the log.
    No imports (other than sys and builtins) inside this function, until after
    the import to be diagnosed has finished!"""
    # Set up the import hook and the associated log data structures
    global import_root_frame
    import_root_frame = sys._getframe()
    if sys.version_info[0] > 2:
        import builtins
    else:
        import __builtin__ as builtins
    global _builtin_import
    _builtin_import = builtins.__import__
    builtins.__import__ = _import_wrapper
    global _active_log_records
    _active_log_records = []
    global _active_import_stack_frames
    _active_import_stack_frames = []
    global import_log
    # Run the imports
    try:
        for (i, module) in enumerate(imports, start=1):
            import_log.append(ImportLogRecord(None, None, i))
            __import__(module)
    finally:
        # Restore the import hook
        builtins.__import__ = _builtin_import

def _import_wrapper(module, globals_=globals(), locals_=[], fromlist=None, level=-1):
    """This is a replacement for __import__ that traces imports.
    It recordes the current stack state in import_log and calls the original
    __import__ (as saved in _builtin_import).
    """

    # This function will find stack frames, from current frame outwards:
    # 0) _import_wrapper itself
    # 1) import statement (or __import__ call)
    # 2) zero or more normal (non-importing) function calls
    # 3) Restart with (0)
    # This function can assume that the outermost frame is not of type (0).

    # Stack frames are unhashable, as their f_lineno attributes are modified
    # during the course of Python execution.

    # The import log is a tree that
    # - does not containy any (0) frames
    # - contains an "import statement" entry for each (1) frame
    # - contains a "subroutine call" entry for each (2) frame
    #   these will have None in the imported_module and imported_names
    #   attributes

    # We expect the stack record to
    # - be empty (if we're in the outermost _import_wrapper call), or
    # - end with an "import statement" record (if we're in a nested
    #   _import_wrapper call).

    current_frame = sys._getframe()
    frames = []
    global __name__
    global import_root_frame
    #print('Scanning stack')
    while current_frame != import_root_frame:
        def get_from_dict(d, key):
            not_found = object()
            value=d.get(key, not_found)
            if value is not_found:
                result = '[Not present]'
            elif value is None:
                result = 'None'
            else:
                result= repr(value)
            return result
        #frame_globals = current_frame.f_globals
        #print('package: %s' % get_from_dict(frame_globals,  '__package__'))
        #print('   name: %s' % get_from_dict(frame_globals, '__name__'))
        #print('   file: %s' % get_from_dict(frame_globals, '__file__'))
        #print('   line: %d' % current_frame.f_lineno)
        if __name__ != current_frame.f_globals['__name__']:
            # not same module as this function, hence not type (0)
            frames.append(current_frame)
        current_frame = current_frame.f_back
    global import_log
    current_log = import_log
    for frame in reversed(frames):
        current_log = current_log[-1].children
        frame_file = frame.f_code.co_filename
        frame_module = frame.f_globals['__name__']
        frame_line = frame.f_lineno
        if (not current_log
            or current_log[-1].file_name != frame_file
            or current_log[-1].module != frame_module
            or current_log[-1].line != frame_line
        ):
            current_log.append(ImportLogRecord(
                file_name=frame_file, module=frame_module, line=frame_line))
    log_entry = current_log[-1]
    # Run the import and record the outcome
    result = None
    try:
        result = _builtin_import(module, globals_, locals_, fromlist, level)
        # Import worked, now record the actual module name
        log_entry.imported_module = result.__name__
    except BaseException as e:
        log_entry.imported_module = \
            'importing %s failed with with exception %s' % (module, e)
    finally:
        # Fill in imported_names
        if fromlist != None:
            # fromlist == None means we have "import <module>"
            # instead of "from <module> import <name>, ..."
            # These do not add to the module's namespace.
            if '*' in fromlist:
                # Starred import, equivalent to importing all names that the
                # module defines. These are the keys of result.__dict__.
                log_entry.imported_names.extend(result.__dict__.iterkeys())
            else:
                log_entry.imported_names.extend(fromlist)
            global origins
            for name in log_entry.imported_names:
                origins.setdefault((log_entry.file_name, name), []).append(log_entry)
    # Return the outcome from the original __import__ in _builtin_import
    return result

def printable_filename(file_name, dirs):
    """file_name, with the first matching prefix from dirs replaced with a
    dollar sign and the corresponding name.

    >>> printable_filename(
    ...     '/home/john/sympy/install/sympy.py',
    ...     [['SYMPY': '/home/john/sympy/install']]
    $SYMPY_DIR/sympy.py
    """
    if file_name == None:
        return None
    import os.path
    for (dir_shorthand, dir_name) in dirs:
        if dir_name == os.path.commonprefix([dir_name, file_name]):
            return os.path.join('$' + dir_shorthand, os.path.relpath(file_name, dir_name))
    return file_name

def filtered(log_entry, options):
    '''Should log_entry be filtered according to options?'''
    filtered = None
    module_name = log_entry.imported_module
    if module_name != None:
        # It's an import, not a call
        while filtered != None:
            (_, _, module_name) = module_name.rpartition('.')
            if module_name == module_name:
                break
            if module_name in options.includes:
                filtered = False
            if module_name in options.excludes:
                filtered = True
    return filtered

def nop(*args):
    """The no-operation function. Does nothing except returning True."""
    return True

def iterate_log(preorder, postorder, options):
    """Iterate over the log entries, applying the functions before resp.
    after visiting subentries.

    Specify preorder(log_entry, options) to do something before subentries are
    visited (options is the SubprocessOptions object).
    log is the current ImportLogRecord.
    Let preorder return False to prevent iterate_log from visiting subentries.

    Specify postorder(log_entry, options) to do something after subentries have
    been visited.
    Parameters are as for preorder.
    The return value is ignored."""
    global import_log
    _iterate_log(preorder, postorder, import_log, options)

def _iterate_log(preorder, postorder, list_of_log_entries, options):
    for log_entry in list_of_log_entries:
        if filtered(log_entry, options):
            # Visit the subentries in case they are included
            _iterate_log(preorder, postorder, log_entry.children, options)
        else:
            # Sandwich subnodes between preorder() and postorder() calls
            if preorder(log_entry, options):
                _iterate_log(preorder, postorder, log_entry.children, options)
            postorder(log_entry, options)

def print_import(options, log_entry, depth):
    indent = '%*s' % (options.indent * depth, '')
    if log_entry.module == None:
        location = 'Command-line import request #%d' % log_entry.line
    else:
        location = '%s (%s) line %d' % (
            log_entry.module,
            printable_filename(log_entry.file_name, options.dirs),
            log_entry.line)
    if log_entry.imported_module == None:
        msg = ' (function call)'
    else:
        msg = ': '
        if log_entry.imported_names:
            msg += 'from %s import %s' % (
                log_entry.imported_module,
                ', '.join(log_entry.imported_names))
        else:
            msg += 'import %s' % log_entry.imported_module
    print('%s%s%s' % (indent, location, msg))
    if (log_entry.children):
        print_import_log(options, log_entry.children, depth + 1)

def print_import_log(options, log_entry_list=import_log, depth=0):
    for log_entry in log_entry_list:
        print_import(options, log_entry, depth)

def generate_output(options):
    global origins
    # Useful to check the data structures when debugging
    # print(origins) # >200 MB due to repetitions in output
    print_import_log(options)
    global import_log; print(import_log) # 1-2 MB
    # cd Projekte/sympy-workspace/sympy-project/sympy/
    # export PYTHONPATH=.
    # bin/diagnose_imports sympy --dir SYMPY /home/jo/Projekte/sympy-workspace/sympy-project/sympy/sympy --dir PYTON /usr/lib/python2.7 --origins 2>&1 | less -SN

    return;

    # All other logging options
    # These visit all included nodes, apply whatever checking is requested,
    # collect the output lines in a map where the key establishes the requested
    # sort order, then print the map.
    # TODO
    # As report option:
    # --duplicate
    # --redefinition
    # --indirect
    # --origins
    # As output formatter/mangler:
    # --sort no|importer|imported
    # --report global|local|calls
    from operator import attrgetter
    if options.sort == 'no':
        key_fn = attrgetter('seq')
    elif options.sort == 'importer':
        key_fn = attrgetter('file_name', 'seq')
    elif options.sort == 'imported':
        key_fn = attrgetter('module', 'seq')

    # --origins processing
    if options.origins:
        import collections.Counter

        for key in sorted(origins.keys()):
            for log_entry in origins[key]:
                if not filtered(log_entry, options):
                    (importer, name) = key
                    print('%s.%s from %s (%s line %d)' % (
                        log_entry.module,
                        name,
                        log_entry.imported_module,
                        printable_filename(log_entry.file_name, options.dirs),
                        log_entry.line))

    # --trace processing
    if options.trace:
        options.depth = 0
        def trace_pre(log_entry, options):
            if log_entry.file_name == None:
                output_line = "%*sModule #%d" % (
                    options.indent * options.depth, '',
                    log_entry.line)
            else:
                output_line = "%*s%s line %d" % (
                    options.indent * options.depth, '',
                    printable_filename(log_entry.file_name, options.dirs),
                    log_entry.line)
            if log_entry.imported_module != None:
                output_line += ': '
                if log_entry.imported_names:
                    output_line += 'from %s import %s' % (
                        log_entry.imported_module,
                        ', '.join(log_entry.imported_names))
                else:
                    output_line += 'import %s' % log_entry.imported_module
            print(output_line)
            options.depth = options.depth + 1
            return True
        def trace_post(log_entry, options):
            options.depth = options.depth - 1
        iterate_log(trace_pre, trace_post, options)

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--!':
        options = SubprocessOptions(sys.argv[2:])
        run_import(options.options)
        generate_output(options)
    else:
        prepare_diagnostics()
