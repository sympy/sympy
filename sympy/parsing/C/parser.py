from __future__ import unicode_literals, print_function
import argparse
import os
import re
import sys

from clang.cindex import (
    CursorKind,
    Index,
    StorageClass,
    TranslationUnit,
    TypeKind,
    UnaryOperator
)


class CodeWriter(object):
    def __init__(self, out, preamble=None):
        self.out = out
        self.line_cleared = True
        self.blank_lines = 2
        self.depth = 0
        self.empty = True

        if preamble:
            self.out.write(preamble)

    def write(self, content):
        if not self.empty:
            for i in range(0, self.blank_lines):
                self.out.write('\n')
        self.blank_lines = 0
        if content:
            if self.line_cleared:
                self.out.write('    ' * self.depth)
            self.out.write(content)
            self.empty = False
            self.line_cleared = False

    def clear_line(self):
        if not self.line_cleared:
            self.out.write('\n')
            self.line_cleared = True
            self.blank_lines = 0

    def clear_minor_block(self):
        if not self.line_cleared:
            self.out.write('\n')
            self.line_cleared = True
        while self.blank_lines < 1:
            self.blank_lines += 1

    def clear_major_block(self):
        if not self.line_cleared:
            self.out.write('\n')
            self.line_cleared = True
        while self.blank_lines < max(1, 2 - self.depth):
            self.blank_lines += 1

    def start_block(self):
        self.empty = True
        self.depth += 1
        self.blank_lines = 2

    def end_block(self):
        self.depth -= 1

class BaseParser(object):
    def __init__(self):
        self.index = Index.create()

    def diagnostics(self, out):
        for diag in self.tu.diagnostics:
            print('%s %s (line %s, col %s) %s' % (
                    {
                        4: 'FATAL',
                        3: 'ERROR',
                        2: 'WARNING',
                        1: 'NOTE',
                        0: 'IGNORED',
                    }[diag.severity],
                    diag.location.file,
                    diag.location.line,
                    diag.location.column,
                    diag.spelling
                ), file=out)

class CodeConverter(BaseParser):
    def __init__(self, name, verbosity=0):
        super(CodeConverter, self).__init__()
        # Tools for debugging.
        self.verbosity = verbosity
        self._depth = 0

        self.root_module = Module(name)
        self.filenames = set()
        self.macros = {}
        self.instantiated_macros = {}

        self.ignored_files = set()
        self.last_decl = []

        self.namespace = self.root_module

    def output(self, module, out):
        module_path = module.split('.')

        mod = None
        for i, mod_name in enumerate(module_path):
            if mod is None:
                mod = self.root_module
            else:
                mod = mod.submodules[mod_name]

            if mod_name != mod.name:
                raise Exception("Unknown module '%s'" % '.'.join(module_path[:i+1]))

        if mod:
            mod.output(CodeWriter(out))
        else:
            raise Exception('No module name specified')

    def _output_module(self, mod, out):
        out.write('===== %s.py ==================================================\n' % mod.full_name)
        mod.output(CodeWriter(out))
        for submodule in mod.submodules.values():
            self._output_module(submodule, out)

    def output_all(self, out):
        self._output_module(self.root_module, out)

    def parse(self, filenames, flags):
        abs_filenames = [os.path.abspath(f) for f in filenames]
        self.filenames.update(abs_filenames)

        for filename in abs_filenames:
            if os.path.splitext(filename)[1] != '.h':
                self.tu = self.index.parse(
                    filename,
                    args=flags,
                    options=TranslationUnit.PARSE_DETAILED_PROCESSING_RECORD
                )
                self.handle(self.tu.cursor, self.root_module)

    def parse_text(self, content, flags):
        for f, c in content:
            abs_filename = os.path.abspath(f)
            self.filenames.add(abs_filename)

            self.tu = self.index.parse(
                f,
                args=flags,
                unsaved_files=[(f, c)],
                options=TranslationUnit.PARSE_DETAILED_PROCESSING_RECORD
            )
            self.handle(self.tu.cursor, self.root_module)



def convert_file(src_in, src_out):
    opts = argparse.ArgumentParser(
            description='Convert C++ code to Python.',
        )

    opts.add_argument(
            'filename',
            metavar='file.cpp',
            help='The file(s) to compile.',
            nargs="+"
        )
    opts.add_argument(
            '-o', '--output',
            metavar='module',
            help='The name of the output module to write.',
        )

    # args = opts.parse_args()
    args.filename = src_in
    args.output = src_out
    converter = CodeConverter('output')
    converter.parse(args.filename)

    converter.diagnostics(sys.stderr)

    if args.output:
        with open('%s.py' % args.output, 'w') as out:
            converter.output('%s.py' % args.output, out)
    else:
        print("Not Implemented")


convert_file("incode.cpp","outcode")
f = open("outcode.py", "r")
print(f.read())
