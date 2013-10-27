"""
Checks that running "import sympy" does not trigger a package import from any
SymPy module (module imports are fine).

Importing sympy packages from inside a SymPy package's __init__ is fragile
because the imported package may itself assume that other packages are already
loaded and fully initialized, causing subtle ordering dependencies.
Reporting package imports as errors does not prevent these problems, but it
removes the indirection through the package import and makes them easier to
diagnose.

Implementation note: Forcing Python into actually unloading already-imported
sympy submodules is a tricky and partly undocumented process. To avoid these
issues, the code simply runs an import of the sympy module in a separate,
pristine Python process.
"""

from __future__ import print_function

import subprocess
import sys
import inspect
import __builtin__

def test_module_imports_are_direct():
    process = subprocess.Popen(
        [sys.executable,__file__],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=-1)
    output, _ = process.communicate()
    assert output == '', "There are import problems:\n" + output

if __name__ == '__main__':

    builtin_import = __builtin__.__import__

    class Definition(object):
        "A symbol with its value."
        def __init__(self, name, value):
            self.name = name
            self.value = value
        def __hash__(self):
            return hash(self.name)
        def __eq__(self, other):
            return self.name == other.name and self.value == other.value
        def __ne__(self, other):
            return not (self == other)
        def __repr__(self):
            return 'Definition(' + repr(self.name) + ',' + repr(self.value) + ')'

    symbol_definers = {} # Maps each function/variable to name of module to define it

    def in_module(what, where):
        "Is what the same or a submodule of where?"
        return what == where or what != None and what.startswith(where + '.')

    def relevant(module):
        return in_module(module, "sympy") and not in_module(module, "sympy.mpmath")

    def tracking_import(module, globals=globals(), locals=[], fromlist=None, level=-1):
        caller_frame = inspect.getframeinfo(sys._getframe(1))
        importer_filename = caller_frame.filename
        importer_module = globals['__name__']
        if importer_filename == caller_frame.filename:
            importer_reference = importer_filename + " line " + str(caller_frame.lineno)
        else:
            importer_reference = importer_filename
        result = builtin_import(module, globals, locals, fromlist, level)
        importee_module = result.__name__
        # We're only interested if importer and importee are in SymPy
        if relevant(importer_module) and relevant(importee_module):
            global symbol_definers
            for symbol in result.__dict__.iterkeys():
                definition = Definition(symbol, result.__dict__[symbol])
                # Note that in nested imports, the innermost one is filled in first.
                # When an __init__ re-exports a function, the initial import will
                # already have been recorded as definer and won't be overwritten.
                if not symbol_definers.has_key(definition):
                    symbol_definers[definition] = importee_module
            if hasattr(result, 'path'):
                print("Error: %s contains package import %s" % (importer_reference, module))
            if fromlist != None:
                symbol_list = fromlist
                if '*' in symbol_list:
                    if (importer_filename.endswith('__init__.py')
                        or importer_filename.endswith('__init__.pyc')
                        or importer_filename.endswith('__init__.pyo')):
                        # We do not check starred imports inside __init__
                        # That's the normal "please copy over its imports to my namespace"
                        symbol_list = []
                    else:
                        symbol_list = result.__dict__.iterkeys()
                for symbol in symbol_list:
                    if not symbol in result.__dict__:
                        print("Error: %s imports %s from %s, which does not define it (yet)"
                              % (importer_reference, symbol, importee_module))
                    else:
                        definition = Definition(symbol, result.__dict__[symbol])
                        symbol_definer = symbol_definers[definition]
                        if symbol_definer != importee_module:
                            print("Error: %s imports %s from %s, but should import it from its definer %s"
                                  % (importer_reference, symbol, importee_module, symbol_definer))
        return result

    __builtin__.__import__ = tracking_import

    __import__('sympy')
