#!/usr/bin/env python
"""Runner script

Runs all the known working examples.

 Usage:

 When all examples run:
   $ ./all > out
   $

 When some examples fail:
   $ ./all > out
   Traceback (most recent call last):
     File "./limits_examples.py", line 17, in ?
   [...]
   $

 Obviously, we want to achieve the first result.
"""

import imp
import os
import sys
import traceback

WORKING_EXAMPLES = [
    "beginner.basic",
    "beginner.differentiation",
    "beginner.expansion",
    "beginner.functions",
    "beginner.limits_examples",
    #"beginner.plotting_nice_plot",
    "beginner.precision",
    "beginner.print_pretty",
    "beginner.series",
    "beginner.substitution",
    "beginner.expansion",
    "intermediate.differential_equations",
    #"intermediate.mplot2d",
    #"intermediate.mplot3d",
    #"intermediate.print_gtk",
    "intermediate.trees",
    "intermediate.vandermonde",
    "advanced.fem",
    "advanced.gibbs_phenomenon",
    "advanced.pidigits",
    #"advanced.plotting",
    "advanced.qft",
    "advanced.relativity",
    ]

example_dir = os.path.dirname(__file__)
example_modules = []

def __import__(name, globals=None, locals=None, fromlist=None):
    """An alternative to the import function so that we can import
    modules defined as strings.

    This code was taken from: http://docs.python.org/lib/examples-imp.html
    """
    # Fast path: see if the module has already been imported.
    try:
        return sys.modules[name]
    except KeyError:
        pass

    # If any of the following calls raises an exception,
    # there's a problem we can't handle -- let the caller handle it.
    module_name = name.split('.')[-1]
    module_path = os.path.join(example_dir, *name.split('.')[:-1])

    fp, pathname, description = imp.find_module(module_name, [module_path])

    try:
        return imp.load_module(module_name, fp, pathname, description)
    finally:
        # Since we may exit via an exception, close fp explicitly.
        if fp:
            fp.close()

def load_example_modules ():
    """Loads modules based upon the given package name"""
    global example_modules

    for entry in WORKING_EXAMPLES:
        mod = __import__(entry)
        example_modules.append(mod)

def setup_path ():
    """Put example directories in the path to load easily"""
    sys.path.insert(0,example_dir)

def run_examples():
    success = []
    fail = []
    for mod in example_modules:
        print "="*79
        print "Running: ", mod.__name__
        try:
            mod.main()
            success.append(mod.__name__)
        except:
            traceback.print_exc()
            fail.append(mod.__name__)
    print "SUCCESS: ", success
    print "FAIL: ", fail

def main (*args, **kws):
    setup_path()
    load_example_modules()
    run_examples()

if __name__ == "__main__":
    main(*sys.argv[1:])
