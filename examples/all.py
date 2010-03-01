#!/usr/bin/env python
"""all.py

Runs all the examples for testing purposes and reports success and failure
to stderr.  An example is marked successful if the running thread does not
throw an exception, for threaded examples, such as plotting, one needs to
check the stderr messages as well.

   $ ./all.py [-hw]

Options:
    -h     print this help message and exit
    -w     Also run examples requiring windowed environment.

Example Usage:
   When no examples fail:
     $ ./all.py > out
     SUCCESSFUL:
       - beginner.basic
       [...]
     NO FAILED EXAMPLES
     $

   When examples fail:
     $ ./all.py -w > out
     Traceback (most recent call last):
       File "./all.py", line 111, in run_examples
     [...]
     SUCCESSFUL:
       - beginner.basic
       [...]
     FAILED:
       - intermediate.mplot2D
       [...]
     $

   Obviously, we want to achieve the first result.
"""

import imp
import os
import sys
import traceback
import getopt

TERMINAL_EXAMPLES = [
    "beginner.basic",
    "beginner.differentiation",
    "beginner.expansion",
    "beginner.functions",
    "beginner.limits_examples",
    "beginner.precision",
    "beginner.print_pretty",
    "beginner.series",
    "beginner.substitution",
    "beginner.expansion",
    "intermediate.coupled_cluster",
    "intermediate.differential_equations",
    "intermediate.partial_differential_eqs",
    "intermediate.trees",
    "intermediate.vandermonde",
    "advanced.fem",
    "advanced.gibbs_phenomenon",
    "advanced.pidigits",
    "advanced.qft",
    "advanced.relativity",
    "advanced.curvilinear_coordinates",
    ]

WINDOWED_EXAMPLES = [
    "beginner.plotting_nice_plot",
    "intermediate.print_gtk",
    "intermediate.mplot2d",
    "intermediate.mplot3d",
    "advanced.plotting",
    ]

EXAMPLE_DIR = os.path.dirname(__file__)

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
    module_path = os.path.join(EXAMPLE_DIR, *name.split('.')[:-1])

    fp, pathname, description = imp.find_module(module_name, [module_path])

    this_file = os.path.abspath(__file__)
    sympy_dir = os.path.join(os.path.dirname(this_file), "..")
    sympy_dir = os.path.normpath(sympy_dir)
    sys.path.insert(0, sympy_dir)

    try:
        return imp.load_module(module_name, fp, pathname, description)
    finally:
        # Since we may exit via an exception, close fp explicitly.
        if fp:
            fp.close()


def load_example_module(example):
    """Loads modules based upon the given package name"""
    mod = __import__(example)
    return mod


def run_examples(windowed=False):
    """Run example in list of modules"""
    success = []
    fail = []
    examples = TERMINAL_EXAMPLES
    if windowed:
        examples += WINDOWED_EXAMPLES
    for example in examples:
        print "="*79
        print "Running: ", example
        try:
            mod = load_example_module(example)
            mod.main()
            success.append(example)
        except:
            traceback.print_exc()
            fail.append(example)
    if success:
        print >> sys.stderr, "SUCCESSFUL: "
        for example in success:
            print >> sys.stderr, "  -", example
    else:
        print >> sys.stderr, "NO SUCCESSFUL EXAMPLES"
    if fail:
        print >> sys.stderr, "FAILED: "
        for example in fail:
            print >> sys.stderr, "  -", example
    else:
        print >> sys.stderr, "NO FAILED EXAMPLES"


def main (*args, **kws):
    """Main script runner"""

    use_windowed = False
    try:
        opts, remainder = getopt.getopt(args, "hw")
        for opt_key, opt_val in opts:
            if opt_key == '-w':
                use_windowed = True
            elif opt_key == "-h":
                print __doc__
                sys.exit(0)
            else:
                raise getopt.GetoptError, "option %s not processed" % opt_key
    except getopt.GetoptError, message:
        print >> sys.stderr, message
        print >> sys.stderr, "Use -h option for usage.\n"
        sys.exit(1)

    run_examples(use_windowed)


if __name__ == "__main__":
    main(*sys.argv[1:])
