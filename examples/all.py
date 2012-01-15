#!/usr/bin/env python

"""all.py

Runs all the examples for testing purposes and reports successes and failures
to stderr.  An example is marked successful if the running thread does not
throw an exception, for threaded examples, such as plotting, one needs to
check the stderr messages as well.

   $ ./all.py [-hqsw]

Options:
    -h     Print this help message and exit.
    -q     Runs examples in quiet mode.  This will suppress example output and
           show simple status messages.
    -s     Hides the summary at the end of testing the examples.
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
    "intermediate.coupled_cluster",
    "intermediate.differential_equations",
    "intermediate.infinite_1d_box",
    "intermediate.partial_differential_eqs",
    "intermediate.trees",
    "intermediate.vandermonde",
    "advanced.curvilinear_coordinates",
    "advanced.fem",
    "advanced.gibbs_phenomenon",
    "advanced.grover_example",
    "advanced.pidigits",
    "advanced.qft",
    "advanced.relativity",
    ]

WINDOWED_EXAMPLES = [
    "beginner.plotting_nice_plot",
    "intermediate.print_gtk",
    "intermediate.mplot2d",
    "intermediate.mplot3d",
    "advanced.autowrap_integrators",
    "advanced.autowrap_ufuncify",
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


def run_examples(windowed=False, quiet=False, summary=True):
    """Run all examples in the list of modules.

    Returns a boolean value indicating whether all the examples were
    successful.
    """
    successes = []
    failures = []
    examples = TERMINAL_EXAMPLES
    if windowed:
        examples += WINDOWED_EXAMPLES

    for example in examples:
        if run_example(example, quiet):
            successes.append(example)
        else:
            failures.append(example)

    if summary:
        show_summary(successes, failures, quiet)

    return len(failures) == 0


def run_example(example, quiet=False):
    """Run a specific example.

    Returns a boolean value indicating whether the example was successful.
    """
    if quiet:
        print example + " " * (72 - len(example)),
    else:
        print "=" * 79
        print "Running: ", example

    try:
        mod = load_example_module(example)
        if quiet:
            suppress_output(mod.main)
            print "[PASS]"
        else:
            mod.main()
        return True
    except:
        if quiet:
            print "[FAIL]"
        traceback.print_exc()
        return False


class DummyFile(object):
    def write(self, x): pass


def suppress_output(fn):
    """Suppresses the output of fn on sys.stdout."""
    save_stdout = sys.stdout
    try:
        sys.stdout = DummyFile()
        fn()
    finally:
        sys.stdout = save_stdout


def show_summary(successes, failures, quiet=False):
    """Shows a summary detailing which examples were successful and which failed."""
    if quiet:
        print "-" * 79
        if failures:
            print "FAILED:"
            for example in failures:
                print "  " + example
        else:
            print "ALL EXAMPLES PASSED"
    else:
        if successes:
            print >> sys.stderr, "SUCCESSFUL: "
            for example in successes:
                print >> sys.stderr, "  -", example
        else:
            print >> sys.stderr, "NO SUCCESSFUL EXAMPLES"

        if failures:
            print >> sys.stderr, "FAILED: "
            for example in failures:
                print >> sys.stderr, "  -", example
        else:
            print >> sys.stderr, "NO FAILED EXAMPLES"


def main(*args, **kws):
    """Main script runner"""
    windowed = False
    summary = True
    quiet = False

    try:
        opts, remainder = getopt.getopt(args, "hqsw")
        for opt_key, opt_val in opts:
            if opt_key == '-w':
                windowed = True
            elif opt_key == "-h":
                print __doc__
                sys.exit(0)
            elif opt_key == "-s":
                summary = False
            elif opt_key == "-q":
                quiet = True
            else:
                raise getopt.GetoptError, "option %s not processed" % opt_key
    except getopt.GetoptError, message:
        print >> sys.stderr, message
        print >> sys.stderr, "Use -h option for usage.\n"
        sys.exit(1)

    return 0 if run_examples(windowed, quiet, summary) else 1


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
