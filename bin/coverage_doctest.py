#!/usr/bin/env python

"""
Program to test that all methods/functions have at least one example
doctest.


Usage:

bin/coverage_doctest.py sympy/core

or

bin/coverage_doctest.py sympy/core/basic.py

"""

from __future__ import with_statement

import os
import re
import sys
import string
import inspect
from optparse import OptionParser

def _is_skip(module_path, name, obj):
    """ Given a module object and the type object, the function
    determines if the objected should be ignored on following
    remarks:
      * It is a builtin
      * Starts with an underscore (is private)
    """

    skip = False
    filename = None

    # Removes the built-in types, a bit hacky
    try:
        filename = inspect.getfile(obj)
    except TypeError as (strerror):
        skip = True

    # If import imported something other than our module
    if inspect.getmodule(obj).__name__ != module_path:
        skip = True

    if skip or name.startswith('_'):
            skip = True

    return skip

def print_coverage(c, c_md, c_mdt, f, f_md, f_mdt):

    # Print routines (can be deported?)

    print
    print 'CLASSES'
    print '*'*10

    if not c:
        print
        print 'No classes found!\n'
    else:
        if c_md:
            print
            print 'Missing docstrings'
            print '-'*15
            for md in c_md:
               print md[0]
        if c_mdt:
            print
            print 'Missing doctests'
            print '-'*15
            for md in c_mdt:
                print md[0]

    print
    print 'DEFINITIONS'
    print '*'*10
    if not f:
        print
        print 'No functions found!\n'
    else:
        if f_md:
            print
            print 'Missing docstrings'
            print '-'*15
            for md in f_md:
               print md[0]
        if f_mdt:
            print
            print 'Missing doctests'
            print '-'*15
            for md in f_mdt:
                print md[0]

def coverage(module_path, verbose=False):

    """ Given a module path, builds an index of all classes and functions
    contained. It then goes through each of the classes/functions to get
    the docstring and doctest coverage of the module. """

    # Import the package and find membmers
    m = None
    try:
        __import__(module_path)
        m = sys.modules[module_path]
    except:
        # Most likely cause, absence of __init__
        print module_path + ' could not be loaded!'
        return 0, 0

    # Get the list of members (currently everything possible)
    m_members = inspect.getmembers(m)

    # Create a list of class and definitions
    m_classes = filter(lambda x: inspect.isclass(x[1]), m_members)
    m_functions = filter(lambda x: inspect.isfunction(x[1]), m_members)


    print
    print '='*70
    print module_path
    print '='*70

    # Iterate over classes
    c_skipped = []
    c_md = []
    c_mdt = []
    c_has_doctest = []
    c_indirect_doctest = []
    classes = 0

    for c in m_classes:
        class_name, class_obj = c[0], c[1]
        # Check if the class needs to be skipped
        skip = _is_skip(module_path, class_name, class_obj)
        if skip:  c_skipped.append(c)
        else:   classes = classes + 1
        # Check if the class has docstrings
        if not skip and not class_obj.__doc__:
             c_md.append(c)
        # Check for a docstring
        if not skip and class_obj.__doc__ and \
            not '>>>' in class_obj.__doc__:
                c_mdt.append(c)

    # Iterate over functions
    f_skipped = []
    f_md = []
    f_mdt = []
    f_has_doctest = []
    f_indirect_doctest = []
    functions = 0


    for f in m_functions:
        f_name, f_obj = f[0], f[1]
        # Check if the class needs to be skipped
        skip = _is_skip(module_path, f_name, f_obj)
        if skip:  f_skipped.append(f)
        else:   functions = functions + 1
        # Check if the class has docstrings
        if not skip and not f_obj.__doc__:
             f_md.append(f)
        # Check for a docstring
        if not skip and f_obj.__doc__ and \
            not '>>>' in f_obj.__doc__:
                f_mdt.append(f)

    print_coverage(classes, c_md, c_mdt, functions, f_md, f_mdt)

    return 0, 0

def go(file, verbose=False, exact=True):

    if os.path.isdir(file):
        doctests, num_functions = 0, 0
        for F in os.listdir(file):
            _doctests, _num_functions = go('%s/%s'%(file,F), verbose, exact=False)
            doctests += _doctests
            num_functions += _num_functions
        return doctests, num_functions
    if not (file.endswith('.py') or file.endswith('.pyx')) or \
        file.endswith('__init__.py') or \
        not exact and ('test_' in file or 'bench_' in file):
            return 0, 0
    if not os.path.exists(file):
        print "File %s does not exist."%file
        sys.exit(1)
    return coverage(string.replace(file,'/','.')[:-3], verbose)

if __name__ == "__main__":

    bintest_dir = os.path.abspath(os.path.dirname(__file__))   # bin/cover...
    sympy_top  = os.path.split(bintest_dir)[0]      # ../
    sympy_dir  = os.path.join(sympy_top, 'sympy')  # ../sympy/
    if os.path.isdir(sympy_dir):
        sys.path.insert(0, sympy_top)

    parser = OptionParser()
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
            default=False)

    options, args = parser.parse_args()

    if len(args) == 0:
        parser.print_help()
    else:
        for file in args:
            doctests, num_functions = go(file, options.verbose)
            if num_functions == 0:
                score = 100
            else:
                score = 100 * float(doctests) / num_functions
                score = int(score)
            print
            print '='*70
            print "TOTAL SCORE for %s: %s%% (%s of %s)" % \
                (file, score, doctests, num_functions)
            print
