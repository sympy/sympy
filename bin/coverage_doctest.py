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


def print_header(name, underline=None, overline=None):
  print
  print name
  if underline: print underline*len(name)

def print_coverage(module_path, c, c_md, c_mdt, f, f_md, f_mdt, score, total_doctests, total_members, verbose=False):


    if verbose:
        print '\n'+'-'*70

    print "%s: %s%% (%s of %s)" % (module_path, score, total_doctests, total_members)

    if verbose:
        print '-'*70


    if verbose:
        print_header('CLASSES', '*')
        if not c:
            print_header('No classes found!')

        else:
            if c_md:
                print_header('Missing docstrings','-')
                for md in c_md:
                   print md
            if c_mdt:
                print_header('Missing doctests','-')
                for md in c_mdt:
                    print md

        print_header('FUNCTIONS','*')
        if not f:
            print_header('No functions found!')
        else:
            if f_md:
                print_header('Missing docstrings', '-')
                for md in f_md:
                   print md
            if f_mdt:
                print_header('Missing doctests', '-')
                for md in f_mdt:
                    print md

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
    m_members = dir(m)


    c_skipped = []
    c_md = []
    c_mdt = []
    c_has_doctest = []
    c_indirect_doctest = []
    classes = 0
    c_doctests = 0

    f_skipped = []
    f_md = []
    f_mdt = []
    f_has_doctest = []
    f_indirect_doctest = []
    functions = 0
    f_doctests = 0

    for member in m_members:
        # Identify if the member (class/def) a part of this module
        obj = getattr(m, member)
        obj_mod = inspect.getmodule(obj)

        # Function not a part of this module
        if not obj_mod or not obj_mod.__name__ == module_path:
          continue

        # If it's a function, simply process it
        if inspect.isfunction(obj) or inspect.ismethod(obj):
            # Various scenarios
            if member.startswith('_'): f_skipped.append(member)
            else:
                if not obj.__doc__: f_md.append(member)
                elif not '>>>' in obj.__doc__: f_mdt.append(member)
                else: f_doctests = f_doctests + 1
                functions = functions + 1

        # If it's a class, look at it's methods too
        elif inspect.isclass(obj):
            classes = classes + 1
            # Process the class first
            if not obj.__doc__: c_md.append(member)
            elif not '>>>' in obj.__doc__: c_mdt.append(member)
            else: c_doctests = c_doctests + 1

            # Iterate through it's members
            for class_m in dir(obj):

                # Check if the method is a part of this module
                class_m_mod = None
                class_m_obj = None

                # Gutsy hack; need to expand reasons
                try:
                    class_m_obj = getattr(obj, class_m)
                    class_m_mod = inspect.getmodule(class_m_obj)
                except:
                    continue

                if not class_m_mod or not class_m_obj or \
                    not class_m_mod.__name__ == module_path:
                    continue

                # Check function for various categories
                full_name = member + '.' + class_m
                if class_m.startswith('_'): f_skipped.append(full_name)
                else:
                    if not class_m_obj.__doc__: f_md.append(full_name)
                    elif not '>>>' in class_m_obj.__doc__: f_mdt.append(full_name)
                    else: f_doctests = f_doctests + 1
                    functions = functions + 1

    total_doctests = c_doctests + f_doctests
    total_members = classes + functions
    if total_members: score = 100 * float(total_doctests) / (total_members)
    else: score = 0
    score = int(score)

    print_coverage(module_path, classes, c_md, c_mdt, functions, f_md, f_mdt, score, total_doctests, total_members, verbose)


    return total_doctests, total_members

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

    # Remove the file extension
    file, ign = os.path.splitext(file)

    # Replace separators by . for module path
    file_module = ""
    h, t = os.path.split(file)
    while h or t:
        if t: file_module = t + '.' + file_module
        h, t = os.path.split(h)

    return coverage(file_module[:-1], verbose)

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
            file = os.path.normpath(file)
            print 'DOCTEST for %s' % (file)
            print '='*70
            print
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
