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

def print_coverage(module_path, c, c_md, c_mdt, c_idt, f, f_md, f_mdt, f_idt, score, total_doctests, total_members, verbose=False):

    """ Prints details (depending on verbose) of a module """

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
                   print '\t* '+md
            if c_mdt:
                print_header('Missing doctests','-')
                for md in c_mdt:
                    print '\t* '+md
            if c_idt:
                print_header('Indirect doctests', '-')
                for md in c_idt:
                    print '\t* '+md


        print_header('FUNCTIONS','*')
        if not f:
            print_header('No functions found!')
        else:
            if f_md:
                print_header('Missing docstrings', '-')
                for md in f_md:
                   print '\t* '+md
            if f_mdt:
                print_header('Missing doctests', '-')
                for md in f_mdt:
                    print '\t* '+md
            if f_idt:
                print_header('Indirect doctests', '-')
                for md in f_idt:
                    print '\t* '+md

def _is_indirect(member, doc):

  """ Given string repr of doc and member checks if the member
  contains indirect documentation """

  d = member in doc
  e = 'indirect doctest' in doc
  if not d and not e:
      return True
  else:
      return False

def _get_arg_list(name, fobj):

    """ Given a function object, constructs a list of arguments
    and their defaults. Takes care of varargs and kwargs """

    trunc = 20 # Sometimes argument length can be huge

    argspec = inspect.getargspec(fobj)

    arg_list = []

    if argspec.args:
        for arg in argspec.args: arg_list.append(str(arg))

    arg_list.reverse()

    # Now add the defaults
    if argspec.defaults:
        rev_defaults = list(argspec.defaults).reverse()
        for i in range(len(argspec.defaults)):
            arg_list[i] = str(arg_list[i]) + '=' + str(argspec.defaults[-i])

    # Get the list in right order
    arg_list.reverse()

    # Add var args
    if argspec.varargs:
          arg_list.append(argspec.varargs)
    if argspec.keywords:
          arg_list.append(argspec.keywords)

    # Truncate long arguments
    arg_list = map(lambda x: x[:trunc], arg_list)

    # Construct the parameter string (enclosed in brackets)
    str_param = "%s(%s)" % (name, ', '.join(arg_list))

    return str_param

def process_function(name, c_name, b_obj, mod_path, f_sk, f_md, f_mdt, f_idt, f_has_doctest, sk_list):

    """ Processes a function to get information regarding documentation.
    It is assume that the function calling this subrouting has already
    verified that it is a valid module function """

    if name in sk_list: return False, False

    # We add in the end, as inspect.getsourcelines is slow
    add_md = False
    add_mdt = False
    add_idt = False
    f_doctest = False
    function = False

    if inspect.isclass(b_obj):
        obj = getattr(b_obj, name)
    else:
        obj = b_obj

    # Check function for various categories
    if inspect.isclass(b_obj):
        full_name = _get_arg_list(c_name + '.' + name, obj)
    else:
        full_name = _get_arg_list(name, obj)
    if name.startswith('_'): f_sk.append(full_name)
    else:
        if not obj.__doc__:
            add_md = True
        elif not '>>>' in obj.__doc__:
            add_mdt = True
        else:
            # Indirect doctest
            if _is_indirect(name, obj.__doc__):
                add_idt = True
            f_doctest = True
        function = True

    if add_md or add_mdt or add_idt:

        try:
            line_no = inspect.getsourcelines(obj)[1]
        except IOError:
            # Raised when source does not exist
            # which means the function is not there.
            return False, False

        full_name = "LINE %d: %s" % (line_no, full_name)
        if add_md: f_md.append(full_name)
        elif add_mdt: f_mdt.append(full_name)
        elif add_idt: f_idt.append(full_name)

    return f_doctest, function


def process_class(c_name, obj, c_md, c_mdt, c_idt, c_has_doctest):

    """ Extracts information about the class regarding documentation.
    It is assumed that the function calling this subroutine has already
    checked that the class is valid. """

    c = False
    c_dt = False
    # Get the line number of class
    try:
      line_no = inspect.getsourcelines(obj)[1]
    except IOError:
      # Raised when source does not exist
      # which means the class is not there.
      return c_dt, c

    c = True
    full_name = "LINE %d: %s" % (line_no, c_name)
    if not obj.__doc__: c_md.append(full_name)
    elif not '>>>' in obj.__doc__: c_mdt.append(full_name)
    else:
        c_dt =  True
        c_has_doctest.append(full_name)
        # indirect doctest
        if _is_indirect(c_name, obj.__doc__):
            c_idt.append(full_name)
    return c_dt, c

def coverage(module_path, verbose=False):

    """ Given a module path, builds an index of all classes and functions
    contained. It then goes through each of the classes/functions to get
    the docstring and doctest coverage of the module. """

    # Import the package and find membmers
    m = None
    try:
        __import__(module_path)
        m = sys.modules[module_path]
    except Exception, a:
        # Most likely cause, absence of __init__
        print module_path + ' could not be loaded due to, \"' + a.args[0] + '\"'
        return 0, 0


    c_skipped = []
    c_md = []
    c_mdt = []
    c_has_doctest = []
    c_idt = []
    classes = 0
    c_doctests = 0

    f_skipped = []
    f_md = []
    f_mdt = []
    f_has_doctest = []
    f_idt = []
    functions = 0
    f_doctests = 0

    skip_members = ['__abstractmethods__']

    # Get the list of members
    m_members = dir(m)
    for member in m_members:

        # Check for skipped functions first, they throw nasty errors
        # when combined with getattr
        if member in skip_members: continue

        # Identify if the member (class/def) a part of this module
        obj = getattr(m, member)
        obj_mod = inspect.getmodule(obj)

        # Function not a part of this module
        if not obj_mod or not obj_mod.__name__ == module_path:
          continue

        # If it's a function
        if inspect.isfunction(obj) or inspect.ismethod(obj):

            f_dt, f = process_function(member, '',  obj, module_path, f_skipped, f_md, f_mdt, f_idt, f_has_doctest, skip_members)
            if f: functions += 1
            if f_dt: f_doctests += 1

        # If it's a class, look at it's methods too
        elif inspect.isclass(obj):

            # Process the class first
            c_dt, c = process_class(member, obj, c_md, c_mdt, c_idt, c_has_doctest)
            if c: classes += 1
            if c_dt: c_doctests += 1

            # Iterate through it's members
            for f_name in dir(obj):

                if f_name in skip_members: continue

                # Identify the module of the current class member
                f_obj = getattr(obj, f_name)
                obj_mod = inspect.getmodule(f_obj)

                # Function not a part of this module
                if not obj_mod or not obj_mod.__name__ == module_path:
                  continue

                # If it's a function
                if inspect.isfunction(f_obj) or inspect.ismethod(f_obj):

                  f_dt, f = process_function(f_name, member, obj, module_path, f_skipped, f_md, f_mdt, f_idt, f_has_doctest, skip_members)
                  if f: functions += 1
                  if f_dt: f_doctests += 1


    # Evaluate the percent coverage
    total_doctests = c_doctests + f_doctests
    total_members = classes + functions
    if total_members: score = 100 * float(total_doctests) / (total_members)
    else: score = 0
    score = int(score)

    print_coverage(module_path, classes, c_md, c_mdt, c_idt, functions, f_md, f_mdt, f_idt, score, total_doctests, total_members, verbose)


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

    usage = "usage: ./bin/doctest_coverage.py sympy/core"

    parser = OptionParser(
        description = __doc__,
        usage = usage,
    )

    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
            default=False)

    options, args = parser.parse_args()

    if len(args) == 0:
        parser.print_help()
    else:
        for file in args:
            file = os.path.normpath(file)
            print 'DOCTEST COVERAGE for %s' % (file)
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
