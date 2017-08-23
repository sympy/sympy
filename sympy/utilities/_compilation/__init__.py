# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)
""" This sub-module is private, i.e. external code should not depend on it.

These functions are used by tests run as part of continuous integration.
Once the implementation is mature (it should support the major
platforms: Windows, OS X & Linux) it may become official API which
 may be relied upon by downstream libraries. Until then API may break
without prior notice.

TODO:
- (optionally) clean up after tempfile.mkdtemp()
- cross-platform testing
- caching of compiler choice and intermediate files

"""

from .compilation import compile_link_import_strings, compile_run_strings


def has_fortran():
    if not hasattr(has_fortran, 'result'):
        try:
            (stdout, stderr), info = compile_run_strings(
                [('main.f90', (
                    'program foo\n'
                    'print *, "hello world"\n'
                    'end program'
                ))]
            )
            assert 'hello world' in stdout
            assert stderr == ''
        except:
            has_fortran.result = False
        else:
            has_fortran.result = True
    return has_fortran.result


def has_c():
    if not hasattr(has_c, 'result'):
        try:
            (stdout, stderr), info = compile_run_strings(
                [('main.c', (
                    '#include <stdio.h>\n'
                    'int main(){\n'
                    'printf("hello world\\n");\n'
                    '}'
                ))]
            )
            assert 'hello world' in stdout
            assert stderr == ''
        except:
            has_c.result = False
        else:
            has_c.result = True
    return has_c.result


def has_cxx():
    if not hasattr(has_cxx, 'result'):
        try:
            (stdout, stderr), info = compile_run_strings(
                [('main.cxx', (
                    '#include <iostream>\n'
                    'int main(){\n'
                    'std::cout << "hello world" << std::endl;\n'
                    '}'
                ))]
            )
            assert 'hello world' in stdout
            assert stderr == ''
        except:
            has_cxx.result = False
        else:
            has_cxx.result = True
    return has_cxx.result
