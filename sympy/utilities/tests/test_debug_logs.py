"""
A way to check if debugging works.
"""
from __future__ import print_function
import subprocess
import sys
import os
import inspect
from os.path import abspath, dirname, join, normpath
from contextlib import contextmanager

@contextmanager
def environ(env):
    """
    set environment variables
    :param env: environment variable to modify
    """
    original_environ = os.environ.copy()
    os.environ.update(env)
    yield
    os.environ = original_environ # reset the environment variable to prevent side effect

from sympy.utilities.pytest import XFAIL


def test_debug_logs():
    my_filename = abspath(inspect.getfile(inspect.currentframe()))
    my_dirname = dirname(my_filename)
    diagnose_debug_logs_filename = join(my_dirname, 'diagnose_debug_logs.py')
    diagnose_debug_logs_filename = normpath(diagnose_debug_logs_filename)
    with environ({'SYMPY_DEBUG': 'True'}):
        output = subprocess.check_output(['python', diagnose_debug_logs_filename])
    assert check_debug_content(output), "There are problems or debug flag is not set \n"

def check_debug_content(log):
    """
    check for token that is representative of our log
    :param log: delog log as an input
    return: True if it is debug, False otherwise
    """
    lower_case_log = log.lower()
    if b"debug" in lower_case_log or b"rule" in lower_case_log:
        return True
    return False
