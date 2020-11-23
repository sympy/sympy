#!/usr/bin/env python

from subprocess import check_call

def run(*cmdline, cwd=None, env=None):
    """
    Run command in subprocess and get lines of output
    """
    return check_call(cmdline, encoding='utf-8', cwd=cwd, env=env)
