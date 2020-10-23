#!/usr/bin/env python

from subprocess import check_output

def run(*cmdline, cwd=None, env=None):
    """
    Run command in subprocess and get lines of output
    """
    return check_output(cmdline, encoding='utf-8', cwd=cwd, env=env).splitlines()
