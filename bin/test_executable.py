#!/usr/bin/env python
"""
Test that only executable files have an executable bit set
"""
from __future__ import print_function

import os
import sys

from get_sympy import path_hack
base_dir = path_hack()

def test_executable(path):
    if not os.path.isdir(path):
        if os.access(path, os.X_OK):
            with open(path, 'r') as f:
                if f.readline()[:2] != "#!":
                    exn_msg = "File at " + path + " either should not be executable or should have a shebang line"
                    raise SystemError(exn_msg)
    else:
        for file in os.listdir(path):
            test_executable(os.path.join(path, file))

test_executable(base_dir)
