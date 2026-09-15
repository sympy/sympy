#!/usr/bin/env python
"""
Test that the installed modules in setup.py are up-to-date.

If this test fails, run

python bin/generate_test_list.py

and

python bin/generate_module_list.py

to generate the up-to-date test and modules list to put in setup.py.

"""
import generate_test_list
import generate_module_list

from get_sympy import path_hack
path_hack()

import setup

module_list = generate_module_list.generate_module_list()
test_list = generate_test_list.generate_test_list()

assert setup.modules == module_list, set(setup.modules).symmetric_difference(set(module_list))
assert setup.tests == test_list, set(setup.tests).symmetric_difference(set(test_list))
print("setup.py modules and tests are OK")
