"""
This module is a copy of python's random module but using
a Random instance for all supported methods. It is needed
for replicating results when using seeds.
"""

import random as std_random
import sys

_inst = std_random.Random()
_thismodule = sys.modules[__name__]

for method_name in std_random.__all__:
    in_inst = hasattr(_inst, method_name)
    parent = _inst if in_inst else std_random
    method = getattr(parent, method_name)
    setattr(_thismodule, method_name, method)
