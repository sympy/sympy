
import sys

from iterables import make_list, flatten
if sys.version_info[1] < 5:
    from iterables import any, all
else:
    any = any
    all = all
