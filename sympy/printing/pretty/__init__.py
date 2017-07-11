"""ASCII-ART 2D pretty-printer"""

from .pretty import pager_print, pprint, pprint_try_use_unicode, \
    pprint_use_unicode, pretty, pretty_print

# if unicode output is available -- let's use it
pprint_try_use_unicode()
