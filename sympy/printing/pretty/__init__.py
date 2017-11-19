"""ASCII-ART 2D pretty-printer"""

__all__ = []

from .pretty import (
    pager_print, pretty, pretty_print, pprint,
    pprint_use_unicode, pprint_try_use_unicode
)
__all__ += [
    "pager_print", "pretty", "pretty_print", "pprint",
    "pprint_use_unicode", "pprint_try_use_unicode"
]

# if unicode output is available -- let's use it
pprint_try_use_unicode()
