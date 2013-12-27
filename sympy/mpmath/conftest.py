# The py library is part of the "py.test" testing suite (python-codespeak-lib
# on Debian), see http://codespeak.net/py/

import py

#this makes py.test put mpath directory into the sys.path, so that we can
#"import mpmath" from tests nicely
rootdir = py.magic.autopath().dirpath()
