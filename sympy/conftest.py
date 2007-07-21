import py
#this makes py.test put sympy directory into the sys.path, so that we can
#"import sympy" from tests nicely
rootdir = py.magic.autopath().dirpath()
