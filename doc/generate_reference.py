#! /usr/bin/python2.5

# See the README in this directory how to generate the documentation.

# This script generates the reference, however, it currently doesn't produce a
# nicely polished documentation.

# You need to run python2.5 with this, because python2.4 has some weird
# Exception hierarchy classes that causes Exceptions to be ignored by this
# script

import types
import os
import sys

def isclass(x):
    from inspect import isclass as _isclass
    return _isclass(x)

def ismethod(x):
    from inspect import ismethod as _ismethod
    if _ismethod(x) or isinstance(x, (types.MethodType, types.FunctionType,
                    types.UnboundMethodType)) or str(type(x)) in [
                    "<type 'classmethod'>", "<type 'staticmethod'>"]:
        return True
    else:
        return False

def get_method_name(x):
    if hasattr(x, "__name__"):
        return x.__name__
    # This is a static method, don't know how to read the name.
    return None

def get_method_args(x):
    from inspect import getsource
    try:
        s = getsource(x)
    except TypeError:
        return ""
    s = s[s.find("("):s.find(":")]
    assert s is not None
    return s

def getdoc(x):
    from inspect import getdoc as _getdoc
    s = _getdoc(x)
    return s

class Parser(object):

    def __init__(self):
        self.modules = {}

    def generate_doc(self, module, outdir="/tmp/", importpath = ""):
        """
        Takes the "module" string and generates a rst for this module.

        modules is just the name, like "sympy", i.e. without a path.
        """

        print "Generating API documentation ... (this may take a while)"
        sys.path.insert(0, importpath)
        m = __import__(module)
        self.main_module = module
        self.handle_module(m)
        print "  saving..."
        for x in self.modules:
            if x.startswith(module):
                f = open(outdir+x+".txt", "w")
                f.write(self.modules[x])
        print "API generated in %s." % (outdir)


    def handle_module(self, mod):
        mname = mod.__name__
        if mname in self.modules:
            # we already handled this module
            return
        self.modules[mname] = ""

        # if you want to get the top level modules without "sympy.", uncomment
        # this:
        #if mname.startswith(self.main_module):
        #    #strip the "sympy.":
        #    s = ".. module:: %s\n\n" % mname[len(self.main_module)+1:]
        #else:
        #    s = ".. module:: %s\n\n" % mname

        s = "=" * len(mname) + "\n"
        s += mname + "\n"
        s += "=" * len(mname) + "\n" + "\n"

        s += ".. module:: %s\n\n" % mname
        if hasattr(mod, __file__):
            s += "filename: %s\n" % mod.__file__
        for x in mod.__dict__.values():
            if isinstance(x, types.ModuleType):
                self.handle_module(x)
            elif x is None:
                pass
            elif isinstance(x, (int, float, str, list, dict)):
                # skip these
                pass
            elif str(type(x)) == "<type 'classobj'>":
                # old style classes
                pass
            elif hasattr(x, "__class__"):
                s += self.handle_class(x)
            else:
                print "  Ignored:", type(x), x
        self.modules[mod.__name__] = s

    def handle_class(self, cls):
        if hasattr(cls, "__name__"):
            s = "\n.. class:: %s\n" % cls.__name__
        else:
            return ""

        # Uncomment this to generate class docstrings too:
        # unfortunately, sphinx fails to read them, so we need to fix sympy
        # first.
        #doc = getdoc(cls)
        #if doc is not None:
        #    s += doc
        #    s += "\n"

        if hasattr(cls, "__dict__"):
            for x in cls.__dict__.values():
                if isinstance(x, types.ModuleType):
                    self.handle_module(x)
                elif str(type(x)) == "<type 'classobj'>":
                    # old style classes
                    pass
                elif x is None:
                    pass
                elif ismethod(x):
                    s += self.handle_method(x)
                elif str(type(x)) == "<type 'property'>":
                    pass
                elif isinstance(x, (int, float, str, list, tuple, dict)):
                    # skip these
                    pass
                elif hasattr(x, "__class__"):
                    # ignore nested classes
                    pass
                else:
                    print "    Ignored in class:", type(x), x

        return s

    def handle_method(self, m):
        mname = get_method_name(m)
        if mname is None:
            s = ""
        else:
            s = "\n.. method:: %s%s\n\n" % (mname, get_method_args(m))
            doc = getdoc(m)
            if doc is not None:
                s += doc
                s += "\n"
        return s

Parser().generate_doc("sympy", importpath = "..", outdir="api/modules/")
