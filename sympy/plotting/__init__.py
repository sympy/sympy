try:
    try:
        from ctypes import *
    except:
        raise ImportError("ctypes is required for plotting.\n'easy_install ctypes' or visit http://sourceforge.net/projects/ctypes/")

    import sys
    import os
    try:
        libdir = os.path.abspath(os.path.dirname(__file__))
        sys.path.insert(0, libdir)
    except:
        pass

    from plot import Plot

except Exception, e:
    class Plot(object):
        def __init__(*args, **kwargs):
            raise e
