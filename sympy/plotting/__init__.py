"""Plotting module that can plot 2D and 3D functions
"""
try:
    try:
        from ctypes import *
    except:
        raise ImportError("ctypes is required for plotting.\n'easy_install ctypes' or visit http://sourceforge.net/projects/ctypes/")

    def Plot(*args, **kwargs):
        import plot
        return plot.Plot(*args, **kwargs)

except Exception, e:
    def Plot(*args, **kwargs):
        raise e

from textplot import textplot
