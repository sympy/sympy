from threading import Lock
from arg_parsing import parse_plot_args, parse_function_args
from plot_object import PlotObject
from plot_function import PlotFunction
from plot_window import PlotWindow
from bounding_box import BoundingBox
from grid_plane import GridPlane

class Plot(object):
    """
    Flexible interface for plotting SymPy functions.
    """

    default_width = 600
    default_height = 600

    def __init__(self, *args, **kwargs):
        self._render_object_lock = Lock() # see lock_begin/lock_end

        self.clear() # initialize _functions list
        self._clear_plotobjects() # initialize _plotobjects list

        self.width = kwargs.get('width', self.default_width)
        self.height = kwargs.get('height', self.default_height)
        self.wireframe = kwargs.get('wireframe', False)
        self.antialiasing = kwargs.get('antialiasing', True)
        self.vsync = kwargs.get('vsync', False)
        self.ortho = kwargs.get('ortho', False)

        self.bounding_box = BoundingBox()
        self.bounding_box.visible = kwargs.get('bounding_box', False)
        self._append_plotobject(self.bounding_box)

        self.grid = kwargs.get('grid', None)
        if isinstance(self.grid, str):
            self.grid = (self.grid,)
        if self.grid != None:
            self.grid = GridPlane(*self.grid)
            self._append_plotobject(self.grid)

        self._calculations_in_progress = 0        

        for f in parse_plot_args(*args):
            self.append(f)     

        self.window = None
        if kwargs.get('show', True):
            self.show()

    def show(self):
        """
        Displays a UI window representing the Plot.
        """
        if self.window != None and not self.window.has_exit and self.window.context != None:
            if self.window.visible:
                self.window.activate()
        else:
            self.window = None
            self.window = PlotWindow(self,
                                     wireframe=self.wireframe,
                                     antialiasing=self.antialiasing,
                                     width=self.width,
                                     height=self.height,
                                     vsync=self.vsync,
                                     ortho=self.ortho)

    def close(self):
        self.window.close()

    def saveimage(self, filepath, **kwargs):
        raise NotImplementedError()

    ### PlotFunction List Interfaces (for end-users) ###

    def clear(self):
        self.lock_begin()
        self._functions = {}
        self.lock_end()

    def __getitem__(self, i):
        return self._functions[i]

    def __setitem__(self, i, args):
        if not (isinstance(i, int) and i > 0):
            raise ValueError("Function index must be a positive integer.")
        if isinstance(args, PlotFunction):
            f = args
        else:
            if not isinstance(args, (tuple, list)): args = [args]
            self._calculations_in_progress += 1
            #print "Calculating..."
            f = parse_function_args(*args)
            #print "Finished."
            self._calculations_in_progress -= 1
            assert isinstance(f, PlotFunction)
        self.bounding_box.consider_function(f)
        self.lock_begin()
        self._functions[i] = f
        self.lock_end()

    def __delitem__(self, i):
        self.lock_begin()
        del self._functions[i] 
        self.lock_end()

    def firstavailableindex(self):
        i=1
        self.lock_begin()
        while i in self._functions: i += 1
        self.lock_end()
        return i

    def append(self, args):
        # synchronization handled in __setitem__
        self[self.firstavailableindex()] = args

    def __len__(self):
        return len(self._functions)

    def __iter__(self):
        return self._functions.itervalues()

    def __str__(self):
        s = ""
        if len(self._functions) == 0:
            s += "<empty>"
        else:
            self.lock_begin()
            s += "\n".join(["%s[%i]: %s" % ("", i, str(self._functions[i]))
                              for i in self._functions])
            self.lock_end()
        return s

    ### PlotObject List Interfaces (for internal use and intrepid hackery) ##

    def _clear_plotobjects(self):
        self.lock_begin()
        self._plotobjects = []
        self.lock_end()

    def _append_plotobject(self, o):
        assert isinstance(o, PlotObject)
        self.lock_begin()
        self._plotobjects.append(o)
        self.lock_end()

    def _remove_plotobject(self, o):
        self.lock_begin()
        self._plotobjects.remove(o)
        self.lock_end()

    ### Thread Synchronization ###

    def lock_begin(self):
        self._render_object_lock.acquire()

    def lock_end(self):
        self._render_object_lock.release()

