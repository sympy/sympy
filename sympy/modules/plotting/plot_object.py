class PlotObject(object):
    """
    Base class for all objects which can be displayed on a Plot.
    This includes basis vectors, boundary boxes, and most importantly
    PlotFunctions.
    """
    visible = True

    def do_render(self):
        """
        Do not override in child class.
        """
        if self.visible:
            self.render()

    def render(self):
        """
        OpenGL code to display this plot object.
        """
        raise Exception("PlotObject class is abstract. Override render in subclass.")
