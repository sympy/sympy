class Geometry3dEntity(tuple):
    """The base class for any three-dimensional geometrical entity."""
    
    def __new__(cls, *args, **kwargs):
        return tuple.__new__(cls, args)
        
    def __getnewargs__(self):
        return tuple(self)

