class Geometry3dEntity(tuple):

    def __new__(cls, *args, **kwargs):
        return tuple.__new__(cls, args)

    def __getnewargs__(self):
        return tuple(self)
