from sympy.core.decorators import deprecated
from sympy.core.core import BasicMeta, Registry, all_classes


class ClassRegistry(Registry):
    """
    Namespace for SymPy classes

    This is needed to avoid problems with cyclic imports.
    To get a SymPy class, use `C.<class_name>` e.g. `C.Rational`, `C.Add`.

    For performance reasons, this is coupled with a set `all_classes` holding
    the classes, which should not be modified directly.
    """
    __slots__ = []

    def __setattr__(self, name, cls):
        Registry.__setattr__(self, name, cls)
        all_classes.add(cls)

    def __delattr__(self, name):
        cls = getattr(self, name)
        Registry.__delattr__(self, name)
        # The same class could have different names, so make sure
        # it's really gone from C before removing it from all_classes.
        if cls not in self.__class__.__dict__.itervalues():
            all_classes.remove(cls)

    @deprecated(
        feature='C, including its class ClassRegistry,',
        last_supported_version='0.7.7',
        useinstead='direct imports from the defining module',
        issue=7124,
        deprecated_since_version='0.7.7')
    def __getattr__(self, name):
        return any(cls.__name__ == name for cls in all_classes)


C = ClassRegistry()
C.BasicMeta = BasicMeta
