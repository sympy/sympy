
class BasicType(type):

    classnamespace = dict()
    def __init__(cls,*args,**kws):
        if not cls.undefined_Function:
            # make Basic subclasses available as attributes
            # set is_<classname> = True and False for other classes
            n = cls.__name__
            c = BasicType.classnamespace.get(n)
            if c is None:
                setattr(cls, 'is_' + n, True)
                for k, v in BasicType.classnamespace.items():
                    setattr(v, 'is_' + n, False)
                BasicType.classnamespace[n] = cls
            else:
                print 'Ignoring redefinition of %s: %s defined earlier than %s' % (n, c, cls)
        type.__init__(cls, *args, **kws)

    def __getattr__(cls, name):
        try: return BasicType.classnamespace[name]
        except KeyError: pass
        raise AttributeError("'%s' object has no attribute '%s'"%
                             (cls.__name__, name))


class Basic(object):

    __metaclass__ = BasicType
    undefined_Function = False

    @staticmethod
    def sympify(a):
        if isinstance(a, Basic):
            return a
        elif isinstance(a, bool):
            raise NotImplementedError("bool support")
        elif isinstance(a, (int, long)):
            return Basic.Integer(a)
        raise ValueError("%s is NOT a valid SymPy expression" % `a`)

    def __repr__(self):
        if isinstance(self, type):
            return self.__class__.torepr(self)
        return self.torepr()


class Atom(Basic):

    canonical = lambda self: self

    def torepr(self):
        return '%s()' % (self.__class__.__name__)

class Composite(Basic):

    def torepr(self):
        return '%s(%s)' % (self.__class__.__name__,', '.join(map(repr, self)))

class MutableCompositeDict(Composite, dict):

    # constructor methods
    def __new__(cls, *args, **options):
        """
        To make MutableClass immutable, execute
          obj.__class__ = Class
        """
        obj = dict.__new__(cls)
        [obj.update(a) for a in args]
        return obj

    def __init__(self, *args, **options):
        pass

    # representation methods
    def torepr(self):
        return '%s(%s)' % (self.__class__.__name__, dict(self))
