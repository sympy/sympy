# $Id:$

'''Minimal Windows COM interface.

Allows pyglet to use COM interfaces on Windows without comtypes.  Unlike
comtypes, this module does not provide property interfaces, read typelibs,
nice-ify return values or permit Python implementations of COM interfaces.  We
don't need anything that sophisticated to work with DirectX.

All interfaces should derive from IUnknown (defined in this module).  The
Python COM interfaces are actually pointers to the implementation (take note
when translating methods that take an interface as argument).

Interfaces can define methods::

    class IDirectSound8(com.IUnknown):
        _methods_ = [
            ('CreateSoundBuffer', com.STDMETHOD()),
            ('GetCaps', com.STDMETHOD(LPDSCAPS)),
            ...
        ]

Only use STDMETHOD or METHOD for the method types (not ordinary ctypes
function types).  The 'this' pointer is bound automatically... e.g., call::

    device = IDirectSound8()
    DirectSoundCreate8(None, ctypes.byref(device), None)

    caps = DSCAPS()
    device.GetCaps(caps)

Because STDMETHODs use HRESULT as the return type, there is no need to check
the return value.

Don't forget to manually manage memory... call Release() when you're done with
an interface.
'''

import ctypes

class GUID(ctypes.Structure):
    _fields_ = [
        ('Data1', ctypes.c_ulong),
        ('Data2', ctypes.c_ushort),
        ('Data3', ctypes.c_ushort),
        ('Data4', ctypes.c_ubyte * 8)
    ]

    def __init__(self, l, w1, w2, b1, b2, b3, b4, b5, b6, b7, b8):
        self.Data1 = l
        self.Data2 = w1
        self.Data3 = w2
        self.Data4[:] = (b1, b2, b3, b4, b5, b6, b7, b8)

LPGUID = ctypes.POINTER(GUID)
IID = GUID
REFIID = ctypes.POINTER(IID)

class METHOD(object):
    '''COM method.'''
    def __init__(self, restype, *args):
        self.restype = restype
        self.argtypes = args

    def get_field(self):
        return ctypes.WINFUNCTYPE(self.restype, *self.argtypes)

class STDMETHOD(METHOD):
    '''COM method with HRESULT return value.'''
    def __init__(self, *args):
        super(STDMETHOD, self).__init__(ctypes.HRESULT, *args)

class COMMethodInstance(object):
    '''Binds a COM interface method.'''
    def __init__(self, name, i, method):
        self.name = name
        self.i = i
        self.method = method

    def __get__(self, obj, tp):
        if obj is not None:
            return lambda *args: \
                self.method.get_field()(self.i, self.name)(obj, *args)
        raise AttributeError()

class COMInterface(ctypes.Structure):
    '''Dummy struct to serve as the type of all COM pointers.'''
    _fields_ = [
        ('lpVtbl', ctypes.c_void_p),
    ]

class InterfaceMetaclass(type(ctypes.POINTER(COMInterface))):
    '''Creates COM interface pointers.'''
    def __new__(cls, name, bases, dct):
        methods = []
        for base in bases[::-1]:
            methods.extend(base.__dict__.get('_methods_', ()))
        methods.extend(dct.get('_methods_', ()))

        for i, (n, method) in enumerate(methods):
            dct[n] = COMMethodInstance(n, i, method)

        dct['_type_'] = COMInterface

        return super(InterfaceMetaclass, cls).__new__(cls, name, bases, dct)

class Interface(ctypes.POINTER(COMInterface)):
    '''Base COM interface pointer.'''
    __metaclass__ = InterfaceMetaclass

class IUnknown(Interface):
    _methods_ = [
        ('QueryInterface', STDMETHOD(REFIID, ctypes.c_void_p)),
        ('AddRef', METHOD(ctypes.c_int)),
        ('Release', METHOD(ctypes.c_int))
    ]

