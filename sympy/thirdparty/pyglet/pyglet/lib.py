'''Functions for loading dynamic libraries.

These extend and correct ctypes functions.
'''

__docformat__ = 'restructuredtext'
__version__ = '$Id: $'

import os
import sys

import ctypes
import ctypes.util

class LibraryLoader(object):
    def load_library(self, *names, **kwargs):
        '''Find and load a library.  
        
        More than one name can be specified, they will be tried in order.
        Platform-specific library names (given as kwargs) are tried first.

        Raises ImportError if library is not found.
        '''
        if 'framework' in kwargs and self.platform == 'darwin':
            return self.load_framework(kwargs['framework'])
        
        platform_names = kwargs.get(self.platform, [])
        if type(platform_names) in (str, unicode):
            platform_names = [platform_names]
        elif type(platform_names) is tuple:
            platform_names = list(platform_names)

        if self.platform == 'linux2':
            platform_names.extend(['lib%s.so' % n for n in names])

        platform_names.extend(names)
        for name in platform_names:
            try:
                return ctypes.cdll.LoadLibrary(name)
            except OSError:
                path = self.find_library(name)
                if path:
                    return ctypes.cdll.LoadLibrary(path)
        raise ImportError('Library "%s" not found.' % names[0])

    find_library = lambda self, name: ctypes.util.find_library(name)

    platform = sys.platform
    if platform == 'cygwin':
        platform = 'win32'

    def load_framework(self, path):
        raise RuntimeError("Can't load framework on this platform.")

class MachOLibraryLoader(LibraryLoader):
    def __init__(self):
        if 'LD_LIBRARY_PATH' in os.environ:
            self.ld_library_path = os.environ['LD_LIBRARY_PATH'].split(':')
        else:
            self.ld_library_path = []

        if 'DYLD_LIBRARY_PATH' in os.environ:
            self.dyld_library_path = os.environ['DYLD_LIBRARY_PATH'].split(':')
        else:
            self.dyld_library_path = []

        if 'DYLD_FALLBACK_LIBRARY_PATH' in os.environ:
            self.dyld_fallback_library_path = \
                os.environ['DYLD_FALLBACK_LIBRARY_PATH'].split(':')
        else:
            self.dyld_fallback_library_path = [
                os.path.expanduser('~/lib'),
                '/usr/local/lib',
                '/usr/lib']
 
    def find_library(self, path):
        '''Implements the dylib search as specified in Apple documentation:
        
        http://developer.apple.com/documentation/DeveloperTools/Conceptual/DynamicLibraries/Articles/DynamicLibraryUsageGuidelines.html
        '''

        libname = os.path.basename(path)
        if '/' in path:
            search_path = (
                [os.path.join(p, libname) \
                    for p in self.dyld_library_path] +
                [path] + 
                [os.path.join(p, libname) \
                    for p in self.dyld_fallback_library_path])
        else:
            search_path = (
                [os.path.join(p, libname) \
                    for p in self.ld_library_path] +
                [os.path.join(p, libname) \
                    for p in self.dyld_library_path] +
                [path] + 
                [os.path.join(p, libname) \
                    for p in self.dyld_fallback_library_path])

        for path in search_path:
            if os.path.exists(path):
                return path

        return None

    def find_framework(self, path):
        '''Implement runtime framework search as described by:

        http://developer.apple.com/documentation/MacOSX/Conceptual/BPFrameworks/Concepts/FrameworkBinding.html
        '''

        # e.g. path == '/System/Library/Frameworks/OpenGL.framework'
        #      name == 'OpenGL'
        # return '/System/Library/Frameworks/OpenGL.framework/OpenGL'
        name = os.path.splitext(os.path.split(path)[1])[0]

        realpath = os.path.join(path, name)
        if os.path.exists(realpath):
            return realpath

        for dir in ('/Library/Frameworks',
                    '/System/Library/Frameworks'):
            realpath = os.path.join(dir, '%s.framework' % name, name)
            if os.path.exists(realpath):
                return realpath

        return None

    def load_framework(self, path):
        realpath = self.find_framework(path)
        if realpath:
            return ctypes.cdll.LoadLibrary(realpath)

        raise ImportError("Can't find framework %s." % path)

if sys.platform == 'darwin':
    loader = MachOLibraryLoader()
else:
    loader = LibraryLoader()
load_library = loader.load_library
