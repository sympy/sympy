# ----------------------------------------------------------------------------
# pyglet
# Copyright (c) 2006-2007 Alex Holkner
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions 
# are met:
#
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright 
#    notice, this list of conditions and the following disclaimer in
#    the documentation and/or other materials provided with the
#    distribution.
#  * Neither the name of the pyglet nor the names of its
#    contributors may be used to endorse or promote products
#    derived from this software without specific prior written
#    permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
# ----------------------------------------------------------------------------

'''
'''

__docformat__ = 'restructuredtext'
__version__ = '$Id: lib.py 1322 2007-10-23 12:58:03Z Alex.Holkner $'

import sys
import ctypes

__all__ = ['link_GL', 'link_GLU', 'link_AGL', 'link_GLX', 'link_WGL']

class MissingFunctionException(Exception):
    def __init__(self, name, requires=None, suggestions=None):
        msg = '%s is not exported by the available OpenGL driver.' % name
        if requires:
            msg += '  %s is required for this functionality.' % requires
        if suggestions:
            msg += '  Consider alternative(s) %s.' % ', '.join(suggestions)
        Exception.__init__(self, msg)

def missing_function(name, requires=None, suggestions=None):
    def MissingFunction(*args, **kwargs):
        raise MissingFunctionException(name, requires, suggestions)
    return MissingFunction

_int_types = (ctypes.c_int16, ctypes.c_int32)
if hasattr(ctypes, 'c_int64'):
    # Some builds of ctypes apparently do not have c_int64
    # defined; it's a pretty good bet that these builds do not
    # have 64-bit pointers.
    _int_types += (ctypes.c_int64,)
for t in _int_types:
    if ctypes.sizeof(t) == ctypes.sizeof(ctypes.c_size_t):
        c_ptrdiff_t = t

class c_void(ctypes.Structure):
    # c_void_p is a buggy return type, converting to int, so
    # POINTER(None) == c_void_p is actually written as
    # POINTER(c_void), so it can be treated as a real pointer.
    _fields_ = [('dummy', ctypes.c_int)]

class GLException(Exception):
    pass

def errcheck(result, func, arguments):
    from pyglet.gl import get_current_context
    context = get_current_context()
    if not context:
        raise GLException('No GL context; create a Window first')
    if not context._gl_begin:
        from pyglet.gl import glGetError, gluErrorString
        error = glGetError()
        if error:
            message = ctypes.cast(gluErrorString(error), ctypes.c_char_p).value
            raise GLException(message)
        return result

def errcheck_glbegin(result, func, arguments):
    from pyglet.gl import get_current_context
    context = get_current_context()
    if not context:
        raise GLException('No GL context; create a Window first')
    context._gl_begin = True
    return result

def errcheck_glend(result, func, arguments):
    from pyglet.gl import get_current_context
    context = get_current_context()
    if not context:
        raise GLException('No GL context; create a Window first')
    context._gl_begin = False
    return errcheck(result, func, arguments)

def decorate_function(func, name):
    from pyglet import options
    if options['debug_gl']:
        if name == 'glBegin':
            func.errcheck = errcheck_glbegin
        elif name == 'glEnd':
            func.errcheck = errcheck_glend
        elif name not in ('glGetError', 'gluErrorString') and \
             name[:3] not in ('glX', 'agl', 'wgl'):
            func.errcheck = errcheck

link_AGL = None
link_GLX = None
link_WGL = None

if sys.platform in ('win32', 'cygwin'):
    from pyglet.gl.lib_wgl import link_GL, link_GLU, link_WGL
elif sys.platform == 'darwin':
    from pyglet.gl.lib_agl import link_GL, link_GLU, link_AGL
else:
    from pyglet.gl.lib_glx import link_GL, link_GLU, link_GLX

