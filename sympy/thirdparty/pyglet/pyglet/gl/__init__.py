# ----------------------------------------------------------------------------
# pyglet
# Copyright (c) 2006-2008 Alex Holkner
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
#  * Neither the name of pyglet nor the names of its
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

'''OpenGL and GLU interface.

This package imports all OpenGL, GLU and registered OpenGL extension
functions.  Functions have identical signatures to their C counterparts.  For
example::

    from pyglet.gl import *
    
    # [...omitted: set up a GL context and framebuffer]
    glBegin(GL_QUADS)
    glVertex3f(0, 0, 0)
    glVertex3f(0.1, 0.2, 0.3)
    glVertex3f(0.1, 0.2, 0.3)
    glEnd()

OpenGL is documented in full at the `OpenGL Reference Pages`_.  

The `OpenGL Programming Guide`_ is a popular reference manual organised by
topic.  The free online version documents only OpenGL 1.1.  `Later editions`_
cover more recent versions of the API and can be purchased from a book store.

.. _OpenGL Reference Pages: http://www.opengl.org/documentation/red_book/
.. _OpenGL Programming Guide: http://fly.cc.fer.hr/~unreal/theredbook/
.. _Later editions: http://www.opengl.org/documentation/red_book/

The following subpackages are imported into this "mega" package already (and
so are available by importing ``pyglet.gl``):

``pyglet.gl.gl``
    OpenGL
``pyglet.gl.glu``
    GLU
``pyglet.gl.gl.glext_arb``
    ARB registered OpenGL extension functions
``pyglet.gl.gl.glext_missing``
    ARB registered OpenGL extension functions not included in the ARB C header

These subpackages are also available, but are not imported into this namespace
by default:

``pyglet.gl.glext_nv``
    nVidia OpenGL extension functions
``pyglet.gl.agl``
    AGL (Mac OS X OpenGL context functions)
``pyglet.gl.glx``
    GLX (Linux OpenGL context functions)
``pyglet.gl.glxext_arb``
    ARB registered GLX extension functions
``pyglet.gl.glxext_nv``
    nvidia GLX extension functions
``pyglet.gl.wgl``
    WGL (Windows OpenGL context functions)
``pyglet.gl.wglext_arb``
    ARB registered WGL extension functions
``pyglet.gl.wglext_nv``
    nvidia WGL extension functions

The information modules are provided for convenience, and are documented
below.
'''

__docformat__ = 'restructuredtext'
__version__ = '$Id: __init__.py 2175 2008-08-15 04:29:16Z Alex.Holkner $'

from pyglet.gl.lib import GLException
from pyglet.gl.gl import *
from pyglet.gl.glu import *
from pyglet.gl.glext_arb import *
from pyglet.gl.glext_missing import *
from pyglet.gl import gl_info

import sys as _sys
_is_epydoc = hasattr(_sys, 'is_epydoc') and _sys.is_epydoc

# List of contexts currently in use, so we can create new contexts that
# share objects with.  Remember to remove from this list when context is
# destroyed.
_contexts = []

#: The active OpenGL context.
#:
#: You can change the current context by calling `Context.set_current`; do not
#: modify this global.
#:
#: :type: `Context`
#:
#: :since: pyglet 1.1
current_context = None

def get_current_context():
    '''Return the active OpenGL context.

    You can change the current context by calling `Context.set_current`.

    :deprecated: Use `current_context`

    :rtype: `Context`
    :return: the context to which OpenGL commands are directed, or None
        if there is no selected context.
    '''
    return current_context

class Config(object):
    '''Graphics configuration.

    A GLConfig stores the preferences for OpenGL attributes such as the
    number of auxilliary buffers, size of the colour and depth buffers,
    double buffering, stencilling, multi- and super-sampling, and so on.

    Different platforms support a different set of attributes, so these
    are set with a string key and a value which is integer or boolean.

    See also `pyglet.window.Screen.get_best_config` and 
    `pyglet.window.Screen.get_matching_configs`.

    :Ivariables:
        `double_buffer` : bool
            Specify the presence of a back-buffer for every color buffer.
        `stereo` : bool
            Specify the presence of separate left and right buffer sets.
        `buffer_size` : int
            Total bits per sample per color buffer.
        `aux_buffers` : int
            The number of auxilliary color buffers.
        `sample_buffers` : int
            The number of multisample buffers.
        `samples` : int
            The number of samples per pixel, or 0 if there are no multisample
            buffers.
        `red_size` : int
            Bits per sample per buffer devoted to the red component.
        `green_size` : int
            Bits per sample per buffer devoted to the green component.
        `blue_size` : int
            Bits per sample per buffer devoted to the blue component.
        `alpha_size` : int
            Bits per sample per buffer devoted to the alpha component.
        `depth_size` : int
            Bits per sample in the depth buffer.
        `stencil_size` : int
            Bits per sample in the stencil buffer.
        `accum_red_size` : int
            Bits per pixel devoted to the red component in the accumulation
            buffer.
        `accum_green_size` : int
            Bits per pixel devoted to the green component in the accumulation
            buffer.
        `accum_blue_size` : int
            Bits per pixel devoted to the blue component in the accumulation
            buffer.
        `accum_alpha_size` : int
            Bits per pixel devoted to the alpha component in the accumulation
            buffer.
    '''

    _attribute_names = [
        'double_buffer',
        'stereo',
        'buffer_size',
        'aux_buffers',
        'sample_buffers',
        'samples',
        'red_size',
        'green_size',
        'blue_size',
        'alpha_size',
        'depth_size',
        'stencil_size',
        'accum_red_size',
        'accum_green_size',
        'accum_blue_size',
        'accum_alpha_size',
    ]

    def __init__(self, **kwargs):
        '''Create a template config with the given attributes.

        Specify attributes as keyword arguments, for example::

            template = Config(double_buffer=True)

        '''
        for name in self._attribute_names:
            if name in kwargs:
                setattr(self, name, kwargs[name])
            else:
                setattr(self, name, None)

    def get_gl_attributes(self):
        '''Return a list of attributes set on this config.

        :rtype: list of tuple (name, value)
        :return: All attributes, with unset attributes having a value of
            ``None``.
        '''
        return [(name, getattr(self, name)) for name in self._attribute_names]

    def create_context(self, share):
        '''Create a GL context that satisifies this configuration.

        :Parameters:
            `share` : `Context`
                If not None, a context with which to share objects with.

        :rtype: `Context`
        :return: The new context.
        '''
        raise ConfigException(
            'This config is not complete.  Use Screen.get_matching_configs')

    def is_complete(self):
        '''Determine if this config is complete and able to create a context.

        Configs created directly are not complete, they can only serve
        as templates for retrieving a supported config from the system.
        For example, `pyglet.window.Screen.get_matching_configs` returns
        complete configs.

        :rtype: bool
        :return: True if the config is complete and can create a context.
        '''
        return False

    def __repr__(self):
        import pprint
        return '%s(%s)' % (self.__class__.__name__, 
                           pprint.pformat(self.get_gl_attributes()))
                                          

class ObjectSpace(object):
    def __init__(self):
        # Textures and buffers scheduled for deletion the next time this
        # object space is active.
        self._doomed_textures = []
        self._doomed_buffers = []

class Context(object):
    '''OpenGL context for drawing.

    Windows in pyglet each have their own GL context.  This class boxes
    the context in a platform-independent manner.  Applications will have
    no need to deal with contexts directly.

    :Ivariables:
        `object_space` : `ObjectSpace`
            An object which is shared between all contexts that share
            GL objects.

    '''

    #: Context share behaviour indicating that objects should not be
    #: shared with existing contexts.
    CONTEXT_SHARE_NONE = None

    #: Context share behaviour indicating that objects are shared with
    #: the most recently created context (the default).
    CONTEXT_SHARE_EXISTING = 1
    
    # Used for error checking, True if currently within a glBegin/End block.
    # Ignored if error checking is disabled.
    _gl_begin = False

    # gl_info.GLInfo instance, filled in on first set_current
    _info = None

    # List of (attr, check) for each driver/device-specific workaround that is
    # implemented.  The `attr` attribute on this context is set to the result
    # of evaluating `check(gl_info)` the first time this context is used.
    _workaround_checks = [
        # GDI Generic renderer on Windows does not implement
        # GL_UNPACK_ROW_LENGTH correctly.
        ('_workaround_unpack_row_length',
         lambda info: info.get_renderer() == 'GDI Generic'),

        # Reportedly segfaults in text_input.py example.
        ('_workaround_vbo',
         lambda info: info.get_renderer() == 'ATI Radeon X1600 OpenGL Engine'),

        # Some ATI cards on OS X start drawing from a VBO before it's written
        # to.  In these cases pyglet needs to call glFinish() to flush the
        # pipeline after updating a buffer but before rendering.
        ('_workaround_vbo_finish',
         lambda info: ('ATI' in info.get_renderer() and 
                       info.have_version(1, 5) and
                       _sys.platform == 'darwin')),
    ]

    def __init__(self, context_share=None):
        self.window = None
        _contexts.append(self)
        if context_share:
            assert context_share in _contexts
            self.object_space = context_share.object_space
        else:
            self.object_space = ObjectSpace()
    
    def __repr__(self):
        return '%s()' % self.__class__.__name__

    def set_current(self):
        global current_context
        assert self in _contexts
        current_context = self

        # Implement workarounds
        if not self._info:
            self._info = gl_info.GLInfo()
            self._info.set_active_context()
            for attr, check in self._workaround_checks:
                setattr(self, attr, check(self._info))

        # Release textures on this context scheduled for deletion
        if self.object_space._doomed_textures:
            textures = self.object_space._doomed_textures
            textures = (GLuint * len(textures))(*textures)
            glDeleteTextures(len(textures), textures)
            self.object_space._doomed_textures = []

        # Release buffers on this context scheduled for deletion
        if self.object_space._doomed_buffers:
            buffers = self.object_space._doomed_buffers
            buffers = (GLuint * len(buffers))(*buffers)
            glDeleteBuffers(len(buffers), buffers)
            self.object_space._doomed_buffers = []

    def destroy(self):
        '''Release the context.

        The context will not be useable after being destroyed.  Each platform
        has its own convention for releasing the context and the buffer(s)
        that depend on it in the correct order; this should never be called
        by an application.
        '''
        global current_context
        if current_context is self:
            current_context = None
            gl_info.remove_active_context()

            # Switch back to shadow context.
            if _shadow_window is not None:
                _shadow_window.switch_to()
        _contexts.remove(self)

    def delete_texture(self, texture_id):
        '''Safely delete a texture belonging to this context.

        Usually, the texture is released immediately using
        ``glDeleteTextures``, however if another context that does not share
        this context's object space is currently active, the deletion will
        be deferred until an appropriate context is activated.

        :Parameters:
            `texture_id` : int
                The OpenGL name of the texture to delete.

        '''
        if self.object_space is current_context.object_space:
            id = GLuint(texture_id)
            glDeleteTextures(1, id)
        else:
            self.object_space._doomed_textures.append(texture_id)

    def delete_buffer(self, buffer_id):
        '''Safely delete a buffer object belonging to this context.

        This method behaves similarly to `delete_texture`, though for
        ``glDeleteBuffers`` instead of ``glDeleteTextures``.

        :Parameters:
            `buffer_id` : int
                The OpenGL name of the buffer to delete.

        :since: pyglet 1.1
        '''
        if self.object_space is current_context.object_space and False:
            id = GLuint(buffer_id)
            glDeleteBuffers(1, id)
        else:
            self.object_space._doomed_buffers.append(buffer_id)

class ContextException(Exception):
    pass

class ConfigException(Exception):
    pass

import pyglet as _pyglet

if _pyglet.options['debug_texture']:
    _debug_texture_total = 0
    _debug_texture_sizes = {}
    _debug_texture = None

    def _debug_texture_alloc(texture, size):
        global _debug_texture_total

        _debug_texture_sizes[texture] = size
        _debug_texture_total += size

        print '%d (+%d)' % (_debug_texture_total, size)

    def _debug_texture_dealloc(texture):
        global _debug_texture_total

        size = _debug_texture_sizes[texture]
        del _debug_texture_sizes[texture]
        _debug_texture_total -= size

        print '%d (-%d)' % (_debug_texture_total, size)

    _glBindTexture = glBindTexture
    def glBindTexture(target, texture):
        global _debug_texture
        _debug_texture = texture
        return _glBindTexture(target, texture)

    _glTexImage2D = glTexImage2D
    def glTexImage2D(target, level, internalformat, width, height, border,
                     format, type, pixels):
        try:
            _debug_texture_dealloc(_debug_texture)
        except KeyError:
            pass

        if internalformat in (1, GL_ALPHA, GL_INTENSITY, GL_LUMINANCE):
            depth = 1
        elif internalformat in (2, GL_RGB16, GL_RGBA16):
            depth = 2
        elif internalformat in (3, GL_RGB):
            depth = 3
        else:
            depth = 4 # Pretty crap assumption
        size = (width + 2 * border) * (height + 2 * border) * depth
        _debug_texture_alloc(_debug_texture, size)

        return _glTexImage2D(target, level, internalformat, width, height,
                             border, format, type, pixels)

    _glDeleteTextures = glDeleteTextures
    def glDeleteTextures(n, textures):
        if not hasattr(textures, '__len__'):
            _debug_texture_dealloc(textures.value)
        else:
            for i in range(n):
                _debug_texture_dealloc(textures[i].value)
        
        return _glDeleteTextures(n, textures)

def _create_shadow_window():
    global _shadow_window

    import pyglet
    if not pyglet.options['shadow_window'] or _is_epydoc:
        return
    
    from pyglet.window import Window
    _shadow_window = Window(width=1, height=1, visible=False)
    _shadow_window.switch_to()

    from pyglet import app
    app.windows.remove(_shadow_window)


_shadow_window = None

# Import pyglet.window now if it isn't currently being imported (this creates
# the shadow window).
if (not _is_epydoc and
    'pyglet.window' not in _sys.modules and 
    _pyglet.options['shadow_window']):
    # trickery is for circular import 
    _pyglet.gl = _sys.modules[__name__]
    import pyglet.window
