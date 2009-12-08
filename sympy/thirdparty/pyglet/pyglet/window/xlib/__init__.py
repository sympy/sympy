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

__docformat__ = 'restructuredtext'
__version__ = '$Id: __init__.py 2262 2008-09-16 14:33:09Z Alex.Holkner $'

from ctypes import *
import unicodedata
import warnings

import pyglet
from pyglet.window import WindowException, NoSuchDisplayException, \
    MouseCursorException, Platform, Display, Screen, MouseCursor, \
    DefaultMouseCursor, ImageMouseCursor, BaseWindow, _PlatformEventHandler
from pyglet.window import event
from pyglet.window import key
from pyglet.window import mouse
from pyglet.event import EventDispatcher

from pyglet import gl
from pyglet.gl import gl_info
from pyglet.gl import glu_info
from pyglet.gl import glx
from pyglet.gl import glxext_arb
from pyglet.gl import glxext_mesa
from pyglet.gl import glx_info

import pyglet.window.xlib.xlib
from pyglet.window.xlib import cursorfont
try:
    import pyglet.window.xlib.xinerama
    _have_xinerama = True
except:
    _have_xinerama = False

try:
    import pyglet.window.xlib.xsync
    _have_xsync = True
except:
    _have_xsync = False

class mwmhints_t(Structure):
    _fields_ = [
        ('flags', c_uint32),
        ('functions', c_uint32),
        ('decorations', c_uint32),
        ('input_mode', c_int32),
        ('status', c_uint32)
    ]

XA_CARDINAL = 6 # Xatom.h:14

# Do we have the November 2000 UTF8 extension?
_have_utf8 = hasattr(xlib._lib, 'Xutf8TextListToTextProperty')

# symbol,ctrl -> motion mapping
_motion_map = {
    (key.UP, False):        key.MOTION_UP,
    (key.RIGHT, False):     key.MOTION_RIGHT,
    (key.DOWN, False):      key.MOTION_DOWN,
    (key.LEFT, False):      key.MOTION_LEFT,
    (key.RIGHT, True):      key.MOTION_NEXT_WORD,
    (key.LEFT, True):       key.MOTION_PREVIOUS_WORD,
    (key.HOME, False):      key.MOTION_BEGINNING_OF_LINE,
    (key.END, False):       key.MOTION_END_OF_LINE,
    (key.PAGEUP, False):    key.MOTION_PREVIOUS_PAGE,
    (key.PAGEDOWN, False):  key.MOTION_NEXT_PAGE,
    (key.HOME, True):       key.MOTION_BEGINNING_OF_FILE,
    (key.END, True):        key.MOTION_END_OF_FILE,
    (key.BACKSPACE, False): key.MOTION_BACKSPACE,
    (key.DELETE, False):    key.MOTION_DELETE,
}

# Set up error handler
def _error_handler(display, event):
    # By default, all errors are silently ignored: this has a better chance
    # of working than the default behaviour of quitting ;-)
    #
    # We've actually never seen an error that was our fault; they're always
    # driver bugs (and so the reports are useless).  Nevertheless, set
    # environment variable PYGLET_DEBUG_X11 to 1 to get dumps of the error
    # and a traceback (execution will continue).
    if pyglet.options['debug_x11']:
        event = event.contents
        buf = c_buffer(1024)
        xlib.XGetErrorText(display, event.error_code, buf, len(buf))
        print 'X11 error:', buf.value
        print '   serial:', event.serial
        print '  request:', event.request_code
        print '    minor:', event.minor_code
        print ' resource:', event.resourceid

        import traceback
        print 'Python stack trace (innermost last):'
        traceback.print_stack()
    return 0
_error_handler_ptr = xlib.XErrorHandler(_error_handler)
xlib.XSetErrorHandler(_error_handler_ptr)

class XlibException(WindowException):
    '''An X11-specific exception.  This exception is probably a programming
    error in pyglet.'''
    pass

class XlibMouseCursor(MouseCursor):
    drawable = False

    def __init__(self, cursor):
        self.cursor = cursor

class XlibPlatform(Platform):
    def __init__(self):
        self._displays = {}

    def get_display(self, name):
        if name not in self._displays:
            self._displays[name] = XlibDisplayDevice(name)
        return self._displays[name]

    def get_default_display(self):
        return self.get_display('')

class XlibDisplayDevice(Display):
    _display = None     # POINTER(xlib.Display)

    _x_im = None        # X input method
                        # TODO close _x_im when display connection closed.
    _enable_xsync = False

    def __init__(self, name):
        super(XlibDisplayDevice, self).__init__()

        self._display = xlib.XOpenDisplay(name)
        if not self._display:
            raise NoSuchDisplayException('Cannot connect to "%s"' % name)
        self.info = glx_info.GLXInfo(self._display)

        # Also set the default GLX display for future info queries
        glx_info.set_display(self._display.contents)

        self._fileno = xlib.XConnectionNumber(self._display)
        self._window_map = {}

        # Initialise XSync
        if _have_xsync:
            event_base = c_int()
            error_base = c_int()
            if xsync.XSyncQueryExtension(self._display, 
                                         byref(event_base),
                                         byref(error_base)):
                major_version = c_int()
                minor_version = c_int()
                if xsync.XSyncInitialize(self._display,
                                         byref(major_version),
                                         byref(minor_version)):
                    self._enable_xsync = True

    def fileno(self):
        return self._fileno

    def get_screens(self):
        x_screen = xlib.XDefaultScreen(self._display)
        if _have_xinerama and xinerama.XineramaIsActive(self._display):
            number = c_int()
            infos = xinerama.XineramaQueryScreens(self._display, 
                                                  byref(number))
            infos = cast(infos, 
                 POINTER(xinerama.XineramaScreenInfo * number.value)).contents
            result = []
            for info in infos:
                result.append(XlibScreen(self,
                                         x_screen,
                                         info.x_org,
                                         info.y_org,
                                         info.width,
                                         info.height,
                                         True))
            xlib.XFree(infos)
            return result
        else:
            # No xinerama
            screen_count = xlib.XScreenCount(self._display)
            result = []
            for i in range(screen_count):
                screen = xlib.XScreenOfDisplay(self._display, i)
                result.append(XlibScreen(self,
                                         i, 
                                         0, 0,
                                         screen.contents.width,
                                         screen.contents.height,
                                         False))
            # Move default screen to be first in list.
            s = result.pop(x_screen)
            result.insert(0, s)
            return result

class XlibScreen(Screen):
    def __init__(self, display, x_screen_id, x, y, width, height, xinerama):
        super(XlibScreen, self).__init__(x, y, width, height)
        self.display = display
        self._x_screen_id = x_screen_id
        self._xinerama = xinerama
 
    def get_matching_configs(self, template):
        x_display = self.display._display

        have_13 = self.display.info.have_version(1, 3)
        if have_13:
            config_class = XlibGLConfig13
        else:
            if 'ATI' in self.display.info.get_client_vendor():
                config_class = XlibGLConfig10ATI
            else:
                config_class = XlibGLConfig10
        
        # Construct array of attributes
        attrs = []
        for name, value in template.get_gl_attributes():
            attr = config_class.attribute_ids.get(name, None)
            if attr and value is not None:
                attrs.extend([attr, int(value)])

        if have_13:
            attrs.extend([glx.GLX_X_RENDERABLE, True])
        else:
            attrs.extend([glx.GLX_RGBA, True])

        if len(attrs):
            attrs.extend([0, 0])
            attrib_list = (c_int * len(attrs))(*attrs)
        else:
            attrib_list = None

        if have_13:
            elements = c_int()
            configs = glx.glXChooseFBConfig(x_display, self._x_screen_id,
                attrib_list, byref(elements))
            if not configs:
                return []

            configs = cast(configs, 
                           POINTER(glx.GLXFBConfig * elements.value)).contents

            result = [config_class(self, c) for c in configs]

            # Can't free array until all XlibGLConfig13's are GC'd.  Too much
            # hassle, live with leak. XXX
            #xlib.XFree(configs)

            return result
        else:
            try:
                return [config_class(self, attrib_list)]
            except gl.ContextException:
                return []

    def __repr__(self):
        return 'XlibScreen(screen=%d, x=%d, y=%d, ' \
               'width=%d, height=%d, xinerama=%d)' % \
            (self._x_screen_id, self.x, self.y, self.width, self.height,
             self._xinerama)

class XlibGLConfig(gl.Config):
    attribute_ids = {
        'buffer_size': glx.GLX_BUFFER_SIZE,
        'level': glx.GLX_LEVEL,     # Not supported
        'double_buffer': glx.GLX_DOUBLEBUFFER,
        'stereo': glx.GLX_STEREO,
        'aux_buffers': glx.GLX_AUX_BUFFERS,
        'red_size': glx.GLX_RED_SIZE,
        'green_size': glx.GLX_GREEN_SIZE,
        'blue_size': glx.GLX_BLUE_SIZE,
        'alpha_size': glx.GLX_ALPHA_SIZE,
        'depth_size': glx.GLX_DEPTH_SIZE,
        'stencil_size': glx.GLX_STENCIL_SIZE,
        'accum_red_size': glx.GLX_ACCUM_RED_SIZE,
        'accum_green_size': glx.GLX_ACCUM_GREEN_SIZE,
        'accum_blue_size': glx.GLX_ACCUM_BLUE_SIZE,
        'accum_alpha_size': glx.GLX_ACCUM_ALPHA_SIZE,
    }

    def create_context(self, share):
        context = self._create_glx_context(share)

        if context == glx.GLX_BAD_CONTEXT:
            raise gl.ContextException('Invalid context share')
        elif context == glx.GLXBadFBConfig:
            raise gl.ContextException('Invalid GL configuration')
        elif context < 0:
            raise gl.ContextException('Could not create GL context') 

        return XlibGLContext(self, context, share)

    def _create_glx_context(self, share):
        raise NotImplementedError('abstract')

    def is_complete(self):
        return True

    def get_visual_info(self):
        raise NotImplementedError('abstract')

class XlibGLConfig10(XlibGLConfig):
    def __init__(self, screen, attrib_list):
        self.screen = screen
        self._display = screen.display._display
        self._visual_info = glx.glXChooseVisual(self._display,
            screen._x_screen_id, attrib_list)
        if not self._visual_info:
            raise gl.ContextException('No conforming visual exists')

        for name, attr in self.attribute_ids.items():
            value = c_int()
            result = glx.glXGetConfig(self._display,
                self._visual_info, attr, byref(value))
            if result >= 0:
                setattr(self, name, value.value)
        self.sample_buffers = 0
        self.samples = 0

    def get_visual_info(self):
        return self._visual_info.contents

    def _create_glx_context(self, share):
        if share:
            return glx.glXCreateContext(self._display, self._visual_info,
                share._context, True)
        else:
            return glx.glXCreateContext(self._display, self._visual_info,
                None, True)

class XlibGLConfig10ATI(XlibGLConfig10):
    attribute_ids = XlibGLConfig.attribute_ids.copy()
    del attribute_ids['stereo']
    stereo = False

class XlibGLConfig13(XlibGLConfig):
    attribute_ids = XlibGLConfig.attribute_ids.copy()
    attribute_ids.update({
        'sample_buffers': glx.GLX_SAMPLE_BUFFERS,
        'samples': glx.GLX_SAMPLES,

        # Not supported in current pyglet API:
        'render_type': glx.GLX_RENDER_TYPE,
        'config_caveat': glx.GLX_CONFIG_CAVEAT,
        'transparent_type': glx.GLX_TRANSPARENT_TYPE,
        'transparent_index_value': glx.GLX_TRANSPARENT_INDEX_VALUE,
        'transparent_red_value': glx.GLX_TRANSPARENT_RED_VALUE,
        'transparent_green_value': glx.GLX_TRANSPARENT_GREEN_VALUE,
        'transparent_blue_value': glx.GLX_TRANSPARENT_BLUE_VALUE,
        'transparent_alpha_value': glx.GLX_TRANSPARENT_ALPHA_VALUE,

        # Used internally
        'x_renderable': glx.GLX_X_RENDERABLE,
    })

    def __init__(self, screen, fbconfig):
        super(XlibGLConfig13, self).__init__()
        self.screen = screen
        self._display = screen.display._display
        self._fbconfig = fbconfig
        for name, attr in self.attribute_ids.items():
            value = c_int()
            result = glx.glXGetFBConfigAttrib(
                self._display, self._fbconfig, attr, byref(value))
            if result >= 0:
                setattr(self, name, value.value)

    def get_visual_info(self):
        return glx.glXGetVisualFromFBConfig(
            self._display, self._fbconfig).contents

    def _create_glx_context(self, share):
        if share:
            return glx.glXCreateNewContext(self._display, self._fbconfig,
                glx.GLX_RGBA_TYPE, share._context, True)
        else:
            return glx.glXCreateNewContext(self._display, self._fbconfig,
                glx.GLX_RGBA_TYPE, None, True)

class XlibGLContext(gl.Context):
    def __init__(self, config, context, share):
        super(XlibGLContext, self).__init__(share)
        self.config = config
        self._context = context
        self._x_display = config.screen.display._display

    def destroy(self):
        super(XlibGLContext, self).destroy()
        glx.glXDestroyContext(self._x_display, self._context)

    def is_direct(self):
        return glx.glXIsDirect(self._x_display, self._context)

# Platform event data is single item, so use platform event handler directly.
XlibEventHandler = _PlatformEventHandler

class XlibWindow(BaseWindow):
    _x_display = None       # X display connection
    _x_screen_id = None     # X screen index
    _x_ic = None            # X input context
    _glx_context = None     # GLX context handle
    _glx_window = None      # GLX window handle
    _window = None          # Xlib window handle
    _minimum_size = None
    _maximum_size = None

    _x = 0
    _y = 0                  # Last known window position
    _width = 0
    _height = 0             # Last known window size
    _mouse_exclusive_client = None  # x,y of "real" mouse during exclusive
    _mouse_buttons = [False] * 6 # State of each xlib button
    _keyboard_exclusive = False
    _active = True
    _applied_mouse_exclusive = False
    _applied_keyboard_exclusive = False
    _mapped = False
    _lost_context = False
    _lost_context_state = False

    _enable_xsync = False
    _current_sync_value = None
    _current_sync_valid = False

    _needs_resize = False   # True when resize event has been received but not
                            # dispatched

    _default_event_mask = (0x1ffffff 
        & ~xlib.PointerMotionHintMask
        & ~xlib.ResizeRedirectMask)

    def __init__(self, *args, **kwargs):
        # Bind event handlers
        self._event_handlers = {}
        for name in self._platform_event_names:
            if not hasattr(self, name):
                continue
            func = getattr(self, name)
            for message in func._platform_event_data:
                self._event_handlers[message] = func 

        super(XlibWindow, self).__init__(*args, **kwargs)
        
    def _recreate(self, changes):
        # If flipping to/from fullscreen and using override_redirect (we
        # always are, _NET_WM_FULLSCREEN doesn't work), need to recreate the
        # window.
        #
        # A possible improvement could be to just hide the top window,
        # destroy the GLX window, and reshow it again when leaving fullscreen.
        # This would prevent the floating window from being moved by the
        # WM.
        if 'fullscreen' in changes or 'resizable' in changes:
            # clear out the GLX context
            self.switch_to()
            gl.glFlush()
            glx.glXMakeCurrent(self._x_display, 0, None)
            if self._glx_window:
                glx.glXDestroyWindow(self._x_display, self._glx_window)
            xlib.XDestroyWindow(self._x_display, self._window)
            self._glx_window = None
            del self.display._window_map[self._window]
            self._window = None 
            self._mapped = False

        # TODO: detect state loss only by examining context share.
        if 'context' in changes:
            self._lost_context = True
            self._lost_context_state = True

        self._create()

    def _create(self):
        # Unmap existing window if necessary while we fiddle with it.
        if self._window and self._mapped:
            self._unmap()

        self.context.window = self
        self._x_display = self.config._display
        self._x_screen_id = self.screen._x_screen_id
        self._glx_context = self.context._context

        self._glx_1_3 = self.display.info.have_version(1, 3)
        self._have_SGI_video_sync = \
            self.display.info.have_extension('GLX_SGI_video_sync')
        self._have_SGI_swap_control = \
            self.display.info.have_extension('GLX_SGI_swap_control')
        self._have_MESA_swap_control = \
            self.display.info.have_extension('GLX_MESA_swap_control')

        # In order of preference:
        # 1. GLX_MESA_swap_control (more likely to work where video_sync will
        #    not)
        # 2. GLX_SGI_video_sync (does not work on Intel 945GM, but that has
        #    MESA)
        # 3. GLX_SGI_swap_control (cannot be disabled once enabled).
        self._use_video_sync = (self._have_SGI_video_sync and 
                                not self._have_MESA_swap_control)

        # Create X window if not already existing.
        if not self._window:
            root = xlib.XRootWindow(self._x_display, self._x_screen_id)

            visual_info = self.config.get_visual_info()

            visual = visual_info.visual
            visual_id = xlib.XVisualIDFromVisual(visual)
            default_visual = xlib.XDefaultVisual(
                self._x_display, self._x_screen_id)
            default_visual_id = xlib.XVisualIDFromVisual(default_visual)
            window_attributes = xlib.XSetWindowAttributes()
            if visual_id != default_visual_id:
                window_attributes.colormap = xlib.XCreateColormap(
                    self._x_display, root, visual, xlib.AllocNone)
            else:
                window_attributes.colormap = xlib.XDefaultColormap(
                    self._x_display, self._x_screen_id)
            window_attributes.bit_gravity = xlib.NorthWestGravity

            # Issue 287: Compiz on Intel/Mesa doesn't draw window decoration
            #            unless CWBackPixel is given in mask.  Should have
            #            no effect on other systems, so it's set
            #            unconditionally.
            mask = xlib.CWColormap | xlib.CWBitGravity | xlib.CWBackPixel

            self._window = xlib.XCreateWindow(self._x_display, root,
                0, 0, self._width, self._height, 0, visual_info.depth,
                xlib.InputOutput, visual, mask,
                byref(window_attributes))
            self.display._window_map[self._window] = self

            # Setting null background pixmap disables drawing the background,
            # preventing flicker while resizing (in theory).
            #
            # Issue 287: Compiz on Intel/Mesa doesn't draw window decoration if
            #            this is called.  As it doesn't seem to have any
            #            effect anyway, it's just commented out.
            #xlib.XSetWindowBackgroundPixmap(self._x_display, self._window, 0)

            self._enable_xsync = (pyglet.options['xsync'] and
                                  self.display._enable_xsync and
                                  self.config.double_buffer)

            # Set supported protocols
            protocols = []
            protocols.append(xlib.XInternAtom(self._x_display,
                                              'WM_DELETE_WINDOW', False))
            if self._enable_xsync:
                protocols.append(xlib.XInternAtom(self._x_display,
                                                  '_NET_WM_SYNC_REQUEST',
                                                  False))
            protocols = (c_ulong * len(protocols))(*protocols)
            xlib.XSetWMProtocols(self._x_display, self._window,
                                 protocols, len(protocols))

            # Create window resize sync counter
            if self._enable_xsync:
                value = xsync.XSyncValue()
                self._sync_counter = xlib.XID(
                    xsync.XSyncCreateCounter(self._x_display, value))
                atom = xlib.XInternAtom(self._x_display,
                                        '_NET_WM_SYNC_REQUEST_COUNTER', False)
                ptr = pointer(self._sync_counter)

                xlib.XChangeProperty(self._x_display, self._window,
                                     atom, XA_CARDINAL, 32,
                                     xlib.PropModeReplace,
                                     cast(ptr, POINTER(c_ubyte)), 1)
        # Set window attributes
        attributes = xlib.XSetWindowAttributes()
        attributes_mask = 0

        # Bypass the window manager in fullscreen.  This is the only reliable
        # technique (over _NET_WM_STATE_FULLSCREEN, Motif, KDE and Gnome
        # hints) that is pretty much guaranteed to work.  Unfortunately
        # we run into window activation and focus problems that require
        # attention.  Search for "override_redirect" for all occurences.
        attributes.override_redirect = self._fullscreen
        attributes_mask |= xlib.CWOverrideRedirect

        if self._fullscreen:
            xlib.XMoveResizeWindow(self._x_display, self._window, 
                self.screen.x, self.screen.y, 
                self.screen.width, self.screen.height)
        else:
            xlib.XResizeWindow(self._x_display, self._window, 
                self._width, self._height)

        xlib.XChangeWindowAttributes(self._x_display, self._window, 
            attributes_mask, byref(attributes))

        # Set style
        styles = {
            self.WINDOW_STYLE_DEFAULT: '_NET_WM_WINDOW_TYPE_NORMAL',
            self.WINDOW_STYLE_DIALOG: '_NET_WM_WINDOW_TYPE_DIALOG',
            self.WINDOW_STYLE_TOOL: '_NET_WM_WINDOW_TYPE_UTILITY',
        }
        if self._style in styles:
            self._set_atoms_property('_NET_WM_WINDOW_TYPE', 
                                     (styles[self._style],))
        elif self._style == self.WINDOW_STYLE_BORDERLESS:
            MWM_HINTS_DECORATIONS = 1 << 1
            PROP_MWM_HINTS_ELEMENTS = 5
            mwmhints = mwmhints_t()
            mwmhints.flags = MWM_HINTS_DECORATIONS
            mwmhints.decorations = 0
            name = xlib.XInternAtom(self._x_display, '_MOTIF_WM_HINTS', False)
            xlib.XChangeProperty(self._x_display, self._window,
                name, name, 32, xlib.PropModeReplace, 
                cast(pointer(mwmhints), POINTER(c_ubyte)),
                PROP_MWM_HINTS_ELEMENTS)

        # Set resizeable
        if not self._resizable:
            self.set_minimum_size(self._width, self._height)
            self.set_maximum_size(self._width, self._height)

        # Set caption
        self.set_caption(self._caption)

        # Create input context.  A good but very outdated reference for this
        # is http://www.sbin.org/doc/Xlib/chapt_11.html
        if _have_utf8 and not self._x_ic:
            if not self.display._x_im:
                xlib.XSetLocaleModifiers('@im=none')
                self.display._x_im = \
                    xlib.XOpenIM(self._x_display, None, None, None)

            xlib.XFlush(self._x_display);

            # Need to set argtypes on this function because it's vararg,
            # and ctypes guesses wrong.
            xlib.XCreateIC.argtypes = [xlib.XIM,    
                                       c_char_p, c_int,
                                       c_char_p, xlib.Window,
                                       c_char_p, xlib.Window,
                                       c_void_p]
            self._x_ic = xlib.XCreateIC(self.display._x_im, 
                'inputStyle', xlib.XIMPreeditNothing|xlib.XIMStatusNothing,
                'clientWindow', self._window,
                'focusWindow', self._window,
                None)

            filter_events = c_ulong()
            xlib.XGetICValues(self._x_ic,
                              'filterEvents', byref(filter_events),
                              None)
            self._default_event_mask |= filter_events.value
            xlib.XSetICFocus(self._x_ic)

        self.switch_to()
        if self._visible:
            self.set_visible(True)

        self.set_mouse_platform_visible()

    def _map(self):
        if self._mapped:
            return

        # Map the window, wait for map event before continuing.
        xlib.XSelectInput(
            self._x_display, self._window, xlib.StructureNotifyMask)
        xlib.XMapRaised(self._x_display, self._window)
        e = xlib.XEvent()
        while True:
            xlib.XNextEvent(self._x_display, e)
            if e.type == xlib.MapNotify:
                break
        xlib.XSelectInput(
            self._x_display, self._window, self._default_event_mask)
        self._mapped = True

        if self._fullscreen:
            # Possibly an override_redirect issue.
            self.activate()

        self.dispatch_event('on_resize', self._width, self._height)
        self.dispatch_event('on_show')
        self.dispatch_event('on_expose')

    def _unmap(self):
        if not self._mapped:
            return

        xlib.XSelectInput(
            self._x_display, self._window, xlib.StructureNotifyMask)
        xlib.XUnmapWindow(self._x_display, self._window)
        e = xlib.XEvent()
        while True:
            xlib.XNextEvent(self._x_display, e)
            if e.type == xlib.UnmapNotify:
                break

        xlib.XSelectInput(
            self._x_display, self._window, self._default_event_mask)
        self._mapped = False

    def _get_root(self):
        attributes = xlib.XWindowAttributes()
        xlib.XGetWindowAttributes(self._x_display, self._window,
                                  byref(attributes))
        return attributes.root

    def close(self):
        if not self._window:
            return

        # clear out the GLX context.  Can fail if current context already
        # destroyed (on exit, say).
        try:
            gl.glFlush()
        except gl.GLException:
            pass

        if self._glx_1_3:
            glx.glXMakeContextCurrent(self._x_display, 0, 0, None)
        else:
            glx.glXMakeCurrent(self._x_display, 0, None)

        self._unmap()
        if self._glx_window:
            glx.glXDestroyWindow(self._x_display, self._glx_window)
        if self._window:
            xlib.XDestroyWindow(self._x_display, self._window)

        del self.display._window_map[self._window]
        self._window = None
        self._glx_window = None

        if _have_utf8:
            xlib.XDestroyIC(self._x_ic)
            self._x_ic = None

        super(XlibWindow, self).close()

    def switch_to(self):
        if self._glx_1_3:
            if not self._glx_window:
                self._glx_window = glx.glXCreateWindow(self._x_display,
                    self._config._fbconfig, self._window, None)
            glx.glXMakeContextCurrent(self._x_display,
                self._glx_window, self._glx_window, self._glx_context)
        else:
            glx.glXMakeCurrent(self._x_display, self._window, self._glx_context)

        self.set_vsync(self._vsync)

        self._context.set_current()
        gl_info.set_active_context()
        glu_info.set_active_context()

    def flip(self):
        self.draw_mouse_cursor()

        if self._vsync and self._have_SGI_video_sync and self._use_video_sync:
            count = c_uint()
            glxext_arb.glXGetVideoSyncSGI(byref(count))
            glxext_arb.glXWaitVideoSyncSGI(
                2, (count.value + 1) % 2, byref(count))

        if self._glx_1_3:
            if not self._glx_window:
                self._glx_window = glx.glXCreateWindow(self._x_display,
                    self._config._fbconfig, self._window, None)
            glx.glXSwapBuffers(self._x_display, self._glx_window)
        else:
            glx.glXSwapBuffers(self._x_display, self._window)

        self._sync_resize()

    def set_vsync(self, vsync):
        if pyglet.options['vsync'] is not None:
            vsync = pyglet.options['vsync']
        self._vsync = vsync
        if not self._use_video_sync:
            interval = vsync and 1 or 0
            if self._have_MESA_swap_control:
                glxext_mesa.glXSwapIntervalMESA(interval)
            elif self._have_SGI_swap_control and interval:
                # SGI_swap_control interval cannot be set to 0
                glxext_arb.glXSwapIntervalSGI(interval)

    def set_caption(self, caption):
        if caption is None:
            caption = ''
        self._caption = caption
        self._set_text_property('WM_NAME', caption, allow_utf8=False)
        self._set_text_property('WM_ICON_NAME', caption, allow_utf8=False)
        self._set_text_property('_NET_WM_NAME', caption)
        self._set_text_property('_NET_WM_ICON_NAME', caption)

    def get_caption(self):
        return self._caption

    def set_size(self, width, height):
        if self._fullscreen:
            raise WindowException('Cannot set size of fullscreen window.')
        self._width = width
        self._height = height
        if not self._resizable:
            self.set_minimum_size(width, height)
            self.set_maximum_size(width, height)
        xlib.XResizeWindow(self._x_display, self._window, width, height)
        self.dispatch_event('on_resize', width, height)

    def get_size(self):
        # XGetGeometry and XWindowAttributes seem to always return the
        # original size of the window, which is wrong after the user
        # has resized it.
        # XXX this is probably fixed now, with fix of resize.
        return self._width, self._height

    def set_location(self, x, y):
        # Assume the window manager has reparented our top-level window
        # only once, in which case attributes.x/y give the offset from
        # the frame to the content window.  Better solution would be
        # to use _NET_FRAME_EXTENTS, where supported.
        attributes = xlib.XWindowAttributes()
        xlib.XGetWindowAttributes(self._x_display, self._window,
                                  byref(attributes))
        # XXX at least under KDE's WM these attrs are both 0
        x -= attributes.x
        y -= attributes.y
        xlib.XMoveWindow(self._x_display, self._window, x, y)

    def get_location(self):
        child = xlib.Window()
        x = c_int()
        y = c_int()
        xlib.XTranslateCoordinates(self._x_display,
                                   self._window,
                                   self._get_root(),
                                   0, 0,
                                   byref(x),
                                   byref(y),
                                   byref(child))
        return x.value, y.value

    def activate(self):
        xlib.XSetInputFocus(self._x_display, self._window,
            xlib.RevertToParent, xlib.CurrentTime)

    def set_visible(self, visible=True):
        if visible:
            self._map()
        else:
            self._unmap()
        self._visible = visible

    def set_minimum_size(self, width, height):
        self._minimum_size = width, height
        self._set_wm_normal_hints()

    def set_maximum_size(self, width, height):
        self._maximum_size = width, height
        self._set_wm_normal_hints()

    def minimize(self):
        xlib.XIconifyWindow(self._x_display, self._window, self._x_screen_id)

    def maximize(self):
        self._set_wm_state('_NET_WM_STATE_MAXIMIZED_HORZ',
                           '_NET_WM_STATE_MAXIMIZED_VERT')

    def set_mouse_platform_visible(self, platform_visible=None):
        if platform_visible is None:
            platform_visible = self._mouse_visible and \
                               not self._mouse_cursor.drawable

        if not platform_visible:
            # Hide pointer by creating an empty cursor
            black = xlib.XBlackPixel(self._x_display, self._x_screen_id)
            black = xlib.XColor()
            bmp = xlib.XCreateBitmapFromData(self._x_display, self._window,
                c_buffer(8), 8, 8)
            cursor = xlib.XCreatePixmapCursor(self._x_display, bmp, bmp,
                black, black, 0, 0)
            xlib.XDefineCursor(self._x_display, self._window, cursor)
            xlib.XFreeCursor(self._x_display, cursor)
            xlib.XFreePixmap(self._x_display, bmp)
        else:
            # Restore cursor
            if isinstance(self._mouse_cursor, XlibMouseCursor):
                xlib.XDefineCursor(self._x_display, self._window, 
                                   self._mouse_cursor.cursor)
            else:
                xlib.XUndefineCursor(self._x_display, self._window)

    def _update_exclusivity(self):
        mouse_exclusive = self._active and self._mouse_exclusive
        keyboard_exclusive = self._active and self._keyboard_exclusive

        if mouse_exclusive != self._applied_mouse_exclusive:
            if mouse_exclusive:
                self.set_mouse_platform_visible(False)

                # Restrict to client area
                xlib.XGrabPointer(self._x_display, self._window,
                    True,
                    0,
                    xlib.GrabModeAsync,
                    xlib.GrabModeAsync,
                    self._window,
                    0,
                    xlib.CurrentTime)

                # Move pointer to center of window
                x = self._width / 2
                y = self._height / 2
                self._mouse_exclusive_client = x, y
                xlib.XWarpPointer(self._x_display,
                    0,              # src window
                    self._window,   # dst window
                    0, 0,           # src x, y
                    0, 0,           # src w, h
                    x, y)
            else:
                # Unclip
                xlib.XUngrabPointer(self._x_display, xlib.CurrentTime)
                self.set_mouse_platform_visible()

            self._applied_mouse_exclusive = mouse_exclusive

        if keyboard_exclusive != self._applied_keyboard_exclusive:
            if keyboard_exclusive:
                xlib.XGrabKeyboard(self._x_display,
                    self._window,
                    False,
                    xlib.GrabModeAsync,
                    xlib.GrabModeAsync,
                    xlib.CurrentTime)
            else:
                 xlib.XUngrabKeyboard(self._x_display, xlib.CurrentTime)
            self._applied_keyboard_exclusive = keyboard_exclusive

    def set_exclusive_mouse(self, exclusive=True):
        if exclusive == self._mouse_exclusive:
            return

        self._mouse_exclusive = exclusive
        self._update_exclusivity()

    def set_exclusive_keyboard(self, exclusive=True):
        if exclusive == self._keyboard_exclusive:
            return
        
        self._keyboard_exclusive = exclusive
        self._update_exclusivity()

    def get_system_mouse_cursor(self, name):
        if name == self.CURSOR_DEFAULT:
            return DefaultMouseCursor()

        # NQR means default shape is not pretty... surely there is another
        # cursor font?
        cursor_shapes = {
            self.CURSOR_CROSSHAIR:       cursorfont.XC_crosshair,
            self.CURSOR_HAND:            cursorfont.XC_hand2,
            self.CURSOR_HELP:            cursorfont.XC_question_arrow,  # NQR
            self.CURSOR_NO:              cursorfont.XC_pirate,          # NQR
            self.CURSOR_SIZE:            cursorfont.XC_fleur,
            self.CURSOR_SIZE_UP:         cursorfont.XC_top_side,
            self.CURSOR_SIZE_UP_RIGHT:   cursorfont.XC_top_right_corner,
            self.CURSOR_SIZE_RIGHT:      cursorfont.XC_right_side,
            self.CURSOR_SIZE_DOWN_RIGHT: cursorfont.XC_bottom_right_corner,
            self.CURSOR_SIZE_DOWN:       cursorfont.XC_bottom_side,
            self.CURSOR_SIZE_DOWN_LEFT:  cursorfont.XC_bottom_left_corner,
            self.CURSOR_SIZE_LEFT:       cursorfont.XC_left_side,
            self.CURSOR_SIZE_UP_LEFT:    cursorfont.XC_top_left_corner,
            self.CURSOR_SIZE_UP_DOWN:    cursorfont.XC_sb_v_double_arrow,
            self.CURSOR_SIZE_LEFT_RIGHT: cursorfont.XC_sb_h_double_arrow,
            self.CURSOR_TEXT:            cursorfont.XC_xterm,
            self.CURSOR_WAIT:            cursorfont.XC_watch,
            self.CURSOR_WAIT_ARROW:      cursorfont.XC_watch,           # NQR
        }
        if name not in cursor_shapes:
            raise MouseCursorException('Unknown cursor name "%s"' % name)
        cursor = xlib.XCreateFontCursor(self._x_display, cursor_shapes[name])
        return XlibMouseCursor(cursor)

    def set_icon(self, *images):
        # Careful!  XChangeProperty takes an array of long when data type
        # is 32-bit (but long can be 64 bit!), so pad high bytes of format if
        # necessary.

        import sys
        format = {
            ('little', 4): 'BGRA',
            ('little', 8): 'BGRAAAAA',
            ('big', 4):    'ARGB',
            ('big', 8):    'AAAAARGB'
        }[(sys.byteorder, sizeof(c_ulong))]

        data = ''
        for image in images:
            image = image.get_image_data()
            pitch = -(image.width * len(format))
            s = c_buffer(sizeof(c_ulong) * 2)
            memmove(s, cast((c_ulong * 2)(image.width, image.height), 
                            POINTER(c_ubyte)), len(s))
            data += s.raw + image.get_data(format, pitch)
        buffer = (c_ubyte * len(data))()
        memmove(buffer, data, len(data))
        atom = xlib.XInternAtom(self._x_display, '_NET_WM_ICON', False)
        xlib.XChangeProperty(self._x_display, self._window, atom, XA_CARDINAL,
            32, xlib.PropModeReplace, buffer, len(data)/sizeof(c_ulong))

    # Private utility

    def _set_wm_normal_hints(self):
        hints = xlib.XAllocSizeHints().contents
        if self._minimum_size:
            hints.flags |= xlib.PMinSize
            hints.min_width, hints.min_height = self._minimum_size
        if self._maximum_size:
            hints.flags |= xlib.PMaxSize
            hints.max_width, hints.max_height = self._maximum_size
        xlib.XSetWMNormalHints(self._x_display, self._window, byref(hints))

    def _set_text_property(self, name, value, allow_utf8=True):
        atom = xlib.XInternAtom(self._x_display, name, False)
        if not atom:
            raise XlibException('Undefined atom "%s"' % name)
        assert type(value) in (str, unicode)
        property = xlib.XTextProperty()
        if _have_utf8 and allow_utf8:
            buf = create_string_buffer(value.encode('utf8'))
            result = xlib.Xutf8TextListToTextProperty(self._x_display,
                cast(pointer(buf), c_char_p), 1, xlib.XUTF8StringStyle, 
                byref(property))
            if result < 0:
                raise XlibException('Could not create UTF8 text property')
        else:
            buf = create_string_buffer(value.encode('ascii', 'ignore'))
            result = xlib.XStringListToTextProperty(
                cast(pointer(buf), c_char_p), 1, byref(property))
            if result < 0:
                raise XlibException('Could not create text property')
        xlib.XSetTextProperty(self._x_display,
            self._window, byref(property), atom)
        # XXX <rj> Xlib doesn't like us freeing this
        #xlib.XFree(property.value)

    def _set_atoms_property(self, name, values, mode=xlib.PropModeReplace):
        name_atom = xlib.XInternAtom(self._x_display, name, False)
        atoms = []
        for value in values:
            atoms.append(xlib.XInternAtom(self._x_display, value, False))
        atom_type = xlib.XInternAtom(self._x_display, 'ATOM', False)
        if len(atoms):
            atoms_ar = (xlib.Atom * len(atoms))(*atoms)
            xlib.XChangeProperty(self._x_display, self._window,
                name_atom, atom_type, 32, mode,
                cast(pointer(atoms_ar), POINTER(c_ubyte)), len(atoms))
        else:
            xlib.XDeleteProperty(self._x_display, self._window, net_wm_state)

    def _set_wm_state(self, *states):
        # Set property
        net_wm_state = xlib.XInternAtom(self._x_display, '_NET_WM_STATE', False)
        atoms = []
        for state in states:
            atoms.append(xlib.XInternAtom(self._x_display, state, False))
        atom_type = xlib.XInternAtom(self._x_display, 'ATOM', False)
        if len(atoms):
            atoms_ar = (xlib.Atom * len(atoms))(*atoms)
            xlib.XChangeProperty(self._x_display, self._window,
                net_wm_state, atom_type, 32, xlib.PropModePrepend,
                cast(pointer(atoms_ar), POINTER(c_ubyte)), len(atoms))
        else:
            xlib.XDeleteProperty(self._x_display, self._window, net_wm_state)

        # Nudge the WM
        e = xlib.XEvent()
        e.xclient.type = xlib.ClientMessage
        e.xclient.message_type = net_wm_state
        e.xclient.display = cast(self._x_display, POINTER(xlib.Display))
        e.xclient.window = self._window
        e.xclient.format = 32
        e.xclient.data.l[0] = xlib.PropModePrepend
        for i, atom in enumerate(atoms):
            e.xclient.data.l[i + 1] = atom
        xlib.XSendEvent(self._x_display, self._get_root(),
            False, xlib.SubstructureRedirectMask, byref(e))

    # Event handling

    def dispatch_events(self):
        self.dispatch_pending_events()

        self._allow_dispatch_event = True

        e = xlib.XEvent()

        # Cache these in case window is closed from an event handler
        _x_display = self._x_display
        _window = self._window

        # Check for the events specific to this window
        while xlib.XCheckWindowEvent(_x_display, _window,
                                     0x1ffffff, byref(e)):
            # Key events are filtered by the xlib window event
            # handler so they get a shot at the prefiltered event.
            if e.xany.type not in (xlib.KeyPress, xlib.KeyRelease):
                if xlib.XFilterEvent(e, 0):
                    continue
            self.dispatch_platform_event(e)

        # Generic events for this window (the window close event).
        while xlib.XCheckTypedWindowEvent(_x_display, _window, 
                xlib.ClientMessage, byref(e)):
            self.dispatch_platform_event(e) 
        
        if self._needs_resize:
            self.dispatch_event('on_resize', self._width, self._height)
            self.dispatch_event('on_expose')
            self._needs_resize = False

        self._allow_dispatch_event = False

    def dispatch_pending_events(self):
        while self._event_queue:
            EventDispatcher.dispatch_event(self, *self._event_queue.pop(0))

        # Dispatch any context-related events
        if self._lost_context:
            self._lost_context = False
            EventDispatcher.dispatch_event(self, 'on_context_lost')
        if self._lost_context_state:
            self._lost_context_state = False
            EventDispatcher.dispatch_event(self, 'on_context_state_lost')

    def dispatch_platform_event(self, e):
        event_handler = self._event_handlers.get(e.type)
        if event_handler:
            event_handler(e)

    @staticmethod
    def _translate_modifiers(state):
        modifiers = 0
        if state & xlib.ShiftMask:
            modifiers |= key.MOD_SHIFT
        if state & xlib.ControlMask:
            modifiers |= key.MOD_CTRL
        if state & xlib.LockMask:
            modifiers |= key.MOD_CAPSLOCK
        if state & xlib.Mod1Mask:
            modifiers |= key.MOD_ALT
        if state & xlib.Mod2Mask:
            modifiers |= key.MOD_NUMLOCK
        if state & xlib.Mod4Mask:
            modifiers |= key.MOD_WINDOWS
        if state & xlib.Mod5Mask:
            modifiers |= key.MOD_SCROLLLOCK
        return modifiers

    # Event handlers
    '''
    def _event_symbol(self, event):
        # pyglet.self.key keysymbols are identical to X11 keysymbols, no
        # need to map the keysymbol.
        symbol = xlib.XKeycodeToKeysym(self._x_display, event.xkey.keycode, 0)
        if symbol == 0:
            # XIM event
            return None
        elif symbol not in key._key_names.keys():
            symbol = key.user_key(event.xkey.keycode)
        return symbol
    '''

    def _event_text_symbol(self, ev):
        text = None
        symbol = xlib.KeySym()
        buffer = create_string_buffer(128)

        # Look up raw keysym before XIM filters it (default for keypress and
        # keyrelease)
        count = xlib.XLookupString(ev.xkey,
                                   buffer, len(buffer) - 1,
                                   byref(symbol), None)

        # Give XIM a shot
        filtered = xlib.XFilterEvent(ev, ev.xany.window)

        if ev.type == xlib.KeyPress and not filtered:
            status = c_int()
            if _have_utf8:
                encoding = 'utf8'
                count = xlib.Xutf8LookupString(self._x_ic,
                                               ev.xkey,
                                               buffer, len(buffer) - 1,
                                               byref(symbol), byref(status))
                if status.value == xlib.XBufferOverflow:
                    raise NotImplementedError('TODO: XIM buffer resize')

            else:
                encoding = 'ascii'
                count = xlib.XLookupString(ev.xkey,
                                           buffer, len(buffer) - 1,
                                           byref(symbol), None)
                if count:
                    status.value = xlib.XLookupBoth

            if status.value & (xlib.XLookupChars | xlib.XLookupBoth):
                text = buffer.value[:count].decode(encoding)

            # Don't treat Unicode command codepoints as text, except Return.
            if text and unicodedata.category(text) == 'Cc' and text != '\r':
                text = None

        symbol = symbol.value

        # If the event is a XIM filtered event, the keysym will be virtual
        # (e.g., aacute instead of A after a dead key).  Drop it, we don't
        # want these kind of key events.
        if ev.xkey.keycode == 0 and not filtered:
            symbol = None

        # pyglet.self.key keysymbols are identical to X11 keysymbols, no
        # need to map the keysymbol.  For keysyms outside the pyglet set, map
        # raw key code to a user key.
        if symbol and symbol not in key._key_names and ev.xkey.keycode:
            symbol = key.user_key(ev.xkey.keycode)

        if filtered:
            # The event was filtered, text must be ignored, but the symbol is
            # still good.
            return None, symbol

        return text, symbol

    def _event_text_motion(self, symbol, modifiers):
        if modifiers & key.MOD_ALT:
            return None
        ctrl = modifiers & key.MOD_CTRL != 0
        return _motion_map.get((symbol, ctrl), None)

    @XlibEventHandler(xlib.KeyPress)
    @XlibEventHandler(xlib.KeyRelease)
    def _event_key(self, ev):
        if ev.type == xlib.KeyRelease:
            # Look in the queue for a matching KeyPress with same timestamp,
            # indicating an auto-repeat rather than actual key event.
            saved = []
            while True:
                auto_event = xlib.XEvent()
                result = xlib.XCheckWindowEvent(self._x_display,
                    self._window, xlib.KeyPress|xlib.KeyRelease, 
                    byref(auto_event))
                if not result:
                    break
                saved.append(auto_event)
                if auto_event.type == xlib.KeyRelease:
                    # just save this off for restoration back to the queue
                    continue
                if ev.xkey.keycode == auto_event.xkey.keycode:
                    # Found a key repeat: dispatch EVENT_TEXT* event
                    text, symbol = self._event_text_symbol(auto_event)
                    modifiers = self._translate_modifiers(ev.xkey.state)
                    modifiers_ctrl = modifiers & (key.MOD_CTRL | key.MOD_ALT)
                    motion = self._event_text_motion(symbol, modifiers)
                    if motion:
                        if modifiers & key.MOD_SHIFT:
                            self.dispatch_event(
                                'on_text_motion_select', motion)
                        else:
                            self.dispatch_event('on_text_motion', motion)
                    elif text and not modifiers_ctrl:
                        self.dispatch_event('on_text', text)

                    ditched = saved.pop()
                    for auto_event in reversed(saved):
                        xlib.XPutBackEvent(self._x_display, byref(auto_event))
                    return
                else:
                    # Key code of press did not match, therefore no repeating
                    # is going on, stop searching.
                    break
            # Whoops, put the events back, it's for real.
            for auto_event in reversed(saved):
                xlib.XPutBackEvent(self._x_display, byref(auto_event))

        text, symbol = self._event_text_symbol(ev)
        modifiers = self._translate_modifiers(ev.xkey.state)
        modifiers_ctrl = modifiers & (key.MOD_CTRL | key.MOD_ALT)
        motion = self._event_text_motion(symbol, modifiers)

        if ev.type == xlib.KeyPress:
            if symbol:
                self.dispatch_event('on_key_press', symbol, modifiers)
            if motion:
                if modifiers & key.MOD_SHIFT:
                    self.dispatch_event('on_text_motion_select', motion)
                else:
                    self.dispatch_event('on_text_motion', motion)
            elif text and not modifiers_ctrl:
                self.dispatch_event('on_text', text)
        elif ev.type == xlib.KeyRelease:
            if symbol:
                self.dispatch_event('on_key_release', symbol, modifiers)

    @XlibEventHandler(xlib.MotionNotify)
    def _event_motionnotify(self, ev):
        x = ev.xmotion.x
        y = self.height - ev.xmotion.y

        dx = x - self._mouse_x
        dy = y - self._mouse_y

        if self._applied_mouse_exclusive and \
           (ev.xmotion.x, ev.xmotion.y) == self._mouse_exclusive_client:
            # Ignore events caused by XWarpPointer
            self._mouse_x = x
            self._mouse_y = y
            return

        if self._applied_mouse_exclusive:
            # Reset pointer position
            ex, ey = self._mouse_exclusive_client
            xlib.XWarpPointer(self._x_display,
                0,
                self._window,
                0, 0,
                0, 0,
                ex, ey)

        self._mouse_x = x
        self._mouse_y = y
        self._mouse_in_window = True

        buttons = 0
        if ev.xmotion.state & xlib.Button1MotionMask:
            buttons |= mouse.LEFT
        if ev.xmotion.state & xlib.Button2MotionMask:
            buttons |= mouse.MIDDLE
        if ev.xmotion.state & xlib.Button3MotionMask:
            buttons |= mouse.RIGHT

        if buttons:
            # Drag event
            modifiers = self._translate_modifiers(ev.xmotion.state)
            self.dispatch_event('on_mouse_drag',
                x, y, dx, dy, buttons, modifiers)
        else:
            # Motion event
            self.dispatch_event('on_mouse_motion', x, y, dx, dy)

    @XlibEventHandler(xlib.ClientMessage)
    def _event_clientmessage(self, ev):
        atom = ev.xclient.data.l[0]
        if atom == xlib.XInternAtom(ev.xclient.display,
                                    'WM_DELETE_WINDOW', False):
            self.dispatch_event('on_close')
        elif (self._enable_xsync and 
              atom == xlib.XInternAtom(ev.xclient.display, 
                                       '_NET_WM_SYNC_REQUEST', False)):
            lo = ev.xclient.data.l[2]
            hi = ev.xclient.data.l[3]
            self._current_sync_value = xsync.XSyncValue(hi, lo)

    def _sync_resize(self):
        if self._enable_xsync and self._current_sync_valid:
            if xsync.XSyncValueIsZero(self._current_sync_value):
                self._current_sync_valid = False
                return
            xsync.XSyncSetCounter(self._x_display, 
                                  self._sync_counter,
                                  self._current_sync_value)
            self._current_sync_value = None
            self._current_sync_valid = False

    @XlibEventHandler(xlib.ButtonPress)
    @XlibEventHandler(xlib.ButtonRelease)
    def _event_button(self, ev):
        x = ev.xbutton.x
        y = self.height - ev.xbutton.y
        button = 1 << (ev.xbutton.button - 1)  # 1, 2, 3 -> 1, 2, 4
        modifiers = self._translate_modifiers(ev.xbutton.state)
        if ev.type == xlib.ButtonPress:
            # override_redirect issue: manually activate this window if
            # fullscreen.
            if self._fullscreen and not self._active:
                self.activate()

            if ev.xbutton.button == 4:
                self.dispatch_event('on_mouse_scroll', x, y, 0, 1)
            elif ev.xbutton.button == 5:
                self.dispatch_event('on_mouse_scroll', x, y, 0, -1)
            elif ev.xbutton.button < len(self._mouse_buttons):
                self._mouse_buttons[ev.xbutton.button] = True
                self.dispatch_event('on_mouse_press', 
                    x, y, button, modifiers)
        else:
            if ev.xbutton.button < 4:
                self._mouse_buttons[ev.xbutton.button] = False
                self.dispatch_event('on_mouse_release', 
                    x, y, button, modifiers)

    @XlibEventHandler(xlib.Expose)
    def _event_expose(self, ev):
        # Ignore all expose events except the last one. We could be told
        # about exposure rects - but I don't see the point since we're
        # working with OpenGL and we'll just redraw the whole scene.
        if ev.xexpose.count > 0: return
        self.dispatch_event('on_expose')

    @XlibEventHandler(xlib.EnterNotify)
    def _event_enternotify(self, ev):
        # figure active mouse buttons
        # XXX ignore modifier state?
        state = ev.xcrossing.state
        self._mouse_buttons[1] = state & xlib.Button1Mask
        self._mouse_buttons[2] = state & xlib.Button2Mask
        self._mouse_buttons[3] = state & xlib.Button3Mask
        self._mouse_buttons[4] = state & xlib.Button4Mask
        self._mouse_buttons[5] = state & xlib.Button5Mask

        # mouse position
        x = self._mouse_x = ev.xcrossing.x
        y = self._mouse_y = self.height - ev.xcrossing.y
        self._mouse_in_window = True

        # XXX there may be more we could do here
        self.dispatch_event('on_mouse_enter', x, y)

    @XlibEventHandler(xlib.LeaveNotify)
    def _event_leavenotify(self, ev):
        x = self._mouse_x = ev.xcrossing.x
        y = self._mouse_y = self.height - ev.xcrossing.y
        self._mouse_in_window = False
        self.dispatch_event('on_mouse_leave', x, y)

    @XlibEventHandler(xlib.ConfigureNotify)
    def _event_configurenotify(self, ev):
        if self._enable_xsync and self._current_sync_value:
            self._current_sync_valid = True

        self.switch_to()

        w, h = ev.xconfigure.width, ev.xconfigure.height
        x, y = ev.xconfigure.x, ev.xconfigure.y
        if self._width != w or self._height != h:
            self._width = w
            self._height = h
            self._needs_resize = True
        if self._x != x or self._y != y:
            self.dispatch_event('on_move', x, y)
            self._x = x
            self._y = y

    @XlibEventHandler(xlib.FocusIn)
    def _event_focusin(self, ev):
        self._active = True
        self._update_exclusivity()
        self.dispatch_event('on_activate')
        xlib.XSetICFocus(self._x_ic)

    @XlibEventHandler(xlib.FocusOut)
    def _event_focusout(self, ev):
        self._active = False
        self._update_exclusivity()
        self.dispatch_event('on_deactivate')
        xlib.XUnsetICFocus(self._x_ic)

    @XlibEventHandler(xlib.MapNotify)
    def _event_mapnotify(self, ev):
        self._mapped = True
        self.dispatch_event('on_show')

    @XlibEventHandler(xlib.UnmapNotify)
    def _event_unmapnotify(self, ev):
        self._mapped = False
        self.dispatch_event('on_hide')
