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
__version__ = '$Id: $'

from ctypes import *
import os.path
import unicodedata
import warnings

import pyglet
from pyglet.window import WindowException, Platform, Display, Screen, \
    BaseWindow, MouseCursor, DefaultMouseCursor, _PlatformEventHandler
from pyglet.window import key
from pyglet.window import mouse
from pyglet.window import event
from pyglet.window.carbon.constants import *
from pyglet.window.carbon.types import *
from pyglet.window.carbon.quartzkey import keymap

import pyglet.lib
from pyglet import gl
from pyglet.gl import agl
from pyglet.gl import gl_info
from pyglet.gl import glu_info
from pyglet.event import EventDispatcher

class CarbonException(WindowException):
    pass

carbon = pyglet.lib.load_library(
    framework='/System/Library/Frameworks/Carbon.framework')
quicktime = pyglet.lib.load_library(
    framework='/System/Library/Frameworks/QuickTime.framework')

carbon.GetEventDispatcherTarget.restype = EventTargetRef
carbon.ReceiveNextEvent.argtypes = \
    [c_uint32, c_void_p, c_double, c_ubyte, POINTER(EventRef)]
carbon.GetWindowPort.restype = agl.AGLDrawable
EventHandlerProcPtr = CFUNCTYPE(c_int, c_int, c_void_p, c_void_p)
carbon.NewEventHandlerUPP.restype = c_void_p
carbon.GetCurrentKeyModifiers = c_uint32
carbon.NewRgn.restype = RgnHandle
carbon.CGDisplayBounds.argtypes = [c_void_p]
carbon.CGDisplayBounds.restype = CGRect

# Map symbol,modifiers -> motion
# Determined by experiment with TextEdit.app
_motion_map = {
    (key.UP, False):                    key.MOTION_UP,
    (key.RIGHT, False):                 key.MOTION_RIGHT,
    (key.DOWN, False):                  key.MOTION_DOWN,
    (key.LEFT, False):                  key.MOTION_LEFT,
    (key.LEFT, key.MOD_OPTION):         key.MOTION_PREVIOUS_WORD,
    (key.RIGHT, key.MOD_OPTION):        key.MOTION_NEXT_WORD,
    (key.LEFT, key.MOD_COMMAND):        key.MOTION_BEGINNING_OF_LINE,
    (key.RIGHT, key.MOD_COMMAND):       key.MOTION_END_OF_LINE,
    (key.PAGEUP, False):                key.MOTION_PREVIOUS_PAGE,
    (key.PAGEDOWN, False):              key.MOTION_NEXT_PAGE,
    (key.HOME, False):                  key.MOTION_BEGINNING_OF_FILE,
    (key.END, False):                   key.MOTION_END_OF_FILE,
    (key.UP, key.MOD_COMMAND):          key.MOTION_BEGINNING_OF_FILE,
    (key.DOWN, key.MOD_COMMAND):        key.MOTION_END_OF_FILE,
    (key.BACKSPACE, False):             key.MOTION_BACKSPACE,
    (key.DELETE, False):                key.MOTION_DELETE,
}

class CarbonPlatform(Platform):
    _display = None

    def get_default_display(self):
        if not self._display:
            self._display = CarbonDisplay()
        return self._display

class CarbonDisplay(Display):
    # TODO: CarbonDisplay could be per display device, which would make
    # reporting of screens and available configs more accurate.  The number of
    # Macs with more than one video card is probably small, though.
    def __init__(self):
        super(CarbonDisplay, self).__init__()

        import MacOS
        if not MacOS.WMAvailable():
            raise CarbonException('Window manager is not available.  ' \
                                  'Ensure you run "pythonw", not "python"')

        self._install_application_event_handlers()
        
    def get_screens(self):
        count = CGDisplayCount()
        carbon.CGGetActiveDisplayList(0, None, byref(count))
        displays = (CGDirectDisplayID * count.value)()
        carbon.CGGetActiveDisplayList(count.value, displays, byref(count))
        return [CarbonScreen(self, id) for id in displays]

    def _install_application_event_handlers(self):
        self._carbon_event_handlers = []
        self._carbon_event_handler_refs = []

        target = carbon.GetApplicationEventTarget()

        # TODO something with a metaclass or hacky like CarbonWindow
        # to make this list extensible
        handlers = [
            (self._on_mouse_down, kEventClassMouse, kEventMouseDown),
            (self._on_apple_event, kEventClassAppleEvent, kEventAppleEvent),
            (self._on_command, kEventClassCommand, kEventProcessCommand),
        ]

        ae_handlers = [
            (self._on_ae_quit, kCoreEventClass, kAEQuitApplication),
        ]

        # Install the application-wide handlers
        for method, cls, event in handlers:
            proc = EventHandlerProcPtr(method)
            self._carbon_event_handlers.append(proc)
            upp = carbon.NewEventHandlerUPP(proc)
            types = EventTypeSpec()
            types.eventClass = cls
            types.eventKind = event
            handler_ref = EventHandlerRef()
            carbon.InstallEventHandler(
                target,
                upp,
                1,
                byref(types),
                c_void_p(),
                byref(handler_ref))
            self._carbon_event_handler_refs.append(handler_ref)

        # Install Apple event handlers
        for method, cls, event in ae_handlers:
            proc = EventHandlerProcPtr(method)
            self._carbon_event_handlers.append(proc)
            upp = carbon.NewAEEventHandlerUPP(proc)
            carbon.AEInstallEventHandler(
                cls,
                event,
                upp,
                0,
                False)

    def _on_command(self, next_handler, ev, data):
        command = HICommand()
        carbon.GetEventParameter(ev, kEventParamDirectObject,
            typeHICommand, c_void_p(), sizeof(command), c_void_p(),
            byref(command))

        if command.commandID == kHICommandQuit:
            self._on_quit()

        return noErr

    def _on_mouse_down(self, next_handler, ev, data):
        # Check for menubar hit
        position = Point()
        carbon.GetEventParameter(ev, kEventParamMouseLocation,
            typeQDPoint, c_void_p(), sizeof(position), c_void_p(),
            byref(position))
        if carbon.FindWindow(position, None) == inMenuBar:
            # Mouse down in menu bar.  MenuSelect() takes care of all
            # menu tracking and blocks until the menu is dismissed.
            # Use command events to handle actual menu item invokations.
            carbon.MenuSelect(position)
            carbon.HiliteMenu(0)

        carbon.CallNextEventHandler(next_handler, ev)
        return noErr

    def _on_apple_event(self, next_handler, ev, data):
        # Somewhat involved way of redispatching Apple event contained
        # within a Carbon event, described in
        # http://developer.apple.com/documentation/AppleScript/
        #  Conceptual/AppleEvents/dispatch_aes_aepg/chapter_4_section_3.html

        release = False
        if carbon.IsEventInQueue(carbon.GetMainEventQueue(), ev):
            carbon.RetainEvent(ev)
            release = True
            carbon.RemoveEventFromQueue(carbon.GetMainEventQueue(), ev)

        ev_record = EventRecord()
        carbon.ConvertEventRefToEventRecord(ev, byref(ev_record))
        carbon.AEProcessAppleEvent(byref(ev_record))

        if release:
            carbon.ReleaseEvent(ev)
        
        return noErr

    def _on_ae_quit(self, ae, reply, refcon):
        self._on_quit()
        return noErr

    def _on_quit(self):
        '''Called when the user tries to quit the application.

        This is not an actual event handler, it is called in response
        to Command+Q, the Quit menu item, and the Dock context menu's Quit
        item.

        The default implementation sets `has_exit` to true on all open
        windows.
        '''
        for window in self.get_windows():
            window.has_exit = True
    
class CarbonScreen(Screen):
    def __init__(self, display, id):
        self.display = display
        rect = carbon.CGDisplayBounds(id)
        super(CarbonScreen, self).__init__(
            int(rect.origin.x), int(rect.origin.y),
            int(rect.size.width), int(rect.size.height))
        self.id = id

    def get_gdevice(self):
        gdevice = GDHandle()
        r = carbon.DMGetGDeviceByDisplayID(self.id, byref(gdevice), False)
        _oscheck(r)
        return gdevice

    def get_matching_configs(self, template):
        # Construct array of attributes for aglChoosePixelFormat
        attrs = []
        for name, value in template.get_gl_attributes():
            attr = CarbonGLConfig._attribute_ids.get(name, None)
            if not attr or not value:
                continue
            attrs.append(attr)
            if attr not in CarbonGLConfig._boolean_attributes:
                attrs.append(int(value))

        # Support for RAGE-II, which is not compliant
        attrs.append(agl.AGL_ALL_RENDERERS)

        # Force selection policy and RGBA
        attrs.append(agl.AGL_MAXIMUM_POLICY)
        attrs.append(agl.AGL_RGBA)

        # In 10.3 and later, AGL_FULLSCREEN is specified so the window can
        # be toggled to/from fullscreen without losing context.  pyglet
        # no longer supports earlier versions of OS X, so we always supply it.
        attrs.append(agl.AGL_FULLSCREEN)

        # Terminate the list.
        attrs.append(agl.AGL_NONE)
        attrib_list = (c_int * len(attrs))(*attrs)

        device = self.get_gdevice()
        pformat = agl.aglChoosePixelFormat(device, 1, attrib_list)
        _aglcheck()

        if not pformat:
            return []
        else:
            return [CarbonGLConfig(self, pformat)]

class CarbonGLConfig(gl.Config):
    # Valid names for GL attributes, and their corresponding AGL constant. 
    _attribute_ids = {
        'double_buffer': agl.AGL_DOUBLEBUFFER,
        'stereo': agl.AGL_STEREO,
        'buffer_size': agl.AGL_BUFFER_SIZE, 
        'sample_buffers': agl.AGL_SAMPLE_BUFFERS_ARB,
        'samples': agl.AGL_SAMPLES_ARB,
        'aux_buffers': agl.AGL_AUX_BUFFERS,
        'red_size': agl.AGL_RED_SIZE,
        'green_size': agl.AGL_GREEN_SIZE,
        'blue_size': agl.AGL_BLUE_SIZE,
        'alpha_size': agl.AGL_ALPHA_SIZE,
        'depth_size': agl.AGL_DEPTH_SIZE,
        'stencil_size': agl.AGL_STENCIL_SIZE,
        'accum_red_size': agl.AGL_ACCUM_RED_SIZE,
        'accum_green_size': agl.AGL_ACCUM_GREEN_SIZE,
        'accum_blue_size': agl.AGL_ACCUM_BLUE_SIZE,
        'accum_alpha_size': agl.AGL_ACCUM_ALPHA_SIZE,

        # Not exposed by pyglet API (set internally)
        'all_renderers': agl.AGL_ALL_RENDERERS,
        'rgba': agl.AGL_RGBA,
        'fullscreen': agl.AGL_FULLSCREEN,
        'minimum_policy': agl.AGL_MINIMUM_POLICY,
        'maximum_policy': agl.AGL_MAXIMUM_POLICY,

        # Not supported in current pyglet API
        'level': agl.AGL_LEVEL, 
        'pixel_size': agl.AGL_PIXEL_SIZE,   # == buffer_size
        'aux_depth_stencil': agl.AGL_AUX_DEPTH_STENCIL,
        'color_float': agl.AGL_COLOR_FLOAT,
        'offscreen': agl.AGL_OFFSCREEN,
        'sample_alpha': agl.AGL_SAMPLE_ALPHA,
        'multisample': agl.AGL_MULTISAMPLE,
        'supersample': agl.AGL_SUPERSAMPLE,
    }

    # AGL constants which do not require a value.
    _boolean_attributes = \
        (agl.AGL_ALL_RENDERERS, 
         agl.AGL_RGBA,
         agl.AGL_DOUBLEBUFFER,
         agl.AGL_STEREO,
         agl.AGL_MINIMUM_POLICY,
         agl.AGL_MAXIMUM_POLICY,
         agl.AGL_OFFSCREEN,
         agl.AGL_FULLSCREEN,
         agl.AGL_AUX_DEPTH_STENCIL,
         agl.AGL_COLOR_FLOAT,
         agl.AGL_MULTISAMPLE,
         agl.AGL_SUPERSAMPLE,
         agl.AGL_SAMPLE_ALPHA)

    def __init__(self, screen, pformat):
        super(CarbonGLConfig, self).__init__()
        self.screen = screen
        self._pformat = pformat
        self._attributes = {}

        for name, attr in self._attribute_ids.items():
            value = c_int()
            result = agl.aglDescribePixelFormat(pformat, attr, byref(value))
            if result:
                setattr(self, name, value.value)
 
    def create_context(self, share):
        if share:
            context = agl.aglCreateContext(self._pformat, share._context)
        else:
            context = agl.aglCreateContext(self._pformat, None)
        _aglcheck()
        return CarbonGLContext(self, context, share, self._pformat)

class CarbonGLContext(gl.Context):
    def __init__(self, config, context, share, pixelformat):
        super(CarbonGLContext, self).__init__(share)
        self.config = config
        self._context = context
        self._pixelformat = pixelformat

    def destroy(self):
        super(CarbonGLContext, self).destroy()
        agl.aglDestroyContext(self._context)

class CarbonMouseCursor(MouseCursor):
    drawable = False
    def __init__(self, theme):
        self.theme = theme

def CarbonEventHandler(event_class, event_kind):
    return _PlatformEventHandler((event_class, event_kind))

class CarbonWindow(BaseWindow):
    _window = None                  # Carbon WindowRef
    _agl_context = None             # AGL context ID
    _recreate_deferred = None

    # Window properties
    _minimum_size = None
    _maximum_size = None
    _fullscreen_restore = None
    _event_dispatcher = None
    _current_modifiers = 0
    _mapped_modifers = 0
    _carbon_event_handlers = []
    _carbon_event_handler_refs = []
    _track_ref = 0
    _track_region = None

    _mouse_exclusive = False
    _mouse_platform_visible = True

    def _recreate(self, changes):
        # We can't destroy the window while event handlers are active,
        # otherwise the (OS X) event dispatcher gets lost and segfaults.
        #
        # Defer actual recreation until dispatch_events next finishes.
        self._recreate_deferred = changes

    def _recreate_immediate(self):
        # The actual _recreate function.
        changes = self._recreate_deferred
        self._recreate_deferred = None

        if ('context' in changes):
            agl.aglSetDrawable(self._agl_context, None)

        if ('fullscreen' in changes and
            not self._fullscreen and
            self._fullscreen_restore):

            # Leaving fullscreen -- destroy everything before the window.
            self._remove_track_region()
            self._remove_event_handlers()
            agl.aglSetDrawable(self._agl_context, None)
            # EndFullScreen disposes _window.
            quicktime.EndFullScreen(self._fullscreen_restore, 0)
            self._window = None

        self._create()

    def _create(self):
        self._agl_context = self.context._context

        if self._window:
            # The window is about to be recreated; destroy everything
            # associated with the old window, then the window itself.
            self._remove_track_region()
            self._remove_event_handlers()
            agl.aglSetDrawable(self._agl_context, None)
            carbon.DisposeWindow(self._window)
            self._window = None

        self._window = WindowRef()

        if self._fullscreen:
            # Switch to fullscreen mode with QuickTime
            fs_width = c_short(0)
            fs_height = c_short(0)
            self._fullscreen_restore = c_void_p()
            quicktime.BeginFullScreen(byref(self._fullscreen_restore),
                                      self.screen.get_gdevice(),
                                      byref(fs_width),
                                      byref(fs_height),
                                      byref(self._window),
                                      None,
                                      0)
            # the following may be used for debugging if you have a second
            # monitor - only the main monitor will go fullscreen
            agl.aglEnable(self._agl_context, agl.AGL_FS_CAPTURE_SINGLE)
            self._width = fs_width.value
            self._height = fs_height.value
            #self._width = self.screen.width
            #self._height = self.screen.height
            agl.aglSetFullScreen(self._agl_context, 
                                 self._width, self._height, 0, 0)
            self._mouse_in_window = True
            self.dispatch_event('on_resize', self._width, self._height)
            self.dispatch_event('on_show')
            self.dispatch_event('on_expose')
        else:
            # Create floating window
            rect = Rect()
            location = None # TODO
            if location is not None:
                rect.left = location[0]
                rect.top = location[1]
            else:
                rect.top = rect.left = 0
            rect.right = rect.left + self._width
            rect.bottom = rect.top + self._height

            styles = {
                self.WINDOW_STYLE_DEFAULT:  (kDocumentWindowClass,
                                             kWindowCloseBoxAttribute |
                                             kWindowCollapseBoxAttribute),
                self.WINDOW_STYLE_DIALOG:   (kDocumentWindowClass,
                                             kWindowCloseBoxAttribute),
                self.WINDOW_STYLE_TOOL:     (kUtilityWindowClass,
                                             kWindowCloseBoxAttribute),
                self.WINDOW_STYLE_BORDERLESS:    (kSimpleWindowClass,
                                                  kWindowNoAttributes)
            }
            window_class, window_attributes = \
                styles.get(self._style, kDocumentWindowClass)

            if self._resizable:
                window_attributes |= (kWindowFullZoomAttribute |
                                      kWindowResizableAttribute)

            r = carbon.CreateNewWindow(window_class,
                                       window_attributes,
                                       byref(rect),
                                       byref(self._window))
            _oscheck(r)

            if location is None:
                carbon.RepositionWindow(self._window, c_void_p(),
                    kWindowCascadeOnMainScreen)

            agl.aglSetDrawable(self._agl_context,
                carbon.GetWindowPort(self._window))

        _aglcheck()

        self.set_caption(self._caption)

        # Get initial state
        self._event_dispatcher = carbon.GetEventDispatcherTarget()
        self._current_modifiers = carbon.GetCurrentKeyModifiers().value
        self._mapped_modifiers = self._map_modifiers(self._current_modifiers)

        # (re)install Carbon event handlers 
        self._install_event_handlers()

        self._create_track_region()

        self.switch_to()
        self.set_vsync(self._vsync)

        if self._visible:
            self.set_visible(True)

    def _create_track_region(self):
        self._remove_track_region()

        # Create a tracking region for the content part of the window
        # to receive enter/leave events.
        track_id = MouseTrackingRegionID()
        track_id.signature = DEFAULT_CREATOR_CODE
        track_id.id = 1
        self._track_ref = MouseTrackingRef()
        self._track_region = carbon.NewRgn()
        carbon.GetWindowRegion(self._window, 
            kWindowContentRgn, self._track_region)
        carbon.CreateMouseTrackingRegion(self._window,  
            self._track_region, None, kMouseTrackingOptionsGlobalClip,
            track_id, None, None,
            byref(self._track_ref))

    def _remove_track_region(self):
        if self._track_region:
            carbon.ReleaseMouseTrackingRegion(self._track_region)
            self._track_region = None

    def close(self):
        super(CarbonWindow, self).close()
        self._agl_context = None
        self._remove_event_handlers()
        self._remove_track_region()

        # Restore cursor visibility
        self.set_mouse_platform_visible(True)
        self.set_exclusive_mouse(False)

        if self._fullscreen:
            quicktime.EndFullScreen(self._fullscreen_restore, 0)
        else:
            carbon.DisposeWindow(self._window)
        self._window = None

    def switch_to(self):
        agl.aglSetCurrentContext(self._agl_context)
        self._context.set_current()
        _aglcheck()
        gl_info.set_active_context()
        glu_info.set_active_context()

    def flip(self):
        self.draw_mouse_cursor()
        agl.aglSwapBuffers(self._agl_context)
        _aglcheck()

    def _get_vsync(self):
        swap = c_long()
        agl.aglGetInteger(self._agl_context, agl.AGL_SWAP_INTERVAL, byref(swap))
        return bool(swap.value)
    vsync = property(_get_vsync) # overrides BaseWindow property

    def set_vsync(self, vsync):
        if pyglet.options['vsync'] is not None:
            vsync = pyglet.options['vsync']
        self._vsync = vsync # _recreate depends on this
        swap = c_long(int(vsync))
        agl.aglSetInteger(self._agl_context, agl.AGL_SWAP_INTERVAL, byref(swap))

    def dispatch_events(self):
        if self._recreate_deferred:
            self._recreate_immediate()

        self._allow_dispatch_event = True
        while self._event_queue:
            EventDispatcher.dispatch_event(self, *self._event_queue.pop(0))

        e = EventRef()
        result = carbon.ReceiveNextEvent(0, c_void_p(), 0, True, byref(e))
        while result == noErr:
            carbon.SendEventToEventTarget(e, self._event_dispatcher)
            carbon.ReleaseEvent(e)

            if self._recreate_deferred:
                self._recreate_immediate()
            result = carbon.ReceiveNextEvent(0, c_void_p(), 0, True, byref(e))

        self._allow_dispatch_event = False

        # Return value from ReceiveNextEvent can be ignored if not
        # noErr; we check here only to look for new bugs.
        # eventLoopQuitErr: the inner event loop was quit, see
        # http://lists.apple.com/archives/Carbon-dev/2006/Jun/msg00850.html
        # Can occur when mixing with other toolkits, e.g. Tk.
        # Fixes issue 180.
        if result not in (eventLoopTimedOutErr, eventLoopQuitErr):
            raise 'Error %d' % result

    def set_caption(self, caption):
        self._caption = caption
        s = _create_cfstring(caption)
        carbon.SetWindowTitleWithCFString(self._window, s)
        carbon.CFRelease(s)

    def set_location(self, x, y):
        rect = Rect()
        carbon.GetWindowBounds(self._window, kWindowContentRgn, byref(rect))
        rect.right += x - rect.left
        rect.bottom += y - rect.top
        rect.left = x
        rect.top = y
        carbon.SetWindowBounds(self._window, kWindowContentRgn, byref(rect))

    def get_location(self):
        rect = Rect()
        carbon.GetWindowBounds(self._window, kWindowContentRgn, byref(rect))
        return rect.left, rect.top

    def set_size(self, width, height):
        if self._fullscreen:
            raise WindowException('Cannot set size of fullscreen window.')
        rect = Rect()
        carbon.GetWindowBounds(self._window, kWindowContentRgn, byref(rect))
        rect.right = rect.left + width
        rect.bottom = rect.top + height
        carbon.SetWindowBounds(self._window, kWindowContentRgn, byref(rect))

        self._width = width
        self._height = height
        self.dispatch_event('on_resize', width, height)
        self.dispatch_event('on_expose')

    def get_size(self):
        if self._fullscreen:
            return self._width, self._height
        rect = Rect()
        carbon.GetWindowBounds(self._window, kWindowContentRgn, byref(rect))
        return rect.right - rect.left, rect.bottom - rect.top

    def set_minimum_size(self, width, height):
        self._minimum_size = (width, height)
        minimum = HISize()
        minimum.width = width
        minimum.height = height
        if self._maximum_size:
            maximum = HISize()
            maximum.width, maximum.height = self._maximum_size
            maximum = byref(maximum)
        else:
            maximum = None
        carbon.SetWindowResizeLimits(self._window, 
            byref(minimum), maximum)

    def set_maximum_size(self, width, height):
        self._maximum_size = (width, height)
        maximum = HISize()
        maximum.width = width
        maximum.height = height
        if self._minimum_size:
            minimum = HISize()
            minimum.width, minimum.height = self._minimum_size
            minimum = byref(minimum)
        else:
            minimum = None
        carbon.SetWindowResizeLimits(self._window, 
            minimum, byref(maximum))

    def activate(self):
        carbon.ActivateWindow(self._window, 1)

        # Also make the application the "front" application.  TODO
        # maybe don't bring forward all of the application's windows?
        psn = ProcessSerialNumber()
        psn.highLongOfPSN = 0
        psn.lowLongOfPSN = kCurrentProcess
        carbon.SetFrontProcess(byref(psn))

    def set_visible(self, visible=True):
        self._visible = visible
        if visible:
            self.dispatch_event('on_resize', self._width, self._height)
            self.dispatch_event('on_show')
            carbon.ShowWindow(self._window)
        else:
            carbon.HideWindow(self._window)

    def minimize(self):
        carbon.CollapseWindow(self._window, True)

    def maximize(self):
        # Maximum "safe" value, gets trimmed to screen size automatically.
        p = Point()
        p.v, p.h = 16000,16000 
        if not carbon.IsWindowInStandardState(self._window, byref(p), None):
            carbon.ZoomWindowIdeal(self._window, inZoomOut, byref(p))

    def set_mouse_platform_visible(self, platform_visible=None):
        if platform_visible is None:
            platform_visible = self._mouse_visible and \
                               not self._mouse_exclusive and \
                               not self._mouse_cursor.drawable
        if not self._mouse_in_window:
            platform_visible = True

        if self._mouse_in_window and \
           isinstance(self._mouse_cursor, CarbonMouseCursor):
            carbon.SetThemeCursor(self._mouse_cursor.theme)
        else:
            carbon.SetThemeCursor(kThemeArrowCursor)

        if self._mouse_platform_visible == platform_visible:
            return

        if platform_visible:
            carbon.ShowCursor()
        else:
            carbon.HideCursor()
        self._mouse_platform_visible = platform_visible

    def set_exclusive_mouse(self, exclusive=True):
        self._mouse_exclusive = exclusive
        if exclusive:
            # Move mouse to center of window
            rect = Rect()
            carbon.GetWindowBounds(self._window, kWindowContentRgn, byref(rect))
            point = CGPoint()
            point.x = (rect.right + rect.left) / 2
            point.y = (rect.bottom + rect.top) / 2
            carbon.CGWarpMouseCursorPosition(point)
            carbon.CGAssociateMouseAndMouseCursorPosition(False)
        else:
            carbon.CGAssociateMouseAndMouseCursorPosition(True)
        self.set_mouse_platform_visible()

    def set_exclusive_keyboard(self, exclusive=True):
        if exclusive:
            # Note: power switch can also be disabled, with
            # kUIOptionDisableSessionTerminate.  That seems
            # a little extreme though.
            carbon.SetSystemUIMode(kUIModeAllHidden,
                (kUIOptionDisableAppleMenu |
                 kUIOptionDisableProcessSwitch |
                 kUIOptionDisableForceQuit |
                 kUIOptionDisableHide))
        else:
            carbon.SetSystemUIMode(kUIModeNormal, 0)

    def get_system_mouse_cursor(self, name):
        if name == self.CURSOR_DEFAULT:
            return DefaultMouseCursor()

        themes = {
            self.CURSOR_CROSSHAIR:       kThemeCrossCursor,
            self.CURSOR_HAND:            kThemePointingHandCursor,
            self.CURSOR_HELP:            kThemeArrowCursor,
            self.CURSOR_NO:              kThemeNotAllowedCursor,
            self.CURSOR_SIZE:            kThemeArrowCursor,
            self.CURSOR_SIZE_UP:         kThemeResizeUpCursor,
            self.CURSOR_SIZE_UP_RIGHT:   kThemeArrowCursor,
            self.CURSOR_SIZE_RIGHT:      kThemeResizeRightCursor,
            self.CURSOR_SIZE_DOWN_RIGHT: kThemeArrowCursor,
            self.CURSOR_SIZE_DOWN:       kThemeResizeDownCursor,
            self.CURSOR_SIZE_DOWN_LEFT:  kThemeArrowCursor,
            self.CURSOR_SIZE_LEFT:       kThemeResizeLeftCursor,
            self.CURSOR_SIZE_UP_LEFT:    kThemeArrowCursor,
            self.CURSOR_SIZE_UP_DOWN:    kThemeResizeUpDownCursor,
            self.CURSOR_SIZE_LEFT_RIGHT: kThemeResizeLeftRightCursor,
            self.CURSOR_TEXT:            kThemeIBeamCursor,
            self.CURSOR_WAIT:            kThemeWatchCursor,
            self.CURSOR_WAIT_ARROW:      kThemeWatchCursor,
        }
        if name not in themes:
            raise CarbonException('Unknown cursor name "%s"' % name)
        return CarbonMouseCursor(themes[name])

    def set_icon(self, *images):
        # Only use the biggest image
        image = images[0]
        size = image.width * image.height
        for img in images:
            if img.width * img.height > size:
                size = img.width * img.height
                image = img

        image = image.image_data
        image.format = 'ARGB'
        image.pitch = -len(image.format) * image.width

        provider = carbon.CGDataProviderCreateWithData(
            None, image.data, len(image.data), None)

        colorspace = carbon.CGColorSpaceCreateDeviceRGB()

        cgi = carbon.CGImageCreate(
            image.width, image.height, 8, 32, -image.pitch,
            colorspace,
            kCGImageAlphaFirst,
            provider,
            None,
            True,
            kCGRenderingIntentDefault)

        carbon.SetApplicationDockTileImage(cgi)

        carbon.CGDataProviderRelease(provider)
        carbon.CGColorSpaceRelease(colorspace)

    # Non-public utilities

    def _update_drawable(self):
        # We can get there after context has been disposed, in which case
        # just do nothing.
        if not self._agl_context:
            return

        agl.aglUpdateContext(self._agl_context)
        _aglcheck()

        # Need a redraw
        self.dispatch_event('on_expose')

    def _update_track_region(self):
        carbon.GetWindowRegion(self._window, 
            kWindowContentRgn, self._track_region)
        carbon.ChangeMouseTrackingRegion(self._track_ref,
            self._track_region, None)

    def _install_event_handlers(self):
        self._remove_event_handlers()

        if self._fullscreen:
            target = carbon.GetApplicationEventTarget()
        else:
            target = carbon.GetWindowEventTarget(self._window)
        carbon.InstallStandardEventHandler(target)

        self._carbon_event_handlers = []
        self._carbon_event_handler_refs = []

        for func_name in self._platform_event_names:
            if not hasattr(self, func_name):
                continue

            func = getattr(self, func_name)
            for event_class, event_kind in func._platform_event_data:
                # TODO: could just build up array of class/kind
                proc = EventHandlerProcPtr(func)
                self._carbon_event_handlers.append(proc)
                upp = carbon.NewEventHandlerUPP(proc)
                types = EventTypeSpec()
                types.eventClass = event_class
                types.eventKind = event_kind
                handler_ref = EventHandlerRef()
                carbon.InstallEventHandler(
                    target,
                    upp,
                    1,
                    byref(types),
                    c_void_p(),
                    byref(handler_ref))
                self._carbon_event_handler_refs.append(handler_ref)

    def _remove_event_handlers(self):
        for ref in self._carbon_event_handler_refs:
            carbon.RemoveEventHandler(ref)
        self._carbon_event_handler_refs = []
        self._carbon_event_handlers = []

    # Carbon event handlers

    @CarbonEventHandler(kEventClassTextInput, kEventTextInputUnicodeForKeyEvent)
    def _on_text_input(self, next_handler, ev, data):
        size = c_uint32()
        carbon.GetEventParameter(ev, kEventParamTextInputSendText,
            typeUTF8Text, c_void_p(), 0, byref(size), c_void_p())
        text = create_string_buffer(size.value)
        carbon.GetEventParameter(ev, kEventParamTextInputSendText,
            typeUTF8Text, c_void_p(), size.value, c_void_p(), byref(text))
        text = text.value.decode('utf8')

        raw_event = EventRef()
        carbon.GetEventParameter(ev, kEventParamTextInputSendKeyboardEvent,
            typeEventRef, c_void_p(), sizeof(raw_event), c_void_p(),
            byref(raw_event))
        symbol, modifiers = self._get_symbol_and_modifiers(raw_event)

        motion_modifiers = modifiers & \
            (key.MOD_COMMAND | key.MOD_CTRL | key.MOD_OPTION)
        if (symbol, motion_modifiers) in _motion_map:
            motion = _motion_map[symbol, motion_modifiers]
            if modifiers & key.MOD_SHIFT:
                self.dispatch_event('on_text_motion_select', motion)
            else:
                self.dispatch_event('on_text_motion', motion)
        elif ((unicodedata.category(text[0]) != 'Cc' or text == u'\r') and
            not (modifiers & key.MOD_COMMAND)):
            self.dispatch_event('on_text', text)
        return noErr

    @CarbonEventHandler(kEventClassKeyboard, kEventRawKeyUp)
    def _on_key_up(self, next_handler, ev, data):
        symbol, modifiers = self._get_symbol_and_modifiers(ev)
        if symbol:
            self.dispatch_event('on_key_release', symbol, modifiers)
        carbon.CallNextEventHandler(next_handler, ev)
        return noErr

    @CarbonEventHandler(kEventClassKeyboard, kEventRawKeyDown)
    def _on_key_down(self, next_handler, ev, data):
        symbol, modifiers = self._get_symbol_and_modifiers(ev)
        if symbol:
            self.dispatch_event('on_key_press', symbol, modifiers)
        carbon.CallNextEventHandler(next_handler, ev)
        return noErr

    @staticmethod
    def _get_symbol_and_modifiers(ev):
        sym = c_uint32()
        carbon.GetEventParameter(ev, kEventParamKeyCode,
            typeUInt32, c_void_p(), sizeof(sym), c_void_p(), byref(sym))
        modifiers = c_uint32()
        carbon.GetEventParameter(ev, kEventParamKeyModifiers,
            typeUInt32, c_void_p(), sizeof(modifiers), c_void_p(),
            byref(modifiers))

        symbol = keymap.get(sym.value, None)
        if symbol is None:
            symbol = key.user_key(sym.value)

        return (symbol, CarbonWindow._map_modifiers(modifiers.value))

    @staticmethod
    def _map_modifiers(modifiers):
        mapped_modifiers = 0
        if modifiers & (shiftKey | rightShiftKey):
            mapped_modifiers |= key.MOD_SHIFT
        if modifiers & (controlKey | rightControlKey):
            mapped_modifiers |= key.MOD_CTRL
        if modifiers & (optionKey | rightOptionKey):
            mapped_modifiers |= key.MOD_OPTION
        if modifiers & alphaLock:
            mapped_modifiers |= key.MOD_CAPSLOCK
        if modifiers & cmdKey:
            mapped_modifiers |= key.MOD_COMMAND

        return mapped_modifiers

    @CarbonEventHandler(kEventClassKeyboard, kEventRawKeyModifiersChanged)
    def _on_modifiers_changed(self, next_handler, ev, data):
        modifiers = c_uint32()
        carbon.GetEventParameter(ev, kEventParamKeyModifiers,
            typeUInt32, c_void_p(), sizeof(modifiers), c_void_p(),
            byref(modifiers))
        modifiers = modifiers.value
        deltas = modifiers ^ self._current_modifiers
        for mask, k in [
            (controlKey, key.LCTRL),
            (shiftKey, key.LSHIFT),
            (cmdKey, key.LCOMMAND),
            (optionKey, key.LOPTION),
            (rightShiftKey, key.RSHIFT),
            (rightOptionKey, key.ROPTION),
            (rightControlKey, key.RCTRL),
            (alphaLock, key.CAPSLOCK),
            (numLock, key.NUMLOCK)]:
            if deltas & mask:
                if modifiers & mask:
                    self.dispatch_event('on_key_press', 
                        k, self._mapped_modifiers)
                else:
                    self.dispatch_event('on_key_release',
                        k, self._mapped_modifiers)
        carbon.CallNextEventHandler(next_handler, ev)

        self._mapped_modifiers = self._map_modifiers(modifiers)
        self._current_modifiers = modifiers
        return noErr

    def _get_mouse_position(self, ev):
        position = HIPoint()
        carbon.GetEventParameter(ev, kEventParamMouseLocation,
            typeHIPoint, c_void_p(), sizeof(position), c_void_p(),
            byref(position))

        bounds = Rect()
        carbon.GetWindowBounds(self._window, kWindowContentRgn, byref(bounds))
        return int(position.x - bounds.left), int(position.y - bounds.top)

    @staticmethod
    def _get_mouse_button_and_modifiers(ev):
        button = EventMouseButton()
        carbon.GetEventParameter(ev, kEventParamMouseButton,
            typeMouseButton, c_void_p(), sizeof(button), c_void_p(),
            byref(button))
        
        if button.value == 1: 
            button = mouse.LEFT
        elif button.value == 2: 
            button = mouse.RIGHT
        elif button.value == 3: 
            button = mouse.MIDDLE

        modifiers = c_uint32()
        carbon.GetEventParameter(ev, kEventParamKeyModifiers,
            typeUInt32, c_void_p(), sizeof(modifiers), c_void_p(),
            byref(modifiers))

        return button, CarbonWindow._map_modifiers(modifiers.value)

    @staticmethod
    def _get_mouse_in_content(ev):
        position = Point()
        carbon.GetEventParameter(ev, kEventParamMouseLocation,
            typeQDPoint, c_void_p(), sizeof(position), c_void_p(),
            byref(position)) 
        return carbon.FindWindow(position, None) == inContent 

    @CarbonEventHandler(kEventClassMouse, kEventMouseDown)
    def _on_mouse_down(self, next_handler, ev, data):
        if self._fullscreen or self._get_mouse_in_content(ev):
            button, modifiers = self._get_mouse_button_and_modifiers(ev)
            x, y = self._get_mouse_position(ev)
            y = self.height - y
            self.dispatch_event('on_mouse_press', x, y, button, modifiers)

        carbon.CallNextEventHandler(next_handler, ev)
        return noErr

    @CarbonEventHandler(kEventClassMouse, kEventMouseUp)
    def _on_mouse_up(self, next_handler, ev, data):
        # Always report mouse up, even out of content area, because it's
        # probably after a drag gesture.
        button, modifiers = self._get_mouse_button_and_modifiers(ev)
        x, y = self._get_mouse_position(ev)
        y = self.height - y
        self.dispatch_event('on_mouse_release', x, y, button, modifiers)

        carbon.CallNextEventHandler(next_handler, ev)
        return noErr

    @CarbonEventHandler(kEventClassMouse, kEventMouseMoved)
    def _on_mouse_moved(self, next_handler, ev, data):
        if self._fullscreen or self._get_mouse_in_content(ev):
            x, y = self._get_mouse_position(ev)
            y = self.height - y

            self._mouse_x = x
            self._mouse_y = y

            delta = HIPoint()
            carbon.GetEventParameter(ev, kEventParamMouseDelta,
                typeHIPoint, c_void_p(), sizeof(delta), c_void_p(),
                byref(delta))

            # Motion event
            self.dispatch_event('on_mouse_motion', 
                x, y, delta.x, -delta.y)

        carbon.CallNextEventHandler(next_handler, ev)
        return noErr

    @CarbonEventHandler(kEventClassMouse, kEventMouseDragged)
    def _on_mouse_dragged(self, next_handler, ev, data):
        button, modifiers = self._get_mouse_button_and_modifiers(ev)
        x, y = self._get_mouse_position(ev)
        y = self.height - y

        self._mouse_x = x
        self._mouse_y = y

        delta = HIPoint()
        carbon.GetEventParameter(ev, kEventParamMouseDelta,
            typeHIPoint, c_void_p(), sizeof(delta), c_void_p(),
            byref(delta))

        # Drag event
        self.dispatch_event('on_mouse_drag',
            x, y, delta.x, -delta.y, button, modifiers)

        carbon.CallNextEventHandler(next_handler, ev)
        return noErr

    @CarbonEventHandler(kEventClassMouse, kEventMouseEntered)
    def _on_mouse_entered(self, next_handler, ev, data):
        x, y = self._get_mouse_position(ev)
        y = self.height - y

        self._mouse_x = x
        self._mouse_y = y
        self._mouse_in_window = True
        self.set_mouse_platform_visible()

        self.dispatch_event('on_mouse_enter', x, y)

        carbon.CallNextEventHandler(next_handler, ev)
        return noErr

    @CarbonEventHandler(kEventClassMouse, kEventMouseExited)
    def _on_mouse_exited(self, next_handler, ev, data):
        if not self._fullscreen:
            x, y = self._get_mouse_position(ev)
            y = self.height - y

            self._mouse_in_window = False
            self.set_mouse_platform_visible()

            self.dispatch_event('on_mouse_leave', x, y)

        carbon.CallNextEventHandler(next_handler, ev)
        return noErr

    @CarbonEventHandler(kEventClassMouse, kEventMouseWheelMoved)
    def _on_mouse_wheel_moved(self, next_handler, ev, data):

        x, y = self._get_mouse_position(ev)
        y = self.height - y

        axis = EventMouseWheelAxis()
        carbon.GetEventParameter(ev, kEventParamMouseWheelAxis,
            typeMouseWheelAxis, c_void_p(), sizeof(axis), c_void_p(),
            byref(axis))
        delta = c_long()
        carbon.GetEventParameter(ev, kEventParamMouseWheelDelta,
            typeSInt32, c_void_p(), sizeof(delta), c_void_p(),
            byref(delta))
        if axis.value == kEventMouseWheelAxisX:
            self.dispatch_event('on_mouse_scroll', 
                x, y, delta.value, 0)
        else:
            self.dispatch_event('on_mouse_scroll', 
                x, y, 0, delta.value)
                
        # _Don't_ call the next handler, which is application, as this then
        # calls our window handler again.
        #carbon.CallNextEventHandler(next_handler, ev)
        return noErr

    @CarbonEventHandler(kEventClassWindow, kEventWindowClose)
    def _on_window_close(self, next_handler, ev, data):
        self.dispatch_event('on_close')

        # Presumably the next event handler is the one that closes
        # the window; don't do that here. 
        #carbon.CallNextEventHandler(next_handler, ev)
        return noErr

    @CarbonEventHandler(kEventClassWindow, kEventWindowResizeCompleted)
    def _on_window_resize_completed(self, next_handler, ev, data):
        rect = Rect()
        carbon.GetWindowBounds(self._window, kWindowContentRgn, byref(rect))
        width = rect.right - rect.left
        height = rect.bottom - rect.top

        self.dispatch_event('on_resize', width, height)
        self.dispatch_event('on_expose')

        carbon.CallNextEventHandler(next_handler, ev)
        return noErr

    @CarbonEventHandler(kEventClassWindow, kEventWindowDragCompleted)
    def _on_window_drag_completed(self, next_handler, ev, data):
        rect = Rect()
        carbon.GetWindowBounds(self._window, kWindowContentRgn, byref(rect))

        self.dispatch_event('on_move', rect.left, rect.top)

        carbon.CallNextEventHandler(next_handler, ev)
        return noErr

    @CarbonEventHandler(kEventClassWindow, kEventWindowBoundsChanged)
    def _on_window_bounds_change(self, next_handler, ev, data):
        self._update_track_region()
        self._update_drawable()

        carbon.CallNextEventHandler(next_handler, ev)
        return noErr

    @CarbonEventHandler(kEventClassWindow, kEventWindowZoomed)
    def _on_window_zoomed(self, next_handler, ev, data):
        rect = Rect()
        carbon.GetWindowBounds(self._window, kWindowContentRgn, byref(rect))
        width = rect.right - rect.left
        height = rect.bottom - rect.top

        self.dispatch_event('on_move', rect.left, rect.top)
        self.dispatch_event('on_resize', width, height)

        carbon.CallNextEventHandler(next_handler, ev)
        return noErr

    @CarbonEventHandler(kEventClassWindow, kEventWindowActivated)
    def _on_window_activated(self, next_handler, ev, data):
        self.dispatch_event('on_activate')

        carbon.CallNextEventHandler(next_handler, ev)
        return noErr

    @CarbonEventHandler(kEventClassWindow, kEventWindowDeactivated)
    def _on_window_deactivated(self, next_handler, ev, data):
        self.dispatch_event('on_deactivate')

        carbon.CallNextEventHandler(next_handler, ev)
        return noErr
        
    @CarbonEventHandler(kEventClassWindow, kEventWindowShown)
    @CarbonEventHandler(kEventClassWindow, kEventWindowExpanded)
    def _on_window_shown(self, next_handler, ev, data):
        self._update_drawable()
        self.dispatch_event('on_show')

        carbon.CallNextEventHandler(next_handler, ev)
        return noErr

    @CarbonEventHandler(kEventClassWindow, kEventWindowHidden)
    @CarbonEventHandler(kEventClassWindow, kEventWindowCollapsed)
    def _on_window_hidden(self, next_handler, ev, data):
        self.dispatch_event('on_hide')

        carbon.CallNextEventHandler(next_handler, ev)
        return noErr

    @CarbonEventHandler(kEventClassWindow, kEventWindowDrawContent)
    def _on_window_draw_content(self, next_handler, ev, data):
        self.dispatch_event('on_expose')

        carbon.CallNextEventHandler(next_handler, ev)
        return noErr
        

       
def _create_cfstring(text):
    return carbon.CFStringCreateWithCString(c_void_p(), 
                                            text.encode('utf8'),
                                            kCFStringEncodingUTF8)

def _oscheck(result):
    if result != noErr:
        raise 'Carbon error %d' % result
    return result

def _aglcheck():
    err = agl.aglGetError()
    if err != agl.AGL_NO_ERROR:
        raise CarbonException(cast(agl.aglErrorString(err), c_char_p).value)

