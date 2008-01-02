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
import unicodedata
import warnings
import sys

if sys.platform not in ('cygwin', 'win32'):
    raise ImportError('Not a win32 platform.')

import pyglet
from pyglet.window import Platform, Display, Screen, BaseWindow, \
    WindowException, MouseCursor, DefaultMouseCursor, _PlatformEventHandler
from pyglet.window import event
from pyglet.event import EventDispatcher
from pyglet.window import key
from pyglet.window import mouse
from pyglet.window.win32.constants import *
from pyglet.window.win32.winkey import *
from pyglet.window.win32.types import *

from pyglet import gl
from pyglet.gl import gl_info
from pyglet.gl import glu_info
from pyglet.gl import wgl
from pyglet.gl import wglext_arb
from pyglet.gl import wgl_info

_debug_win32 = pyglet.options['debug_win32']

if _debug_win32:
    import traceback
    _GetLastError = windll.kernel32.GetLastError
    _SetLastError = windll.kernel32.SetLastError
    _FormatMessageA = windll.kernel32.FormatMessageA

    _log_win32 = open('debug_win32.log', 'w')
    
    def format_error(err):
        msg = create_string_buffer(256)
        _FormatMessageA(FORMAT_MESSAGE_FROM_SYSTEM,
                          c_void_p(),
                          err,
                          0,
                          msg,
                          len(msg),
                          c_void_p())
        return msg.value
    
    class DebugLibrary(object):
        def __init__(self, lib):
            self.lib = lib

        def __getattr__(self, name):
            fn = getattr(self.lib, name)
            def f(*args):
                _SetLastError(0)
                result = fn(*args)
                err = _GetLastError()
                if err != 0:
                    map(_log_win32.write,
                        traceback.format_list(traceback.extract_stack()[:-1]))
                    print >> _log_win32, format_error(err)
                return result
            return f
else:
    DebugLibrary = lambda lib: lib
            
_gdi32 = DebugLibrary(windll.gdi32)
_kernel32 = DebugLibrary(windll.kernel32)
_user32 = DebugLibrary(windll.user32)


_user32.GetKeyState.restype = c_short
_gdi32.CreateDIBitmap.argtypes = [HDC, POINTER(BITMAPINFOHEADER), DWORD,
    c_void_p, POINTER(BITMAPINFO), c_uint]

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
    (key.PAGEUP, True):     key.MOTION_BEGINNING_OF_FILE,
    (key.PAGEDOWN, True):   key.MOTION_END_OF_FILE,
    (key.BACKSPACE, False): key.MOTION_BACKSPACE,
    (key.DELETE, False):    key.MOTION_DELETE,
}

class Win32Exception(WindowException):
    pass

class Win32Platform(Platform):
    _display = None

    def get_default_display(self):
        if not self._display:
            self._display = Win32Display()
        return self._display
    
class Win32Display(Display):
    def get_screens(self):
        screens = []
        def enum_proc(hMonitor, hdcMonitor, lprcMonitor, dwData):
            r = lprcMonitor.contents
            width = r.right - r.left
            height = r.bottom - r.top
            screens.append(
                Win32Screen(self, hMonitor, r.left, r.top, width, height))
            return True
        enum_proc_type = WINFUNCTYPE(BOOL, HMONITOR, HDC, POINTER(RECT), LPARAM)
        enum_proc_ptr = enum_proc_type(enum_proc)
        _user32.EnumDisplayMonitors(NULL, NULL, enum_proc_ptr, 0)
        return screens

class Win32Screen(Screen):
    def __init__(self, display, handle, x, y, width, height):
        super(Win32Screen, self).__init__(x, y, width, height)
        self.display = display        
        self._handle = handle

    def get_matching_configs(self, template):
        # Determine which technique should be used for finding matching configs.        
        # Use the builtin PIXELFORMATDESCRIPTOR if possible, otherwise resort
        # to the WGL_ARB_pixel_format extension.
        need_pixel_format_arb = False
        if template.sample_buffers or template.samples:
            need_pixel_format_arb = True
            
        if need_pixel_format_arb:
            # Need a GL context before we can query WGL extensions.
            dummy_window = None
            if not gl_info.have_context():
                # Create a dummy context
                config = self.get_best_config()
                context = config.create_context(None)
                dummy_window = Win32Window(visible=False, context=context)
            
            try:            
                # Check for required extensions
                if not wgl_info.have_extension('WGL_ARB_pixel_format'):
                    return []
                return self._get_arb_pixel_format_matching_configs(template)
            finally:
                if dummy_window:
                    dummy_window.close()

        return self._get_pixel_format_descriptor_matching_configs(template)

    def _get_pixel_format_descriptor_matching_configs(self, template):
        '''Get matching configs using standard PIXELFORMATDESCRIPTOR
        technique.'''
        pfd = PIXELFORMATDESCRIPTOR()
        pfd.nSize = sizeof(PIXELFORMATDESCRIPTOR)
        pfd.nVersion = 1
        pfd.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL

        if template.double_buffer:
            pfd.dwFlags |= PFD_DOUBLEBUFFER
        else:
            pfd.dwFlags |= PFD_DOUBLEBUFFER_DONTCARE

        if template.stereo:
            pfd.dwFlags |= PFD_STEREO
        else:
            pfd.dwFlags |= PFD_STEREO_DONTCARE

        '''Not supported in pyglet API        
        if attributes.get('swap_copy', False):
            pfd.dwFlags |= PFD_SWAP_COPY
        if attributes.get('swap_exchange', False):
            pfd.dwFlags |= PFD_SWAP_EXCHANGE
        '''

        if not template.depth_size:
            pfd.dwFlags |= PFD_DEPTH_DONTCARE

        pfd.iPixelType = PFD_TYPE_RGBA
        pfd.cColorBits = template.buffer_size or 0
        pfd.cRedBits = template.red_size or 0
        pfd.cGreenBits = template.green_size or 0
        pfd.cBlueBits = template.blue_size or 0
        pfd.cAlphaBits = template.alpha_size or 0
        pfd.cAccumRedBits = template.accum_red_size or 0
        pfd.cAccumGreenBits = template.accum_green_size or 0
        pfd.cAccumBlueBits = template.accum_blue_size or 0
        pfd.cAccumAlphaBits = template.accum_alpha_size or 0
        pfd.cDepthBits = template.depth_size or 0
        pfd.cStencilBits = template.stencil_size or 0
        pfd.cAuxBuffers = template.aux_buffers or 0

        # No window created yet, so lets create a config based on
        # the DC of the entire screen.
        hdc = _user32.GetDC(0)

        pf = _gdi32.ChoosePixelFormat(hdc, byref(pfd))
        if pf:
            return [Win32Config(self, hdc, pf)]
        else:
            return []                    
                    
    def _get_arb_pixel_format_matching_configs(self, template):
        '''Get configs using the WGL_ARB_pixel_format extension.
        This method assumes a (dummy) GL context is already created.'''
        
        # Check for required extensions        
        if template.sample_buffers or template.samples:
            if not gl_info.have_extension('GL_ARB_multisample'):
                return []

        # Construct array of attributes
        attrs = []
        for name, value in template.get_gl_attributes():
            attr = Win32ConfigARB.attribute_ids.get(name, None)
            if attr and value is not None:
                attrs.extend([attr, int(value)])
        attrs.append(0)        
        attrs = (c_int * len(attrs))(*attrs)

        hdc = _user32.GetDC(0)     
        
        pformats = (c_int * 16)()
        nformats = c_uint(16)
        wglext_arb.wglChoosePixelFormatARB(hdc, attrs, None, 
                                           nformats, pformats, nformats)

        formats = [Win32ConfigARB(self, hdc, pf) \
                   for pf in pformats[:nformats.value]]
        return formats

class Win32Config(gl.Config):
    def __init__(self, screen, hdc, pf):
        self.screen = screen
        self._hdc = hdc
        self._pf = pf
        self._pfd = PIXELFORMATDESCRIPTOR()
        _gdi32.DescribePixelFormat(self._hdc, 
            self._pf, sizeof(PIXELFORMATDESCRIPTOR), byref(self._pfd))

        self.double_buffer = bool(self._pfd.dwFlags & PFD_DOUBLEBUFFER)
        self.sample_buffers = 0
        self.samples = 0
        self.stereo = bool(self._pfd.dwFlags & PFD_STEREO)
        self.buffer_size = self._pfd.cColorBits
        self.red_size = self._pfd.cRedBits
        self.green_size = self._pfd.cGreenBits
        self.blue_size = self._pfd.cBlueBits
        self.alpha_size = self._pfd.cAlphaBits
        self.accum_red_size = self._pfd.cAccumRedBits
        self.accum_green_size = self._pfd.cAccumGreenBits
        self.accum_blue_size = self._pfd.cAccumBlueBits
        self.accum_alpha_size = self._pfd.cAccumAlphaBits
        self.depth_size = self._pfd.cDepthBits
        self.stencil_size = self._pfd.cStencilBits
        self.aux_buffers = self._pfd.cAuxBuffers

    def create_context(self, share):
        # The context can't be created until we have the DC of the
        # window.  It's _possible_ that this could screw things up
        # (for example, destroying the share context before the new
        # window is created), but these are unlikely and not in the
        # ordinary workflow.
        return Win32Context(self, share)

    def is_complete(self):
        return True

class Win32ConfigARB(Win32Config):    
    attribute_ids = {
        'double_buffer': wglext_arb.WGL_DOUBLE_BUFFER_ARB,
        'stereo': wglext_arb.WGL_STEREO_ARB,
        'buffer_size': wglext_arb.WGL_COLOR_BITS_ARB,
        'aux_buffers': wglext_arb.WGL_AUX_BUFFERS_ARB,
        'sample_buffers': wglext_arb.WGL_SAMPLE_BUFFERS_ARB,
        'samples': wglext_arb.WGL_SAMPLES_ARB,
        'red_size': wglext_arb.WGL_RED_BITS_ARB,
        'green_size': wglext_arb.WGL_GREEN_BITS_ARB,
        'blue_size': wglext_arb.WGL_BLUE_BITS_ARB,
        'alpha_size': wglext_arb.WGL_ALPHA_BITS_ARB,
        'depth_size': wglext_arb.WGL_DEPTH_BITS_ARB,
        'stencil_size': wglext_arb.WGL_STENCIL_BITS_ARB,
        'accum_red_size': wglext_arb.WGL_ACCUM_RED_BITS_ARB,
        'accum_green_size': wglext_arb.WGL_ACCUM_GREEN_BITS_ARB,
        'accum_blue_size': wglext_arb.WGL_ACCUM_BLUE_BITS_ARB,
        'accum_alpha_size': wglext_arb.WGL_ACCUM_ALPHA_BITS_ARB,
    }
    def __init__(self, screen, hdc, pf):
        self.screen = screen
        self._hdc = hdc
        self._pf = pf
        
        names, attrs = map(None, *self.attribute_ids.items())
        attrs = (c_int * len(attrs))(*attrs)
        values = (c_int * len(attrs))()
        
        result = wglext_arb.wglGetPixelFormatAttribivARB(hdc,
            pf, 0, len(attrs), attrs, values)

        for name, value in zip(names, values):
            setattr(self, name, value)

    def create_context(self, share):
        return Win32ContextARB(self, share)

class Win32Context(gl.Context):
    _context = None
    def __init__(self, config, share):
        super(Win32Context, self).__init__(share)
        self.config = config
        self._share = share

    def _set_window(self, window):
        assert self._context is None
        _gdi32.SetPixelFormat(
            window._dc, self.config._pf, byref(self.config._pfd))
        self._context = wgl.wglCreateContext(window._dc)
        if self._share:
            assert self._share._context is not None
            wgl.wglShareLists(self._share._context, self._context)

    def destroy(self):
        super(Win32Context, self).destroy()
        wgl.wglDeleteContext(self._context)

class Win32ContextARB(Win32Context):
    def _set_window(self, window):
        assert self._context is None
        _gdi32.SetPixelFormat(window._dc, self.config._pf, None)
        self._context = wgl.wglCreateContext(window._dc)
        if self._share:
            assert self._share._context is not None
            wgl.wglShareLists(self._share._context, self._context)

class Win32MouseCursor(MouseCursor):
    drawable = False
    def __init__(self, cursor):
        self.cursor = cursor

# This is global state, we have to be careful not to set the same state twice,
# which will throw off the ShowCursor counter.
_win32_cursor_visible = True

Win32EventHandler = _PlatformEventHandler

class Win32Window(BaseWindow):
    _window_class = None
    _hwnd = None
    _dc = None
    _wgl_context = None
    _tracking = False
    _hidden = False
    _has_focus = False

    _exclusive_keyboard = False
    _exclusive_keyboard_focus = True
    _exclusive_mouse = False
    _exclusive_mouse_focus = True
    _exclusive_mouse_screen = None
    _exclusive_mouse_client = None
    _mouse_platform_visible = True

    _ws_style = 0
    _ex_ws_style = 0
    _minimum_size = None
    _maximum_size = None

    def __init__(self, *args, **kwargs):
        # Bind event handlers
        self._event_handlers = {}
        for func_name in self._platform_event_names:
            if not hasattr(self, func_name):
                continue
            func = getattr(self, func_name)
            for message in func._platform_event_data:
                self._event_handlers[message] = func

        super(Win32Window, self).__init__(*args, **kwargs)
        
    def _recreate(self, changes):
        if 'context' in changes:
            self._wgl_context = None
        
        self._create()

    def _create(self):
        # Ensure style is set before determining width/height.
        if self._fullscreen:
            self._ws_style = WS_POPUP
            self._ex_ws_style = 0 # WS_EX_TOPMOST
        else:
            styles = {
                self.WINDOW_STYLE_DEFAULT: (WS_OVERLAPPEDWINDOW, 0),
                self.WINDOW_STYLE_DIALOG:  (WS_OVERLAPPED|WS_CAPTION|WS_SYSMENU,
                                            WS_EX_DLGMODALFRAME),
                self.WINDOW_STYLE_TOOL:    (WS_OVERLAPPED|WS_CAPTION|WS_SYSMENU,
                                            WS_EX_TOOLWINDOW),
                self.WINDOW_STYLE_BORDERLESS: (WS_POPUP, 0),
            }
            self._ws_style, self._ex_ws_style = styles[self._style]

        if self._resizable and not self._fullscreen:
            self._ws_style |= WS_THICKFRAME
        else:
            self._ws_style &= ~(WS_THICKFRAME|WS_MAXIMIZEBOX)

        width, height = self._client_to_window_size(self._width, self._height)

        if not self._window_class:
            module = _kernel32.GetModuleHandleW(None)
            white = _gdi32.GetStockObject(WHITE_BRUSH)
            self._window_class = WNDCLASS()
            self._window_class.lpszClassName = u'GenericAppClass%d' % id(self)
            self._window_class.lpfnWndProc = WNDPROC(self._wnd_proc)
            self._window_class.style = CS_VREDRAW | CS_HREDRAW
            self._window_class.hInstance = 0
            self._window_class.hIcon = _user32.LoadIconW(module, 1)
            self._window_class.hbrBackground = white
            self._window_class.lpszMenuName = None
            self._window_class.cbClsExtra = 0
            self._window_class.cbWndExtra = 0
            _user32.RegisterClassW(byref(self._window_class))
        
        if not self._hwnd:
            self._hwnd = _user32.CreateWindowExW(
                self._ex_ws_style,
                self._window_class.lpszClassName,
                u'',
                self._ws_style,
                CW_USEDEFAULT,
                CW_USEDEFAULT,
                width,
                height,
                0,
                0,
                self._window_class.hInstance,
                0)

            self._dc = _user32.GetDC(self._hwnd)
        else:
            # Window already exists, update it with new style

            # We need to hide window here, otherwise Windows forgets
            # to redraw the whole screen after leaving fullscreen.
            _user32.ShowWindow(self._hwnd, SW_HIDE)

            _user32.SetWindowLongW(self._hwnd,
                GWL_STYLE,
                self._ws_style)
            _user32.SetWindowLongW(self._hwnd,
                GWL_EXSTYLE,
                self._ex_ws_style)

        if self._fullscreen:
            hwnd_after = HWND_TOPMOST
        else:
            hwnd_after = HWND_NOTOPMOST

        # Position and size window
        if self._fullscreen:
            _user32.SetWindowPos(self._hwnd, hwnd_after,
                self._screen.x, self._screen.y, width, height, SWP_FRAMECHANGED)
        elif False: # TODO location not in pyglet API
            x, y = self._client_to_window_pos(*factory.get_location())
            _user32.SetWindowPos(self._hwnd, hwnd_after,
                x, y, width, height, SWP_FRAMECHANGED)
        else:
            _user32.SetWindowPos(self._hwnd, hwnd_after,
                0, 0, width, height, SWP_NOMOVE | SWP_FRAMECHANGED)

        # Context must be created after window is created.
        if not self._wgl_context:
            self.context._set_window(self)
            self._wgl_context = self.context._context

        self.set_caption(self._caption)

        self.switch_to()
        self.set_vsync(self._vsync)

        if self._visible:
            self.set_visible()
            self.dispatch_event('on_expose')

    def close(self):
        super(Win32Window, self).close()
        _user32.DestroyWindow(self._hwnd)
        _user32.UnregisterClassW(self._window_class.lpszClassName, 0)
        self.set_mouse_platform_visible(True)
        self._hwnd = None
        self._dc = None
        self._wgl_context = None

    def _get_vsync(self):
        if wgl_info.have_extension('WGL_EXT_swap_control'):
            return bool(wglext_arb.wglGetSwapIntervalEXT())
    vsync = property(_get_vsync) # overrides BaseWindow property

    def set_vsync(self, vsync):
        if pyglet.options['vsync'] is not None:
            vsync = pyglet.options['vsync']
        if wgl_info.have_extension('WGL_EXT_swap_control'):
            wglext_arb.wglSwapIntervalEXT(int(vsync))
        else:
            warnings.warn('Could not set vsync; unsupported extension.')

    def switch_to(self):
        wgl.wglMakeCurrent(self._dc, self._wgl_context)
        self._context.set_current()
        gl_info.set_active_context()
        glu_info.set_active_context()

    def flip(self):
        self.draw_mouse_cursor()
        wgl.wglSwapLayerBuffers(self._dc, wgl.WGL_SWAP_MAIN_PLANE)

    def set_location(self, x, y):
        x, y = self._client_to_window_pos(x, y)
        _user32.SetWindowPos(self._hwnd, 0, x, y, 0, 0, 
            (SWP_NOZORDER |
             SWP_NOSIZE |
             SWP_NOOWNERZORDER))

    def get_location(self):
        rect = RECT()
        _user32.GetClientRect(self._hwnd, byref(rect))
        _user32.ClientToScreen(self._hwnd, byref(rect))
        return rect.left, rect.top

    def set_size(self, width, height):
        if self._fullscreen:
            raise WindowException('Cannot set size of fullscreen window.')
        width, height = self._client_to_window_size(width, height)
        _user32.SetWindowPos(self._hwnd, 0, 0, 0, width, height,
            (SWP_NOZORDER |
             SWP_NOMOVE |
             SWP_NOOWNERZORDER))

    def get_size(self):
        rect = RECT()
        _user32.GetClientRect(self._hwnd, byref(rect))
        return rect.right - rect.left, rect.bottom - rect.top

    def set_minimum_size(self, width, height):
        self._minimum_size = width, height

    def set_maximum_size(self, width, height):
        self._maximum_size = width, height

    def activate(self):
        _user32.SetForegroundWindow(self._hwnd)

    def set_visible(self, visible=True):
        if visible:
            if self._fullscreen:
                _user32.SetWindowPos(self._hwnd, HWND_TOPMOST, 0, 0, 0, 0,
                    SWP_NOMOVE | SWP_NOSIZE | SWP_SHOWWINDOW)
            else:
                _user32.ShowWindow(self._hwnd, SW_SHOW)
            self.dispatch_event('on_show')
            self.activate()
        else:
            _user32.ShowWindow(self._hwnd, SW_HIDE)
            self.dispatch_event('on_hide')
        self._visible = visible
        self.set_mouse_platform_visible()

    def minimize(self):
        _user32.ShowWindow(self._hwnd, SW_MINIMIZE)

    def maximize(self):
        _user32.ShowWindow(self._hwnd, SW_MAXIMIZE)

    def set_caption(self, caption):
        self._caption = caption
        _user32.SetWindowTextW(self._hwnd, c_wchar_p(caption))

    def set_mouse_platform_visible(self, platform_visible=None):
        if platform_visible is None:
            platform_visible = (self._mouse_visible and
                                not self._exclusive_mouse and
                                not self._mouse_cursor.drawable) or \
                               (not self._mouse_in_window or 
                                not self._has_focus)

        if platform_visible and not self._mouse_cursor.drawable:
            if isinstance(self._mouse_cursor, Win32MouseCursor):
                cursor = self._mouse_cursor.cursor
            else:
                cursor = _user32.LoadCursorW(None, IDC_ARROW)
            _user32.SetClassLongW(self._hwnd, GCL_HCURSOR, cursor)
            _user32.SetCursor(cursor)

        if platform_visible == self._mouse_platform_visible:
            return

        # Avoid calling ShowCursor with the current visibility (which would
        # push the counter too far away from zero).
        global _win32_cursor_visible
        if _win32_cursor_visible != platform_visible:
            _user32.ShowCursor(platform_visible)
            _win32_cursor_visible = platform_visible

        self._mouse_platform_visible = platform_visible

    def _reset_exclusive_mouse_screen(self):
        '''Recalculate screen coords of mouse warp point for exclusive
        mouse.'''
        p = POINT()
        rect = RECT()
        _user32.GetClientRect(self._hwnd, byref(rect))
        _user32.MapWindowPoints(self._hwnd, HWND_DESKTOP, byref(rect), 2)
        p.x = (rect.left + rect.right) / 2
        p.y = (rect.top + rect.bottom) / 2

        # This is the point the mouse will be kept at while in exclusive
        # mode.
        self._exclusive_mouse_screen = p.x, p.y
        self._exclusive_mouse_client = p.x - rect.left, p.y - rect.top

    def set_exclusive_mouse(self, exclusive=True):
        if self._exclusive_mouse == exclusive and \
           self._exclusive_mouse_focus == self._has_focus:
            return
    
        if exclusive and self._has_focus:
            # Move mouse to the center of the window.
            self._reset_exclusive_mouse_screen()
            _user32.SetCursorPos(*self._exclusive_mouse_screen)

            # Clip to client area, to prevent large mouse movements taking
            # it outside the client area.
            rect = RECT()
            _user32.GetClientRect(self._hwnd, byref(rect))
            _user32.MapWindowPoints(self._hwnd, HWND_DESKTOP, byref(rect), 2)
            _user32.ClipCursor(byref(rect))
        else:
            # Release clip
            _user32.ClipCursor(c_void_p())

        self._exclusive_mouse = exclusive
        self._exclusive_mouse_focus = self._has_focus
        self.set_mouse_platform_visible()

    def set_exclusive_keyboard(self, exclusive=True):
        if self._exclusive_keyboard == exclusive and \
           self._exclusive_keyboard_focus == self._has_focus:
            return

        if exclusive and self._has_focus:
            _user32.RegisterHotKey(self._hwnd, 0, WIN32_MOD_ALT, VK_TAB)
        else:
            _user32.UnregisterHotKey(self._hwnd, 0)

        self._exclusive_keyboard = exclusive
        self._exclusive_keyboard_focus = self._has_focus

    def get_system_mouse_cursor(self, name):
        if name == self.CURSOR_DEFAULT:
            return DefaultMouseCursor()

        names = {
            self.CURSOR_CROSSHAIR:       IDC_CROSS,
            self.CURSOR_HAND:            IDC_HAND,
            self.CURSOR_HELP:            IDC_HELP,
            self.CURSOR_NO:              IDC_NO,
            self.CURSOR_SIZE:            IDC_SIZEALL,
            self.CURSOR_SIZE_UP:         IDC_SIZENS,
            self.CURSOR_SIZE_UP_RIGHT:   IDC_SIZENESW,
            self.CURSOR_SIZE_RIGHT:      IDC_SIZEWE,
            self.CURSOR_SIZE_DOWN_RIGHT: IDC_SIZENWSE,
            self.CURSOR_SIZE_DOWN:       IDC_SIZENS,
            self.CURSOR_SIZE_DOWN_LEFT:  IDC_SIZENESW,
            self.CURSOR_SIZE_LEFT:       IDC_SIZEWE,
            self.CURSOR_SIZE_UP_LEFT:    IDC_SIZENWSE,
            self.CURSOR_SIZE_UP_DOWN:    IDC_SIZENS,
            self.CURSOR_SIZE_LEFT_RIGHT: IDC_SIZEWE,
            self.CURSOR_TEXT:            IDC_IBEAM,
            self.CURSOR_WAIT:            IDC_WAIT,
            self.CURSOR_WAIT_ARROW:      IDC_APPSTARTING,
        }
        if name not in names:
            raise Win32Exception('Unknown cursor name "%s"' % name)
        cursor = _user32.LoadCursorW(None, names[name])
        return Win32MouseCursor(cursor)

    def set_icon(self, *images):
        # XXX Undocumented AFAICT, but XP seems happy to resize an image
        # of any size, so no scaling necessary.

        def best_image(width, height):
            # A heuristic for finding closest sized image to required size.
            image = images[0]
            for img in images:
                if img.width == width and img.height == height:
                    # Exact match always used
                    return img
                elif img.width >= width and \
                     img.width * img.height > image.width * image.height:
                    # At least wide enough, and largest area
                    image = img
            return image

        def get_icon(image):
            # Alpha-blended icon: see http://support.microsoft.com/kb/318876
            image.format = 'BGRA'
            image.pitch = len(image.format) * image.width

            header = BITMAPV5HEADER()
            header.bV5Size = sizeof(header)
            header.bV5Width = image.width
            header.bV5Height = image.height
            header.bV5Planes = 1
            header.bV5BitCount = 32
            header.bV5Compression = BI_BITFIELDS
            header.bV5RedMask = 0x00ff0000
            header.bV5GreenMask = 0x0000ff00
            header.bV5BlueMask = 0x000000ff
            header.bV5AlphaMask = 0xff000000

            hdc = _user32.GetDC(None)
            dataptr = c_void_p()
            bitmap = _gdi32.CreateDIBSection(hdc, byref(header), DIB_RGB_COLORS,
                byref(dataptr), None, 0)
            _user32.ReleaseDC(None, hdc)

            memmove(dataptr, image.data, len(image.data))

            mask = _gdi32.CreateBitmap(image.width, image.height, 1, 1, None)

            iconinfo = ICONINFO()
            iconinfo.fIcon = True
            iconinfo.hbmMask = mask
            iconinfo.hbmColor = bitmap
            icon = _user32.CreateIconIndirect(byref(iconinfo))

            _gdi32.DeleteObject(mask)
            _gdi32.DeleteObject(bitmap)
            
            return icon

        # Set large icon
        image = best_image(_user32.GetSystemMetrics(SM_CXICON),
                           _user32.GetSystemMetrics(SM_CYICON))
        icon = get_icon(image)
        _user32.SetClassLongW(self._hwnd, GCL_HICON, icon)

        # Set small icon
        image = best_image(_user32.GetSystemMetrics(SM_CXSMICON),
                           _user32.GetSystemMetrics(SM_CYSMICON))
        icon = get_icon(image)
        _user32.SetClassLongW(self._hwnd, GCL_HICONSM, icon)

    # Private util

    def _client_to_window_size(self, width, height):
        rect = RECT()
        rect.left = 0
        rect.top = 0
        rect.right = width
        rect.bottom = height
        _user32.AdjustWindowRectEx(byref(rect),
            self._ws_style, False, self._ex_ws_style)
        return rect.right - rect.left, rect.bottom - rect.top

    def _client_to_window_pos(self, x, y):
        rect = RECT()
        rect.left = x
        rect.top = y
        _user32.AdjustWindowRectEx(byref(rect),
            self._ws_style, False, self._ex_ws_style)
        return rect.left, rect.top

    # Event dispatching

    def dispatch_events(self):
        self._allow_dispatch_event = True
        while self._event_queue:
            event = self._event_queue.pop(0)
            if type(event[0]) is str:
                # pyglet event
                EventDispatcher.dispatch_event(self, *event)
            else:
                # win32 event
                event[0](*event[1:])

        msg = MSG()
        while _user32.PeekMessageW(byref(msg), 0, 0, 0, PM_REMOVE):
            _user32.TranslateMessage(byref(msg))
            _user32.DispatchMessageW(byref(msg))
        self._allow_dispatch_event = False

    def _wnd_proc(self, hwnd, msg, wParam, lParam):
        event_handler = self._event_handlers.get(msg, None)
        result = 0
        if event_handler:
            if self._allow_dispatch_event:
                result = event_handler(msg, wParam, lParam)
            else:
                self._event_queue.append((event_handler, msg, wParam, lParam))
                result = 0
        if not result and msg != WM_CLOSE:
            result = _user32.DefWindowProcW(c_int(hwnd), c_int(msg),
                c_int(wParam), c_int(lParam)) 
        return result

    # Event handlers

    def _get_modifiers(self, key_lParam=0):
        modifiers = 0
        if _user32.GetKeyState(VK_SHIFT) & 0xff00:
            modifiers |= key.MOD_SHIFT
        if _user32.GetKeyState(VK_CONTROL) & 0xff00:
            modifiers |= key.MOD_CTRL
        if _user32.GetKeyState(VK_LWIN) & 0xff00:
            modifiers |= key.MOD_WINDOWS
        if _user32.GetKeyState(VK_CAPITAL) & 0x00ff:    # toggle
            modifiers |= key.MOD_CAPSLOCK
        if _user32.GetKeyState(VK_NUMLOCK) & 0x00ff:    # toggle
            modifiers |= key.MOD_NUMLOCK
        if _user32.GetKeyState(VK_SCROLL) & 0x00ff:    # toggle
            modifiers |= key.MOD_SCROLLLOCK
        if key_lParam:
            if key_lParam & (1 << 29):
                modifiers |= key.MOD_ALT
        elif _user32.GetKeyState(VK_MENU) < 0:
            modifiers |= key.MOD_ALT
        return modifiers

    @staticmethod
    def _get_location(lParam):
        x = c_int16(lParam & 0xffff).value
        y = c_int16(lParam >> 16).value
        return x, y

    @Win32EventHandler(WM_KEYDOWN)
    @Win32EventHandler(WM_KEYUP)
    @Win32EventHandler(WM_SYSKEYDOWN)
    @Win32EventHandler(WM_SYSKEYUP)
    def _event_key(self, msg, wParam, lParam):
        repeat = False
        if lParam & (1 << 30):
            if msg not in (WM_KEYUP, WM_SYSKEYUP):
                repeat = True
            ev = 'on_key_release'
        else:
            ev = 'on_key_press'

        symbol = keymap.get(wParam, None)
        if symbol is None:
            ch = _user32.MapVirtualKeyW(wParam, MAPVK_VK_TO_CHAR)
            symbol = chmap.get(ch)

        if symbol is None:
            symbol = key.user_key(wParam)
        elif symbol == key.LCTRL and lParam & (1 << 24):
            symbol = key.RCTRL
        elif symbol == key.LALT and lParam & (1 << 24):
            symbol = key.RALT
        elif symbol == key.LSHIFT:
            pass # TODO: some magic with getstate to find out if it's the
                 # right or left shift key. 

        modifiers = self._get_modifiers(lParam)
        
        if not repeat:
            self.dispatch_event(ev, symbol, modifiers)

        ctrl = modifiers & key.MOD_CTRL != 0
        if (symbol, ctrl) in _motion_map and msg not in (WM_KEYUP, WM_SYSKEYUP):
            motion = _motion_map[symbol, ctrl]
            if modifiers & key.MOD_SHIFT:
                self.dispatch_event('on_text_motion_select', motion)
            else:
                self.dispatch_event('on_text_motion', motion)

        # Send on to DefWindowProc if not exclusive.
        if self._exclusive_keyboard:
            return 0
        else:
            return None

    @Win32EventHandler(WM_CHAR)
    def _event_char(self, msg, wParam, lParam):
        text = unichr(wParam)
        if unicodedata.category(text) != 'Cc' or text == '\r':
            self.dispatch_event('on_text', text)
        return 0

    @Win32EventHandler(WM_MOUSEMOVE)
    def _event_mousemove(self, msg, wParam, lParam):
        x, y = self._get_location(lParam)

        if (x, y) == self._exclusive_mouse_client:
            # Ignore the event caused by SetCursorPos
            self._mouse_x = x
            self._mouse_y = y
            return 0

        y = self.height - y

        if self._exclusive_mouse and self._has_focus:
            # Reset mouse position (so we don't hit the edge of the screen).
            _user32.SetCursorPos(*self._exclusive_mouse_screen)
            
        dx = x - self._mouse_x
        dy = y - self._mouse_y

        if not self._tracking:
            # There is no WM_MOUSEENTER message (!), so fake it from the
            # first WM_MOUSEMOVE event after leaving.  Use self._tracking
            # to determine when to recreate the tracking structure after
            # re-entering (to track the next WM_MOUSELEAVE).
            self._mouse_in_window = True
            self.set_mouse_platform_visible()
            self.dispatch_event('on_mouse_enter', x, y)
            self._tracking = True
            track = TRACKMOUSEEVENT()
            track.cbSize = sizeof(track)
            track.dwFlags = TME_LEAVE
            track.hwndTrack = self._hwnd
            _user32.TrackMouseEvent(byref(track))

        self._mouse_x = x
        self._mouse_y = y
        
        buttons = 0
        if wParam & MK_LBUTTON:
            buttons |= mouse.LEFT
        if wParam & MK_MBUTTON:
            buttons |= mouse.MIDDLE
        if wParam & MK_RBUTTON:
            buttons |= mouse.RIGHT

        if buttons:
            # Drag event
            modifiers = self._get_modifiers()
            self.dispatch_event('on_mouse_drag', 
                x, y, dx, dy, buttons, modifiers)
        else:
            # Motion event
            self.dispatch_event('on_mouse_motion', x, y, dx, dy)
        return 0

    @Win32EventHandler(WM_MOUSELEAVE)
    def _event_mouseleave(self, msg, wParam, lParam):
        point = POINT()
        _user32.GetCursorPos(byref(point))
        _user32.ScreenToClient(self._hwnd, byref(point))
        x = point.x
        y = self.height - point.y
        self._tracking = False
        self._mouse_in_window = False
        self.set_mouse_platform_visible()
        self.dispatch_event('on_mouse_leave', x, y)
        return 0

    def _event_mousebutton(self, ev, button, lParam):
        if ev == 'on_mouse_press':
            _user32.SetCapture(self._hwnd)
        else:
            _user32.ReleaseCapture()
        x, y = self._get_location(lParam)
        y = self.height - y
        self.dispatch_event(ev, x, y, button, self._get_modifiers())
        return 0

    @Win32EventHandler(WM_LBUTTONDOWN)
    def _event_lbuttondown(self, msg, wParam, lParam):
        return self._event_mousebutton(
            'on_mouse_press', mouse.LEFT, lParam)

    @Win32EventHandler(WM_LBUTTONUP)
    def _event_lbuttonup(self, msg, wParam, lParam):
        return self._event_mousebutton(
            'on_mouse_release', mouse.LEFT, lParam)

    @Win32EventHandler(WM_MBUTTONDOWN)
    def _event_mbuttondown(self, msg, wParam, lParam):
        return self._event_mousebutton(
            'on_mouse_press', mouse.MIDDLE, lParam)

    @Win32EventHandler(WM_MBUTTONUP)
    def _event_mbuttonup(self, msg, wParam, lParam):
        return self._event_mousebutton(
            'on_mouse_release', mouse.MIDDLE, lParam)

    @Win32EventHandler(WM_RBUTTONDOWN)
    def _event_rbuttondown(self, msg, wParam, lParam):
        return self._event_mousebutton(
            'on_mouse_press', mouse.RIGHT, lParam)

    @Win32EventHandler(WM_RBUTTONUP)
    def _event_rbuttonup(self, msg, wParam, lParam):
        return self._event_mousebutton(
            'on_mouse_release', mouse.RIGHT, lParam)

    @Win32EventHandler(WM_MOUSEWHEEL)
    def _event_mousewheel(self, msg, wParam, lParam):
        delta = c_short(wParam >> 16).value
        self.dispatch_event('on_mouse_scroll', 
            self._mouse_x, self._mouse_y, 0, delta / WHEEL_DELTA)
        return 0

    @Win32EventHandler(WM_CLOSE)
    def _event_close(self, msg, wParam, lParam):
        self.dispatch_event('on_close')
        return 0

    @Win32EventHandler(WM_PAINT)
    def _event_paint(self, msg, wParam, lParam):
        self.dispatch_event('on_expose')
        # Validating the window using ValidateRect or ValidateRgn
        # doesn't clear the paint message when more than one window
        # is open [why?]; defer to DefWindowProc instead.
        return None

    @Win32EventHandler(WM_SIZE)
    def _event_size(self, msg, wParam, lParam):
        if not self._dc:
            # Ignore window creation size event (appears for fullscreen
            # only) -- we haven't got DC or HWND yet.
            return None

        if wParam == SIZE_MINIMIZED:
            # Minimized, not resized.
            self._hidden = True
            self.dispatch_event('on_hide')
            return 0
        if self._hidden:
            # Restored
            self._hidden = False
            self.dispatch_event('on_show')
        w, h = self._get_location(lParam)
        self._reset_exclusive_mouse_screen()
        self.dispatch_event('on_resize', w, h)
        return 0

    @Win32EventHandler(WM_MOVE)
    def _event_move(self, msg, wParam, lParam):
        x, y = self._get_location(lParam)
        self._reset_exclusive_mouse_screen()
        self.dispatch_event('on_move', x, y)
        return 0

    '''
    # Alternative to using WM_SETFOCUS and WM_KILLFOCUS.  Which
    # is better?

    @Win32EventHandler(WM_ACTIVATE)
    def _event_activate(self, msg, wParam, lParam):
        if wParam & 0xffff == WA_INACTIVE:
            self.dispatch_event('on_deactivate')
        else:
            self.dispatch_event('on_activate')
            _user32.SetFocus(self._hwnd)
        return 0
    '''

    @Win32EventHandler(WM_SETFOCUS)
    def _event_setfocus(self, msg, wParam, lParam):
        self.dispatch_event('on_activate')
        self._has_focus = True
        self.set_exclusive_keyboard(self._exclusive_keyboard)
        self.set_exclusive_mouse(self._exclusive_mouse)
        return 0

    @Win32EventHandler(WM_KILLFOCUS)
    def _event_killfocus(self, msg, wParam, lParam):
        self.dispatch_event('on_deactivate')
        self._has_focus = False
        self.set_exclusive_keyboard(self._exclusive_keyboard)
        self.set_exclusive_mouse(self._exclusive_mouse)
        return 0

    @Win32EventHandler(WM_GETMINMAXINFO)
    def _event_getminmaxinfo(self, msg, wParam, lParam):
        info = MINMAXINFO.from_address(lParam)
        if self._minimum_size:
            info.ptMinTrackSize.x, info.ptMinTrackSize.y = \
                self._client_to_window_size(*self._minimum_size)
        if self._maximum_size:
            info.ptMaxTrackSize.x, info.ptMaxTrackSize.y = \
                self._client_to_window_size(*self._maximum_size)
        return 0

    @Win32EventHandler(WM_ERASEBKGND)
    def _event_erasebkgnd(self, msg, wParam, lParam):
        # Prevent flicker during resize.
        return 0
