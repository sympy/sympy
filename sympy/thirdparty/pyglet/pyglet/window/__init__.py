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

'''Windowing and user-interface events.

This module allows applications to create and display windows with an
OpenGL context.  Windows can be created with a variety of border styles 
or set fullscreen.

You can register event handlers for keyboard, mouse and window events.
For games and kiosks you can also restrict the input to your windows,
for example disabling users from switching away from the application
with certain key combinations or capturing and hiding the mouse.

Getting started
---------------

Call the Window constructor to create a new window::

    from pyglet.window import Window
    win = Window(width=640, height=480)

Attach your own event handlers::

    @win.event
    def on_key_press(symbol, modifiers):
        # ... handle this event ...

Place drawing code for the window within the `Window.on_draw` event handler::

    @win.event
    def on_draw():
        # ... drawing code ...

Call `pyglet.app.run` to enter the main event loop (by default, this
returns when all open windows are closed)::

    from pyglet import app
    app.run()

Creating a game window
----------------------

Use `Window.set_exclusive_mouse` to hide the mouse cursor and receive relative
mouse movement events.  Specify ``fullscreen=True`` as a keyword argument to
the `Window` constructor to render to the entire screen rather than opening a
window::

    win = Window(fullscreen=True)
    win.set_exclusive_mouse()

Working with multiple screens
-----------------------------

By default, fullscreen windows are opened on the primary display (typically
set by the user in their operating system settings).  You can retrieve a list
of attached screens and select one manually if you prefer.  This is useful for
opening a fullscreen window on each screen::

    display = window.get_platform().get_default_display()
    screens = display.get_screens()
    windows = []
    for screen in screens:
        windows.append(window.Window(fullscreen=True, screen=screen))

Specifying a screen has no effect if the window is not fullscreen.

Specifying the OpenGL context properties
----------------------------------------

Each window has its own context which is created when the window is created.
You can specify the properties of the context before it is created
by creating a "template" configuration::

    from pyglet import gl
    # Create template config
    config = gl.Config()
    config.stencil_size = 8
    config.aux_buffers = 4
    # Create a window using this config
    win = window.Window(config=config)

To determine if a given configuration is supported, query the screen (see
above, "Working with multiple screens")::

    configs = screen.get_matching_configs(config)
    if not configs:
        # ... config is not supported
    else:
        win = window.Window(config=configs[0])

'''

__docformat__ = 'restructuredtext'
__version__ = '$Id: __init__.py 2215 2008-08-24 04:53:34Z Alex.Holkner $'

import pprint
import sys

import pyglet
from pyglet import gl
from pyglet.gl import gl_info
from pyglet.event import EventDispatcher
import pyglet.window.key

class WindowException(Exception):
    '''The root exception for all window-related errors.'''
    pass

class NoSuchDisplayException(WindowException):
    '''An exception indicating the requested display is not available.'''
    pass

class NoSuchConfigException(WindowException):
    '''An exception indicating the requested configuration is not
    available.'''
    pass

class MouseCursorException(WindowException):
    '''The root exception for all mouse cursor-related errors.'''
    pass

class Platform(object):
    '''Operating-system-level functionality.

    The platform instance can only be obtained with `get_platform`.  Use
    the platform to obtain a `Display` instance.
    '''
    def get_display(self, name):
        '''Get a display device by name.

        This is meaningful only under X11, where the `name` is a
        string including the host name and display number; for example
        ``"localhost:1"``.

        On platforms other than X11, `name` is ignored and the default
        display is returned.  pyglet does not support multiple multiple
        video devices on Windows or OS X.  If more than one device is
        attached, they will appear as a single virtual device comprising
        all the attached screens.

        :Parameters:
            `name` : str
                The name of the display to connect to.

        :rtype: `Display`
        '''
        return get_default_display()

    def get_default_display(self):
        '''Get the default display device.

        :rtype: `Display`
        '''
        raise NotImplementedError('abstract')

class Display(object):
    '''A display device supporting one or more screens.

    Use `Platform.get_display` or `Platform.get_default_display` to obtain
    an instance of this class.  Use a display to obtain `Screen` instances.
    '''
    def __init__(self):
        from pyglet import app
        app.displays.add(self)

    def get_screens(self):
        '''Get the available screens.

        A typical multi-monitor workstation comprises one `Display` with
        multiple `Screen` s.  This method returns a list of screens which
        can be enumerated to select one for full-screen display.

        For the purposes of creating an OpenGL config, the default screen
        will suffice.

        :rtype: list of `Screen`
        '''
        raise NotImplementedError('abstract')    

    def get_default_screen(self):
        '''Get the default screen as specified by the user's operating system
        preferences.

        :rtype: `Screen`
        '''
        return self.get_screens()[0]

    def get_windows(self):
        '''Get the windows currently attached to this display.

        :rtype: sequence of `Window`
        '''
        from pyglet import app
        return [window for window in app.windows if window.display is self]

class Screen(object):
    '''A virtual monitor that supports fullscreen windows.

    Screens typically map onto a physical display such as a
    monitor, television or projector.  Selecting a screen for a window
    has no effect unless the window is made fullscreen, in which case
    the window will fill only that particular virtual screen.

    The `width` and `height` attributes of a screen give the current
    resolution of the screen.  The `x` and `y` attributes give the global
    location of the top-left corner of the screen.  This is useful for
    determining if screens arranged above or next to one another.  
    
    You cannot always rely on the origin to give the placement of monitors.
    For example, an X server with two displays without Xinerama enabled
    will present two logically separate screens with no relation to each
    other.

    Use `Display.get_screens` or `Display.get_default_screen` to obtain an
    instance of this class.

    :Ivariables:
        `x` : int
            Left edge of the screen on the virtual desktop.
        `y` : int
            Top edge of the screen on the virtual desktop.
        `width` : int
            Width of the screen, in pixels.
        `height` : int
            Height of the screen, in pixels.

    '''

    def __init__(self, x, y, width, height):
        self.x = x
        self.y = y
        self.width = width
        self.height = height

    def __repr__(self):
        return '%s(x=%d, y=%d, width=%d, height=%d)' % \
            (self.__class__.__name__, self.x, self.y, self.width, self.height)

    def get_best_config(self, template=None):
        '''Get the best available GL config.

        Any required attributes can be specified in `template`.  If
        no configuration matches the template, `NoSuchConfigException` will
        be raised.

        :Parameters:
            `template` : `pyglet.gl.Config`
                A configuration with desired attributes filled in.

        :rtype: `pyglet.gl.Config`
        :return: A configuration supported by the platform that best
            fulfils the needs described by the template.
        '''
        if template is None:
            template = gl.Config()
        configs = self.get_matching_configs(template)
        if not configs:
            raise NoSuchConfigException()
        return configs[0]

    def get_matching_configs(self, template):
        '''Get a list of configs that match a specification.

        Any attributes specified in `template` will have values equal
        to or greater in each returned config.  If no configs satisfy
        the template, an empty list is returned.

        :Parameters:
            `template` : `pyglet.gl.Config`
                A configuration with desired attributes filled in.

        :rtype: list of `pyglet.gl.Config`
        :return: A list of matching configs.
        '''
        raise NotImplementedError('abstract')

class MouseCursor(object):
    '''An abstract mouse cursor.'''

    #: Indicates if the cursor is drawn using OpenGL.  This is True
    #: for all mouse cursors except system cursors.
    drawable = True

    def draw(self, x, y):
        '''Abstract render method.

        The cursor should be drawn with the "hot" spot at the given
        coordinates.  The projection is set to the pyglet default (i.e., 
        orthographic in window-space), however no other aspects of the 
        state can be assumed.

        :Parameters:
            `x` : int
                X coordinate of the mouse pointer's hot spot.
            `y` : int
                Y coordinate of the mouse pointer's hot spot.

        '''
        raise NotImplementedError('abstract')

class DefaultMouseCursor(MouseCursor):
    '''The default mouse cursor used by the operating system.'''
    drawable = False

class ImageMouseCursor(MouseCursor):
    '''A user-defined mouse cursor created from an image.

    Use this class to create your own mouse cursors and assign them
    to windows.  There are no constraints on the image size or format.
    '''
    drawable = True

    def __init__(self, image, hot_x=0, hot_y=0):
        '''Create a mouse cursor from an image.

        :Parameters:
            `image` : `pyglet.image.AbstractImage`
                Image to use for the mouse cursor.  It must have a
                valid ``texture`` attribute.
            `hot_x` : int
                X coordinate of the "hot" spot in the image relative to the
                image's anchor.
            `hot_y` : int
                Y coordinate of the "hot" spot in the image, relative to the
                image's anchor.
        '''
        self.texture = image.get_texture()
        self.hot_x = hot_x
        self.hot_y = hot_y

    def draw(self, x, y):
        gl.glPushAttrib(gl.GL_ENABLE_BIT | gl.GL_CURRENT_BIT)
        gl.glColor4f(1, 1, 1, 1)
        gl.glEnable(gl.GL_BLEND)
        gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
        self.texture.blit(x - self.hot_x, y - self.hot_y, 0)
        gl.glPopAttrib()

def _PlatformEventHandler(data):
    '''Decorator for platform event handlers.  
    
    Apply giving the platform-specific data needed by the window to associate
    the method with an event.  See platform-specific subclasses of this
    decorator for examples.

    The following attributes are set on the function, which is returned
    otherwise unchanged:

    _platform_event
        True
    _platform_event_data
        List of data applied to the function (permitting multiple decorators
        on the same method).
    '''
    def _event_wrapper(f):
        f._platform_event = True
        if not hasattr(f, '_platform_event_data'):
            f._platform_event_data = []
        f._platform_event_data.append(data)
        return f
    return _event_wrapper

class _WindowMetaclass(type):
    '''Sets the _platform_event_names class variable on the window
    subclass.
    '''
    def __init__(cls, name, bases, dict):
        cls._platform_event_names = set()
        for base in bases:
            if hasattr(base, '_platform_event_names'):
                cls._platform_event_names.update(base._platform_event_names)
        for name, func in dict.items():
            if hasattr(func, '_platform_event'):
                cls._platform_event_names.add(name)
        super(_WindowMetaclass, cls).__init__(name, bases, dict)

class BaseWindow(EventDispatcher):
    '''Platform-independent application window.

    A window is a "heavyweight" object occupying operating system resources.
    The "client" or "content" area of a window is filled entirely with
    an OpenGL viewport.  Applications have no access to operating system
    widgets or controls; all rendering must be done via OpenGL.

    Windows may appear as floating regions or can be set to fill an entire
    screen (fullscreen).  When floating, windows may appear borderless or
    decorated with a platform-specific frame (including, for example, the
    title bar, minimize and close buttons, resize handles, and so on).

    While it is possible to set the location of a window, it is recommended
    that applications allow the platform to place it according to local
    conventions.  This will ensure it is not obscured by other windows,
    and appears on an appropriate screen for the user.

    To render into a window, you must first call `switch_to`, to make
    it the current OpenGL context.  If you use only one window in the
    application, there is no need to do this.

    :Ivariables:
        `has_exit` : bool
            True if the user has attempted to close the window.

            :deprecated: Windows are closed immediately by the default
                `on_close` handler when `pyglet.app.event_loop` is being
                used.

    '''
    __metaclass__ = _WindowMetaclass

    # Filled in by metaclass with the names of all methods on this (sub)class
    # that are platform event handlers.
    _platform_event_names = set()

    #: The default window style.
    WINDOW_STYLE_DEFAULT = None
    #: The window style for pop-up dialogs.
    WINDOW_STYLE_DIALOG = 'dialog'
    #: The window style for tool windows.
    WINDOW_STYLE_TOOL = 'tool'
    #: A window style without any decoration.
    WINDOW_STYLE_BORDERLESS = 'borderless' 

    #: The default mouse cursor.
    CURSOR_DEFAULT = None
    #: A crosshair mouse cursor.
    CURSOR_CROSSHAIR = 'crosshair'
    #: A pointing hand mouse cursor.
    CURSOR_HAND = 'hand'
    #: A "help" mouse cursor; typically a question mark and an arrow.
    CURSOR_HELP = 'help'
    #: A mouse cursor indicating that the selected operation is not permitted.
    CURSOR_NO = 'no'
    #: A mouse cursor indicating the element can be resized.
    CURSOR_SIZE = 'size'
    #: A mouse cursor indicating the element can be resized from the top
    #: border.
    CURSOR_SIZE_UP = 'size_up'
    #: A mouse cursor indicating the element can be resized from the
    #: upper-right corner.
    CURSOR_SIZE_UP_RIGHT = 'size_up_right'
    #: A mouse cursor indicating the element can be resized from the right 
    #: border.
    CURSOR_SIZE_RIGHT = 'size_right'
    #: A mouse cursor indicating the element can be resized from the lower-right
    #: corner.
    CURSOR_SIZE_DOWN_RIGHT = 'size_down_right'
    #: A mouse cursor indicating the element can be resized from the bottom 
    #: border.
    CURSOR_SIZE_DOWN = 'size_down'
    #: A mouse cursor indicating the element can be resized from the lower-left
    #: corner.
    CURSOR_SIZE_DOWN_LEFT = 'size_down_left'
    #: A mouse cursor indicating the element can be resized from the left 
    #: border.
    CURSOR_SIZE_LEFT = 'size_left'
    #: A mouse cursor indicating the element can be resized from the upper-left
    #: corner.
    CURSOR_SIZE_UP_LEFT = 'size_up_left'
    #: A mouse cursor indicating the element can be resized vertically.
    CURSOR_SIZE_UP_DOWN = 'size_up_down'
    #: A mouse cursor indicating the element can be resized horizontally.
    CURSOR_SIZE_LEFT_RIGHT = 'size_left_right'
    #: A text input mouse cursor (I-beam).
    CURSOR_TEXT = 'text'
    #: A "wait" mouse cursor; typically an hourglass or watch.
    CURSOR_WAIT = 'wait'
    #: The "wait" mouse cursor combined with an arrow.
    CURSOR_WAIT_ARROW = 'wait_arrow'

    has_exit = False

    #: Window display contents validity.  The `pyglet.app` event loop
    #: examines every window each iteration and only dispatches the `on_draw`
    #: event to windows that have `invalid` set.  By default, windows always
    #: have `invalid` set to ``True``.
    #:
    #: You can prevent redundant redraws by setting this variable to ``False``
    #: in the window's `on_draw` handler, and setting it to True again in
    #: response to any events that actually do require a window contents
    #: update.
    #:
    #: :type: bool
    #: :since: pyglet 1.1
    invalid = True

    # Instance variables accessible only via properties

    _width = None
    _height = None
    _caption = None
    _resizable = False
    _style = WINDOW_STYLE_DEFAULT
    _fullscreen = False
    _visible = False
    _vsync = False
    _screen = None
    _config = None
    _context = None

    # Used to restore window size and position after fullscreen
    _windowed_size = None
    _windowed_location = None

    # Subclasses should update these after relevant events
    _mouse_cursor = DefaultMouseCursor()
    _mouse_x = 0
    _mouse_y = 0
    _mouse_visible = True
    _mouse_exclusive = False
    _mouse_in_window = True

    _event_queue = None
    _enable_event_queue = True    # overridden by EventLoop.
    _allow_dispatch_event = False # controlled by dispatch_events stack frame

    # Class attributes

    _default_width = 640
    _default_height = 480

    def __init__(self, 
                 width=None,
                 height=None,
                 caption=None,
                 resizable=False,
                 style=WINDOW_STYLE_DEFAULT,
                 fullscreen=False,
                 visible=True,
                 vsync=True,
                 display=None,
                 screen=None,
                 config=None,
                 context=None):
        '''Create a window.

        All parameters are optional, and reasonable defaults are assumed
        where they are not specified.

        The `display`, `screen`, `config` and `context` parameters form
        a hierarchy of control: there is no need to specify more than 
        one of these.  For example, if you specify `screen` the `display`
        will be inferred, and a default `config` and `context` will be
        created.

        `config` is a special case; it can be a template created by the
        user specifying the attributes desired, or it can be a complete
        `config` as returned from `Screen.get_matching_configs` or similar.

        The context will be active as soon as the window is created, as if
        `switch_to` was just called.

        :Parameters:
            `width` : int
                Width of the window, in pixels.  Not valid if `fullscreen`
                is True.  Defaults to 640.
            `height` : int
                Height of the window, in pixels.  Not valid if `fullscreen`
                is True.  Defaults to 480.
            `caption` : str or unicode
                Initial caption (title) of the window.  Defaults to
                ``sys.argv[0]``.
            `resizable` : bool
                If True, the window will be resizable.  Defaults to False.
            `style` : int
                One of the ``WINDOW_STYLE_*`` constants specifying the
                border style of the window.
            `fullscreen` : bool
                If True, the window will cover the entire screen rather
                than floating.  Defaults to False.
            `visible` : bool
                Determines if the window is visible immediately after
                creation.  Defaults to True.  Set this to False if you
                would like to change attributes of the window before
                having it appear to the user.
            `vsync` : bool
                If True, buffer flips are synchronised to the primary screen's
                vertical retrace, eliminating flicker.
            `display` : `Display`
                The display device to use.  Useful only under X11.
            `screen` : `Screen`
                The screen to use, if in fullscreen.
            `config` : `pyglet.gl.Config`
                Either a template from which to create a complete config,
                or a complete config.
            `context` : `pyglet.gl.Context`
                The context to attach to this window.  The context must
                not already be attached to another window.

        '''
        EventDispatcher.__init__(self)
        self._event_queue = []

        if not display:
            display = get_platform().get_default_display()

        if not screen:
            screen = display.get_default_screen()

        if not config:
            for template_config in [
                gl.Config(double_buffer=True, depth_size=24),
                gl.Config(double_buffer=True, depth_size=16)]:
                try:
                    config = screen.get_best_config(template_config)
                    break
                except NoSuchConfigException:
                    pass
            if not config:
                raise NoSuchConfigException('No standard config is available.')

        if not config.is_complete():
            config = screen.get_best_config(config)

        if not context:
            context = config.create_context(gl.current_context)

        # Set these in reverse order to above, to ensure we get user
        # preference
        self._context = context
        self._config = self._context.config
        self._screen = self._config.screen
        self._display = self._screen.display

        if fullscreen:
            if width is not None or height is not None:
                raise WindowException(
                    'Width and height cannot be specified with fullscreen.')
            self._windowed_size = self._default_width, self._default_height
            width = self._screen.width
            height = self._screen.height
        else:
            if width is None:
                width = self._default_width
            if height is None:
                height = self._default_height

        self._width = width
        self._height = height
        self._resizable = resizable
        self._fullscreen = fullscreen
        self._style = style
        if pyglet.options['vsync'] is not None:
            self._vsync = pyglet.options['vsync']
        else:
            self._vsync = vsync


        if caption is None:
            caption = sys.argv[0]
        self._caption = caption

        from pyglet import app
        app.windows.add(self)
        self._create()

        self.switch_to()
        if visible:
            self.set_visible(True)
            self.activate()

    def _create(self):
        raise NotImplementedError('abstract')

    def _recreate(self, changes):
        '''Recreate the window with current attributes.

        :Parameters:
            `changes` : list of str
                List of attribute names that were changed since the last
                `_create` or `_recreate`.  For example, ``['fullscreen']``
                is given if the window is to be toggled to or from fullscreen. 
        '''
        raise NotImplementedError('abstract')

    def flip(self):
        '''Swap the OpenGL front and back buffers.

        Call this method on a double-buffered window to update the
        visible display with the back buffer.  The contents of the back buffer
        is undefined after this operation.

        Windows are double-buffered by default.  This method is called
        automatically by `EventLoop` after the `on_draw` event.
        '''
        raise NotImplementedError('abstract')

    def switch_to(self):
        '''Make this window the current OpenGL rendering context.

        Only one OpenGL context can be active at a time.  This method sets
        the current window's context to be current.  You should use this
        method in preference to `pyglet.gl.Context.set_current`, as it may
        perform additional initialisation functions.
        '''
        raise NotImplementedError('abstract')

    def set_fullscreen(self, fullscreen=True, screen=None):
        '''Toggle to or from fullscreen.

        After toggling fullscreen, the GL context should have retained its
        state and objects, however the buffers will need to be cleared and
        redrawn.

        :Parameters:
            `fullscreen` : bool
                True if the window should be made fullscreen, False if it
                should be windowed.
            `screen` : Screen
                If not None and fullscreen is True, the window is moved to the
                given screen.  The screen must belong to the same display as
                the window.

        '''
        if fullscreen == self._fullscreen and screen is None:
            return

        if not self._fullscreen:
            # Save windowed size
            self._windowed_size = self.get_size()
            self._windowed_location = self.get_location()

        if fullscreen and screen is not None:
            assert screen.display is self.display
            self._screen = screen

        self._fullscreen = fullscreen
        if self._fullscreen:
            self._width = self.screen.width
            self._height = self.screen.height
        else:
            self._width, self._height = self._windowed_size

        self._recreate(['fullscreen'])

        if not self._fullscreen and self._windowed_location:
            # Restore windowed location -- no effect on OS X because of
            # deferred recreate.  Move into platform _create? XXX
            self.set_location(*self._windowed_location)

    def on_resize(self, width, height):
        '''A default resize event handler.

        This default handler updates the GL viewport to cover the entire
        window and sets the ``GL_PROJECTION`` matrix to be orthagonal in
        window space.  The bottom-left corner is (0, 0) and the top-right
        corner is the width and height of the window in pixels.

        Override this event handler with your own to create another
        projection, for example in perspective.
        '''
        gl.glViewport(0, 0, width, height)
        gl.glMatrixMode(gl.GL_PROJECTION)
        gl.glLoadIdentity()
        gl.glOrtho(0, width, 0, height, -1, 1)
        gl.glMatrixMode(gl.GL_MODELVIEW)

    def on_close(self):
        '''Default on_close handler.'''
        self.has_exit = True
        from pyglet import app
        if app.event_loop is not None:
            self.close()

    def on_key_press(self, symbol, modifiers):
        '''Default on_key_press handler.'''
        if symbol == key.ESCAPE:
            self.dispatch_event('on_close')

    def close(self):
        '''Close the window.

        After closing the window, the GL context will be invalid.  The
        window instance cannot be reused once closed (see also `set_visible`).

        The `pyglet.app.EventLoop.on_window_close` event is dispatched on
        `pyglet.app.event_loop` when this method is called.
        '''
        from pyglet import app
        if not self._context:
            return
        app.windows.remove(self)
        self._context.destroy()
        self._config = None
        self._context = None
        if app.event_loop:
            app.event_loop.dispatch_event('on_window_close', self)

    def draw_mouse_cursor(self):
        '''Draw the custom mouse cursor.

        If the current mouse cursor has ``drawable`` set, this method
        is called before the buffers are flipped to render it.  
        
        This method always leaves the ``GL_MODELVIEW`` matrix as current,
        regardless of what it was set to previously.  No other GL state
        is affected.

        There is little need to override this method; instead, subclass
        ``MouseCursor`` and provide your own ``draw`` method.
        '''
        # Draw mouse cursor if set and visible.
        # XXX leaves state in modelview regardless of starting state
        if (self._mouse_cursor.drawable and 
            self._mouse_visible and 
            self._mouse_in_window):
            gl.glMatrixMode(gl.GL_PROJECTION)
            gl.glPushMatrix()
            gl.glLoadIdentity()
            gl.glOrtho(0, self.width, 0, self.height, -1, 1)

            gl.glMatrixMode(gl.GL_MODELVIEW)
            gl.glPushMatrix()
            gl.glLoadIdentity()

            self._mouse_cursor.draw(self._mouse_x, self._mouse_y)

            gl.glMatrixMode(gl.GL_PROJECTION)
            gl.glPopMatrix()

            gl.glMatrixMode(gl.GL_MODELVIEW)
            gl.glPopMatrix()

    # Properties provide read-only access to instance variables.  Use
    # set_* methods to change them if applicable.

    caption = property(lambda self: self._caption,
        doc='''The window caption (title).  Read-only.

        :type: str
        ''')
    resizable = property(lambda self: self._resizable,
        doc='''True if the window is resizeable.  Read-only.

        :type: bool
        ''')
    style = property(lambda self: self._style,
        doc='''The window style; one of the ``WINDOW_STYLE_*`` constants.
        Read-only.
        
        :type: int
        ''')
    fullscreen = property(lambda self: self._fullscreen,
        doc='''True if the window is currently fullscreen.  Read-only.
        
        :type: bool
        ''')
    visible = property(lambda self: self._visible,
        doc='''True if the window is currently visible.  Read-only.
        
        :type: bool
        ''')
    vsync = property(lambda self: self._vsync,
        doc='''True if buffer flips are synchronised to the screen's vertical
        retrace.  Read-only.
        
        :type: bool
        ''')
    display = property(lambda self: self._display,
        doc='''The display this window belongs to.  Read-only.

        :type: `Display`
        ''')
    screen = property(lambda self: self._screen,
        doc='''The screen this window is fullscreen in.  Read-only.
        
        :type: `Screen`
        ''')
    config = property(lambda self: self._config,
        doc='''A GL config describing the context of this window.  Read-only.
        
        :type: `pyglet.gl.Config`
        ''')
    context = property(lambda self: self._context,
        doc='''The OpenGL context attached to this window.  Read-only.
        
        :type: `pyglet.gl.Context`
        ''')

    # These are the only properties that can be set
    width = property(lambda self: self.get_size()[0],
                     lambda self, width: self.set_size(width, self.height),
         doc='''The width of the window, in pixels.  Read-write.
         
         :type: int
         ''')

    height = property(lambda self: self.get_size()[1],
                      lambda self, height: self.set_size(self.width, height),
         doc='''The height of the window, in pixels.  Read-write.
         
         :type: int
         ''')

    def set_caption(self, caption):
        '''Set the window's caption.

        The caption appears in the titlebar of the window, if it has one,
        and in the taskbar on Windows and many X11 window managers.

        :Parameters:
            `caption` : str or unicode
                The caption to set.

        '''
        raise NotImplementedError('abstract')

    def set_minimum_size(self, width, height):
        '''Set the minimum size of the window.

        Once set, the user will not be able to resize the window smaller
        than the given dimensions.  There is no way to remove the
        minimum size constraint on a window (but you could set it to 0,0).

        The behaviour is undefined if the minimum size is set larger than
        the current size of the window.

        The window size does not include the border or title bar.

        :Parameters:
            `width` : int
                Minimum width of the window, in pixels.
            `height` : int
                Minimum height of the window, in pixels.

        '''
        raise NotImplementedError('abstract')

    def set_maximum_size(self, width, height):
        '''Set the maximum size of the window.

        Once set, the user will not be able to resize the window larger
        than the given dimensions.  There is no way to remove the
        maximum size constraint on a window (but you could set it to a large
        value).

        The behaviour is undefined if the maximum size is set smaller than
        the current size of the window.

        The window size does not include the border or title bar.

        :Parameters:
            `width` : int
                Maximum width of the window, in pixels.
            `height` : int
                Maximum height of the window, in pixels.

        '''
        raise NotImplementedError('abstract')

    def set_size(self, width, height):
        '''Resize the window.
        
        The behaviour is undefined if the window is not resizable, or if
        it is currently fullscreen.

        The window size does not include the border or title bar.

        :Parameters:
            `width` : int
                New width of the window, in pixels.
            `height` : int
                New height of the window, in pixels.

        '''
        raise NotImplementedError('abstract')

    def get_size(self):
        '''Return the current size of the window.

        The window size does not include the border or title bar.

        :rtype: (int, int)
        :return: The width and height of the window, in pixels.
        '''
        raise NotImplementedError('abstract')

    def set_location(self, x, y):
        '''Set the position of the window.

        :Parameters:
            `x` : int
                Distance of the left edge of the window from the left edge
                of the virtual desktop, in pixels.
            `y` : int
                Distance of the top edge of the window from the top edge of
                the virtual desktop, in pixels.

        '''
        raise NotImplementedError('abstract')

    def get_location(self):
        '''Return the current position of the window.

        :rtype: (int, int)
        :return: The distances of the left and top edges from their respective
            edges on the virtual desktop, in pixels.
        '''
        raise NotImplementedError('abstract')

    def activate(self):
        '''Attempt to restore keyboard focus to the window.

        Depending on the window manager or operating system, this may not
        be successful.  For example, on Windows XP an application is not
        allowed to "steal" focus from another application.  Instead, the
        window's taskbar icon will flash, indicating it requires attention.
        '''
        raise NotImplementedError('abstract')

    def set_visible(self, visible=True):    
        '''Show or hide the window.

        :Parameters:
            `visible` : bool
                If True, the window will be shown; otherwise it will be
                hidden.

        '''
        raise NotImplementedError('abstract')

    def minimize(self):
        '''Minimize the window.
        '''
        raise NotImplementedError('abstract')

    def maximize(self):
        '''Maximize the window.

        The behaviour of this method is somewhat dependent on the user's
        display setup.  On a multi-monitor system, the window may maximize
        to either a single screen or the entire virtual desktop.
        '''
        raise NotImplementedError('abstract')

    def set_vsync(self, vsync):
        '''Enable or disable vertical sync control.

        When enabled, this option ensures flips from the back to the front
        buffer are performed only during the vertical retrace period of the
        primary display.  This can prevent "tearing" or flickering when
        the buffer is updated in the middle of a video scan.

        Note that LCD monitors have an analagous time in which they are not
        reading from the video buffer; while it does not correspond to
        a vertical retrace it has the same effect.

        With multi-monitor systems the secondary monitor cannot be
        synchronised to, so tearing and flicker cannot be avoided when the
        window is positioned outside of the primary display.  In this case
        it may be advisable to forcibly reduce the framerate (for example,
        using `pyglet.clock.set_fps_limit`).

        :Parameters:
            `vsync` : bool
                If True, vsync is enabled, otherwise it is disabled.

        '''
        raise NotImplementedError('abstract')

    def set_mouse_visible(self, visible=True):
        '''Show or hide the mouse cursor.

        The mouse cursor will only be hidden while it is positioned within
        this window.  Mouse events will still be processed as usual.

        :Parameters:
            `visible` : bool
                If True, the mouse cursor will be visible, otherwise it
                will be hidden.

        '''
        self._mouse_visible = visible
        self.set_mouse_platform_visible()

    def set_mouse_platform_visible(self, platform_visible=None):
        '''Set the platform-drawn mouse cursor visibility.  This is called
        automatically after changing the mouse cursor or exclusive mode.

        Applications should not normally need to call this method, see
        `set_mouse_visible` instead.

        :Parameters:
            `platform_visible` : bool or None
                If None, sets platform visibility to the required visibility
                for the current exclusive mode and cursor type.  Otherwise,
                a bool value will override and force a visibility.

        '''
        raise NotImplementedError()

    def set_mouse_cursor(self, cursor=None):
        '''Change the appearance of the mouse cursor.

        The appearance of the mouse cursor is only changed while it is
        within this window.

        :Parameters:
            `cursor` : `MouseCursor`
                The cursor to set, or None to restore the default cursor.

        '''
        if cursor is None:
            cursor = DefaultMouseCursor()
        self._mouse_cursor = cursor
        self.set_mouse_platform_visible()

    def set_exclusive_mouse(self, exclusive=True):
        '''Hide the mouse cursor and direct all mouse events to this
        window.

        When enabled, this feature prevents the mouse leaving the window.  It
        is useful for certain styles of games that require complete control of
        the mouse.  The position of the mouse as reported in subsequent events
        is meaningless when exclusive mouse is enabled; you should only use
        the relative motion parameters ``dx`` and ``dy``.

        :Parameters:
            `exclusive` : bool
                If True, exclusive mouse is enabled, otherwise it is disabled.

        '''
        raise NotImplementedError('abstract')

    def set_exclusive_keyboard(self, exclusive=True):
        '''Prevent the user from switching away from this window using
        keyboard accelerators.

        When enabled, this feature disables certain operating-system specific
        key combinations such as Alt+Tab (Command+Tab on OS X).  This can be
        useful in certain kiosk applications, it should be avoided in general
        applications or games.

        :Parameters:
            `exclusive` : bool
                If True, exclusive keyboard is enabled, otherwise it is
                disabled.

        '''
        raise NotImplementedError('abstract')

    def get_system_mouse_cursor(self, name):
        '''Obtain a system mouse cursor.

        Use `set_mouse_cursor` to make the cursor returned by this method
        active.  The names accepted by this method are the ``CURSOR_*``
        constants defined on this class.

        :Parameters:
            `name` : str
                Name describing the mouse cursor to return.  For example,
                ``CURSOR_WAIT``, ``CURSOR_HELP``, etc.

        :rtype: `MouseCursor`
        :return: A mouse cursor which can be used with `set_mouse_cursor`.
        '''
        raise NotImplementedError()

    def set_icon(self, *images):
        '''Set the window icon.

        If multiple images are provided, one with an appropriate size 
        will be selected (if the correct size is not provided, the image
        will be scaled).

        Useful sizes to provide are 16x16, 32x32, 64x64 (Mac only) and
        128x128 (Mac only).

        :Parameters:
            `images` : sequence of `pyglet.image.AbstractImage`
                List of images to use for the window icon.
        
        '''
        pass

    def clear(self):
        '''Clear the window.

        This is a convenience method for clearing the color and depth
        buffer.  The window must be the active context (see `switch_to`).
        '''
        gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
    
    def dispatch_event(self, *args):
        if not self._enable_event_queue or self._allow_dispatch_event:
            EventDispatcher.dispatch_event(self, *args)
        else:
            self._event_queue.append(args)

    def dispatch_events(self):
        '''Poll the operating system event queue for new events and call
        attached event handlers.

        This method is provided for legacy applications targeting pyglet 1.0,
        and advanced applications that must integrate their event loop
        into another framework.

        Typical applications should use `pyglet.app.run`.
        '''
        raise NotImplementedError('abstract')

    # If documenting, show the event methods.  Otherwise, leave them out
    # as they are not really methods.
    if hasattr(sys, 'is_epydoc') and sys.is_epydoc:
        def on_key_press(symbol, modifiers):
            '''A key on the keyboard was pressed (and held down).

            In pyglet 1.0 the default handler sets `has_exit` to ``True`` if
            the ``ESC`` key is pressed.

            In pyglet 1.1 the default handler dispatches the `on_close`
            event if the ``ESC`` key is pressed.

            :Parameters:
                `symbol` : int
                    The key symbol pressed.
                `modifiers` : int
                    Bitwise combination of the key modifiers active.
            
            :event:
            '''

        def on_key_release(symbol, modifiers):
            '''A key on the keyboard was released.

            :Parameters:
                `symbol` : int
                    The key symbol pressed.
                `modifiers` : int
                    Bitwise combination of the key modifiers active.

            :event:
            '''

        def on_text(text):
            '''The user input some text.

            Typically this is called after `on_key_press` and before
            `on_key_release`, but may also be called multiple times if the key
            is held down (key repeating); or called without key presses if
            another input method was used (e.g., a pen input).

            You should always use this method for interpreting text, as the
            key symbols often have complex mappings to their unicode
            representation which this event takes care of.

            :Parameters:
                `text` : unicode
                    The text entered by the user.

            :event:
            '''

        def on_text_motion(motion):
            '''The user moved the text input cursor.

            Typically this is called after `on_key_press` and before
            `on_key_release`, but may also be called multiple times if the key
            is help down (key repeating).

            You should always use this method for moving the text input cursor
            (caret), as different platforms have different default keyboard
            mappings, and key repeats are handled correctly.

            The values that `motion` can take are defined in
            `pyglet.window.key`:

            * MOTION_UP
            * MOTION_RIGHT
            * MOTION_DOWN
            * MOTION_LEFT
            * MOTION_NEXT_WORD
            * MOTION_PREVIOUS_WORD
            * MOTION_BEGINNING_OF_LINE
            * MOTION_END_OF_LINE
            * MOTION_NEXT_PAGE
            * MOTION_PREVIOUS_PAGE
            * MOTION_BEGINNING_OF_FILE
            * MOTION_END_OF_FILE
            * MOTION_BACKSPACE
            * MOTION_DELETE

            :Parameters:
                `motion` : int
                    The direction of motion; see remarks.

            :event:
            '''

        def on_text_motion_select(motion):
            '''The user moved the text input cursor while extending the
            selection.

            Typically this is called after `on_key_press` and before
            `on_key_release`, but may also be called multiple times if the key
            is help down (key repeating).

            You should always use this method for responding to text selection
            events rather than the raw `on_key_press`, as different platforms
            have different default keyboard mappings, and key repeats are
            handled correctly.

            The values that `motion` can take are defined in `pyglet.window.key`:

            * MOTION_UP
            * MOTION_RIGHT
            * MOTION_DOWN
            * MOTION_LEFT
            * MOTION_NEXT_WORD
            * MOTION_PREVIOUS_WORD
            * MOTION_BEGINNING_OF_LINE
            * MOTION_END_OF_LINE
            * MOTION_NEXT_PAGE
            * MOTION_PREVIOUS_PAGE
            * MOTION_BEGINNING_OF_FILE
            * MOTION_END_OF_FILE

            :Parameters:
                `motion` : int
                    The direction of selection motion; see remarks.

            :event:
            '''

        def on_mouse_motion(x, y, dx, dy):
            '''The mouse was moved with no buttons held down.

            :Parameters:
                `x` : int
                    Distance in pixels from the left edge of the window.
                `y` : int
                    Distance in pixels from the bottom edge of the window.
                `dx` : int
                    Relative X position from the previous mouse position.
                `dy` : int
                    Relative Y position from the previous mouse position.

            :event:
            '''

        def on_mouse_drag(x, y, dx, dy, buttons, modifiers):
            '''The mouse was moved with one or more mouse buttons pressed.

            This event will continue to be fired even if the mouse leaves
            the window, so long as the drag buttons are continuously held down.

            :Parameters:
                `x` : int
                    Distance in pixels from the left edge of the window.
                `y` : int
                    Distance in pixels from the bottom edge of the window.
                `dx` : int
                    Relative X position from the previous mouse position.
                `dy` : int
                    Relative Y position from the previous mouse position.
                `buttons` : int
                    Bitwise combination of the mouse buttons currently pressed.
                `modifiers` : int
                    Bitwise combination of any keyboard modifiers currently
                    active.

            :event:
            '''

        def on_mouse_press(x, y, button, modifiers):
            '''A mouse button was pressed (and held down).

            :Parameters:
                `x` : int
                    Distance in pixels from the left edge of the window.
                `y` : int
                    Distance in pixels from the bottom edge of the window.
                `button` : int
                    The mouse button that was pressed.
                `modifiers` : int
                    Bitwise combination of any keyboard modifiers currently
                    active.
                
            :event:
            '''

        def on_mouse_release(x, y, button, modifiers):
            '''A mouse button was released.

            :Parameters:
                `x` : int
                    Distance in pixels from the left edge of the window.
                `y` : int
                    Distance in pixels from the bottom edge of the window.
                `button` : int
                    The mouse button that was released.
                `modifiers` : int
                    Bitwise combination of any keyboard modifiers currently
                    active.

            :event:
            '''
                
        def on_mouse_scroll(x, y, scroll_x, scroll_y):
            '''The mouse wheel was scrolled.

            Note that most mice have only a vertical scroll wheel, so
            `scroll_x` is usually 0.  An exception to this is the Apple Mighty
            Mouse, which has a mouse ball in place of the wheel which allows
            both `scroll_x` and `scroll_y` movement.

            :Parameters:
                `x` : int
                    Distance in pixels from the left edge of the window.
                `y` : int
                    Distance in pixels from the bottom edge of the window.
                `scroll_x` : int
                    Number of "clicks" towards the right (left if negative).
                `scroll_y` : int
                    Number of "clicks" upwards (downwards if negative).

            :event:
            '''

        def on_close():
            '''The user attempted to close the window.

            This event can be triggered by clicking on the "X" control box in
            the window title bar, or by some other platform-dependent manner.

            The default handler sets `has_exit` to ``True``.  In pyglet 1.1, if
            `pyglet.app.event_loop` is being used, `close` is also called,
            closing the window immediately.

            :event:
            '''

        def on_mouse_enter(x, y):
            '''The mouse was moved into the window.

            This event will not be trigged if the mouse is currently being
            dragged.

            :Parameters:
                `x` : int
                    Distance in pixels from the left edge of the window.
                `y` : int
                    Distance in pixels from the bottom edge of the window.

            :event:
            '''

        def on_mouse_leave(x, y):
            '''The mouse was moved outside of the window.

            This event will not be trigged if the mouse is currently being
            dragged.  Note that the coordinates of the mouse pointer will be
            outside of the window rectangle.

            :Parameters:
                `x` : int
                    Distance in pixels from the left edge of the window.
                `y` : int
                    Distance in pixels from the bottom edge of the window.

            :event:
            '''

        def on_expose():
            '''A portion of the window needs to be redrawn.

            This event is triggered when the window first appears, and any time
            the contents of the window is invalidated due to another window
            obscuring it.

            There is no way to determine which portion of the window needs
            redrawing.  Note that the use of this method is becoming
            increasingly uncommon, as newer window managers composite windows
            automatically and keep a backing store of the window contents.

            :event:
            '''

        def on_resize(width, height):
            '''The window was resized.

            The window will have the GL context when this event is dispatched;
            there is no need to call `switch_to` in this handler.

            :Parameters:
                `width` : int
                    The new width of the window, in pixels.
                `height` : int
                    The new height of the window, in pixels.

            :event:
            '''

        def on_move(x, y):
            '''The window was moved.

            :Parameters:
                `x` : int
                    Distance from the left edge of the screen to the left edge
                    of the window.
                `y` : int
                    Distance from the top edge of the screen to the top edge of
                    the window.  Note that this is one of few methods in pyglet
                    which use a Y-down coordinate system.

            :event:
            '''

        def on_activate():
            '''The window was activated.

            This event can be triggered by clicking on the title bar, bringing
            it to the foreground; or by some platform-specific method.

            When a window is "active" it has the keyboard focus.

            :event:
            '''

        def on_deactivate():
            '''The window was deactivated.

            This event can be triggered by clicking on another application
            window.  When a window is deactivated it no longer has the
            keyboard focus.

            :event:
            '''

        def on_show():
            '''The window was shown.

            This event is triggered when a window is restored after being
            minimised, or after being displayed for the first time.

            :event:
            '''

        def on_hide():
            '''The window was hidden.

            This event is triggered when a window is minimised or (on Mac OS X)
            hidden by the user.

            :event:
            '''

        def on_context_lost():
            '''The window's GL context was lost.
            
            When the context is lost no more GL methods can be called until it
            is recreated.  This is a rare event, triggered perhaps by the user
            switching to an incompatible video mode.  When it occurs, an
            application will need to reload all objects (display lists, texture
            objects, shaders) as well as restore the GL state.

            :event:
            '''

        def on_context_state_lost():
            '''The state of the window's GL context was lost.

            pyglet may sometimes need to recreate the window's GL context if
            the window is moved to another video device, or between fullscreen
            or windowed mode.  In this case it will try to share the objects
            (display lists, texture objects, shaders) between the old and new
            contexts.  If this is possible, only the current state of the GL
            context is lost, and the application should simply restore state.

            :event:
            '''

        def on_draw():
            '''The window contents must be redrawn.

            The `EventLoop` will dispatch this event when the window
            should be redrawn.  This will happen during idle time after
            any window events and after any scheduled functions were called.

            The window will already have the GL context, so there is no
            need to call `switch_to`.  The window's `flip` method will
            be called after this event, so your event handler should not.

            You should make no assumptions about the window contents when
            this event is triggered; a resize or expose event may have
            invalidated the framebuffer since the last time it was drawn.

            :since: pyglet 1.1

            :event:
            '''

BaseWindow.register_event_type('on_key_press')
BaseWindow.register_event_type('on_key_release')
BaseWindow.register_event_type('on_text')
BaseWindow.register_event_type('on_text_motion')
BaseWindow.register_event_type('on_text_motion_select')
BaseWindow.register_event_type('on_mouse_motion')
BaseWindow.register_event_type('on_mouse_drag')
BaseWindow.register_event_type('on_mouse_press')
BaseWindow.register_event_type('on_mouse_release')
BaseWindow.register_event_type('on_mouse_scroll')
BaseWindow.register_event_type('on_mouse_enter')
BaseWindow.register_event_type('on_mouse_leave')
BaseWindow.register_event_type('on_close')
BaseWindow.register_event_type('on_expose')
BaseWindow.register_event_type('on_resize')
BaseWindow.register_event_type('on_move')
BaseWindow.register_event_type('on_activate')
BaseWindow.register_event_type('on_deactivate')
BaseWindow.register_event_type('on_show')
BaseWindow.register_event_type('on_hide')
BaseWindow.register_event_type('on_context_lost')
BaseWindow.register_event_type('on_context_state_lost')
BaseWindow.register_event_type('on_draw')

def get_platform():
    '''Get an instance of the Platform most appropriate for this
    system.

    :rtype: `Platform`
    :return: The platform instance.
    '''
    return _platform

_is_epydoc = hasattr(sys, 'is_epydoc') and sys.is_epydoc

if _is_epydoc:
    # We are building documentation
    Window = BaseWindow
    Window.__name__ = 'Window'
    del BaseWindow
else:
    # Try to determine which platform to use.
    if sys.platform == 'darwin':
        from pyglet.window.carbon import CarbonPlatform, CarbonWindow
        _platform = CarbonPlatform()
        Window = CarbonWindow
    elif sys.platform in ('win32', 'cygwin'):
        from pyglet.window.win32 import Win32Platform, Win32Window
        _platform = Win32Platform()
        Window = Win32Window
    else:
        from pyglet.window.xlib import XlibPlatform, XlibWindow
        _platform = XlibPlatform()
        Window = XlibWindow

# Create shadow window. (trickery is for circular import)
if not _is_epydoc:
    pyglet.window = sys.modules[__name__]
    gl._create_shadow_window()
