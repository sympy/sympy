#!/usr/bin/env python
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

Within your main run loop, you must call `Window.dispatch_events` regularly.
Windows are double-buffered by default, so you must call `Window.flip` to
update the display::

    while not win.has_exit:
        win.dispatch_events()
        # ... drawing commands ...
        win.flip()

Creating a game window
----------------------

Use `Window.set_exclusive_mouse` to hide the mouse cursor and receive relative
mouse movement events.  Specify ``fullscreen=True`` as a keyword argument to
the `Window` constructor to render to the entire screen rather than opening a
window::

    win = Window(fullscreen=True)
    win.set_mouse_exclusive()

Working with multiple windows
-----------------------------

You can open any number of windows and render to them individually.  Each
window must have the event handlers set on it that you are interested in
(i.e., each window will have its own mouse event handler).  

You must call `Window.dispatch_events` for each window.  Before rendering
to a window, you must call `Window.switch_to` to set the active GL context.
Here is an example run loop for a list of windows::

    windows = # list of Window instances
    while windows:
        for win in windows:
            win.dispatch_events()
            if win.has_exit:
                win.close()
        windows = [w for w in windows if not w.has_exit]
   
        for win in windows:
            win.switch_to()
            # ... drawing commands for this window ...
            win.flip()

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
__version__ = '$Id: __init__.py 1036 2007-07-14 07:08:07Z Alex.Holkner $'

import pprint
import sys

from pyglet import gl
from pyglet.gl import gl_info
from pyglet.window.event import WindowEventDispatcher, WindowExitHandler
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

    def __init__(self, image, hot_x, hot_y):
        '''Create a mouse cursor from an image.

        :Parameters:
            `image` : `pyglet.image.AbstractImage`
                Image to use for the mouse cursor.  It must have a
                valid `texture` attribute.
            `hot_x` : int
                X coordinate of the "hot" spot in the image.
            `hot_y` : int
                Y coordinate of the "hot" spot in the image, measured
                from the bottom.
        '''
        self.texture = image.texture
        self.hot_x = hot_x
        self.hot_y = hot_y

    def draw(self, x, y):
        gl.glPushAttrib(gl.GL_ENABLE_BIT)
        gl.glEnable(gl.GL_BLEND)
        gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
        self.texture.blit(x - self.hot_x, y - self.hot_y, 0)
        gl.glPopAttrib()

class BaseWindow(WindowEventDispatcher, WindowExitHandler):
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
    '''
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

    def __init__(self, 
                 width=640,
                 height=480,
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
                Width of the window, in pixels.  Ignored if `fullscreen`
                is True.  Defaults to 640.
            `height` : int
                Height of the window, in pixels.  Ignored if `fullscreen`
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
        WindowEventDispatcher.__init__(self)

        if not display:
            display = get_platform().get_default_display()

        if not screen:
            screen = display.get_default_screen()

        if not config:
            config = gl.Config(double_buffer=True,
                               depth_size=24)

        if not config.is_complete():
            config = screen.get_best_config(config)

        if not context:
            context = config.create_context(gl.get_current_context())

        if fullscreen:
            self._windowed_size = width, height
            width = screen.width
            height = screen.height

        self._width = width
        self._height = height
        self._resizable = resizable
        self._fullscreen = fullscreen
        self._style = style
        self._vsync = vsync

        # Set these in reverse order to above, to ensure we get user
        # preference
        self._context = context
        self._config = self._context.config
        self._screen = self._config.screen
        self._display = self._screen.display

        if caption is None:
            caption = sys.argv[0]
        self._caption = caption

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

    def set_fullscreen(self, fullscreen=True):
        '''Toggle to or from fullscreen.

        After toggling fullscreen, the GL context should have retained its
        state and objects, however the buffers will need to be cleared and
        redrawn.

        :Parameters:
            `fullscreen` : bool
                True if the window should be made fullscreen, False if it
                should be windowed.

        '''
        if fullscreen == self._fullscreen:
            return

        if not self._fullscreen:
            # Save windowed size
            self._windowed_size = self.get_size()
            self._windowed_location = self.get_location()

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
        self.switch_to()
        gl.glViewport(0, 0, width, height)
        gl.glMatrixMode(gl.GL_PROJECTION)
        gl.glLoadIdentity()
        gl.glOrtho(0, width, 0, height, -1, 1)
        gl.glMatrixMode(gl.GL_MODELVIEW)

    def close(self):
        '''Close the window.

        Windows are closed automatically when the process exits, so this
        method need only be called when multiple windows or console input
        are being used.

        After closing the window, the GL context will be invalid.  The
        window instance cannot be reused once closed (see also `set_visible`).
        '''
        self._context.destroy()
        self._config = None
        self._context = None

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

    def dispatch_events(self):
        '''Process the operating system event queue and call attached
        event handlers.
        '''
        raise NotImplementedError('abstract')

def get_platform():
    '''Get an instance of the Platform most appropriate for this
    system.

    :rtype: `Platform`
    :return: The platform instance.
    '''
    return _platform

if hasattr(sys, 'is_epydoc') and sys.is_epydoc:
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

