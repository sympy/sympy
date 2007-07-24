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

'''Events for `pyglet.window`.

See `WindowEventHandler` for a description of the window event types.
'''

__docformat__ = 'restructuredtext'
__version__ = '$Id: event.py 979 2007-06-29 11:04:36Z Alex.Holkner $'

import sys

from pyglet.window import key
from pyglet.window import mouse
from pyglet.event import EventDispatcher

# If documenting, document the event routines, otherwise the class is empty.

if hasattr(sys, 'is_epydoc') and sys.is_epydoc:
    class WindowEventDispatcher(EventDispatcher):
        '''The abstract window event dispatcher.
        '''

        def on_key_press(symbol, modifiers):
            '''A key on the keyboard was pressed (and held down).

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
            `on_key_release`, but may also be called multiple times if the key is
            held down (key repeating); or called without key presses if another
            input method was used (e.g., a pen input).

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
            `on_key_release`, but may also be called multiple times if the key is
            help down (key repeating).

            You should always use this method for moving the text input cursor
            (caret), as different platforms have different default keyboard
            mappings, and key repeats are handled correctly.

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
            * MOTION_BACKSPACE
            * MOTION_DELETE

            :Parameters:
                `motion` : int
                    The direction of motion; see remarks.

            :event:
            '''

        def on_text_motion_select(motion):
            '''The user moved the text input cursor while extending the selection.

            Typically this is called after `on_key_press` and before
            `on_key_release`, but may also be called multiple times if the key is
            help down (key repeating).

            You should always use this method for responding to text selection
            events rather than the raw `on_key_press`, as different platforms have
            different default keyboard mappings, and key repeats are handled
            correctly.

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
                `x` : float
                    Distance in pixels from the left edge of the window.
                `y` : float
                    Distance in pixels from the bottom edge of the window.
                `dx` : float
                    Relative X position from the previous mouse position.
                `dy` : float
                    Relative Y position from the previous mouse position.

            :event:
            '''

        def on_mouse_drag(x, y, dx, dy, buttons, modifiers):
            '''The mouse was moved with one or more mouse buttons pressed.

            This event will continue to be fired even if the mouse leaves
            the window, so long as the drag buttons are continuously held down.

            :Parameters:
                `x` : float
                    Distance in pixels from the left edge of the window.
                `y` : float
                    Distance in pixels from the bottom edge of the window.
                `dx` : float
                    Relative X position from the previous mouse position.
                `dy` : float
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
                `x` : float
                    Distance in pixels from the left edge of the window.
                `y` : float
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
                `x` : float
                    Distance in pixels from the left edge of the window.
                `y` : float
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
                `x` : float
                    Distance in pixels from the left edge of the window.
                `y` : float
                    Distance in pixels from the bottom edge of the window.
                `scroll_x` : int
                    Number of "clicks" towards the right (left if negative).
                `scroll_y` : int
                    Number of "clicks" upwards (downards if negative).

            :event:
            '''

        def on_close():
            '''The user attempted to close the window.

            This event can be triggered by clicking on the "X" control box in
            the window title bar, or by some other platform-dependent manner.

            :event:
            '''

        def on_mouse_enter(x, y):
            '''The mouse was moved into the window.

            This event will not be trigged if the mouse is currently being
            dragged.

            :Parameters:
                `x` : float
                    Distance in pixels from the left edge of the window.
                `y` : float
                    Distance in pixels from the bottom edge of the window.

            :event:
            '''

        def on_mouse_leave(x, y):
            '''The mouse was moved outside of the window.

            This event will not be trigged if the mouse is currently being
            dragged.  Note that the coordinates of the mouse pointer will be
            outside of the window rectangle.

            :Parameters:
                `x` : float
                    Distance in pixels from the left edge of the window.
                `y` : float
                    Distance in pixels from the bottom edge of the window.

            :event:
            '''

        def on_expose():
            '''A portion of the window needs to be redrawn.

            This event is triggered when the window first appears, and any time
            the contents of the window is invalidated due to another window
            obscuring it.

            There is no way to determine which portion of the window needs
            redrawing.  Note that the use of this method is becoming increasingly
            uncommon, as newer window managers composite windows automatically and
            keep a backing store of the window contents.

            :event:
            '''

        def on_resize(width, height):
            '''The window was resized.

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

            This event can be triggered by clicking on another application window.
            When a window is deactivated it no longer has the keyboard focus.

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
            
            When the context is lost no more GL methods can be called until it is
            recreated.  This is a rare event, triggered perhaps by the user
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

else:
    class WindowEventDispatcher(EventDispatcher):
        pass

class WindowExitHandler(object):
    '''Determine if the window should be closed.

    This event handler watches for the ESC key or the window close event
    and sets `self.has_exit` to True when either is pressed.  An instance
    of this class is automatically attached to all new `pyglet.window.Window`
    objects.

    :Ivariables:
        `has_exit` : bool
            True if the user wants to close the window.

    '''
    has_exit = False

    def on_close(self):
        self.has_exit = True

    def on_key_press(self, symbol, modifiers):
        if symbol == key.ESCAPE:
            self.has_exit = True

class WindowEventLogger(object):
    '''Print all events to a file.

    When this event handler is added to a window it prints out all events
    and their parameters; useful for debugging or discovering which events
    you need to handle.

    Example::

        win = window.Window()
        win.push_handlers(WindowEventLogger())

    '''
    def __init__(self, logfile=None):
        '''Create a `WindowEventLogger` which writes to `logfile`.

        :Parameters:
            `logfile` : file-like object
                The file to write to.  If unspecified, stdout will be used.

        '''
        if logfile is None:
            import sys
            logfile = sys.stdout
        self.file = logfile

    def on_key_press(self, symbol, modifiers):
        print >> self.file, 'on_key_press(symbol=%s, modifiers=%s)' % (
            key.symbol_string(symbol), key.modifiers_string(modifiers))

    def on_key_release(self, symbol, modifiers):
        print >> self.file, 'on_key_release(symbol=%s, modifiers=%s)' % (
            key.symbol_string(symbol), key.modifiers_string(modifiers))

    def on_text(self, text):
        print >> self.file, 'on_text(text=%r)' % text

    def on_text_motion(self, motion):
        print >> self.file, 'on_text_motion(motion=%s)' % (
            key.motion_string(motion))

    def on_text_motion_select(self, motion):
        print >> self.file, 'on_text_motion_select(motion=%s)' % (
            key.motion_string(motion))

    def on_mouse_motion(self, x, y, dx, dy):
        print >> self.file, 'on_mouse_motion(x=%d, y=%d, dx=%d, dy=%d)' % (
            x, y, dx, dy)

    def on_mouse_drag(self, x, y, dx, dy, buttons, modifiers):
        print >> self.file, 'on_mouse_drag(x=%d, y=%d, dx=%d, dy=%d, '\
                            'buttons=%s, modifiers=%s)' % (
              x, y, dx, dy, 
              mouse.buttons_string(buttons), key.modifiers_string(modifiers))

    def on_mouse_press(self, x, y, button, modifiers):
        print >> self.file, 'on_mouse_press(x=%d, y=%d, button=%r, '\
                            'modifiers=%s)' % (x, y,
            mouse.buttons_string(button), key.modifiers_string(modifiers))

    def on_mouse_release(self, x, y, button, modifiers):
        print >> self.file, 'on_mouse_release(x=%d, y=%d, button=%r, '\
                            'modifiers=%s)' % (x, y, 
            mouse.buttons_string(button), key.modifiers_string(modifiers))

    def on_mouse_scroll(self, x, y, dx, dy):
        print >> self.file, 'on_mouse_scroll(x=%f, y=%f, dx=%f, dy=%f)' % (
            x, y, dx, dy)

    def on_close(self):
        print >> self.file, 'on_close()'

    def on_mouse_enter(self, x, y):
        print >> self.file, 'on_mouse_enter(x=%d, y=%d)' % (x, y)

    def on_mouse_leave(self, x, y):
        print >> self.file, 'on_mouse_leave(x=%d, y=%d)' % (x, y)

    def on_expose(self):
        print >> self.file, 'on_expose()'

    def on_resize(self, width, height):
        print >> self.file, 'on_resize(width=%d, height=%d)' % (width, height)

    def on_move(self, x, y):
        print >> self.file, 'on_move(x=%d, y=%d)' % (x, y)

    def on_activate(self):
        print >> self.file, 'on_activate()'

    def on_deactivate(self):
        print >> self.file, 'on_deactivate()'

    def on_show(self):
        print >> self.file, 'on_show()'

    def on_hide(self):
        print >> self.file, 'on_hide()'

    def on_context_lost(self):
        print >> self.file, 'on_context_lost()'

    def on_context_state_lost(self):
        print >> self.file, 'on_context_state_lost()'


# symbolic names for the window events
EVENT_KEY_PRESS = WindowEventDispatcher.register_event_type('on_key_press')
EVENT_KEY_RELEASE = WindowEventDispatcher.register_event_type('on_key_release')
EVENT_TEXT = WindowEventDispatcher.register_event_type('on_text')
EVENT_TEXT_MOTION = WindowEventDispatcher.register_event_type('on_text_motion')
EVENT_TEXT_MOTION_SELECT = \
    WindowEventDispatcher.register_event_type('on_text_motion_select')
EVENT_MOUSE_MOTION = \
    WindowEventDispatcher.register_event_type('on_mouse_motion')
EVENT_MOUSE_DRAG = WindowEventDispatcher.register_event_type('on_mouse_drag')
EVENT_MOUSE_PRESS = WindowEventDispatcher.register_event_type('on_mouse_press')
EVENT_MOUSE_RELEASE = \
    WindowEventDispatcher.register_event_type('on_mouse_release')
EVENT_MOUSE_SCROLL = \
    WindowEventDispatcher.register_event_type('on_mouse_scroll')
EVENT_MOUSE_ENTER = WindowEventDispatcher.register_event_type('on_mouse_enter')
EVENT_MOUSE_LEAVE = WindowEventDispatcher.register_event_type('on_mouse_leave')
EVENT_CLOSE = WindowEventDispatcher.register_event_type('on_close')
EVENT_EXPOSE = WindowEventDispatcher.register_event_type('on_expose')
EVENT_RESIZE = WindowEventDispatcher.register_event_type('on_resize')
EVENT_MOVE = WindowEventDispatcher.register_event_type('on_move')
EVENT_ACTIVATE = WindowEventDispatcher.register_event_type('on_activate')
EVENT_DEACTIVATE = WindowEventDispatcher.register_event_type('on_deactivate')
EVENT_SHOW = WindowEventDispatcher.register_event_type('on_show')
EVENT_HIDE = WindowEventDispatcher.register_event_type('on_hide')
EVENT_CONTEXT_LOST = \
    WindowEventDispatcher.register_event_type('on_context_lost')
EVENT_CONTEXT_STATE_LOST = \
    WindowEventDispatcher.register_event_type('on_context_state_lost')
