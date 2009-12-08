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

'''Application-wide functionality.

Most applications need only call `run` after creating one or more windows
to begin processing events.  For example, a simple application consisting of
one window is::

    from pyglet import app
    from pyglet import window

    win = window.Window()
    app.run()

To handle events on the main event loop, instantiate it manually.  The
following example exits the application as soon as any window is closed (the
default policy is to wait until all windows are closed)::

    event_loop = app.EventLoop()

    @event_loop.event
    def on_window_close(window):
        event_loop.exit()

:since: pyglet 1.1
'''

__docformat__ = 'restructuredtext'
__version__ = '$Id: __init__.py 2140 2008-07-27 04:15:52Z Alex.Holkner $'

import sys
import weakref

from pyglet import clock
from pyglet import event

_is_epydoc = hasattr(sys, 'is_epydoc') and sys.is_epydoc

class WeakSet(object):
    '''Set of objects, referenced weakly.

    Adding an object to this set does not prevent it from being garbage
    collected.  Upon being garbage collected, the object is automatically
    removed from the set.
    '''
    def __init__(self):
        self._dict = weakref.WeakKeyDictionary()

    def add(self, value):
        self._dict[value] = True

    def remove(self, value):
        del self._dict[value]

    def __iter__(self):
        for key in self._dict.keys():
            yield key

    def __contains__(self, other):
        return other in self._dict

    def __len__(self):
        return len(self._dict)

#: Set of all open displays.  Instances of `Display` are automatically added
#: to this set upon construction.  The set uses weak references, so displays
#: are removed from the set when they are no longer referenced.
#:
#: :type: `WeakSet`
displays = WeakSet()

#: Set of all open windows (including invisible windows).  Instances of
#: `Window` are automatically added to this set upon construction.  The set
#: uses weak references, so windows are removed from the set when they are no
#: longer referenced or are closed explicitly.
#:
#: :type: `WeakSet`
windows = WeakSet()


class BaseEventLoop(event.EventDispatcher):
    '''The main run loop of the application.

    Calling `run` begins the application event loop, which processes
    operating system events, calls `pyglet.clock.tick` to call scheduled
    functions and calls `pyglet.window.Window.on_draw` and
    `pyglet.window.Window.flip` to update window contents.

    Applications can subclass `EventLoop` and override certain methods
    to integrate another framework's run loop, or to customise processing
    in some other way.  You should not in general override `run`, as
    this method contains platform-specific code that ensures the application
    remains responsive to the user while keeping CPU usage to a minimum.
    '''

    #: Flag indicating if the event loop will exit in the next iteration.
    #: This is
    has_exit = False

    def run(self):
        '''Begin processing events, scheduled functions and window updates.

        This method returns when `has_exit` is set to True.

        Developers are discouraged from overriding this method, as the
        implementation is platform-specific.
        '''
        raise NotImplementedError('abstract')

    def _setup(self):
        global event_loop
        event_loop = self

        # Disable event queuing for dispatch_events
        from pyglet.window import Window
        Window._enable_event_queue = False
        
        # Dispatch pending events
        for window in windows:
            window.switch_to()
            window.dispatch_pending_events()

    def _idle_chance(self):
        '''If timeout has expired, manually force an idle loop.

        Called by window that have blocked the event loop (e.g. during
        resizing).
        '''

    def idle(self):
        '''Called during each iteration of the event loop.

        The method is called immediately after any window events (i.e., after
        any user input).  The method can return a duration after which
        the idle method will be called again.  The method may be called
        earlier if the user creates more input events.  The method
        can return `None` to only wait for user events.

        For example, return ``1.0`` to have the idle method called every
        second, or immediately after any user events.

        The default implementation dispatches the
        `pyglet.window.Window.on_draw` event for all windows and uses
        `pyglet.clock.tick` and `pyglet.clock.get_sleep_time` on the default
        clock to determine the return value.

        This method should be overridden by advanced users only.  To have
        code execute at regular intervals, use the
        `pyglet.clock.schedule` methods.

        :rtype: float
        :return: The number of seconds before the idle method should
            be called again, or `None` to block for user input.
        '''
        dt = clock.tick(True)

        # Redraw all windows
        for window in windows:
            if window.invalid:
                window.switch_to()
                window.dispatch_event('on_draw')
                window.flip()

        # Update timout
        return clock.get_sleep_time(True)

    def exit(self):
        '''Safely exit the event loop at the end of the current iteration.

        This method is convenience for setting `has_exit` to ``True``.
        '''
        self.has_exit = True

    def on_window_close(self, window):
        '''Default window close handler.'''
        if not windows:
            self.exit()

    if _is_epydoc:
        def on_window_close(window):
            '''A window was closed.

            This event is dispatched when a window is closed.  It is not
            dispatched if the window's close button was pressed but the
            window did not close.

            The default handler calls `exit` if no more windows are open.  You
            can override this handler to base your application exit on some
            other policy.

            :event:
            '''

        def on_enter():
            '''The event loop is about to begin.

            This is dispatched when the event loop is prepared to enter
            the main run loop, and represents the last chance for an 
            application to initialise itself.

            :event:
            '''

        def on_exit():
            '''The event loop is about to exit.

            After dispatching this event, the `run` method returns (the
            application may not actually exit if you have more code
            following the `run` invocation).

            :event:
            '''

BaseEventLoop.register_event_type('on_window_close')
BaseEventLoop.register_event_type('on_enter')
BaseEventLoop.register_event_type('on_exit')

#: The global event loop.  Set to the correct instance when an `EventLoop` is
#: started.
#:
#: :type: `EventLoop`
event_loop = None

def run():
    '''Begin processing events, scheduled functions and window updates.

    This is a convenience function, equivalent to::

        EventLoop().run()

    '''
    EventLoop().run()

def exit():
    '''Exit the application event loop.

    Causes the application event loop to finish, if an event loop is currently
    running.  The application may not necessarily exit (for example, there may
    be additional code following the `run` invocation).

    This is a convenience function, equivalent to::

        event_loop.exit()

    '''
    if event_loop:
        event_loop.exit()

if _is_epydoc:
    EventLoop = BaseEventLoop
    EventLoop.__name__ = 'EventLoop'
    del BaseEventLoop
else:
    # Permit cyclic import.
    import pyglet
    pyglet.app = sys.modules[__name__]

    if sys.platform == 'darwin':
        from pyglet.app.carbon import CarbonEventLoop as EventLoop
    elif sys.platform in ('win32', 'cygwin'):
        from pyglet.app.win32 import Win32EventLoop as EventLoop
    else:
        from pyglet.app.xlib import XlibEventLoop as EventLoop


