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

'''
'''

__docformat__ = 'restructuredtext'
__version__ = '$Id: xlib.py 2045 2008-05-01 16:59:49Z Alex.Holkner $'

import select
import weakref

from pyglet.app import displays, windows, BaseEventLoop
from pyglet.window.xlib import xlib

class XlibEventLoop(BaseEventLoop):
    def run(self):
        self._setup()

        e = xlib.XEvent()
        t = 0
        sleep_time = 0.

        self.dispatch_event('on_enter')

        while not self.has_exit:
            # Check for already pending events
            for display in displays:
                if xlib.XPending(display._display):
                    pending_displays = (display,)
                    break
            else:
                # None found; select on all file descriptors or timeout
                iwtd = self.get_select_files()
                pending_displays, _, _ = select.select(iwtd, (), (), sleep_time)

            # Dispatch platform events
            for display in pending_displays:
                while xlib.XPending(display._display):
                    xlib.XNextEvent(display._display, e)

                    # Key events are filtered by the xlib window event
                    # handler so they get a shot at the prefiltered event.
                    if e.xany.type not in (xlib.KeyPress, xlib.KeyRelease):
                        if xlib.XFilterEvent(e, e.xany.window):
                            continue
                    try:
                        window = display._window_map[e.xany.window]
                    except KeyError:
                        continue

                    window.dispatch_platform_event(e)

            # Dispatch resize events
            for window in windows:
                if window._needs_resize:
                    window.switch_to()
                    window.dispatch_event('on_resize', 
                                          window._width, window._height)
                    window.dispatch_event('on_expose')
                    window._needs_resize = False

            sleep_time = self.idle()

        self.dispatch_event('on_exit')

    def get_select_files(self):
        return list(displays)
