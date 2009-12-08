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
# $Id:$

__docformat__ = 'restructuredtext'
__version__ = '$Id: $'

import ctypes
import time

from pyglet.app import windows, BaseEventLoop
from pyglet.window.win32 import _user32, types, constants

class Win32EventLoop(BaseEventLoop):
    def run(self):
        self._setup()

        self._timer_proc = types.TIMERPROC(self._timer_func)
        self._timer = timer = _user32.SetTimer(0, 0, 0, self._timer_proc)
        self._polling = False
        self._allow_polling = True
        msg = types.MSG()
        
        self.dispatch_event('on_enter')

        while not self.has_exit:
            if self._polling:
                while _user32.PeekMessageW(ctypes.byref(msg), 
                                           0, 0, 0, constants.PM_REMOVE):
                    _user32.TranslateMessage(ctypes.byref(msg))
                    _user32.DispatchMessageW(ctypes.byref(msg))
                self._timer_func(0, 0, timer, 0)
            else:
                _user32.GetMessageW(ctypes.byref(msg), 0, 0, 0)
                _user32.TranslateMessage(ctypes.byref(msg))
                _user32.DispatchMessageW(ctypes.byref(msg))
            
                # Manual idle event
                msg_types = \
                    _user32.GetQueueStatus(constants.QS_ALLINPUT) & 0xffff0000
                if (msg.message != constants.WM_TIMER and
                    not msg_types & ~(constants.QS_TIMER<<16)):
                    self._timer_func(0, 0, timer, 0)

        self.dispatch_event('on_exit')

    def _idle_chance(self):
        if (self._next_idle_time is not None and 
            self._next_idle_time <= time.time()):
            self._timer_func(0, 0, self._timer, 0)

    def _timer_func(self, hwnd, msg, timer, t):
        sleep_time = self.idle()

        if sleep_time is None:
            # Block indefinitely
            millis = constants.USER_TIMER_MAXIMUM
            self._next_idle_time = None
            self._polling = False
            _user32.SetTimer(0, timer, millis, self._timer_proc)
        elif sleep_time < 0.01 and self._allow_polling:
            # Degenerate to polling
            millis = constants.USER_TIMER_MAXIMUM
            self._next_idle_time = 0.
            if not self._polling:
                self._polling = True
                _user32.SetTimer(0, timer, millis, self._timer_proc)
        else:
            # Block until timer
            # XXX hack to avoid oversleep; needs to be api
            sleep_time = max(sleep_time - 0.01, 0) 
            millis = int(sleep_time * 1000)
            self._next_idle_time = time.time() + sleep_time
            self._polling = False
            _user32.SetTimer(0, timer, millis, self._timer_proc)
