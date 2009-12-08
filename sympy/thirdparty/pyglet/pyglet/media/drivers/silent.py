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

'''Fallback driver producing no audio.
'''

__docformat__ = 'restructuredtext'
__version__ = '$Id: silent.py 1680 2008-01-27 09:13:50Z Alex.Holkner $'

import time

from pyglet.media import AudioPlayer, Listener, AudioData
from pyglet.media import MediaException

class SilentAudioPlayer(AudioPlayer):
    UPDATE_PERIOD = 0.1
    
    def __init__(self, audio_format):
        super(SilentAudioPlayer, self).__init__(audio_format)

        self._playing = False
        self._eos_count = 0

        self._audio_data_list = []
        self._head_time = 0.0
        self._head_timestamp = 0.0
        self._head_system_time = time.time()

    def get_write_size(self):
        bytes = int(self.audio_format.bytes_per_second * self.UPDATE_PERIOD)
        return max(0, bytes - sum(
            [a.length for a in self._audio_data_list if a is not None]))

    def write(self, audio_data):
        if not self._audio_data_list:
            self._head_time = 0.0
            self._head_timestamp = audio_data.timestamp
            self._head_system_time = time.time()
        self._audio_data_list.append(
            AudioData(None, 
                      audio_data.length, 
                      audio_data.timestamp,
                      audio_data.duration))
        audio_data.consume(audio_data.length, self.audio_format)

    def write_eos(self):
        if self._audio_data_list:
            self._audio_data_list.append(None)

    def write_end(self):
        pass

    def play(self):
        self._playing = True
        self._head_system_time = time.time()

    def stop(self):
        self._playing = False
        self._head_time = time.time() - self._head_system_time

    def clear(self):
        self._audio_data_list = []
        self._head_time = 0.0
        self._head_system_time = time.time()
        self._eos_count = 0

    def pump(self):
        if not self._playing:
            return
        system_time = time.time()
        head_time = system_time - self._head_system_time
        try:
            while head_time >= self._audio_data_list[0].duration:
                head_time -= self._audio_data_list[0].duration
                self._audio_data_list.pop(0)
                while self._audio_data_list[0] is None:
                    self._eos_count += 1
                    self._audio_data_list.pop(0)
            self._head_timestamp = self._audio_data_list[0].timestamp
            self._head_system_time = system_time - head_time
        except IndexError:
            pass

    def get_time(self):
        if not self._audio_data_list:
            return time.time() - self._head_system_time + self._head_timestamp

        if self._playing:
            system_time = time.time()
            head_time = system_time - self._head_system_time
            return head_time + self._audio_data_list[0].timestamp 
        else:
            return self._audio_data_list[0].timestamp + self._head_time

    def clear_eos(self):
        if self._eos_count:
            self._eos_count -= 1 
            return True
        return False

class SilentListener(Listener):
    def _set_volume(self, volume):
        self._volume = volume

    def _set_position(self, position):
        self._position = position

    def _set_velocity(self, velocity):
        self._velocity = velocity

    def _set_forward_orientation(self, orientation):
        self._forward_orientation = orientation

    def _set_up_orientation(self, orientation):
        self._up_orientation = orientation

    def _set_doppler_factor(self, factor):
        self._doppler_factor = factor

    def _set_speed_of_sound(self, speed_of_sound):
        self._speed_of_sound = speed_of_sound

def driver_init():
    pass

driver_listener = SilentListener()
driver_audio_player_class = SilentAudioPlayer
