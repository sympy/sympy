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

'''Linux ALSA audio implementation.
'''

__docformat__ = 'restructuredtext'
__version__ = '$Id: __init__.py 2084 2008-05-27 12:42:19Z Alex.Holkner $'

import ctypes

import pyglet
from pyglet.media import AudioPlayer, Listener, MediaException

from pyglet.media.drivers.alsa import asound

alsa_debug = None
if pyglet.options['debug_media']:
    alsa_debug = 'alsa.log'

class ALSAException(MediaException):
    pass

if alsa_debug:
    def check(err):
        if err < 0:
            raise ALSAException(asound.snd_strerror(err))
        return err
else:
    check = lambda v: v

class ALSAAudioPlayer(AudioPlayer):
    _device_name = 'default'
    _buffer_time = 0.3
    _min_write_bytes = 10000
    
    def __init__(self, audio_format):
        super(ALSAAudioPlayer, self).__init__(audio_format)

        format = {
            8:  asound.SND_PCM_FORMAT_U8,
            16: asound.SND_PCM_FORMAT_S16,
            24: asound.SND_PCM_FORMAT_S24,  # probably won't work
            32: asound.SND_PCM_FORMAT_S32
        }.get(audio_format.sample_size)
        if format is None:
            raise ALSAException('Unsupported audio format.')
            
        self.pcm = ctypes.POINTER(asound.snd_pcm_t)()
        self.hwparams = ctypes.POINTER(asound.snd_pcm_hw_params_t)()
        self.swparams = ctypes.POINTER(asound.snd_pcm_sw_params_t)() 

        check(asound.snd_pcm_open(ctypes.byref(self.pcm),
                                  self._device_name,
                                  asound.SND_PCM_STREAM_PLAYBACK,
                                  asound.SND_PCM_NONBLOCK))
        check(asound.snd_pcm_hw_params_malloc(ctypes.byref(self.hwparams)))
        check(asound.snd_pcm_sw_params_malloc(ctypes.byref(self.swparams)))
        check(asound.snd_pcm_hw_params_any(self.pcm, self.hwparams))

        check(asound.snd_pcm_hw_params_set_access(self.pcm, self.hwparams,
            asound.SND_PCM_ACCESS_RW_INTERLEAVED))
        check(asound.snd_pcm_hw_params_set_format(self.pcm, self.hwparams,
            format))
        check(asound.snd_pcm_hw_params_set_channels(self.pcm, self.hwparams,
            audio_format.channels))

        rate = ctypes.c_uint(audio_format.sample_rate)
        dir = ctypes.c_int(0)
        check(asound.snd_pcm_hw_params_set_rate_near(self.pcm, self.hwparams,
            rate, dir))
        # Note actual sample rate is in rate.value.  Ignored for now because
        # difference is negligible.

        buffer_size = int(self._buffer_time * audio_format.sample_rate)
        bs = asound.snd_pcm_uframes_t(buffer_size)
        check(asound.snd_pcm_hw_params_set_buffer_size_near(self.pcm,
            self.hwparams, bs))
        # Actual buffer size is in bs.value.  Ignored for now.

        check(asound.snd_pcm_hw_params(self.pcm, self.hwparams))

        if alsa_debug:
            asound.snd_output_printf(debug_output, 
                'New device: %s\n' % self._device_name)
            check(asound.snd_pcm_dump(self.pcm, debug_output))

        self.can_pause = asound.snd_pcm_hw_params_can_pause(self.hwparams)

        # List of (alsatime, timestamp)
        self._timestamps = []
        self._stop_alsatime = None
        self._eos_count = 0
        self._playing = False

    def __del__(self):
        try:
            check(asound.snd_pcm_close(self.pcm))
        except:
            pass

    def get_write_size(self):
        state = asound.snd_pcm_state(self.pcm)
        if state == asound.SND_PCM_STATE_PAUSED:
            return 0

        avail = max(0, asound.snd_pcm_avail_update(self.pcm))
        bytes = avail * self.audio_format.bytes_per_sample
        if bytes < self._min_write_bytes:
            return 0
        return bytes

    def write(self, audio_data):
        samples = audio_data.length // self.audio_format.bytes_per_sample
        samples_out = asound.snd_pcm_writei(self.pcm, audio_data.data,
                                            samples)
        if samples_out < 0:
            if samples_out == -11: # EAGAIN
                return
            elif samples_out == -32: # EPIPE (xrun)
                check(asound.snd_pcm_prepare(self.pcm))
                return
            else:
                raise ALSAException(asound.snd_strerror(samples_out))

        delay = asound.snd_pcm_sframes_t()
        check(asound.snd_pcm_delay(self.pcm, delay))
        alsatime = self._get_asound_time() + \
            delay.value / float(self.audio_format.sample_rate)

        self._timestamps.append((alsatime, audio_data.timestamp))

        audio_data.consume(samples_out * self.audio_format.bytes_per_sample,
                           self.audio_format)

    def write_eos(self):
        if self._timestamps:
            self._timestamps.append((None, None))

    def write_end(self):
        pass

    def play(self):
        if self._playing:
            return

        state = asound.snd_pcm_state(self.pcm)
        if self.can_pause and state == asound.SND_PCM_STATE_PAUSED:
            check(asound.snd_pcm_pause(self.pcm, 0))
        elif state not in (asound.SND_PCM_STATE_RUNNING,
                           asound.SND_PCM_STATE_PREPARED):
            check(asound.snd_pcm_prepare(self.pcm))
        self._playing = True

        if self._stop_alsatime is not None:
            diff = self._get_asound_time() - self._stop_alsatime
            self._timestamps = [(a + diff, t) for a, t in self._timestamps]
            self._stop_alsatime = None

    def stop(self):
        if not self._playing:
            return

        if self.can_pause and self._playing:
            check(asound.snd_pcm_pause(self.pcm, 1))
            self._stop_alsatime = self._get_asound_time()
        else:
            # Hardware can't pause, so we'll just drop everything that's
            # been buffered.  Improvement could be to rebuffer data that
            # wasn't played.
            self.clear()
        self._playing = False

    def clear(self):
        check(asound.snd_pcm_drop(self.pcm))
        self._stop_alsatime = None
        self._timestamps = []

    def _get_asound_time(self):
        status = ctypes.POINTER(asound.snd_pcm_status_t)()
        timestamp = asound.snd_timestamp_t()

        check(asound.snd_pcm_status_malloc(ctypes.byref(status)))
        check(asound.snd_pcm_status(self.pcm, status))
        asound.snd_pcm_status_get_tstamp(status, ctypes.byref(timestamp))
        asound.snd_pcm_status_free(status)
        return timestamp.tv_sec + timestamp.tv_usec * 0.000001

    def pump(self):
        underrun = False

        if self._stop_alsatime is not None:
            return underrun

        # Check that ALSA's still playing
        if self._playing:
            state = asound.snd_pcm_state(self.pcm)
            if state not in (asound.SND_PCM_STATE_RUNNING,
                             asound.SND_PCM_STATE_PREPARED):
                # Underrun!
                check(asound.snd_pcm_prepare(self.pcm))
                underrun = True

        alsatime = self._get_asound_time()
        try:
            while self._timestamps[0][0] < alsatime:
                self._timestamps.pop(0)
                while self._timestamps[0][0] is None:
                    self._eos_count += 1
                    self._timestamps.pop(0)
        except IndexError:
            pass

        return underrun

    def get_time(self):
        if self._stop_alsatime is None:
            alsatime = self._get_asound_time()
        else:
            alsatime = self._stop_alsatime

        if not self._timestamps:
            self._playing = False
            return 0.0

        alsats, ts = self._timestamps[0]
        t =  alsatime - alsats + ts
        return t

    def clear_eos(self):
        if self._eos_count > 0:
            self._eos_count -= 1
            return True
        return False

class ALSAListener(Listener):
    def _set_volume(self, volume):
        # TODO master volume
        self._volume = volume

    # All other properties are silently ignored.

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
    global debug_output
    debug_output = ctypes.POINTER(asound.snd_output_t)()
    if alsa_debug:
        check(asound.snd_output_stdio_open(ctypes.byref(debug_output),
                                           alsa_debug,
                                           'w'))

driver_listener = ALSAListener()
driver_audio_player_class = ALSAAudioPlayer
