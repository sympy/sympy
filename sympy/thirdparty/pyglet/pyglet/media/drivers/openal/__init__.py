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
# $Id: __init__.py 1352 2007-11-04 03:35:06Z Alex.Holkner $

import ctypes
import sys
import time

from pyglet.media import AudioPlayer, Listener, MediaException

from pyglet.media.drivers.openal import lib_openal as al
from pyglet.media.drivers.openal import lib_alc as alc

class OpenALException(MediaException):
    pass

def _split_nul_strings(s):
    # NUL-separated list of strings, double-NUL-terminated.
    nul = False
    i = 0
    while True:
        if s[i] == '\0':
            if nul:
                break
            else:
                nul = True
        else:
            nul = False
        i += 1
    s = s[:i - 1]
    return s.split('\0')

def get_version():
    major = alc.ALCint()
    minor = alc.ALCint()
    alc.alcGetIntegerv(_device, alc.ALC_MAJOR_VERSION, 
                       ctypes.sizeof(major), major)
    alc.alcGetIntegerv(_device, alc.ALC_MINOR_VERSION, 
                       ctypes.sizeof(minor), minor)
    return major.value, minor.value

def have_version(major, minor):
    return (major, minor) <= get_version()

def get_extensions():
    extensions = alc.alcGetString(_device, alc.ALC_EXTENSIONS)
    if sys.platform == 'darwin':
        return ctypes.cast(extensions, ctypes.c_char_p).value.split(' ')
    else:
        return _split_nul_strings(extensions)

def have_extension(extension):
    return extension in get_extensions()

format_map = {
    (1,  8): al.AL_FORMAT_MONO8,
    (1, 16): al.AL_FORMAT_MONO16,
    (2,  8): al.AL_FORMAT_STEREO8,
    (2, 16): al.AL_FORMAT_STEREO16,
}

class OpenALAudioPlayer(AudioPlayer):
    #: Seconds ahead to buffer audio.  Keep small for low latency, but large
    #: enough to avoid underruns. (0.05 is the minimum for my 2.2 GHz Linux)
    _update_buffer_time = 0.2

    #: Minimum size of an OpenAL buffer worth bothering with
    _min_buffer_size = 512

    #: Maximum size of an OpenAL buffer, in bytes.  TODO: use OpenAL maximum
    _max_buffer_size = 65536

    def __init__(self, audio_format):
        super(OpenALAudioPlayer, self).__init__(audio_format)

        try:
            self._al_format = format_map[(audio_format.channels,
                                          audio_format.sample_size)]
        except KeyError:
            raise OpenALException('Unsupported audio format.')

        self._al_source = al.ALuint()
        al.alGenSources(1, self._al_source)

        # Seconds of audio currently queued not processed (estimate)
        self._buffered_time = 0.0

        # Seconds of audio into current (head) buffer
        self._current_buffer_time = 0.0

        # List of (timestamp, duration) corresponding to currently queued AL
        # buffers
        self._timestamps = []

        # OpenAL 1.0 timestamp interpolation
        self._timestamp_system_time = 0.0

        # Desired play state (True even if stopped due to underrun)
        self._playing = False

        # Timestamp when paused
        self._pause_timestamp = 0.0

        self._eos_count = 0

    def __del__(self):
        try:
            al.alDeleteSources(1, self._al_source)
        except:
            pass

    def get_write_size(self):
        t = self._buffered_time - self._current_buffer_time
        size = int(max(0, self._update_buffer_time - t) * \
            self.audio_format.bytes_per_second)
        if size < self._min_buffer_size:
            size = 0
        return size

    def write(self, audio_data):
        buffer = al.ALuint()
        al.alGenBuffers(1, buffer)
        al.alBufferData(buffer, 
                        self._al_format,
                        audio_data.data,
                        audio_data.length,
                        self.audio_format.sample_rate)
        al.alSourceQueueBuffers(self._al_source, 1, ctypes.byref(buffer)) 

        self._buffered_time += audio_data.duration
        self._timestamps.append((audio_data.timestamp, audio_data.duration))
        audio_data.consume(audio_data.length, self.audio_format)

    def write_eos(self):
        if self._timestamps:
            self._timestamps.append((None, None))

    def write_end(self):
        pass

    def play(self):
        if self._playing:
            return

        self._playing = True
        self._al_play()
        if not _have_1_1:
            self._timestamp_system_time = time.time()

    def _al_play(self):
        if not self._timestamps:
            return
        state = al.ALint()
        al.alGetSourcei(self._al_source, al.AL_SOURCE_STATE, state)
        if state.value != al.AL_PLAYING:
            al.alSourcePlay(self._al_source)

    def stop(self):
        if not self._playing:
            return

        self._pause_timestamp = self.get_time()
        al.alSourcePause(self._al_source)
        self._playing = False

    def clear(self):
        al.alSourceStop(self._al_source)
        self._playing = False

        processed = al.ALint()
        al.alGetSourcei(self._al_source, al.AL_BUFFERS_PROCESSED, processed)
        if processed.value:
            buffers = (al.ALuint * processed.value)()
            al.alSourceUnqueueBuffers(self._al_source, len(buffers), buffers)
            al.alDeleteBuffers(len(buffers), buffers)

        self._pause_timestamp = 0.0
        self._buffered_time = 0.0
        self._timestamps = []

    def pump(self):
        # Release spent buffers
        processed = al.ALint()
        al.alGetSourcei(self._al_source, al.AL_BUFFERS_PROCESSED, processed)
        processed = processed.value
        if processed:
            buffers = (al.ALuint * processed)()
            al.alSourceUnqueueBuffers(self._al_source, len(buffers), buffers)
            al.alDeleteBuffers(len(buffers), buffers)

        # Pop timestamps and check for eos markers
        try:
            while processed:
                if not _have_1_1:
                    self._timestamp_system_time = time.time()
                _, duration = self._timestamps.pop(0)
                self._buffered_time -= duration
                while self._timestamps[0][0] is None:
                    self._eos_count += 1
                    self._timestamps.pop(0)
                processed -= 1
        except IndexError:
            pass

        if _have_1_1:
            samples = al.ALint()
            al.alGetSourcei(self._al_source, al.AL_SAMPLE_OFFSET, samples)
            self._current_buffer_time = samples.value / \
                float(self.audio_format.sample_rate)
        else:
            # Interpolate system time past buffer timestamp
            self._current_buffer_time = time.time() - \
                self._timestamp_system_time

        # Begin playing if underrun previously
        if self._playing:
            self._al_play()

    def get_time(self):
        state = al.ALint()
        al.alGetSourcei(self._al_source, al.AL_SOURCE_STATE, state)
        if state.value != al.AL_PLAYING:
            return self._pause_timestamp

        if not self._timestamps:
            return self._pause_timestamp

        ts, _ = self._timestamps[0]

        return ts + self._current_buffer_time

    def clear_eos(self):
        while self._eos_count > 0:
            self._eos_count -= 1
            return True
        return False

    def set_volume(self, volume):
        al.alSourcef(self._al_source, al.AL_GAIN, max(0, volume))

    def set_position(self, position):
        x, y, z = position
        al.alSource3f(self._al_source, al.AL_POSITION, x, y, z)

    def set_min_distance(self, min_distance):
        al.alSourcef(self._al_source, al.AL_REFERENCE_DISTANCE, min_distance)

    def set_max_distance(self, max_distance):
        al.alSourcef(self._al_source, al.AL_MAX_DISTANCE, max_distance)

    def set_pitch(self, pitch):
        al.alSourcef(self._al_source, al.AL_PITCH, max(0, pitch))

    def set_cone_orientation(self, cone_orientation):
        x, y, z = cone_orientation
        al.alSource3f(self._al_source, al.AL_DIRECTION, x, y, z)

    def set_cone_inner_angle(self, cone_inner_angle):
        al.alSourcef(self._al_source, al.AL_CONE_INNER_ANGLE, cone_inner_angle)

    def set_cone_outer_angle(self, cone_outer_angle):
        al.alSourcef(self._al_source, al.AL_CONE_OUTER_ANGLE, cone_outer_angle)

    def set_cone_outer_gain(self, cone_outer_gain):
        al.alSourcef(self._al_source, al.AL_CONE_OUTER_GAIN, cone_outer_gain)

class OpenALListener(Listener):
    def _set_volume(self, volume):
        al.alListenerf(al.AL_GAIN, volume)
        self._volume = volume

    def _set_position(self, position):
        x, y, z = position
        al.alListener3f(al.AL_POSITION, x, y, z)
        self._position = position 

    def _set_forward_orientation(self, orientation):
        val = (al.ALfloat * 6)(*(orientation + self._up_orientation))
        al.alListenerfv(al.AL_ORIENTATION, val)
        self._forward_orientation = orientation

    def _set_up_orientation(self, orientation):
        val = (al.ALfloat * 6)(*(self._forward_orientation + orientation))
        al.alListenerfv(al.AL_ORIENTATION, val)
        self._up_orientation = orientation

_device = None
_have_1_1 = False

def driver_init(device_name = None):
    global _device
    global _have_1_1

    # TODO devices must be enumerated on Windows, otherwise 1.0 context is
    # returned.

    _device = alc.alcOpenDevice(device_name)
    if not _device:
        raise Exception('No OpenAL device.')

    alcontext = alc.alcCreateContext(_device, None)
    alc.alcMakeContextCurrent(alcontext)

    if have_version(1, 1):
        # Good version info to cache
        _have_1_1 = True

    # See issue #163.
    import sys
    if sys.platform in ('win32', 'cygwin'):
        from pyglet import clock
        clock.Clock._force_sleep = True

driver_listener = OpenALListener()
driver_audio_player_class = OpenALAudioPlayer

