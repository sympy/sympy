#!/usr/bin/python
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
# $Id: openal.py 1045 2007-07-17 02:40:18Z r1chardj0n3s $

import ctypes
import sys
import time

from pyglet.media import Sound, Listener, EVENT_FINISHED
from pyglet.media import lib_openal as al
from pyglet.media import lib_alc as alc

_device = None
_is_init = False
_have_1_1 = False
def init(device_name = None):
    global _device
    global _is_init
    global _have_1_1

    _device = alc.alcOpenDevice(device_name)
    if not _device:
        raise Exception('No OpenAL device.')

    alcontext = alc.alcCreateContext(_device, None)
    alc.alcMakeContextCurrent(alcontext)
    _is_init = True

    if have_version(1, 1):
        # Good version info to cache
        _have_1_1 = True

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

class BufferPool(list):
    def __init__(self):
        self.timestamps = {}

    def get(self, timestamp):
        if not self:
            buffer = al.ALuint()
            al.alGenBuffers(1, buffer)
        else:
            buffer = al.ALuint(self.pop(0))
        self.timestamps[buffer.value] = timestamp
        return buffer

    def timestamp(self, buffer):
        return self.timestamps[buffer]

    def replace(self, buffers):
        self.extend(buffers)

    def __del__(self, al=al):
        if al and al.ALuint and al.alDeleteBuffers:
            buffers = (al.ALuint * len(self))(*self)
            al.alDeleteBuffers(len(self), buffers)

buffer_pool = BufferPool()

_format_map = {
    (1,  8): al.AL_FORMAT_MONO8,
    (1, 16): al.AL_FORMAT_MONO16,
    (2,  8): al.AL_FORMAT_STEREO8,
    (2, 16): al.AL_FORMAT_STEREO16,
}
def get_format(channels, depth):
    return _format_map[channels, depth]

class OpenALSound(Sound):
    _processed_buffers = 0
    _queued_buffers = 0

    def __init__(self):
        super(OpenALSound, self).__init__()

        self.source = al.ALuint()
        al.alGenSources(1, self.source)
        self.play_when_buffered = False

        # OpenAL on Linux lacks the time functions, so this is a stab at
        # interpolating time between the known buffer timestamps.  When a
        # timestamp is read from the active buffer, the current system time is
        # stored in self._last_buffer_time, and the buffer that was used is
        # stored in self._last_buffer. 
        #
        # If the same buffer is in use the next time the time is requested,
        # the elapsed system time is added to the buffer timestamp.
        #
        # The (completely invalid) assumption is that the buffer is just
        # beginning when _last_buffer_time is stored.  This is more likely
        # than not, at least, if the buffers are long.  If the buffers are
        # short, hopefully not much interpolation will be needed and the
        # difference won't be noticeable.
        #
        # The alternative -- reusing the same timestamp without adding system
        # time -- results in slightly jumpy video due to lost frames.
        #
        # When OpenAL 1.1 is present (i.e., on OS X), we add the buffer's
        # timestamp to the retrieved sample offset within the current buffer.
        #
        # XXX Need special consideration when pausing (Linux only).
        self._last_buffer = 0
        self._last_buffer_time = 0

    def __del__(self, al=al):
        if al and al.alDeleteSources:
            al.alDeleteSources(1, self.source)

    def _get_time(self):
        buffer = al.ALint()
        al.alGetSourcei(self.source, al.AL_BUFFER, buffer)
        if not buffer:
            return 0.

        # The playback position at the start of the current buffer
        buffer_timestamp = buffer_pool.timestamp(buffer.value)

        if _have_1_1:
            # Add buffer timestamp to sample offset
            # XXX this occasionally goes backwards
            buffer_samples = al.ALint()
            al.alGetSourcei(self.source, al.AL_SAMPLE_OFFSET, buffer_samples)
            sample_rate = al.ALint()
            al.alGetBufferi(buffer.value, al.AL_FREQUENCY, sample_rate)
            buffer_time = buffer_samples.value / float(sample_rate.value)
            return buffer_timestamp + buffer_time
        else:
            # Interpolate system time past buffer timestamp
            if not self.playing:
                return self._last_buffer_time + buffer_timestamp
            elif buffer.value == self._last_buffer:
                return time.time() - self._last_buffer_time + buffer_timestamp
            else:
                self._last_buffer = buffer.value
                self._last_buffer_time = time.time()
                return buffer_timestamp

    def play(self):
        self._openal_play()

    def _openal_play(self):
        if self.playing:
            return

        buffers = al.ALint()
        al.alGetSourcei(self.source, al.AL_BUFFERS_QUEUED, buffers)
        if buffers.value == 0:
            self.play_when_buffered = True
        else:
            al.alSourcePlay(self.source)
            self.play_when_buffered = False
            self.playing = True
            self._last_buffer_time = time.time() - self._last_buffer_time

    def pause(self):
        self.playing = False
        self._last_buffer_time = time.time() - self._last_buffer_time
        al.alSourcePause(self.source)

    def stop(self):
        self.playing = False
        self.finished = True
        al.alSourceStop(self.source)

    def _set_volume(self, volume):
        al.alSourcef(self.source, al.AL_GAIN, max(0, volume))
        self._volume = volume

    def _set_min_gain(self, min_gain):
        al.alSourcef(self.source, al.AL_MIN_GAIN, max(0, min_gain))
        self._min_gain = min_gain

    def _set_max_gain(self, max_gain):
        al.alSourcef(self.source, al.AL_MAX_GAIN, max(0, max_gain))
        self._max_gain = max_gain

    def _set_position(self, position):
        x, y, z = position
        al.alSource3f(self.source, al.AL_POSITION, x, y, z)
        self._position = position

    def _set_velocity(self, velocity):
        x, y, z = velocity
        al.alSource3f(self.source, al.AL_VELOCITY, x, y, z)
        self._velocity = velocity

    def _set_pitch(self, pitch):
        al.alSourcef(self.source, al.AL_PITCH, max(0, pitch))
        self._pitch = pitch

    def _set_cone_orientation(self, cone_orientation):
        x, y, z = cone_orientation
        al.alSource3f(self.source, al.AL_DIRECTION, x, y, z)
        self._cone_orientation = cone_orientation

    def _set_cone_inner_angle(self, cone_inner_angle):
        al.alSourcef(self.source, al.AL_CONE_INNER_ANGLE, cone_inner_angle)
        self._cone_inner_angle = cone_inner_angle

    def _set_cone_outer_angle(self, cone_outer_angle):
        al.alSourcef(self.source, al.AL_CONE_OUTER_ANGLE, cone_outer_angle)
        self._cone_outer_angle = cone_outer_angle

    def _set_cone_outer_gain(self, cone_outer_gain):
        al.alSourcef(self.source, al.AL_CONE_OUTER_GAIN, cone_outer_gain)
        self._cone_outer_gain = cone_outer_gain

    def dispatch_events(self):
        queued = al.ALint()
        processed = al.ALint()
        al.alGetSourcei(self.source, al.AL_BUFFERS_QUEUED, queued)
        al.alGetSourcei(self.source, al.AL_BUFFERS_PROCESSED, processed)
        if processed.value == queued.value:
            self.finished = True
            self.playing = False
            self.dispatch_event(EVENT_FINISHED)
        self._processed_buffers = processed.value
        self._queued_buffers = queued.value

        if self.play_when_buffered and queued.value:
            self._openal_play()

class OpenALStreamingSound(OpenALSound):
    def __del__(self):
        try:
            self.stop()
        except (ValueError, TypeError, NameError):
            # we're being __del__'ed during interpreter exit
            # (ValueError would come from trying to unschedule us from the
            # event instances list, NameError from trying to use already-
            # collected modules and TypeError from trying to invoke methods
            # on collected objects)
            pass

    def stop(self):
        super(OpenALStreamingSound, self).stop()

        # Release all buffers
        queued = al.ALint()
        al.alGetSourcei(self.source, al.AL_BUFFERS_QUEUED, queued)
        self._release_buffers(queued.value)

    def dispatch_events(self):
        super(OpenALStreamingSound, self).dispatch_events()

        # Release spent buffers
        if self._processed_buffers:
            self._release_buffers(self._processed_buffers)

    def _release_buffers(self, num_buffers):
        discard_buffers = (al.ALuint * num_buffers)()
        al.alSourceUnqueueBuffers(
            self.source, len(discard_buffers), discard_buffers)
        buffer_pool.replace(discard_buffers)


class OpenALStaticSound(OpenALSound):
    def __init__(self, medium):
        super(OpenALStaticSound, self).__init__()

        # Keep a reference to the medium to avoid premature release of
        # buffers.
        self.medium = medium

    def stop(self):
        super(OpenALSound, self).stop()

        self.medium = None

class OpenALListener(Listener):
    def _set_position(self, position):
        x, y, z = position
        al.alListener3f(al.AL_POSITION, x, y, z)
        self._position = position 

    def _set_velocity(self, velocity):
        x, y, z = velocity
        al.alListener3f(al.AL_VELOCITY, x, y, z)
        self._velocity = velocity 

    def _set_forward_orientation(self, orientation):
        val = (ALfloat * 6)(*(orientation + self._up_orientation))
        al.alListenerfv(al.AL_ORIENTATION, val)
        self._forward_orientation = orientation

    def _set_up_orientation(self, orientation):
        val = (ALfloat * 6)(*(self._forward_orientation + orientation))
        al.alListenerfv(al.AL_ORIENTATION, val)
        self._up_orientation = orientation

    def _set_doppler_factor(self, factor):
        al.alDopplerFactor(factor)
        self._doppler_factor = factor

    def _set_speed_of_sound(self, speed_of_sound):
        al.alSpeedOfSound(speed_of_sound)
        self._speed_of_sound = speed_of_sound
