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

from pyglet.media import Source, AudioFormat, AudioData

import ctypes
import os
import math

class ProceduralSource(Source):
    def __init__(self, duration, sample_rate=44800, sample_size=16):
        self._duration = float(duration)
        self.audio_format = AudioFormat(
            channels=1,
            sample_size=sample_size,
            sample_rate=sample_rate)

        self._offset = 0
        self._bytes_per_sample = sample_size >> 3
        self._bytes_per_second = self._bytes_per_sample * sample_rate
        self._max_offset = int(self._bytes_per_second * self._duration)
        
        if self._bytes_per_sample == 2:
            self._max_offset &= 0xfffffffe

    def _get_audio_data(self, bytes):
        bytes = min(bytes, self._max_offset - self._offset)
        if bytes <= 0:
            return None
        
        timestamp = float(self._offset) / self._bytes_per_second
        duration = float(bytes) / self._bytes_per_second
        data = self._generate_data(bytes, self._offset)
        self._offset += bytes
        is_eos = self._offset >= self._max_offset

        return AudioData(data,
                         bytes,
                         timestamp,
                         duration)

    def _generate_data(self, bytes, offset):
        '''Generate `bytes` bytes of data.

        Return data as ctypes array or string.
        '''
        raise NotImplementedError('abstract')

    def _seek(self, timestamp):
        self._offset = int(timestamp * self._bytes_per_second)

        # Bound within duration
        self._offset = min(max(self._offset, 0), self._max_offset)

        # Align to sample
        if self._bytes_per_sample == 2:
            self._offset &= 0xfffffffe

class Silence(ProceduralSource):
    def _generate_data(self, bytes, offset):
        if self._bytes_per_sample == 1:
            return '\127' * bytes
        else:
            return '\0' * bytes

class WhiteNoise(ProceduralSource):
    def _generate_data(self, bytes, offset):
        return os.urandom(bytes)

class Sine(ProceduralSource):
    def __init__(self, duration, frequency=440, **kwargs):
        super(Sine, self).__init__(duration, **kwargs)
        self.frequency = frequency
        
    def _generate_data(self, bytes, offset):
        if self._bytes_per_sample == 1:
            start = offset
            samples = bytes
            bias = 127
            amplitude = 127
            data = (ctypes.c_ubyte * samples)()
        else:
            start = offset >> 1
            samples = bytes >> 1
            bias = 0
            amplitude = 32767
            data = (ctypes.c_short * samples)()
        step = self.frequency * (math.pi * 2) / self.audio_format.sample_rate
        for i in range(samples):
            data[i] = int(math.sin(step * (i + start)) * amplitude + bias)
        return data

class Saw(ProceduralSource):
    def __init__(self, duration, frequency=440, **kwargs):
        super(Saw, self).__init__(duration, **kwargs)
        self.frequency = frequency
        
    def _generate_data(self, bytes, offset):
        # XXX TODO consider offset
        if self._bytes_per_sample == 1:
            samples = bytes
            value = 127
            max = 255
            min = 0
            data = (ctypes.c_ubyte * samples)()
        else:
            samples = bytes >> 1
            value = 0
            max = 32767
            min = -32768
            data = (ctypes.c_short * samples)()
        step = (max - min) * 2 * self.frequency / self.audio_format.sample_rate
        for i in range(samples):
            value += step
            if value > max:
                value = max - (value - max)
                step = -step
            if value < min:
                value = min - (value - min)
                step = -step
            data[i] = value
        return data

class Square(ProceduralSource):
    def __init__(self, duration, frequency=440, **kwargs):
        super(Square, self).__init__(duration, **kwargs)
        self.frequency = frequency
        
    def _generate_data(self, bytes, offset):
        # XXX TODO consider offset
        if self._bytes_per_sample == 1:
            samples = bytes
            value = 0
            amplitude = 255
            data = (ctypes.c_ubyte * samples)()
        else:
            samples = bytes >> 1
            value = -32768
            amplitude = 65535
            data = (ctypes.c_short * samples)()
        period = self.audio_format.sample_rate / self.frequency / 2
        count = 0
        for i in range(samples):
            count += 1
            if count == period:
                value = amplitude - value
                count = 0
            data[i] = value
        return data
