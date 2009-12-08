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
# $Id: __init__.py 2154 2008-08-05 22:53:37Z Alex.Holkner $

'''Audio and video playback.

pyglet can play WAV files, and if AVbin is installed, many other audio and
video formats.

Playback is handled by the `Player` class, which reads raw data from `Source`
objects and provides methods for pausing, seeking, adjusting the volume, and
so on.  The `Player` class implements a the best available audio device
(currently, only OpenAL is supported)::

    player = Player()

A `Source` is used to decode arbitrary audio and video files.  It is
associated with a single player by "queuing" it::

    source = load('background_music.mp3')
    player.queue(source)

Use the `Player` to control playback.  

If the source contains video, the `Source.video_format` attribute will be
non-None, and the `Player.texture` attribute will contain the current video
image synchronised to the audio.

Decoding sounds can be processor-intensive and may introduce latency,
particularly for short sounds that must be played quickly, such as bullets or
explosions.  You can force such sounds to be decoded and retained in memory
rather than streamed from disk by wrapping the source in a `StaticSource`::

    bullet_sound = StaticSource(load('bullet.wav'))

The other advantage of a `StaticSource` is that it can be queued on any number
of players, and so played many times simultaneously.

'''

__docformat__ = 'restructuredtext'
__version__ = '$Id: __init__.py 2154 2008-08-05 22:53:37Z Alex.Holkner $'

import ctypes
import sys
import time
import StringIO

import pyglet
from pyglet import clock
from pyglet import event

_debug_media = pyglet.options['debug_media']

class MediaException(Exception):
    pass

class MediaFormatException(MediaException):
    pass

class CannotSeekException(MediaException):
    pass

class AudioFormat(object):
    '''Audio details.

    An instance of this class is provided by sources with audio tracks.  You
    should not modify the fields, as they are used internally to describe the
    format of data provided by the source.

    :Ivariables:
        `channels` : int
            The number of channels: 1 for mono or 2 for stereo (pyglet does
            not yet support surround-sound sources).
        `sample_size` : int
            Bits per sample; only 8 or 16 are supported.
        `sample_rate` : int
            Samples per second (in Hertz).

    '''

    def __init__(self, channels, sample_size, sample_rate):
        self.channels = channels
        self.sample_size = sample_size
        self.sample_rate = sample_rate
        
        # Convenience
        self.bytes_per_sample = (sample_size >> 3) * channels
        self.bytes_per_second = self.bytes_per_sample * sample_rate

    def __eq__(self, other):
        return (self.channels == other.channels and 
                self.sample_size == other.sample_size and
                self.sample_rate == other.sample_rate)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return '%s(channels=%d, sample_size=%d, sample_rate=%d)' % (
            self.__class__.__name__, self.channels, self.sample_size,
            self.sample_rate)

class VideoFormat(object):
    '''Video details.

    An instance of this class is provided by sources with a video track.  You
    should not modify the fields.

    Note that the sample aspect has no relation to the aspect ratio of the
    video image.  For example, a video image of 640x480 with sample aspect 2.0
    should be displayed at 1280x480.  It is the responsibility of the
    application to perform this scaling.

    :Ivariables:
        `width` : int
            Width of video image, in pixels.
        `height` : int
            Height of video image, in pixels.
        `sample_aspect` : float
            Aspect ratio (width over height) of a single video pixel.

    '''
    
    def __init__(self, width, height, sample_aspect=1.0):
        self.width = width
        self.height = height
        self.sample_aspect = sample_aspect

class AudioData(object):
    '''A single packet of audio data.

    This class is used internally by pyglet.

    :Ivariables:
        `data` : str or ctypes array or pointer
            Sample data.
        `length` : int
            Size of sample data, in bytes.
        `timestamp` : float
            Time of the first sample, in seconds.
        `duration` : float
            Total data duration, in seconds.

    '''
    def __init__(self, data, length, timestamp, duration):
        self.data = data
        self.length = length
        self.timestamp = timestamp
        self.duration = duration

    def consume(self, bytes, audio_format):
        '''Remove some data from beginning of packet.'''
        if bytes == self.length:
            self.data = None
            self.length = 0
            self.timestamp += self.duration
            self.duration = 0.
            return
        elif bytes == 0:
            return

        if not isinstance(self.data, str):
            # XXX Create a string buffer for the whole packet then
            #     chop it up.  Could do some pointer arith here and
            #     save a bit of data pushing, but my guess is this is
            #     faster than fudging aruond with ctypes (and easier).
            data = ctypes.create_string_buffer(self.length)
            ctypes.memmove(data, self.data, self.length)
            self.data = data
        self.data = self.data[bytes:]
        self.length -= bytes
        self.duration -= bytes / float(audio_format.bytes_per_second)
        self.timestamp += bytes / float(audio_format.bytes_per_second)

    def get_string_data(self):
        '''Return data as a string.'''
        if type(self.data) is str:
            return self.data

        buf = ctypes.create_string_buffer(self.length)
        ctypes.memmove(buf, self.data, self.length)
        return buf.raw

class AudioPlayer(object):
    '''Abstract low-level interface for playing audio.

    AudioPlayer has no knowledge of sources or eos behaviour.  Once
    created, its audio format cannot be modified.  The player will attempt
    to recover automatically from a buffer underrun (but this is not
    guaranteed).
    
    Applications should not use this class directly, but instead use `Player`.

    :Ivariables:
        `audio_format` : `AudioFormat`
            The player's audio format (read-only).

    '''

    UPDATE_PERIOD = 0.15
    
    def __init__(self, audio_format):
        '''Create a new audio player for the given audio format.

        :Parameters:
            `audio_format` : `AudioFormat`
                Audio format parameters.

        '''
        self.audio_format = audio_format

    def get_write_size(self):
        '''Return the maximum number of bytes that can be written.  
        
        This is used as a hint for preparing data for `write`, not as a strict
        contract.
        
        :rtype: int
        '''
        raise NotImplementedError('abstract')

    def write(self, audio_data):
        '''Write audio_data to the stream.

        This method calls `AudioData.consume` to remove data actually written.
        
        :Parameters:
            `audio_data` : `AudioData`
                Data to write.

        '''
        raise NotImplementedError('abstract')

    def write_eos(self):
        '''Write an EOS marker to the stream at the current write point.'''
        raise NotImplementedError('abstract')

    def write_end(self):
        '''Mark that there will be no more audio data past the current write
        point.
        
        This does not produce an EOS, but is required to prevent data
        underrun artifacts.
        '''
        raise NotImplementedError('abstract')

    def play(self):
        '''Begin playback.'''
        raise NotImplementedError('abstract')

    def stop(self):
        '''Stop playback.'''
        raise NotImplementedError('abstract')

    def clear(self):
        '''Clear all buffered data and prepare for replacement data.

        The player should be stopped before calling this method.
        '''
        raise NotImplementedError('abstract')

    def pump(self):
        '''Called once per loop iteration before checking for eos
        triggers.'''
        raise NotImplementedError('abstract')

    def get_time(self):
        '''Return best guess of current playback time.  The time is relative
        to the timestamps provided in the data supplied to `write`.  The time
        is meaningless unless proper care has been taken to clear EOS markers.

        :rtype: float
        :return: current play cursor time, in seconds.
        '''
        raise NotImplementedError('abstract')

    def clear_eos(self):
        '''Check if an EOS marker has been passed, and clear it.
        
        This method should be called repeatedly to clear all pending EOS
        markers.

        :rtype: bool
        :return: True if an EOS marker was cleared.
        '''
        raise NotImplementedError('abstract')

    def set_volume(self, volume):
        '''See `Player.volume`.'''
        pass

    def set_position(self, position):
        '''See `Player.position`.'''
        pass

    def set_min_distance(self, min_distance):
        '''See `Player.min_distance`.'''
        pass

    def set_max_distance(self, max_distance):
        '''See `Player.max_distance`.'''
        pass

    def set_pitch(self, pitch):
        '''See `Player.pitch`.'''
        pass

    def set_cone_orientation(self, cone_orientation):
        '''See `Player.cone_orientation`.'''
        pass

    def set_cone_inner_angle(self, cone_inner_angle):
        '''See `Player.cone_inner_angle`.'''
        pass

    def set_cone_outer_angle(self, cone_outer_angle):
        '''See `Player.cone_outer_angle`.'''
        pass

    def set_cone_outer_gain(self, cone_outer_gain):
        '''See `Player.cone_outer_gain`.'''
        pass

class Source(object):
    '''An audio and/or video source.

    :Ivariables:
        `audio_format` : `AudioFormat`
            Format of the audio in this source, or None if the source is
            silent.
        `video_format` : `VideoFormat`
            Format of the video in this source, or None if there is no
            video.
    '''

    _duration = None
    
    audio_format = None
    video_format = None

    def _get_duration(self):
        return self._duration

    duration = property(lambda self: self._get_duration(),
                        doc='''The length of the source, in seconds.

        Not all source durations can be determined; in this case the value
        is None.

        Read-only.

        :type: float
        ''')

    def play(self):
        '''Play the source.

        This is a convenience method which creates a ManagedSoundPlayer for
        this source and plays it immediately.

        :rtype: `ManagedSoundPlayer`
        '''
        player = ManagedSoundPlayer()
        player.queue(self)
        player.play()
        return player

    def get_animation(self):
        '''Import all video frames into memory as an `Animation`.

        An empty animation will be returned if the source has no video.
        Otherwise, the animation will contain all unplayed video frames (the
        entire source, if it has not been queued on a player).  After creating
        the animation, the source will be at EOS.

        This method is unsuitable for videos running longer than a
        few seconds.

        :since: pyglet 1.1

        :rtype: `pyglet.image.Animation`
        '''
        from pyglet.image import Animation, AnimationFrame
        if not self.video_format:
            return Animation([])
        else:
            # Create a dummy player for the source to push its textures onto.
            frames = []
            last_ts = 0
            next_ts = self.get_next_video_timestamp()
            while next_ts is not None:
                image = self.get_next_video_frame()
                assert image is not None
                delay = next_ts - last_ts
                frames.append(AnimationFrame(image, delay))
                last_ts = next_ts
                next_ts = self.get_next_video_timestamp()
            return Animation(frames)

    def get_next_video_timestamp(self):
        '''Get the timestamp of the next video frame.

        :since: pyglet 1.1

        :rtype: float
        :return: The next timestamp, or ``None`` if there are no more video
            frames.
        '''
        pass

    def get_next_video_frame(self):
        '''Get the next video frame.

        Video frames may share memory: the previous frame may be invalidated
        or corrupted when this method is called unless the application has
        made a copy of it.

        :since: pyglet 1.1

        :rtype: `pyglet.image.AbstractImage`
        :return: The next video frame image, or ``None`` if there are no more
            video frames.
        '''
        pass

    # Internal methods that Players call on the source:

    def _play(self):
        '''Begin decoding in real-time.'''
        pass

    def _pause(self):
        '''Pause decoding, but remain prerolled.'''
        pass

    def _stop(self):
        '''Stop forever and clean up.'''
        pass

    def _seek(self, timestamp):
        '''Seek to given timestamp.'''
        raise CannotSeekException()

    def _get_queue_source(self):
        '''Return the `Source` to be used as the queue source for a player.

        Default implementation returns self.'''
        return self

    def _get_audio_data(self, bytes):
        '''Get next packet of audio data.

        :Parameters:
            `bytes` : int
                Maximum number of bytes of data to return.

        :rtype: `AudioData`
        :return: Next packet of audio data, or None if there is no (more)
            data.
        '''
        return None

    def _init_texture(self, player):
        '''Create the player's texture.'''
        pass
    
    def _update_texture(self, player, timestamp):
        '''Update the texture on player.'''
        pass

    def _release_texture(self, player):
        '''Release the player's texture.'''
        pass

class StreamingSource(Source):
    '''A source that is decoded as it is being played, and can only be
    queued once.
    '''
    
    _is_queued = False

    is_queued = property(lambda self: self._is_queued,
                         doc='''Determine if this source has been queued
        on a `Player` yet.

        Read-only.

        :type: bool
        ''')

    def _get_queue_source(self):
        '''Return the `Source` to be used as the queue source for a player.

        Default implementation returns self.'''
        if self._is_queued:
            raise MediaException('This source is already queued on a player.')
        self._is_queued = True
        return self

class StaticSource(Source):
    '''A source that has been completely decoded in memory.  This source can
    be queued onto multiple players any number of times.
    '''
    
    def __init__(self, source):
        '''Construct a `StaticSource` for the data in `source`.

        :Parameters:
            `source` : `Source`
                The source to read and decode audio and video data from.

        '''
        source = source._get_queue_source()
        if source.video_format:
            raise NotImplementedException(
                'Static sources not supported for video yet.')

        self.audio_format = source.audio_format
        if not self.audio_format:
            return

        # TODO enable time-insensitive playback 
        source._play()

        # Arbitrary: number of bytes to request at a time.
        buffer_size = 1 << 20 # 1 MB

        # Naive implementation.  Driver-specific implementations may override
        # to load static audio data into device (or at least driver) memory. 
        data = StringIO.StringIO()
        while True:
            audio_data = source._get_audio_data(buffer_size)
            if not audio_data:
                break
            data.write(audio_data.get_string_data())
        self._data = data.getvalue()

    def _get_queue_source(self):
        return StaticMemorySource(self._data, self.audio_format)

    def _get_audio_data(self, bytes):
        raise RuntimeError('StaticSource cannot be queued.')

class StaticMemorySource(StaticSource):
    '''Helper class for default implementation of `StaticSource`.  Do not use
    directly.'''

    def __init__(self, data, audio_format):
        '''Construct a memory source over the given data buffer.
        '''
        self._file = StringIO.StringIO(data)
        self._max_offset = len(data)
        self.audio_format = audio_format
        self._duration = len(data) / float(audio_format.bytes_per_second)

    def _seek(self, timestamp):
        offset = int(timestamp * self.audio_format.bytes_per_second)

        # Align to sample
        if self.audio_format.bytes_per_sample == 2:
            offset &= 0xfffffffe
        elif self.audio_format.bytes_per_sample == 4:
            offset &= 0xfffffffc

        self._file.seek(offset)

    def _get_audio_data(self, bytes):
        offset = self._file.tell()
        timestamp = float(offset) / self.audio_format.bytes_per_second

        # Align to sample size
        if self.audio_format.bytes_per_sample == 2:
            bytes &= 0xfffffffe
        elif self.audio_format.bytes_per_sample == 4:
            bytes &= 0xfffffffc

        data = self._file.read(bytes)
        if not len(data):
            return None

        duration = float(len(data)) / self.audio_format.bytes_per_second
        return AudioData(data, len(data), timestamp, duration)

class Player(event.EventDispatcher):
    '''A sound and/or video player.

    Queue sources on this player to play them.
    '''

    #: The player will pause when it reaches the end of the stream.
    EOS_PAUSE = 'pause'
    #: The player will loop the current stream continuosly.
    EOS_LOOP = 'loop'
    #: The player will move on to the next queued stream when it reaches the
    #: end of the current source.  If there is no source queued, the player
    #: will pause.
    EOS_NEXT = 'next'
    #: The player will stop entirely; valid only for ManagedSoundPlayer.
    EOS_STOP = 'stop'

    # Source and queuing attributes
    _source_read_index = 0
    _eos_action = EOS_NEXT
    _playing = False

    # If True and _playing is False, user is currently seeking while paused;
    # should refrain from filling the audio buffer.
    _pause_seek = False

    # Override audio timestamp for seeking and silent video
    _timestamp = None

    # Used to track timestamp for silent sources
    _last_system_time = 0.

    # Audio attributes
    _audio = None
    _audio_finished = False
    _next_audio_data = None

    # Video attributes
    _texture = None

    # Spacialisation attributes, preserved between audio players
    _volume = 1.0
    _min_distance = 1.0
    _max_distance = 100000000.

    _position = (0, 0, 0)
    _pitch = 1.0

    _cone_orientation = (0, 0, 1)
    _cone_inner_angle = 360.
    _cone_outer_angle = 360.
    _cone_outer_gain = 1.

    def __init__(self):
        self._sources = []

    def _create_audio(self):
        '''Create _audio for sources[0].  
        
        Reuses existing _audio if it exists and is compatible.
        '''
        if not self._sources:
            return

        source = self._sources[0]
        if not source.audio_format:
            self._audio = None
            return

        if self._audio:
            self._audio_finished = False
            if self._audio.audio_format == source.audio_format:
                return
            else:
                self._audio = None

        self._audio = audio_player_class(source.audio_format)
        self._audio.set_volume(self._volume)
        self._audio.set_min_distance(self._min_distance)
        self._audio.set_max_distance(self._max_distance)
        self._audio.set_position(self._position)
        self._audio.set_pitch(self._pitch)
        self._audio.set_cone_orientation(self._cone_orientation)
        self._audio.set_cone_inner_angle(self._cone_inner_angle)
        self._audio.set_cone_outer_angle(self._cone_outer_angle)
        self._audio.set_cone_outer_gain(self._cone_outer_gain)

    def _fill_audio(self):
        '''Ensure _audio is full.'''
        if not self._audio or self._audio_finished:
            return

        write_size = self._audio.get_write_size()
        if not write_size:
            return

        for audio_data, audio_format in self._get_audio_data(write_size):
            if audio_data == 'eos':
                self._audio.write_eos()
                continue
            elif audio_data == 'end':
                self._audio.write_end()
                self._audio_finished = True
                return
            if audio_format != self._audio.audio_format:
                self._next_audio_data = audio_format, audio_data
                return
            length = audio_data.length
            self._audio.write(audio_data)
            if audio_data.length:
                self._next_audio_data = audio_format, audio_data
                return

            write_size -= length
            if write_size <= 0:
                return

    def _get_audio_data(self, bytes):
        '''Yields pairs of (audio_data, audio_format).'''
        if self._next_audio_data:
            audio_format, audio_data = self._next_audio_data
            self._next_audio_data = None
            bytes -= audio_data.length
            yield audio_data, audio_format

        try:
            source = self._sources[self._source_read_index]
        except IndexError:
            source = None

        while source and bytes > 4: # bytes > 4 compensates for alignment loss
            audio_data = source._get_audio_data(bytes)
            if audio_data:
                bytes -= audio_data.length
                yield audio_data, source.audio_format
            else:
                yield 'eos', source.audio_format
                if self._eos_action == self.EOS_NEXT:
                    self._source_read_index += 1
                    try:
                        source = self._sources[self._source_read_index]
                        source._play()
                    except IndexError:
                        source = None
                elif self._eos_action == self.EOS_LOOP:
                    source._seek(0)
                elif self._eos_action == self.EOS_PAUSE:
                    source = None
                elif self._eos_action == self.EOS_STOP:
                    source = None
                else:
                    assert False, 'Invalid eos_action'
                    source = None

        if not source:
            yield 'end', None

    def _update_schedule(self):
        clock.unschedule(self.dispatch_events)
        if self._playing and self._sources:
            interval = 1000.
            if self._sources[0].video_format:
                interval = min(interval, 1/24.)
            if self._audio:
                interval = min(interval, self._audio.UPDATE_PERIOD)
            clock.schedule_interval_soft(self.dispatch_events, interval)

    def queue(self, source):
        '''Queue the source on this player.

        If the player has no source, the player will be paused immediately
        on this source.

        :Parameters:
            `source` : Source
                The source to queue.

        '''
        self._sources.append(source._get_queue_source())
        if len(self._sources) == 1:
            self._source_read_index = 0
            self._begin_source()

    def play(self):
        '''Begin playing the current source.

        This has no effect if the player is already playing.
        '''
        self._playing = True
        self._pause_seek = False

        if self._audio:
            self._timestamp = None
            self._audio.play()
        else:
            self._last_system_time = time.time()

        self.dispatch_events()
        self._update_schedule()

    def pause(self):
        '''Pause playback of the current source.

        This has no effect if the player is already paused.
        '''
        self._playing = False
        self._pause_seek = False
        
        if self._audio:
            self._audio.stop()
        self._update_schedule()

    def seek(self, timestamp):
        '''Seek for playback to the indicated timestamp in seconds on the
        current source.  If the timestamp is outside the duration of the
        source, it will be clamped to the end.

        :Parameters:
            `timestamp` : float
                Timestamp to seek to.
        '''
        if not self._sources:
            return

        if not self._playing:
            self._pause_seek = True

        self._audio_finished = False
        source = self._sources[0]
        self._source_read_index = 0
        self._next_audio_data = None
        source._seek(timestamp)
        self._timestamp = timestamp

        if self._audio:
            self._audio.stop()
            self._audio.clear()
        else:
            self._last_system_time = time.time()

        self.dispatch_events()
        
    def next(self):
        '''Move immediately to the next queued source.

        There may be a gap in playback while the audio buffer is refilled.
        '''
        if not self._sources:
            return

        if self._audio:
            self._audio.stop()
            self._audio.clear()
        else:
            self._last_system_time = time.time()
            self._timestamp = 0.

        self._next_source()

    def _next_source(self):
        if not self._sources:
            self._update_schedule()
            return

        self._source_read_index = max(0, self._source_read_index - 1)
        source = self._sources.pop(0)
        source._release_texture(self)
        source._stop()
        self._begin_source()

    def _begin_source(self):
        if not self._sources:
            return

        source = self._sources[0]
        source._init_texture(self)
        self._create_audio()
        self._fill_audio()

        if not self._audio:
            self._timestamp = 0.

        if self._playing:
            self.play()
            self._update_schedule()

    def _on_eos(self):
        '''Internal method when EOS is encountered.  Returns False if
        dispatch_events should be immediately aborted.'''
        if self._eos_action == self.EOS_NEXT:
            self._next_source()
        elif self._eos_action == self.EOS_PAUSE:
            self._playing = False
            self._timestamp = self._sources[0].duration
        elif self._eos_action == self.EOS_STOP:
            self.stop()
            self._sources = []
            return False
        self.dispatch_event('on_eos')
        return True

    def dispatch_events(self, dt=None):
        '''Dispatch any pending events and perform regular heartbeat functions
        to maintain playback.

        :Parameters:
            `dt` : None
                Ignored (for compatibility with `pyglet.clock.schedule`)

        :deprecated: Since pyglet 1.1, Player objects schedule themselves on
            the default clock automatically.  Applications should not call
            this method.

        '''
        if not self._sources:
            return

        if not self._pause_seek:
            self._fill_audio()

        if self._audio:
            underrun = self._audio.pump()
            while self._audio.clear_eos():
                if not self._on_eos():
                    return
            if underrun:
                self._audio.UPDATE_PERIOD *= 0.75
                self._audio.__class__.UPDATE_PERIOD *= 0.75
                self._update_schedule()
                if _debug_media:
                    print '%r underrun: reducing update period to %.2f' % \
                        (self._audio, self._audio.UPDATE_PERIOD)
        else:
            if self._playing:
                t = time.time()
                self._timestamp += t - self._last_system_time
                self._last_system_time = t
                while self._timestamp > self._sources[0].duration:
                    if not self._on_eos():
                        return
                    if self._eos_action == self.EOS_LOOP:
                        self._timestamp -= self._sources[0].duration

        if self._texture:
            self._sources[0]._update_texture(self, self._get_time())

    def _get_time(self):
        if self._timestamp is not None:
            return self._timestamp
        elif self._audio:
            return self._audio.get_time()

    time = property(lambda self: self._get_time(),
                    doc='''Retrieve the current playback time of the current
         source.
                    
         The playback time is a float expressed in seconds, with 0.0 being
         the beginning of the sound.  The playback time returned represents
         the time encoded in the source, and may not reflect actual time
         passed due to pitch shifting or pausing.

         Read-only.

         :type: float
         ''')

    def _get_source(self):
        if self._sources:
            return self._sources[0]

    source = property(lambda self: self._get_source(),
                      doc='''Return the current source.

         Read-only.

         :type: Source
         ''')

    def _set_eos_action(self, action):
        self._eos_action = action

    eos_action = property(lambda self: self._eos_action,
                          _set_eos_action,
                          doc='''Set the behaviour of the player when it
        reaches the end of the current source.

        This must be one of the constants `EOS_NEXT`, `EOS_PAUSE` or
        `EOS_LOOP`.

        :type: str
        ''')

    playing = property(lambda self: self._playing,
                       doc='''Determine if the player state is playing.

        The `playing` property is irrespective of whether or not there is
        actually a source to play.  If `playing` is True and a source is
        queued, it will begin playing immediately.  If `playing` is False, 
        it is implied that the player is paused.  There is no other possible
        state.

        Read-only.

        :type: bool
        ''')

    def _set_volume(self, volume):
        self._volume = volume
        if self._audio:
            self._audio.set_volume(volume)

    volume = property(lambda self: self._volume,
                      lambda self, volume: self._set_volume(volume),
                      doc='''The volume level of sound playback.

         The nominal level is 1.0, and 0.0 is silence.

         The volume level is affected by the distance from the listener (if
         positioned).
         
         :type: float
         ''')

    def _set_position(self, position):
        self._position = position
        if self._audio:
            self._audio.set_position(position)

    position = property(lambda self: self._position,
                        lambda self, position: self._set_position(position),
                        doc='''The position of the sound in 3D space.

        The position is given as a tuple of floats (x, y, z).  The unit
        defaults to meters, but can be modified with the listener
        properties.
        
        :type: 3-tuple of float
        ''')

    def _set_min_distance(self, min_distance):
        self._min_distance = min_distance
        if self._audio:
            self._audio.set_min_distance(min_distance)

    min_distance = property(lambda self: self._min_distance,
                            lambda self, v: self._set_min_distance(v),
                            doc='''The distance beyond which the sound volume
        drops by half, and within which no attenuation is applied.

        The minimum distance controls how quickly a sound is attenuated
        as it moves away from the listener.  The gain is clamped at the
        nominal value within the min distance.  By default the value is
        1.0.
        
        The unit defaults to meters, but can be modified with the listener
        properties.
        
        :type: float
        ''')

    def _set_max_distance(self, max_distance):
        self._max_distance = max_distance
        if self._audio:
            self._audio.set_max_distance(max_distance)

    max_distance = property(lambda self: self._max_distance,
                            lambda self, v: self._set_max_distance(v),
                            doc='''The distance at which no further attenuation
        is applied.

        When the distance from the listener to the player is greater than 
        this value, attenuation is calculated as if the distance
        were value.  By default the maximum distance is infinity.
        
        The unit defaults to meters, but can be modified with the listener
        properties.
        
        :type: float
        ''')

    def _set_pitch(self, pitch):
        self._pitch = pitch
        if self._audio:
            self._audio.set_pitch(pitch)

    pitch = property(lambda self: self._pitch,
                     lambda self, pitch: self._set_pitch(pitch),
                     doc='''The pitch shift to apply to the sound.

        The nominal pitch is 1.0.  A pitch of 2.0 will sound one octave
        higher, and play twice as fast.  A pitch of 0.5 will sound one octave
        lower, and play twice as slow.  A pitch of 0.0 is not permitted.
        
        :type: float
        ''')

    def _set_cone_orientation(self, cone_orientation):
        self._cone_orientation = cone_orientation
        if self._audio:
            self._audio.set_cone_orientation(cone_orientation)

    cone_orientation = property(lambda self: self._cone_orientation,
                                lambda self, c: self._set_cone_orientation(c),
                                doc='''The direction of the sound in 3D space.
                                
        The direction is specified as a tuple of floats (x, y, z), and has no
        unit.  The default direction is (0, 0, -1).  Directional effects are
        only noticeable if the other cone properties are changed from their
        default values.
        
        :type: 3-tuple of float
        ''')

    def _set_cone_inner_angle(self, cone_inner_angle):
        self._cone_inner_angle = cone_inner_angle
        if self._audio:
            self._audio.set_cone_inner_angle(cone_inner_angle)

    cone_inner_angle = property(lambda self: self._cone_inner_angle,
                                lambda self, a: self._set_cone_inner_angle(a),
                                doc='''The interior angle of the inner cone.
                                
        The angle is given in degrees, and defaults to 360.  When the listener
        is positioned within the volume defined by the inner cone, the sound
        is played at normal gain (see `volume`).
        
        :type: float
        ''')

    def _set_cone_outer_angle(self, cone_outer_angle):
        self._cone_outer_angle = cone_outer_angle
        if self._audio:
            self._audio.set_cone_outer_angle(cone_outer_angle)

    cone_outer_angle = property(lambda self: self._cone_outer_angle,
                                lambda self, a: self._set_cone_outer_angle(a),
                                doc='''The interior angle of the outer cone.
                                
        The angle is given in degrees, and defaults to 360.  When the listener
        is positioned within the volume defined by the outer cone, but outside
        the volume defined by the inner cone, the gain applied is a smooth
        interpolation between `volume` and `cone_outer_gain`.
        
        :type: float
        ''')

    def _set_cone_outer_gain(self, cone_outer_gain):
        self._cone_outer_gain = cone_outer_gain
        if self._audio:
            self._audio.set_cone_outer_gain(cone_outer_gain)

    cone_outer_gain = property(lambda self: self._cone_outer_gain,
                                lambda self, g: self._set_cone_outer_gain(g),
                                doc='''The gain applied outside the cone.
                                
        When the listener is positioned outside the volume defined by the
        outer cone, this gain is applied instead of `volume`.
        
        :type: float
        ''')

    def get_texture(self):
        '''Get the texture for the current video frame.

        You should call this method every time you display a frame
        of video, as multiple textures might be used.  The return value will
        be `None` if there is no video in the current source.

        :since: pyglet 1.1

        :rtype: `pyglet.image.Texture`
        '''
        return self._texture

    texture = property(lambda self: self._texture,
                       doc='''The video texture.

        You should rerequest this property every time you display a frame
        of video, as multiple textures might be used.  This property will
        be `None` if there is no video in the current source.

        :deprecated: Use `get_texture`.

        :type: `pyglet.image.Texture`
        ''')

    if getattr(sys, 'is_epydoc', False):
        def on_eos():
            '''The player has reached the end of the current source.

            This event is dispatched regardless of the EOS action.  You
            can alter the EOS action in this event handler, however playback
            may stutter as the media device will not have enough time to
            decode and buffer the new data in advance.

            :event:
            '''
Player.register_event_type('on_eos')

class ManagedSoundPlayer(Player):
    '''A player which takes care of updating its own audio buffers.

    This player will continue playing the sound until the sound is
    finished, even if the application discards the player early.

    Only one source can be queued on the player; the player will be 
    discarded when the source finishes.
    '''

    #: The only possible end of stream action for a managed player.
    EOS_STOP = 'stop'

    _eos_action = EOS_STOP
    eos_action = property(lambda self: EOS_STOP,
                          doc='''The fixed eos_action is `EOS_STOP`,
        in which the player is discarded as soon as the source has
        finished.

        Read-only.
        
        :type: str
        ''')

    def __init__(self):
        super(ManagedSoundPlayer, self).__init__()
        managed_players.append(self)

    def stop(self):
        self._timestamp = 0.
        clock.unschedule(self.dispatch_events)
        managed_players.remove(self)

class Listener(object):
    '''The listener properties for positional audio.

    You can obtain the singleton instance of this class as
    `pyglet.media.listener`.
    '''
    _volume = 1.0
    _position = (0, 0, 0)
    _forward_orientation = (0, 0, -1)
    _up_orientation = (0, 1, 0)

    def _set_volume(self, volume):
        raise NotImplementedError('abstract')

    volume = property(lambda self: self._volume,
                      lambda self, volume: self._set_volume(volume),
                      doc='''The master volume for sound playback.

        All sound volumes are multiplied by this master volume before being
        played.  A value of 0 will silence playback (but still consume
        resources).  The nominal volume is 1.0.
        
        :type: float
        ''')

    def _set_position(self, position):
        raise NotImplementedError('abstract')

    position = property(lambda self: self._position,
                        lambda self, position: self._set_position(position),
                        doc='''The position of the listener in 3D space.

        The position is given as a tuple of floats (x, y, z).  The unit
        defaults to meters, but can be modified with the listener
        properties.
        
        :type: 3-tuple of float
        ''')

    def _set_forward_orientation(self, orientation):
        raise NotImplementedError('abstract')

    forward_orientation = property(lambda self: self._forward_orientation,
                               lambda self, o: self._set_forward_orientation(o),
                               doc='''A vector giving the direction the
        listener is facing.

        The orientation is given as a tuple of floats (x, y, z), and has
        no unit.  The forward orientation should be orthagonal to the
        up orientation.
        
        :type: 3-tuple of float
        ''')

    def _set_up_orientation(self, orientation):
        raise NotImplementedError('abstract')

    up_orientation = property(lambda self: self._up_orientation,
                              lambda self, o: self._set_up_orientation(o),
                              doc='''A vector giving the "up" orientation
        of the listener.

        The orientation is given as a tuple of floats (x, y, z), and has
        no unit.  The up orientation should be orthagonal to the
        forward orientation.
        
        :type: 3-tuple of float
        ''')


if getattr(sys, 'is_epydoc', False):
    #: The singleton listener.
    #:
    #: :type: `Listener`
    listener = Listener()

    #: Indication of the presence of AVbin.  When `have_avbin` is ``True``
    #: pyglet will be able to play back compressed media streams such as
    #: MP3, OGG and various video formats.  If ``False`` only uncompressed
    #: Wave files can be loaded.
    #:
    #: :type: bool
    have_avbin = False
else:
    # Find best available sound driver according to user preference
    import pyglet
    driver = None
    for driver_name in pyglet.options['audio']:
        try:
            driver_name = 'pyglet.media.drivers.' + driver_name
            __import__(driver_name)
            driver = sys.modules[driver_name]
            driver.driver_init()
            break
        except (ImportError, AttributeError, MediaException):
            pass

    if not driver:
        raise ImportError('No suitable audio driver could be loaded.')

    audio_player_class = driver.driver_audio_player_class
    listener = driver.driver_listener

    # Find best available source loader
    have_avbin = False
    try:
        from pyglet.media import avbin
        _source_class = avbin.AVbinSource
        have_avbin = True
    except ImportError:
        from pyglet.media import riff
        _source_class = riff.WaveSource

# Pretend to import some common audio drivers so that py2exe/py2app
# are fooled into packagin them.
if False:
    import pyglet.media.drivers.silent
    import pyglet.media.drivers.openal
    import pyglet.media.drivers.directsound
    import pyglet.media.drivers.alsa

def load(filename, file=None, streaming=True):
    '''Load a source from a file.

    Currently the `file` argument is not supported; media files must exist
    as real paths.

    :Parameters:
        `filename` : str
            Filename of the media file to load.
        `file` : file-like object
            Not yet supported.
        `streaming` : bool
            If False, a `StaticSource` will be returned; otherwise (default) a
            `StreamingSource` is created.

    :rtype: `Source`
    '''
    source = _source_class(filename, file)
    if not streaming:
        source = StaticSource(source)
    return source

managed_players = []
def dispatch_events():
    '''Process managed audio events.

    You must call this function regularly (typically once per run loop
    iteration) in order to keep audio buffers of managed players full.

    :deprecated: Since pyglet 1.1, Player objects schedule themselves on
        the default clock automatically.  Applications should not call this
        method.

    '''
    for player in managed_players:
        player.dispatch_events()
