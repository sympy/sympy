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
# $Id: __init__.py 1222 2007-09-01 11:00:40Z Alex.Holkner $

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

Use the `Player` to control playback.  Within your main run loop, you must
periodically call `dispatch_events` to ensure the audio buffers are refilled::

    player.play()
    while player.source:    # While the source hasn't finished
        player.dispatch_events()

If the source contains video, its `video_format` attribute will be non-None,
and the player's `texture` attribute will contain the current video image
synchronised to the audio.

Decoding sounds can be processor-intensive and may introduce latency,
particularly for short sounds that must be played quickly, such as bullets or
explosions.  You can force such sounds to be decoded and retained in memory
rather than streamed from disk by wrapping the source in a `StaticSource`::

    bullet_sound = StaticSource(load('bullet.wav'))

The other advantage of a `StaticSound` is that it can be queued on any number
of players, and so played many times simultaneously.

'''

__docformat__ = 'restructuredtext'
__version__ = '$Id: __init__.py 1222 2007-09-01 11:00:40Z Alex.Holkner $'

import ctypes
import sys
import StringIO

from pyglet import event

class MediaException(Exception):
    pass

class MediaFormatException(Exception):
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
            Bits per sample; typically 8 or 16.
        `sample_rate` : int
            Samples per second (in Herz).

    '''

    def __init__(self, channels, sample_size, sample_rate):
        self.channels = channels
        self.sample_size = sample_size
        self.sample_rate = sample_rate
        
        # Convenience
        self.bytes_per_sample = (sample_size >> 3) * channels
        self.bytes_per_second = self.bytes_per_sample * sample_rate

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
        `is_eos` : bool
            If True, this is the last audio packet in the source.

    '''
    def __init__(self, data, length, timestamp, duration, is_eos=False):
        self.data = data
        self.length = length
        self.timestamp = timestamp
        self.duration = duration
        self.is_eos = is_eos

    def consume(self, bytes, audio_format):
        '''Remove some data from beginning of packet.'''
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
        player.eos_action = player.EOS_STOP
        player.queue(self)
        player.play()
        return player

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
        if source.video_format:
            raise NotImplementedException(
                'Static sources not supported for video yet.')

        self.audio_format = source.audio_format
        if not self.audio_format:
            return

        # TODO enable time-insensitive playback 
        source.play()

        # Arbitrary: number of bytes to request at a time.
        buffer_size = 1 << 20 # 1 MB

        # Naive implementation.  Driver-specific implementations may override
        # to load static audio data into device (or at least driver) memory. 
        data = StringIO.StringIO()
        while True:
            audio_data = source._get_audio_data(buffer_size)
            if not audio_data:
                break
            data.write(audio_data.data)
        self._data = data.getvalue()

    def _get_queue_source(self):
        return StaticMemorySource(self._data, self.audio_format)

    def _get_audio_data(self, bytes):
        raise RuntimeError('StaticSource cannot be queued.')

class StaticMemorySource(StaticSource):
    '''Helper class for default implementation of `StaticSource`.  Do not use
    directly.'''

    def __init__(self, data, audio_format):
        self._file = StringIO.StringIO(data)
        self._max_offset = len(data)
        self.audio_format = audio_format
        self._duration = len(data) / float(audio_format.bytes_per_second)

    def _seek(self, timestamp):
        offset = int(timestamp * self.audio_format.bytes_per_second)

        # Align to sample
        if self.audio_format.bytes_per_sample == 2:
            offset &= 0xfffffffe
        elif self.audio_foramt.bytes_per_sample == 4:
            offset &= 0xfffffffc

        self._file.seek(offset)

    def _get_audio_data(self, bytes):
        offset = self._file.tell()
        timestamp = float(offset) / self.audio_format.bytes_per_second

        data = self._file.read(bytes)
        if not data:
            return None

        duration = float(len(data)) / self.audio_format.bytes_per_second
        is_eos = self._file.tell() == self._max_offset
        return AudioData(data,
                         len(data),
                         timestamp,
                         duration,
                         is_eos)

class BasePlayer(event.EventDispatcher):
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
    _source = None
    _eos_action = EOS_NEXT
    _playing = False

    # Sound and spacialisation attributes
    _volume = 1.0
    _max_gain = 1.0
    _min_gain = 0.0

    _position = (0, 0, 0)
    _velocity = (0, 0, 0)
    _pitch = 1.0

    _cone_orientation = (0, 0, 0)
    _cone_inner_angle = 360.
    _cone_outer_angle = 360.
    _cone_outer_gain = 1.

    # Video attributes
    _texture = None

    def queue(self, source):
        '''Queue the source on this player.

        If the player has no source, the player will be paused immediately
        on this source.

        :Parameters:
            `source` : Source
                The source to queue.

        '''

    def play(self):
        '''Begin playing the current source.

        This has no effect if the player is already playing.
        '''
        raise NotImplementedError('abstract')

    def pause(self):
        '''Pause playback of the current source.

        This has no effect if the player is already paused.
        '''
        raise NotImplementedError('abstract')

    def seek(self, timestamp):
        '''Seek for playback to the indicated timestamp in seconds on the
        current source.  If the timestamp is outside the duration of the
        source, it will be clamped to the end.

        :Parameters:
            `timestamp` : float
                Timestamp to seek to.
        '''
        raise NotImplementedError('abstract')

    def next(self):
        '''Move immediately to the next queued source.

        If the `eos_action` of this player is `EOS_NEXT`, and the source has
        been queued for long enough, there will be no gap in the audio or
        video playback.  Otherwise, there may be some delay as the next source
        is prerolled and the first frames decoded and buffered.
        '''
        raise NotImplementedError('abstract')

    def dispatch_events(self):
        '''Dispatch any pending events and perform regular heartbeat functions
        to maintain playback.
        '''
        pass

    def _get_time(self):
        raise NotImplementedError('abstract')

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
        return self._source

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
        raise NotImplementedError('abstract')

    volume = property(lambda self: self._volume,
                      lambda self, volume: self._set_volume(volume),
                      doc='''The volume level of sound playback.

         The nominal level is 1.0, and 0.0 is silence.

         The volume level is affected by factors such as the distance from the
         listener (if positioned), and is clamped (after distance attenuation)
         to the range [min_gain, max_gain].
         
         :type: float
         ''')

    def _set_min_gain(self, min_gain):
        raise NotImplementedError('abstract')

    min_gain = property(lambda self: self._min_gain,
                        lambda self, min_gain: self._set_min_gain(min_gain),
                        doc='''The minimum gain to apply to the sound, even

         The gain is clamped after distance attenuation.  The default value
         is 0.0.
         
         :type: float
         ''')

    def _set_max_gain(self, max_gain):
        raise NotImplementedError('abstract')

    max_gain = property(lambda self: self._max_gain,
                        lambda self, max_gain: self._set_max_gain(max_gain),
                        doc='''The maximum gain to apply to the sound.
         
         The gain is clamped after distance attenuation.  The default value
         is 1.0.
         
         :type: float
         ''')

    def _set_position(self, position):
        raise NotImplementedError('abstract')

    position = property(lambda self: self._position,
                        lambda self, position: self._set_position(position),
                        doc='''The position of the sound in 3D space.

        The position is given as a tuple of floats (x, y, z).  The unit
        defaults to meters, but can be modified with the listener
        properties.
        
        :type: 3-tuple of float
        ''')

    def _set_velocity(self, velocity):
        raise NotImplementedError('abstract')

    velocity = property(lambda self: self._velocity,
                        lambda self, velocity: self._set_velocity(velocity),
                        doc='''The velocity of the sound in 3D space.

        The velocity is given as a tuple of floats (x, y, z).  The unit
        defaults to meters per second, but can be modified with the listener
        properties.
        
        :type: 3-tuple of float
        ''')

    def _set_pitch(self, pitch):
        raise NotImplementedError('abstract')

    pitch = property(lambda self: self._pitch,
                     lambda self, pitch: self._set_pitch(pitch),
                     doc='''The pitch shift to apply to the sound.

        The nominal pitch is 1.0.  A pitch of 2.0 will sound one octave
        higher, and play twice as fast.  A pitch of 0.5 will sound one octave
        lower, and play twice as slow.  A pitch of 0.0 is not permitted.
        
        The pitch shift is applied to the source before doppler effects.
        
        :type: float
        ''')

    def _set_cone_orientation(self, cone_orientation):
        raise NotImplementedError('abstract')

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
        raise NotImplementedError('abstract')

    cone_inner_angle = property(lambda self: self._cone_inner_angle,
                                lambda self, a: self._set_cone_inner_angle(a),
                                doc='''The interior angle of the inner cone.
                                
        The angle is given in degrees, and defaults to 360.  When the listener
        is positioned within the volume defined by the inner cone, the sound
        is played at normal gain (see `volume`).
        
        :type: float
        ''')

    def _set_cone_outer_angle(self, cone_outer_angle):
        raise NotImplementedError('abstract')

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
        raise NotImplementedError('abstract')

    cone_outer_gain = property(lambda self: self._cone_outer_gain,
                                lambda self, g: self._set_cone_outer_gain(g),
                                doc='''The gain applied outside the cone.
                                
        When the listener is positioned outside the volume defined by the
        outer cone, this gain is applied instead of `volume`.
        
        :type: float
        ''')

    texture = property(lambda self: self._texture,
                       doc='''The video texture.

        You should rerequest this property every time you display a frame
        of video, as multiple textures might be used.  This property will
        be `None` if there is no video in the current source.

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
BasePlayer.register_event_type('on_eos')

class ManagedSoundPlayerMixIn(object):
    def __init__(self):
        super(ManagedSoundPlayerMixIn, self).__init__()
        managed_players.append(self)

    def stop(self):
        managed_players.remove(self)

class Listener(object):
    '''The listener properties for positional audio.

    You can obtain the singleton instance of this class as
    `pyglet.media.listener`.
    '''
    _volume = 1.0
    _position = (0, 0, 0)
    _velocity = (0, 0, 0)
    _forward_orientation = (0, 0, -1)
    _up_orientation = (0, 1, 0)

    _doppler_factor = 1.
    _speed_of_sound = 343.3

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

    def _set_velocity(self, velocity):
        raise NotImplementedError('abstract')

    velocity = property(lambda self: self._velocity,
                        lambda self, velocity: self._set_velocity(velocity),
                        doc='''The velocity of the listener in 3D space.

        The velocity is given as a tuple of floats (x, y, z).  The unit
        defaults to meters per second, but can be modified with the listener
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

    def _set_doppler_factor(self, factor):
        raise NotImplementedError('abstract')

    doppler_factor = property(lambda self: self._doppler_factor,
                              lambda self, f: self._set_doppler_factor(f),
                              doc='''The emphasis to apply to the doppler
        effect for sounds that move relative to the listener.

        The default value is 1.0, which results in a physically-based
        calculation.  The effect can be enhanced by using a higher factor,
        or subdued using a fractional factor (negative factors are
        ignored).
        
        :type: float
        ''')

    def _set_speed_of_sound(self, speed_of_sound):
        raise NotImplementedError('abstract')

    speed_of_sound = property(lambda self: self._speed_of_sound,
                              lambda self, s: self._set_speed_of_sound(s),
                              doc='''The speed of sound, in units per second.

        The default value is 343.3, a typical result at sea-level on a mild
        day, using meters as the distance unit.

        The speed of sound only affects the calculation of pitch shift to 
        apply due to doppler effects; in particular, no propogation delay
        or relative phase adjustment is applied (in current implementations
        of audio devices).

        :type: float
        ''')

if getattr(sys, 'is_epydoc', False):
    #: The singleton listener.
    #:
    #: :type: `Listener`
    listener = Listener()

    # Document imaginary Player class
    Player = BasePlayer
    Player.__name__ = 'Player'
    del BasePlayer

    # Document imaginary ManagedSoundPlayer class.  (Actually implemented
    # by ManagedSoundPlayerMixIn).
    class ManagedSoundPlayer(Player):
        '''A player which takes care of updating its own audio buffers.

        This player will continue playing the sound until the sound is
        finished, even if the application discards the player early.
        There is no need to call `Player.dispatch_events` on this player,
        though you must call `pyglet.media.dispatch_events`.
        '''

        #: The only possible end of stream action for a managed player.
        EOS_STOP = 'stop'

        eos_action = property(lambda self: EOS_STOP,
                              doc='''The fixed eos_action is `EOS_STOP`,
            in which the player is discarded as soon as the source has
            finished.

            Read-only.
            
            :type: str
            ''')

else:
    # Find best available sound driver according to user preference
    import pyglet
    driver = None
    for driver_name in pyglet.options['audio_driver']:
        try:
            driver_name = 'pyglet.media.drivers.' + driver_name
            __import__(driver_name)
            driver = sys.modules[driver_name]
            break
        except ImportError:
            pass

    if not driver:
        raise ImportError('No suitable audio driver could be loaded.')

    driver.driver_init()
    Player = driver.DriverPlayer
    ManagedSoundPlayer = driver.DriverManagedSoundPlayer
    listener = driver.driver_listener

    # Find best available source loader
    try:
        from pyglet.media import avbin
        _source_class = avbin.AVbinSource
    except ImportError:
        from pyglet.media import riff
        _source_class = riff.WaveSource

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
    '''
    for player in managed_players:
        player.dispatch_events()
