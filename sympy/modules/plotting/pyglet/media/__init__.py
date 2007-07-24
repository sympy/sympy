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
# $Id: __init__.py 912 2007-06-18 01:21:19Z r1chardj0n3s $

'''Media playback interface.

Only basic functionality is described here; for full reference see the
accompanying documentation.

To load some media::

    from pyglet import media
    sound = media.load('sound.mp3')
    audio = medium.get_audio()

    movie = media.load('movie.mp4')
    video = medium.get_video()

The supported media file types include WAV, MP3, and many more,
depending on the operating system.

Both audio and video support the same API implemeted in `MediumInstance`.

To have media actually play, you will need to invoke
``media.dispatch_events()`` in your application's event loop.

'''

__docformat__ = 'restructuredtext'
__version__ = '$Id: __init__.py 912 2007-06-18 01:21:19Z r1chardj0n3s $'

import sys

from pyglet import event

class MediaException(Exception):
    pass

class InvalidMediumException(MediaException):
    pass

class Medium(object):
    '''An audio and/or video medium that can be played.

    A medium cannot itself be played, but it can provide a sound instance
    which represents an instance of a playing sample.  You can retrieve
    an instance of the sound, which will be queued up as soon as possible,
    with `get_sound`.  
    
    For convenience, the `play` method will return a sound and begin playing
    it as soon as possible (for a static Medium, this will be almost
    immediately).

    :Ivariables:
        `has_audio` : bool
            If True, there is an audio track in the medium, and `get_sound`
            can be called to use it.
        `has_video` : bool
            If True, there is a video track in the medium, and `get_video`
            can be called to use it.

    '''

    _duration = None
    
    has_audio = False
    has_video = False

    def _get_duration(self):
        return self._duration

    duration = property(lambda self: self._get_duration(),
                        doc='''The length of the medium, in seconds.

        Not all media can determine their duration; in this case the value
        is None.

        Read-only.

        :type: float
        ''')

    def get_sound(self):
        '''Create and return a new sound that can be played.

        Each call to this method creates a new sound instance; the multiple
        sounds can be played, paused, and otherwise manipulated independently
        and simultaneously.

        :rtype: `Sound`
        '''
        raise NotImplementedError('abstract')

    def get_video(self):
        '''Create and return a new video that can be played.

        Each call to this method creates a new video instance, which can be
        played and paused independently of any other videos.

        This call can fail if the medium has no video data (i.e., it is a
        sound file), or if the medium was loaded statically instead of
        streaming.

        :rtype: `Video`
        '''

    def play(self):
        '''Play the sound.

        This is a convenience method which creates a sound and plays it
        immediately.

        :rtype: `Sound`
        '''
        sound = self.get_sound()
        sound.play()
        return sound

class MediumInstance(event.EventDispatcher):
    '''An instance of a sound or video.

    :Ivariables:
        `playing` : bool
            If True, the sound is currently playing.  Even after calling
            `play`, the sound may not begin playing until enough audio has
            been buffered.  If False, the sound is either buffering and
            about to play, is explicitly paused, or has finished.
            This variable is read-only.
        `finished` : bool
            If True, the sound has finished playing.  This variable is
            read-only.
    '''

    playing = False
    finished = False

    def play(self):
        '''Begin playing the instance.

        This has no effect if the instance is already playing.
        '''
        raise NotImplementedError('abstract')

    def pause(self):
        '''Pause playback of the instance.

        This has no effect if the instance is already paused.
        '''
        raise NotImplementedError('abstract')

    def stop(self):
        '''Stop playback of the instance and release all resources.

        Once an instance has been stopped, it cannot be started again.
        '''
        raise NotImplementedError('abstract')

    def seek(self, timestamp):
        '''Seek for playback to the indicated timestamp in seconds.
        '''
        raise NotImplementedError('abstract')

    def _get_time(self):
        raise NotImplementedError('abstract')

    time = property(lambda self: self._get_time(),
                    doc='''Retrieve the current playback time of the instance.
                    
         The playback time is a float expressed in seconds, with 0.0 being
         the beginning of the sound.  The playback time returned represents
         the time encoded in the media, and may not reflect actual time
         passed due to pitch shifting or pausing.

         Read-only.

         :type: float
         ''')
 
    def dispatch_events(self):
        '''Dispatch any pending events and perform regular heartbeat functions
        to maintain playback.

        This method is called automatically by `pyglet.media.dispatch_events`,
        there is no need to call this from an application.
        '''
        pass

    def unschedule(self):
        '''Stop event processing for this instance.
        
        This will prevent any further calls to `dispatch_events`.
        '''
        pass


EVENT_FINISHED = MediumInstance.register_event_type('on_finished')

class Sound(MediumInstance):
    '''An instance of a sound, either currently playing or ready to be played.

    By default, monaural sounds are played at nominal volume equally among
    all speakers, however they may also be positioned in 3D space.  Stereo
    sounds are not positionable.

    :Ivariables:
        `depth` : int
            The number of bits per sample per channel (usually 8 or 16).
            The value is None if the audio properties have not yet been
            determined.
        `channels` : int
            The number of audio channels provided: 1 for monoaural sound, 2
            for stereo, or more for multi-channel sound. The value is None if
            the audio properties have not yet been determined.
        `sample_rate` : float
            The audio sample rate, in Hz.  The sound may be resampled
            to match the audio device's sample rate; this value gives
            the original sample rate.  The value is None if the audio
            properties have not yet been determined.

    '''
        
    depth = None
    sample_rate = None
    channels = None

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


class Video(MediumInstance):
    '''A video that can be played.

    :Ivariables:
        `sound` : `Sound`
            Reference to the sound instance that accompanies this video.
        `texture` : `pyglet.image.Texture`
            Reference to the texture object that holds the current frame of
            video.
        `width` : int
            Width of the video, in pixels.  None if unknown.
        `height` : int
            Height of the video, in pixels.  None if unknown.

    '''
    sound = None
    texture = None

    width = None
    height = None

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

    def load(filename, file=None, streaming=None):
        '''Load a medium.

        :Parameters:
            `filename` : str
                Filename to load.
            `file` : file-like object
                File to load data from.  If unspecified, the filename will be
                opened.
            `streaming` : bool
                If True, the medium will be decoded as it is played; otherwise
                it will be decoded immediately and stored in memory.

        :rtype: `Medium`
        '''

    def dispatch_events():
        '''Process audio events.

        You must call this function regularly (typically once per run loop
        iteration) in order to keep audio buffers full and video textures
        up-to-date.
        '''

    def cleanup():
        '''Release all media resources.

        You should call this function before your application exits, if it has
        imported this module.
        '''
else:
    if sys.platform == 'linux2':
        from pyglet.media import gst_openal
        _device = gst_openal
    elif sys.platform == 'darwin':
        from pyglet.media import quicktime
        _device = quicktime
    elif sys.platform in ('win32', 'cygwin'):
        from pyglet.media import directshow
        _device = directshow
    else:
        raise ImportError('pyglet.media not yet supported on %s' % sys.platform)

    load = _device.load
    dispatch_events = _device.dispatch_events
    cleanup = _device.cleanup
    listener = _device.listener
    _device.init()
