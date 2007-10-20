
'''Fallback driver producing no audio.
'''

__docformat__ = 'restructuredtext'
__version__ = '$Id: silent.py 1213 2007-08-30 10:47:05Z Alex.Holkner $'

import time

from pyglet.media import BasePlayer, ManagedSoundPlayerMixIn, Listener
from pyglet.media import MediaException

class SilentPlayer(BasePlayer):
    def __init__(self):
        super(SilentPlayer, self).__init__()

        self._sources = []
        self._playing = False
        self._timestamp = 0.
        self._timestamp_time = None

    def queue(self, source):
        source = source._get_queue_source()

        if not self._sources:
            source._init_texture(self)
        self._sources.append(source)

    def next(self):
        if self._sources:
            old_source = self._sources.pop(0)
            old_source._release_texture(self)
            old_source._stop()

        if self._sources:
            self._sources[0]._init_texture(self)

    def dispatch_events(self):
        if not self._sources:
            return

        if self._playing:
            now = time.time()
            self._timestamp += now - self._timestamp_time
            self._timestamp_time = now

        if self._texture:
            self._sources[0]._update_texture(self, self._timestamp)

        if self._timestamp > self._sources[0].duration:
            self.next()
            self._timestamp = 0.

    def _get_time(self):
        if not self._playing:
            return self._timestamp

        return self._timestamp + time.time() - self._timestamp_time

    def play(self):
        if self._playing:
            return

        self._playing = True

        if not self._sources:
            return

        self._timestamp_time = time.time()

    def pause(self):
        self._playing = False

        if not self._sources:
            return

    def seek(self, timestamp):
        if self._sources:
            self._sources[0]._seek(timestamp)
            self._timestamp = timestamp
            self._timestamp_time = time.time()

    def _get_source(self):
        if self._sources:
            return self._sources[0]
        return None

    def _set_volume(self, volume):
        self._volume = volume

    def _set_min_gain(self, min_gain):
        self._min_gain = min_gain

    def _set_max_gain(self, max_gain):
        self._max_gain = max_gain

    def _set_position(self, position):
        self._position = position

    def _set_velocity(self, velocity):
        self._velocity = velocity

    def _set_pitch(self, pitch):
        self._pitch = pitch

    def _set_cone_orientation(self, cone_orientation):
        self._cone_orientation = cone_orientation

    def _set_cone_inner_angle(self, cone_inner_angle):
        self._cone_inner_angle = cone_inner_angle

    def _set_cone_outer_gain(self, cone_outer_gain):
        self._cone_outer_gain = cone_outer_gain

class SilentManagedSoundPlayer(SilentPlayer, ManagedSoundPlayerMixIn):
    pass

class SilentListener(Listener):
    def set_volume(self, volume):
        self._volume = volume

    def set_position(self, position):
        self._position = position

    def set_velocity(self, velocity):
        self._velocity = velocity

    def set_forward_orientation(self, orientation):
        self._forward_orientation = orientation

    def set_up_orientation(self, orientation):
        self._up_orientation = orientation

    def set_doppler_factor(self, factor):
        self._doppler_factor = factor

    def set_speed_of_sound(self, speed_of_sound):
        self._speed_of_sound = speed_of_sound

def driver_init():
    pass

driver_listener = SilentListener()
DriverPlayer = SilentPlayer
DriverManagedSoundPlayer = SilentManagedSoundPlayer
