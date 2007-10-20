
'''
'''

__docformat__ = 'restructuredtext'
__version__ = '$Id: __init__.py 1223 2007-09-01 11:33:40Z Alex.Holkner $'

import ctypes

from pyglet.media import BasePlayer, ManagedSoundPlayerMixIn, Listener
from pyglet.media import MediaException

from pyglet.media.drivers.alsa import asound

alsa_debug = 'alsa.log'

class ALSAException(MediaException):
    pass

def check(err):
    if err < 0:
        raise ALSAException(asound.snd_strerror(err))
    return err

class Device(object):
    _ring_buffer_time = 0.3
    
    def __init__(self, name):
        self.name = name
        self.pcm = ctypes.POINTER(asound.snd_pcm_t)()
        self.hwparams = ctypes.POINTER(asound.snd_pcm_hw_params_t)()
        self.swparams = ctypes.POINTER(asound.snd_pcm_sw_params_t)() 

        check(asound.snd_pcm_open(ctypes.byref(self.pcm),
                                  name,
                                  asound.SND_PCM_STREAM_PLAYBACK,
                                  asound.SND_PCM_NONBLOCK))
        check(asound.snd_pcm_hw_params_malloc(ctypes.byref(self.hwparams)))
        check(asound.snd_pcm_sw_params_malloc(ctypes.byref(self.swparams)))
        check(asound.snd_pcm_hw_params_any(self.pcm, self.hwparams))

        if alsa_debug:
            asound.snd_output_printf(debug_output, 'New device: %s\n' % name)


    def __del__(self):
        try:
            check(asound.snd_pcm_close(self.pcm))
        except (NameError, AttributeError):
            pass

    def prepare(self, source):
        if not source.audio_format:
            # TODO avoid creating in this case.
            return

        format = {
            8:  asound.SND_PCM_FORMAT_U8,
            16: asound.SND_PCM_FORMAT_S16,
            24: asound.SND_PCM_FORMAT_S24,  # probably won't work
            32: asound.SND_PCM_FORMAT_S32
        }.get(source.audio_format.sample_size)
        if format is None:
            raise ALSAException('Unsupported audio format.')

        check(asound.snd_pcm_hw_params_set_access(self.pcm, self.hwparams,
            asound.SND_PCM_ACCESS_RW_INTERLEAVED))
        check(asound.snd_pcm_hw_params_set_format(self.pcm, self.hwparams,
            format))
        check(asound.snd_pcm_hw_params_set_channels(self.pcm, self.hwparams,
            source.audio_format.channels))
        check(asound.snd_pcm_hw_params_set_rate(self.pcm, self.hwparams,
            source.audio_format.sample_rate, 0))
        check(asound.snd_pcm_hw_params_set_buffer_size(self.pcm, self.hwparams,
            int(self._ring_buffer_time * source.audio_format.sample_rate)))
        check(asound.snd_pcm_hw_params(self.pcm, self.hwparams))
        if alsa_debug:
            check(asound.snd_pcm_dump(self.pcm, debug_output))

        self.can_pause = asound.snd_pcm_hw_params_can_pause(self.hwparams)

class ALSAPlayer(BasePlayer):
    _min_buffer_time = 0.1
    _max_buffer_size = 65536

    def __init__(self):
        super(ALSAPlayer, self).__init__()

        self._sources = []
        self._source_read_index = 0
        self._playing = False
        self._device = None

        # timestamp of last packet buffered (must be of _source_read_index)
        self._current_buffer_time = 0.

        # cumulative duration of all sources completely buffered.
        self._cumulative_buffer_time = 0.

        # asound time for current source begin time
        self._start_time = None

        # saved timestamp when paused
        self._timestamp = 0.

        self._queue_audio_data = None

    def queue(self, source):
        source = source._get_queue_source()

        if not self._sources:
            self._source_read_index = 0
            source._init_texture(self)
        self._sources.append(source)

    def next(self):
        if self._sources:
            old_source = self._sources.pop(0)
            old_source._release_texture(self)
            old_source._stop()
            self._source_read_index -= 1

        if self._sources:
            self._sources[0]._init_texture(self)

    def dispatch_events(self):
        if not self._sources:
            return

        if not self._playing:
            # If paused, just update the video texture.
            if self._texture:
                self._sources[0]._update_texture(self, self.time)
            return

        # Create a device if there isn't one.  TODO only if source and source
        # has audio
        if not self._device:
            self._device = Device('plug:front')
            self._device.prepare(self._sources[0])

        self_time = self.time 

        # Passed EOS?
        source = self._sources[0]
        while source and source.duration < self_time:
            if self._eos_action == self.EOS_NEXT:
                self.next()
            elif self._eos_action == self.EOS_STOP:
                self.stop()
                self._sources = []
                return
            self.dispatch_event('on_eos')

            self_time -= source.duration
            self._cumulative_buffer_time -= source.duration
            assert self._cumulative_buffer_time >= -0.001 # some float err ok
            try:
                source = self._sources[0]
                self._set_start_time(self._sources[0], self_time)
            except IndexError:
                source = None
                self._start_time = None

        # Ensure device buffer is full
        try:
            source = self._sources[self._source_read_index]
        except IndexError:
            source = None
        while (source and 
               self._cumulative_buffer_time + self._current_buffer_time - self_time
                  < self._min_buffer_time):
            if self._queue_audio_data:
                audio_data = self._queue_audio_data
                self._queue_audio_data = None
            else:
                max_bytes = int(
                  self._min_buffer_time * source.audio_format.bytes_per_second)
                max_bytes = min(max_bytes, self._max_buffer_size)
                audio_data = source._get_audio_data(max_bytes)

            if audio_data: 
                samples = \
                    audio_data.length // source.audio_format.bytes_per_sample
                samples_out = asound.snd_pcm_writei(self._device.pcm, 
                    audio_data.data, samples)
                if samples_out < 0:
                    if samples_out == -11: # EAGAIN
                        self._queue_audio_data = audio_data
                    elif samples_out == -32: # EPIPE
                        # xrun recovery
                        check(asound.snd_pcm_prepare(self._device.pcm))
                        self._queue_audio_data = audio_data
                    else:
                        raise ALSAException(asound.snd_strerror(samples_out))
                elif samples_out < samples:
                    audio_data.consume(
                        samples_out * source.audio_format.bytes_per_sample,
                        source.audio_format)
                    self._current_buffer_time = audio_data.timestamp
                    self._queue_audio_data = audio_data
                else:
                    self._current_buffer_time = \
                        audio_data.timestamp + audio_data.duration
                    
                if self._start_time is None:
                    # XXX start playback
                    self._set_start_time(source, audio_data.timestamp)
                    
            else:
                # EOS on read source
                self._cumulative_buffer_time += source.duration
                self._current_buffer_time = 0.
                if self._eos_action == self.EOS_NEXT:
                    self._source_read_index += 1
                    try:
                        # preroll
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

        # Update video texture
        if self._texture:
            self._sources[0]._update_texture(self, self_time)

    def _get_time(self):
        if self._start_time is None or not self._playing:
            return self._timestamp
        return self._get_asound_time() - self._start_time

    def _get_asound_time(self):
        status = ctypes.POINTER(asound.snd_pcm_status_t)()
        timestamp = asound.snd_timestamp_t()

        check(asound.snd_pcm_status_malloc(ctypes.byref(status)))
        check(asound.snd_pcm_status(self._device.pcm, status))
        asound.snd_pcm_status_get_tstamp(status, ctypes.byref(timestamp))
        asound.snd_pcm_status_free(status)
        return timestamp.tv_sec + timestamp.tv_usec * 0.000001

    def _clear_ring_buffer(self):
        check(asound.snd_pcm_drop(self._device.pcm))
        check(asound.snd_pcm_prepare(self._device.pcm))
        self._current_buffer_time = 0.
        self._cumulative_buffer_time = 0.
        self._source_read_index = 0
        self._queue_audio_data = None
        self._start_time = None

    def _set_start_time(self, source, timestamp):
        delay = asound.snd_pcm_sframes_t()
        check(asound.snd_pcm_delay(self._device.pcm, ctypes.byref(delay)))
        delay = (delay.value / 
                 source.audio_format.bytes_per_sample /
                 float(source.audio_format.bytes_per_second))
        self._start_time = self._get_asound_time() - timestamp - delay

    def play(self):
        if self._playing:
            return

        self._playing = True

        if not self._sources:
            return

        if self._device:
            state = asound.snd_pcm_state(self._device.pcm)
            if self._device.can_pause and state == asound.SND_PCM_STATE_PAUSED:
                check(asound.snd_pcm_pause(self._device.pcm, 0))
                self._set_start_time(self._sources[0], self._timestamp)

    def pause(self):
        if not self._playing:
            return

        self._playing = False

        if not self._sources:
            return

        if self._device:
            self._timestamp = self._get_asound_time() - self._start_time
            if self._device.can_pause:
                check(asound.snd_pcm_pause(self._device.pcm, 1))
            else:
                # Hardware can't pause, so we'll just drop everything that's
                # been buffered.  Improvement could be to rebuffer data that
                # wasn't played.
                self._clear_ring_buffer()

    def seek(self, timestamp):
        if self._sources:
            self._sources[0]._seek(timestamp)
            self._timestamp = timestamp
            self._start_time = None
            
            if self._device:
                self._clear_ring_buffer()

    def _get_source(self):
        if self._sources:
            return self._sources[0]
        return None

    def _stop(self):
        raise RuntimeError('Invalid eos_action for this player.')

    def _set_volume(self, volume):
        self._volume = volume
        # TODO apply to device

    # All other properties are silently ignored.

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

class ALSAManagedSoundPlayer(ALSAPlayer, ManagedSoundPlayerMixIn):
    pass

class ALSAListener(Listener):
    def set_volume(self, volume):
        # TODO master volume
        self._volume = volume

    # All other properties are silently ignored.

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
    global debug_output
    print asound.snd_asoundlib_version()
    debug_output = ctypes.POINTER(asound.snd_output_t)()
    if alsa_debug:
        check(asound.snd_output_stdio_open(ctypes.byref(debug_output),
                                           alsa_debug,
                                           'w'))

driver_listener = ALSAListener()
DriverPlayer = ALSAPlayer
DriverManagedSoundPlayer = ALSAManagedSoundPlayer
