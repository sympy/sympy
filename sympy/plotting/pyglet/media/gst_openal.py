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
# $Id: gst_openal.py 1033 2007-07-13 03:38:16Z Alex.Holkner $

'''Provide a media device that decodes with Gstreamer and plays back with
OpenAL.
'''

# There are a LOT of threads going on here, all of them managed by Gstreamer.
# If pyglet ever needs to run under a Python that doesn't have a GIL, some
# locks will need to be introduced to prevent concurrency catastrophes.  
#
# At the moment, no locks are used because we assume only one thread is
# executing Python code at a time.  Some semaphores are used to block and wake
# up the main thread when needed, these are all instances of
# threading.Semaphore.  Note that these don't represent any kind of
# thread-safety.

import ctypes
import time
import threading

import pyglet.lib
from pyglet.gl import *
from pyglet import image
from pyglet.media import Medium, MediaException, InvalidMediumException, Video
from pyglet.media import gstreamer
from pyglet.media import lib_openal as al
from pyglet.media import openal
from pyglet.media.gstreamer import gst, gstbase, gobject

# Minimum time (secs) to allow for decoding another frame or audio buffer, in
# case sleep time is exactly equal to buffered time.  Also compensates for
# late return from sleep.
SLEEP_UNDERSHOOT = 0.1

DecodebinNewDecodedPad = ctypes.CFUNCTYPE(None,
    ctypes.POINTER(gstreamer.GstElement), ctypes.POINTER(gstreamer.GstPad),
    ctypes.c_int, ctypes.c_void_p)
DecodebinNoMorePads = ctypes.CFUNCTYPE(None,
    ctypes.POINTER(gstreamer.GstElement), ctypes.c_void_p)

GstBusFunc = ctypes.CFUNCTYPE(ctypes.c_int,
    ctypes.c_void_p, ctypes.POINTER(gstreamer.GstMessage), ctypes.c_void_p)

class PythonFileObjectSrcPad(gstreamer.Pad):
    name = 'src'
    direction = gstreamer.GST_PAD_SRC
    caps = 'ANY'

class PythonFileObjectSrc(gstreamer.Element):
    name = 'pyfileobjectsrc'
    klass = 'Source/PythonFileObject'
    description = 'Read from a Python file object'
    author = 'pyglet.org'

    class_struct = gstreamer.GstBaseSrcClass
    instance_struct = gstreamer.GstBaseSrc
    parent_type = gstbase.gst_base_src_get_type

    vmethods = ('create', 'start', 'stop', 'is_seekable', 'do_seek')

    pad_templates = (PythonFileObjectSrcPad,)
    pad_instances = [
        # Instance is created from template by superclass (GstBaseSrcClass)
    ]

    def init(self, file):
        # Not a vmethod; this is the Python interface
        self.file = file
        self.offset =  0

    def start(self):
        return True

    def stop(self):
        return True

    def is_seekable(self):
        return False
        # TODO return hasattr(self.file, 'seek')

    def do_seek(self, segment):
        # TODO
        return True

    def create(self, offset, size, bufp):
        if offset != self.offset:
            self.file.seek(offset)
            self.offset = offset

        try:
            data = self.file.read(size)
            if not data and size > 0:
                # EOF
                return gstreamer.GST_FLOW_UNEXPECTED
        except IOError:
            return gstreamer.GST_FLOW_ERROR

        gst.gst_buffer_new_and_alloc.restype = \
            ctypes.POINTER(gstreamer.GstBuffer)
        buffer = gst.gst_buffer_new_and_alloc(size)

        size = len(data)
        ctypes.memmove(buffer.contents.data, data, size)

        buffer.contents.size = size
        buffer.contents.offset = offset
        buffer.contents.offset_end = offset + size
        buffer.contents.timestamp = gstreamer.GST_CLOCK_TIME_NONE

        self.offset += size

        bufp.contents.contents = buffer.contents

        return gstreamer.GST_FLOW_OK


class GstreamerDecoder(object):
    _decodebin_pads_semaphore = None
    _have_pads = False

    def __init__(self):
        # Lock a semaphore to block until all pads are decoded.
        self._decodebin_pads_semaphore = threading.Semaphore(0)

    def _create_decoder_pipeline(self, filename, file):
        if file is None:
            # Hooray, can just use filename
            src = gst.gst_element_factory_make('filesrc', 'src')
            gobject.g_object_set(src, 'location', filename, None)
        else:
            # Boo, use PythonFileObjectSrc
            src = gst.gst_element_factory_make('pyfileobjectsrc', 'src')
            pysrc = PythonFileObjectSrc.get_instance(src)
            pysrc.init(file)
            
            # Don't GC the source
            self._file_source = pysrc

        # Create pipeline to manage pause/play state.
        pipeline = gst.gst_pipeline_new('pipeline')

        # Create decodebin to find audio and video streams in file
        decoder = gst.gst_element_factory_make('decodebin', 'decoder')
        self._new_decoded_pad_func = DecodebinNewDecodedPad(
            self._new_decoded_pad)
        self._no_more_pads_func = DecodebinNoMorePads(self._no_more_pads)
        gobject.g_signal_connect_data(decoder, 'new-decoded-pad', 
            self._new_decoded_pad_func, None, None, 0)
        gobject.g_signal_connect_data(decoder, 'no-more-pads',
            self._no_more_pads_func, None, None, 0)
        gst.gst_bin_add_many(pipeline, src, decoder, None)
        gst.gst_element_link(src, decoder)

        # Create audioconvert
        self.convert = gst.gst_element_factory_make('audioconvert', 'convert')
        gst.gst_bin_add(pipeline, self.convert)

        return pipeline

    def _destroy_pipeline(self, pipeline):
        gst.gst_element_set_state(pipeline, gstreamer.GST_STATE_NULL)
        gst.gst_object_unref(pipeline)
        del self._no_more_pads_func
        del self._new_decoded_pad_func

    def _new_decoded_pad(self, decodebin, pad, last, data):
        '''Called by decodebin element when a source pad is created.'''
        caps = gst.gst_pad_get_caps(pad)
        struct = gst.gst_caps_get_structure(caps, 0)
        gst.gst_structure_get_name.restype = ctypes.c_char_p
        name = gst.gst_structure_get_name(struct)

        if name.startswith('audio/x-raw'):
            channels = ctypes.c_int()
            gst.gst_structure_get_int(
                struct, 'channels', ctypes.byref(channels))
            depth = sample_rate = 0 # TODO
            self._new_audio_pad(pad, channels.value, depth, sample_rate)
        elif name.startswith('video/x-raw'):
            self._new_video_pad(pad)

    def _new_audio_pad(self, pad, channels, depth, sample_rate):
        pass

    def _new_video_pad(self, pad):
        pass

    def _no_more_pads(self, decodebin, data):
        self._decodebin_pads_semaphore.release()

    def _no_more_pads_deferred(self):
        '''Callback from main application thread when all pads have been
        decoded.  In contrast, _no_more_pads is called from the streaming
        thread.'''
        pass

    def _wait_no_more_pads(self):
        self._decodebin_pads_semaphore.acquire()
        self._decodebin_pads_semaphore.release()
        if not self._have_pads:
            self._have_pads = True
            self._no_more_pads_deferred()


class OpenALSinkPad(gstreamer.Pad):
    name = 'sink'
    direction = gstreamer.GST_PAD_SINK
    caps = '''audio/x-raw-int,
                width = (int) 16,
                depth = (int) 16,
                signed = (boolean) TRUE,
                endianness = (int) BYTE_ORDER,
                channels = (int) { 1, 2 },
                rate = (int) [ 1, MAX ];
              audio/x-raw-int,
                width = (int) 8,
                depth = (int) 8,
                signed = (boolean) FALSE,
                channels = (int) { 1, 2 },
                rate = (int) [ 1, MAX ]
        '''

    def setcaps(self, this, caps):
        gst.gst_caps_get_structure.restype = ctypes.c_void_p
        structure = gst.gst_caps_get_structure(caps, 0)

        rate = ctypes.c_int()
        gst.gst_structure_get_int(structure, 'rate', ctypes.byref(rate))
        self.rate = rate.value

        channels = ctypes.c_int()
        gst.gst_structure_get_int(structure, 'channels', ctypes.byref(channels))

        depth = ctypes.c_int()
        gst.gst_structure_get_int(structure, 'depth', ctypes.byref(depth))

        self.format = openal.get_format(channels.value, depth.value)
        self.element.channels = channels.value
        self.element.rate = rate.value

        self.bytes_per_second = \
            float(channels.value * depth.value / 8 * rate.value)

        return True

    def chain(self, this, buffer):
        timestamp = buffer.contents.timestamp * 0.000000001
        albuffer = openal.buffer_pool.get(timestamp)
        size = buffer.contents.size
        al.alBufferData(albuffer, self.format, 
                        buffer.contents.data, size,
                        self.rate)
        gst.gst_mini_object_unref(buffer)

        self.element._add_buffer(albuffer, size / self.bytes_per_second)
        return gstreamer.GST_FLOW_OK

    def event(self, this, event):
        if event.contents.type == gstreamer.GST_EVENT_EOS:
            self.element._finished_buffering()
            return True
        return False

class TextureSinkPad(gstreamer.Pad):
    name = 'texturesink'
    direction = gstreamer.GST_PAD_SINK
    caps = '''video/x-raw-rgb,
                bpp = (int) 32,
                depth = (int) 24,
                endianness = (int) BIG_ENDIAN,
                red_mask = (int) 0xFF000000,
                green_mask = (int) 0x00FF0000,
                blue_mask = (int) 0x0000FF00
        '''

    def setcaps(self, this, caps):
        gst.gst_caps_get_structure.restype = ctypes.c_void_p
        structure = gst.gst_caps_get_structure(caps, 0)

        width = ctypes.c_int()
        gst.gst_structure_get_int(structure, 'width', ctypes.byref(width))

        height = ctypes.c_int()
        gst.gst_structure_get_int(structure, 'height', ctypes.byref(height))

        self.element._set_frame_format(width.value, height.value)

        return True

    def chain(self, this, buffer):
        self.element._add_frame(buffer.contents.data, 
                                buffer.contents.size,
                                buffer.contents.timestamp * 0.000000001)
        gst.gst_mini_object_unref(buffer)
        return gstreamer.GST_FLOW_OK

    def event(self, this, event):
        return False

# Gstreamer thread management
# -----------------------------------------------------------------------------

class GstreamerMedium(Medium, GstreamerDecoder):
    # The have_video and have_sound properties have to block until all pads
    # have been discovered.
    _has_audio = False
    _has_video = False

    def _get_has_audio(self):
        self._wait_no_more_pads()
        return self._has_audio
    has_audio = property(_get_has_audio)

    def _get_has_video(self):
        self._wait_no_more_pads()
        return self._has_video
    has_video = property(_get_has_video)

# OpenAL streaming
# -----------------------------------------------------------------------------

class OpenALStreamingSinkElement(gstreamer.Element):
    name = 'openalstreamingsink'
    klass = 'audio/openal-streaming-sink'
    description = 'Sink to streaming OpenAL buffers'
    author = 'pyglet.org'

    pad_templates = (OpenALSinkPad,TextureSinkPad)
    pad_instances = [
        ('sink', OpenALSinkPad),
        ('texturesink', TextureSinkPad),
    ]

    def init(self, sound, video=None):
        self.sound = sound
        self.video = video

    def _add_buffer(self, buffer, buffer_time):
        al.alSourceQueueBuffers(self.sound.source, 1, buffer)
        if buffer_time - self.sound.time > self.sound._buffer_time:
            time.sleep(self.sound._buffer_time - SLEEP_UNDERSHOOT)

    def _set_frame_format(self, width, height):
        self.video._set_frame_format(width, height)

    def _add_frame(self, data, length, timestamp):
        self.video._add_frame(data, length, timestamp)

    def _finished_buffering(self):
        self.sound._on_eos()

class GstreamerOpenALStreamingSound(openal.OpenALStreamingSound, 
                                    GstreamerDecoder):
    _buffer_time = .5  # seconds ahead to buffer

    def __init__(self, filename, file):
        super(GstreamerOpenALStreamingSound, self).__init__()
        self._pipeline = self._create_decoder_pipeline(filename, file)

        self.sink = gst.gst_element_factory_make('openalstreamingsink', 'sink')
        gst.gst_bin_add(self._pipeline, self.sink)
        gst.gst_element_link(self.convert, self.sink)

        gst.gst_element_set_state(self._pipeline, gstreamer.GST_STATE_PAUSED)

    def play(self):
        instances.append(self)
        super(GstreamerOpenALStreamingSound, self).play()
        gst.gst_element_set_state(self._pipeline, gstreamer.GST_STATE_PLAYING)

    def pause(self):
        super(GstreamerOpenALStreamingSound, self).pause()
        gst.gst_element_set_state(self._pipeline, gstreamer.GST_STATE_PAUSED)

    def stop(self):
        if self in instances:
            instances.remove(self)
        super(GstreamerOpenALStreamingSound, self).stop()
        self._destroy_pipeline(self._pipeline)
        self._pipeline = None

    def _new_audio_pad(self, pad, channels, depth, sample_rate):
        '''Create and connect the sink for the given source pad.'''
        convertpad = gst.gst_element_get_pad(self.convert, 'sink')
        gst.gst_pad_link(pad, convertpad)
        gst.gst_object_unref(convertpad)

        sink = self.sink
        pysink = OpenALStreamingSinkElement.get_instance(self.sink)
        pysink.init(self)

    def _on_eos(self):
        pass

class GstreamerOpenALStreamingVideo(Video,
                                    GstreamerDecoder):
    _buffer_time = .5  # seconds ahead to buffer

    _sound = None
    _texture = None
    
    _width = None
    _height = None
    _frame_format_semaphore = None

    # For videos with no sound, we use the system clock as a timer (and no
    # pitch shifting is possible).  _start_time gives the system time the
    # video was started (and is updated to reflect pause/resume).  When
    # paused, _start_time has the video position time instead.
    _start_time = 0

    def __init__(self, filename, file):
        super(GstreamerOpenALStreamingVideo, self).__init__()

        # Block frame format until width/height arrives
        self._frame_format_semaphore = threading.Semaphore(0)

        self._pipeline = self._create_decoder_pipeline(filename, file)

        self.sink = gst.gst_element_factory_make('openalstreamingsink', 'sink')
        gst.gst_bin_add(self._pipeline, self.sink)
        gst.gst_element_link(self.convert, self.sink)

        self.videoconvert = \
            gst.gst_element_factory_make('ffmpegcolorspace', 'videoconvert')
        gst.gst_bin_add(self._pipeline, self.videoconvert)
        gst.gst_element_link(self.videoconvert, self.sink)

        gst.gst_element_set_state(self._pipeline, gstreamer.GST_STATE_PAUSED)

        self._frames = [] # queued upcoming video frames
        self._last_frame_timestamp = 0.

    def _get_sound(self):
        self._wait_no_more_pads()
        return self._sound
    sound = property(_get_sound)

    def _wait_frame_format(self):
        self._frame_format_semaphore.acquire()
        self._frame_format_semaphore.release()

    def _get_width(self):
        self._wait_frame_format()
        return self._width
    width = property(_get_width)

    def _get_height(self):
        self._wait_frame_format()
        return self._height
    height = property(_get_height)

    def _get_texture(self):
        self._wait_frame_format()
        if not self._texture:
            glEnable(GL_TEXTURE_2D)
            texture = image.Texture.create_for_size(GL_TEXTURE_2D,
                self.width, self.height, GL_RGB)
            if texture.width != self.width or texture.height != self.height:
                self._texture = texture.get_region(
                    0, 0, self.width, self.height)
            else:
                self._texture = texture

            # Flip texture coords (good enough for simple apps).
            bl, br, tr, tl = self._texture.tex_coords
            self._texture.tex_coords = tl, tr, br, bl
        return self._texture
    texture = property(_get_texture)

    def play(self):
        instances.append(self)
        self._wait_no_more_pads()
        gst.gst_element_set_state(self._pipeline, gstreamer.GST_STATE_PLAYING)
        if self.sound:
            self.sound.play()
        else:
            self._start_time = time.time() - self._start_time
        self.playing = True

    def pause(self):
        self._wait_no_more_pads()
        gst.gst_element_set_state(self._pipeline, gstreamer.GST_STATE_PAUSED)
        if self.sound:
            self.sound.pause()
        else:
            self._start_time = time.time() - self._start_time
        self.playing = False

    def stop(self):
        instances.remove(self)
        if self.sound:
            self.sound.stop()
        self._destroy_pipeline(self._pipeline)
        self._pipeline = None
        self.playing = False
        self.finished = True

    def _get_time(self):
        self._wait_no_more_pads()
        if self.sound:
            return self.sound.time
        else:
            if self.playing:
                return time.time() - self._start_time
            else:
                return self._start_time

    def _new_audio_pad(self, pad, channels, depth, sample_rate):
        '''Create and connect the sink for the given source pad.'''
        convertpad = gst.gst_element_get_pad(self.convert, 'sink')
        gst.gst_pad_link(pad, convertpad)
        gst.gst_object_unref(convertpad)

        self._sound = openal.OpenALStreamingSound()
        self._sound._buffer_time = self._buffer_time

    def _new_video_pad(self, pad):
        sinkpad = gst.gst_element_get_pad(self.videoconvert, 'sink')
        gst.gst_pad_link(pad, sinkpad)
        gst.gst_object_unref(sinkpad)

    def _set_frame_format(self, width, height):
        self._width = width
        self._height = height
        self._frame_format_semaphore.release()
    
    def _add_frame(self, buffer, length, timestamp):
        # Make a copy of the frame and queue it up with its timestamp
        frame_buffer = (ctypes.c_byte * length)()
        ctypes.memmove(frame_buffer, buffer, length)
        self._frames.append((frame_buffer, timestamp))

        # If we've queued enough, sleep (otherwise chain will just run hot
        # and exhaust cpu and memory).
        if timestamp > self._frames[0][1] + self._buffer_time:
            time.sleep(self._buffer_time - SLEEP_UNDERSHOOT)

    def _no_more_pads(self, decodebin, data):
        super(GstreamerOpenALStreamingVideo, self)._no_more_pads(
            decodebin, data)
        sink = self.sink
        pysink = OpenALStreamingSinkElement.get_instance(self.sink)
        pysink.init(self.sound, self)

    def _on_eos(self):
        pass

    def dispatch_events(self):
        if self.sound:
            self.sound.dispatch_events()

        time = self.time
        frame = None
        timestamp = self._last_frame_timestamp
        while self._frames and time >= timestamp:
            frame, timestamp = self._frames.pop(0)
        self._last_frame_timestamp = timestamp

        if frame:
            texture = self.texture
            glBindTexture(texture.target, texture.id)
            glTexSubImage2D(texture.target,
                            texture.level,
                            0, 0,
                            self.width, self.height,
                            GL_RGBA, # TODO
                            GL_UNSIGNED_BYTE,
                            frame)

class GstreamerOpenALStreamingMedium(GstreamerMedium):
    def __init__(self, filename, file):
        super(GstreamerOpenALStreamingMedium, self).__init__()
        self.filename = filename
        self.file = file

        # Set up a decoder pipeline just enough to determine if there's
        # sound/video in the medium.  The pipeline is otherwise unused.
        self._pipeline = self._create_decoder_pipeline(filename, file)
        gst.gst_element_set_state(self._pipeline, gstreamer.GST_STATE_PLAYING)

    def _new_audio_pad(self, pad, channels, depth, sample_rate):
        self._has_audio = True

    def _new_video_pad(self, pad):
        self._has_video = True

    def _no_more_pads(self, decodebin, data):
        super(GstreamerOpenALStreamingMedium, self)._no_more_pads(
            decodebin, data)

    def _no_more_pads_deferred(self):
        gst.gst_element_set_state(self._pipeline, gstreamer.GST_STATE_NULL)
        gst.gst_object_unref(self._pipeline)

    def get_sound(self):
        if not self.has_audio:
            raise InvalidMediumException('No audio in medium')
        return GstreamerOpenALStreamingSound(self.filename, self.file)

    def get_video(self):
        if not self.has_video:
            raise InvalidMediumException('No video in medium')
        return GstreamerOpenALStreamingVideo(self.filename, self.file)

# OpenAL static
# -----------------------------------------------------------------------------

class OpenALStaticSinkElement(gstreamer.Element):
    name = 'openalstaticsink'
    klass = 'audio/openal-static-sink'
    description = 'Sink to static OpenAL buffers'
    author = 'pyglet.org'

    pad_templates = (OpenALSinkPad,)
    pad_instances = [
        ('sink', OpenALSinkPad),
    ]

    def init(self):
        # List of buffers, locked until completed decoding.
        self.buffers = []
        self.buffers_semaphore = threading.Semaphore(0)

    def _add_buffer(self, buffer, buffer_time):
        self.buffers.append(buffer)

    def _finished_buffering(self):
        # Convert list of buffers to an array of buffers for fast queueing
        self.buffers = (al.ALuint * len(self.buffers))(*self.buffers)
        self.buffers_semaphore.release()

class GstreamerOpenALSound(openal.OpenALSound):
    '''Subclass to allow scheduling and unscheduling of event dispatch.
    '''
    def play(self):
        instances.append(self)
        super(GstreamerOpenALSound, self).play()

    def stop(self):
        instances.remove(self)
        super(GstreamerOpenALSound, self).stop()
       
class GstreamerOpenALStaticMedium(GstreamerMedium):
    '''A Medium which decodes all data into a list of OpenAL buffers which
    are shared by many OpenALSound instances.

    The pipeline created during decoding is::

        filesrc ! decodebin ! audioconvert ! openalstaticsink

    OpenALStaticSinkElement implements ``openalstaticsink`` and contains
    a reference to this instance.
    '''
    _buffers = None

    # XXX The __del__ method is not called if get_sound was never called, as
    # there is a circular reference with bound callback methods from
    # GstreamerDecoder.  This can cause streaming thread to still be active
    # while Python is cleaning up.  Presents as:
    #
    #   Error in sys.excepthook:
    #   Original exception was:
    #
    # which is ugly but harmless.
    #
    # If it _were_ called, we would call self._destroy_pipeline, to
    # prevent this problem.  Since the pipeline is already destroyed if
    # get_sound is called, there is no point in a destructor with this flaw.

    def __init__(self, filename, file=None):
        super(GstreamerOpenALStaticMedium, self).__init__()

        self._pipeline = self._create_decoder_pipeline(filename, file)
        self._sink = gst.gst_element_factory_make('openalstaticsink', 'sink')
        gst.gst_bin_add(self._pipeline, self._sink)
        gst.gst_element_link(self.convert, self._sink)

        self._element = OpenALStaticSinkElement.get_instance(self._sink)
        self._element.init()

        gst.gst_element_set_state(self._pipeline, gstreamer.GST_STATE_PLAYING)

        
    def _new_audio_pad(self, pad, channels, depth, sample_rate):
        # Connect the decoded audio pad to the audioconvert
        convertpad = gst.gst_element_get_pad(self.convert, 'sink')
        gst.gst_pad_link(pad, convertpad)
        gst.gst_object_unref(convertpad)

        self._has_audio = True

    def get_sound(self):
        if not self.has_audio:
            raise InvalidMediumException('No audio in medium')

        if self._buffers is None:
            # Wait until sound has finished decoding
            self._element.buffers_semaphore.acquire()
            self._element.buffers_semaphore.release()

            # Save buffers here, then throw away the decoding pipeline
            self._buffers = self._element.buffers
            self._destroy_pipeline(self._pipeline)
            self._pipeline = None

        sound = GstreamerOpenALSound()
        al.alSourceQueueBuffers(sound.source, len(self._buffers), self._buffers)
        return sound

 
class OpenALPlugin(gstreamer.Plugin):
    '''A static plugin to hold the ``openalstreamingsink`` and
    ``openalstaticsink`` elements.
    
    This is equivalent to the GST_PLUGIN_DEFINE macro.
    '''

    name = 'openal-plugin'
    description = 'OpenAL plugin'
    version = '1.0'
    license = 'LGPL'
    package = 'pyglet'
    origin = 'http://www.pyglet.org'

    elements = (
        PythonFileObjectSrc,
        OpenALStreamingSinkElement, 
        OpenALStaticSinkElement)


# Device interface
# --------------------------------------------------------------------------

def load(filename, file=None, streaming=None):
    if streaming is None:
        # Gstreamer can't tell us the duration of a file, so don't stream
        # unless asked to.  This gives best runtime performance at the
        # expense of higher CPU and memory usage for long music clips.
        streaming = False
    if streaming:
        return GstreamerOpenALStreamingMedium(filename, file)
    else:
        return GstreamerOpenALStaticMedium(filename, file)

def dispatch_events():
    global instances

    gstreamer.heartbeat()
    for instance in instances:
        instance.dispatch_events()
    instances = [instance for instance in instances if not instance.finished]

def init():
    gthread = pyglet.lib.load_library('gthread-2.0')
    thread_supported = ctypes.cast(gthread.g_threads_got_initialized,
                                   ctypes.POINTER(ctypes.c_int)).contents.value
    if not thread_supported:
        gthread.g_thread_init(None)

    openal.init()
    gstreamer.init()
    openal_plugin = OpenALPlugin()
    openal_plugin.register()

def cleanup():
    while instances:
        instance = instances.pop()
        instance.stop()
    import gc
    gc.collect()

listener = openal.OpenALListener()

# Active instances
instances = []


