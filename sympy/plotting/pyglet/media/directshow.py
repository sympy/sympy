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
# $Id:$

import atexit
import ctypes
import comtypes
from comtypes import client, GUID

from pyglet.gl import *
import pyglet.window.win32.types as wintypes
from pyglet import image

from pyglet.media import Sound, Video, Medium, Listener, MediaException

_qedit = client.GetModule('qedit.dll')
_quartz = client.GetModule('quartz.dll')
_dsound = client.GetModule('dx8vb.dll')

_qedit_c = ctypes.windll.LoadLibrary('qedit.dll')
_dsound_c = ctypes.windll.LoadLibrary('dsound.dll')

CLSID_FilterGraph =   '{e436ebb3-524f-11ce-9f53-0020af0ba770}'
CLSID_SampleGrabber = '{c1f400a0-3f08-11d3-9f0b-006008039e37}'
CLSID_NullRenderer =  '{c1f400a4-3f08-11d3-9f0b-006008039e37}'
CLSID_DSoundAudioRenderer = '{79376820-07D0-11CF-A24D-0020AFD79767}'

MEDIATYPE_Audio =  GUID('{73647561-0000-0010-8000-00AA00389B71}')
MEDIATYPE_Video =  GUID('{73646976-0000-0010-8000-00AA00389B71}')
MEDIASUBTYPE_PCM = GUID('{00000001-0000-0010-8000-00AA00389B71}')
MEDIASUBTYPE_RGB24 = GUID('{e436eb7d-524f-11ce-9f53-0020af0ba770}')
MEDIASUBTYPE_RGB32 = GUID('{e436eb7e-524f-11ce-9f53-0020af0ba770}')

FORMAT_VideoInfo = GUID('{05589F80-C356-11CE-BF01-00AA0055595A}')
FORMAT_VideoInfo2 = GUID('{F72A76A0-EB0A-11D0-ACE4-0000C0CC16BA}')

INFINITE = 0xFFFFFFFF

AM_MEDIA_TYPE = _qedit._AMMediaType

PINDIR_INPUT = 0
PINDIR_OUTPUT = 1

DBSVOLUME_MIN = -10000
DBSVOLUME_MAX = 0

class VIDEOINFOHEADER(ctypes.Structure):
    _fields_ = [
        ('rcSource', wintypes.RECT),
        ('rcTarget', wintypes.RECT),
        ('dwBitRate', wintypes.DWORD),
        ('dwBitErrorRate', wintypes.DWORD),
        ('AvgTimePerFrame', ctypes.c_int64),
        ('bmiHeader', wintypes.BITMAPINFOHEADER),
    ]

def _connect_filters(graph, source, sink):
    source = source.QueryInterface(_qedit.IBaseFilter)
    sink = sink.QueryInterface(_qedit.IBaseFilter)
    source_pin = _get_filter_pin(source, PINDIR_OUTPUT)
    sink_pin = _get_filter_pin(sink, PINDIR_INPUT)

    graph_builder = graph.QueryInterface(_qedit.IGraphBuilder)
    graph_builder.Connect(source_pin, sink_pin)
    
def _get_filter_pin(filter, direction):
    enum = filter.EnumPins()
    pin, count = enum.Next(1)
    while pin:
        if pin.QueryDirection() == direction:
            return pin
        pin, count = enum.Next(1)
        
class DirectShow_BufferGrabber(comtypes.COMObject):
    _com_interfaces_ = [_qedit.ISampleGrabberCB]

    def __init__(self, buffers):
        self.buffers = buffers

    def ISampleGrabberCB_BufferCB(self, this, sample_time, buffer, len):
        # TODO store sample_time
        dest = (ctypes.c_byte * len)()
        ctypes.memmove(dest, buffer, len)
        self.buffers.append(dest)
        return 0

    def ISampleGrabberCB_SampleCB(self, this, sample_time, sample):
        raise NotImplementedError('Use BufferCB instead')
        return 0

class DirectSound3DBuffer(Sound):
    '''Volume and positioning methods for any sound buffer.'''

    def _set_filter_graph(self, filter_graph):
        '''Locate and set the sound buffer from the filter graph.'''
        enum_filters = filter_graph.EnumFilters()
        filter, _ = enum_filters.Next(1)
        while filter:
            name = filter.QueryFilterInfo().achName
            name = u''.join([unichr(c) for c in name]).strip(u'\0')
            if name == 'Default DirectSound Device':
                self._set_buffer(filter)
                break
            filter, _ = enum_filters.Next(1)            

    def _set_buffer(self, buffer):
        self._basic_audio = buffer.QueryInterface(_quartz.IBasicAudio)
    
    def _set_volume(self, volume):
        self._volume = volume
        # Scale volume to log between -10000 and 0
        volume = int(-10 ** (2 * (1. - volume))  * 100 + 1)
        volume = min(max(-10000, volume), 0)
        self._basic_audio.Volume = volume

class DirectShowStreamingSound(DirectSound3DBuffer):
    def __init__(self, filename):
        filter_graph = client.CreateObject(
            CLSID_FilterGraph, interface=_qedit.IFilterGraph)
        
        filter_builder = filter_graph.QueryInterface(_qedit.IGraphBuilder)
        filter_builder.RenderFile(filename, None)
        self._position = filter_graph.QueryInterface(_quartz.IMediaPosition)
        self._control = filter_graph.QueryInterface(_quartz.IMediaControl)

        self._control.Pause()

        self._stop_time = self._position.StopTime

        self._set_filter_graph(filter_graph)

    def play(self):
        instances.append(self)
        self._control.Run()
        self.playing = True

    def pause(self):
        self._control.Pause()
        self.playing = False

    def stop(self):
        instances.remove(self)
        self._control.Stop()
        self.playing = False

    def dispatch_events(self):
        position = self._position.CurrentPosition
          
        if position >= self._stop_time:
            self.finished = True
            self.playing = False

    def _get_time(self):
        return self._position.CurrentPosition

class DirectShowStreamingVideo(Video):
    _texture = None
    
    def __init__(self, filename):
        # Create filter graph with filename at source
        filter_graph = client.CreateObject(
            CLSID_FilterGraph, interface=_qedit.IFilterGraph)

        graph_builder = filter_graph.QueryInterface(_qedit.IGraphBuilder)
        source_filter = graph_builder.AddSourceFilter(filename, filename)

        # Create a sample grabber for RGB24 video
        self._sample_grabber = client.CreateObject(
            CLSID_SampleGrabber, interface=_qedit.ISampleGrabber)
        filter_graph.AddFilter(self._sample_grabber, None)

        media_type = AM_MEDIA_TYPE()
        media_type.majortype = MEDIATYPE_Video
        media_type.subtype = MEDIASUBTYPE_RGB24
        self._sample_grabber.SetMediaType(media_type)
        self._sample_grabber.SetBufferSamples(True)

        # Samples will be added to self.buffers list by
        # DirectShow_BufferGrabber
        self._buffers = []
        self._buffer_grabber = DirectShow_BufferGrabber(self._buffers)
        self._sample_grabber.SetCallback(self._buffer_grabber, 1)

        _connect_filters(filter_graph, source_filter, self._sample_grabber)
        
        # Connect output of samplegrabber to null to prevent video window
        # showing up.
        null_renderer = client.CreateObject(
            CLSID_NullRenderer, interface=_qedit.IBaseFilter)
        filter_graph.AddFilter(null_renderer, None)
        _connect_filters(filter_graph, self._sample_grabber, null_renderer)
        
        self._control = filter_graph.QueryInterface(_quartz.IMediaControl)
        self._position = filter_graph.QueryInterface(_quartz.IMediaPosition)
        self._stop_time = self._position.StopTime

        # Enumerate the graph and look for an unconnected output pin.
        # Hopefully this will be the audio stream.
        enum_filters = filter_graph.EnumFilters()
        filter, _ = enum_filters.Next(1)
        audio_source_pin = None
        while filter:
            enum_pins = filter.EnumPins()
            pin, _ = enum_pins.Next(1)
            while pin:
                if pin.QueryDirection() == PINDIR_OUTPUT:
                    try:
                        pin.ConnectedTo()
                    except comtypes.COMError:
                        audio_source_pin = pin
                pin, _ = enum_pins.Next(1)
            filter, _ = enum_filters.Next(1)

        # Create a DirectSound Audio Renderer and connect it to what we think
        # is the audio source
        if audio_source_pin:
            audio_renderer = client.CreateObject(
                CLSID_DSoundAudioRenderer, interface=_qedit.IBaseFilter)
            filter_graph.AddFilter(audio_renderer, None)
            audio_renderer_pin = _get_filter_pin(audio_renderer, PINDIR_INPUT)
            graph_builder.Connect(audio_source_pin, audio_renderer_pin)

            # Create public sound interface 
            self.sound = DirectSound3DBuffer()
            self.sound._set_buffer(audio_renderer)

        # Cue for playing, this should fill up buffers enough to determine
        # format.
        self._control.Pause()

        # Find out video format
        self._sample_grabber.GetConnectedMediaType(media_type)
        if media_type.formattype == FORMAT_VideoInfo:
            format = ctypes.cast(media_type.pbFormat,
                                 ctypes.POINTER(VIDEOINFOHEADER)).contents
            self.width = format.bmiHeader.biWidth
            self.height = format.bmiHeader.biHeight
            self._pitch = (self.width * 3 + 3) & ~0x3 # align to 4 bytes
        else:
            raise MediaException('Unsupported video format type')
            
    def play(self):
        instances.append(self)
        self._control.Run()
        self.playing = True

    def pause(self):
        self._control.Pause()
        self.playing = False

    def stop(self):
        instances.remove(self)
        self._control.Stop()
        self.playing = False

        # Wait until it transitions (otherwise we may GC the sample grabber
        # while it's still in use).
        self._control.GetState(INFINITE)
        
        # Clean up in this order to prevent a crash
        del self.sound._basic_audio # yuk, but necessarily like this
        del self._control
        del self._position
        del self._sample_grabber
        del self._buffer_grabber

    def _get_time(self):
        return self._position.CurrentPosition

    def _get_texture(self):
        if not self._texture:
            glEnable(GL_TEXTURE_2D)
            texture = image.Texture.create_for_size(GL_TEXTURE_2D,
                self.width, self.height, GL_RGB)
            if texture.width != self.width or texture.height != self.height:
                self._texture = texture.get_region(
                    0, 0, self.width, self.height)
            else:
                self._texture = texture
        return self._texture
    texture = property(_get_texture) 

    def dispatch_events(self):
        position = self._position.CurrentPosition
          
        if position >= self._stop_time:
            self.finished = True
            self.playing = False
        
        # TODO look at timestamps
        buffer = None
        while self._buffers:
            buffer = self._buffers.pop(0)

        if buffer:
            texture = self.texture
            glBindTexture(texture.target, texture.id)

            # Create an intermediate ImageData and use it to update the
            # texture, as it does a good job of swapping BGR to RGB in the
            # most efficient way.
            imagedata = image.ImageData(
                self.width, self.height, 'BGR', buffer, pitch=self._pitch)
            imagedata.blit_to_texture(texture.target, 0, 0, 0, 0)

    def unschedule(self):
        instances.remove(self)
            

class DirectShowStreamingMedium(Medium):
    def __init__(self, filename, file=None):
        if file is not None:
            raise NotImplementedError('TODO file objects')

        self.filename = filename

        self._check_file()

    def _check_file(self):
        # Open and render a filtergraph for the file to find out if it has
        # audio and/or video, and set self.has_(audio|video).
        
        filter_graph = client.CreateObject(CLSID_FilterGraph, 
                                           interface=_qedit.IFilterGraph)
        
        filter_builder = filter_graph.QueryInterface(_qedit.IGraphBuilder)
        try:
            filter_builder.RenderFile(self.filename, None)

            # Enumerate filters in the graph, see if pins have audio or video
            # on them.
            enum_filters = filter_graph.EnumFilters()
            filter, _ = enum_filters.Next(1)
            while filter:
                enum_pins = filter.EnumPins()
                pin, _ = enum_pins.Next(1)
                while pin:
                    media_type = pin.ConnectionMediaType()
                    self.has_audio |= media_type.majortype == MEDIATYPE_Audio
                    self.has_video |= media_type.majortype == MEDIATYPE_Video
                    pin, _ = enum_pins.Next(1)
                filter, _ = enum_filters.Next(1)
        except comtypes.COMError:
            # Can't render file, then has_audio and has_video are both False
            pass

    def get_sound(self):
        # TODO if has_video, sink it to null so it doesn't pop up in it's own
        # window.
        return DirectShowStreamingSound(self.filename)

    def get_video(self):
        return DirectShowStreamingVideo(self.filename)

class DirectShowStaticMedium(Medium):
    def __init__(self, filename, file=None):
        filter_graph = client.CreateObject(
            CLSID_FilterGraph, interface=_qedit.IFilterGraph)

        graph_builder = filter_graph.QueryInterface(_qedit.IGraphBuilder)
        source_filter = graph_builder.AddSourceFilter(filename, filename)

        sample_grabber = client.CreateObject(
            CLSID_SampleGrabber, interface=_qedit.ISampleGrabber)
        filter_graph.AddFilter(sample_grabber, None)

        media_type = AM_MEDIA_TYPE()
        media_type.majortype = MEDIATYPE_Audio
        media_type.subtype = MEDIASUBTYPE_PCM
        sample_grabber.SetMediaType(media_type)
        sample_grabber.SetBufferSamples(True)

        self.buffers = []
        buffer_grabber = DirectShow_BufferGrabber(self.buffers)
        sample_grabber.SetCallback(buffer_grabber, 1)

        try:
            _connect_filters(filter_graph, source_filter, sample_grabber)
        except comtypes.COMError:
            # Couldn't connect sample grabber... means nothing could provide
            # an audio pin.  Give up now, leave has_audio = False.
            return

        null_renderer = client.CreateObject(
            CLSID_NullRenderer, interface=_qedit.IBaseFilter)
        filter_graph.AddFilter(null_renderer, None)
        _connect_filters(filter_graph, sample_grabber, null_renderer)
        
        media_filter = filter_graph.QueryInterface(_qedit.IMediaFilter)
        media_filter.SetSyncSource(None)

        control = filter_graph.QueryInterface(_quartz.IMediaControl)
        control.Run()

        media_event = filter_graph.QueryInterface(_quartz.IMediaEvent)
        media_event.WaitForCompletion(INFINITE)

        sample_grabber.GetConnectedMediaType(media_type)
        self.format = ctypes.cast(media_type.pbFormat,
                                  ctypes.POINTER(_dsound.WAVEFORMATEX)).contents

        self.sample_rate = self.format.lSamplesPerSec
        self.channels = self.format.nChannels

        self.buffer = []
        for buffer in self.buffers:
            self.buffer += buffer[:]
        self.buffer = (ctypes.c_byte * len(self.buffer))(*self.buffer)

        self.has_audio = True

        # TODO video window will probably pop up if there's a video stream in
        # medium; need to suppress it.

    def get_sound(self):
        return DirectShowStaticSound(self.format, self.buffer)


class DirectShowStaticSound(Sound):
    def __init__(self, format, buffer):
        desc = _dsound.DSBUFFERDESC()
        desc.lFlags = _dsound.DSBCAPS_CTRL3D
        desc.lBufferBytes = len(buffer)
        desc.fxFormat = format
        desc.guid3DAlgorithm = _dsound.GUID_DS3DALG_NO_VIRTUALIZATION
        
        try:
            self.sound_buffer = directsound.CreateSoundBuffer(desc)
            self.sound_buffer.WriteBuffer(0, 0, buffer,
                                          _dsound.DSBLOCK_ENTIREBUFFER)
        except comtypes.COMError:
            raise MediaException('Sound is too long for static buffer')

        self._bytes_per_second = float(format.lAvgBytesPerSec)
        self._last_time = 0

    def play(self):
        instances.append(self)
        self.sound_buffer.Play(_dsound.DSBPLAY_DEFAULT)
        self.playing = True

    def pause(self):
        self.sound_buffer.Stop()
        self.playing = False

    def stop(self):
        instances.remove(self)
        self.sound_buffer.Stop()
        self.playing = False

    def _get_time(self):
        position = self.sound_buffer.GetCurrentPosition()
        return position.lPlay / self._bytes_per_second

    def dispatch_events(self):
        time = self.time
        if time < self._last_time:
            self.finished = True
            self.playing = False
        self._last_time = time

    def unschedule(self):
        instances.remove(self)

# Device interface
# -----------------------------------------------------------------------------

def load(filename, file=None, streaming=None):
    if streaming is None:
        streaming = True

    if streaming:
        return DirectShowStreamingMedium(filename, file)
    else:
        return DirectShowStaticMedium(filename, file)

def dispatch_events():
    global instances
    for instance in instances:
        instance.dispatch_events()
    instances = [instance for instance in instances if not instance.finished]

directx = None
directsound = None

# TODO TODO TODO TODO 
# grab any available window and SetCooperativeLevel when possible.
from pyglet.window import Window
w = Window(visible=False)
# TODO TODO TODO TODO

def init():
    global directx
    global directsound
    directx = client.CreateObject(_dsound.DirectX8._reg_clsid_,
                                  interface=_dsound.IDirectX8)
    directsound = directx.DirectSoundCreate(None)
    directsound.SetCooperativeLevel(w._hwnd, _dsound.DSSCL_PRIORITY)

def cleanup():
    for instance in instances[:]:
        instance.stop()
    
    global directx
    global directsound
    del directsound
    del directx

# XXX temporary, this is critical on Windows to avoid a crash on exit.
# should probably do same/similar on other platforms and remove public
# cleanup() call.
atexit.register(cleanup)

# XXX temporary
listener = Listener()

instances = []

'''
Grabbing real-time samples (e.g., for video):

size = ctypes.c_long(0)
# Need to use raw_func here because func has argtypes mangled too
# much.
GetCurrentBuffer = \
    sample_grabber._ISampleGrabber__com_GetCurrentBuffer
GetCurrentBuffer(size, None)
buffer = (ctypes.c_byte * size.value)()
GetCurrentBuffer(size, 
                 ctypes.cast(buffer, ctypes.POINTER(ctypes.c_long)))
print size, buffer[:]
'''

'''
Valid fields in format

desc.fxFormat.nSize = ctypes.sizeof(desc.fxFormat)
desc.fxFormat.nFormatTag = _dsound.WAVE_FORMAT_PCM
desc.fxFormat.nChannels = 2
desc.fxFormat.lSamplesPerSec = 22050
desc.fxFormat.nBitsPerSample = 16
desc.fxFormat.nBlockAlign = \
    (desc.fxFormat.nChannels * desc.fxFormat.nBitsPerSample) / 8
desc.fxFormat.nAvgBytesPerSec =  \
    (desc.fxFormat.lSamplesPerSec * desc.fxFormat.nBlockAlign)
'''
