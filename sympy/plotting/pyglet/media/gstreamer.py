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
# $Id: gstreamer.py 1033 2007-07-13 03:38:16Z Alex.Holkner $

'''High-level Gstreamer interface.
'''

import _ctypes
import ctypes
from ctypes import util
import time
import sys

import pyglet.lib

glib = pyglet.lib.load_library('glib-2.0')
gobject = pyglet.lib.load_library('gobject-2.0')
gst = pyglet.lib.load_library('gstreamer-0.10')
gstbase = pyglet.lib.load_library('gstbase-0.10')

GST_VERSION_MAJOR = 0
GST_VERSION_MINOR = 10
GST_VERSION_BUILD = 11

GST_STATE_VOID_PENDING = 0
GST_STATE_NULL = 1
GST_STATE_READY = 2
GST_STATE_PAUSED = 3
GST_STATE_PLAYING = 4

GST_FORMAT_TIME = 3


GstPluginInitFunc = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p)

GST_PADDING = 4
GST_PADDING_LARGE = 20



class GstPluginDesc(ctypes.Structure):
    _fields_ = [
        ('major_version', ctypes.c_int),
        ('minor_version', ctypes.c_int),
        ('name', ctypes.c_char_p),
        ('description', ctypes.c_char_p),
        ('plugin_init', GstPluginInitFunc),
        ('version', ctypes.c_char_p),
        ('license', ctypes.c_char_p),
        ('source', ctypes.c_char_p),
        ('package', ctypes.c_char_p),
        ('origin', ctypes.c_char_p),
        ('_gst_reserved', ctypes.c_void_p * GST_PADDING)
    ]

GType = ctypes.c_ulong
GBaseInitFunc = ctypes.CFUNCTYPE(None, ctypes.c_void_p)
GClassInitFunc = ctypes.CFUNCTYPE(None, ctypes.c_void_p, ctypes.c_void_p)
GInstanceInitFunc = ctypes.CFUNCTYPE(None, ctypes.c_void_p, ctypes.c_void_p)

class GTypeInfo(ctypes.Structure):
    _fields_ = [
        ('class_size', ctypes.c_uint16),
        ('base_init', GBaseInitFunc),
        ('base_finalize', ctypes.c_void_p),
        ('class_init', GClassInitFunc),
        ('class_finalize', ctypes.c_void_p),
        ('class_data', ctypes.c_void_p),
        ('instance_size', ctypes.c_uint16),
        ('n_preallocs', ctypes.c_uint16),
        ('instance_init', GInstanceInitFunc),
        ('value_table', ctypes.c_void_p),
    ]

class GTypeClass(ctypes.Structure):
    _fields_ = [
        ('g_type', GType),
    ]

class GTypeInstance(ctypes.Structure):
    _fields_ = [
        ('g_class', ctypes.POINTER(GTypeClass)),
    ]

class GObject(ctypes.Structure):
    _fields_ = [
        ('g_type_instance', GTypeInstance),
        ('ref_count', ctypes.c_uint),
        ('qdata', ctypes.c_void_p),
    ]

class GObjectClass(ctypes.Structure):
    _fields_ = [
        ('g_type_class', GTypeClass),
        ('construct_properties', ctypes.c_void_p),
        ('constructor', ctypes.c_void_p),
        ('set_property', ctypes.c_void_p),
        ('get_property', ctypes.c_void_p),
        ('dispose', ctypes.c_void_p),
        ('finalize', ctypes.c_void_p),
        ('dispatch_properties_changed', ctypes.c_void_p),
        ('notify', ctypes.c_void_p),
        ('pdummy', ctypes.c_void_p * 8),
    ]

class GstCaps(ctypes.Structure):
    _fields_ = [
        ('type', GType),
        ('refcount', ctypes.c_int),
        ('flags', ctypes.c_int),
        ('structs', ctypes.c_void_p),
        ('_gst_reserved', ctypes.c_void_p * GST_PADDING),
    ]

class GstStaticCaps(ctypes.Structure):
    _fields_ = [
        ('caps', GstCaps),
        ('string', ctypes.c_char_p),
        ('_gst_reserved', ctypes.c_void_p * GST_PADDING),
    ]

def GST_STATIC_CAPS(string):
    return GstStaticCaps(GstCaps(), string)

GST_PAD_ALWAYS = 0
GST_PAD_SOMETIMES = 1
GST_PAD_REQUEST = 2

GST_PAD_UNKNOWN = 0
GST_PAD_SRC = 1
GST_PAD_SINK = 2

class GstStaticPadTemplate(ctypes.Structure):
    _fields_ = [
        ('name_template', ctypes.c_char_p),
        ('direction', ctypes.c_int),
        ('presence', ctypes.c_int),
        ('static_caps', GstStaticCaps),
    ]


class GstObject(ctypes.Structure):
    _fields_ = [
        ('object', GObject),
        ('refcount', ctypes.c_int),
        ('lock', ctypes.c_void_p),
        ('name', ctypes.c_char_p),
        ('name_prefix', ctypes.c_char_p),
        ('parent', ctypes.c_void_p),
        ('flags', ctypes.c_uint32),
        ('_gst_reserved', ctypes.c_void_p),
    ]

class GstObjectClass(ctypes.Structure):
    _fields_ = [
        ('parent_class', GObjectClass),
        ('path_string_separator', ctypes.c_char),
        ('signal_object', ctypes.c_void_p),
        ('lock', ctypes.c_void_p),
        ('parent_set', ctypes.c_void_p),
        ('parent_unset', ctypes.c_void_p),
        ('object_saved', ctypes.c_void_p),
        ('deep_notify', ctypes.c_void_p),
        ('save_thyself', ctypes.c_void_p),
        ('restore_thyself', ctypes.c_void_p),
        ('_gst_reserved', ctypes.c_void_p * GST_PADDING)
    ]

class GstElementDetails(ctypes.Structure):
    _fields_ = [
        ('longname', ctypes.c_char_p),
        ('klass', ctypes.c_char_p),
        ('description', ctypes.c_char_p),
        ('author', ctypes.c_char_p),
        ('_gst_reserved', ctypes.c_void_p * GST_PADDING),
    ]

GstState = ctypes.c_int
GstStateChangeReturn = ctypes.c_int
GstClockTimeDiff = ctypes.c_int64
GstClockID = ctypes.c_void_p
GstClockTime = ctypes.c_uint64
GstClockTimeDiff = ctypes.c_int64

GST_CLOCK_TIME_NONE = -1

class GstMiniObject(ctypes.Structure):
    _fields_ = [
        ('instance', GTypeInstance),
        ('refcount', ctypes.c_int),
        ('flags', ctypes.c_uint),
        ('_gst_reserved', ctypes.c_void_p),
    ]

GST_EVENT_TYPE_SHIFT = 4
GST_EVENT_TYPE_UPSTREAM = 1 << 0
GST_EVENT_TYPE_DOWNSTREAM = 1 << 1
GST_EVENT_TYPE_SERIALIZED = 1 << 2

def GST_EVENT_MAKE_TYPE(num, flags):
    return (num << GST_EVENT_TYPE_SHIFT) | flags

GstEventType = ctypes.c_int
GST_EVENT_EOS = GST_EVENT_MAKE_TYPE(5, 
    GST_EVENT_TYPE_DOWNSTREAM | GST_EVENT_TYPE_SERIALIZED)

class GstEvent(ctypes.Structure):
    _fields_ = [
        ('mini_object', GstMiniObject),
        ('type', GstEventType),
        ('timestamp', ctypes.c_uint64),
        ('src', ctypes.c_void_p),
        ('structure', ctypes.c_void_p),
        ('_gst_reserved', ctypes.c_void_p),
    ]


class GstBuffer(ctypes.Structure):
    _fields_ = [
        ('mini_object', GstMiniObject),
        ('data', ctypes.c_void_p),
        ('size', ctypes.c_uint),
        ('timestamp', GstClockTime),
        ('duration', GstClockTime),
        ('caps', ctypes.POINTER(GstCaps)),
        ('offset', ctypes.c_uint64),
        ('offset_end', ctypes.c_uint64),
        ('malloc_data', ctypes.c_void_p),
        ('_gst_reserved', ctypes.c_void_p * GST_PADDING),
    ]

class GstSegment(ctypes.Structure):
    _fields_ = [
        ('rate', ctypes.c_double),
        ('abs_rate', ctypes.c_double),
        ('format', ctypes.c_int),
        ('flags', ctypes.c_int),
        ('start', ctypes.c_int64),
        ('stop', ctypes.c_int64),
        ('time', ctypes.c_int64),
        ('accum', ctypes.c_int64),
        ('last_stop', ctypes.c_int64),
        ('duration', ctypes.c_int64),
        #('applied_rate', ctypes.c_double),  ABI added 0.10.6
        ('_gst_reserved', ctypes.c_void_p * GST_PADDING),
    ]

#TEMP
GstPad = ctypes.c_int

GstFlowReturn = ctypes.c_int
GstPadChainFunction = ctypes.CFUNCTYPE(GstFlowReturn, 
    ctypes.POINTER(GstPad), ctypes.POINTER(GstBuffer))
GstPadSetCapsFunction = ctypes.CFUNCTYPE(ctypes.c_int,
    ctypes.POINTER(GstPad), ctypes.POINTER(GstCaps))
GstPadEventFunction = ctypes.CFUNCTYPE(ctypes.c_int,
    ctypes.POINTER(GstPad), ctypes.POINTER(GstEvent))
GstPadGetRangeFunction = ctypes.CFUNCTYPE(ctypes.c_int,
    ctypes.POINTER(GstPad), ctypes.c_uint64, ctypes.c_uint,
    ctypes.POINTER(ctypes.POINTER(GstBuffer)))


class GstElement(ctypes.Structure):
    _fields_ = [
        ('object', GstObject),
        ('state_lock', ctypes.c_void_p),
        ('state_cond', ctypes.c_void_p),
        ('state_cookie', ctypes.c_uint32),
        ('current_state', GstState),
        ('next_state', GstState),
        ('pending_state', GstState),
        ('last_return', GstStateChangeReturn),
        ('bus', ctypes.c_void_p),
        ('clock', ctypes.c_void_p),
        ('base_time', GstClockTimeDiff),
        ('numpads', ctypes.c_uint16),
        ('pads', ctypes.c_void_p),
        ('numsrcpads', ctypes.c_uint16),
        ('srcpads', ctypes.c_void_p),
        ('numsinkpads', ctypes.c_uint16),
        ('sinkpads', ctypes.c_void_p),
        ('pads_cookie', ctypes.c_uint32),
        ('_gst_reserved', ctypes.c_void_p * GST_PADDING),
    ]

class GstElementClass(ctypes.Structure):
    _fields_ = [
        ('parent_class', GstObjectClass),
        ('details', GstElementDetails),
        ('elementfactory', ctypes.c_void_p),
        ('padtemplates', ctypes.c_void_p),
        ('numpadtemplates', ctypes.c_int),
        ('pad_templ_cookie', ctypes.c_uint32),
        ('pad_added', ctypes.c_void_p),
        ('pad_removed', ctypes.c_void_p),
        ('no_more_pads', ctypes.c_void_p),
        ('request_new_pad', ctypes.c_void_p),
        ('release_pad', ctypes.c_void_p),
        ('get_state', ctypes.c_void_p),
        ('set_state', ctypes.c_void_p),
        ('change_state', ctypes.c_void_p),
        ('set_bus', ctypes.c_void_p),
        ('provide_clock', ctypes.c_void_p),
        ('set_clock', ctypes.c_void_p),
        ('get_index', ctypes.c_void_p),
        ('set_index', ctypes.c_void_p),
        ('send_event', ctypes.c_void_p),
        ('get_query_types', ctypes.c_void_p),
        ('query', ctypes.c_void_p),
        ('_gst_reserved', ctypes.c_void_p * GST_PADDING),
    ]

class GstBaseSrc(ctypes.Structure):
    _fields_ = [
        ('element', GstElement),
        ('srcpad', ctypes.POINTER(GstPad)),
        ('live_lock', ctypes.c_void_p),
        ('live_cond', ctypes.c_void_p),
        ('is_live', ctypes.c_int),
        ('live_running', ctypes.c_int),
        ('blocksize', ctypes.c_int),
        ('can_activate_push', ctypes.c_int),
        ('pad_mode', ctypes.c_int),
        ('seekable', ctypes.c_int),
        ('random_access', ctypes.c_int),
        ('clock_id', GstClockID),
        ('end_time', GstClockTime),
        ('segment', GstSegment),
        ('need_newsegment', ctypes.c_int),
        ('offset', ctypes.c_uint64),
        ('size', ctypes.c_uint64),
        ('num_buffers', ctypes.c_int),
        ('num_buffers_left', ctypes.c_int),
        ('_gst_reserved', ctypes.c_void_p * (GST_PADDING_LARGE - 1)),
        ('priv', ctypes.c_void_p),
    ]

class GstBaseSrcClass(ctypes.Structure):
    _fields_ = [
        ('parent_class', GstElementClass),
        ('get_caps', ctypes.CFUNCTYPE(ctypes.c_void_p,
            ctypes.POINTER(GstBaseSrc))),
        ('set_caps', ctypes.CFUNCTYPE(ctypes.c_int,
            ctypes.POINTER(GstBaseSrc), ctypes.c_void_p)),
        ('negotiate', ctypes.CFUNCTYPE(ctypes.c_int,
            ctypes.POINTER(GstBaseSrc))),
        ('newsegment', ctypes.CFUNCTYPE(ctypes.c_int,
            ctypes.POINTER(GstBaseSrc))),
        ('start', ctypes.CFUNCTYPE(ctypes.c_int,
            ctypes.POINTER(GstBaseSrc))),
        ('stop', ctypes.CFUNCTYPE(ctypes.c_int,
            ctypes.POINTER(GstBaseSrc))),
        ('get_times', ctypes.CFUNCTYPE(None,
            ctypes.POINTER(GstBaseSrc), ctypes.POINTER(GstBuffer),
            ctypes.POINTER(GstClockTime), ctypes.POINTER(GstClockTime))),
        ('get_size', ctypes.CFUNCTYPE(ctypes.c_int,
            ctypes.POINTER(GstBaseSrc), ctypes.POINTER(ctypes.c_uint64))),
        ('is_seekable', ctypes.CFUNCTYPE(ctypes.c_int,
            ctypes.POINTER(GstBaseSrc))),
        ('unlock', ctypes.CFUNCTYPE(ctypes.c_int,
            ctypes.POINTER(GstBaseSrc))),
        ('event', ctypes.CFUNCTYPE(ctypes.c_int,
            ctypes.POINTER(GstBaseSrc), ctypes.POINTER(GstEvent))),
        ('create', ctypes.CFUNCTYPE(GstFlowReturn,
            ctypes.POINTER(GstBaseSrc), ctypes.c_uint64, ctypes.c_uint, 
            ctypes.POINTER(ctypes.POINTER(GstBuffer)))),
        ('do_seek', ctypes.CFUNCTYPE(ctypes.c_int,
            ctypes.POINTER(GstBaseSrc), ctypes.POINTER(GstSegment))),
        ('query', ctypes.c_void_p),
        ('check_get_range', ctypes.c_void_p),
        ('fixate', ctypes.c_void_p),
        ('_gst_reserved', (ctypes.c_void_p * (GST_PADDING_LARGE - 4))),
    ]


class GstBaseSink(ctypes.Structure):
    _fields_ = [
        ('element', GstElement),
        ('sinkpad', ctypes.c_void_p),
        ('pad_mode', ctypes.c_int),
        ('offset', ctypes.c_uint64),
        ('can_activate_pull', ctypes.c_int),
        ('can_activate_push', ctypes.c_int),
        ('preroll_queue', ctypes.c_void_p),
        ('preroll_queue_max_len', ctypes.c_int),
        ('preroll_queued', ctypes.c_int),
        ('buffers_queued', ctypes.c_int),
        ('events_queued', ctypes.c_int),
        ('eos', ctypes.c_int),
        ('eos_queued', ctypes.c_int),
        ('need_preroll', ctypes.c_int),
        ('have_preroll', ctypes.c_int),
        ('playing_async', ctypes.c_int),
        ('have_newsegment', ctypes.c_int),
        ('segment', GstSegment),
        ('clock_id', ctypes.c_void_p),
        ('end_time', ctypes.c_uint64),
        ('sync', ctypes.c_int),
        ('flushing', ctypes.c_int),
        ('_gst_reserved', ctypes.c_void_p * (GST_PADDING_LARGE - 1)),
        ('priv', ctypes.c_void_p),
    ]

class GstBaseSinkClass(ctypes.Structure):
    _fields_ = [
        ('parent_class', GstElementClass),
        ('get_caps', ctypes.c_void_p),
        ('set_caps', ctypes.c_void_p),
        ('buffer_alloc', ctypes.c_void_p),
        ('get_times', ctypes.c_void_p),
        ('start', ctypes.c_void_p),
        ('stop', ctypes.c_void_p),
        ('unlock', ctypes.c_void_p),
        ('event', ctypes.c_void_p),
        ('preroll', ctypes.c_void_p),
        ('render', ctypes.c_void_p),
        ('async_play', ctypes.c_void_p),
        ('_gst_reserved', ctypes.c_void_p * (GST_PADDING_LARGE - 1)),
    ]

class GstBaseAudioSink(ctypes.Structure):
    _fields_ = [
        ('element', GstBaseSink),
        ('ringbuffer', ctypes.c_void_p),
        ('buffer_time', ctypes.c_uint64),
        ('latency_time', ctypes.c_uint64),
        ('next_sample', ctypes.c_uint64),
        ('provide_clock', ctypes.c_int),
        ('provided_clock', ctypes.c_void_p),
        ('_gst_reserved', ctypes.c_void_p * GST_PADDING),
    ]

class GstBaseAudioSinkClass(ctypes.Structure):
    _fields_ = [
        ('parent_class', GstBaseSinkClass),
        ('create_ringbuffer', ctypes.c_void_p),
        ('_gst_reserved', ctypes.c_void_p * GST_PADDING),
    ]

class GstAudioSink(ctypes.Structure):
    _fields_ = [
        ('element', GstBaseAudioSink),
        ('thread', ctypes.c_void_p),
        ('_gst_reserved', ctypes.c_void_p * GST_PADDING),
    ]

GstBufferFormatType = ctypes.c_int
GST_BUFFER_LINEAR = 0
GST_BUFFER_FLOAT = 1
GST_BUFFER_MU_LAW = 2
GST_BUFFER_A_LAW = 3
GST_BUFFER_IMA_ADPCM = 4
GST_BUFFER_MPEG = 5
GST_BUFFER_GSM = 6

GstBufferFormat = ctypes.c_int

class GstRingBufferSpec(ctypes.Structure):
    _fields_ = [
        ('caps', ctypes.POINTER(GstCaps)),
        ('type', GstBufferFormatType),
        ('format', GstBufferFormat),
        ('sign', ctypes.c_int),
        ('bigend', ctypes.c_int),
        ('width', ctypes.c_int),
        ('depth', ctypes.c_int),
        ('rate', ctypes.c_int),
        ('channels', ctypes.c_int),
        ('latency_time', ctypes.c_uint64),
        ('buffer_time', ctypes.c_uint64),
        ('segsize', ctypes.c_int),
        ('segtotal', ctypes.c_int),
        ('bytes_per_sample', ctypes.c_int),
        ('silence_sample', ctypes.c_uint8 * 32),
        ('_gst_reserved', ctypes.c_void_p * GST_PADDING),
    ]

class GstAudioSinkClass(ctypes.Structure):
    _fields_ = [
        ('parent_class', GstBaseAudioSinkClass),
        ('open', ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p)),
        ('prepare', ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p,
                                     ctypes.POINTER(GstRingBufferSpec))),
        ('unprepare', ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p)),
        ('close', ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p)),
        ('write', ctypes.CFUNCTYPE(ctypes.c_uint, ctypes.c_void_p,
                                   ctypes.POINTER(ctypes.c_byte), ctypes.c_uint)),
        ('delay', ctypes.CFUNCTYPE(ctypes.c_uint, ctypes.c_void_p)),
        ('reset', ctypes.CFUNCTYPE(None, ctypes.c_void_p)),
        ('_gst_reserved', ctypes.c_void_p * GST_PADDING),
    ]

def py_derived_element(base):
    '''Derive a GstElement with a 'pyobject' member.'''
    class PyGstElement(ctypes.Structure):
        _fields_ = [
            ('element', base),
            ('pyobject', ctypes.py_object),
        ]
    return PyGstElement

GstMessageType = ctypes.c_int
GST_MESSAGE_UNKNOWN           = 0
GST_MESSAGE_EOS               = (1 << 0)
GST_MESSAGE_ERROR             = (1 << 1)
GST_MESSAGE_WARNING           = (1 << 2)
GST_MESSAGE_INFO              = (1 << 3)
GST_MESSAGE_TAG               = (1 << 4)
GST_MESSAGE_BUFFERING         = (1 << 5)
GST_MESSAGE_STATE_CHANGED     = (1 << 6)
GST_MESSAGE_STATE_DIRTY       = (1 << 7)
GST_MESSAGE_STEP_DONE         = (1 << 8)
GST_MESSAGE_CLOCK_PROVIDE     = (1 << 9)
GST_MESSAGE_CLOCK_LOST        = (1 << 10)
GST_MESSAGE_NEW_CLOCK         = (1 << 11)
GST_MESSAGE_STRUCTURE_CHANGE  = (1 << 12)
GST_MESSAGE_STREAM_STATUS     = (1 << 13)
GST_MESSAGE_APPLICATION       = (1 << 14)
GST_MESSAGE_ELEMENT           = (1 << 15)
GST_MESSAGE_SEGMENT_START     = (1 << 16)
GST_MESSAGE_SEGMENT_DONE      = (1 << 17)
GST_MESSAGE_DURATION          = (1 << 18)
GST_MESSAGE_LATENCY           = (1 << 19)
GST_MESSAGE_ANY               = ~0

class GstMessage(ctypes.Structure):
    _fields_ = [
        ('mini_object', GstMiniObject),
        ('lock', ctypes.c_void_p),
        ('cond', ctypes.c_void_p),
        ('type', GstMessageType),
        ('timestamp', ctypes.c_uint64),
        ('src', ctypes.c_void_p),
        ('structure', ctypes.c_void_p),
        ('_gst_reserved', ctypes.c_void_p * GST_PADDING),
    ]

gst.gst_bus_poll.restype = ctypes.POINTER(GstMessage)

# Python high-level classes

class Plugin(object):
    name = ''
    description = ''
    version = ''
    license = ''
    source = ''
    package = ''
    origin = ''

    elements = ()

    def register(self):
        self._func_plugin_init = GstPluginInitFunc(self.plugin_init)
        desc = GstPluginDesc(GST_VERSION_MAJOR,
                             GST_VERSION_MINOR,
                             self.name,
                             self.description,
                             self._func_plugin_init,
                             self.version,
                             self.license,
                             self.source,
                             self.package,
                             self.origin)
        gst._gst_plugin_register_static(ctypes.byref(desc))

    def plugin_init(self, plugin):
        for element_class in self.elements:
            gst.gst_element_register(plugin, 
                                     element_class.name, 
                                     element_class.rank,
                                     element_class.get_type())
        return True

def virtual_method(cls, name, ctype):
    def f(this, *args):
        this = ctypes.cast(this, ctypes.POINTER(cls._py_instance_struct))
        self = this.contents.pyobject
        if hasattr(self, name):
            return getattr(self, name)(*args)
    return ctype(f)

class Element(object):
    # You must fill these in as class attributes
    name = 'Unknown Element'
    klass = 'Element/Unknown'
    description = 'An unknown element'
    author = 'nobody@nowhere.com'

    class_struct = GstElementClass
    instance_struct = GstElement
    parent_type = gst.gst_element_get_type

    rank = 0
    pad_templates = ()
    pad_instances = ()
    vmethods = ()

    # Filled in on class
    _object_type = None
    _py_instance_struct = None
    _funcptrs = []

    # Filled in upon instantiation
    pads = None

    @classmethod
    def get_type(cls):
        if cls._object_type is None:
            cls._py_instance_struct = py_derived_element(cls.instance_struct)
            cls._func_baseinit = GBaseInitFunc(cls.base_init)
            cls._func_classinit = GClassInitFunc(cls.class_init)
            cls._func_instanceinit = GInstanceInitFunc(cls.instance_init)
            object_info = GTypeInfo(ctypes.sizeof(cls.class_struct),
                                    cls._func_baseinit,
                                    None,
                                    cls._func_classinit,
                                    None,
                                    None,
                                    ctypes.sizeof(cls._py_instance_struct),
                                    0,
                                    cls._func_instanceinit)
            cls._object_type = gobject.g_type_register_static(
                    cls.parent_type(),
                    cls.__name__, 
                    ctypes.byref(object_info), 
                    0)
        return cls._object_type

    @classmethod
    def base_init(cls, _class):
        details = GstElementDetails(cls.name,
                                    cls.klass,
                                    cls.description,
                                    cls.author)
        gst.gst_element_class_set_details(_class, ctypes.byref(details))

        # register pad templates
        gst.gst_static_pad_template_get.restype = ctypes.c_void_p
        for pad_template in cls.pad_templates:
            template_struct = GstStaticPadTemplate(
                pad_template.name,
                pad_template.direction,
                GST_PAD_ALWAYS,
                GST_STATIC_CAPS(pad_template.caps))
            template = gst.gst_static_pad_template_get(
                ctypes.byref(template_struct))
            gst.gst_element_class_add_pad_template(_class, template)

    @classmethod
    def class_init(cls, _class, data):
        # install properties, signals, set up vmethods.
        _class = ctypes.cast(_class, ctypes.POINTER(cls.class_struct))

        for method in cls.vmethods:
            c = _class.contents
            while not hasattr(c, method):
               c = c.parent_class
            
            for name, ftype in c._fields_:
                if name == method:
                    break
            funcptr = virtual_method(cls, method, ftype)
            cls._funcptrs.append(funcptr)

            setattr(c, method, funcptr)

    @classmethod
    def instance_init(cls, _instance, _class):
        self = cls()

        _pyinstance = ctypes.cast(_instance, 
                                  ctypes.POINTER(cls._py_instance_struct))
        _pyinstance.contents.pyobject = self
        _ctypes.Py_INCREF(_pyinstance.contents.pyobject)

        self.pads = []
        for name, pad_class in cls.pad_instances:
            pad = gst.gst_pad_new_from_template(
                gst.gst_element_class_get_pad_template(
                    _class, pad_class.name), 
                name)

            pypad = pad_class(pad, self)
            self.pads.append(pypad)

            gst.gst_element_add_pad(_instance, pad)

    @classmethod
    def get_instance(cls, element):
        _pyinstance = ctypes.cast(element,
                                  ctypes.POINTER(cls._py_instance_struct))
        return _pyinstance.contents.pyobject

    def __del__(self):
        print 'kill me!'

GST_FLOW_RESENT = 1
GST_FLOW_OK = 0
GST_FLOW_NOT_LINKED = -1
GST_FLOW_WRONG_STATE = -2
GST_FLOW_UNEXPECTED = -3
GST_FLOW_NOT_NEGOTIATED = -4
GST_FLOW_ERROR = -5
GST_FLOW_NOT_SUPPORTED = -6

class Pad(object):
    # Fill these in as class attributes
    name = ''
    direction = GST_PAD_UNKNOWN
    caps = 'ANY'

    # Filled in on instance
    this = None
    element = None
    funcptrs = [] # for refcounting

    def __init__(self, this, element):
        self.this = this
        self.element = element

        if hasattr(self, 'setcaps'):
            setcaps = GstPadSetCapsFunction(self.setcaps)
            gst.gst_pad_set_setcaps_function(this, setcaps)
            self.funcptrs.append(setcaps)

        if hasattr(self, 'chain'):
            chain = GstPadChainFunction(self.chain)
            gst.gst_pad_set_chain_function(this, chain)
            self.funcptrs.append(chain)

        if hasattr(self, 'event'):
            event = GstPadEventFunction(self.event)
            gst.gst_pad_set_event_function(this, event)
            self.funcptrs.append(event)

        if hasattr(self, 'get_range'):
            get_range = GstPadGetRangeFunction(self.get_range)
            gst.gst_pad_set_get_range_function(this, get_range)
            self.funcptrs.append(get_range)

    #def setcaps(self, this, caps):
    #    return True

    #def chain(self, this, buf):
    #    return GST_FLOW_OK

    #def event(self, this, event):
    #    return False

mainloop = None
mainloop_context = None

def init():
    global mainloop, mainloop_context
    glib.g_main_loop_new.restype = ctypes.c_void_p
    mainloop = glib.g_main_loop_new(None, True)
    glib.g_main_loop_get_context.restype = ctypes.c_void_p
    mainloop_context = glib.g_main_loop_get_context(mainloop)

    argc = ctypes.c_int(0)
    argv = ctypes.c_void_p()
    gst.gst_init(ctypes.byref(argc), ctypes.byref(argv))

def play(uri):
    playbin = gst.gst_element_factory_make('playbin', 'play')
    gobject.g_object_set(playbin, 'uri', uri, None)

    gst.gst_element_set_state(playbin, GST_STATE_PLAYING)

def heartbeat():
    glib.g_main_context_iteration(mainloop_context, False)


