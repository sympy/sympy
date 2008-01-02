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
'''Wrapper for asound

Generated with:
tools/wraptypes/wrap.py -o asound.py -lasound /usr/include/alsa/alisp.h /usr/include/alsa/asoundef.h /usr/include/alsa/asoundlib.h /usr/include/alsa/conf.h /usr/include/alsa/control.h /usr/include/alsa/control_external.h /usr/include/alsa/conv.h /usr/include/alsa/error.h /usr/include/alsa/global.h /usr/include/alsa/hwdep.h /usr/include/alsa/iatomic.h /usr/include/alsa/input.h /usr/include/alsa/instr.h /usr/include/alsa/mixer.h /usr/include/alsa/mixer_abst.h /usr/include/alsa/pcm.h /usr/include/alsa/pcm_external.h /usr/include/alsa/pcm_extplug.h /usr/include/alsa/pcm_ioplug.h /usr/include/alsa/pcm_old.h /usr/include/alsa/pcm_plugin.h /usr/include/alsa/pcm_rate.h /usr/include/alsa/rawmidi.h /usr/include/alsa/seq.h /usr/include/alsa/seq_event.h /usr/include/alsa/seq_midi_event.h /usr/include/alsa/seqmid.h /usr/include/alsa/timer.h /usr/include/alsa/version.h

 -- And then hacked to work with libasound.so, grep for XXX  

Do not regenerate this file.
'''

__docformat__ =  'restructuredtext'
__version__ = '$Id: asound.py 1286 2007-09-29 02:33:49Z Alex.Holkner $'

import ctypes
from ctypes import *

import pyglet.lib

_lib = pyglet.lib.load_library('asound')

_int_types = (c_int16, c_int32)
if hasattr(ctypes, 'c_int64'):
    # Some builds of ctypes apparently do not have c_int64
    # defined; it's a pretty good bet that these builds do not
    # have 64-bit pointers.
    _int_types += (ctypes.c_int64,)
for t in _int_types:
    if sizeof(t) == sizeof(c_size_t):
        c_ptrdiff_t = t

class c_void(Structure):
    # c_void_p is a buggy return type, converting to int, so
    # POINTER(None) == c_void_p is actually written as
    # POINTER(c_void), so it can be treated as a real pointer.
    _fields_ = [('dummy', c_int)]



class struct_alisp_cfg(Structure):
    __slots__ = [
    ]
struct_alisp_cfg._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/alisp.h:39
alsa_lisp_default_cfg_free = _lib.alsa_lisp_default_cfg_free
alsa_lisp_default_cfg_free.restype = None
alsa_lisp_default_cfg_free.argtypes = [POINTER(struct_alisp_cfg)]

class struct_alisp_cfg(Structure):
    __slots__ = [
    ]
struct_alisp_cfg._fields_ = [
    ('_opaque_struct', c_int)
]

class struct_alisp_instance(Structure):
    __slots__ = [
    ]
struct_alisp_instance._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/alisp.h:40
alsa_lisp = _lib.alsa_lisp
alsa_lisp.restype = c_int
alsa_lisp.argtypes = [POINTER(struct_alisp_cfg), POINTER(POINTER(struct_alisp_instance))]

class struct_alisp_instance(Structure):
    __slots__ = [
    ]
struct_alisp_instance._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/alisp.h:41
alsa_lisp_free = _lib.alsa_lisp_free
alsa_lisp_free.restype = None
alsa_lisp_free.argtypes = [POINTER(struct_alisp_instance)]

class struct_alisp_instance(Structure):
    __slots__ = [
    ]
struct_alisp_instance._fields_ = [
    ('_opaque_struct', c_int)
]

class struct_alisp_seq_iterator(Structure):
    __slots__ = [
    ]
struct_alisp_seq_iterator._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/alisp.h:48
alsa_lisp_result_free = _lib.alsa_lisp_result_free
alsa_lisp_result_free.restype = None
alsa_lisp_result_free.argtypes = [POINTER(struct_alisp_instance), POINTER(struct_alisp_seq_iterator)]

class struct_alisp_instance(Structure):
    __slots__ = [
    ]
struct_alisp_instance._fields_ = [
    ('_opaque_struct', c_int)
]

class struct_alisp_seq_iterator(Structure):
    __slots__ = [
    ]
struct_alisp_seq_iterator._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/alisp.h:50
alsa_lisp_seq_first = _lib.alsa_lisp_seq_first
alsa_lisp_seq_first.restype = c_int
alsa_lisp_seq_first.argtypes = [POINTER(struct_alisp_instance), c_char_p, POINTER(POINTER(struct_alisp_seq_iterator))]

class struct_alisp_seq_iterator(Structure):
    __slots__ = [
    ]
struct_alisp_seq_iterator._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/alisp.h:52
alsa_lisp_seq_next = _lib.alsa_lisp_seq_next
alsa_lisp_seq_next.restype = c_int
alsa_lisp_seq_next.argtypes = [POINTER(POINTER(struct_alisp_seq_iterator))]

class struct_alisp_seq_iterator(Structure):
    __slots__ = [
    ]
struct_alisp_seq_iterator._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/alisp.h:53
alsa_lisp_seq_count = _lib.alsa_lisp_seq_count
alsa_lisp_seq_count.restype = c_int
alsa_lisp_seq_count.argtypes = [POINTER(struct_alisp_seq_iterator)]

class struct_alisp_seq_iterator(Structure):
    __slots__ = [
    ]
struct_alisp_seq_iterator._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/alisp.h:54
alsa_lisp_seq_integer = _lib.alsa_lisp_seq_integer
alsa_lisp_seq_integer.restype = c_int
alsa_lisp_seq_integer.argtypes = [POINTER(struct_alisp_seq_iterator), POINTER(c_long)]

class struct_alisp_seq_iterator(Structure):
    __slots__ = [
    ]
struct_alisp_seq_iterator._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/alisp.h:55
alsa_lisp_seq_pointer = _lib.alsa_lisp_seq_pointer
alsa_lisp_seq_pointer.restype = c_int
alsa_lisp_seq_pointer.argtypes = [POINTER(struct_alisp_seq_iterator), c_char_p, POINTER(POINTER(None))]

IEC958_AES0_PROFESSIONAL = 1 	# /usr/include/alsa/asoundef.h:96
IEC958_AES0_NONAUDIO = 2 	# /usr/include/alsa/asoundef.h:97
IEC958_AES0_PRO_EMPHASIS = 28 	# /usr/include/alsa/asoundef.h:98
IEC958_AES0_PRO_EMPHASIS_NOTID = 0 	# /usr/include/alsa/asoundef.h:99
IEC958_AES0_PRO_EMPHASIS_NONE = 4 	# /usr/include/alsa/asoundef.h:100
IEC958_AES0_PRO_EMPHASIS_5015 = 12 	# /usr/include/alsa/asoundef.h:101
IEC958_AES0_PRO_EMPHASIS_CCITT = 28 	# /usr/include/alsa/asoundef.h:102
IEC958_AES0_PRO_FREQ_UNLOCKED = 32 	# /usr/include/alsa/asoundef.h:103
IEC958_AES0_PRO_FS = 192 	# /usr/include/alsa/asoundef.h:104
IEC958_AES0_PRO_FS_NOTID = 0 	# /usr/include/alsa/asoundef.h:105
IEC958_AES0_PRO_FS_44100 = 64 	# /usr/include/alsa/asoundef.h:106
IEC958_AES0_PRO_FS_48000 = 128 	# /usr/include/alsa/asoundef.h:107
IEC958_AES0_PRO_FS_32000 = 192 	# /usr/include/alsa/asoundef.h:108
IEC958_AES0_CON_NOT_COPYRIGHT = 4 	# /usr/include/alsa/asoundef.h:109
IEC958_AES0_CON_EMPHASIS = 56 	# /usr/include/alsa/asoundef.h:110
IEC958_AES0_CON_EMPHASIS_NONE = 0 	# /usr/include/alsa/asoundef.h:111
IEC958_AES0_CON_EMPHASIS_5015 = 8 	# /usr/include/alsa/asoundef.h:112
IEC958_AES0_CON_MODE = 192 	# /usr/include/alsa/asoundef.h:113
IEC958_AES1_PRO_MODE = 15 	# /usr/include/alsa/asoundef.h:114
IEC958_AES1_PRO_MODE_NOTID = 0 	# /usr/include/alsa/asoundef.h:115
IEC958_AES1_PRO_MODE_STEREOPHONIC = 2 	# /usr/include/alsa/asoundef.h:116
IEC958_AES1_PRO_MODE_SINGLE = 4 	# /usr/include/alsa/asoundef.h:117
IEC958_AES1_PRO_MODE_TWO = 8 	# /usr/include/alsa/asoundef.h:118
IEC958_AES1_PRO_MODE_PRIMARY = 12 	# /usr/include/alsa/asoundef.h:119
IEC958_AES1_PRO_MODE_BYTE3 = 15 	# /usr/include/alsa/asoundef.h:120
IEC958_AES1_PRO_USERBITS = 240 	# /usr/include/alsa/asoundef.h:121
IEC958_AES1_PRO_USERBITS_NOTID = 0 	# /usr/include/alsa/asoundef.h:122
IEC958_AES1_PRO_USERBITS_192 = 128 	# /usr/include/alsa/asoundef.h:123
IEC958_AES1_PRO_USERBITS_UDEF = 192 	# /usr/include/alsa/asoundef.h:124
IEC958_AES1_CON_CATEGORY = 127 	# /usr/include/alsa/asoundef.h:125
IEC958_AES1_CON_GENERAL = 0 	# /usr/include/alsa/asoundef.h:126
IEC958_AES1_CON_EXPERIMENTAL = 64 	# /usr/include/alsa/asoundef.h:127
IEC958_AES1_CON_SOLIDMEM_MASK = 15 	# /usr/include/alsa/asoundef.h:128
IEC958_AES1_CON_SOLIDMEM_ID = 8 	# /usr/include/alsa/asoundef.h:129
IEC958_AES1_CON_BROADCAST1_MASK = 7 	# /usr/include/alsa/asoundef.h:130
IEC958_AES1_CON_BROADCAST1_ID = 4 	# /usr/include/alsa/asoundef.h:131
IEC958_AES1_CON_DIGDIGCONV_MASK = 7 	# /usr/include/alsa/asoundef.h:132
IEC958_AES1_CON_DIGDIGCONV_ID = 2 	# /usr/include/alsa/asoundef.h:133
IEC958_AES1_CON_ADC_COPYRIGHT_MASK = 31 	# /usr/include/alsa/asoundef.h:134
IEC958_AES1_CON_ADC_COPYRIGHT_ID = 6 	# /usr/include/alsa/asoundef.h:135
IEC958_AES1_CON_ADC_MASK = 31 	# /usr/include/alsa/asoundef.h:136
IEC958_AES1_CON_ADC_ID = 22 	# /usr/include/alsa/asoundef.h:137
IEC958_AES1_CON_BROADCAST2_MASK = 15 	# /usr/include/alsa/asoundef.h:138
IEC958_AES1_CON_BROADCAST2_ID = 14 	# /usr/include/alsa/asoundef.h:139
IEC958_AES1_CON_LASEROPT_MASK = 7 	# /usr/include/alsa/asoundef.h:140
IEC958_AES1_CON_LASEROPT_ID = 1 	# /usr/include/alsa/asoundef.h:141
IEC958_AES1_CON_MUSICAL_MASK = 7 	# /usr/include/alsa/asoundef.h:142
IEC958_AES1_CON_MUSICAL_ID = 5 	# /usr/include/alsa/asoundef.h:143
IEC958_AES1_CON_MAGNETIC_MASK = 7 	# /usr/include/alsa/asoundef.h:144
IEC958_AES1_CON_MAGNETIC_ID = 3 	# /usr/include/alsa/asoundef.h:145
IEC958_AES1_CON_IEC908_CD = 1 	# /usr/include/alsa/asoundef.h:146
IEC958_AES1_CON_NON_IEC908_CD = 9 	# /usr/include/alsa/asoundef.h:147
IEC958_AES1_CON_PCM_CODER = 2 	# /usr/include/alsa/asoundef.h:148
IEC958_AES1_CON_SAMPLER = 34 	# /usr/include/alsa/asoundef.h:149
IEC958_AES1_CON_MIXER = 18 	# /usr/include/alsa/asoundef.h:150
IEC958_AES1_CON_RATE_CONVERTER = 26 	# /usr/include/alsa/asoundef.h:151
IEC958_AES1_CON_SYNTHESIZER = 5 	# /usr/include/alsa/asoundef.h:152
IEC958_AES1_CON_MICROPHONE = 13 	# /usr/include/alsa/asoundef.h:153
IEC958_AES1_CON_DAT = 3 	# /usr/include/alsa/asoundef.h:154
IEC958_AES1_CON_VCR = 11 	# /usr/include/alsa/asoundef.h:155
IEC958_AES1_CON_ORIGINAL = 128 	# /usr/include/alsa/asoundef.h:156
IEC958_AES2_PRO_SBITS = 7 	# /usr/include/alsa/asoundef.h:157
IEC958_AES2_PRO_SBITS_20 = 2 	# /usr/include/alsa/asoundef.h:158
IEC958_AES2_PRO_SBITS_24 = 4 	# /usr/include/alsa/asoundef.h:159
IEC958_AES2_PRO_SBITS_UDEF = 6 	# /usr/include/alsa/asoundef.h:160
IEC958_AES2_PRO_WORDLEN = 56 	# /usr/include/alsa/asoundef.h:161
IEC958_AES2_PRO_WORDLEN_NOTID = 0 	# /usr/include/alsa/asoundef.h:162
IEC958_AES2_PRO_WORDLEN_22_18 = 16 	# /usr/include/alsa/asoundef.h:163
IEC958_AES2_PRO_WORDLEN_23_19 = 32 	# /usr/include/alsa/asoundef.h:164
IEC958_AES2_PRO_WORDLEN_24_20 = 40 	# /usr/include/alsa/asoundef.h:165
IEC958_AES2_PRO_WORDLEN_20_16 = 48 	# /usr/include/alsa/asoundef.h:166
IEC958_AES2_CON_SOURCE = 15 	# /usr/include/alsa/asoundef.h:167
IEC958_AES2_CON_SOURCE_UNSPEC = 0 	# /usr/include/alsa/asoundef.h:168
IEC958_AES2_CON_CHANNEL = 240 	# /usr/include/alsa/asoundef.h:169
IEC958_AES2_CON_CHANNEL_UNSPEC = 0 	# /usr/include/alsa/asoundef.h:170
IEC958_AES3_CON_FS = 15 	# /usr/include/alsa/asoundef.h:171
IEC958_AES3_CON_FS_44100 = 0 	# /usr/include/alsa/asoundef.h:172
IEC958_AES3_CON_FS_48000 = 2 	# /usr/include/alsa/asoundef.h:173
IEC958_AES3_CON_FS_32000 = 3 	# /usr/include/alsa/asoundef.h:174
IEC958_AES3_CON_CLOCK = 48 	# /usr/include/alsa/asoundef.h:175
IEC958_AES3_CON_CLOCK_1000PPM = 0 	# /usr/include/alsa/asoundef.h:176
IEC958_AES3_CON_CLOCK_50PPM = 16 	# /usr/include/alsa/asoundef.h:177
IEC958_AES3_CON_CLOCK_VARIABLE = 32 	# /usr/include/alsa/asoundef.h:178
MIDI_CHANNELS = 16 	# /usr/include/alsa/asoundef.h:188
MIDI_GM_DRUM_CHANNEL = 9 	# /usr/include/alsa/asoundef.h:189
MIDI_CMD_NOTE_OFF = 128 	# /usr/include/alsa/asoundef.h:197
MIDI_CMD_NOTE_ON = 144 	# /usr/include/alsa/asoundef.h:198
MIDI_CMD_NOTE_PRESSURE = 160 	# /usr/include/alsa/asoundef.h:199
MIDI_CMD_CONTROL = 176 	# /usr/include/alsa/asoundef.h:200
MIDI_CMD_PGM_CHANGE = 192 	# /usr/include/alsa/asoundef.h:201
MIDI_CMD_CHANNEL_PRESSURE = 208 	# /usr/include/alsa/asoundef.h:202
MIDI_CMD_BENDER = 224 	# /usr/include/alsa/asoundef.h:203
MIDI_CMD_COMMON_SYSEX = 240 	# /usr/include/alsa/asoundef.h:205
MIDI_CMD_COMMON_MTC_QUARTER = 241 	# /usr/include/alsa/asoundef.h:206
MIDI_CMD_COMMON_SONG_POS = 242 	# /usr/include/alsa/asoundef.h:207
MIDI_CMD_COMMON_SONG_SELECT = 243 	# /usr/include/alsa/asoundef.h:208
MIDI_CMD_COMMON_TUNE_REQUEST = 246 	# /usr/include/alsa/asoundef.h:209
MIDI_CMD_COMMON_SYSEX_END = 247 	# /usr/include/alsa/asoundef.h:210
MIDI_CMD_COMMON_CLOCK = 248 	# /usr/include/alsa/asoundef.h:211
MIDI_CMD_COMMON_START = 250 	# /usr/include/alsa/asoundef.h:212
MIDI_CMD_COMMON_CONTINUE = 251 	# /usr/include/alsa/asoundef.h:213
MIDI_CMD_COMMON_STOP = 252 	# /usr/include/alsa/asoundef.h:214
MIDI_CMD_COMMON_SENSING = 254 	# /usr/include/alsa/asoundef.h:215
MIDI_CMD_COMMON_RESET = 255 	# /usr/include/alsa/asoundef.h:216
MIDI_CTL_MSB_BANK = 0 	# /usr/include/alsa/asoundef.h:226
MIDI_CTL_MSB_MODWHEEL = 1 	# /usr/include/alsa/asoundef.h:227
MIDI_CTL_MSB_BREATH = 2 	# /usr/include/alsa/asoundef.h:228
MIDI_CTL_MSB_FOOT = 4 	# /usr/include/alsa/asoundef.h:229
MIDI_CTL_MSB_PORTAMENTO_TIME = 5 	# /usr/include/alsa/asoundef.h:230
MIDI_CTL_MSB_DATA_ENTRY = 6 	# /usr/include/alsa/asoundef.h:231
MIDI_CTL_MSB_MAIN_VOLUME = 7 	# /usr/include/alsa/asoundef.h:232
MIDI_CTL_MSB_BALANCE = 8 	# /usr/include/alsa/asoundef.h:233
MIDI_CTL_MSB_PAN = 10 	# /usr/include/alsa/asoundef.h:234
MIDI_CTL_MSB_EXPRESSION = 11 	# /usr/include/alsa/asoundef.h:235
MIDI_CTL_MSB_EFFECT1 = 12 	# /usr/include/alsa/asoundef.h:236
MIDI_CTL_MSB_EFFECT2 = 13 	# /usr/include/alsa/asoundef.h:237
MIDI_CTL_MSB_GENERAL_PURPOSE1 = 16 	# /usr/include/alsa/asoundef.h:238
MIDI_CTL_MSB_GENERAL_PURPOSE2 = 17 	# /usr/include/alsa/asoundef.h:239
MIDI_CTL_MSB_GENERAL_PURPOSE3 = 18 	# /usr/include/alsa/asoundef.h:240
MIDI_CTL_MSB_GENERAL_PURPOSE4 = 19 	# /usr/include/alsa/asoundef.h:241
MIDI_CTL_LSB_BANK = 32 	# /usr/include/alsa/asoundef.h:242
MIDI_CTL_LSB_MODWHEEL = 33 	# /usr/include/alsa/asoundef.h:243
MIDI_CTL_LSB_BREATH = 34 	# /usr/include/alsa/asoundef.h:244
MIDI_CTL_LSB_FOOT = 36 	# /usr/include/alsa/asoundef.h:245
MIDI_CTL_LSB_PORTAMENTO_TIME = 37 	# /usr/include/alsa/asoundef.h:246
MIDI_CTL_LSB_DATA_ENTRY = 38 	# /usr/include/alsa/asoundef.h:247
MIDI_CTL_LSB_MAIN_VOLUME = 39 	# /usr/include/alsa/asoundef.h:248
MIDI_CTL_LSB_BALANCE = 40 	# /usr/include/alsa/asoundef.h:249
MIDI_CTL_LSB_PAN = 42 	# /usr/include/alsa/asoundef.h:250
MIDI_CTL_LSB_EXPRESSION = 43 	# /usr/include/alsa/asoundef.h:251
MIDI_CTL_LSB_EFFECT1 = 44 	# /usr/include/alsa/asoundef.h:252
MIDI_CTL_LSB_EFFECT2 = 45 	# /usr/include/alsa/asoundef.h:253
MIDI_CTL_LSB_GENERAL_PURPOSE1 = 48 	# /usr/include/alsa/asoundef.h:254
MIDI_CTL_LSB_GENERAL_PURPOSE2 = 49 	# /usr/include/alsa/asoundef.h:255
MIDI_CTL_LSB_GENERAL_PURPOSE3 = 50 	# /usr/include/alsa/asoundef.h:256
MIDI_CTL_LSB_GENERAL_PURPOSE4 = 51 	# /usr/include/alsa/asoundef.h:257
MIDI_CTL_SUSTAIN = 64 	# /usr/include/alsa/asoundef.h:258
MIDI_CTL_PORTAMENTO = 65 	# /usr/include/alsa/asoundef.h:259
MIDI_CTL_SOSTENUTO = 66 	# /usr/include/alsa/asoundef.h:260
MIDI_CTL_SUSTENUTO = 66 	# /usr/include/alsa/asoundef.h:261
MIDI_CTL_SOFT_PEDAL = 67 	# /usr/include/alsa/asoundef.h:262
MIDI_CTL_LEGATO_FOOTSWITCH = 68 	# /usr/include/alsa/asoundef.h:263
MIDI_CTL_HOLD2 = 69 	# /usr/include/alsa/asoundef.h:264
MIDI_CTL_SC1_SOUND_VARIATION = 70 	# /usr/include/alsa/asoundef.h:265
MIDI_CTL_SC2_TIMBRE = 71 	# /usr/include/alsa/asoundef.h:266
MIDI_CTL_SC3_RELEASE_TIME = 72 	# /usr/include/alsa/asoundef.h:267
MIDI_CTL_SC4_ATTACK_TIME = 73 	# /usr/include/alsa/asoundef.h:268
MIDI_CTL_SC5_BRIGHTNESS = 74 	# /usr/include/alsa/asoundef.h:269
MIDI_CTL_SC6 = 75 	# /usr/include/alsa/asoundef.h:270
MIDI_CTL_SC7 = 76 	# /usr/include/alsa/asoundef.h:271
MIDI_CTL_SC8 = 77 	# /usr/include/alsa/asoundef.h:272
MIDI_CTL_SC9 = 78 	# /usr/include/alsa/asoundef.h:273
MIDI_CTL_SC10 = 79 	# /usr/include/alsa/asoundef.h:274
MIDI_CTL_GENERAL_PURPOSE5 = 80 	# /usr/include/alsa/asoundef.h:275
MIDI_CTL_GENERAL_PURPOSE6 = 81 	# /usr/include/alsa/asoundef.h:276
MIDI_CTL_GENERAL_PURPOSE7 = 82 	# /usr/include/alsa/asoundef.h:277
MIDI_CTL_GENERAL_PURPOSE8 = 83 	# /usr/include/alsa/asoundef.h:278
MIDI_CTL_PORTAMENTO_CONTROL = 84 	# /usr/include/alsa/asoundef.h:279
MIDI_CTL_E1_REVERB_DEPTH = 91 	# /usr/include/alsa/asoundef.h:280
MIDI_CTL_E2_TREMOLO_DEPTH = 92 	# /usr/include/alsa/asoundef.h:281
MIDI_CTL_E3_CHORUS_DEPTH = 93 	# /usr/include/alsa/asoundef.h:282
MIDI_CTL_E4_DETUNE_DEPTH = 94 	# /usr/include/alsa/asoundef.h:283
MIDI_CTL_E5_PHASER_DEPTH = 95 	# /usr/include/alsa/asoundef.h:284
MIDI_CTL_DATA_INCREMENT = 96 	# /usr/include/alsa/asoundef.h:285
MIDI_CTL_DATA_DECREMENT = 97 	# /usr/include/alsa/asoundef.h:286
MIDI_CTL_NONREG_PARM_NUM_LSB = 98 	# /usr/include/alsa/asoundef.h:287
MIDI_CTL_NONREG_PARM_NUM_MSB = 99 	# /usr/include/alsa/asoundef.h:288
MIDI_CTL_REGIST_PARM_NUM_LSB = 100 	# /usr/include/alsa/asoundef.h:289
MIDI_CTL_REGIST_PARM_NUM_MSB = 101 	# /usr/include/alsa/asoundef.h:290
MIDI_CTL_ALL_SOUNDS_OFF = 120 	# /usr/include/alsa/asoundef.h:291
MIDI_CTL_RESET_CONTROLLERS = 121 	# /usr/include/alsa/asoundef.h:292
MIDI_CTL_LOCAL_CONTROL_SWITCH = 122 	# /usr/include/alsa/asoundef.h:293
MIDI_CTL_ALL_NOTES_OFF = 123 	# /usr/include/alsa/asoundef.h:294
MIDI_CTL_OMNI_OFF = 124 	# /usr/include/alsa/asoundef.h:295
MIDI_CTL_OMNI_ON = 125 	# /usr/include/alsa/asoundef.h:296
MIDI_CTL_MONO1 = 126 	# /usr/include/alsa/asoundef.h:297
MIDI_CTL_MONO2 = 127 	# /usr/include/alsa/asoundef.h:298
IEC958_AES0_PROFESSIONAL = 1 	# /usr/include/alsa/asoundef.h:41
IEC958_AES0_NONAUDIO = 2 	# /usr/include/alsa/asoundef.h:42
IEC958_AES0_PRO_EMPHASIS = 28 	# /usr/include/alsa/asoundef.h:43
IEC958_AES0_PRO_EMPHASIS_NOTID = 0 	# /usr/include/alsa/asoundef.h:44
IEC958_AES0_PRO_EMPHASIS_NONE = 4 	# /usr/include/alsa/asoundef.h:45
IEC958_AES0_PRO_EMPHASIS_5015 = 12 	# /usr/include/alsa/asoundef.h:46
IEC958_AES0_PRO_EMPHASIS_CCITT = 28 	# /usr/include/alsa/asoundef.h:47
IEC958_AES0_PRO_FREQ_UNLOCKED = 32 	# /usr/include/alsa/asoundef.h:48
IEC958_AES0_PRO_FS = 192 	# /usr/include/alsa/asoundef.h:49
IEC958_AES0_PRO_FS_NOTID = 0 	# /usr/include/alsa/asoundef.h:50
IEC958_AES0_PRO_FS_44100 = 64 	# /usr/include/alsa/asoundef.h:51
IEC958_AES0_PRO_FS_48000 = 128 	# /usr/include/alsa/asoundef.h:52
IEC958_AES0_PRO_FS_32000 = 192 	# /usr/include/alsa/asoundef.h:53
IEC958_AES0_CON_NOT_COPYRIGHT = 4 	# /usr/include/alsa/asoundef.h:54
IEC958_AES0_CON_EMPHASIS = 56 	# /usr/include/alsa/asoundef.h:55
IEC958_AES0_CON_EMPHASIS_NONE = 0 	# /usr/include/alsa/asoundef.h:56
IEC958_AES0_CON_EMPHASIS_5015 = 8 	# /usr/include/alsa/asoundef.h:57
IEC958_AES0_CON_MODE = 192 	# /usr/include/alsa/asoundef.h:58
IEC958_AES1_PRO_MODE = 15 	# /usr/include/alsa/asoundef.h:59
IEC958_AES1_PRO_MODE_NOTID = 0 	# /usr/include/alsa/asoundef.h:60
IEC958_AES1_PRO_MODE_STEREOPHONIC = 2 	# /usr/include/alsa/asoundef.h:61
IEC958_AES1_PRO_MODE_SINGLE = 4 	# /usr/include/alsa/asoundef.h:62
IEC958_AES1_PRO_MODE_TWO = 8 	# /usr/include/alsa/asoundef.h:63
IEC958_AES1_PRO_MODE_PRIMARY = 12 	# /usr/include/alsa/asoundef.h:64
IEC958_AES1_PRO_MODE_BYTE3 = 15 	# /usr/include/alsa/asoundef.h:65
IEC958_AES1_PRO_USERBITS = 240 	# /usr/include/alsa/asoundef.h:66
IEC958_AES1_PRO_USERBITS_NOTID = 0 	# /usr/include/alsa/asoundef.h:67
IEC958_AES1_PRO_USERBITS_192 = 128 	# /usr/include/alsa/asoundef.h:68
IEC958_AES1_PRO_USERBITS_UDEF = 192 	# /usr/include/alsa/asoundef.h:69
IEC958_AES1_CON_CATEGORY = 127 	# /usr/include/alsa/asoundef.h:70
IEC958_AES1_CON_GENERAL = 0 	# /usr/include/alsa/asoundef.h:71
IEC958_AES1_CON_EXPERIMENTAL = 64 	# /usr/include/alsa/asoundef.h:72
IEC958_AES1_CON_SOLIDMEM_MASK = 15 	# /usr/include/alsa/asoundef.h:73
IEC958_AES1_CON_SOLIDMEM_ID = 8 	# /usr/include/alsa/asoundef.h:74
IEC958_AES1_CON_BROADCAST1_MASK = 7 	# /usr/include/alsa/asoundef.h:75
IEC958_AES1_CON_BROADCAST1_ID = 4 	# /usr/include/alsa/asoundef.h:76
IEC958_AES1_CON_DIGDIGCONV_MASK = 7 	# /usr/include/alsa/asoundef.h:77
IEC958_AES1_CON_DIGDIGCONV_ID = 2 	# /usr/include/alsa/asoundef.h:78
IEC958_AES1_CON_ADC_COPYRIGHT_MASK = 31 	# /usr/include/alsa/asoundef.h:79
IEC958_AES1_CON_ADC_COPYRIGHT_ID = 6 	# /usr/include/alsa/asoundef.h:80
IEC958_AES1_CON_ADC_MASK = 31 	# /usr/include/alsa/asoundef.h:81
IEC958_AES1_CON_ADC_ID = 22 	# /usr/include/alsa/asoundef.h:82
IEC958_AES1_CON_BROADCAST2_MASK = 15 	# /usr/include/alsa/asoundef.h:83
IEC958_AES1_CON_BROADCAST2_ID = 14 	# /usr/include/alsa/asoundef.h:84
IEC958_AES1_CON_LASEROPT_MASK = 7 	# /usr/include/alsa/asoundef.h:85
IEC958_AES1_CON_LASEROPT_ID = 1 	# /usr/include/alsa/asoundef.h:86
IEC958_AES1_CON_MUSICAL_MASK = 7 	# /usr/include/alsa/asoundef.h:87
IEC958_AES1_CON_MUSICAL_ID = 5 	# /usr/include/alsa/asoundef.h:88
IEC958_AES1_CON_MAGNETIC_MASK = 7 	# /usr/include/alsa/asoundef.h:89
IEC958_AES1_CON_MAGNETIC_ID = 3 	# /usr/include/alsa/asoundef.h:90
IEC958_AES1_CON_IEC908_CD = 1 	# /usr/include/alsa/asoundef.h:91
IEC958_AES1_CON_NON_IEC908_CD = 9 	# /usr/include/alsa/asoundef.h:92
IEC958_AES1_CON_PCM_CODER = 2 	# /usr/include/alsa/asoundef.h:93
IEC958_AES1_CON_SAMPLER = 34 	# /usr/include/alsa/asoundef.h:94
IEC958_AES1_CON_MIXER = 18 	# /usr/include/alsa/asoundef.h:95
IEC958_AES1_CON_RATE_CONVERTER = 26 	# /usr/include/alsa/asoundef.h:96
IEC958_AES1_CON_SYNTHESIZER = 5 	# /usr/include/alsa/asoundef.h:97
IEC958_AES1_CON_MICROPHONE = 13 	# /usr/include/alsa/asoundef.h:98
IEC958_AES1_CON_DAT = 3 	# /usr/include/alsa/asoundef.h:99
IEC958_AES1_CON_VCR = 11 	# /usr/include/alsa/asoundef.h:100
IEC958_AES1_CON_ORIGINAL = 128 	# /usr/include/alsa/asoundef.h:101
IEC958_AES2_PRO_SBITS = 7 	# /usr/include/alsa/asoundef.h:102
IEC958_AES2_PRO_SBITS_20 = 2 	# /usr/include/alsa/asoundef.h:103
IEC958_AES2_PRO_SBITS_24 = 4 	# /usr/include/alsa/asoundef.h:104
IEC958_AES2_PRO_SBITS_UDEF = 6 	# /usr/include/alsa/asoundef.h:105
IEC958_AES2_PRO_WORDLEN = 56 	# /usr/include/alsa/asoundef.h:106
IEC958_AES2_PRO_WORDLEN_NOTID = 0 	# /usr/include/alsa/asoundef.h:107
IEC958_AES2_PRO_WORDLEN_22_18 = 16 	# /usr/include/alsa/asoundef.h:108
IEC958_AES2_PRO_WORDLEN_23_19 = 32 	# /usr/include/alsa/asoundef.h:109
IEC958_AES2_PRO_WORDLEN_24_20 = 40 	# /usr/include/alsa/asoundef.h:110
IEC958_AES2_PRO_WORDLEN_20_16 = 48 	# /usr/include/alsa/asoundef.h:111
IEC958_AES2_CON_SOURCE = 15 	# /usr/include/alsa/asoundef.h:112
IEC958_AES2_CON_SOURCE_UNSPEC = 0 	# /usr/include/alsa/asoundef.h:113
IEC958_AES2_CON_CHANNEL = 240 	# /usr/include/alsa/asoundef.h:114
IEC958_AES2_CON_CHANNEL_UNSPEC = 0 	# /usr/include/alsa/asoundef.h:115
IEC958_AES3_CON_FS = 15 	# /usr/include/alsa/asoundef.h:116
IEC958_AES3_CON_FS_44100 = 0 	# /usr/include/alsa/asoundef.h:117
IEC958_AES3_CON_FS_48000 = 2 	# /usr/include/alsa/asoundef.h:118
IEC958_AES3_CON_FS_32000 = 3 	# /usr/include/alsa/asoundef.h:119
IEC958_AES3_CON_CLOCK = 48 	# /usr/include/alsa/asoundef.h:120
IEC958_AES3_CON_CLOCK_1000PPM = 0 	# /usr/include/alsa/asoundef.h:121
IEC958_AES3_CON_CLOCK_50PPM = 16 	# /usr/include/alsa/asoundef.h:122
IEC958_AES3_CON_CLOCK_VARIABLE = 32 	# /usr/include/alsa/asoundef.h:123
MIDI_CHANNELS = 16 	# /usr/include/alsa/asoundef.h:133
MIDI_GM_DRUM_CHANNEL = 9 	# /usr/include/alsa/asoundef.h:134
MIDI_CMD_NOTE_OFF = 128 	# /usr/include/alsa/asoundef.h:142
MIDI_CMD_NOTE_ON = 144 	# /usr/include/alsa/asoundef.h:143
MIDI_CMD_NOTE_PRESSURE = 160 	# /usr/include/alsa/asoundef.h:144
MIDI_CMD_CONTROL = 176 	# /usr/include/alsa/asoundef.h:145
MIDI_CMD_PGM_CHANGE = 192 	# /usr/include/alsa/asoundef.h:146
MIDI_CMD_CHANNEL_PRESSURE = 208 	# /usr/include/alsa/asoundef.h:147
MIDI_CMD_BENDER = 224 	# /usr/include/alsa/asoundef.h:148
MIDI_CMD_COMMON_SYSEX = 240 	# /usr/include/alsa/asoundef.h:150
MIDI_CMD_COMMON_MTC_QUARTER = 241 	# /usr/include/alsa/asoundef.h:151
MIDI_CMD_COMMON_SONG_POS = 242 	# /usr/include/alsa/asoundef.h:152
MIDI_CMD_COMMON_SONG_SELECT = 243 	# /usr/include/alsa/asoundef.h:153
MIDI_CMD_COMMON_TUNE_REQUEST = 246 	# /usr/include/alsa/asoundef.h:154
MIDI_CMD_COMMON_SYSEX_END = 247 	# /usr/include/alsa/asoundef.h:155
MIDI_CMD_COMMON_CLOCK = 248 	# /usr/include/alsa/asoundef.h:156
MIDI_CMD_COMMON_START = 250 	# /usr/include/alsa/asoundef.h:157
MIDI_CMD_COMMON_CONTINUE = 251 	# /usr/include/alsa/asoundef.h:158
MIDI_CMD_COMMON_STOP = 252 	# /usr/include/alsa/asoundef.h:159
MIDI_CMD_COMMON_SENSING = 254 	# /usr/include/alsa/asoundef.h:160
MIDI_CMD_COMMON_RESET = 255 	# /usr/include/alsa/asoundef.h:161
MIDI_CTL_MSB_BANK = 0 	# /usr/include/alsa/asoundef.h:171
MIDI_CTL_MSB_MODWHEEL = 1 	# /usr/include/alsa/asoundef.h:172
MIDI_CTL_MSB_BREATH = 2 	# /usr/include/alsa/asoundef.h:173
MIDI_CTL_MSB_FOOT = 4 	# /usr/include/alsa/asoundef.h:174
MIDI_CTL_MSB_PORTAMENTO_TIME = 5 	# /usr/include/alsa/asoundef.h:175
MIDI_CTL_MSB_DATA_ENTRY = 6 	# /usr/include/alsa/asoundef.h:176
MIDI_CTL_MSB_MAIN_VOLUME = 7 	# /usr/include/alsa/asoundef.h:177
MIDI_CTL_MSB_BALANCE = 8 	# /usr/include/alsa/asoundef.h:178
MIDI_CTL_MSB_PAN = 10 	# /usr/include/alsa/asoundef.h:179
MIDI_CTL_MSB_EXPRESSION = 11 	# /usr/include/alsa/asoundef.h:180
MIDI_CTL_MSB_EFFECT1 = 12 	# /usr/include/alsa/asoundef.h:181
MIDI_CTL_MSB_EFFECT2 = 13 	# /usr/include/alsa/asoundef.h:182
MIDI_CTL_MSB_GENERAL_PURPOSE1 = 16 	# /usr/include/alsa/asoundef.h:183
MIDI_CTL_MSB_GENERAL_PURPOSE2 = 17 	# /usr/include/alsa/asoundef.h:184
MIDI_CTL_MSB_GENERAL_PURPOSE3 = 18 	# /usr/include/alsa/asoundef.h:185
MIDI_CTL_MSB_GENERAL_PURPOSE4 = 19 	# /usr/include/alsa/asoundef.h:186
MIDI_CTL_LSB_BANK = 32 	# /usr/include/alsa/asoundef.h:187
MIDI_CTL_LSB_MODWHEEL = 33 	# /usr/include/alsa/asoundef.h:188
MIDI_CTL_LSB_BREATH = 34 	# /usr/include/alsa/asoundef.h:189
MIDI_CTL_LSB_FOOT = 36 	# /usr/include/alsa/asoundef.h:190
MIDI_CTL_LSB_PORTAMENTO_TIME = 37 	# /usr/include/alsa/asoundef.h:191
MIDI_CTL_LSB_DATA_ENTRY = 38 	# /usr/include/alsa/asoundef.h:192
MIDI_CTL_LSB_MAIN_VOLUME = 39 	# /usr/include/alsa/asoundef.h:193
MIDI_CTL_LSB_BALANCE = 40 	# /usr/include/alsa/asoundef.h:194
MIDI_CTL_LSB_PAN = 42 	# /usr/include/alsa/asoundef.h:195
MIDI_CTL_LSB_EXPRESSION = 43 	# /usr/include/alsa/asoundef.h:196
MIDI_CTL_LSB_EFFECT1 = 44 	# /usr/include/alsa/asoundef.h:197
MIDI_CTL_LSB_EFFECT2 = 45 	# /usr/include/alsa/asoundef.h:198
MIDI_CTL_LSB_GENERAL_PURPOSE1 = 48 	# /usr/include/alsa/asoundef.h:199
MIDI_CTL_LSB_GENERAL_PURPOSE2 = 49 	# /usr/include/alsa/asoundef.h:200
MIDI_CTL_LSB_GENERAL_PURPOSE3 = 50 	# /usr/include/alsa/asoundef.h:201
MIDI_CTL_LSB_GENERAL_PURPOSE4 = 51 	# /usr/include/alsa/asoundef.h:202
MIDI_CTL_SUSTAIN = 64 	# /usr/include/alsa/asoundef.h:203
MIDI_CTL_PORTAMENTO = 65 	# /usr/include/alsa/asoundef.h:204
MIDI_CTL_SOSTENUTO = 66 	# /usr/include/alsa/asoundef.h:205
MIDI_CTL_SUSTENUTO = 66 	# /usr/include/alsa/asoundef.h:206
MIDI_CTL_SOFT_PEDAL = 67 	# /usr/include/alsa/asoundef.h:207
MIDI_CTL_LEGATO_FOOTSWITCH = 68 	# /usr/include/alsa/asoundef.h:208
MIDI_CTL_HOLD2 = 69 	# /usr/include/alsa/asoundef.h:209
MIDI_CTL_SC1_SOUND_VARIATION = 70 	# /usr/include/alsa/asoundef.h:210
MIDI_CTL_SC2_TIMBRE = 71 	# /usr/include/alsa/asoundef.h:211
MIDI_CTL_SC3_RELEASE_TIME = 72 	# /usr/include/alsa/asoundef.h:212
MIDI_CTL_SC4_ATTACK_TIME = 73 	# /usr/include/alsa/asoundef.h:213
MIDI_CTL_SC5_BRIGHTNESS = 74 	# /usr/include/alsa/asoundef.h:214
MIDI_CTL_SC6 = 75 	# /usr/include/alsa/asoundef.h:215
MIDI_CTL_SC7 = 76 	# /usr/include/alsa/asoundef.h:216
MIDI_CTL_SC8 = 77 	# /usr/include/alsa/asoundef.h:217
MIDI_CTL_SC9 = 78 	# /usr/include/alsa/asoundef.h:218
MIDI_CTL_SC10 = 79 	# /usr/include/alsa/asoundef.h:219
MIDI_CTL_GENERAL_PURPOSE5 = 80 	# /usr/include/alsa/asoundef.h:220
MIDI_CTL_GENERAL_PURPOSE6 = 81 	# /usr/include/alsa/asoundef.h:221
MIDI_CTL_GENERAL_PURPOSE7 = 82 	# /usr/include/alsa/asoundef.h:222
MIDI_CTL_GENERAL_PURPOSE8 = 83 	# /usr/include/alsa/asoundef.h:223
MIDI_CTL_PORTAMENTO_CONTROL = 84 	# /usr/include/alsa/asoundef.h:224
MIDI_CTL_E1_REVERB_DEPTH = 91 	# /usr/include/alsa/asoundef.h:225
MIDI_CTL_E2_TREMOLO_DEPTH = 92 	# /usr/include/alsa/asoundef.h:226
MIDI_CTL_E3_CHORUS_DEPTH = 93 	# /usr/include/alsa/asoundef.h:227
MIDI_CTL_E4_DETUNE_DEPTH = 94 	# /usr/include/alsa/asoundef.h:228
MIDI_CTL_E5_PHASER_DEPTH = 95 	# /usr/include/alsa/asoundef.h:229
MIDI_CTL_DATA_INCREMENT = 96 	# /usr/include/alsa/asoundef.h:230
MIDI_CTL_DATA_DECREMENT = 97 	# /usr/include/alsa/asoundef.h:231
MIDI_CTL_NONREG_PARM_NUM_LSB = 98 	# /usr/include/alsa/asoundef.h:232
MIDI_CTL_NONREG_PARM_NUM_MSB = 99 	# /usr/include/alsa/asoundef.h:233
MIDI_CTL_REGIST_PARM_NUM_LSB = 100 	# /usr/include/alsa/asoundef.h:234
MIDI_CTL_REGIST_PARM_NUM_MSB = 101 	# /usr/include/alsa/asoundef.h:235
MIDI_CTL_ALL_SOUNDS_OFF = 120 	# /usr/include/alsa/asoundef.h:236
MIDI_CTL_RESET_CONTROLLERS = 121 	# /usr/include/alsa/asoundef.h:237
MIDI_CTL_LOCAL_CONTROL_SWITCH = 122 	# /usr/include/alsa/asoundef.h:238
MIDI_CTL_ALL_NOTES_OFF = 123 	# /usr/include/alsa/asoundef.h:239
MIDI_CTL_OMNI_OFF = 124 	# /usr/include/alsa/asoundef.h:240
MIDI_CTL_OMNI_ON = 125 	# /usr/include/alsa/asoundef.h:241
MIDI_CTL_MONO1 = 126 	# /usr/include/alsa/asoundef.h:242
MIDI_CTL_MONO2 = 127 	# /usr/include/alsa/asoundef.h:243
SND_LIB_MAJOR = 1 	# /usr/include/alsa/version.h:5
SND_LIB_MINOR = 0 	# /usr/include/alsa/version.h:6
SND_LIB_SUBMINOR = 14 	# /usr/include/alsa/version.h:7
SND_LIB_EXTRAVER = 100003 	# /usr/include/alsa/version.h:8
SND_LIB_VERSION = 65550 	# /usr/include/alsa/version.h:10
# /usr/include/alsa/global.h:47
snd_asoundlib_version = _lib.snd_asoundlib_version
snd_asoundlib_version.restype = c_char_p
snd_asoundlib_version.argtypes = []

# /usr/include/alsa/global.h:100
snd_dlopen = _lib.snd_dlopen
snd_dlopen.restype = POINTER(c_void)
snd_dlopen.argtypes = [c_char_p, c_int]

# /usr/include/alsa/global.h:101
snd_dlsym = _lib.snd_dlsym
snd_dlsym.restype = POINTER(c_void)
snd_dlsym.argtypes = [POINTER(None), c_char_p, c_char_p]

# /usr/include/alsa/global.h:102
snd_dlclose = _lib.snd_dlclose
snd_dlclose.restype = c_int
snd_dlclose.argtypes = [POINTER(None)]

class struct__snd_async_handler(Structure):
    __slots__ = [
    ]
struct__snd_async_handler._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_async_handler(Structure):
    __slots__ = [
    ]
struct__snd_async_handler._fields_ = [
    ('_opaque_struct', c_int)
]

snd_async_handler_t = struct__snd_async_handler 	# /usr/include/alsa/global.h:111
snd_async_callback_t = CFUNCTYPE(None, POINTER(snd_async_handler_t)) 	# /usr/include/alsa/global.h:118
# /usr/include/alsa/global.h:120
snd_async_add_handler = _lib.snd_async_add_handler
snd_async_add_handler.restype = c_int
snd_async_add_handler.argtypes = [POINTER(POINTER(snd_async_handler_t)), c_int, snd_async_callback_t, POINTER(None)]

# /usr/include/alsa/global.h:122
snd_async_del_handler = _lib.snd_async_del_handler
snd_async_del_handler.restype = c_int
snd_async_del_handler.argtypes = [POINTER(snd_async_handler_t)]

# /usr/include/alsa/global.h:123
snd_async_handler_get_fd = _lib.snd_async_handler_get_fd
snd_async_handler_get_fd.restype = c_int
snd_async_handler_get_fd.argtypes = [POINTER(snd_async_handler_t)]

# /usr/include/alsa/global.h:124
snd_async_handler_get_signo = _lib.snd_async_handler_get_signo
snd_async_handler_get_signo.restype = c_int
snd_async_handler_get_signo.argtypes = [POINTER(snd_async_handler_t)]

# /usr/include/alsa/global.h:125
snd_async_handler_get_callback_private = _lib.snd_async_handler_get_callback_private
snd_async_handler_get_callback_private.restype = POINTER(c_void)
snd_async_handler_get_callback_private.argtypes = [POINTER(snd_async_handler_t)]

class struct_snd_shm_area(Structure):
    __slots__ = [
    ]
struct_snd_shm_area._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/global.h:127
snd_shm_area_create = _lib.snd_shm_area_create
snd_shm_area_create.restype = POINTER(struct_snd_shm_area)
snd_shm_area_create.argtypes = [c_int, POINTER(None)]

class struct_snd_shm_area(Structure):
    __slots__ = [
    ]
struct_snd_shm_area._fields_ = [
    ('_opaque_struct', c_int)
]

class struct_snd_shm_area(Structure):
    __slots__ = [
    ]
struct_snd_shm_area._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/global.h:128
snd_shm_area_share = _lib.snd_shm_area_share
snd_shm_area_share.restype = POINTER(struct_snd_shm_area)
snd_shm_area_share.argtypes = [POINTER(struct_snd_shm_area)]

class struct_snd_shm_area(Structure):
    __slots__ = [
    ]
struct_snd_shm_area._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/global.h:129
snd_shm_area_destroy = _lib.snd_shm_area_destroy
snd_shm_area_destroy.restype = c_int
snd_shm_area_destroy.argtypes = [POINTER(struct_snd_shm_area)]

# /usr/include/alsa/global.h:131
snd_user_file = _lib.snd_user_file
snd_user_file.restype = c_int
snd_user_file.argtypes = [c_char_p, POINTER(c_char_p)]


# XXX from `man gettimeofday`
class struct_timeval(Structure):
    _fields_ = [
        ('tv_sec', c_long),
        ('tv_usec', c_long)
    ]

snd_timestamp_t = struct_timeval 	# /usr/include/alsa/global.h:146


# XXX wrong, but not used ATM
class struct_timespec(Structure):
    __slots__ = [
    ]
struct_timespec._fields_ = [
    ('_opaque_struct', c_int)
]

class struct_timespec(Structure):
    __slots__ = [
    ]
struct_timespec._fields_ = [
    ('_opaque_struct', c_int)
]

snd_htimestamp_t = struct_timespec 	# /usr/include/alsa/global.h:148
class struct__snd_input(Structure):
    __slots__ = [
    ]
struct__snd_input._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_input(Structure):
    __slots__ = [
    ]
struct__snd_input._fields_ = [
    ('_opaque_struct', c_int)
]

snd_input_t = struct__snd_input 	# /usr/include/alsa/input.h:54
enum__snd_input_type = c_int
SND_INPUT_STDIO = 1
SND_INPUT_BUFFER = 2
snd_input_type_t = enum__snd_input_type 	# /usr/include/alsa/input.h:62
# /usr/include/alsa/input.h:64
snd_input_stdio_open = _lib.snd_input_stdio_open
snd_input_stdio_open.restype = c_int
snd_input_stdio_open.argtypes = [POINTER(POINTER(snd_input_t)), c_char_p, c_char_p]

class struct__IO_FILE(Structure):
    __slots__ = [
    ]
struct__IO_FILE._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__IO_FILE(Structure):
    __slots__ = [
    ]
struct__IO_FILE._fields_ = [
    ('_opaque_struct', c_int)
]

FILE = struct__IO_FILE 	# /usr/include/gentoo-multilib/amd64/stdio.h:46
# /usr/include/alsa/input.h:65
snd_input_stdio_attach = _lib.snd_input_stdio_attach
snd_input_stdio_attach.restype = c_int
snd_input_stdio_attach.argtypes = [POINTER(POINTER(snd_input_t)), POINTER(FILE), c_int]

__ssize_t = c_long 	# /usr/include/gentoo-multilib/amd64/bits/types.h:183
ssize_t = __ssize_t 	# /usr/include/gentoo-multilib/amd64/unistd.h:189
# /usr/include/alsa/input.h:66
snd_input_buffer_open = _lib.snd_input_buffer_open
snd_input_buffer_open.restype = c_int
snd_input_buffer_open.argtypes = [POINTER(POINTER(snd_input_t)), c_char_p, ssize_t]

# /usr/include/alsa/input.h:67
snd_input_close = _lib.snd_input_close
snd_input_close.restype = c_int
snd_input_close.argtypes = [POINTER(snd_input_t)]

# /usr/include/alsa/input.h:68
snd_input_scanf = _lib.snd_input_scanf
snd_input_scanf.restype = c_int
snd_input_scanf.argtypes = [POINTER(snd_input_t), c_char_p]

# /usr/include/alsa/input.h:73
snd_input_gets = _lib.snd_input_gets
snd_input_gets.restype = c_char_p
snd_input_gets.argtypes = [POINTER(snd_input_t), c_char_p, c_size_t]

# /usr/include/alsa/input.h:74
snd_input_getc = _lib.snd_input_getc
snd_input_getc.restype = c_int
snd_input_getc.argtypes = [POINTER(snd_input_t)]

# /usr/include/alsa/input.h:75
snd_input_ungetc = _lib.snd_input_ungetc
snd_input_ungetc.restype = c_int
snd_input_ungetc.argtypes = [POINTER(snd_input_t), c_int]

SND_ERROR_BEGIN = 500000 	# /usr/include/alsa/error.h:41
SND_ERROR_INCOMPATIBLE_VERSION = 500000 	# /usr/include/alsa/error.h:42
SND_ERROR_ALISP_NIL = 500001 	# /usr/include/alsa/error.h:43
# /usr/include/alsa/error.h:45
snd_strerror = _lib.snd_strerror
snd_strerror.restype = c_char_p
snd_strerror.argtypes = [c_int]

snd_lib_error_handler_t = CFUNCTYPE(None, c_char_p, c_int, c_char_p, c_int, c_char_p) 	# /usr/include/alsa/error.h:59
# /usr/include/alsa/error.h:61
snd_lib_error_set_handler = _lib.snd_lib_error_set_handler
snd_lib_error_set_handler.restype = c_int
snd_lib_error_set_handler.argtypes = [snd_lib_error_handler_t]

SND_CONFIG_DLSYM_VERSION_EVALUATE = 0 	# /usr/include/alsa/conf.h:43
SND_CONFIG_DLSYM_VERSION_HOOK = 0 	# /usr/include/alsa/conf.h:45
enum__snd_config_type = c_int
SND_CONFIG_TYPE_INTEGER = 1
SND_CONFIG_TYPE_INTEGER64 = 2
SND_CONFIG_TYPE_REAL = 3
SND_CONFIG_TYPE_STRING = 4
SND_CONFIG_TYPE_POINTER = 5
SND_CONFIG_TYPE_COMPOUND = 1024
snd_config_type_t = enum__snd_config_type 	# /usr/include/alsa/conf.h:61
class struct__snd_config(Structure):
    __slots__ = [
    ]
struct__snd_config._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_config(Structure):
    __slots__ = [
    ]
struct__snd_config._fields_ = [
    ('_opaque_struct', c_int)
]

snd_config_t = struct__snd_config 	# /usr/include/alsa/conf.h:69
class struct__snd_config_iterator(Structure):
    __slots__ = [
    ]
struct__snd_config_iterator._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_config_iterator(Structure):
    __slots__ = [
    ]
struct__snd_config_iterator._fields_ = [
    ('_opaque_struct', c_int)
]

snd_config_iterator_t = POINTER(struct__snd_config_iterator) 	# /usr/include/alsa/conf.h:77
class struct__snd_config_update(Structure):
    __slots__ = [
    ]
struct__snd_config_update._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_config_update(Structure):
    __slots__ = [
    ]
struct__snd_config_update._fields_ = [
    ('_opaque_struct', c_int)
]

snd_config_update_t = struct__snd_config_update 	# /usr/include/alsa/conf.h:83
# /usr/include/alsa/conf.h:87
snd_config_top = _lib.snd_config_top
snd_config_top.restype = c_int
snd_config_top.argtypes = [POINTER(POINTER(snd_config_t))]

# /usr/include/alsa/conf.h:89
snd_config_load = _lib.snd_config_load
snd_config_load.restype = c_int
snd_config_load.argtypes = [POINTER(snd_config_t), POINTER(snd_input_t)]

# /usr/include/alsa/conf.h:90
snd_config_load_override = _lib.snd_config_load_override
snd_config_load_override.restype = c_int
snd_config_load_override.argtypes = [POINTER(snd_config_t), POINTER(snd_input_t)]

class struct__snd_output(Structure):
    __slots__ = [
    ]
struct__snd_output._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_output(Structure):
    __slots__ = [
    ]
struct__snd_output._fields_ = [
    ('_opaque_struct', c_int)
]

snd_output_t = struct__snd_output 	# /usr/include/alsa/output.h:54

# XXX output.h was not generated because of varargs, but we need this..
snd_output_stdio_open = _lib.snd_output_stdio_open
snd_output_stdio_open.restype = c_int
snd_output_stdio_open.argtypes = [POINTER(POINTER(snd_output_t)), c_char_p, c_char_p]

# XXX no args for varargs function (python can do formatting)
snd_output_printf = _lib.snd_output_printf
snd_output_printf.restype = c_int
snd_output_printf.argtypes = [POINTER(snd_output_t), c_char_p]

# /usr/include/alsa/conf.h:91
snd_config_save = _lib.snd_config_save
snd_config_save.restype = c_int
snd_config_save.argtypes = [POINTER(snd_config_t), POINTER(snd_output_t)]

# /usr/include/alsa/conf.h:92
snd_config_update = _lib.snd_config_update
snd_config_update.restype = c_int
snd_config_update.argtypes = []

# /usr/include/alsa/conf.h:93
snd_config_update_r = _lib.snd_config_update_r
snd_config_update_r.restype = c_int
snd_config_update_r.argtypes = [POINTER(POINTER(snd_config_t)), POINTER(POINTER(snd_config_update_t)), c_char_p]

# /usr/include/alsa/conf.h:94
snd_config_update_free = _lib.snd_config_update_free
snd_config_update_free.restype = c_int
snd_config_update_free.argtypes = [POINTER(snd_config_update_t)]

# /usr/include/alsa/conf.h:95
snd_config_update_free_global = _lib.snd_config_update_free_global
snd_config_update_free_global.restype = c_int
snd_config_update_free_global.argtypes = []

# /usr/include/alsa/conf.h:97
snd_config_search = _lib.snd_config_search
snd_config_search.restype = c_int
snd_config_search.argtypes = [POINTER(snd_config_t), c_char_p, POINTER(POINTER(snd_config_t))]

# /usr/include/alsa/conf.h:99
snd_config_searchv = _lib.snd_config_searchv
snd_config_searchv.restype = c_int
snd_config_searchv.argtypes = [POINTER(snd_config_t), POINTER(POINTER(snd_config_t))]

# /usr/include/alsa/conf.h:101
snd_config_search_definition = _lib.snd_config_search_definition
snd_config_search_definition.restype = c_int
snd_config_search_definition.argtypes = [POINTER(snd_config_t), c_char_p, c_char_p, POINTER(POINTER(snd_config_t))]

# /usr/include/alsa/conf.h:105
snd_config_expand = _lib.snd_config_expand
snd_config_expand.restype = c_int
snd_config_expand.argtypes = [POINTER(snd_config_t), POINTER(snd_config_t), c_char_p, POINTER(snd_config_t), POINTER(POINTER(snd_config_t))]

# /usr/include/alsa/conf.h:108
snd_config_evaluate = _lib.snd_config_evaluate
snd_config_evaluate.restype = c_int
snd_config_evaluate.argtypes = [POINTER(snd_config_t), POINTER(snd_config_t), POINTER(snd_config_t), POINTER(POINTER(snd_config_t))]

# /usr/include/alsa/conf.h:111
snd_config_add = _lib.snd_config_add
snd_config_add.restype = c_int
snd_config_add.argtypes = [POINTER(snd_config_t), POINTER(snd_config_t)]

# /usr/include/alsa/conf.h:112
snd_config_delete = _lib.snd_config_delete
snd_config_delete.restype = c_int
snd_config_delete.argtypes = [POINTER(snd_config_t)]

# /usr/include/alsa/conf.h:113
snd_config_delete_compound_members = _lib.snd_config_delete_compound_members
snd_config_delete_compound_members.restype = c_int
snd_config_delete_compound_members.argtypes = [POINTER(snd_config_t)]

# /usr/include/alsa/conf.h:114
snd_config_copy = _lib.snd_config_copy
snd_config_copy.restype = c_int
snd_config_copy.argtypes = [POINTER(POINTER(snd_config_t)), POINTER(snd_config_t)]

# /usr/include/alsa/conf.h:116
snd_config_make = _lib.snd_config_make
snd_config_make.restype = c_int
snd_config_make.argtypes = [POINTER(POINTER(snd_config_t)), c_char_p, snd_config_type_t]

# /usr/include/alsa/conf.h:118
snd_config_make_integer = _lib.snd_config_make_integer
snd_config_make_integer.restype = c_int
snd_config_make_integer.argtypes = [POINTER(POINTER(snd_config_t)), c_char_p]

# /usr/include/alsa/conf.h:119
snd_config_make_integer64 = _lib.snd_config_make_integer64
snd_config_make_integer64.restype = c_int
snd_config_make_integer64.argtypes = [POINTER(POINTER(snd_config_t)), c_char_p]

# /usr/include/alsa/conf.h:120
snd_config_make_real = _lib.snd_config_make_real
snd_config_make_real.restype = c_int
snd_config_make_real.argtypes = [POINTER(POINTER(snd_config_t)), c_char_p]

# /usr/include/alsa/conf.h:121
snd_config_make_string = _lib.snd_config_make_string
snd_config_make_string.restype = c_int
snd_config_make_string.argtypes = [POINTER(POINTER(snd_config_t)), c_char_p]

# /usr/include/alsa/conf.h:122
snd_config_make_pointer = _lib.snd_config_make_pointer
snd_config_make_pointer.restype = c_int
snd_config_make_pointer.argtypes = [POINTER(POINTER(snd_config_t)), c_char_p]

# /usr/include/alsa/conf.h:123
snd_config_make_compound = _lib.snd_config_make_compound
snd_config_make_compound.restype = c_int
snd_config_make_compound.argtypes = [POINTER(POINTER(snd_config_t)), c_char_p, c_int]

# /usr/include/alsa/conf.h:125
snd_config_imake_integer = _lib.snd_config_imake_integer
snd_config_imake_integer.restype = c_int
snd_config_imake_integer.argtypes = [POINTER(POINTER(snd_config_t)), c_char_p, c_long]

# /usr/include/alsa/conf.h:126
snd_config_imake_integer64 = _lib.snd_config_imake_integer64
snd_config_imake_integer64.restype = c_int
snd_config_imake_integer64.argtypes = [POINTER(POINTER(snd_config_t)), c_char_p, c_longlong]

# /usr/include/alsa/conf.h:127
snd_config_imake_real = _lib.snd_config_imake_real
snd_config_imake_real.restype = c_int
snd_config_imake_real.argtypes = [POINTER(POINTER(snd_config_t)), c_char_p, c_double]

# /usr/include/alsa/conf.h:128
snd_config_imake_string = _lib.snd_config_imake_string
snd_config_imake_string.restype = c_int
snd_config_imake_string.argtypes = [POINTER(POINTER(snd_config_t)), c_char_p, c_char_p]

# /usr/include/alsa/conf.h:129
snd_config_imake_pointer = _lib.snd_config_imake_pointer
snd_config_imake_pointer.restype = c_int
snd_config_imake_pointer.argtypes = [POINTER(POINTER(snd_config_t)), c_char_p, POINTER(None)]

# /usr/include/alsa/conf.h:131
snd_config_get_type = _lib.snd_config_get_type
snd_config_get_type.restype = snd_config_type_t
snd_config_get_type.argtypes = [POINTER(snd_config_t)]

# /usr/include/alsa/conf.h:133
snd_config_set_id = _lib.snd_config_set_id
snd_config_set_id.restype = c_int
snd_config_set_id.argtypes = [POINTER(snd_config_t), c_char_p]

# /usr/include/alsa/conf.h:134
snd_config_set_integer = _lib.snd_config_set_integer
snd_config_set_integer.restype = c_int
snd_config_set_integer.argtypes = [POINTER(snd_config_t), c_long]

# /usr/include/alsa/conf.h:135
snd_config_set_integer64 = _lib.snd_config_set_integer64
snd_config_set_integer64.restype = c_int
snd_config_set_integer64.argtypes = [POINTER(snd_config_t), c_longlong]

# /usr/include/alsa/conf.h:136
snd_config_set_real = _lib.snd_config_set_real
snd_config_set_real.restype = c_int
snd_config_set_real.argtypes = [POINTER(snd_config_t), c_double]

# /usr/include/alsa/conf.h:137
snd_config_set_string = _lib.snd_config_set_string
snd_config_set_string.restype = c_int
snd_config_set_string.argtypes = [POINTER(snd_config_t), c_char_p]

# /usr/include/alsa/conf.h:138
snd_config_set_ascii = _lib.snd_config_set_ascii
snd_config_set_ascii.restype = c_int
snd_config_set_ascii.argtypes = [POINTER(snd_config_t), c_char_p]

# /usr/include/alsa/conf.h:139
snd_config_set_pointer = _lib.snd_config_set_pointer
snd_config_set_pointer.restype = c_int
snd_config_set_pointer.argtypes = [POINTER(snd_config_t), POINTER(None)]

# /usr/include/alsa/conf.h:140
snd_config_get_id = _lib.snd_config_get_id
snd_config_get_id.restype = c_int
snd_config_get_id.argtypes = [POINTER(snd_config_t), POINTER(c_char_p)]

# /usr/include/alsa/conf.h:141
snd_config_get_integer = _lib.snd_config_get_integer
snd_config_get_integer.restype = c_int
snd_config_get_integer.argtypes = [POINTER(snd_config_t), POINTER(c_long)]

# /usr/include/alsa/conf.h:142
snd_config_get_integer64 = _lib.snd_config_get_integer64
snd_config_get_integer64.restype = c_int
snd_config_get_integer64.argtypes = [POINTER(snd_config_t), POINTER(c_longlong)]

# /usr/include/alsa/conf.h:143
snd_config_get_real = _lib.snd_config_get_real
snd_config_get_real.restype = c_int
snd_config_get_real.argtypes = [POINTER(snd_config_t), POINTER(c_double)]

# /usr/include/alsa/conf.h:144
snd_config_get_ireal = _lib.snd_config_get_ireal
snd_config_get_ireal.restype = c_int
snd_config_get_ireal.argtypes = [POINTER(snd_config_t), POINTER(c_double)]

# /usr/include/alsa/conf.h:145
snd_config_get_string = _lib.snd_config_get_string
snd_config_get_string.restype = c_int
snd_config_get_string.argtypes = [POINTER(snd_config_t), POINTER(c_char_p)]

# /usr/include/alsa/conf.h:146
snd_config_get_ascii = _lib.snd_config_get_ascii
snd_config_get_ascii.restype = c_int
snd_config_get_ascii.argtypes = [POINTER(snd_config_t), POINTER(c_char_p)]

# /usr/include/alsa/conf.h:147
snd_config_get_pointer = _lib.snd_config_get_pointer
snd_config_get_pointer.restype = c_int
snd_config_get_pointer.argtypes = [POINTER(snd_config_t), POINTER(POINTER(None))]

# /usr/include/alsa/conf.h:148
snd_config_test_id = _lib.snd_config_test_id
snd_config_test_id.restype = c_int
snd_config_test_id.argtypes = [POINTER(snd_config_t), c_char_p]

# /usr/include/alsa/conf.h:150
snd_config_iterator_first = _lib.snd_config_iterator_first
snd_config_iterator_first.restype = snd_config_iterator_t
snd_config_iterator_first.argtypes = [POINTER(snd_config_t)]

# /usr/include/alsa/conf.h:151
snd_config_iterator_next = _lib.snd_config_iterator_next
snd_config_iterator_next.restype = snd_config_iterator_t
snd_config_iterator_next.argtypes = [snd_config_iterator_t]

# /usr/include/alsa/conf.h:152
snd_config_iterator_end = _lib.snd_config_iterator_end
snd_config_iterator_end.restype = snd_config_iterator_t
snd_config_iterator_end.argtypes = [POINTER(snd_config_t)]

# /usr/include/alsa/conf.h:153
snd_config_iterator_entry = _lib.snd_config_iterator_entry
snd_config_iterator_entry.restype = POINTER(snd_config_t)
snd_config_iterator_entry.argtypes = [snd_config_iterator_t]

# /usr/include/alsa/conf.h:168
snd_config_get_bool_ascii = _lib.snd_config_get_bool_ascii
snd_config_get_bool_ascii.restype = c_int
snd_config_get_bool_ascii.argtypes = [c_char_p]

# /usr/include/alsa/conf.h:169
snd_config_get_bool = _lib.snd_config_get_bool
snd_config_get_bool.restype = c_int
snd_config_get_bool.argtypes = [POINTER(snd_config_t)]

# /usr/include/alsa/conf.h:170
snd_config_get_ctl_iface_ascii = _lib.snd_config_get_ctl_iface_ascii
snd_config_get_ctl_iface_ascii.restype = c_int
snd_config_get_ctl_iface_ascii.argtypes = [c_char_p]

# /usr/include/alsa/conf.h:171
snd_config_get_ctl_iface = _lib.snd_config_get_ctl_iface
snd_config_get_ctl_iface.restype = c_int
snd_config_get_ctl_iface.argtypes = [POINTER(snd_config_t)]

class struct_snd_devname(Structure):
    __slots__ = [
    ]
struct_snd_devname._fields_ = [
    ('_opaque_struct', c_int)
]

class struct_snd_devname(Structure):
    __slots__ = [
    ]
struct_snd_devname._fields_ = [
    ('_opaque_struct', c_int)
]

snd_devname_t = struct_snd_devname 	# /usr/include/alsa/conf.h:178
# /usr/include/alsa/conf.h:189
snd_names_list = _lib.snd_names_list
snd_names_list.restype = c_int
snd_names_list.argtypes = [c_char_p, POINTER(POINTER(snd_devname_t))]

# /usr/include/alsa/conf.h:190
snd_names_list_free = _lib.snd_names_list_free
snd_names_list_free.restype = None
snd_names_list_free.argtypes = [POINTER(snd_devname_t)]

SND_PCM_DLSYM_VERSION = 0 	# /usr/include/alsa/pcm.h:43
class struct__snd_pcm_info(Structure):
    __slots__ = [
    ]
struct__snd_pcm_info._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_pcm_info(Structure):
    __slots__ = [
    ]
struct__snd_pcm_info._fields_ = [
    ('_opaque_struct', c_int)
]

snd_pcm_info_t = struct__snd_pcm_info 	# /usr/include/alsa/pcm.h:46
class struct__snd_pcm_hw_params(Structure):
    __slots__ = [
    ]
struct__snd_pcm_hw_params._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_pcm_hw_params(Structure):
    __slots__ = [
    ]
struct__snd_pcm_hw_params._fields_ = [
    ('_opaque_struct', c_int)
]

snd_pcm_hw_params_t = struct__snd_pcm_hw_params 	# /usr/include/alsa/pcm.h:48
class struct__snd_pcm_sw_params(Structure):
    __slots__ = [
    ]
struct__snd_pcm_sw_params._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_pcm_sw_params(Structure):
    __slots__ = [
    ]
struct__snd_pcm_sw_params._fields_ = [
    ('_opaque_struct', c_int)
]

snd_pcm_sw_params_t = struct__snd_pcm_sw_params 	# /usr/include/alsa/pcm.h:50
class struct__snd_pcm_status(Structure):
    __slots__ = [
    ]
struct__snd_pcm_status._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_pcm_status(Structure):
    __slots__ = [
    ]
struct__snd_pcm_status._fields_ = [
    ('_opaque_struct', c_int)
]

snd_pcm_status_t = struct__snd_pcm_status 	# /usr/include/alsa/pcm.h:52
class struct__snd_pcm_access_mask(Structure):
    __slots__ = [
    ]
struct__snd_pcm_access_mask._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_pcm_access_mask(Structure):
    __slots__ = [
    ]
struct__snd_pcm_access_mask._fields_ = [
    ('_opaque_struct', c_int)
]

snd_pcm_access_mask_t = struct__snd_pcm_access_mask 	# /usr/include/alsa/pcm.h:54
class struct__snd_pcm_format_mask(Structure):
    __slots__ = [
    ]
struct__snd_pcm_format_mask._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_pcm_format_mask(Structure):
    __slots__ = [
    ]
struct__snd_pcm_format_mask._fields_ = [
    ('_opaque_struct', c_int)
]

snd_pcm_format_mask_t = struct__snd_pcm_format_mask 	# /usr/include/alsa/pcm.h:56
class struct__snd_pcm_subformat_mask(Structure):
    __slots__ = [
    ]
struct__snd_pcm_subformat_mask._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_pcm_subformat_mask(Structure):
    __slots__ = [
    ]
struct__snd_pcm_subformat_mask._fields_ = [
    ('_opaque_struct', c_int)
]

snd_pcm_subformat_mask_t = struct__snd_pcm_subformat_mask 	# /usr/include/alsa/pcm.h:58
enum__snd_pcm_class = c_int
SND_PCM_CLASS_GENERIC = 0
SND_PCM_CLASS_MULTI = 1
SND_PCM_CLASS_MODEM = 2
SND_PCM_CLASS_DIGITIZER = 3
SND_PCM_CLASS_LAST = 0
snd_pcm_class_t = enum__snd_pcm_class 	# /usr/include/alsa/pcm.h:72
enum__snd_pcm_subclass = c_int
SND_PCM_SUBCLASS_GENERIC_MIX = 0
SND_PCM_SUBCLASS_MULTI_MIX = 1
SND_PCM_SUBCLASS_LAST = 0
snd_pcm_subclass_t = enum__snd_pcm_subclass 	# /usr/include/alsa/pcm.h:81
enum__snd_pcm_stream = c_int
SND_PCM_STREAM_PLAYBACK = 0
SND_PCM_STREAM_CAPTURE = 1
SND_PCM_STREAM_LAST = 0
snd_pcm_stream_t = enum__snd_pcm_stream 	# /usr/include/alsa/pcm.h:90
enum__snd_pcm_access = c_int
SND_PCM_ACCESS_MMAP_INTERLEAVED = 0
SND_PCM_ACCESS_MMAP_NONINTERLEAVED = 1
SND_PCM_ACCESS_MMAP_COMPLEX = 2
SND_PCM_ACCESS_RW_INTERLEAVED = 3
SND_PCM_ACCESS_RW_NONINTERLEAVED = 4
SND_PCM_ACCESS_LAST = 0
snd_pcm_access_t = enum__snd_pcm_access 	# /usr/include/alsa/pcm.h:105
enum__snd_pcm_format = c_int
SND_PCM_FORMAT_UNKNOWN = 1
SND_PCM_FORMAT_S8 = 0
SND_PCM_FORMAT_U8 = 1
SND_PCM_FORMAT_S16_LE = 2
SND_PCM_FORMAT_S16_BE = 3
SND_PCM_FORMAT_U16_LE = 4
SND_PCM_FORMAT_U16_BE = 5
SND_PCM_FORMAT_S24_LE = 6
SND_PCM_FORMAT_S24_BE = 7
SND_PCM_FORMAT_U24_LE = 8
SND_PCM_FORMAT_U24_BE = 9
SND_PCM_FORMAT_S32_LE = 10
SND_PCM_FORMAT_S32_BE = 11
SND_PCM_FORMAT_U32_LE = 12
SND_PCM_FORMAT_U32_BE = 13
SND_PCM_FORMAT_FLOAT_LE = 14
SND_PCM_FORMAT_FLOAT_BE = 15
SND_PCM_FORMAT_FLOAT64_LE = 16
SND_PCM_FORMAT_FLOAT64_BE = 17
SND_PCM_FORMAT_IEC958_SUBFRAME_LE = 18
SND_PCM_FORMAT_IEC958_SUBFRAME_BE = 19
SND_PCM_FORMAT_MU_LAW = 20
SND_PCM_FORMAT_A_LAW = 21
SND_PCM_FORMAT_IMA_ADPCM = 22
SND_PCM_FORMAT_MPEG = 23
SND_PCM_FORMAT_GSM = 24
SND_PCM_FORMAT_SPECIAL = 31
SND_PCM_FORMAT_S24_3LE = 32
SND_PCM_FORMAT_S24_3BE = 33
SND_PCM_FORMAT_U24_3LE = 34
SND_PCM_FORMAT_U24_3BE = 35
SND_PCM_FORMAT_S20_3LE = 36
SND_PCM_FORMAT_S20_3BE = 37
SND_PCM_FORMAT_U20_3LE = 38
SND_PCM_FORMAT_U20_3BE = 39
SND_PCM_FORMAT_S18_3LE = 40
SND_PCM_FORMAT_S18_3BE = 41
SND_PCM_FORMAT_U18_3LE = 42
SND_PCM_FORMAT_U18_3BE = 43
SND_PCM_FORMAT_LAST = 0

# XXX wraptypes didn't pick up byte order detection
import sys
if sys.byteorder == 'little':
    SND_PCM_FORMAT_S16 = SND_PCM_FORMAT_S16_LE
    SND_PCM_FORMAT_U16 = SND_PCM_FORMAT_U16_LE
    SND_PCM_FORMAT_S24 = SND_PCM_FORMAT_S24_LE
    SND_PCM_FORMAT_U24 = SND_PCM_FORMAT_U24_LE
    SND_PCM_FORMAT_S32 = SND_PCM_FORMAT_S32_LE
    SND_PCM_FORMAT_U32 = SND_PCM_FORMAT_U32_LE
    SND_PCM_FORMAT_FLOAT = SND_PCM_FORMAT_FLOAT_LE
    SND_PCM_FORMAT_FLOAT64 = SND_PCM_FORMAT_FLOAT64_LE
    SND_PCM_FORMAT_IEC958_SUBFRAME = SND_PCM_FORMAT_IEC958_SUBFRAME_LE
else:
    SND_PCM_FORMAT_S16 = SND_PCM_FORMAT_S16_BE
    SND_PCM_FORMAT_U16 = SND_PCM_FORMAT_U16_BE
    SND_PCM_FORMAT_S24 = SND_PCM_FORMAT_S24_BE
    SND_PCM_FORMAT_U24 = SND_PCM_FORMAT_U24_BE
    SND_PCM_FORMAT_S32 = SND_PCM_FORMAT_S32_BE
    SND_PCM_FORMAT_U32 = SND_PCM_FORMAT_U32_BE
    SND_PCM_FORMAT_FLOAT = SND_PCM_FORMAT_FLOAT_BE
    SND_PCM_FORMAT_FLOAT64 = SND_PCM_FORMAT_FLOAT64_BE
    SND_PCM_FORMAT_IEC958_SUBFRAME = SND_PCM_FORMAT_IEC958_SUBFRAME_BE

snd_pcm_format_t = enum__snd_pcm_format 	# /usr/include/alsa/pcm.h:230
enum__snd_pcm_subformat = c_int
SND_PCM_SUBFORMAT_STD = 0
SND_PCM_SUBFORMAT_LAST = 0
snd_pcm_subformat_t = enum__snd_pcm_subformat 	# /usr/include/alsa/pcm.h:237
enum__snd_pcm_state = c_int
SND_PCM_STATE_OPEN = 0
SND_PCM_STATE_SETUP = 1
SND_PCM_STATE_PREPARED = 2
SND_PCM_STATE_RUNNING = 3
SND_PCM_STATE_XRUN = 4
SND_PCM_STATE_DRAINING = 5
SND_PCM_STATE_PAUSED = 6
SND_PCM_STATE_SUSPENDED = 7
SND_PCM_STATE_DISCONNECTED = 8
SND_PCM_STATE_LAST = 0
snd_pcm_state_t = enum__snd_pcm_state 	# /usr/include/alsa/pcm.h:260
enum__snd_pcm_start = c_int
SND_PCM_START_DATA = 0
SND_PCM_START_EXPLICIT = 1
SND_PCM_START_LAST = 0
snd_pcm_start_t = enum__snd_pcm_start 	# /usr/include/alsa/pcm.h:269
enum__snd_pcm_xrun = c_int
SND_PCM_XRUN_NONE = 0
SND_PCM_XRUN_STOP = 1
SND_PCM_XRUN_LAST = 0
snd_pcm_xrun_t = enum__snd_pcm_xrun 	# /usr/include/alsa/pcm.h:278
enum__snd_pcm_tstamp = c_int
SND_PCM_TSTAMP_NONE = 0
SND_PCM_TSTAMP_MMAP = 1
SND_PCM_TSTAMP_LAST = 0
snd_pcm_tstamp_t = enum__snd_pcm_tstamp 	# /usr/include/alsa/pcm.h:287
snd_pcm_uframes_t = c_ulong 	# /usr/include/alsa/pcm.h:290
snd_pcm_sframes_t = c_long 	# /usr/include/alsa/pcm.h:292
SND_PCM_NONBLOCK = 1 	# /usr/include/alsa/pcm.h:295
SND_PCM_ASYNC = 2 	# /usr/include/alsa/pcm.h:297
class struct__snd_pcm(Structure):
    __slots__ = [
    ]
struct__snd_pcm._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_pcm(Structure):
    __slots__ = [
    ]
struct__snd_pcm._fields_ = [
    ('_opaque_struct', c_int)
]

snd_pcm_t = struct__snd_pcm 	# /usr/include/alsa/pcm.h:300
enum__snd_pcm_type = c_int
snd_pcm_type_t = enum__snd_pcm_type 	# /usr/include/alsa/pcm.h:369
class struct__snd_pcm_channel_area(Structure):
    __slots__ = [
        'addr',
        'first',
        'step',
    ]
struct__snd_pcm_channel_area._fields_ = [
    ('addr', POINTER(None)),
    ('first', c_uint),
    ('step', c_uint),
]

snd_pcm_channel_area_t = struct__snd_pcm_channel_area 	# /usr/include/alsa/pcm.h:379
class struct__snd_pcm_sync_id(Union):
    __slots__ = [
        'id',
        'id16',
        'id32',
    ]
struct__snd_pcm_sync_id._fields_ = [
    ('id', c_ubyte * 16),
    ('id16', c_ushort * 8),
    ('id32', c_uint * 4),
]

snd_pcm_sync_id_t = struct__snd_pcm_sync_id 	# /usr/include/alsa/pcm.h:389
class struct__snd_pcm_scope(Structure):
    __slots__ = [
    ]
struct__snd_pcm_scope._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_pcm_scope(Structure):
    __slots__ = [
    ]
struct__snd_pcm_scope._fields_ = [
    ('_opaque_struct', c_int)
]

snd_pcm_scope_t = struct__snd_pcm_scope 	# /usr/include/alsa/pcm.h:392
# /usr/include/alsa/pcm.h:394
snd_pcm_open = _lib.snd_pcm_open
snd_pcm_open.restype = c_int
snd_pcm_open.argtypes = [POINTER(POINTER(snd_pcm_t)), c_char_p, snd_pcm_stream_t, c_int]

# /usr/include/alsa/pcm.h:396
snd_pcm_open_lconf = _lib.snd_pcm_open_lconf
snd_pcm_open_lconf.restype = c_int
snd_pcm_open_lconf.argtypes = [POINTER(POINTER(snd_pcm_t)), c_char_p, snd_pcm_stream_t, c_int, POINTER(snd_config_t)]

# /usr/include/alsa/pcm.h:400
snd_pcm_close = _lib.snd_pcm_close
snd_pcm_close.restype = c_int
snd_pcm_close.argtypes = [POINTER(snd_pcm_t)]

# /usr/include/alsa/pcm.h:401
snd_pcm_name = _lib.snd_pcm_name
snd_pcm_name.restype = c_char_p
snd_pcm_name.argtypes = [POINTER(snd_pcm_t)]

# /usr/include/alsa/pcm.h:402
snd_pcm_type = _lib.snd_pcm_type
snd_pcm_type.restype = snd_pcm_type_t
snd_pcm_type.argtypes = [POINTER(snd_pcm_t)]

# /usr/include/alsa/pcm.h:403
snd_pcm_stream = _lib.snd_pcm_stream
snd_pcm_stream.restype = snd_pcm_stream_t
snd_pcm_stream.argtypes = [POINTER(snd_pcm_t)]

# /usr/include/alsa/pcm.h:404
snd_pcm_poll_descriptors_count = _lib.snd_pcm_poll_descriptors_count
snd_pcm_poll_descriptors_count.restype = c_int
snd_pcm_poll_descriptors_count.argtypes = [POINTER(snd_pcm_t)]

class struct_pollfd(Structure):
    __slots__ = [
    ]
struct_pollfd._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/pcm.h:405
snd_pcm_poll_descriptors = _lib.snd_pcm_poll_descriptors
snd_pcm_poll_descriptors.restype = c_int
snd_pcm_poll_descriptors.argtypes = [POINTER(snd_pcm_t), POINTER(struct_pollfd), c_uint]

class struct_pollfd(Structure):
    __slots__ = [
    ]
struct_pollfd._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/pcm.h:406
snd_pcm_poll_descriptors_revents = _lib.snd_pcm_poll_descriptors_revents
snd_pcm_poll_descriptors_revents.restype = c_int
snd_pcm_poll_descriptors_revents.argtypes = [POINTER(snd_pcm_t), POINTER(struct_pollfd), c_uint, POINTER(c_ushort)]

# /usr/include/alsa/pcm.h:407
snd_pcm_nonblock = _lib.snd_pcm_nonblock
snd_pcm_nonblock.restype = c_int
snd_pcm_nonblock.argtypes = [POINTER(snd_pcm_t), c_int]

# /usr/include/alsa/pcm.h:408
snd_async_add_pcm_handler = _lib.snd_async_add_pcm_handler
snd_async_add_pcm_handler.restype = c_int
snd_async_add_pcm_handler.argtypes = [POINTER(POINTER(snd_async_handler_t)), POINTER(snd_pcm_t), snd_async_callback_t, POINTER(None)]

# /usr/include/alsa/pcm.h:410
snd_async_handler_get_pcm = _lib.snd_async_handler_get_pcm
snd_async_handler_get_pcm.restype = POINTER(snd_pcm_t)
snd_async_handler_get_pcm.argtypes = [POINTER(snd_async_handler_t)]

# /usr/include/alsa/pcm.h:411
snd_pcm_info = _lib.snd_pcm_info
snd_pcm_info.restype = c_int
snd_pcm_info.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_info_t)]

# /usr/include/alsa/pcm.h:412
snd_pcm_hw_params_current = _lib.snd_pcm_hw_params_current
snd_pcm_hw_params_current.restype = c_int
snd_pcm_hw_params_current.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t)]

# /usr/include/alsa/pcm.h:413
snd_pcm_hw_params = _lib.snd_pcm_hw_params
snd_pcm_hw_params.restype = c_int
snd_pcm_hw_params.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t)]

# /usr/include/alsa/pcm.h:414
snd_pcm_hw_free = _lib.snd_pcm_hw_free
snd_pcm_hw_free.restype = c_int
snd_pcm_hw_free.argtypes = [POINTER(snd_pcm_t)]

# /usr/include/alsa/pcm.h:415
snd_pcm_sw_params_current = _lib.snd_pcm_sw_params_current
snd_pcm_sw_params_current.restype = c_int
snd_pcm_sw_params_current.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_sw_params_t)]

# /usr/include/alsa/pcm.h:416
snd_pcm_sw_params = _lib.snd_pcm_sw_params
snd_pcm_sw_params.restype = c_int
snd_pcm_sw_params.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_sw_params_t)]

# /usr/include/alsa/pcm.h:417
snd_pcm_prepare = _lib.snd_pcm_prepare
snd_pcm_prepare.restype = c_int
snd_pcm_prepare.argtypes = [POINTER(snd_pcm_t)]

# /usr/include/alsa/pcm.h:418
snd_pcm_reset = _lib.snd_pcm_reset
snd_pcm_reset.restype = c_int
snd_pcm_reset.argtypes = [POINTER(snd_pcm_t)]

# /usr/include/alsa/pcm.h:419
snd_pcm_status = _lib.snd_pcm_status
snd_pcm_status.restype = c_int
snd_pcm_status.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_status_t)]

# /usr/include/alsa/pcm.h:420
snd_pcm_start = _lib.snd_pcm_start
snd_pcm_start.restype = c_int
snd_pcm_start.argtypes = [POINTER(snd_pcm_t)]

# /usr/include/alsa/pcm.h:421
snd_pcm_drop = _lib.snd_pcm_drop
snd_pcm_drop.restype = c_int
snd_pcm_drop.argtypes = [POINTER(snd_pcm_t)]

# /usr/include/alsa/pcm.h:422
snd_pcm_drain = _lib.snd_pcm_drain
snd_pcm_drain.restype = c_int
snd_pcm_drain.argtypes = [POINTER(snd_pcm_t)]

# /usr/include/alsa/pcm.h:423
snd_pcm_pause = _lib.snd_pcm_pause
snd_pcm_pause.restype = c_int
snd_pcm_pause.argtypes = [POINTER(snd_pcm_t), c_int]

# /usr/include/alsa/pcm.h:424
snd_pcm_state = _lib.snd_pcm_state
snd_pcm_state.restype = snd_pcm_state_t
snd_pcm_state.argtypes = [POINTER(snd_pcm_t)]

# /usr/include/alsa/pcm.h:425
snd_pcm_hwsync = _lib.snd_pcm_hwsync
snd_pcm_hwsync.restype = c_int
snd_pcm_hwsync.argtypes = [POINTER(snd_pcm_t)]

# /usr/include/alsa/pcm.h:426
snd_pcm_delay = _lib.snd_pcm_delay
snd_pcm_delay.restype = c_int
snd_pcm_delay.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_sframes_t)]

# /usr/include/alsa/pcm.h:427
snd_pcm_resume = _lib.snd_pcm_resume
snd_pcm_resume.restype = c_int
snd_pcm_resume.argtypes = [POINTER(snd_pcm_t)]

# /usr/include/alsa/pcm.h:428
snd_pcm_avail_update = _lib.snd_pcm_avail_update
snd_pcm_avail_update.restype = snd_pcm_sframes_t
snd_pcm_avail_update.argtypes = [POINTER(snd_pcm_t)]

# /usr/include/alsa/pcm.h:429
snd_pcm_rewind = _lib.snd_pcm_rewind
snd_pcm_rewind.restype = snd_pcm_sframes_t
snd_pcm_rewind.argtypes = [POINTER(snd_pcm_t), snd_pcm_uframes_t]

# /usr/include/alsa/pcm.h:430
snd_pcm_forward = _lib.snd_pcm_forward
snd_pcm_forward.restype = snd_pcm_sframes_t
snd_pcm_forward.argtypes = [POINTER(snd_pcm_t), snd_pcm_uframes_t]

# /usr/include/alsa/pcm.h:431
snd_pcm_writei = _lib.snd_pcm_writei
snd_pcm_writei.restype = snd_pcm_sframes_t
snd_pcm_writei.argtypes = [POINTER(snd_pcm_t), POINTER(None), snd_pcm_uframes_t]

# /usr/include/alsa/pcm.h:432
snd_pcm_readi = _lib.snd_pcm_readi
snd_pcm_readi.restype = snd_pcm_sframes_t
snd_pcm_readi.argtypes = [POINTER(snd_pcm_t), POINTER(None), snd_pcm_uframes_t]

# /usr/include/alsa/pcm.h:433
snd_pcm_writen = _lib.snd_pcm_writen
snd_pcm_writen.restype = snd_pcm_sframes_t
snd_pcm_writen.argtypes = [POINTER(snd_pcm_t), POINTER(POINTER(None)), snd_pcm_uframes_t]

# /usr/include/alsa/pcm.h:434
snd_pcm_readn = _lib.snd_pcm_readn
snd_pcm_readn.restype = snd_pcm_sframes_t
snd_pcm_readn.argtypes = [POINTER(snd_pcm_t), POINTER(POINTER(None)), snd_pcm_uframes_t]

# /usr/include/alsa/pcm.h:435
snd_pcm_wait = _lib.snd_pcm_wait
snd_pcm_wait.restype = c_int
snd_pcm_wait.argtypes = [POINTER(snd_pcm_t), c_int]

# /usr/include/alsa/pcm.h:437
snd_pcm_link = _lib.snd_pcm_link
snd_pcm_link.restype = c_int
snd_pcm_link.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_t)]

# /usr/include/alsa/pcm.h:438
snd_pcm_unlink = _lib.snd_pcm_unlink
snd_pcm_unlink.restype = c_int
snd_pcm_unlink.argtypes = [POINTER(snd_pcm_t)]

# /usr/include/alsa/pcm.h:445
snd_pcm_recover = _lib.snd_pcm_recover
snd_pcm_recover.restype = c_int
snd_pcm_recover.argtypes = [POINTER(snd_pcm_t), c_int, c_int]

# /usr/include/alsa/pcm.h:446
snd_pcm_set_params = _lib.snd_pcm_set_params
snd_pcm_set_params.restype = c_int
snd_pcm_set_params.argtypes = [POINTER(snd_pcm_t), snd_pcm_format_t, snd_pcm_access_t, c_uint, c_uint, c_int, c_uint]

# /usr/include/alsa/pcm.h:453
snd_pcm_get_params = _lib.snd_pcm_get_params
snd_pcm_get_params.restype = c_int
snd_pcm_get_params.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_uframes_t), POINTER(snd_pcm_uframes_t)]

# /usr/include/alsa/pcm.h:466
snd_pcm_info_sizeof = _lib.snd_pcm_info_sizeof
snd_pcm_info_sizeof.restype = c_size_t
snd_pcm_info_sizeof.argtypes = []

# /usr/include/alsa/pcm.h:472
snd_pcm_info_malloc = _lib.snd_pcm_info_malloc
snd_pcm_info_malloc.restype = c_int
snd_pcm_info_malloc.argtypes = [POINTER(POINTER(snd_pcm_info_t))]

# /usr/include/alsa/pcm.h:473
snd_pcm_info_free = _lib.snd_pcm_info_free
snd_pcm_info_free.restype = None
snd_pcm_info_free.argtypes = [POINTER(snd_pcm_info_t)]

# /usr/include/alsa/pcm.h:474
snd_pcm_info_copy = _lib.snd_pcm_info_copy
snd_pcm_info_copy.restype = None
snd_pcm_info_copy.argtypes = [POINTER(snd_pcm_info_t), POINTER(snd_pcm_info_t)]

# /usr/include/alsa/pcm.h:475
snd_pcm_info_get_device = _lib.snd_pcm_info_get_device
snd_pcm_info_get_device.restype = c_uint
snd_pcm_info_get_device.argtypes = [POINTER(snd_pcm_info_t)]

# /usr/include/alsa/pcm.h:476
snd_pcm_info_get_subdevice = _lib.snd_pcm_info_get_subdevice
snd_pcm_info_get_subdevice.restype = c_uint
snd_pcm_info_get_subdevice.argtypes = [POINTER(snd_pcm_info_t)]

# /usr/include/alsa/pcm.h:477
snd_pcm_info_get_stream = _lib.snd_pcm_info_get_stream
snd_pcm_info_get_stream.restype = snd_pcm_stream_t
snd_pcm_info_get_stream.argtypes = [POINTER(snd_pcm_info_t)]

# /usr/include/alsa/pcm.h:478
snd_pcm_info_get_card = _lib.snd_pcm_info_get_card
snd_pcm_info_get_card.restype = c_int
snd_pcm_info_get_card.argtypes = [POINTER(snd_pcm_info_t)]

# /usr/include/alsa/pcm.h:479
snd_pcm_info_get_id = _lib.snd_pcm_info_get_id
snd_pcm_info_get_id.restype = c_char_p
snd_pcm_info_get_id.argtypes = [POINTER(snd_pcm_info_t)]

# /usr/include/alsa/pcm.h:480
snd_pcm_info_get_name = _lib.snd_pcm_info_get_name
snd_pcm_info_get_name.restype = c_char_p
snd_pcm_info_get_name.argtypes = [POINTER(snd_pcm_info_t)]

# /usr/include/alsa/pcm.h:481
snd_pcm_info_get_subdevice_name = _lib.snd_pcm_info_get_subdevice_name
snd_pcm_info_get_subdevice_name.restype = c_char_p
snd_pcm_info_get_subdevice_name.argtypes = [POINTER(snd_pcm_info_t)]

# /usr/include/alsa/pcm.h:482
snd_pcm_info_get_class = _lib.snd_pcm_info_get_class
snd_pcm_info_get_class.restype = snd_pcm_class_t
snd_pcm_info_get_class.argtypes = [POINTER(snd_pcm_info_t)]

# /usr/include/alsa/pcm.h:483
snd_pcm_info_get_subclass = _lib.snd_pcm_info_get_subclass
snd_pcm_info_get_subclass.restype = snd_pcm_subclass_t
snd_pcm_info_get_subclass.argtypes = [POINTER(snd_pcm_info_t)]

# /usr/include/alsa/pcm.h:484
snd_pcm_info_get_subdevices_count = _lib.snd_pcm_info_get_subdevices_count
snd_pcm_info_get_subdevices_count.restype = c_uint
snd_pcm_info_get_subdevices_count.argtypes = [POINTER(snd_pcm_info_t)]

# /usr/include/alsa/pcm.h:485
snd_pcm_info_get_subdevices_avail = _lib.snd_pcm_info_get_subdevices_avail
snd_pcm_info_get_subdevices_avail.restype = c_uint
snd_pcm_info_get_subdevices_avail.argtypes = [POINTER(snd_pcm_info_t)]

# /usr/include/alsa/pcm.h:486
snd_pcm_info_get_sync = _lib.snd_pcm_info_get_sync
snd_pcm_info_get_sync.restype = snd_pcm_sync_id_t
snd_pcm_info_get_sync.argtypes = [POINTER(snd_pcm_info_t)]

# /usr/include/alsa/pcm.h:487
snd_pcm_info_set_device = _lib.snd_pcm_info_set_device
snd_pcm_info_set_device.restype = None
snd_pcm_info_set_device.argtypes = [POINTER(snd_pcm_info_t), c_uint]

# /usr/include/alsa/pcm.h:488
snd_pcm_info_set_subdevice = _lib.snd_pcm_info_set_subdevice
snd_pcm_info_set_subdevice.restype = None
snd_pcm_info_set_subdevice.argtypes = [POINTER(snd_pcm_info_t), c_uint]

# /usr/include/alsa/pcm.h:489
snd_pcm_info_set_stream = _lib.snd_pcm_info_set_stream
snd_pcm_info_set_stream.restype = None
snd_pcm_info_set_stream.argtypes = [POINTER(snd_pcm_info_t), snd_pcm_stream_t]

# /usr/include/alsa/pcm.h:500
snd_pcm_hw_params_any = _lib.snd_pcm_hw_params_any
snd_pcm_hw_params_any.restype = c_int
snd_pcm_hw_params_any.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t)]

# /usr/include/alsa/pcm.h:502
snd_pcm_hw_params_can_mmap_sample_resolution = _lib.snd_pcm_hw_params_can_mmap_sample_resolution
snd_pcm_hw_params_can_mmap_sample_resolution.restype = c_int
snd_pcm_hw_params_can_mmap_sample_resolution.argtypes = [POINTER(snd_pcm_hw_params_t)]

# /usr/include/alsa/pcm.h:503
snd_pcm_hw_params_is_double = _lib.snd_pcm_hw_params_is_double
snd_pcm_hw_params_is_double.restype = c_int
snd_pcm_hw_params_is_double.argtypes = [POINTER(snd_pcm_hw_params_t)]

# /usr/include/alsa/pcm.h:504
snd_pcm_hw_params_is_batch = _lib.snd_pcm_hw_params_is_batch
snd_pcm_hw_params_is_batch.restype = c_int
snd_pcm_hw_params_is_batch.argtypes = [POINTER(snd_pcm_hw_params_t)]

# /usr/include/alsa/pcm.h:505
snd_pcm_hw_params_is_block_transfer = _lib.snd_pcm_hw_params_is_block_transfer
snd_pcm_hw_params_is_block_transfer.restype = c_int
snd_pcm_hw_params_is_block_transfer.argtypes = [POINTER(snd_pcm_hw_params_t)]

# /usr/include/alsa/pcm.h:506
snd_pcm_hw_params_can_overrange = _lib.snd_pcm_hw_params_can_overrange
snd_pcm_hw_params_can_overrange.restype = c_int
snd_pcm_hw_params_can_overrange.argtypes = [POINTER(snd_pcm_hw_params_t)]

# /usr/include/alsa/pcm.h:507
snd_pcm_hw_params_can_pause = _lib.snd_pcm_hw_params_can_pause
snd_pcm_hw_params_can_pause.restype = c_int
snd_pcm_hw_params_can_pause.argtypes = [POINTER(snd_pcm_hw_params_t)]

# /usr/include/alsa/pcm.h:508
snd_pcm_hw_params_can_resume = _lib.snd_pcm_hw_params_can_resume
snd_pcm_hw_params_can_resume.restype = c_int
snd_pcm_hw_params_can_resume.argtypes = [POINTER(snd_pcm_hw_params_t)]

# /usr/include/alsa/pcm.h:509
snd_pcm_hw_params_is_half_duplex = _lib.snd_pcm_hw_params_is_half_duplex
snd_pcm_hw_params_is_half_duplex.restype = c_int
snd_pcm_hw_params_is_half_duplex.argtypes = [POINTER(snd_pcm_hw_params_t)]

# /usr/include/alsa/pcm.h:510
snd_pcm_hw_params_is_joint_duplex = _lib.snd_pcm_hw_params_is_joint_duplex
snd_pcm_hw_params_is_joint_duplex.restype = c_int
snd_pcm_hw_params_is_joint_duplex.argtypes = [POINTER(snd_pcm_hw_params_t)]

# /usr/include/alsa/pcm.h:511
snd_pcm_hw_params_can_sync_start = _lib.snd_pcm_hw_params_can_sync_start
snd_pcm_hw_params_can_sync_start.restype = c_int
snd_pcm_hw_params_can_sync_start.argtypes = [POINTER(snd_pcm_hw_params_t)]

# /usr/include/alsa/pcm.h:512
snd_pcm_hw_params_get_rate_numden = _lib.snd_pcm_hw_params_get_rate_numden
snd_pcm_hw_params_get_rate_numden.restype = c_int
snd_pcm_hw_params_get_rate_numden.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_uint)]

# /usr/include/alsa/pcm.h:515
snd_pcm_hw_params_get_sbits = _lib.snd_pcm_hw_params_get_sbits
snd_pcm_hw_params_get_sbits.restype = c_int
snd_pcm_hw_params_get_sbits.argtypes = [POINTER(snd_pcm_hw_params_t)]

# /usr/include/alsa/pcm.h:516
snd_pcm_hw_params_get_fifo_size = _lib.snd_pcm_hw_params_get_fifo_size
snd_pcm_hw_params_get_fifo_size.restype = c_int
snd_pcm_hw_params_get_fifo_size.argtypes = [POINTER(snd_pcm_hw_params_t)]

# /usr/include/alsa/pcm.h:544
snd_pcm_hw_params_sizeof = _lib.snd_pcm_hw_params_sizeof
snd_pcm_hw_params_sizeof.restype = c_size_t
snd_pcm_hw_params_sizeof.argtypes = []

# /usr/include/alsa/pcm.h:550
snd_pcm_hw_params_malloc = _lib.snd_pcm_hw_params_malloc
snd_pcm_hw_params_malloc.restype = c_int
snd_pcm_hw_params_malloc.argtypes = [POINTER(POINTER(snd_pcm_hw_params_t))]

# /usr/include/alsa/pcm.h:551
snd_pcm_hw_params_free = _lib.snd_pcm_hw_params_free
snd_pcm_hw_params_free.restype = None
snd_pcm_hw_params_free.argtypes = [POINTER(snd_pcm_hw_params_t)]

# /usr/include/alsa/pcm.h:552
snd_pcm_hw_params_copy = _lib.snd_pcm_hw_params_copy
snd_pcm_hw_params_copy.restype = None
snd_pcm_hw_params_copy.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_hw_params_t)]

# /usr/include/alsa/pcm.h:556
snd_pcm_hw_params_get_access = _lib.snd_pcm_hw_params_get_access
snd_pcm_hw_params_get_access.restype = c_int
snd_pcm_hw_params_get_access.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_access_t)]

# /usr/include/alsa/pcm.h:557
snd_pcm_hw_params_test_access = _lib.snd_pcm_hw_params_test_access
snd_pcm_hw_params_test_access.restype = c_int
snd_pcm_hw_params_test_access.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), snd_pcm_access_t]

# /usr/include/alsa/pcm.h:558
snd_pcm_hw_params_set_access = _lib.snd_pcm_hw_params_set_access
snd_pcm_hw_params_set_access.restype = c_int
snd_pcm_hw_params_set_access.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), snd_pcm_access_t]

# /usr/include/alsa/pcm.h:559
snd_pcm_hw_params_set_access_first = _lib.snd_pcm_hw_params_set_access_first
snd_pcm_hw_params_set_access_first.restype = c_int
snd_pcm_hw_params_set_access_first.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_access_t)]

# /usr/include/alsa/pcm.h:560
snd_pcm_hw_params_set_access_last = _lib.snd_pcm_hw_params_set_access_last
snd_pcm_hw_params_set_access_last.restype = c_int
snd_pcm_hw_params_set_access_last.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_access_t)]

# /usr/include/alsa/pcm.h:561
snd_pcm_hw_params_set_access_mask = _lib.snd_pcm_hw_params_set_access_mask
snd_pcm_hw_params_set_access_mask.restype = c_int
snd_pcm_hw_params_set_access_mask.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_access_mask_t)]

# /usr/include/alsa/pcm.h:562
snd_pcm_hw_params_get_access_mask = _lib.snd_pcm_hw_params_get_access_mask
snd_pcm_hw_params_get_access_mask.restype = c_int
snd_pcm_hw_params_get_access_mask.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_access_mask_t)]

# /usr/include/alsa/pcm.h:564
snd_pcm_hw_params_get_format = _lib.snd_pcm_hw_params_get_format
snd_pcm_hw_params_get_format.restype = c_int
snd_pcm_hw_params_get_format.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_format_t)]

# /usr/include/alsa/pcm.h:565
snd_pcm_hw_params_test_format = _lib.snd_pcm_hw_params_test_format
snd_pcm_hw_params_test_format.restype = c_int
snd_pcm_hw_params_test_format.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), snd_pcm_format_t]

# /usr/include/alsa/pcm.h:566
snd_pcm_hw_params_set_format = _lib.snd_pcm_hw_params_set_format
snd_pcm_hw_params_set_format.restype = c_int
snd_pcm_hw_params_set_format.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), snd_pcm_format_t]

# /usr/include/alsa/pcm.h:567
snd_pcm_hw_params_set_format_first = _lib.snd_pcm_hw_params_set_format_first
snd_pcm_hw_params_set_format_first.restype = c_int
snd_pcm_hw_params_set_format_first.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_format_t)]

# /usr/include/alsa/pcm.h:568
snd_pcm_hw_params_set_format_last = _lib.snd_pcm_hw_params_set_format_last
snd_pcm_hw_params_set_format_last.restype = c_int
snd_pcm_hw_params_set_format_last.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_format_t)]

# /usr/include/alsa/pcm.h:569
snd_pcm_hw_params_set_format_mask = _lib.snd_pcm_hw_params_set_format_mask
snd_pcm_hw_params_set_format_mask.restype = c_int
snd_pcm_hw_params_set_format_mask.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_format_mask_t)]

# /usr/include/alsa/pcm.h:570
snd_pcm_hw_params_get_format_mask = _lib.snd_pcm_hw_params_get_format_mask
snd_pcm_hw_params_get_format_mask.restype = None
snd_pcm_hw_params_get_format_mask.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_format_mask_t)]

# /usr/include/alsa/pcm.h:572
snd_pcm_hw_params_get_subformat = _lib.snd_pcm_hw_params_get_subformat
snd_pcm_hw_params_get_subformat.restype = c_int
snd_pcm_hw_params_get_subformat.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_subformat_t)]

# /usr/include/alsa/pcm.h:573
snd_pcm_hw_params_test_subformat = _lib.snd_pcm_hw_params_test_subformat
snd_pcm_hw_params_test_subformat.restype = c_int
snd_pcm_hw_params_test_subformat.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), snd_pcm_subformat_t]

# /usr/include/alsa/pcm.h:574
snd_pcm_hw_params_set_subformat = _lib.snd_pcm_hw_params_set_subformat
snd_pcm_hw_params_set_subformat.restype = c_int
snd_pcm_hw_params_set_subformat.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), snd_pcm_subformat_t]

# /usr/include/alsa/pcm.h:575
snd_pcm_hw_params_set_subformat_first = _lib.snd_pcm_hw_params_set_subformat_first
snd_pcm_hw_params_set_subformat_first.restype = c_int
snd_pcm_hw_params_set_subformat_first.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_subformat_t)]

# /usr/include/alsa/pcm.h:576
snd_pcm_hw_params_set_subformat_last = _lib.snd_pcm_hw_params_set_subformat_last
snd_pcm_hw_params_set_subformat_last.restype = c_int
snd_pcm_hw_params_set_subformat_last.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_subformat_t)]

# /usr/include/alsa/pcm.h:577
snd_pcm_hw_params_set_subformat_mask = _lib.snd_pcm_hw_params_set_subformat_mask
snd_pcm_hw_params_set_subformat_mask.restype = c_int
snd_pcm_hw_params_set_subformat_mask.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_subformat_mask_t)]

# /usr/include/alsa/pcm.h:578
snd_pcm_hw_params_get_subformat_mask = _lib.snd_pcm_hw_params_get_subformat_mask
snd_pcm_hw_params_get_subformat_mask.restype = None
snd_pcm_hw_params_get_subformat_mask.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_subformat_mask_t)]

# /usr/include/alsa/pcm.h:580
snd_pcm_hw_params_get_channels = _lib.snd_pcm_hw_params_get_channels
snd_pcm_hw_params_get_channels.restype = c_int
snd_pcm_hw_params_get_channels.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(c_uint)]

# /usr/include/alsa/pcm.h:581
snd_pcm_hw_params_get_channels_min = _lib.snd_pcm_hw_params_get_channels_min
snd_pcm_hw_params_get_channels_min.restype = c_int
snd_pcm_hw_params_get_channels_min.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(c_uint)]

# /usr/include/alsa/pcm.h:582
snd_pcm_hw_params_get_channels_max = _lib.snd_pcm_hw_params_get_channels_max
snd_pcm_hw_params_get_channels_max.restype = c_int
snd_pcm_hw_params_get_channels_max.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(c_uint)]

# /usr/include/alsa/pcm.h:583
snd_pcm_hw_params_test_channels = _lib.snd_pcm_hw_params_test_channels
snd_pcm_hw_params_test_channels.restype = c_int
snd_pcm_hw_params_test_channels.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), c_uint]

# /usr/include/alsa/pcm.h:584
snd_pcm_hw_params_set_channels = _lib.snd_pcm_hw_params_set_channels
snd_pcm_hw_params_set_channels.restype = c_int
snd_pcm_hw_params_set_channels.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), c_uint]

# /usr/include/alsa/pcm.h:585
snd_pcm_hw_params_set_channels_min = _lib.snd_pcm_hw_params_set_channels_min
snd_pcm_hw_params_set_channels_min.restype = c_int
snd_pcm_hw_params_set_channels_min.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint)]

# /usr/include/alsa/pcm.h:586
snd_pcm_hw_params_set_channels_max = _lib.snd_pcm_hw_params_set_channels_max
snd_pcm_hw_params_set_channels_max.restype = c_int
snd_pcm_hw_params_set_channels_max.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint)]

# /usr/include/alsa/pcm.h:587
snd_pcm_hw_params_set_channels_minmax = _lib.snd_pcm_hw_params_set_channels_minmax
snd_pcm_hw_params_set_channels_minmax.restype = c_int
snd_pcm_hw_params_set_channels_minmax.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_uint)]

# /usr/include/alsa/pcm.h:588
snd_pcm_hw_params_set_channels_near = _lib.snd_pcm_hw_params_set_channels_near
snd_pcm_hw_params_set_channels_near.restype = c_int
snd_pcm_hw_params_set_channels_near.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint)]

# /usr/include/alsa/pcm.h:589
snd_pcm_hw_params_set_channels_first = _lib.snd_pcm_hw_params_set_channels_first
snd_pcm_hw_params_set_channels_first.restype = c_int
snd_pcm_hw_params_set_channels_first.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint)]

# /usr/include/alsa/pcm.h:590
snd_pcm_hw_params_set_channels_last = _lib.snd_pcm_hw_params_set_channels_last
snd_pcm_hw_params_set_channels_last.restype = c_int
snd_pcm_hw_params_set_channels_last.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint)]

# /usr/include/alsa/pcm.h:592
snd_pcm_hw_params_get_rate = _lib.snd_pcm_hw_params_get_rate
snd_pcm_hw_params_get_rate.restype = c_int
snd_pcm_hw_params_get_rate.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:593
snd_pcm_hw_params_get_rate_min = _lib.snd_pcm_hw_params_get_rate_min
snd_pcm_hw_params_get_rate_min.restype = c_int
snd_pcm_hw_params_get_rate_min.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:594
snd_pcm_hw_params_get_rate_max = _lib.snd_pcm_hw_params_get_rate_max
snd_pcm_hw_params_get_rate_max.restype = c_int
snd_pcm_hw_params_get_rate_max.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:595
snd_pcm_hw_params_test_rate = _lib.snd_pcm_hw_params_test_rate
snd_pcm_hw_params_test_rate.restype = c_int
snd_pcm_hw_params_test_rate.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), c_uint, c_int]

# /usr/include/alsa/pcm.h:596
snd_pcm_hw_params_set_rate = _lib.snd_pcm_hw_params_set_rate
snd_pcm_hw_params_set_rate.restype = c_int
snd_pcm_hw_params_set_rate.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), c_uint, c_int]

# /usr/include/alsa/pcm.h:597
snd_pcm_hw_params_set_rate_min = _lib.snd_pcm_hw_params_set_rate_min
snd_pcm_hw_params_set_rate_min.restype = c_int
snd_pcm_hw_params_set_rate_min.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:598
snd_pcm_hw_params_set_rate_max = _lib.snd_pcm_hw_params_set_rate_max
snd_pcm_hw_params_set_rate_max.restype = c_int
snd_pcm_hw_params_set_rate_max.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:599
snd_pcm_hw_params_set_rate_minmax = _lib.snd_pcm_hw_params_set_rate_minmax
snd_pcm_hw_params_set_rate_minmax.restype = c_int
snd_pcm_hw_params_set_rate_minmax.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:600
snd_pcm_hw_params_set_rate_near = _lib.snd_pcm_hw_params_set_rate_near
snd_pcm_hw_params_set_rate_near.restype = c_int
snd_pcm_hw_params_set_rate_near.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:601
snd_pcm_hw_params_set_rate_first = _lib.snd_pcm_hw_params_set_rate_first
snd_pcm_hw_params_set_rate_first.restype = c_int
snd_pcm_hw_params_set_rate_first.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:602
snd_pcm_hw_params_set_rate_last = _lib.snd_pcm_hw_params_set_rate_last
snd_pcm_hw_params_set_rate_last.restype = c_int
snd_pcm_hw_params_set_rate_last.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:603
snd_pcm_hw_params_set_rate_resample = _lib.snd_pcm_hw_params_set_rate_resample
snd_pcm_hw_params_set_rate_resample.restype = c_int
snd_pcm_hw_params_set_rate_resample.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), c_uint]

# /usr/include/alsa/pcm.h:604
snd_pcm_hw_params_get_rate_resample = _lib.snd_pcm_hw_params_get_rate_resample
snd_pcm_hw_params_get_rate_resample.restype = c_int
snd_pcm_hw_params_get_rate_resample.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint)]

# /usr/include/alsa/pcm.h:605
snd_pcm_hw_params_set_export_buffer = _lib.snd_pcm_hw_params_set_export_buffer
snd_pcm_hw_params_set_export_buffer.restype = c_int
snd_pcm_hw_params_set_export_buffer.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), c_uint]

# /usr/include/alsa/pcm.h:606
snd_pcm_hw_params_get_export_buffer = _lib.snd_pcm_hw_params_get_export_buffer
snd_pcm_hw_params_get_export_buffer.restype = c_int
snd_pcm_hw_params_get_export_buffer.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint)]

# /usr/include/alsa/pcm.h:608
snd_pcm_hw_params_get_period_time = _lib.snd_pcm_hw_params_get_period_time
snd_pcm_hw_params_get_period_time.restype = c_int
snd_pcm_hw_params_get_period_time.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:609
snd_pcm_hw_params_get_period_time_min = _lib.snd_pcm_hw_params_get_period_time_min
snd_pcm_hw_params_get_period_time_min.restype = c_int
snd_pcm_hw_params_get_period_time_min.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:610
snd_pcm_hw_params_get_period_time_max = _lib.snd_pcm_hw_params_get_period_time_max
snd_pcm_hw_params_get_period_time_max.restype = c_int
snd_pcm_hw_params_get_period_time_max.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:611
snd_pcm_hw_params_test_period_time = _lib.snd_pcm_hw_params_test_period_time
snd_pcm_hw_params_test_period_time.restype = c_int
snd_pcm_hw_params_test_period_time.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), c_uint, c_int]

# /usr/include/alsa/pcm.h:612
snd_pcm_hw_params_set_period_time = _lib.snd_pcm_hw_params_set_period_time
snd_pcm_hw_params_set_period_time.restype = c_int
snd_pcm_hw_params_set_period_time.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), c_uint, c_int]

# /usr/include/alsa/pcm.h:613
snd_pcm_hw_params_set_period_time_min = _lib.snd_pcm_hw_params_set_period_time_min
snd_pcm_hw_params_set_period_time_min.restype = c_int
snd_pcm_hw_params_set_period_time_min.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:614
snd_pcm_hw_params_set_period_time_max = _lib.snd_pcm_hw_params_set_period_time_max
snd_pcm_hw_params_set_period_time_max.restype = c_int
snd_pcm_hw_params_set_period_time_max.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:615
snd_pcm_hw_params_set_period_time_minmax = _lib.snd_pcm_hw_params_set_period_time_minmax
snd_pcm_hw_params_set_period_time_minmax.restype = c_int
snd_pcm_hw_params_set_period_time_minmax.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:616
snd_pcm_hw_params_set_period_time_near = _lib.snd_pcm_hw_params_set_period_time_near
snd_pcm_hw_params_set_period_time_near.restype = c_int
snd_pcm_hw_params_set_period_time_near.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:617
snd_pcm_hw_params_set_period_time_first = _lib.snd_pcm_hw_params_set_period_time_first
snd_pcm_hw_params_set_period_time_first.restype = c_int
snd_pcm_hw_params_set_period_time_first.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:618
snd_pcm_hw_params_set_period_time_last = _lib.snd_pcm_hw_params_set_period_time_last
snd_pcm_hw_params_set_period_time_last.restype = c_int
snd_pcm_hw_params_set_period_time_last.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:620
snd_pcm_hw_params_get_period_size = _lib.snd_pcm_hw_params_get_period_size
snd_pcm_hw_params_get_period_size.restype = c_int
snd_pcm_hw_params_get_period_size.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_uframes_t), POINTER(c_int)]

# /usr/include/alsa/pcm.h:621
snd_pcm_hw_params_get_period_size_min = _lib.snd_pcm_hw_params_get_period_size_min
snd_pcm_hw_params_get_period_size_min.restype = c_int
snd_pcm_hw_params_get_period_size_min.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_uframes_t), POINTER(c_int)]

# /usr/include/alsa/pcm.h:622
snd_pcm_hw_params_get_period_size_max = _lib.snd_pcm_hw_params_get_period_size_max
snd_pcm_hw_params_get_period_size_max.restype = c_int
snd_pcm_hw_params_get_period_size_max.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_uframes_t), POINTER(c_int)]

# /usr/include/alsa/pcm.h:623
snd_pcm_hw_params_test_period_size = _lib.snd_pcm_hw_params_test_period_size
snd_pcm_hw_params_test_period_size.restype = c_int
snd_pcm_hw_params_test_period_size.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), snd_pcm_uframes_t, c_int]

# /usr/include/alsa/pcm.h:624
snd_pcm_hw_params_set_period_size = _lib.snd_pcm_hw_params_set_period_size
snd_pcm_hw_params_set_period_size.restype = c_int
snd_pcm_hw_params_set_period_size.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), snd_pcm_uframes_t, c_int]

# /usr/include/alsa/pcm.h:625
snd_pcm_hw_params_set_period_size_min = _lib.snd_pcm_hw_params_set_period_size_min
snd_pcm_hw_params_set_period_size_min.restype = c_int
snd_pcm_hw_params_set_period_size_min.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_uframes_t), POINTER(c_int)]

# /usr/include/alsa/pcm.h:626
snd_pcm_hw_params_set_period_size_max = _lib.snd_pcm_hw_params_set_period_size_max
snd_pcm_hw_params_set_period_size_max.restype = c_int
snd_pcm_hw_params_set_period_size_max.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_uframes_t), POINTER(c_int)]

# /usr/include/alsa/pcm.h:627
snd_pcm_hw_params_set_period_size_minmax = _lib.snd_pcm_hw_params_set_period_size_minmax
snd_pcm_hw_params_set_period_size_minmax.restype = c_int
snd_pcm_hw_params_set_period_size_minmax.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_uframes_t), POINTER(c_int), POINTER(snd_pcm_uframes_t), POINTER(c_int)]

# /usr/include/alsa/pcm.h:628
snd_pcm_hw_params_set_period_size_near = _lib.snd_pcm_hw_params_set_period_size_near
snd_pcm_hw_params_set_period_size_near.restype = c_int
snd_pcm_hw_params_set_period_size_near.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_uframes_t), POINTER(c_int)]

# /usr/include/alsa/pcm.h:629
snd_pcm_hw_params_set_period_size_first = _lib.snd_pcm_hw_params_set_period_size_first
snd_pcm_hw_params_set_period_size_first.restype = c_int
snd_pcm_hw_params_set_period_size_first.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_uframes_t), POINTER(c_int)]

# /usr/include/alsa/pcm.h:630
snd_pcm_hw_params_set_period_size_last = _lib.snd_pcm_hw_params_set_period_size_last
snd_pcm_hw_params_set_period_size_last.restype = c_int
snd_pcm_hw_params_set_period_size_last.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_uframes_t), POINTER(c_int)]

# /usr/include/alsa/pcm.h:631
snd_pcm_hw_params_set_period_size_integer = _lib.snd_pcm_hw_params_set_period_size_integer
snd_pcm_hw_params_set_period_size_integer.restype = c_int
snd_pcm_hw_params_set_period_size_integer.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t)]

# /usr/include/alsa/pcm.h:633
snd_pcm_hw_params_get_periods = _lib.snd_pcm_hw_params_get_periods
snd_pcm_hw_params_get_periods.restype = c_int
snd_pcm_hw_params_get_periods.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:634
snd_pcm_hw_params_get_periods_min = _lib.snd_pcm_hw_params_get_periods_min
snd_pcm_hw_params_get_periods_min.restype = c_int
snd_pcm_hw_params_get_periods_min.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:635
snd_pcm_hw_params_get_periods_max = _lib.snd_pcm_hw_params_get_periods_max
snd_pcm_hw_params_get_periods_max.restype = c_int
snd_pcm_hw_params_get_periods_max.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:636
snd_pcm_hw_params_test_periods = _lib.snd_pcm_hw_params_test_periods
snd_pcm_hw_params_test_periods.restype = c_int
snd_pcm_hw_params_test_periods.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), c_uint, c_int]

# /usr/include/alsa/pcm.h:637
snd_pcm_hw_params_set_periods = _lib.snd_pcm_hw_params_set_periods
snd_pcm_hw_params_set_periods.restype = c_int
snd_pcm_hw_params_set_periods.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), c_uint, c_int]

# /usr/include/alsa/pcm.h:638
snd_pcm_hw_params_set_periods_min = _lib.snd_pcm_hw_params_set_periods_min
snd_pcm_hw_params_set_periods_min.restype = c_int
snd_pcm_hw_params_set_periods_min.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:639
snd_pcm_hw_params_set_periods_max = _lib.snd_pcm_hw_params_set_periods_max
snd_pcm_hw_params_set_periods_max.restype = c_int
snd_pcm_hw_params_set_periods_max.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:640
snd_pcm_hw_params_set_periods_minmax = _lib.snd_pcm_hw_params_set_periods_minmax
snd_pcm_hw_params_set_periods_minmax.restype = c_int
snd_pcm_hw_params_set_periods_minmax.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:641
snd_pcm_hw_params_set_periods_near = _lib.snd_pcm_hw_params_set_periods_near
snd_pcm_hw_params_set_periods_near.restype = c_int
snd_pcm_hw_params_set_periods_near.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:642
snd_pcm_hw_params_set_periods_first = _lib.snd_pcm_hw_params_set_periods_first
snd_pcm_hw_params_set_periods_first.restype = c_int
snd_pcm_hw_params_set_periods_first.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:643
snd_pcm_hw_params_set_periods_last = _lib.snd_pcm_hw_params_set_periods_last
snd_pcm_hw_params_set_periods_last.restype = c_int
snd_pcm_hw_params_set_periods_last.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:644
snd_pcm_hw_params_set_periods_integer = _lib.snd_pcm_hw_params_set_periods_integer
snd_pcm_hw_params_set_periods_integer.restype = c_int
snd_pcm_hw_params_set_periods_integer.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t)]

# /usr/include/alsa/pcm.h:646
snd_pcm_hw_params_get_buffer_time = _lib.snd_pcm_hw_params_get_buffer_time
snd_pcm_hw_params_get_buffer_time.restype = c_int
snd_pcm_hw_params_get_buffer_time.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:647
snd_pcm_hw_params_get_buffer_time_min = _lib.snd_pcm_hw_params_get_buffer_time_min
snd_pcm_hw_params_get_buffer_time_min.restype = c_int
snd_pcm_hw_params_get_buffer_time_min.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:648
snd_pcm_hw_params_get_buffer_time_max = _lib.snd_pcm_hw_params_get_buffer_time_max
snd_pcm_hw_params_get_buffer_time_max.restype = c_int
snd_pcm_hw_params_get_buffer_time_max.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:649
snd_pcm_hw_params_test_buffer_time = _lib.snd_pcm_hw_params_test_buffer_time
snd_pcm_hw_params_test_buffer_time.restype = c_int
snd_pcm_hw_params_test_buffer_time.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), c_uint, c_int]

# /usr/include/alsa/pcm.h:650
snd_pcm_hw_params_set_buffer_time = _lib.snd_pcm_hw_params_set_buffer_time
snd_pcm_hw_params_set_buffer_time.restype = c_int
snd_pcm_hw_params_set_buffer_time.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), c_uint, c_int]

# /usr/include/alsa/pcm.h:651
snd_pcm_hw_params_set_buffer_time_min = _lib.snd_pcm_hw_params_set_buffer_time_min
snd_pcm_hw_params_set_buffer_time_min.restype = c_int
snd_pcm_hw_params_set_buffer_time_min.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:652
snd_pcm_hw_params_set_buffer_time_max = _lib.snd_pcm_hw_params_set_buffer_time_max
snd_pcm_hw_params_set_buffer_time_max.restype = c_int
snd_pcm_hw_params_set_buffer_time_max.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:653
snd_pcm_hw_params_set_buffer_time_minmax = _lib.snd_pcm_hw_params_set_buffer_time_minmax
snd_pcm_hw_params_set_buffer_time_minmax.restype = c_int
snd_pcm_hw_params_set_buffer_time_minmax.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:654
snd_pcm_hw_params_set_buffer_time_near = _lib.snd_pcm_hw_params_set_buffer_time_near
snd_pcm_hw_params_set_buffer_time_near.restype = c_int
snd_pcm_hw_params_set_buffer_time_near.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:655
snd_pcm_hw_params_set_buffer_time_first = _lib.snd_pcm_hw_params_set_buffer_time_first
snd_pcm_hw_params_set_buffer_time_first.restype = c_int
snd_pcm_hw_params_set_buffer_time_first.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:656
snd_pcm_hw_params_set_buffer_time_last = _lib.snd_pcm_hw_params_set_buffer_time_last
snd_pcm_hw_params_set_buffer_time_last.restype = c_int
snd_pcm_hw_params_set_buffer_time_last.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:658
snd_pcm_hw_params_get_buffer_size = _lib.snd_pcm_hw_params_get_buffer_size
snd_pcm_hw_params_get_buffer_size.restype = c_int
snd_pcm_hw_params_get_buffer_size.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_uframes_t)]

# /usr/include/alsa/pcm.h:659
snd_pcm_hw_params_get_buffer_size_min = _lib.snd_pcm_hw_params_get_buffer_size_min
snd_pcm_hw_params_get_buffer_size_min.restype = c_int
snd_pcm_hw_params_get_buffer_size_min.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_uframes_t)]

# /usr/include/alsa/pcm.h:660
snd_pcm_hw_params_get_buffer_size_max = _lib.snd_pcm_hw_params_get_buffer_size_max
snd_pcm_hw_params_get_buffer_size_max.restype = c_int
snd_pcm_hw_params_get_buffer_size_max.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_uframes_t)]

# /usr/include/alsa/pcm.h:661
snd_pcm_hw_params_test_buffer_size = _lib.snd_pcm_hw_params_test_buffer_size
snd_pcm_hw_params_test_buffer_size.restype = c_int
snd_pcm_hw_params_test_buffer_size.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), snd_pcm_uframes_t]

# /usr/include/alsa/pcm.h:662
snd_pcm_hw_params_set_buffer_size = _lib.snd_pcm_hw_params_set_buffer_size
snd_pcm_hw_params_set_buffer_size.restype = c_int
snd_pcm_hw_params_set_buffer_size.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), snd_pcm_uframes_t]

# /usr/include/alsa/pcm.h:663
snd_pcm_hw_params_set_buffer_size_min = _lib.snd_pcm_hw_params_set_buffer_size_min
snd_pcm_hw_params_set_buffer_size_min.restype = c_int
snd_pcm_hw_params_set_buffer_size_min.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_uframes_t)]

# /usr/include/alsa/pcm.h:664
snd_pcm_hw_params_set_buffer_size_max = _lib.snd_pcm_hw_params_set_buffer_size_max
snd_pcm_hw_params_set_buffer_size_max.restype = c_int
snd_pcm_hw_params_set_buffer_size_max.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_uframes_t)]

# /usr/include/alsa/pcm.h:665
snd_pcm_hw_params_set_buffer_size_minmax = _lib.snd_pcm_hw_params_set_buffer_size_minmax
snd_pcm_hw_params_set_buffer_size_minmax.restype = c_int
snd_pcm_hw_params_set_buffer_size_minmax.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_uframes_t), POINTER(snd_pcm_uframes_t)]

# /usr/include/alsa/pcm.h:666
snd_pcm_hw_params_set_buffer_size_near = _lib.snd_pcm_hw_params_set_buffer_size_near
snd_pcm_hw_params_set_buffer_size_near.restype = c_int
snd_pcm_hw_params_set_buffer_size_near.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_uframes_t)]

# /usr/include/alsa/pcm.h:667
snd_pcm_hw_params_set_buffer_size_first = _lib.snd_pcm_hw_params_set_buffer_size_first
snd_pcm_hw_params_set_buffer_size_first.restype = c_int
snd_pcm_hw_params_set_buffer_size_first.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_uframes_t)]

# /usr/include/alsa/pcm.h:668
snd_pcm_hw_params_set_buffer_size_last = _lib.snd_pcm_hw_params_set_buffer_size_last
snd_pcm_hw_params_set_buffer_size_last.restype = c_int
snd_pcm_hw_params_set_buffer_size_last.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_uframes_t)]

# /usr/include/alsa/pcm.h:670
snd_pcm_hw_params_get_tick_time = _lib.snd_pcm_hw_params_get_tick_time
snd_pcm_hw_params_get_tick_time.restype = c_int
snd_pcm_hw_params_get_tick_time.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:671
snd_pcm_hw_params_get_tick_time_min = _lib.snd_pcm_hw_params_get_tick_time_min
snd_pcm_hw_params_get_tick_time_min.restype = c_int
snd_pcm_hw_params_get_tick_time_min.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:672
snd_pcm_hw_params_get_tick_time_max = _lib.snd_pcm_hw_params_get_tick_time_max
snd_pcm_hw_params_get_tick_time_max.restype = c_int
snd_pcm_hw_params_get_tick_time_max.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:673
snd_pcm_hw_params_test_tick_time = _lib.snd_pcm_hw_params_test_tick_time
snd_pcm_hw_params_test_tick_time.restype = c_int
snd_pcm_hw_params_test_tick_time.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), c_uint, c_int]

# /usr/include/alsa/pcm.h:674
snd_pcm_hw_params_set_tick_time = _lib.snd_pcm_hw_params_set_tick_time
snd_pcm_hw_params_set_tick_time.restype = c_int
snd_pcm_hw_params_set_tick_time.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), c_uint, c_int]

# /usr/include/alsa/pcm.h:675
snd_pcm_hw_params_set_tick_time_min = _lib.snd_pcm_hw_params_set_tick_time_min
snd_pcm_hw_params_set_tick_time_min.restype = c_int
snd_pcm_hw_params_set_tick_time_min.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:676
snd_pcm_hw_params_set_tick_time_max = _lib.snd_pcm_hw_params_set_tick_time_max
snd_pcm_hw_params_set_tick_time_max.restype = c_int
snd_pcm_hw_params_set_tick_time_max.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:677
snd_pcm_hw_params_set_tick_time_minmax = _lib.snd_pcm_hw_params_set_tick_time_minmax
snd_pcm_hw_params_set_tick_time_minmax.restype = c_int
snd_pcm_hw_params_set_tick_time_minmax.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:678
snd_pcm_hw_params_set_tick_time_near = _lib.snd_pcm_hw_params_set_tick_time_near
snd_pcm_hw_params_set_tick_time_near.restype = c_int
snd_pcm_hw_params_set_tick_time_near.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:679
snd_pcm_hw_params_set_tick_time_first = _lib.snd_pcm_hw_params_set_tick_time_first
snd_pcm_hw_params_set_tick_time_first.restype = c_int
snd_pcm_hw_params_set_tick_time_first.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:680
snd_pcm_hw_params_set_tick_time_last = _lib.snd_pcm_hw_params_set_tick_time_last
snd_pcm_hw_params_set_tick_time_last.restype = c_int
snd_pcm_hw_params_set_tick_time_last.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_hw_params_t), POINTER(c_uint), POINTER(c_int)]

# /usr/include/alsa/pcm.h:684
snd_pcm_hw_params_get_min_align = _lib.snd_pcm_hw_params_get_min_align
snd_pcm_hw_params_get_min_align.restype = c_int
snd_pcm_hw_params_get_min_align.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(snd_pcm_uframes_t)]

# /usr/include/alsa/pcm.h:695
snd_pcm_sw_params_sizeof = _lib.snd_pcm_sw_params_sizeof
snd_pcm_sw_params_sizeof.restype = c_size_t
snd_pcm_sw_params_sizeof.argtypes = []

# /usr/include/alsa/pcm.h:701
snd_pcm_sw_params_malloc = _lib.snd_pcm_sw_params_malloc
snd_pcm_sw_params_malloc.restype = c_int
snd_pcm_sw_params_malloc.argtypes = [POINTER(POINTER(snd_pcm_sw_params_t))]

# /usr/include/alsa/pcm.h:702
snd_pcm_sw_params_free = _lib.snd_pcm_sw_params_free
snd_pcm_sw_params_free.restype = None
snd_pcm_sw_params_free.argtypes = [POINTER(snd_pcm_sw_params_t)]

# /usr/include/alsa/pcm.h:703
snd_pcm_sw_params_copy = _lib.snd_pcm_sw_params_copy
snd_pcm_sw_params_copy.restype = None
snd_pcm_sw_params_copy.argtypes = [POINTER(snd_pcm_sw_params_t), POINTER(snd_pcm_sw_params_t)]

# /usr/include/alsa/pcm.h:704
snd_pcm_sw_params_get_boundary = _lib.snd_pcm_sw_params_get_boundary
snd_pcm_sw_params_get_boundary.restype = c_int
snd_pcm_sw_params_get_boundary.argtypes = [POINTER(snd_pcm_sw_params_t), POINTER(snd_pcm_uframes_t)]

# /usr/include/alsa/pcm.h:708
snd_pcm_sw_params_set_tstamp_mode = _lib.snd_pcm_sw_params_set_tstamp_mode
snd_pcm_sw_params_set_tstamp_mode.restype = c_int
snd_pcm_sw_params_set_tstamp_mode.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_sw_params_t), snd_pcm_tstamp_t]

# /usr/include/alsa/pcm.h:709
snd_pcm_sw_params_get_tstamp_mode = _lib.snd_pcm_sw_params_get_tstamp_mode
snd_pcm_sw_params_get_tstamp_mode.restype = c_int
snd_pcm_sw_params_get_tstamp_mode.argtypes = [POINTER(snd_pcm_sw_params_t), POINTER(snd_pcm_tstamp_t)]

# /usr/include/alsa/pcm.h:710
snd_pcm_sw_params_set_sleep_min = _lib.snd_pcm_sw_params_set_sleep_min
snd_pcm_sw_params_set_sleep_min.restype = c_int
snd_pcm_sw_params_set_sleep_min.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_sw_params_t), c_uint]

# /usr/include/alsa/pcm.h:711
snd_pcm_sw_params_get_sleep_min = _lib.snd_pcm_sw_params_get_sleep_min
snd_pcm_sw_params_get_sleep_min.restype = c_int
snd_pcm_sw_params_get_sleep_min.argtypes = [POINTER(snd_pcm_sw_params_t), POINTER(c_uint)]

# /usr/include/alsa/pcm.h:712
snd_pcm_sw_params_set_avail_min = _lib.snd_pcm_sw_params_set_avail_min
snd_pcm_sw_params_set_avail_min.restype = c_int
snd_pcm_sw_params_set_avail_min.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_sw_params_t), snd_pcm_uframes_t]

# /usr/include/alsa/pcm.h:713
snd_pcm_sw_params_get_avail_min = _lib.snd_pcm_sw_params_get_avail_min
snd_pcm_sw_params_get_avail_min.restype = c_int
snd_pcm_sw_params_get_avail_min.argtypes = [POINTER(snd_pcm_sw_params_t), POINTER(snd_pcm_uframes_t)]

# /usr/include/alsa/pcm.h:714
snd_pcm_sw_params_set_xfer_align = _lib.snd_pcm_sw_params_set_xfer_align
snd_pcm_sw_params_set_xfer_align.restype = c_int
snd_pcm_sw_params_set_xfer_align.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_sw_params_t), snd_pcm_uframes_t]

# /usr/include/alsa/pcm.h:715
snd_pcm_sw_params_get_xfer_align = _lib.snd_pcm_sw_params_get_xfer_align
snd_pcm_sw_params_get_xfer_align.restype = c_int
snd_pcm_sw_params_get_xfer_align.argtypes = [POINTER(snd_pcm_sw_params_t), POINTER(snd_pcm_uframes_t)]

# /usr/include/alsa/pcm.h:716
snd_pcm_sw_params_set_start_threshold = _lib.snd_pcm_sw_params_set_start_threshold
snd_pcm_sw_params_set_start_threshold.restype = c_int
snd_pcm_sw_params_set_start_threshold.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_sw_params_t), snd_pcm_uframes_t]

# /usr/include/alsa/pcm.h:717
snd_pcm_sw_params_get_start_threshold = _lib.snd_pcm_sw_params_get_start_threshold
snd_pcm_sw_params_get_start_threshold.restype = c_int
snd_pcm_sw_params_get_start_threshold.argtypes = [POINTER(snd_pcm_sw_params_t), POINTER(snd_pcm_uframes_t)]

# /usr/include/alsa/pcm.h:718
snd_pcm_sw_params_set_stop_threshold = _lib.snd_pcm_sw_params_set_stop_threshold
snd_pcm_sw_params_set_stop_threshold.restype = c_int
snd_pcm_sw_params_set_stop_threshold.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_sw_params_t), snd_pcm_uframes_t]

# /usr/include/alsa/pcm.h:719
snd_pcm_sw_params_get_stop_threshold = _lib.snd_pcm_sw_params_get_stop_threshold
snd_pcm_sw_params_get_stop_threshold.restype = c_int
snd_pcm_sw_params_get_stop_threshold.argtypes = [POINTER(snd_pcm_sw_params_t), POINTER(snd_pcm_uframes_t)]

# /usr/include/alsa/pcm.h:720
snd_pcm_sw_params_set_silence_threshold = _lib.snd_pcm_sw_params_set_silence_threshold
snd_pcm_sw_params_set_silence_threshold.restype = c_int
snd_pcm_sw_params_set_silence_threshold.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_sw_params_t), snd_pcm_uframes_t]

# /usr/include/alsa/pcm.h:721
snd_pcm_sw_params_get_silence_threshold = _lib.snd_pcm_sw_params_get_silence_threshold
snd_pcm_sw_params_get_silence_threshold.restype = c_int
snd_pcm_sw_params_get_silence_threshold.argtypes = [POINTER(snd_pcm_sw_params_t), POINTER(snd_pcm_uframes_t)]

# /usr/include/alsa/pcm.h:722
snd_pcm_sw_params_set_silence_size = _lib.snd_pcm_sw_params_set_silence_size
snd_pcm_sw_params_set_silence_size.restype = c_int
snd_pcm_sw_params_set_silence_size.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_sw_params_t), snd_pcm_uframes_t]

# /usr/include/alsa/pcm.h:723
snd_pcm_sw_params_get_silence_size = _lib.snd_pcm_sw_params_get_silence_size
snd_pcm_sw_params_get_silence_size.restype = c_int
snd_pcm_sw_params_get_silence_size.argtypes = [POINTER(snd_pcm_sw_params_t), POINTER(snd_pcm_uframes_t)]

# /usr/include/alsa/pcm.h:743
snd_pcm_access_mask_sizeof = _lib.snd_pcm_access_mask_sizeof
snd_pcm_access_mask_sizeof.restype = c_size_t
snd_pcm_access_mask_sizeof.argtypes = []

# /usr/include/alsa/pcm.h:749
snd_pcm_access_mask_malloc = _lib.snd_pcm_access_mask_malloc
snd_pcm_access_mask_malloc.restype = c_int
snd_pcm_access_mask_malloc.argtypes = [POINTER(POINTER(snd_pcm_access_mask_t))]

# /usr/include/alsa/pcm.h:750
snd_pcm_access_mask_free = _lib.snd_pcm_access_mask_free
snd_pcm_access_mask_free.restype = None
snd_pcm_access_mask_free.argtypes = [POINTER(snd_pcm_access_mask_t)]

# /usr/include/alsa/pcm.h:751
snd_pcm_access_mask_copy = _lib.snd_pcm_access_mask_copy
snd_pcm_access_mask_copy.restype = None
snd_pcm_access_mask_copy.argtypes = [POINTER(snd_pcm_access_mask_t), POINTER(snd_pcm_access_mask_t)]

# /usr/include/alsa/pcm.h:752
snd_pcm_access_mask_none = _lib.snd_pcm_access_mask_none
snd_pcm_access_mask_none.restype = None
snd_pcm_access_mask_none.argtypes = [POINTER(snd_pcm_access_mask_t)]

# /usr/include/alsa/pcm.h:753
snd_pcm_access_mask_any = _lib.snd_pcm_access_mask_any
snd_pcm_access_mask_any.restype = None
snd_pcm_access_mask_any.argtypes = [POINTER(snd_pcm_access_mask_t)]

# /usr/include/alsa/pcm.h:754
snd_pcm_access_mask_test = _lib.snd_pcm_access_mask_test
snd_pcm_access_mask_test.restype = c_int
snd_pcm_access_mask_test.argtypes = [POINTER(snd_pcm_access_mask_t), snd_pcm_access_t]

# /usr/include/alsa/pcm.h:755
snd_pcm_access_mask_empty = _lib.snd_pcm_access_mask_empty
snd_pcm_access_mask_empty.restype = c_int
snd_pcm_access_mask_empty.argtypes = [POINTER(snd_pcm_access_mask_t)]

# /usr/include/alsa/pcm.h:756
snd_pcm_access_mask_set = _lib.snd_pcm_access_mask_set
snd_pcm_access_mask_set.restype = None
snd_pcm_access_mask_set.argtypes = [POINTER(snd_pcm_access_mask_t), snd_pcm_access_t]

# /usr/include/alsa/pcm.h:757
snd_pcm_access_mask_reset = _lib.snd_pcm_access_mask_reset
snd_pcm_access_mask_reset.restype = None
snd_pcm_access_mask_reset.argtypes = [POINTER(snd_pcm_access_mask_t), snd_pcm_access_t]

# /usr/include/alsa/pcm.h:768
snd_pcm_format_mask_sizeof = _lib.snd_pcm_format_mask_sizeof
snd_pcm_format_mask_sizeof.restype = c_size_t
snd_pcm_format_mask_sizeof.argtypes = []

# /usr/include/alsa/pcm.h:774
snd_pcm_format_mask_malloc = _lib.snd_pcm_format_mask_malloc
snd_pcm_format_mask_malloc.restype = c_int
snd_pcm_format_mask_malloc.argtypes = [POINTER(POINTER(snd_pcm_format_mask_t))]

# /usr/include/alsa/pcm.h:775
snd_pcm_format_mask_free = _lib.snd_pcm_format_mask_free
snd_pcm_format_mask_free.restype = None
snd_pcm_format_mask_free.argtypes = [POINTER(snd_pcm_format_mask_t)]

# /usr/include/alsa/pcm.h:776
snd_pcm_format_mask_copy = _lib.snd_pcm_format_mask_copy
snd_pcm_format_mask_copy.restype = None
snd_pcm_format_mask_copy.argtypes = [POINTER(snd_pcm_format_mask_t), POINTER(snd_pcm_format_mask_t)]

# /usr/include/alsa/pcm.h:777
snd_pcm_format_mask_none = _lib.snd_pcm_format_mask_none
snd_pcm_format_mask_none.restype = None
snd_pcm_format_mask_none.argtypes = [POINTER(snd_pcm_format_mask_t)]

# /usr/include/alsa/pcm.h:778
snd_pcm_format_mask_any = _lib.snd_pcm_format_mask_any
snd_pcm_format_mask_any.restype = None
snd_pcm_format_mask_any.argtypes = [POINTER(snd_pcm_format_mask_t)]

# /usr/include/alsa/pcm.h:779
snd_pcm_format_mask_test = _lib.snd_pcm_format_mask_test
snd_pcm_format_mask_test.restype = c_int
snd_pcm_format_mask_test.argtypes = [POINTER(snd_pcm_format_mask_t), snd_pcm_format_t]

# /usr/include/alsa/pcm.h:780
snd_pcm_format_mask_empty = _lib.snd_pcm_format_mask_empty
snd_pcm_format_mask_empty.restype = c_int
snd_pcm_format_mask_empty.argtypes = [POINTER(snd_pcm_format_mask_t)]

# /usr/include/alsa/pcm.h:781
snd_pcm_format_mask_set = _lib.snd_pcm_format_mask_set
snd_pcm_format_mask_set.restype = None
snd_pcm_format_mask_set.argtypes = [POINTER(snd_pcm_format_mask_t), snd_pcm_format_t]

# /usr/include/alsa/pcm.h:782
snd_pcm_format_mask_reset = _lib.snd_pcm_format_mask_reset
snd_pcm_format_mask_reset.restype = None
snd_pcm_format_mask_reset.argtypes = [POINTER(snd_pcm_format_mask_t), snd_pcm_format_t]

# /usr/include/alsa/pcm.h:793
snd_pcm_subformat_mask_sizeof = _lib.snd_pcm_subformat_mask_sizeof
snd_pcm_subformat_mask_sizeof.restype = c_size_t
snd_pcm_subformat_mask_sizeof.argtypes = []

# /usr/include/alsa/pcm.h:799
snd_pcm_subformat_mask_malloc = _lib.snd_pcm_subformat_mask_malloc
snd_pcm_subformat_mask_malloc.restype = c_int
snd_pcm_subformat_mask_malloc.argtypes = [POINTER(POINTER(snd_pcm_subformat_mask_t))]

# /usr/include/alsa/pcm.h:800
snd_pcm_subformat_mask_free = _lib.snd_pcm_subformat_mask_free
snd_pcm_subformat_mask_free.restype = None
snd_pcm_subformat_mask_free.argtypes = [POINTER(snd_pcm_subformat_mask_t)]

# /usr/include/alsa/pcm.h:801
snd_pcm_subformat_mask_copy = _lib.snd_pcm_subformat_mask_copy
snd_pcm_subformat_mask_copy.restype = None
snd_pcm_subformat_mask_copy.argtypes = [POINTER(snd_pcm_subformat_mask_t), POINTER(snd_pcm_subformat_mask_t)]

# /usr/include/alsa/pcm.h:802
snd_pcm_subformat_mask_none = _lib.snd_pcm_subformat_mask_none
snd_pcm_subformat_mask_none.restype = None
snd_pcm_subformat_mask_none.argtypes = [POINTER(snd_pcm_subformat_mask_t)]

# /usr/include/alsa/pcm.h:803
snd_pcm_subformat_mask_any = _lib.snd_pcm_subformat_mask_any
snd_pcm_subformat_mask_any.restype = None
snd_pcm_subformat_mask_any.argtypes = [POINTER(snd_pcm_subformat_mask_t)]

# /usr/include/alsa/pcm.h:804
snd_pcm_subformat_mask_test = _lib.snd_pcm_subformat_mask_test
snd_pcm_subformat_mask_test.restype = c_int
snd_pcm_subformat_mask_test.argtypes = [POINTER(snd_pcm_subformat_mask_t), snd_pcm_subformat_t]

# /usr/include/alsa/pcm.h:805
snd_pcm_subformat_mask_empty = _lib.snd_pcm_subformat_mask_empty
snd_pcm_subformat_mask_empty.restype = c_int
snd_pcm_subformat_mask_empty.argtypes = [POINTER(snd_pcm_subformat_mask_t)]

# /usr/include/alsa/pcm.h:806
snd_pcm_subformat_mask_set = _lib.snd_pcm_subformat_mask_set
snd_pcm_subformat_mask_set.restype = None
snd_pcm_subformat_mask_set.argtypes = [POINTER(snd_pcm_subformat_mask_t), snd_pcm_subformat_t]

# /usr/include/alsa/pcm.h:807
snd_pcm_subformat_mask_reset = _lib.snd_pcm_subformat_mask_reset
snd_pcm_subformat_mask_reset.restype = None
snd_pcm_subformat_mask_reset.argtypes = [POINTER(snd_pcm_subformat_mask_t), snd_pcm_subformat_t]

# /usr/include/alsa/pcm.h:818
snd_pcm_status_sizeof = _lib.snd_pcm_status_sizeof
snd_pcm_status_sizeof.restype = c_size_t
snd_pcm_status_sizeof.argtypes = []

# /usr/include/alsa/pcm.h:824
snd_pcm_status_malloc = _lib.snd_pcm_status_malloc
snd_pcm_status_malloc.restype = c_int
snd_pcm_status_malloc.argtypes = [POINTER(POINTER(snd_pcm_status_t))]

# /usr/include/alsa/pcm.h:825
snd_pcm_status_free = _lib.snd_pcm_status_free
snd_pcm_status_free.restype = None
snd_pcm_status_free.argtypes = [POINTER(snd_pcm_status_t)]

# /usr/include/alsa/pcm.h:826
snd_pcm_status_copy = _lib.snd_pcm_status_copy
snd_pcm_status_copy.restype = None
snd_pcm_status_copy.argtypes = [POINTER(snd_pcm_status_t), POINTER(snd_pcm_status_t)]

# /usr/include/alsa/pcm.h:827
snd_pcm_status_get_state = _lib.snd_pcm_status_get_state
snd_pcm_status_get_state.restype = snd_pcm_state_t
snd_pcm_status_get_state.argtypes = [POINTER(snd_pcm_status_t)]

# /usr/include/alsa/pcm.h:828
snd_pcm_status_get_trigger_tstamp = _lib.snd_pcm_status_get_trigger_tstamp
snd_pcm_status_get_trigger_tstamp.restype = None
snd_pcm_status_get_trigger_tstamp.argtypes = [POINTER(snd_pcm_status_t), POINTER(snd_timestamp_t)]

# /usr/include/alsa/pcm.h:829
snd_pcm_status_get_trigger_htstamp = _lib.snd_pcm_status_get_trigger_htstamp
snd_pcm_status_get_trigger_htstamp.restype = None
snd_pcm_status_get_trigger_htstamp.argtypes = [POINTER(snd_pcm_status_t), POINTER(snd_htimestamp_t)]

# /usr/include/alsa/pcm.h:830
snd_pcm_status_get_tstamp = _lib.snd_pcm_status_get_tstamp
snd_pcm_status_get_tstamp.restype = None
snd_pcm_status_get_tstamp.argtypes = [POINTER(snd_pcm_status_t), POINTER(snd_timestamp_t)]

# /usr/include/alsa/pcm.h:831
snd_pcm_status_get_htstamp = _lib.snd_pcm_status_get_htstamp
snd_pcm_status_get_htstamp.restype = None
snd_pcm_status_get_htstamp.argtypes = [POINTER(snd_pcm_status_t), POINTER(snd_htimestamp_t)]

# /usr/include/alsa/pcm.h:832
snd_pcm_status_get_delay = _lib.snd_pcm_status_get_delay
snd_pcm_status_get_delay.restype = snd_pcm_sframes_t
snd_pcm_status_get_delay.argtypes = [POINTER(snd_pcm_status_t)]

# /usr/include/alsa/pcm.h:833
snd_pcm_status_get_avail = _lib.snd_pcm_status_get_avail
snd_pcm_status_get_avail.restype = snd_pcm_uframes_t
snd_pcm_status_get_avail.argtypes = [POINTER(snd_pcm_status_t)]

# /usr/include/alsa/pcm.h:834
snd_pcm_status_get_avail_max = _lib.snd_pcm_status_get_avail_max
snd_pcm_status_get_avail_max.restype = snd_pcm_uframes_t
snd_pcm_status_get_avail_max.argtypes = [POINTER(snd_pcm_status_t)]

# /usr/include/alsa/pcm.h:835
snd_pcm_status_get_overrange = _lib.snd_pcm_status_get_overrange
snd_pcm_status_get_overrange.restype = snd_pcm_uframes_t
snd_pcm_status_get_overrange.argtypes = [POINTER(snd_pcm_status_t)]

# /usr/include/alsa/pcm.h:846
snd_pcm_type_name = _lib.snd_pcm_type_name
snd_pcm_type_name.restype = c_char_p
snd_pcm_type_name.argtypes = [snd_pcm_type_t]

# /usr/include/alsa/pcm.h:847
snd_pcm_stream_name = _lib.snd_pcm_stream_name
snd_pcm_stream_name.restype = c_char_p
snd_pcm_stream_name.argtypes = [snd_pcm_stream_t]

# /usr/include/alsa/pcm.h:848
snd_pcm_access_name = _lib.snd_pcm_access_name
snd_pcm_access_name.restype = c_char_p
snd_pcm_access_name.argtypes = [snd_pcm_access_t]

# /usr/include/alsa/pcm.h:849
snd_pcm_format_name = _lib.snd_pcm_format_name
snd_pcm_format_name.restype = c_char_p
snd_pcm_format_name.argtypes = [snd_pcm_format_t]

# /usr/include/alsa/pcm.h:850
snd_pcm_format_description = _lib.snd_pcm_format_description
snd_pcm_format_description.restype = c_char_p
snd_pcm_format_description.argtypes = [snd_pcm_format_t]

# /usr/include/alsa/pcm.h:851
snd_pcm_subformat_name = _lib.snd_pcm_subformat_name
snd_pcm_subformat_name.restype = c_char_p
snd_pcm_subformat_name.argtypes = [snd_pcm_subformat_t]

# /usr/include/alsa/pcm.h:852
snd_pcm_subformat_description = _lib.snd_pcm_subformat_description
snd_pcm_subformat_description.restype = c_char_p
snd_pcm_subformat_description.argtypes = [snd_pcm_subformat_t]

# /usr/include/alsa/pcm.h:853
snd_pcm_format_value = _lib.snd_pcm_format_value
snd_pcm_format_value.restype = snd_pcm_format_t
snd_pcm_format_value.argtypes = [c_char_p]

# /usr/include/alsa/pcm.h:854
snd_pcm_tstamp_mode_name = _lib.snd_pcm_tstamp_mode_name
snd_pcm_tstamp_mode_name.restype = c_char_p
snd_pcm_tstamp_mode_name.argtypes = [snd_pcm_tstamp_t]

# /usr/include/alsa/pcm.h:855
snd_pcm_state_name = _lib.snd_pcm_state_name
snd_pcm_state_name.restype = c_char_p
snd_pcm_state_name.argtypes = [snd_pcm_state_t]

# /usr/include/alsa/pcm.h:866
snd_pcm_dump = _lib.snd_pcm_dump
snd_pcm_dump.restype = c_int
snd_pcm_dump.argtypes = [POINTER(snd_pcm_t), POINTER(snd_output_t)]

# /usr/include/alsa/pcm.h:867
snd_pcm_dump_hw_setup = _lib.snd_pcm_dump_hw_setup
snd_pcm_dump_hw_setup.restype = c_int
snd_pcm_dump_hw_setup.argtypes = [POINTER(snd_pcm_t), POINTER(snd_output_t)]

# /usr/include/alsa/pcm.h:868
snd_pcm_dump_sw_setup = _lib.snd_pcm_dump_sw_setup
snd_pcm_dump_sw_setup.restype = c_int
snd_pcm_dump_sw_setup.argtypes = [POINTER(snd_pcm_t), POINTER(snd_output_t)]

# /usr/include/alsa/pcm.h:869
snd_pcm_dump_setup = _lib.snd_pcm_dump_setup
snd_pcm_dump_setup.restype = c_int
snd_pcm_dump_setup.argtypes = [POINTER(snd_pcm_t), POINTER(snd_output_t)]

# /usr/include/alsa/pcm.h:870
snd_pcm_hw_params_dump = _lib.snd_pcm_hw_params_dump
snd_pcm_hw_params_dump.restype = c_int
snd_pcm_hw_params_dump.argtypes = [POINTER(snd_pcm_hw_params_t), POINTER(snd_output_t)]

# /usr/include/alsa/pcm.h:871
snd_pcm_sw_params_dump = _lib.snd_pcm_sw_params_dump
snd_pcm_sw_params_dump.restype = c_int
snd_pcm_sw_params_dump.argtypes = [POINTER(snd_pcm_sw_params_t), POINTER(snd_output_t)]

# /usr/include/alsa/pcm.h:872
snd_pcm_status_dump = _lib.snd_pcm_status_dump
snd_pcm_status_dump.restype = c_int
snd_pcm_status_dump.argtypes = [POINTER(snd_pcm_status_t), POINTER(snd_output_t)]

# /usr/include/alsa/pcm.h:883
snd_pcm_mmap_begin = _lib.snd_pcm_mmap_begin
snd_pcm_mmap_begin.restype = c_int
snd_pcm_mmap_begin.argtypes = [POINTER(snd_pcm_t), POINTER(POINTER(snd_pcm_channel_area_t)), POINTER(snd_pcm_uframes_t), POINTER(snd_pcm_uframes_t)]

# /usr/include/alsa/pcm.h:887
snd_pcm_mmap_commit = _lib.snd_pcm_mmap_commit
snd_pcm_mmap_commit.restype = snd_pcm_sframes_t
snd_pcm_mmap_commit.argtypes = [POINTER(snd_pcm_t), snd_pcm_uframes_t, snd_pcm_uframes_t]

# /usr/include/alsa/pcm.h:890
snd_pcm_mmap_writei = _lib.snd_pcm_mmap_writei
snd_pcm_mmap_writei.restype = snd_pcm_sframes_t
snd_pcm_mmap_writei.argtypes = [POINTER(snd_pcm_t), POINTER(None), snd_pcm_uframes_t]

# /usr/include/alsa/pcm.h:891
snd_pcm_mmap_readi = _lib.snd_pcm_mmap_readi
snd_pcm_mmap_readi.restype = snd_pcm_sframes_t
snd_pcm_mmap_readi.argtypes = [POINTER(snd_pcm_t), POINTER(None), snd_pcm_uframes_t]

# /usr/include/alsa/pcm.h:892
snd_pcm_mmap_writen = _lib.snd_pcm_mmap_writen
snd_pcm_mmap_writen.restype = snd_pcm_sframes_t
snd_pcm_mmap_writen.argtypes = [POINTER(snd_pcm_t), POINTER(POINTER(None)), snd_pcm_uframes_t]

# /usr/include/alsa/pcm.h:893
snd_pcm_mmap_readn = _lib.snd_pcm_mmap_readn
snd_pcm_mmap_readn.restype = snd_pcm_sframes_t
snd_pcm_mmap_readn.argtypes = [POINTER(snd_pcm_t), POINTER(POINTER(None)), snd_pcm_uframes_t]

# /usr/include/alsa/pcm.h:904
snd_pcm_format_signed = _lib.snd_pcm_format_signed
snd_pcm_format_signed.restype = c_int
snd_pcm_format_signed.argtypes = [snd_pcm_format_t]

# /usr/include/alsa/pcm.h:905
snd_pcm_format_unsigned = _lib.snd_pcm_format_unsigned
snd_pcm_format_unsigned.restype = c_int
snd_pcm_format_unsigned.argtypes = [snd_pcm_format_t]

# /usr/include/alsa/pcm.h:906
snd_pcm_format_linear = _lib.snd_pcm_format_linear
snd_pcm_format_linear.restype = c_int
snd_pcm_format_linear.argtypes = [snd_pcm_format_t]

# /usr/include/alsa/pcm.h:907
snd_pcm_format_float = _lib.snd_pcm_format_float
snd_pcm_format_float.restype = c_int
snd_pcm_format_float.argtypes = [snd_pcm_format_t]

# /usr/include/alsa/pcm.h:908
snd_pcm_format_little_endian = _lib.snd_pcm_format_little_endian
snd_pcm_format_little_endian.restype = c_int
snd_pcm_format_little_endian.argtypes = [snd_pcm_format_t]

# /usr/include/alsa/pcm.h:909
snd_pcm_format_big_endian = _lib.snd_pcm_format_big_endian
snd_pcm_format_big_endian.restype = c_int
snd_pcm_format_big_endian.argtypes = [snd_pcm_format_t]

# /usr/include/alsa/pcm.h:910
snd_pcm_format_cpu_endian = _lib.snd_pcm_format_cpu_endian
snd_pcm_format_cpu_endian.restype = c_int
snd_pcm_format_cpu_endian.argtypes = [snd_pcm_format_t]

# /usr/include/alsa/pcm.h:911
snd_pcm_format_width = _lib.snd_pcm_format_width
snd_pcm_format_width.restype = c_int
snd_pcm_format_width.argtypes = [snd_pcm_format_t]

# /usr/include/alsa/pcm.h:912
snd_pcm_format_physical_width = _lib.snd_pcm_format_physical_width
snd_pcm_format_physical_width.restype = c_int
snd_pcm_format_physical_width.argtypes = [snd_pcm_format_t]

# /usr/include/alsa/pcm.h:913
snd_pcm_build_linear_format = _lib.snd_pcm_build_linear_format
snd_pcm_build_linear_format.restype = snd_pcm_format_t
snd_pcm_build_linear_format.argtypes = [c_int, c_int, c_int, c_int]

# /usr/include/alsa/pcm.h:914
snd_pcm_format_size = _lib.snd_pcm_format_size
snd_pcm_format_size.restype = ssize_t
snd_pcm_format_size.argtypes = [snd_pcm_format_t, c_size_t]

u_int8_t = c_ubyte 	# /usr/include/gentoo-multilib/amd64/sys/types.h:174
# /usr/include/alsa/pcm.h:915
snd_pcm_format_silence = _lib.snd_pcm_format_silence
snd_pcm_format_silence.restype = u_int8_t
snd_pcm_format_silence.argtypes = [snd_pcm_format_t]

u_int16_t = c_uint 	# /usr/include/gentoo-multilib/amd64/sys/types.h:175
# /usr/include/alsa/pcm.h:916
snd_pcm_format_silence_16 = _lib.snd_pcm_format_silence_16
snd_pcm_format_silence_16.restype = u_int16_t
snd_pcm_format_silence_16.argtypes = [snd_pcm_format_t]

u_int32_t = c_uint 	# /usr/include/gentoo-multilib/amd64/sys/types.h:176
# /usr/include/alsa/pcm.h:917
snd_pcm_format_silence_32 = _lib.snd_pcm_format_silence_32
snd_pcm_format_silence_32.restype = u_int32_t
snd_pcm_format_silence_32.argtypes = [snd_pcm_format_t]

u_int64_t = c_ulong 	# /usr/include/gentoo-multilib/amd64/sys/types.h:178
# /usr/include/alsa/pcm.h:918
snd_pcm_format_silence_64 = _lib.snd_pcm_format_silence_64
snd_pcm_format_silence_64.restype = u_int64_t
snd_pcm_format_silence_64.argtypes = [snd_pcm_format_t]

# /usr/include/alsa/pcm.h:919
snd_pcm_format_set_silence = _lib.snd_pcm_format_set_silence
snd_pcm_format_set_silence.restype = c_int
snd_pcm_format_set_silence.argtypes = [snd_pcm_format_t, POINTER(None), c_uint]

# /usr/include/alsa/pcm.h:921
snd_pcm_bytes_to_frames = _lib.snd_pcm_bytes_to_frames
snd_pcm_bytes_to_frames.restype = snd_pcm_sframes_t
snd_pcm_bytes_to_frames.argtypes = [POINTER(snd_pcm_t), ssize_t]

# /usr/include/alsa/pcm.h:922
snd_pcm_frames_to_bytes = _lib.snd_pcm_frames_to_bytes
snd_pcm_frames_to_bytes.restype = ssize_t
snd_pcm_frames_to_bytes.argtypes = [POINTER(snd_pcm_t), snd_pcm_sframes_t]

# /usr/include/alsa/pcm.h:923
snd_pcm_bytes_to_samples = _lib.snd_pcm_bytes_to_samples
snd_pcm_bytes_to_samples.restype = c_long
snd_pcm_bytes_to_samples.argtypes = [POINTER(snd_pcm_t), ssize_t]

# /usr/include/alsa/pcm.h:924
snd_pcm_samples_to_bytes = _lib.snd_pcm_samples_to_bytes
snd_pcm_samples_to_bytes.restype = ssize_t
snd_pcm_samples_to_bytes.argtypes = [POINTER(snd_pcm_t), c_long]

# /usr/include/alsa/pcm.h:926
snd_pcm_area_silence = _lib.snd_pcm_area_silence
snd_pcm_area_silence.restype = c_int
snd_pcm_area_silence.argtypes = [POINTER(snd_pcm_channel_area_t), snd_pcm_uframes_t, c_uint, snd_pcm_format_t]

# /usr/include/alsa/pcm.h:928
snd_pcm_areas_silence = _lib.snd_pcm_areas_silence
snd_pcm_areas_silence.restype = c_int
snd_pcm_areas_silence.argtypes = [POINTER(snd_pcm_channel_area_t), snd_pcm_uframes_t, c_uint, snd_pcm_uframes_t, snd_pcm_format_t]

# /usr/include/alsa/pcm.h:930
snd_pcm_area_copy = _lib.snd_pcm_area_copy
snd_pcm_area_copy.restype = c_int
snd_pcm_area_copy.argtypes = [POINTER(snd_pcm_channel_area_t), snd_pcm_uframes_t, POINTER(snd_pcm_channel_area_t), snd_pcm_uframes_t, c_uint, snd_pcm_format_t]

# /usr/include/alsa/pcm.h:933
snd_pcm_areas_copy = _lib.snd_pcm_areas_copy
snd_pcm_areas_copy.restype = c_int
snd_pcm_areas_copy.argtypes = [POINTER(snd_pcm_channel_area_t), snd_pcm_uframes_t, POINTER(snd_pcm_channel_area_t), snd_pcm_uframes_t, c_uint, snd_pcm_uframes_t, snd_pcm_format_t]

enum__snd_pcm_hook_type = c_int
SND_PCM_HOOK_TYPE_HW_PARAMS = 0
SND_PCM_HOOK_TYPE_HW_FREE = 1
SND_PCM_HOOK_TYPE_CLOSE = 2
SND_PCM_HOOK_TYPE_LAST = 0
snd_pcm_hook_type_t = enum__snd_pcm_hook_type 	# /usr/include/alsa/pcm.h:952
class struct__snd_pcm_hook(Structure):
    __slots__ = [
    ]
struct__snd_pcm_hook._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_pcm_hook(Structure):
    __slots__ = [
    ]
struct__snd_pcm_hook._fields_ = [
    ('_opaque_struct', c_int)
]

snd_pcm_hook_t = struct__snd_pcm_hook 	# /usr/include/alsa/pcm.h:955
snd_pcm_hook_func_t = CFUNCTYPE(c_int, POINTER(snd_pcm_hook_t)) 	# /usr/include/alsa/pcm.h:957
# /usr/include/alsa/pcm.h:958
snd_pcm_hook_get_pcm = _lib.snd_pcm_hook_get_pcm
snd_pcm_hook_get_pcm.restype = POINTER(snd_pcm_t)
snd_pcm_hook_get_pcm.argtypes = [POINTER(snd_pcm_hook_t)]

# /usr/include/alsa/pcm.h:959
snd_pcm_hook_get_private = _lib.snd_pcm_hook_get_private
snd_pcm_hook_get_private.restype = POINTER(c_void)
snd_pcm_hook_get_private.argtypes = [POINTER(snd_pcm_hook_t)]

# /usr/include/alsa/pcm.h:960
snd_pcm_hook_set_private = _lib.snd_pcm_hook_set_private
snd_pcm_hook_set_private.restype = None
snd_pcm_hook_set_private.argtypes = [POINTER(snd_pcm_hook_t), POINTER(None)]

# /usr/include/alsa/pcm.h:961
snd_pcm_hook_add = _lib.snd_pcm_hook_add
snd_pcm_hook_add.restype = c_int
snd_pcm_hook_add.argtypes = [POINTER(POINTER(snd_pcm_hook_t)), POINTER(snd_pcm_t), snd_pcm_hook_type_t, snd_pcm_hook_func_t, POINTER(None)]

# /usr/include/alsa/pcm.h:964
snd_pcm_hook_remove = _lib.snd_pcm_hook_remove
snd_pcm_hook_remove.restype = c_int
snd_pcm_hook_remove.argtypes = [POINTER(snd_pcm_hook_t)]

class struct__snd_pcm_scope_ops(Structure):
    __slots__ = [
        'enable',
        'disable',
        'start',
        'stop',
        'update',
        'reset',
        'close',
    ]
struct__snd_pcm_scope_ops._fields_ = [
    ('enable', POINTER(CFUNCTYPE(c_int, POINTER(snd_pcm_scope_t)))),
    ('disable', POINTER(CFUNCTYPE(None, POINTER(snd_pcm_scope_t)))),
    ('start', POINTER(CFUNCTYPE(None, POINTER(snd_pcm_scope_t)))),
    ('stop', POINTER(CFUNCTYPE(None, POINTER(snd_pcm_scope_t)))),
    ('update', POINTER(CFUNCTYPE(None, POINTER(snd_pcm_scope_t)))),
    ('reset', POINTER(CFUNCTYPE(None, POINTER(snd_pcm_scope_t)))),
    ('close', POINTER(CFUNCTYPE(None, POINTER(snd_pcm_scope_t)))),
]

snd_pcm_scope_ops_t = struct__snd_pcm_scope_ops 	# /usr/include/alsa/pcm.h:1005
# /usr/include/alsa/pcm.h:1007
snd_pcm_meter_get_bufsize = _lib.snd_pcm_meter_get_bufsize
snd_pcm_meter_get_bufsize.restype = snd_pcm_uframes_t
snd_pcm_meter_get_bufsize.argtypes = [POINTER(snd_pcm_t)]

# /usr/include/alsa/pcm.h:1008
snd_pcm_meter_get_channels = _lib.snd_pcm_meter_get_channels
snd_pcm_meter_get_channels.restype = c_uint
snd_pcm_meter_get_channels.argtypes = [POINTER(snd_pcm_t)]

# /usr/include/alsa/pcm.h:1009
snd_pcm_meter_get_rate = _lib.snd_pcm_meter_get_rate
snd_pcm_meter_get_rate.restype = c_uint
snd_pcm_meter_get_rate.argtypes = [POINTER(snd_pcm_t)]

# /usr/include/alsa/pcm.h:1010
snd_pcm_meter_get_now = _lib.snd_pcm_meter_get_now
snd_pcm_meter_get_now.restype = snd_pcm_uframes_t
snd_pcm_meter_get_now.argtypes = [POINTER(snd_pcm_t)]

# /usr/include/alsa/pcm.h:1011
snd_pcm_meter_get_boundary = _lib.snd_pcm_meter_get_boundary
snd_pcm_meter_get_boundary.restype = snd_pcm_uframes_t
snd_pcm_meter_get_boundary.argtypes = [POINTER(snd_pcm_t)]

# /usr/include/alsa/pcm.h:1012
snd_pcm_meter_add_scope = _lib.snd_pcm_meter_add_scope
snd_pcm_meter_add_scope.restype = c_int
snd_pcm_meter_add_scope.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_scope_t)]

# /usr/include/alsa/pcm.h:1013
snd_pcm_meter_search_scope = _lib.snd_pcm_meter_search_scope
snd_pcm_meter_search_scope.restype = POINTER(snd_pcm_scope_t)
snd_pcm_meter_search_scope.argtypes = [POINTER(snd_pcm_t), c_char_p]

# /usr/include/alsa/pcm.h:1014
snd_pcm_scope_malloc = _lib.snd_pcm_scope_malloc
snd_pcm_scope_malloc.restype = c_int
snd_pcm_scope_malloc.argtypes = [POINTER(POINTER(snd_pcm_scope_t))]

# /usr/include/alsa/pcm.h:1015
snd_pcm_scope_set_ops = _lib.snd_pcm_scope_set_ops
snd_pcm_scope_set_ops.restype = None
snd_pcm_scope_set_ops.argtypes = [POINTER(snd_pcm_scope_t), POINTER(snd_pcm_scope_ops_t)]

# /usr/include/alsa/pcm.h:1016
snd_pcm_scope_set_name = _lib.snd_pcm_scope_set_name
snd_pcm_scope_set_name.restype = None
snd_pcm_scope_set_name.argtypes = [POINTER(snd_pcm_scope_t), c_char_p]

# /usr/include/alsa/pcm.h:1017
snd_pcm_scope_get_name = _lib.snd_pcm_scope_get_name
snd_pcm_scope_get_name.restype = c_char_p
snd_pcm_scope_get_name.argtypes = [POINTER(snd_pcm_scope_t)]

# /usr/include/alsa/pcm.h:1018
snd_pcm_scope_get_callback_private = _lib.snd_pcm_scope_get_callback_private
snd_pcm_scope_get_callback_private.restype = POINTER(c_void)
snd_pcm_scope_get_callback_private.argtypes = [POINTER(snd_pcm_scope_t)]

# /usr/include/alsa/pcm.h:1019
snd_pcm_scope_set_callback_private = _lib.snd_pcm_scope_set_callback_private
snd_pcm_scope_set_callback_private.restype = None
snd_pcm_scope_set_callback_private.argtypes = [POINTER(snd_pcm_scope_t), POINTER(None)]

# /usr/include/alsa/pcm.h:1020
snd_pcm_scope_s16_open = _lib.snd_pcm_scope_s16_open
snd_pcm_scope_s16_open.restype = c_int
snd_pcm_scope_s16_open.argtypes = [POINTER(snd_pcm_t), c_char_p, POINTER(POINTER(snd_pcm_scope_t))]

# /usr/include/alsa/pcm.h:1022
snd_pcm_scope_s16_get_channel_buffer = _lib.snd_pcm_scope_s16_get_channel_buffer
snd_pcm_scope_s16_get_channel_buffer.restype = POINTER(c_int16)
snd_pcm_scope_s16_get_channel_buffer.argtypes = [POINTER(snd_pcm_scope_t), c_uint]

enum__snd_spcm_latency = c_int
SND_SPCM_LATENCY_STANDARD = 0
SND_SPCM_LATENCY_MEDIUM = 1
SND_SPCM_LATENCY_REALTIME = 2
snd_spcm_latency_t = enum__snd_spcm_latency 	# /usr/include/alsa/pcm.h:1045
enum__snd_spcm_xrun_type = c_int
SND_SPCM_XRUN_IGNORE = 0
SND_SPCM_XRUN_STOP = 1
snd_spcm_xrun_type_t = enum__snd_spcm_xrun_type 	# /usr/include/alsa/pcm.h:1053
enum__snd_spcm_duplex_type = c_int
SND_SPCM_DUPLEX_LIBERAL = 0
SND_SPCM_DUPLEX_PEDANTIC = 1
snd_spcm_duplex_type_t = enum__snd_spcm_duplex_type 	# /usr/include/alsa/pcm.h:1061
# /usr/include/alsa/pcm.h:1063
snd_spcm_init = _lib.snd_spcm_init
snd_spcm_init.restype = c_int
snd_spcm_init.argtypes = [POINTER(snd_pcm_t), c_uint, c_uint, snd_pcm_format_t, snd_pcm_subformat_t, snd_spcm_latency_t, snd_pcm_access_t, snd_spcm_xrun_type_t]

# /usr/include/alsa/pcm.h:1072
snd_spcm_init_duplex = _lib.snd_spcm_init_duplex
snd_spcm_init_duplex.restype = c_int
snd_spcm_init_duplex.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_t), c_uint, c_uint, snd_pcm_format_t, snd_pcm_subformat_t, snd_spcm_latency_t, snd_pcm_access_t, snd_spcm_xrun_type_t, snd_spcm_duplex_type_t]

# /usr/include/alsa/pcm.h:1083
snd_spcm_init_get_params = _lib.snd_spcm_init_get_params
snd_spcm_init_get_params.restype = c_int
snd_spcm_init_get_params.argtypes = [POINTER(snd_pcm_t), POINTER(c_uint), POINTER(snd_pcm_uframes_t), POINTER(snd_pcm_uframes_t)]

# /usr/include/alsa/pcm.h:1098
snd_pcm_start_mode_name = _lib.snd_pcm_start_mode_name
snd_pcm_start_mode_name.restype = c_char_p
snd_pcm_start_mode_name.argtypes = [snd_pcm_start_t]

# /usr/include/alsa/pcm.h:1099
snd_pcm_xrun_mode_name = _lib.snd_pcm_xrun_mode_name
snd_pcm_xrun_mode_name.restype = c_char_p
snd_pcm_xrun_mode_name.argtypes = [snd_pcm_xrun_t]

# /usr/include/alsa/pcm.h:1100
snd_pcm_sw_params_set_start_mode = _lib.snd_pcm_sw_params_set_start_mode
snd_pcm_sw_params_set_start_mode.restype = c_int
snd_pcm_sw_params_set_start_mode.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_sw_params_t), snd_pcm_start_t]

# /usr/include/alsa/pcm.h:1101
snd_pcm_sw_params_get_start_mode = _lib.snd_pcm_sw_params_get_start_mode
snd_pcm_sw_params_get_start_mode.restype = snd_pcm_start_t
snd_pcm_sw_params_get_start_mode.argtypes = [POINTER(snd_pcm_sw_params_t)]

# /usr/include/alsa/pcm.h:1102
snd_pcm_sw_params_set_xrun_mode = _lib.snd_pcm_sw_params_set_xrun_mode
snd_pcm_sw_params_set_xrun_mode.restype = c_int
snd_pcm_sw_params_set_xrun_mode.argtypes = [POINTER(snd_pcm_t), POINTER(snd_pcm_sw_params_t), snd_pcm_xrun_t]

# /usr/include/alsa/pcm.h:1103
snd_pcm_sw_params_get_xrun_mode = _lib.snd_pcm_sw_params_get_xrun_mode
snd_pcm_sw_params_get_xrun_mode.restype = snd_pcm_xrun_t
snd_pcm_sw_params_get_xrun_mode.argtypes = [POINTER(snd_pcm_sw_params_t)]

SND_RAWMIDI_DLSYM_VERSION = 0 	# /usr/include/alsa/rawmidi.h:42
class struct__snd_rawmidi_info(Structure):
    __slots__ = [
    ]
struct__snd_rawmidi_info._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_rawmidi_info(Structure):
    __slots__ = [
    ]
struct__snd_rawmidi_info._fields_ = [
    ('_opaque_struct', c_int)
]

snd_rawmidi_info_t = struct__snd_rawmidi_info 	# /usr/include/alsa/rawmidi.h:45
class struct__snd_rawmidi_params(Structure):
    __slots__ = [
    ]
struct__snd_rawmidi_params._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_rawmidi_params(Structure):
    __slots__ = [
    ]
struct__snd_rawmidi_params._fields_ = [
    ('_opaque_struct', c_int)
]

snd_rawmidi_params_t = struct__snd_rawmidi_params 	# /usr/include/alsa/rawmidi.h:47
class struct__snd_rawmidi_status(Structure):
    __slots__ = [
    ]
struct__snd_rawmidi_status._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_rawmidi_status(Structure):
    __slots__ = [
    ]
struct__snd_rawmidi_status._fields_ = [
    ('_opaque_struct', c_int)
]

snd_rawmidi_status_t = struct__snd_rawmidi_status 	# /usr/include/alsa/rawmidi.h:49
enum__snd_rawmidi_stream = c_int
SND_RAWMIDI_STREAM_OUTPUT = 0
SND_RAWMIDI_STREAM_INPUT = 1
SND_RAWMIDI_STREAM_LAST = 0
snd_rawmidi_stream_t = enum__snd_rawmidi_stream 	# /usr/include/alsa/rawmidi.h:58
SND_RAWMIDI_APPEND = 1 	# /usr/include/alsa/rawmidi.h:61
SND_RAWMIDI_NONBLOCK = 2 	# /usr/include/alsa/rawmidi.h:63
SND_RAWMIDI_SYNC = 4 	# /usr/include/alsa/rawmidi.h:65
class struct__snd_rawmidi(Structure):
    __slots__ = [
    ]
struct__snd_rawmidi._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_rawmidi(Structure):
    __slots__ = [
    ]
struct__snd_rawmidi._fields_ = [
    ('_opaque_struct', c_int)
]

snd_rawmidi_t = struct__snd_rawmidi 	# /usr/include/alsa/rawmidi.h:68
enum__snd_rawmidi_type = c_int
SND_RAWMIDI_TYPE_HW = 1
SND_RAWMIDI_TYPE_SHM = 2
SND_RAWMIDI_TYPE_INET = 3
SND_RAWMIDI_TYPE_VIRTUAL = 4
snd_rawmidi_type_t = enum__snd_rawmidi_type 	# /usr/include/alsa/rawmidi.h:80
# /usr/include/alsa/rawmidi.h:82
snd_rawmidi_open = _lib.snd_rawmidi_open
snd_rawmidi_open.restype = c_int
snd_rawmidi_open.argtypes = [POINTER(POINTER(snd_rawmidi_t)), POINTER(POINTER(snd_rawmidi_t)), c_char_p, c_int]

# /usr/include/alsa/rawmidi.h:84
snd_rawmidi_open_lconf = _lib.snd_rawmidi_open_lconf
snd_rawmidi_open_lconf.restype = c_int
snd_rawmidi_open_lconf.argtypes = [POINTER(POINTER(snd_rawmidi_t)), POINTER(POINTER(snd_rawmidi_t)), c_char_p, c_int, POINTER(snd_config_t)]

# /usr/include/alsa/rawmidi.h:86
snd_rawmidi_close = _lib.snd_rawmidi_close
snd_rawmidi_close.restype = c_int
snd_rawmidi_close.argtypes = [POINTER(snd_rawmidi_t)]

# /usr/include/alsa/rawmidi.h:87
snd_rawmidi_poll_descriptors_count = _lib.snd_rawmidi_poll_descriptors_count
snd_rawmidi_poll_descriptors_count.restype = c_int
snd_rawmidi_poll_descriptors_count.argtypes = [POINTER(snd_rawmidi_t)]

class struct_pollfd(Structure):
    __slots__ = [
    ]
struct_pollfd._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/rawmidi.h:88
snd_rawmidi_poll_descriptors = _lib.snd_rawmidi_poll_descriptors
snd_rawmidi_poll_descriptors.restype = c_int
snd_rawmidi_poll_descriptors.argtypes = [POINTER(snd_rawmidi_t), POINTER(struct_pollfd), c_uint]

class struct_pollfd(Structure):
    __slots__ = [
    ]
struct_pollfd._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/rawmidi.h:89
snd_rawmidi_poll_descriptors_revents = _lib.snd_rawmidi_poll_descriptors_revents
snd_rawmidi_poll_descriptors_revents.restype = c_int
snd_rawmidi_poll_descriptors_revents.argtypes = [POINTER(snd_rawmidi_t), POINTER(struct_pollfd), c_uint, POINTER(c_ushort)]

# /usr/include/alsa/rawmidi.h:90
snd_rawmidi_nonblock = _lib.snd_rawmidi_nonblock
snd_rawmidi_nonblock.restype = c_int
snd_rawmidi_nonblock.argtypes = [POINTER(snd_rawmidi_t), c_int]

# /usr/include/alsa/rawmidi.h:91
snd_rawmidi_info_sizeof = _lib.snd_rawmidi_info_sizeof
snd_rawmidi_info_sizeof.restype = c_size_t
snd_rawmidi_info_sizeof.argtypes = []

# /usr/include/alsa/rawmidi.h:97
snd_rawmidi_info_malloc = _lib.snd_rawmidi_info_malloc
snd_rawmidi_info_malloc.restype = c_int
snd_rawmidi_info_malloc.argtypes = [POINTER(POINTER(snd_rawmidi_info_t))]

# /usr/include/alsa/rawmidi.h:98
snd_rawmidi_info_free = _lib.snd_rawmidi_info_free
snd_rawmidi_info_free.restype = None
snd_rawmidi_info_free.argtypes = [POINTER(snd_rawmidi_info_t)]

# /usr/include/alsa/rawmidi.h:99
snd_rawmidi_info_copy = _lib.snd_rawmidi_info_copy
snd_rawmidi_info_copy.restype = None
snd_rawmidi_info_copy.argtypes = [POINTER(snd_rawmidi_info_t), POINTER(snd_rawmidi_info_t)]

# /usr/include/alsa/rawmidi.h:100
snd_rawmidi_info_get_device = _lib.snd_rawmidi_info_get_device
snd_rawmidi_info_get_device.restype = c_uint
snd_rawmidi_info_get_device.argtypes = [POINTER(snd_rawmidi_info_t)]

# /usr/include/alsa/rawmidi.h:101
snd_rawmidi_info_get_subdevice = _lib.snd_rawmidi_info_get_subdevice
snd_rawmidi_info_get_subdevice.restype = c_uint
snd_rawmidi_info_get_subdevice.argtypes = [POINTER(snd_rawmidi_info_t)]

# /usr/include/alsa/rawmidi.h:102
snd_rawmidi_info_get_stream = _lib.snd_rawmidi_info_get_stream
snd_rawmidi_info_get_stream.restype = snd_rawmidi_stream_t
snd_rawmidi_info_get_stream.argtypes = [POINTER(snd_rawmidi_info_t)]

# /usr/include/alsa/rawmidi.h:103
snd_rawmidi_info_get_card = _lib.snd_rawmidi_info_get_card
snd_rawmidi_info_get_card.restype = c_int
snd_rawmidi_info_get_card.argtypes = [POINTER(snd_rawmidi_info_t)]

# /usr/include/alsa/rawmidi.h:104
snd_rawmidi_info_get_flags = _lib.snd_rawmidi_info_get_flags
snd_rawmidi_info_get_flags.restype = c_uint
snd_rawmidi_info_get_flags.argtypes = [POINTER(snd_rawmidi_info_t)]

# /usr/include/alsa/rawmidi.h:105
snd_rawmidi_info_get_id = _lib.snd_rawmidi_info_get_id
snd_rawmidi_info_get_id.restype = c_char_p
snd_rawmidi_info_get_id.argtypes = [POINTER(snd_rawmidi_info_t)]

# /usr/include/alsa/rawmidi.h:106
snd_rawmidi_info_get_name = _lib.snd_rawmidi_info_get_name
snd_rawmidi_info_get_name.restype = c_char_p
snd_rawmidi_info_get_name.argtypes = [POINTER(snd_rawmidi_info_t)]

# /usr/include/alsa/rawmidi.h:107
snd_rawmidi_info_get_subdevice_name = _lib.snd_rawmidi_info_get_subdevice_name
snd_rawmidi_info_get_subdevice_name.restype = c_char_p
snd_rawmidi_info_get_subdevice_name.argtypes = [POINTER(snd_rawmidi_info_t)]

# /usr/include/alsa/rawmidi.h:108
snd_rawmidi_info_get_subdevices_count = _lib.snd_rawmidi_info_get_subdevices_count
snd_rawmidi_info_get_subdevices_count.restype = c_uint
snd_rawmidi_info_get_subdevices_count.argtypes = [POINTER(snd_rawmidi_info_t)]

# /usr/include/alsa/rawmidi.h:109
snd_rawmidi_info_get_subdevices_avail = _lib.snd_rawmidi_info_get_subdevices_avail
snd_rawmidi_info_get_subdevices_avail.restype = c_uint
snd_rawmidi_info_get_subdevices_avail.argtypes = [POINTER(snd_rawmidi_info_t)]

# /usr/include/alsa/rawmidi.h:110
snd_rawmidi_info_set_device = _lib.snd_rawmidi_info_set_device
snd_rawmidi_info_set_device.restype = None
snd_rawmidi_info_set_device.argtypes = [POINTER(snd_rawmidi_info_t), c_uint]

# /usr/include/alsa/rawmidi.h:111
snd_rawmidi_info_set_subdevice = _lib.snd_rawmidi_info_set_subdevice
snd_rawmidi_info_set_subdevice.restype = None
snd_rawmidi_info_set_subdevice.argtypes = [POINTER(snd_rawmidi_info_t), c_uint]

# /usr/include/alsa/rawmidi.h:112
snd_rawmidi_info_set_stream = _lib.snd_rawmidi_info_set_stream
snd_rawmidi_info_set_stream.restype = None
snd_rawmidi_info_set_stream.argtypes = [POINTER(snd_rawmidi_info_t), snd_rawmidi_stream_t]

# /usr/include/alsa/rawmidi.h:113
snd_rawmidi_info = _lib.snd_rawmidi_info
snd_rawmidi_info.restype = c_int
snd_rawmidi_info.argtypes = [POINTER(snd_rawmidi_t), POINTER(snd_rawmidi_info_t)]

# /usr/include/alsa/rawmidi.h:114
snd_rawmidi_params_sizeof = _lib.snd_rawmidi_params_sizeof
snd_rawmidi_params_sizeof.restype = c_size_t
snd_rawmidi_params_sizeof.argtypes = []

# /usr/include/alsa/rawmidi.h:120
snd_rawmidi_params_malloc = _lib.snd_rawmidi_params_malloc
snd_rawmidi_params_malloc.restype = c_int
snd_rawmidi_params_malloc.argtypes = [POINTER(POINTER(snd_rawmidi_params_t))]

# /usr/include/alsa/rawmidi.h:121
snd_rawmidi_params_free = _lib.snd_rawmidi_params_free
snd_rawmidi_params_free.restype = None
snd_rawmidi_params_free.argtypes = [POINTER(snd_rawmidi_params_t)]

# /usr/include/alsa/rawmidi.h:122
snd_rawmidi_params_copy = _lib.snd_rawmidi_params_copy
snd_rawmidi_params_copy.restype = None
snd_rawmidi_params_copy.argtypes = [POINTER(snd_rawmidi_params_t), POINTER(snd_rawmidi_params_t)]

# /usr/include/alsa/rawmidi.h:123
snd_rawmidi_params_set_buffer_size = _lib.snd_rawmidi_params_set_buffer_size
snd_rawmidi_params_set_buffer_size.restype = c_int
snd_rawmidi_params_set_buffer_size.argtypes = [POINTER(snd_rawmidi_t), POINTER(snd_rawmidi_params_t), c_size_t]

# /usr/include/alsa/rawmidi.h:124
snd_rawmidi_params_get_buffer_size = _lib.snd_rawmidi_params_get_buffer_size
snd_rawmidi_params_get_buffer_size.restype = c_size_t
snd_rawmidi_params_get_buffer_size.argtypes = [POINTER(snd_rawmidi_params_t)]

# /usr/include/alsa/rawmidi.h:125
snd_rawmidi_params_set_avail_min = _lib.snd_rawmidi_params_set_avail_min
snd_rawmidi_params_set_avail_min.restype = c_int
snd_rawmidi_params_set_avail_min.argtypes = [POINTER(snd_rawmidi_t), POINTER(snd_rawmidi_params_t), c_size_t]

# /usr/include/alsa/rawmidi.h:126
snd_rawmidi_params_get_avail_min = _lib.snd_rawmidi_params_get_avail_min
snd_rawmidi_params_get_avail_min.restype = c_size_t
snd_rawmidi_params_get_avail_min.argtypes = [POINTER(snd_rawmidi_params_t)]

# /usr/include/alsa/rawmidi.h:127
snd_rawmidi_params_set_no_active_sensing = _lib.snd_rawmidi_params_set_no_active_sensing
snd_rawmidi_params_set_no_active_sensing.restype = c_int
snd_rawmidi_params_set_no_active_sensing.argtypes = [POINTER(snd_rawmidi_t), POINTER(snd_rawmidi_params_t), c_int]

# /usr/include/alsa/rawmidi.h:128
snd_rawmidi_params_get_no_active_sensing = _lib.snd_rawmidi_params_get_no_active_sensing
snd_rawmidi_params_get_no_active_sensing.restype = c_int
snd_rawmidi_params_get_no_active_sensing.argtypes = [POINTER(snd_rawmidi_params_t)]

# /usr/include/alsa/rawmidi.h:129
snd_rawmidi_params = _lib.snd_rawmidi_params
snd_rawmidi_params.restype = c_int
snd_rawmidi_params.argtypes = [POINTER(snd_rawmidi_t), POINTER(snd_rawmidi_params_t)]

# /usr/include/alsa/rawmidi.h:130
snd_rawmidi_params_current = _lib.snd_rawmidi_params_current
snd_rawmidi_params_current.restype = c_int
snd_rawmidi_params_current.argtypes = [POINTER(snd_rawmidi_t), POINTER(snd_rawmidi_params_t)]

# /usr/include/alsa/rawmidi.h:131
snd_rawmidi_status_sizeof = _lib.snd_rawmidi_status_sizeof
snd_rawmidi_status_sizeof.restype = c_size_t
snd_rawmidi_status_sizeof.argtypes = []

# /usr/include/alsa/rawmidi.h:137
snd_rawmidi_status_malloc = _lib.snd_rawmidi_status_malloc
snd_rawmidi_status_malloc.restype = c_int
snd_rawmidi_status_malloc.argtypes = [POINTER(POINTER(snd_rawmidi_status_t))]

# /usr/include/alsa/rawmidi.h:138
snd_rawmidi_status_free = _lib.snd_rawmidi_status_free
snd_rawmidi_status_free.restype = None
snd_rawmidi_status_free.argtypes = [POINTER(snd_rawmidi_status_t)]

# /usr/include/alsa/rawmidi.h:139
snd_rawmidi_status_copy = _lib.snd_rawmidi_status_copy
snd_rawmidi_status_copy.restype = None
snd_rawmidi_status_copy.argtypes = [POINTER(snd_rawmidi_status_t), POINTER(snd_rawmidi_status_t)]

# /usr/include/alsa/rawmidi.h:140
snd_rawmidi_status_get_tstamp = _lib.snd_rawmidi_status_get_tstamp
snd_rawmidi_status_get_tstamp.restype = None
snd_rawmidi_status_get_tstamp.argtypes = [POINTER(snd_rawmidi_status_t), POINTER(snd_htimestamp_t)]

# /usr/include/alsa/rawmidi.h:141
snd_rawmidi_status_get_avail = _lib.snd_rawmidi_status_get_avail
snd_rawmidi_status_get_avail.restype = c_size_t
snd_rawmidi_status_get_avail.argtypes = [POINTER(snd_rawmidi_status_t)]

# /usr/include/alsa/rawmidi.h:142
snd_rawmidi_status_get_xruns = _lib.snd_rawmidi_status_get_xruns
snd_rawmidi_status_get_xruns.restype = c_size_t
snd_rawmidi_status_get_xruns.argtypes = [POINTER(snd_rawmidi_status_t)]

# /usr/include/alsa/rawmidi.h:143
snd_rawmidi_status = _lib.snd_rawmidi_status
snd_rawmidi_status.restype = c_int
snd_rawmidi_status.argtypes = [POINTER(snd_rawmidi_t), POINTER(snd_rawmidi_status_t)]

# /usr/include/alsa/rawmidi.h:144
snd_rawmidi_drain = _lib.snd_rawmidi_drain
snd_rawmidi_drain.restype = c_int
snd_rawmidi_drain.argtypes = [POINTER(snd_rawmidi_t)]

# /usr/include/alsa/rawmidi.h:145
snd_rawmidi_drop = _lib.snd_rawmidi_drop
snd_rawmidi_drop.restype = c_int
snd_rawmidi_drop.argtypes = [POINTER(snd_rawmidi_t)]

# /usr/include/alsa/rawmidi.h:146
snd_rawmidi_write = _lib.snd_rawmidi_write
snd_rawmidi_write.restype = ssize_t
snd_rawmidi_write.argtypes = [POINTER(snd_rawmidi_t), POINTER(None), c_size_t]

# /usr/include/alsa/rawmidi.h:147
snd_rawmidi_read = _lib.snd_rawmidi_read
snd_rawmidi_read.restype = ssize_t
snd_rawmidi_read.argtypes = [POINTER(snd_rawmidi_t), POINTER(None), c_size_t]

# /usr/include/alsa/rawmidi.h:148
snd_rawmidi_name = _lib.snd_rawmidi_name
snd_rawmidi_name.restype = c_char_p
snd_rawmidi_name.argtypes = [POINTER(snd_rawmidi_t)]

# /usr/include/alsa/rawmidi.h:149
snd_rawmidi_type = _lib.snd_rawmidi_type
snd_rawmidi_type.restype = snd_rawmidi_type_t
snd_rawmidi_type.argtypes = [POINTER(snd_rawmidi_t)]

# /usr/include/alsa/rawmidi.h:150
snd_rawmidi_stream = _lib.snd_rawmidi_stream
snd_rawmidi_stream.restype = snd_rawmidi_stream_t
snd_rawmidi_stream.argtypes = [POINTER(snd_rawmidi_t)]

SND_TIMER_DLSYM_VERSION = 0 	# /usr/include/alsa/timer.h:42
SND_TIMER_QUERY_DLSYM_VERSION = 0 	# /usr/include/alsa/timer.h:44
class struct__snd_timer_id(Structure):
    __slots__ = [
    ]
struct__snd_timer_id._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_timer_id(Structure):
    __slots__ = [
    ]
struct__snd_timer_id._fields_ = [
    ('_opaque_struct', c_int)
]

snd_timer_id_t = struct__snd_timer_id 	# /usr/include/alsa/timer.h:47
class struct__snd_timer_ginfo(Structure):
    __slots__ = [
    ]
struct__snd_timer_ginfo._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_timer_ginfo(Structure):
    __slots__ = [
    ]
struct__snd_timer_ginfo._fields_ = [
    ('_opaque_struct', c_int)
]

snd_timer_ginfo_t = struct__snd_timer_ginfo 	# /usr/include/alsa/timer.h:49
class struct__snd_timer_gparams(Structure):
    __slots__ = [
    ]
struct__snd_timer_gparams._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_timer_gparams(Structure):
    __slots__ = [
    ]
struct__snd_timer_gparams._fields_ = [
    ('_opaque_struct', c_int)
]

snd_timer_gparams_t = struct__snd_timer_gparams 	# /usr/include/alsa/timer.h:51
class struct__snd_timer_gstatus(Structure):
    __slots__ = [
    ]
struct__snd_timer_gstatus._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_timer_gstatus(Structure):
    __slots__ = [
    ]
struct__snd_timer_gstatus._fields_ = [
    ('_opaque_struct', c_int)
]

snd_timer_gstatus_t = struct__snd_timer_gstatus 	# /usr/include/alsa/timer.h:53
class struct__snd_timer_info(Structure):
    __slots__ = [
    ]
struct__snd_timer_info._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_timer_info(Structure):
    __slots__ = [
    ]
struct__snd_timer_info._fields_ = [
    ('_opaque_struct', c_int)
]

snd_timer_info_t = struct__snd_timer_info 	# /usr/include/alsa/timer.h:55
class struct__snd_timer_params(Structure):
    __slots__ = [
    ]
struct__snd_timer_params._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_timer_params(Structure):
    __slots__ = [
    ]
struct__snd_timer_params._fields_ = [
    ('_opaque_struct', c_int)
]

snd_timer_params_t = struct__snd_timer_params 	# /usr/include/alsa/timer.h:57
class struct__snd_timer_status(Structure):
    __slots__ = [
    ]
struct__snd_timer_status._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_timer_status(Structure):
    __slots__ = [
    ]
struct__snd_timer_status._fields_ = [
    ('_opaque_struct', c_int)
]

snd_timer_status_t = struct__snd_timer_status 	# /usr/include/alsa/timer.h:59
enum__snd_timer_class = c_int
SND_TIMER_CLASS_NONE = 1
SND_TIMER_CLASS_SLAVE = 0
SND_TIMER_CLASS_GLOBAL = 1
SND_TIMER_CLASS_CARD = 2
SND_TIMER_CLASS_PCM = 3
SND_TIMER_CLASS_LAST = 0
snd_timer_class_t = enum__snd_timer_class 	# /usr/include/alsa/timer.h:68
enum__snd_timer_slave_class = c_int
SND_TIMER_SCLASS_NONE = 0
SND_TIMER_SCLASS_APPLICATION = 1
SND_TIMER_SCLASS_SEQUENCER = 2
SND_TIMER_SCLASS_OSS_SEQUENCER = 3
SND_TIMER_SCLASS_LAST = 0
snd_timer_slave_class_t = enum__snd_timer_slave_class 	# /usr/include/alsa/timer.h:77
enum__snd_timer_event = c_int
SND_TIMER_EVENT_RESOLUTION = 0
SND_TIMER_EVENT_TICK = 1
SND_TIMER_EVENT_START = 2
SND_TIMER_EVENT_STOP = 3
SND_TIMER_EVENT_CONTINUE = 4
SND_TIMER_EVENT_PAUSE = 5
SND_TIMER_EVENT_EARLY = 6
SND_TIMER_EVENT_SUSPEND = 7
SND_TIMER_EVENT_RESUME = 8
SND_TIMER_EVENT_MSTART = 9
SND_TIMER_EVENT_MSTOP = 10
SND_TIMER_EVENT_MCONTINUE = 11
SND_TIMER_EVENT_MPAUSE = 12
SND_TIMER_EVENT_MSUSPEND = 13
SND_TIMER_EVENT_MRESUME = 14
snd_timer_event_t = enum__snd_timer_event 	# /usr/include/alsa/timer.h:97
class struct__snd_timer_read(Structure):
    __slots__ = [
        'resolution',
        'ticks',
    ]
struct__snd_timer_read._fields_ = [
    ('resolution', c_uint),
    ('ticks', c_uint),
]

snd_timer_read_t = struct__snd_timer_read 	# /usr/include/alsa/timer.h:103
class struct__snd_timer_tread(Structure):
    __slots__ = [
        'event',
        'tstamp',
        'val',
    ]
struct__snd_timer_tread._fields_ = [
    ('event', snd_timer_event_t),
    ('tstamp', snd_htimestamp_t),
    ('val', c_uint),
]

snd_timer_tread_t = struct__snd_timer_tread 	# /usr/include/alsa/timer.h:110
SND_TIMER_GLOBAL_SYSTEM = 0 	# /usr/include/alsa/timer.h:113
SND_TIMER_GLOBAL_RTC = 1 	# /usr/include/alsa/timer.h:115
SND_TIMER_GLOBAL_HPET = 2 	# /usr/include/alsa/timer.h:117
SND_TIMER_OPEN_NONBLOCK = 1 	# /usr/include/alsa/timer.h:120
SND_TIMER_OPEN_TREAD = 2 	# /usr/include/alsa/timer.h:122
enum__snd_timer_type = c_int
SND_TIMER_TYPE_HW = 0
SND_TIMER_TYPE_SHM = 1
SND_TIMER_TYPE_INET = 2
snd_timer_type_t = enum__snd_timer_type 	# /usr/include/alsa/timer.h:132
class struct__snd_timer_query(Structure):
    __slots__ = [
    ]
struct__snd_timer_query._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_timer_query(Structure):
    __slots__ = [
    ]
struct__snd_timer_query._fields_ = [
    ('_opaque_struct', c_int)
]

snd_timer_query_t = struct__snd_timer_query 	# /usr/include/alsa/timer.h:135
class struct__snd_timer(Structure):
    __slots__ = [
    ]
struct__snd_timer._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_timer(Structure):
    __slots__ = [
    ]
struct__snd_timer._fields_ = [
    ('_opaque_struct', c_int)
]

snd_timer_t = struct__snd_timer 	# /usr/include/alsa/timer.h:137
# /usr/include/alsa/timer.h:140
snd_timer_query_open = _lib.snd_timer_query_open
snd_timer_query_open.restype = c_int
snd_timer_query_open.argtypes = [POINTER(POINTER(snd_timer_query_t)), c_char_p, c_int]

# /usr/include/alsa/timer.h:141
snd_timer_query_open_lconf = _lib.snd_timer_query_open_lconf
snd_timer_query_open_lconf.restype = c_int
snd_timer_query_open_lconf.argtypes = [POINTER(POINTER(snd_timer_query_t)), c_char_p, c_int, POINTER(snd_config_t)]

# /usr/include/alsa/timer.h:142
snd_timer_query_close = _lib.snd_timer_query_close
snd_timer_query_close.restype = c_int
snd_timer_query_close.argtypes = [POINTER(snd_timer_query_t)]

# /usr/include/alsa/timer.h:143
snd_timer_query_next_device = _lib.snd_timer_query_next_device
snd_timer_query_next_device.restype = c_int
snd_timer_query_next_device.argtypes = [POINTER(snd_timer_query_t), POINTER(snd_timer_id_t)]

# /usr/include/alsa/timer.h:144
snd_timer_query_info = _lib.snd_timer_query_info
snd_timer_query_info.restype = c_int
snd_timer_query_info.argtypes = [POINTER(snd_timer_query_t), POINTER(snd_timer_ginfo_t)]

# /usr/include/alsa/timer.h:145
snd_timer_query_params = _lib.snd_timer_query_params
snd_timer_query_params.restype = c_int
snd_timer_query_params.argtypes = [POINTER(snd_timer_query_t), POINTER(snd_timer_gparams_t)]

# /usr/include/alsa/timer.h:146
snd_timer_query_status = _lib.snd_timer_query_status
snd_timer_query_status.restype = c_int
snd_timer_query_status.argtypes = [POINTER(snd_timer_query_t), POINTER(snd_timer_gstatus_t)]

# /usr/include/alsa/timer.h:148
snd_timer_open = _lib.snd_timer_open
snd_timer_open.restype = c_int
snd_timer_open.argtypes = [POINTER(POINTER(snd_timer_t)), c_char_p, c_int]

# /usr/include/alsa/timer.h:149
snd_timer_open_lconf = _lib.snd_timer_open_lconf
snd_timer_open_lconf.restype = c_int
snd_timer_open_lconf.argtypes = [POINTER(POINTER(snd_timer_t)), c_char_p, c_int, POINTER(snd_config_t)]

# /usr/include/alsa/timer.h:150
snd_timer_close = _lib.snd_timer_close
snd_timer_close.restype = c_int
snd_timer_close.argtypes = [POINTER(snd_timer_t)]

# /usr/include/alsa/timer.h:151
snd_async_add_timer_handler = _lib.snd_async_add_timer_handler
snd_async_add_timer_handler.restype = c_int
snd_async_add_timer_handler.argtypes = [POINTER(POINTER(snd_async_handler_t)), POINTER(snd_timer_t), snd_async_callback_t, POINTER(None)]

# /usr/include/alsa/timer.h:153
snd_async_handler_get_timer = _lib.snd_async_handler_get_timer
snd_async_handler_get_timer.restype = POINTER(snd_timer_t)
snd_async_handler_get_timer.argtypes = [POINTER(snd_async_handler_t)]

# /usr/include/alsa/timer.h:154
snd_timer_poll_descriptors_count = _lib.snd_timer_poll_descriptors_count
snd_timer_poll_descriptors_count.restype = c_int
snd_timer_poll_descriptors_count.argtypes = [POINTER(snd_timer_t)]

class struct_pollfd(Structure):
    __slots__ = [
    ]
struct_pollfd._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/timer.h:155
snd_timer_poll_descriptors = _lib.snd_timer_poll_descriptors
snd_timer_poll_descriptors.restype = c_int
snd_timer_poll_descriptors.argtypes = [POINTER(snd_timer_t), POINTER(struct_pollfd), c_uint]

class struct_pollfd(Structure):
    __slots__ = [
    ]
struct_pollfd._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/timer.h:156
snd_timer_poll_descriptors_revents = _lib.snd_timer_poll_descriptors_revents
snd_timer_poll_descriptors_revents.restype = c_int
snd_timer_poll_descriptors_revents.argtypes = [POINTER(snd_timer_t), POINTER(struct_pollfd), c_uint, POINTER(c_ushort)]

# /usr/include/alsa/timer.h:157
snd_timer_info = _lib.snd_timer_info
snd_timer_info.restype = c_int
snd_timer_info.argtypes = [POINTER(snd_timer_t), POINTER(snd_timer_info_t)]

# /usr/include/alsa/timer.h:158
snd_timer_params = _lib.snd_timer_params
snd_timer_params.restype = c_int
snd_timer_params.argtypes = [POINTER(snd_timer_t), POINTER(snd_timer_params_t)]

# /usr/include/alsa/timer.h:159
snd_timer_status = _lib.snd_timer_status
snd_timer_status.restype = c_int
snd_timer_status.argtypes = [POINTER(snd_timer_t), POINTER(snd_timer_status_t)]

# /usr/include/alsa/timer.h:160
snd_timer_start = _lib.snd_timer_start
snd_timer_start.restype = c_int
snd_timer_start.argtypes = [POINTER(snd_timer_t)]

# /usr/include/alsa/timer.h:161
snd_timer_stop = _lib.snd_timer_stop
snd_timer_stop.restype = c_int
snd_timer_stop.argtypes = [POINTER(snd_timer_t)]

# /usr/include/alsa/timer.h:162
snd_timer_continue = _lib.snd_timer_continue
snd_timer_continue.restype = c_int
snd_timer_continue.argtypes = [POINTER(snd_timer_t)]

# /usr/include/alsa/timer.h:163
snd_timer_read = _lib.snd_timer_read
snd_timer_read.restype = ssize_t
snd_timer_read.argtypes = [POINTER(snd_timer_t), POINTER(None), c_size_t]

# /usr/include/alsa/timer.h:165
snd_timer_id_sizeof = _lib.snd_timer_id_sizeof
snd_timer_id_sizeof.restype = c_size_t
snd_timer_id_sizeof.argtypes = []

# /usr/include/alsa/timer.h:168
snd_timer_id_malloc = _lib.snd_timer_id_malloc
snd_timer_id_malloc.restype = c_int
snd_timer_id_malloc.argtypes = [POINTER(POINTER(snd_timer_id_t))]

# /usr/include/alsa/timer.h:169
snd_timer_id_free = _lib.snd_timer_id_free
snd_timer_id_free.restype = None
snd_timer_id_free.argtypes = [POINTER(snd_timer_id_t)]

# /usr/include/alsa/timer.h:170
snd_timer_id_copy = _lib.snd_timer_id_copy
snd_timer_id_copy.restype = None
snd_timer_id_copy.argtypes = [POINTER(snd_timer_id_t), POINTER(snd_timer_id_t)]

# /usr/include/alsa/timer.h:172
snd_timer_id_set_class = _lib.snd_timer_id_set_class
snd_timer_id_set_class.restype = None
snd_timer_id_set_class.argtypes = [POINTER(snd_timer_id_t), c_int]

# /usr/include/alsa/timer.h:173
snd_timer_id_get_class = _lib.snd_timer_id_get_class
snd_timer_id_get_class.restype = c_int
snd_timer_id_get_class.argtypes = [POINTER(snd_timer_id_t)]

# /usr/include/alsa/timer.h:174
snd_timer_id_set_sclass = _lib.snd_timer_id_set_sclass
snd_timer_id_set_sclass.restype = None
snd_timer_id_set_sclass.argtypes = [POINTER(snd_timer_id_t), c_int]

# /usr/include/alsa/timer.h:175
snd_timer_id_get_sclass = _lib.snd_timer_id_get_sclass
snd_timer_id_get_sclass.restype = c_int
snd_timer_id_get_sclass.argtypes = [POINTER(snd_timer_id_t)]

# /usr/include/alsa/timer.h:176
snd_timer_id_set_card = _lib.snd_timer_id_set_card
snd_timer_id_set_card.restype = None
snd_timer_id_set_card.argtypes = [POINTER(snd_timer_id_t), c_int]

# /usr/include/alsa/timer.h:177
snd_timer_id_get_card = _lib.snd_timer_id_get_card
snd_timer_id_get_card.restype = c_int
snd_timer_id_get_card.argtypes = [POINTER(snd_timer_id_t)]

# /usr/include/alsa/timer.h:178
snd_timer_id_set_device = _lib.snd_timer_id_set_device
snd_timer_id_set_device.restype = None
snd_timer_id_set_device.argtypes = [POINTER(snd_timer_id_t), c_int]

# /usr/include/alsa/timer.h:179
snd_timer_id_get_device = _lib.snd_timer_id_get_device
snd_timer_id_get_device.restype = c_int
snd_timer_id_get_device.argtypes = [POINTER(snd_timer_id_t)]

# /usr/include/alsa/timer.h:180
snd_timer_id_set_subdevice = _lib.snd_timer_id_set_subdevice
snd_timer_id_set_subdevice.restype = None
snd_timer_id_set_subdevice.argtypes = [POINTER(snd_timer_id_t), c_int]

# /usr/include/alsa/timer.h:181
snd_timer_id_get_subdevice = _lib.snd_timer_id_get_subdevice
snd_timer_id_get_subdevice.restype = c_int
snd_timer_id_get_subdevice.argtypes = [POINTER(snd_timer_id_t)]

# /usr/include/alsa/timer.h:183
snd_timer_ginfo_sizeof = _lib.snd_timer_ginfo_sizeof
snd_timer_ginfo_sizeof.restype = c_size_t
snd_timer_ginfo_sizeof.argtypes = []

# /usr/include/alsa/timer.h:186
snd_timer_ginfo_malloc = _lib.snd_timer_ginfo_malloc
snd_timer_ginfo_malloc.restype = c_int
snd_timer_ginfo_malloc.argtypes = [POINTER(POINTER(snd_timer_ginfo_t))]

# /usr/include/alsa/timer.h:187
snd_timer_ginfo_free = _lib.snd_timer_ginfo_free
snd_timer_ginfo_free.restype = None
snd_timer_ginfo_free.argtypes = [POINTER(snd_timer_ginfo_t)]

# /usr/include/alsa/timer.h:188
snd_timer_ginfo_copy = _lib.snd_timer_ginfo_copy
snd_timer_ginfo_copy.restype = None
snd_timer_ginfo_copy.argtypes = [POINTER(snd_timer_ginfo_t), POINTER(snd_timer_ginfo_t)]

# /usr/include/alsa/timer.h:190
snd_timer_ginfo_set_tid = _lib.snd_timer_ginfo_set_tid
snd_timer_ginfo_set_tid.restype = c_int
snd_timer_ginfo_set_tid.argtypes = [POINTER(snd_timer_ginfo_t), POINTER(snd_timer_id_t)]

# /usr/include/alsa/timer.h:191
snd_timer_ginfo_get_tid = _lib.snd_timer_ginfo_get_tid
snd_timer_ginfo_get_tid.restype = POINTER(snd_timer_id_t)
snd_timer_ginfo_get_tid.argtypes = [POINTER(snd_timer_ginfo_t)]

# /usr/include/alsa/timer.h:192
snd_timer_ginfo_get_flags = _lib.snd_timer_ginfo_get_flags
snd_timer_ginfo_get_flags.restype = c_uint
snd_timer_ginfo_get_flags.argtypes = [POINTER(snd_timer_ginfo_t)]

# /usr/include/alsa/timer.h:193
snd_timer_ginfo_get_card = _lib.snd_timer_ginfo_get_card
snd_timer_ginfo_get_card.restype = c_int
snd_timer_ginfo_get_card.argtypes = [POINTER(snd_timer_ginfo_t)]

# /usr/include/alsa/timer.h:194
snd_timer_ginfo_get_id = _lib.snd_timer_ginfo_get_id
snd_timer_ginfo_get_id.restype = c_char_p
snd_timer_ginfo_get_id.argtypes = [POINTER(snd_timer_ginfo_t)]

# /usr/include/alsa/timer.h:195
snd_timer_ginfo_get_name = _lib.snd_timer_ginfo_get_name
snd_timer_ginfo_get_name.restype = c_char_p
snd_timer_ginfo_get_name.argtypes = [POINTER(snd_timer_ginfo_t)]

# /usr/include/alsa/timer.h:196
snd_timer_ginfo_get_resolution = _lib.snd_timer_ginfo_get_resolution
snd_timer_ginfo_get_resolution.restype = c_ulong
snd_timer_ginfo_get_resolution.argtypes = [POINTER(snd_timer_ginfo_t)]

# /usr/include/alsa/timer.h:197
snd_timer_ginfo_get_resolution_min = _lib.snd_timer_ginfo_get_resolution_min
snd_timer_ginfo_get_resolution_min.restype = c_ulong
snd_timer_ginfo_get_resolution_min.argtypes = [POINTER(snd_timer_ginfo_t)]

# /usr/include/alsa/timer.h:198
snd_timer_ginfo_get_resolution_max = _lib.snd_timer_ginfo_get_resolution_max
snd_timer_ginfo_get_resolution_max.restype = c_ulong
snd_timer_ginfo_get_resolution_max.argtypes = [POINTER(snd_timer_ginfo_t)]

# /usr/include/alsa/timer.h:199
snd_timer_ginfo_get_clients = _lib.snd_timer_ginfo_get_clients
snd_timer_ginfo_get_clients.restype = c_uint
snd_timer_ginfo_get_clients.argtypes = [POINTER(snd_timer_ginfo_t)]

# /usr/include/alsa/timer.h:201
snd_timer_info_sizeof = _lib.snd_timer_info_sizeof
snd_timer_info_sizeof.restype = c_size_t
snd_timer_info_sizeof.argtypes = []

# /usr/include/alsa/timer.h:204
snd_timer_info_malloc = _lib.snd_timer_info_malloc
snd_timer_info_malloc.restype = c_int
snd_timer_info_malloc.argtypes = [POINTER(POINTER(snd_timer_info_t))]

# /usr/include/alsa/timer.h:205
snd_timer_info_free = _lib.snd_timer_info_free
snd_timer_info_free.restype = None
snd_timer_info_free.argtypes = [POINTER(snd_timer_info_t)]

# /usr/include/alsa/timer.h:206
snd_timer_info_copy = _lib.snd_timer_info_copy
snd_timer_info_copy.restype = None
snd_timer_info_copy.argtypes = [POINTER(snd_timer_info_t), POINTER(snd_timer_info_t)]

# /usr/include/alsa/timer.h:208
snd_timer_info_is_slave = _lib.snd_timer_info_is_slave
snd_timer_info_is_slave.restype = c_int
snd_timer_info_is_slave.argtypes = [POINTER(snd_timer_info_t)]

# /usr/include/alsa/timer.h:209
snd_timer_info_get_card = _lib.snd_timer_info_get_card
snd_timer_info_get_card.restype = c_int
snd_timer_info_get_card.argtypes = [POINTER(snd_timer_info_t)]

# /usr/include/alsa/timer.h:210
snd_timer_info_get_id = _lib.snd_timer_info_get_id
snd_timer_info_get_id.restype = c_char_p
snd_timer_info_get_id.argtypes = [POINTER(snd_timer_info_t)]

# /usr/include/alsa/timer.h:211
snd_timer_info_get_name = _lib.snd_timer_info_get_name
snd_timer_info_get_name.restype = c_char_p
snd_timer_info_get_name.argtypes = [POINTER(snd_timer_info_t)]

# /usr/include/alsa/timer.h:212
snd_timer_info_get_resolution = _lib.snd_timer_info_get_resolution
snd_timer_info_get_resolution.restype = c_long
snd_timer_info_get_resolution.argtypes = [POINTER(snd_timer_info_t)]

# /usr/include/alsa/timer.h:214
snd_timer_params_sizeof = _lib.snd_timer_params_sizeof
snd_timer_params_sizeof.restype = c_size_t
snd_timer_params_sizeof.argtypes = []

# /usr/include/alsa/timer.h:217
snd_timer_params_malloc = _lib.snd_timer_params_malloc
snd_timer_params_malloc.restype = c_int
snd_timer_params_malloc.argtypes = [POINTER(POINTER(snd_timer_params_t))]

# /usr/include/alsa/timer.h:218
snd_timer_params_free = _lib.snd_timer_params_free
snd_timer_params_free.restype = None
snd_timer_params_free.argtypes = [POINTER(snd_timer_params_t)]

# /usr/include/alsa/timer.h:219
snd_timer_params_copy = _lib.snd_timer_params_copy
snd_timer_params_copy.restype = None
snd_timer_params_copy.argtypes = [POINTER(snd_timer_params_t), POINTER(snd_timer_params_t)]

# /usr/include/alsa/timer.h:221
snd_timer_params_set_auto_start = _lib.snd_timer_params_set_auto_start
snd_timer_params_set_auto_start.restype = c_int
snd_timer_params_set_auto_start.argtypes = [POINTER(snd_timer_params_t), c_int]

# /usr/include/alsa/timer.h:222
snd_timer_params_get_auto_start = _lib.snd_timer_params_get_auto_start
snd_timer_params_get_auto_start.restype = c_int
snd_timer_params_get_auto_start.argtypes = [POINTER(snd_timer_params_t)]

# /usr/include/alsa/timer.h:223
snd_timer_params_set_exclusive = _lib.snd_timer_params_set_exclusive
snd_timer_params_set_exclusive.restype = c_int
snd_timer_params_set_exclusive.argtypes = [POINTER(snd_timer_params_t), c_int]

# /usr/include/alsa/timer.h:224
snd_timer_params_get_exclusive = _lib.snd_timer_params_get_exclusive
snd_timer_params_get_exclusive.restype = c_int
snd_timer_params_get_exclusive.argtypes = [POINTER(snd_timer_params_t)]

# /usr/include/alsa/timer.h:225
snd_timer_params_set_early_event = _lib.snd_timer_params_set_early_event
snd_timer_params_set_early_event.restype = c_int
snd_timer_params_set_early_event.argtypes = [POINTER(snd_timer_params_t), c_int]

# /usr/include/alsa/timer.h:226
snd_timer_params_get_early_event = _lib.snd_timer_params_get_early_event
snd_timer_params_get_early_event.restype = c_int
snd_timer_params_get_early_event.argtypes = [POINTER(snd_timer_params_t)]

# /usr/include/alsa/timer.h:227
snd_timer_params_set_ticks = _lib.snd_timer_params_set_ticks
snd_timer_params_set_ticks.restype = None
snd_timer_params_set_ticks.argtypes = [POINTER(snd_timer_params_t), c_long]

# /usr/include/alsa/timer.h:228
snd_timer_params_get_ticks = _lib.snd_timer_params_get_ticks
snd_timer_params_get_ticks.restype = c_long
snd_timer_params_get_ticks.argtypes = [POINTER(snd_timer_params_t)]

# /usr/include/alsa/timer.h:229
snd_timer_params_set_queue_size = _lib.snd_timer_params_set_queue_size
snd_timer_params_set_queue_size.restype = None
snd_timer_params_set_queue_size.argtypes = [POINTER(snd_timer_params_t), c_long]

# /usr/include/alsa/timer.h:230
snd_timer_params_get_queue_size = _lib.snd_timer_params_get_queue_size
snd_timer_params_get_queue_size.restype = c_long
snd_timer_params_get_queue_size.argtypes = [POINTER(snd_timer_params_t)]

# /usr/include/alsa/timer.h:231
snd_timer_params_set_filter = _lib.snd_timer_params_set_filter
snd_timer_params_set_filter.restype = None
snd_timer_params_set_filter.argtypes = [POINTER(snd_timer_params_t), c_uint]

# /usr/include/alsa/timer.h:232
snd_timer_params_get_filter = _lib.snd_timer_params_get_filter
snd_timer_params_get_filter.restype = c_uint
snd_timer_params_get_filter.argtypes = [POINTER(snd_timer_params_t)]

# /usr/include/alsa/timer.h:234
snd_timer_status_sizeof = _lib.snd_timer_status_sizeof
snd_timer_status_sizeof.restype = c_size_t
snd_timer_status_sizeof.argtypes = []

# /usr/include/alsa/timer.h:237
snd_timer_status_malloc = _lib.snd_timer_status_malloc
snd_timer_status_malloc.restype = c_int
snd_timer_status_malloc.argtypes = [POINTER(POINTER(snd_timer_status_t))]

# /usr/include/alsa/timer.h:238
snd_timer_status_free = _lib.snd_timer_status_free
snd_timer_status_free.restype = None
snd_timer_status_free.argtypes = [POINTER(snd_timer_status_t)]

# /usr/include/alsa/timer.h:239
snd_timer_status_copy = _lib.snd_timer_status_copy
snd_timer_status_copy.restype = None
snd_timer_status_copy.argtypes = [POINTER(snd_timer_status_t), POINTER(snd_timer_status_t)]

# /usr/include/alsa/timer.h:241
snd_timer_status_get_timestamp = _lib.snd_timer_status_get_timestamp
snd_timer_status_get_timestamp.restype = snd_htimestamp_t
snd_timer_status_get_timestamp.argtypes = [POINTER(snd_timer_status_t)]

# /usr/include/alsa/timer.h:242
snd_timer_status_get_resolution = _lib.snd_timer_status_get_resolution
snd_timer_status_get_resolution.restype = c_long
snd_timer_status_get_resolution.argtypes = [POINTER(snd_timer_status_t)]

# /usr/include/alsa/timer.h:243
snd_timer_status_get_lost = _lib.snd_timer_status_get_lost
snd_timer_status_get_lost.restype = c_long
snd_timer_status_get_lost.argtypes = [POINTER(snd_timer_status_t)]

# /usr/include/alsa/timer.h:244
snd_timer_status_get_overrun = _lib.snd_timer_status_get_overrun
snd_timer_status_get_overrun.restype = c_long
snd_timer_status_get_overrun.argtypes = [POINTER(snd_timer_status_t)]

# /usr/include/alsa/timer.h:245
snd_timer_status_get_queue = _lib.snd_timer_status_get_queue
snd_timer_status_get_queue.restype = c_long
snd_timer_status_get_queue.argtypes = [POINTER(snd_timer_status_t)]

# /usr/include/alsa/timer.h:248
snd_timer_info_get_ticks = _lib.snd_timer_info_get_ticks
snd_timer_info_get_ticks.restype = c_long
snd_timer_info_get_ticks.argtypes = [POINTER(snd_timer_info_t)]

SND_HWDEP_DLSYM_VERSION = 0 	# /usr/include/alsa/hwdep.h:42
class struct__snd_hwdep_info(Structure):
    __slots__ = [
    ]
struct__snd_hwdep_info._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_hwdep_info(Structure):
    __slots__ = [
    ]
struct__snd_hwdep_info._fields_ = [
    ('_opaque_struct', c_int)
]

snd_hwdep_info_t = struct__snd_hwdep_info 	# /usr/include/alsa/hwdep.h:45
class struct__snd_hwdep_dsp_status(Structure):
    __slots__ = [
    ]
struct__snd_hwdep_dsp_status._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_hwdep_dsp_status(Structure):
    __slots__ = [
    ]
struct__snd_hwdep_dsp_status._fields_ = [
    ('_opaque_struct', c_int)
]

snd_hwdep_dsp_status_t = struct__snd_hwdep_dsp_status 	# /usr/include/alsa/hwdep.h:48
class struct__snd_hwdep_dsp_image(Structure):
    __slots__ = [
    ]
struct__snd_hwdep_dsp_image._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_hwdep_dsp_image(Structure):
    __slots__ = [
    ]
struct__snd_hwdep_dsp_image._fields_ = [
    ('_opaque_struct', c_int)
]

snd_hwdep_dsp_image_t = struct__snd_hwdep_dsp_image 	# /usr/include/alsa/hwdep.h:51
enum__snd_hwdep_iface = c_int
SND_HWDEP_IFACE_OPL2 = 0
SND_HWDEP_IFACE_OPL3 = 1
SND_HWDEP_IFACE_OPL4 = 2
SND_HWDEP_IFACE_SB16CSP = 3
SND_HWDEP_IFACE_EMU10K1 = 4
SND_HWDEP_IFACE_YSS225 = 5
SND_HWDEP_IFACE_ICS2115 = 6
SND_HWDEP_IFACE_SSCAPE = 7
SND_HWDEP_IFACE_VX = 8
SND_HWDEP_IFACE_MIXART = 9
SND_HWDEP_IFACE_USX2Y = 10
SND_HWDEP_IFACE_EMUX_WAVETABLE = 11
SND_HWDEP_IFACE_BLUETOOTH = 12
SND_HWDEP_IFACE_USX2Y_PCM = 13
SND_HWDEP_IFACE_PCXHR = 14
SND_HWDEP_IFACE_SB_RC = 15
SND_HWDEP_IFACE_LAST = 0
snd_hwdep_iface_t = enum__snd_hwdep_iface 	# /usr/include/alsa/hwdep.h:73
SND_HWDEP_OPEN_READ = 0 	# /usr/include/alsa/hwdep.h:76
SND_HWDEP_OPEN_WRITE = 1 	# /usr/include/alsa/hwdep.h:78
SND_HWDEP_OPEN_DUPLEX = 2 	# /usr/include/alsa/hwdep.h:80
SND_HWDEP_OPEN_NONBLOCK = 2048 	# /usr/include/alsa/hwdep.h:82
enum__snd_hwdep_type = c_int
SND_HWDEP_TYPE_HW = 1
SND_HWDEP_TYPE_SHM = 2
SND_HWDEP_TYPE_INET = 3
snd_hwdep_type_t = enum__snd_hwdep_type 	# /usr/include/alsa/hwdep.h:92
class struct__snd_hwdep(Structure):
    __slots__ = [
    ]
struct__snd_hwdep._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_hwdep(Structure):
    __slots__ = [
    ]
struct__snd_hwdep._fields_ = [
    ('_opaque_struct', c_int)
]

snd_hwdep_t = struct__snd_hwdep 	# /usr/include/alsa/hwdep.h:95
# /usr/include/alsa/hwdep.h:97
snd_hwdep_open = _lib.snd_hwdep_open
snd_hwdep_open.restype = c_int
snd_hwdep_open.argtypes = [POINTER(POINTER(snd_hwdep_t)), c_char_p, c_int]

# /usr/include/alsa/hwdep.h:98
snd_hwdep_close = _lib.snd_hwdep_close
snd_hwdep_close.restype = c_int
snd_hwdep_close.argtypes = [POINTER(snd_hwdep_t)]

class struct_pollfd(Structure):
    __slots__ = [
    ]
struct_pollfd._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/hwdep.h:99
snd_hwdep_poll_descriptors = _lib.snd_hwdep_poll_descriptors
snd_hwdep_poll_descriptors.restype = c_int
snd_hwdep_poll_descriptors.argtypes = [POINTER(snd_hwdep_t), POINTER(struct_pollfd), c_uint]

class struct_pollfd(Structure):
    __slots__ = [
    ]
struct_pollfd._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/hwdep.h:100
snd_hwdep_poll_descriptors_revents = _lib.snd_hwdep_poll_descriptors_revents
snd_hwdep_poll_descriptors_revents.restype = c_int
snd_hwdep_poll_descriptors_revents.argtypes = [POINTER(snd_hwdep_t), POINTER(struct_pollfd), c_uint, POINTER(c_ushort)]

# /usr/include/alsa/hwdep.h:101
snd_hwdep_nonblock = _lib.snd_hwdep_nonblock
snd_hwdep_nonblock.restype = c_int
snd_hwdep_nonblock.argtypes = [POINTER(snd_hwdep_t), c_int]

# /usr/include/alsa/hwdep.h:102
snd_hwdep_info = _lib.snd_hwdep_info
snd_hwdep_info.restype = c_int
snd_hwdep_info.argtypes = [POINTER(snd_hwdep_t), POINTER(snd_hwdep_info_t)]

# /usr/include/alsa/hwdep.h:103
snd_hwdep_dsp_status = _lib.snd_hwdep_dsp_status
snd_hwdep_dsp_status.restype = c_int
snd_hwdep_dsp_status.argtypes = [POINTER(snd_hwdep_t), POINTER(snd_hwdep_dsp_status_t)]

# /usr/include/alsa/hwdep.h:104
snd_hwdep_dsp_load = _lib.snd_hwdep_dsp_load
snd_hwdep_dsp_load.restype = c_int
snd_hwdep_dsp_load.argtypes = [POINTER(snd_hwdep_t), POINTER(snd_hwdep_dsp_image_t)]

# /usr/include/alsa/hwdep.h:105
snd_hwdep_ioctl = _lib.snd_hwdep_ioctl
snd_hwdep_ioctl.restype = c_int
snd_hwdep_ioctl.argtypes = [POINTER(snd_hwdep_t), c_uint, POINTER(None)]

# /usr/include/alsa/hwdep.h:106
snd_hwdep_write = _lib.snd_hwdep_write
snd_hwdep_write.restype = ssize_t
snd_hwdep_write.argtypes = [POINTER(snd_hwdep_t), POINTER(None), c_size_t]

# /usr/include/alsa/hwdep.h:107
snd_hwdep_read = _lib.snd_hwdep_read
snd_hwdep_read.restype = ssize_t
snd_hwdep_read.argtypes = [POINTER(snd_hwdep_t), POINTER(None), c_size_t]

# /usr/include/alsa/hwdep.h:109
snd_hwdep_info_sizeof = _lib.snd_hwdep_info_sizeof
snd_hwdep_info_sizeof.restype = c_size_t
snd_hwdep_info_sizeof.argtypes = []

# /usr/include/alsa/hwdep.h:112
snd_hwdep_info_malloc = _lib.snd_hwdep_info_malloc
snd_hwdep_info_malloc.restype = c_int
snd_hwdep_info_malloc.argtypes = [POINTER(POINTER(snd_hwdep_info_t))]

# /usr/include/alsa/hwdep.h:113
snd_hwdep_info_free = _lib.snd_hwdep_info_free
snd_hwdep_info_free.restype = None
snd_hwdep_info_free.argtypes = [POINTER(snd_hwdep_info_t)]

# /usr/include/alsa/hwdep.h:114
snd_hwdep_info_copy = _lib.snd_hwdep_info_copy
snd_hwdep_info_copy.restype = None
snd_hwdep_info_copy.argtypes = [POINTER(snd_hwdep_info_t), POINTER(snd_hwdep_info_t)]

# /usr/include/alsa/hwdep.h:116
snd_hwdep_info_get_device = _lib.snd_hwdep_info_get_device
snd_hwdep_info_get_device.restype = c_uint
snd_hwdep_info_get_device.argtypes = [POINTER(snd_hwdep_info_t)]

# /usr/include/alsa/hwdep.h:117
snd_hwdep_info_get_card = _lib.snd_hwdep_info_get_card
snd_hwdep_info_get_card.restype = c_int
snd_hwdep_info_get_card.argtypes = [POINTER(snd_hwdep_info_t)]

# /usr/include/alsa/hwdep.h:118
snd_hwdep_info_get_id = _lib.snd_hwdep_info_get_id
snd_hwdep_info_get_id.restype = c_char_p
snd_hwdep_info_get_id.argtypes = [POINTER(snd_hwdep_info_t)]

# /usr/include/alsa/hwdep.h:119
snd_hwdep_info_get_name = _lib.snd_hwdep_info_get_name
snd_hwdep_info_get_name.restype = c_char_p
snd_hwdep_info_get_name.argtypes = [POINTER(snd_hwdep_info_t)]

# /usr/include/alsa/hwdep.h:120
snd_hwdep_info_get_iface = _lib.snd_hwdep_info_get_iface
snd_hwdep_info_get_iface.restype = snd_hwdep_iface_t
snd_hwdep_info_get_iface.argtypes = [POINTER(snd_hwdep_info_t)]

# /usr/include/alsa/hwdep.h:121
snd_hwdep_info_set_device = _lib.snd_hwdep_info_set_device
snd_hwdep_info_set_device.restype = None
snd_hwdep_info_set_device.argtypes = [POINTER(snd_hwdep_info_t), c_uint]

# /usr/include/alsa/hwdep.h:123
snd_hwdep_dsp_status_sizeof = _lib.snd_hwdep_dsp_status_sizeof
snd_hwdep_dsp_status_sizeof.restype = c_size_t
snd_hwdep_dsp_status_sizeof.argtypes = []

# /usr/include/alsa/hwdep.h:126
snd_hwdep_dsp_status_malloc = _lib.snd_hwdep_dsp_status_malloc
snd_hwdep_dsp_status_malloc.restype = c_int
snd_hwdep_dsp_status_malloc.argtypes = [POINTER(POINTER(snd_hwdep_dsp_status_t))]

# /usr/include/alsa/hwdep.h:127
snd_hwdep_dsp_status_free = _lib.snd_hwdep_dsp_status_free
snd_hwdep_dsp_status_free.restype = None
snd_hwdep_dsp_status_free.argtypes = [POINTER(snd_hwdep_dsp_status_t)]

# /usr/include/alsa/hwdep.h:128
snd_hwdep_dsp_status_copy = _lib.snd_hwdep_dsp_status_copy
snd_hwdep_dsp_status_copy.restype = None
snd_hwdep_dsp_status_copy.argtypes = [POINTER(snd_hwdep_dsp_status_t), POINTER(snd_hwdep_dsp_status_t)]

# /usr/include/alsa/hwdep.h:130
snd_hwdep_dsp_status_get_version = _lib.snd_hwdep_dsp_status_get_version
snd_hwdep_dsp_status_get_version.restype = c_uint
snd_hwdep_dsp_status_get_version.argtypes = [POINTER(snd_hwdep_dsp_status_t)]

# /usr/include/alsa/hwdep.h:131
snd_hwdep_dsp_status_get_id = _lib.snd_hwdep_dsp_status_get_id
snd_hwdep_dsp_status_get_id.restype = c_char_p
snd_hwdep_dsp_status_get_id.argtypes = [POINTER(snd_hwdep_dsp_status_t)]

# /usr/include/alsa/hwdep.h:132
snd_hwdep_dsp_status_get_num_dsps = _lib.snd_hwdep_dsp_status_get_num_dsps
snd_hwdep_dsp_status_get_num_dsps.restype = c_uint
snd_hwdep_dsp_status_get_num_dsps.argtypes = [POINTER(snd_hwdep_dsp_status_t)]

# /usr/include/alsa/hwdep.h:133
snd_hwdep_dsp_status_get_dsp_loaded = _lib.snd_hwdep_dsp_status_get_dsp_loaded
snd_hwdep_dsp_status_get_dsp_loaded.restype = c_uint
snd_hwdep_dsp_status_get_dsp_loaded.argtypes = [POINTER(snd_hwdep_dsp_status_t)]

# /usr/include/alsa/hwdep.h:134
snd_hwdep_dsp_status_get_chip_ready = _lib.snd_hwdep_dsp_status_get_chip_ready
snd_hwdep_dsp_status_get_chip_ready.restype = c_uint
snd_hwdep_dsp_status_get_chip_ready.argtypes = [POINTER(snd_hwdep_dsp_status_t)]

# /usr/include/alsa/hwdep.h:136
snd_hwdep_dsp_image_sizeof = _lib.snd_hwdep_dsp_image_sizeof
snd_hwdep_dsp_image_sizeof.restype = c_size_t
snd_hwdep_dsp_image_sizeof.argtypes = []

# /usr/include/alsa/hwdep.h:139
snd_hwdep_dsp_image_malloc = _lib.snd_hwdep_dsp_image_malloc
snd_hwdep_dsp_image_malloc.restype = c_int
snd_hwdep_dsp_image_malloc.argtypes = [POINTER(POINTER(snd_hwdep_dsp_image_t))]

# /usr/include/alsa/hwdep.h:140
snd_hwdep_dsp_image_free = _lib.snd_hwdep_dsp_image_free
snd_hwdep_dsp_image_free.restype = None
snd_hwdep_dsp_image_free.argtypes = [POINTER(snd_hwdep_dsp_image_t)]

# /usr/include/alsa/hwdep.h:141
snd_hwdep_dsp_image_copy = _lib.snd_hwdep_dsp_image_copy
snd_hwdep_dsp_image_copy.restype = None
snd_hwdep_dsp_image_copy.argtypes = [POINTER(snd_hwdep_dsp_image_t), POINTER(snd_hwdep_dsp_image_t)]

# /usr/include/alsa/hwdep.h:143
snd_hwdep_dsp_image_get_index = _lib.snd_hwdep_dsp_image_get_index
snd_hwdep_dsp_image_get_index.restype = c_uint
snd_hwdep_dsp_image_get_index.argtypes = [POINTER(snd_hwdep_dsp_image_t)]

# /usr/include/alsa/hwdep.h:144
snd_hwdep_dsp_image_get_name = _lib.snd_hwdep_dsp_image_get_name
snd_hwdep_dsp_image_get_name.restype = c_char_p
snd_hwdep_dsp_image_get_name.argtypes = [POINTER(snd_hwdep_dsp_image_t)]

# /usr/include/alsa/hwdep.h:145
snd_hwdep_dsp_image_get_image = _lib.snd_hwdep_dsp_image_get_image
snd_hwdep_dsp_image_get_image.restype = POINTER(c_void)
snd_hwdep_dsp_image_get_image.argtypes = [POINTER(snd_hwdep_dsp_image_t)]

# /usr/include/alsa/hwdep.h:146
snd_hwdep_dsp_image_get_length = _lib.snd_hwdep_dsp_image_get_length
snd_hwdep_dsp_image_get_length.restype = c_size_t
snd_hwdep_dsp_image_get_length.argtypes = [POINTER(snd_hwdep_dsp_image_t)]

# /usr/include/alsa/hwdep.h:148
snd_hwdep_dsp_image_set_index = _lib.snd_hwdep_dsp_image_set_index
snd_hwdep_dsp_image_set_index.restype = None
snd_hwdep_dsp_image_set_index.argtypes = [POINTER(snd_hwdep_dsp_image_t), c_uint]

# /usr/include/alsa/hwdep.h:149
snd_hwdep_dsp_image_set_name = _lib.snd_hwdep_dsp_image_set_name
snd_hwdep_dsp_image_set_name.restype = None
snd_hwdep_dsp_image_set_name.argtypes = [POINTER(snd_hwdep_dsp_image_t), c_char_p]

# /usr/include/alsa/hwdep.h:150
snd_hwdep_dsp_image_set_image = _lib.snd_hwdep_dsp_image_set_image
snd_hwdep_dsp_image_set_image.restype = None
snd_hwdep_dsp_image_set_image.argtypes = [POINTER(snd_hwdep_dsp_image_t), POINTER(None)]

# /usr/include/alsa/hwdep.h:151
snd_hwdep_dsp_image_set_length = _lib.snd_hwdep_dsp_image_set_length
snd_hwdep_dsp_image_set_length.restype = None
snd_hwdep_dsp_image_set_length.argtypes = [POINTER(snd_hwdep_dsp_image_t), c_size_t]

SND_CONTROL_DLSYM_VERSION = 0 	# /usr/include/alsa/control.h:43
class struct_snd_aes_iec958(Structure):
    __slots__ = [
        'status',
        'subcode',
        'pad',
        'dig_subframe',
    ]
struct_snd_aes_iec958._fields_ = [
    ('status', c_ubyte * 24),
    ('subcode', c_ubyte * 147),
    ('pad', c_ubyte),
    ('dig_subframe', c_ubyte * 4),
]

snd_aes_iec958_t = struct_snd_aes_iec958 	# /usr/include/alsa/control.h:51
class struct__snd_ctl_card_info(Structure):
    __slots__ = [
    ]
struct__snd_ctl_card_info._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_ctl_card_info(Structure):
    __slots__ = [
    ]
struct__snd_ctl_card_info._fields_ = [
    ('_opaque_struct', c_int)
]

snd_ctl_card_info_t = struct__snd_ctl_card_info 	# /usr/include/alsa/control.h:54
class struct__snd_ctl_elem_id(Structure):
    __slots__ = [
    ]
struct__snd_ctl_elem_id._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_ctl_elem_id(Structure):
    __slots__ = [
    ]
struct__snd_ctl_elem_id._fields_ = [
    ('_opaque_struct', c_int)
]

snd_ctl_elem_id_t = struct__snd_ctl_elem_id 	# /usr/include/alsa/control.h:57
class struct__snd_ctl_elem_list(Structure):
    __slots__ = [
    ]
struct__snd_ctl_elem_list._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_ctl_elem_list(Structure):
    __slots__ = [
    ]
struct__snd_ctl_elem_list._fields_ = [
    ('_opaque_struct', c_int)
]

snd_ctl_elem_list_t = struct__snd_ctl_elem_list 	# /usr/include/alsa/control.h:60
class struct__snd_ctl_elem_info(Structure):
    __slots__ = [
    ]
struct__snd_ctl_elem_info._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_ctl_elem_info(Structure):
    __slots__ = [
    ]
struct__snd_ctl_elem_info._fields_ = [
    ('_opaque_struct', c_int)
]

snd_ctl_elem_info_t = struct__snd_ctl_elem_info 	# /usr/include/alsa/control.h:63
class struct__snd_ctl_elem_value(Structure):
    __slots__ = [
    ]
struct__snd_ctl_elem_value._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_ctl_elem_value(Structure):
    __slots__ = [
    ]
struct__snd_ctl_elem_value._fields_ = [
    ('_opaque_struct', c_int)
]

snd_ctl_elem_value_t = struct__snd_ctl_elem_value 	# /usr/include/alsa/control.h:66
class struct__snd_ctl_event(Structure):
    __slots__ = [
    ]
struct__snd_ctl_event._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_ctl_event(Structure):
    __slots__ = [
    ]
struct__snd_ctl_event._fields_ = [
    ('_opaque_struct', c_int)
]

snd_ctl_event_t = struct__snd_ctl_event 	# /usr/include/alsa/control.h:69
enum__snd_ctl_elem_type = c_int
SND_CTL_ELEM_TYPE_NONE = 0
SND_CTL_ELEM_TYPE_BOOLEAN = 1
SND_CTL_ELEM_TYPE_INTEGER = 2
SND_CTL_ELEM_TYPE_ENUMERATED = 3
SND_CTL_ELEM_TYPE_BYTES = 4
SND_CTL_ELEM_TYPE_IEC958 = 5
SND_CTL_ELEM_TYPE_INTEGER64 = 6
SND_CTL_ELEM_TYPE_LAST = 0
snd_ctl_elem_type_t = enum__snd_ctl_elem_type 	# /usr/include/alsa/control.h:88
enum__snd_ctl_elem_iface = c_int
SND_CTL_ELEM_IFACE_CARD = 0
SND_CTL_ELEM_IFACE_HWDEP = 1
SND_CTL_ELEM_IFACE_MIXER = 2
SND_CTL_ELEM_IFACE_PCM = 3
SND_CTL_ELEM_IFACE_RAWMIDI = 4
SND_CTL_ELEM_IFACE_TIMER = 5
SND_CTL_ELEM_IFACE_SEQUENCER = 6
SND_CTL_ELEM_IFACE_LAST = 0
snd_ctl_elem_iface_t = enum__snd_ctl_elem_iface 	# /usr/include/alsa/control.h:107
enum__snd_ctl_event_type = c_int
SND_CTL_EVENT_ELEM = 0
SND_CTL_EVENT_LAST = 0
snd_ctl_event_type_t = enum__snd_ctl_event_type 	# /usr/include/alsa/control.h:114
SND_CTL_EVENT_MASK_REMOVE = -1 	# /usr/include/alsa/control.h:118
SND_CTL_EVENT_MASK_VALUE = 1 	# /usr/include/alsa/control.h:120
SND_CTL_EVENT_MASK_INFO = 2 	# /usr/include/alsa/control.h:122
SND_CTL_EVENT_MASK_ADD = 4 	# /usr/include/alsa/control.h:124
SND_CTL_EVENT_MASK_TLV = 8 	# /usr/include/alsa/control.h:126
SND_CTL_POWER_MASK = 65280 	# /usr/include/alsa/control.h:155
SND_CTL_POWER_D0 = 0 	# /usr/include/alsa/control.h:157
SND_CTL_POWER_D1 = 256 	# /usr/include/alsa/control.h:159
SND_CTL_POWER_D2 = 512 	# /usr/include/alsa/control.h:161
SND_CTL_POWER_D3 = 768 	# /usr/include/alsa/control.h:163
SND_CTL_POWER_D3hot = 768 	# /usr/include/alsa/control.h:165
SND_CTL_POWER_D3cold = 769 	# /usr/include/alsa/control.h:167
SND_CTL_TLVT_CONTAINER = 0 	# /usr/include/alsa/control.h:170
SND_CTL_TLVT_DB_SCALE = 1 	# /usr/include/alsa/control.h:172
SND_CTL_TLVT_DB_LINEAR = 2 	# /usr/include/alsa/control.h:174
SND_CTL_TLVT_DB_RANGE = 3 	# /usr/include/alsa/control.h:176
SND_CTL_TLV_DB_GAIN_MUTE = -9999999 	# /usr/include/alsa/control.h:179
enum__snd_ctl_type = c_int
SND_CTL_TYPE_HW = 1
SND_CTL_TYPE_SHM = 2
SND_CTL_TYPE_INET = 3
SND_CTL_TYPE_EXT = 4
snd_ctl_type_t = enum__snd_ctl_type 	# /usr/include/alsa/control.h:191
SND_CTL_NONBLOCK = 1 	# /usr/include/alsa/control.h:194
SND_CTL_ASYNC = 2 	# /usr/include/alsa/control.h:197
SND_CTL_READONLY = 4 	# /usr/include/alsa/control.h:200
class struct__snd_ctl(Structure):
    __slots__ = [
    ]
struct__snd_ctl._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_ctl(Structure):
    __slots__ = [
    ]
struct__snd_ctl._fields_ = [
    ('_opaque_struct', c_int)
]

snd_ctl_t = struct__snd_ctl 	# /usr/include/alsa/control.h:203
SND_SCTL_NOFREE = 1 	# /usr/include/alsa/control.h:206
class struct__snd_sctl(Structure):
    __slots__ = [
    ]
struct__snd_sctl._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_sctl(Structure):
    __slots__ = [
    ]
struct__snd_sctl._fields_ = [
    ('_opaque_struct', c_int)
]

snd_sctl_t = struct__snd_sctl 	# /usr/include/alsa/control.h:209
# /usr/include/alsa/control.h:211
snd_card_load = _lib.snd_card_load
snd_card_load.restype = c_int
snd_card_load.argtypes = [c_int]

# /usr/include/alsa/control.h:212
snd_card_next = _lib.snd_card_next
snd_card_next.restype = c_int
snd_card_next.argtypes = [POINTER(c_int)]

# /usr/include/alsa/control.h:213
snd_card_get_index = _lib.snd_card_get_index
snd_card_get_index.restype = c_int
snd_card_get_index.argtypes = [c_char_p]

# /usr/include/alsa/control.h:214
snd_card_get_name = _lib.snd_card_get_name
snd_card_get_name.restype = c_int
snd_card_get_name.argtypes = [c_int, POINTER(c_char_p)]

# /usr/include/alsa/control.h:215
snd_card_get_longname = _lib.snd_card_get_longname
snd_card_get_longname.restype = c_int
snd_card_get_longname.argtypes = [c_int, POINTER(c_char_p)]

''' Issue 144: These were added in 1.0.14
# /usr/include/alsa/control.h:217
snd_device_name_hint = _lib.snd_device_name_hint
snd_device_name_hint.restype = c_int
snd_device_name_hint.argtypes = [c_int, c_char_p, POINTER(POINTER(POINTER(None)))]

# /usr/include/alsa/control.h:218
snd_device_name_free_hint = _lib.snd_device_name_free_hint
snd_device_name_free_hint.restype = c_int
snd_device_name_free_hint.argtypes = [POINTER(POINTER(None))]

# /usr/include/alsa/control.h:219
snd_device_name_get_hint = _lib.snd_device_name_get_hint
snd_device_name_get_hint.restype = c_char_p
snd_device_name_get_hint.argtypes = [POINTER(None), c_char_p]
'''

# /usr/include/alsa/control.h:221
snd_ctl_open = _lib.snd_ctl_open
snd_ctl_open.restype = c_int
snd_ctl_open.argtypes = [POINTER(POINTER(snd_ctl_t)), c_char_p, c_int]

# /usr/include/alsa/control.h:222
snd_ctl_open_lconf = _lib.snd_ctl_open_lconf
snd_ctl_open_lconf.restype = c_int
snd_ctl_open_lconf.argtypes = [POINTER(POINTER(snd_ctl_t)), c_char_p, c_int, POINTER(snd_config_t)]

# /usr/include/alsa/control.h:223
snd_ctl_close = _lib.snd_ctl_close
snd_ctl_close.restype = c_int
snd_ctl_close.argtypes = [POINTER(snd_ctl_t)]

# /usr/include/alsa/control.h:224
snd_ctl_nonblock = _lib.snd_ctl_nonblock
snd_ctl_nonblock.restype = c_int
snd_ctl_nonblock.argtypes = [POINTER(snd_ctl_t), c_int]

# /usr/include/alsa/control.h:225
snd_async_add_ctl_handler = _lib.snd_async_add_ctl_handler
snd_async_add_ctl_handler.restype = c_int
snd_async_add_ctl_handler.argtypes = [POINTER(POINTER(snd_async_handler_t)), POINTER(snd_ctl_t), snd_async_callback_t, POINTER(None)]

# /usr/include/alsa/control.h:227
snd_async_handler_get_ctl = _lib.snd_async_handler_get_ctl
snd_async_handler_get_ctl.restype = POINTER(snd_ctl_t)
snd_async_handler_get_ctl.argtypes = [POINTER(snd_async_handler_t)]

# /usr/include/alsa/control.h:228
snd_ctl_poll_descriptors_count = _lib.snd_ctl_poll_descriptors_count
snd_ctl_poll_descriptors_count.restype = c_int
snd_ctl_poll_descriptors_count.argtypes = [POINTER(snd_ctl_t)]

class struct_pollfd(Structure):
    __slots__ = [
    ]
struct_pollfd._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/control.h:229
snd_ctl_poll_descriptors = _lib.snd_ctl_poll_descriptors
snd_ctl_poll_descriptors.restype = c_int
snd_ctl_poll_descriptors.argtypes = [POINTER(snd_ctl_t), POINTER(struct_pollfd), c_uint]

class struct_pollfd(Structure):
    __slots__ = [
    ]
struct_pollfd._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/control.h:230
snd_ctl_poll_descriptors_revents = _lib.snd_ctl_poll_descriptors_revents
snd_ctl_poll_descriptors_revents.restype = c_int
snd_ctl_poll_descriptors_revents.argtypes = [POINTER(snd_ctl_t), POINTER(struct_pollfd), c_uint, POINTER(c_ushort)]

# /usr/include/alsa/control.h:231
snd_ctl_subscribe_events = _lib.snd_ctl_subscribe_events
snd_ctl_subscribe_events.restype = c_int
snd_ctl_subscribe_events.argtypes = [POINTER(snd_ctl_t), c_int]

# /usr/include/alsa/control.h:232
snd_ctl_card_info = _lib.snd_ctl_card_info
snd_ctl_card_info.restype = c_int
snd_ctl_card_info.argtypes = [POINTER(snd_ctl_t), POINTER(snd_ctl_card_info_t)]

# /usr/include/alsa/control.h:233
snd_ctl_elem_list = _lib.snd_ctl_elem_list
snd_ctl_elem_list.restype = c_int
snd_ctl_elem_list.argtypes = [POINTER(snd_ctl_t), POINTER(snd_ctl_elem_list_t)]

# /usr/include/alsa/control.h:234
snd_ctl_elem_info = _lib.snd_ctl_elem_info
snd_ctl_elem_info.restype = c_int
snd_ctl_elem_info.argtypes = [POINTER(snd_ctl_t), POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:235
snd_ctl_elem_read = _lib.snd_ctl_elem_read
snd_ctl_elem_read.restype = c_int
snd_ctl_elem_read.argtypes = [POINTER(snd_ctl_t), POINTER(snd_ctl_elem_value_t)]

# /usr/include/alsa/control.h:236
snd_ctl_elem_write = _lib.snd_ctl_elem_write
snd_ctl_elem_write.restype = c_int
snd_ctl_elem_write.argtypes = [POINTER(snd_ctl_t), POINTER(snd_ctl_elem_value_t)]

# /usr/include/alsa/control.h:237
snd_ctl_elem_lock = _lib.snd_ctl_elem_lock
snd_ctl_elem_lock.restype = c_int
snd_ctl_elem_lock.argtypes = [POINTER(snd_ctl_t), POINTER(snd_ctl_elem_id_t)]

# /usr/include/alsa/control.h:238
snd_ctl_elem_unlock = _lib.snd_ctl_elem_unlock
snd_ctl_elem_unlock.restype = c_int
snd_ctl_elem_unlock.argtypes = [POINTER(snd_ctl_t), POINTER(snd_ctl_elem_id_t)]

# /usr/include/alsa/control.h:239
snd_ctl_elem_tlv_read = _lib.snd_ctl_elem_tlv_read
snd_ctl_elem_tlv_read.restype = c_int
snd_ctl_elem_tlv_read.argtypes = [POINTER(snd_ctl_t), POINTER(snd_ctl_elem_id_t), POINTER(c_uint), c_uint]

# /usr/include/alsa/control.h:241
snd_ctl_elem_tlv_write = _lib.snd_ctl_elem_tlv_write
snd_ctl_elem_tlv_write.restype = c_int
snd_ctl_elem_tlv_write.argtypes = [POINTER(snd_ctl_t), POINTER(snd_ctl_elem_id_t), POINTER(c_uint)]

# /usr/include/alsa/control.h:243
snd_ctl_elem_tlv_command = _lib.snd_ctl_elem_tlv_command
snd_ctl_elem_tlv_command.restype = c_int
snd_ctl_elem_tlv_command.argtypes = [POINTER(snd_ctl_t), POINTER(snd_ctl_elem_id_t), POINTER(c_uint)]

# /usr/include/alsa/control.h:245
snd_ctl_hwdep_next_device = _lib.snd_ctl_hwdep_next_device
snd_ctl_hwdep_next_device.restype = c_int
snd_ctl_hwdep_next_device.argtypes = [POINTER(snd_ctl_t), POINTER(c_int)]

# /usr/include/alsa/control.h:246
snd_ctl_hwdep_info = _lib.snd_ctl_hwdep_info
snd_ctl_hwdep_info.restype = c_int
snd_ctl_hwdep_info.argtypes = [POINTER(snd_ctl_t), POINTER(snd_hwdep_info_t)]

# /usr/include/alsa/control.h:247
snd_ctl_pcm_next_device = _lib.snd_ctl_pcm_next_device
snd_ctl_pcm_next_device.restype = c_int
snd_ctl_pcm_next_device.argtypes = [POINTER(snd_ctl_t), POINTER(c_int)]

# /usr/include/alsa/control.h:248
snd_ctl_pcm_info = _lib.snd_ctl_pcm_info
snd_ctl_pcm_info.restype = c_int
snd_ctl_pcm_info.argtypes = [POINTER(snd_ctl_t), POINTER(snd_pcm_info_t)]

# /usr/include/alsa/control.h:249
snd_ctl_pcm_prefer_subdevice = _lib.snd_ctl_pcm_prefer_subdevice
snd_ctl_pcm_prefer_subdevice.restype = c_int
snd_ctl_pcm_prefer_subdevice.argtypes = [POINTER(snd_ctl_t), c_int]

# /usr/include/alsa/control.h:250
snd_ctl_rawmidi_next_device = _lib.snd_ctl_rawmidi_next_device
snd_ctl_rawmidi_next_device.restype = c_int
snd_ctl_rawmidi_next_device.argtypes = [POINTER(snd_ctl_t), POINTER(c_int)]

# /usr/include/alsa/control.h:251
snd_ctl_rawmidi_info = _lib.snd_ctl_rawmidi_info
snd_ctl_rawmidi_info.restype = c_int
snd_ctl_rawmidi_info.argtypes = [POINTER(snd_ctl_t), POINTER(snd_rawmidi_info_t)]

# /usr/include/alsa/control.h:252
snd_ctl_rawmidi_prefer_subdevice = _lib.snd_ctl_rawmidi_prefer_subdevice
snd_ctl_rawmidi_prefer_subdevice.restype = c_int
snd_ctl_rawmidi_prefer_subdevice.argtypes = [POINTER(snd_ctl_t), c_int]

# /usr/include/alsa/control.h:253
snd_ctl_set_power_state = _lib.snd_ctl_set_power_state
snd_ctl_set_power_state.restype = c_int
snd_ctl_set_power_state.argtypes = [POINTER(snd_ctl_t), c_uint]

# /usr/include/alsa/control.h:254
snd_ctl_get_power_state = _lib.snd_ctl_get_power_state
snd_ctl_get_power_state.restype = c_int
snd_ctl_get_power_state.argtypes = [POINTER(snd_ctl_t), POINTER(c_uint)]

# /usr/include/alsa/control.h:256
snd_ctl_read = _lib.snd_ctl_read
snd_ctl_read.restype = c_int
snd_ctl_read.argtypes = [POINTER(snd_ctl_t), POINTER(snd_ctl_event_t)]

# /usr/include/alsa/control.h:257
snd_ctl_wait = _lib.snd_ctl_wait
snd_ctl_wait.restype = c_int
snd_ctl_wait.argtypes = [POINTER(snd_ctl_t), c_int]

# /usr/include/alsa/control.h:258
snd_ctl_name = _lib.snd_ctl_name
snd_ctl_name.restype = c_char_p
snd_ctl_name.argtypes = [POINTER(snd_ctl_t)]

# /usr/include/alsa/control.h:259
snd_ctl_type = _lib.snd_ctl_type
snd_ctl_type.restype = snd_ctl_type_t
snd_ctl_type.argtypes = [POINTER(snd_ctl_t)]

# /usr/include/alsa/control.h:261
snd_ctl_elem_type_name = _lib.snd_ctl_elem_type_name
snd_ctl_elem_type_name.restype = c_char_p
snd_ctl_elem_type_name.argtypes = [snd_ctl_elem_type_t]

# /usr/include/alsa/control.h:262
snd_ctl_elem_iface_name = _lib.snd_ctl_elem_iface_name
snd_ctl_elem_iface_name.restype = c_char_p
snd_ctl_elem_iface_name.argtypes = [snd_ctl_elem_iface_t]

# /usr/include/alsa/control.h:263
snd_ctl_event_type_name = _lib.snd_ctl_event_type_name
snd_ctl_event_type_name.restype = c_char_p
snd_ctl_event_type_name.argtypes = [snd_ctl_event_type_t]

# /usr/include/alsa/control.h:265
snd_ctl_event_elem_get_mask = _lib.snd_ctl_event_elem_get_mask
snd_ctl_event_elem_get_mask.restype = c_uint
snd_ctl_event_elem_get_mask.argtypes = [POINTER(snd_ctl_event_t)]

# /usr/include/alsa/control.h:266
snd_ctl_event_elem_get_numid = _lib.snd_ctl_event_elem_get_numid
snd_ctl_event_elem_get_numid.restype = c_uint
snd_ctl_event_elem_get_numid.argtypes = [POINTER(snd_ctl_event_t)]

# /usr/include/alsa/control.h:267
snd_ctl_event_elem_get_id = _lib.snd_ctl_event_elem_get_id
snd_ctl_event_elem_get_id.restype = None
snd_ctl_event_elem_get_id.argtypes = [POINTER(snd_ctl_event_t), POINTER(snd_ctl_elem_id_t)]

# /usr/include/alsa/control.h:268
snd_ctl_event_elem_get_interface = _lib.snd_ctl_event_elem_get_interface
snd_ctl_event_elem_get_interface.restype = snd_ctl_elem_iface_t
snd_ctl_event_elem_get_interface.argtypes = [POINTER(snd_ctl_event_t)]

# /usr/include/alsa/control.h:269
snd_ctl_event_elem_get_device = _lib.snd_ctl_event_elem_get_device
snd_ctl_event_elem_get_device.restype = c_uint
snd_ctl_event_elem_get_device.argtypes = [POINTER(snd_ctl_event_t)]

# /usr/include/alsa/control.h:270
snd_ctl_event_elem_get_subdevice = _lib.snd_ctl_event_elem_get_subdevice
snd_ctl_event_elem_get_subdevice.restype = c_uint
snd_ctl_event_elem_get_subdevice.argtypes = [POINTER(snd_ctl_event_t)]

# /usr/include/alsa/control.h:271
snd_ctl_event_elem_get_name = _lib.snd_ctl_event_elem_get_name
snd_ctl_event_elem_get_name.restype = c_char_p
snd_ctl_event_elem_get_name.argtypes = [POINTER(snd_ctl_event_t)]

# /usr/include/alsa/control.h:272
snd_ctl_event_elem_get_index = _lib.snd_ctl_event_elem_get_index
snd_ctl_event_elem_get_index.restype = c_uint
snd_ctl_event_elem_get_index.argtypes = [POINTER(snd_ctl_event_t)]

# /usr/include/alsa/control.h:274
snd_ctl_elem_list_alloc_space = _lib.snd_ctl_elem_list_alloc_space
snd_ctl_elem_list_alloc_space.restype = c_int
snd_ctl_elem_list_alloc_space.argtypes = [POINTER(snd_ctl_elem_list_t), c_uint]

# /usr/include/alsa/control.h:275
snd_ctl_elem_list_free_space = _lib.snd_ctl_elem_list_free_space
snd_ctl_elem_list_free_space.restype = None
snd_ctl_elem_list_free_space.argtypes = [POINTER(snd_ctl_elem_list_t)]

# /usr/include/alsa/control.h:277
snd_ctl_elem_id_sizeof = _lib.snd_ctl_elem_id_sizeof
snd_ctl_elem_id_sizeof.restype = c_size_t
snd_ctl_elem_id_sizeof.argtypes = []

# /usr/include/alsa/control.h:283
snd_ctl_elem_id_malloc = _lib.snd_ctl_elem_id_malloc
snd_ctl_elem_id_malloc.restype = c_int
snd_ctl_elem_id_malloc.argtypes = [POINTER(POINTER(snd_ctl_elem_id_t))]

# /usr/include/alsa/control.h:284
snd_ctl_elem_id_free = _lib.snd_ctl_elem_id_free
snd_ctl_elem_id_free.restype = None
snd_ctl_elem_id_free.argtypes = [POINTER(snd_ctl_elem_id_t)]

# /usr/include/alsa/control.h:285
snd_ctl_elem_id_clear = _lib.snd_ctl_elem_id_clear
snd_ctl_elem_id_clear.restype = None
snd_ctl_elem_id_clear.argtypes = [POINTER(snd_ctl_elem_id_t)]

# /usr/include/alsa/control.h:286
snd_ctl_elem_id_copy = _lib.snd_ctl_elem_id_copy
snd_ctl_elem_id_copy.restype = None
snd_ctl_elem_id_copy.argtypes = [POINTER(snd_ctl_elem_id_t), POINTER(snd_ctl_elem_id_t)]

# /usr/include/alsa/control.h:287
snd_ctl_elem_id_get_numid = _lib.snd_ctl_elem_id_get_numid
snd_ctl_elem_id_get_numid.restype = c_uint
snd_ctl_elem_id_get_numid.argtypes = [POINTER(snd_ctl_elem_id_t)]

# /usr/include/alsa/control.h:288
snd_ctl_elem_id_get_interface = _lib.snd_ctl_elem_id_get_interface
snd_ctl_elem_id_get_interface.restype = snd_ctl_elem_iface_t
snd_ctl_elem_id_get_interface.argtypes = [POINTER(snd_ctl_elem_id_t)]

# /usr/include/alsa/control.h:289
snd_ctl_elem_id_get_device = _lib.snd_ctl_elem_id_get_device
snd_ctl_elem_id_get_device.restype = c_uint
snd_ctl_elem_id_get_device.argtypes = [POINTER(snd_ctl_elem_id_t)]

# /usr/include/alsa/control.h:290
snd_ctl_elem_id_get_subdevice = _lib.snd_ctl_elem_id_get_subdevice
snd_ctl_elem_id_get_subdevice.restype = c_uint
snd_ctl_elem_id_get_subdevice.argtypes = [POINTER(snd_ctl_elem_id_t)]

# /usr/include/alsa/control.h:291
snd_ctl_elem_id_get_name = _lib.snd_ctl_elem_id_get_name
snd_ctl_elem_id_get_name.restype = c_char_p
snd_ctl_elem_id_get_name.argtypes = [POINTER(snd_ctl_elem_id_t)]

# /usr/include/alsa/control.h:292
snd_ctl_elem_id_get_index = _lib.snd_ctl_elem_id_get_index
snd_ctl_elem_id_get_index.restype = c_uint
snd_ctl_elem_id_get_index.argtypes = [POINTER(snd_ctl_elem_id_t)]

# /usr/include/alsa/control.h:293
snd_ctl_elem_id_set_numid = _lib.snd_ctl_elem_id_set_numid
snd_ctl_elem_id_set_numid.restype = None
snd_ctl_elem_id_set_numid.argtypes = [POINTER(snd_ctl_elem_id_t), c_uint]

# /usr/include/alsa/control.h:294
snd_ctl_elem_id_set_interface = _lib.snd_ctl_elem_id_set_interface
snd_ctl_elem_id_set_interface.restype = None
snd_ctl_elem_id_set_interface.argtypes = [POINTER(snd_ctl_elem_id_t), snd_ctl_elem_iface_t]

# /usr/include/alsa/control.h:295
snd_ctl_elem_id_set_device = _lib.snd_ctl_elem_id_set_device
snd_ctl_elem_id_set_device.restype = None
snd_ctl_elem_id_set_device.argtypes = [POINTER(snd_ctl_elem_id_t), c_uint]

# /usr/include/alsa/control.h:296
snd_ctl_elem_id_set_subdevice = _lib.snd_ctl_elem_id_set_subdevice
snd_ctl_elem_id_set_subdevice.restype = None
snd_ctl_elem_id_set_subdevice.argtypes = [POINTER(snd_ctl_elem_id_t), c_uint]

# /usr/include/alsa/control.h:297
snd_ctl_elem_id_set_name = _lib.snd_ctl_elem_id_set_name
snd_ctl_elem_id_set_name.restype = None
snd_ctl_elem_id_set_name.argtypes = [POINTER(snd_ctl_elem_id_t), c_char_p]

# /usr/include/alsa/control.h:298
snd_ctl_elem_id_set_index = _lib.snd_ctl_elem_id_set_index
snd_ctl_elem_id_set_index.restype = None
snd_ctl_elem_id_set_index.argtypes = [POINTER(snd_ctl_elem_id_t), c_uint]

# /usr/include/alsa/control.h:300
snd_ctl_card_info_sizeof = _lib.snd_ctl_card_info_sizeof
snd_ctl_card_info_sizeof.restype = c_size_t
snd_ctl_card_info_sizeof.argtypes = []

# /usr/include/alsa/control.h:306
snd_ctl_card_info_malloc = _lib.snd_ctl_card_info_malloc
snd_ctl_card_info_malloc.restype = c_int
snd_ctl_card_info_malloc.argtypes = [POINTER(POINTER(snd_ctl_card_info_t))]

# /usr/include/alsa/control.h:307
snd_ctl_card_info_free = _lib.snd_ctl_card_info_free
snd_ctl_card_info_free.restype = None
snd_ctl_card_info_free.argtypes = [POINTER(snd_ctl_card_info_t)]

# /usr/include/alsa/control.h:308
snd_ctl_card_info_clear = _lib.snd_ctl_card_info_clear
snd_ctl_card_info_clear.restype = None
snd_ctl_card_info_clear.argtypes = [POINTER(snd_ctl_card_info_t)]

# /usr/include/alsa/control.h:309
snd_ctl_card_info_copy = _lib.snd_ctl_card_info_copy
snd_ctl_card_info_copy.restype = None
snd_ctl_card_info_copy.argtypes = [POINTER(snd_ctl_card_info_t), POINTER(snd_ctl_card_info_t)]

# /usr/include/alsa/control.h:310
snd_ctl_card_info_get_card = _lib.snd_ctl_card_info_get_card
snd_ctl_card_info_get_card.restype = c_int
snd_ctl_card_info_get_card.argtypes = [POINTER(snd_ctl_card_info_t)]

# /usr/include/alsa/control.h:311
snd_ctl_card_info_get_id = _lib.snd_ctl_card_info_get_id
snd_ctl_card_info_get_id.restype = c_char_p
snd_ctl_card_info_get_id.argtypes = [POINTER(snd_ctl_card_info_t)]

# /usr/include/alsa/control.h:312
snd_ctl_card_info_get_driver = _lib.snd_ctl_card_info_get_driver
snd_ctl_card_info_get_driver.restype = c_char_p
snd_ctl_card_info_get_driver.argtypes = [POINTER(snd_ctl_card_info_t)]

# /usr/include/alsa/control.h:313
snd_ctl_card_info_get_name = _lib.snd_ctl_card_info_get_name
snd_ctl_card_info_get_name.restype = c_char_p
snd_ctl_card_info_get_name.argtypes = [POINTER(snd_ctl_card_info_t)]

# /usr/include/alsa/control.h:314
snd_ctl_card_info_get_longname = _lib.snd_ctl_card_info_get_longname
snd_ctl_card_info_get_longname.restype = c_char_p
snd_ctl_card_info_get_longname.argtypes = [POINTER(snd_ctl_card_info_t)]

# /usr/include/alsa/control.h:315
snd_ctl_card_info_get_mixername = _lib.snd_ctl_card_info_get_mixername
snd_ctl_card_info_get_mixername.restype = c_char_p
snd_ctl_card_info_get_mixername.argtypes = [POINTER(snd_ctl_card_info_t)]

# /usr/include/alsa/control.h:316
snd_ctl_card_info_get_components = _lib.snd_ctl_card_info_get_components
snd_ctl_card_info_get_components.restype = c_char_p
snd_ctl_card_info_get_components.argtypes = [POINTER(snd_ctl_card_info_t)]

# /usr/include/alsa/control.h:318
snd_ctl_event_sizeof = _lib.snd_ctl_event_sizeof
snd_ctl_event_sizeof.restype = c_size_t
snd_ctl_event_sizeof.argtypes = []

# /usr/include/alsa/control.h:324
snd_ctl_event_malloc = _lib.snd_ctl_event_malloc
snd_ctl_event_malloc.restype = c_int
snd_ctl_event_malloc.argtypes = [POINTER(POINTER(snd_ctl_event_t))]

# /usr/include/alsa/control.h:325
snd_ctl_event_free = _lib.snd_ctl_event_free
snd_ctl_event_free.restype = None
snd_ctl_event_free.argtypes = [POINTER(snd_ctl_event_t)]

# /usr/include/alsa/control.h:326
snd_ctl_event_clear = _lib.snd_ctl_event_clear
snd_ctl_event_clear.restype = None
snd_ctl_event_clear.argtypes = [POINTER(snd_ctl_event_t)]

# /usr/include/alsa/control.h:327
snd_ctl_event_copy = _lib.snd_ctl_event_copy
snd_ctl_event_copy.restype = None
snd_ctl_event_copy.argtypes = [POINTER(snd_ctl_event_t), POINTER(snd_ctl_event_t)]

# /usr/include/alsa/control.h:328
snd_ctl_event_get_type = _lib.snd_ctl_event_get_type
snd_ctl_event_get_type.restype = snd_ctl_event_type_t
snd_ctl_event_get_type.argtypes = [POINTER(snd_ctl_event_t)]

# /usr/include/alsa/control.h:330
snd_ctl_elem_list_sizeof = _lib.snd_ctl_elem_list_sizeof
snd_ctl_elem_list_sizeof.restype = c_size_t
snd_ctl_elem_list_sizeof.argtypes = []

# /usr/include/alsa/control.h:336
snd_ctl_elem_list_malloc = _lib.snd_ctl_elem_list_malloc
snd_ctl_elem_list_malloc.restype = c_int
snd_ctl_elem_list_malloc.argtypes = [POINTER(POINTER(snd_ctl_elem_list_t))]

# /usr/include/alsa/control.h:337
snd_ctl_elem_list_free = _lib.snd_ctl_elem_list_free
snd_ctl_elem_list_free.restype = None
snd_ctl_elem_list_free.argtypes = [POINTER(snd_ctl_elem_list_t)]

# /usr/include/alsa/control.h:338
snd_ctl_elem_list_clear = _lib.snd_ctl_elem_list_clear
snd_ctl_elem_list_clear.restype = None
snd_ctl_elem_list_clear.argtypes = [POINTER(snd_ctl_elem_list_t)]

# /usr/include/alsa/control.h:339
snd_ctl_elem_list_copy = _lib.snd_ctl_elem_list_copy
snd_ctl_elem_list_copy.restype = None
snd_ctl_elem_list_copy.argtypes = [POINTER(snd_ctl_elem_list_t), POINTER(snd_ctl_elem_list_t)]

# /usr/include/alsa/control.h:340
snd_ctl_elem_list_set_offset = _lib.snd_ctl_elem_list_set_offset
snd_ctl_elem_list_set_offset.restype = None
snd_ctl_elem_list_set_offset.argtypes = [POINTER(snd_ctl_elem_list_t), c_uint]

# /usr/include/alsa/control.h:341
snd_ctl_elem_list_get_used = _lib.snd_ctl_elem_list_get_used
snd_ctl_elem_list_get_used.restype = c_uint
snd_ctl_elem_list_get_used.argtypes = [POINTER(snd_ctl_elem_list_t)]

# /usr/include/alsa/control.h:342
snd_ctl_elem_list_get_count = _lib.snd_ctl_elem_list_get_count
snd_ctl_elem_list_get_count.restype = c_uint
snd_ctl_elem_list_get_count.argtypes = [POINTER(snd_ctl_elem_list_t)]

# /usr/include/alsa/control.h:343
snd_ctl_elem_list_get_id = _lib.snd_ctl_elem_list_get_id
snd_ctl_elem_list_get_id.restype = None
snd_ctl_elem_list_get_id.argtypes = [POINTER(snd_ctl_elem_list_t), c_uint, POINTER(snd_ctl_elem_id_t)]

# /usr/include/alsa/control.h:344
snd_ctl_elem_list_get_numid = _lib.snd_ctl_elem_list_get_numid
snd_ctl_elem_list_get_numid.restype = c_uint
snd_ctl_elem_list_get_numid.argtypes = [POINTER(snd_ctl_elem_list_t), c_uint]

# /usr/include/alsa/control.h:345
snd_ctl_elem_list_get_interface = _lib.snd_ctl_elem_list_get_interface
snd_ctl_elem_list_get_interface.restype = snd_ctl_elem_iface_t
snd_ctl_elem_list_get_interface.argtypes = [POINTER(snd_ctl_elem_list_t), c_uint]

# /usr/include/alsa/control.h:346
snd_ctl_elem_list_get_device = _lib.snd_ctl_elem_list_get_device
snd_ctl_elem_list_get_device.restype = c_uint
snd_ctl_elem_list_get_device.argtypes = [POINTER(snd_ctl_elem_list_t), c_uint]

# /usr/include/alsa/control.h:347
snd_ctl_elem_list_get_subdevice = _lib.snd_ctl_elem_list_get_subdevice
snd_ctl_elem_list_get_subdevice.restype = c_uint
snd_ctl_elem_list_get_subdevice.argtypes = [POINTER(snd_ctl_elem_list_t), c_uint]

# /usr/include/alsa/control.h:348
snd_ctl_elem_list_get_name = _lib.snd_ctl_elem_list_get_name
snd_ctl_elem_list_get_name.restype = c_char_p
snd_ctl_elem_list_get_name.argtypes = [POINTER(snd_ctl_elem_list_t), c_uint]

# /usr/include/alsa/control.h:349
snd_ctl_elem_list_get_index = _lib.snd_ctl_elem_list_get_index
snd_ctl_elem_list_get_index.restype = c_uint
snd_ctl_elem_list_get_index.argtypes = [POINTER(snd_ctl_elem_list_t), c_uint]

# /usr/include/alsa/control.h:351
snd_ctl_elem_info_sizeof = _lib.snd_ctl_elem_info_sizeof
snd_ctl_elem_info_sizeof.restype = c_size_t
snd_ctl_elem_info_sizeof.argtypes = []

# /usr/include/alsa/control.h:357
snd_ctl_elem_info_malloc = _lib.snd_ctl_elem_info_malloc
snd_ctl_elem_info_malloc.restype = c_int
snd_ctl_elem_info_malloc.argtypes = [POINTER(POINTER(snd_ctl_elem_info_t))]

# /usr/include/alsa/control.h:358
snd_ctl_elem_info_free = _lib.snd_ctl_elem_info_free
snd_ctl_elem_info_free.restype = None
snd_ctl_elem_info_free.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:359
snd_ctl_elem_info_clear = _lib.snd_ctl_elem_info_clear
snd_ctl_elem_info_clear.restype = None
snd_ctl_elem_info_clear.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:360
snd_ctl_elem_info_copy = _lib.snd_ctl_elem_info_copy
snd_ctl_elem_info_copy.restype = None
snd_ctl_elem_info_copy.argtypes = [POINTER(snd_ctl_elem_info_t), POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:361
snd_ctl_elem_info_get_type = _lib.snd_ctl_elem_info_get_type
snd_ctl_elem_info_get_type.restype = snd_ctl_elem_type_t
snd_ctl_elem_info_get_type.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:362
snd_ctl_elem_info_is_readable = _lib.snd_ctl_elem_info_is_readable
snd_ctl_elem_info_is_readable.restype = c_int
snd_ctl_elem_info_is_readable.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:363
snd_ctl_elem_info_is_writable = _lib.snd_ctl_elem_info_is_writable
snd_ctl_elem_info_is_writable.restype = c_int
snd_ctl_elem_info_is_writable.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:364
snd_ctl_elem_info_is_volatile = _lib.snd_ctl_elem_info_is_volatile
snd_ctl_elem_info_is_volatile.restype = c_int
snd_ctl_elem_info_is_volatile.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:365
snd_ctl_elem_info_is_inactive = _lib.snd_ctl_elem_info_is_inactive
snd_ctl_elem_info_is_inactive.restype = c_int
snd_ctl_elem_info_is_inactive.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:366
snd_ctl_elem_info_is_locked = _lib.snd_ctl_elem_info_is_locked
snd_ctl_elem_info_is_locked.restype = c_int
snd_ctl_elem_info_is_locked.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:367
snd_ctl_elem_info_is_tlv_readable = _lib.snd_ctl_elem_info_is_tlv_readable
snd_ctl_elem_info_is_tlv_readable.restype = c_int
snd_ctl_elem_info_is_tlv_readable.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:368
snd_ctl_elem_info_is_tlv_writable = _lib.snd_ctl_elem_info_is_tlv_writable
snd_ctl_elem_info_is_tlv_writable.restype = c_int
snd_ctl_elem_info_is_tlv_writable.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:369
snd_ctl_elem_info_is_tlv_commandable = _lib.snd_ctl_elem_info_is_tlv_commandable
snd_ctl_elem_info_is_tlv_commandable.restype = c_int
snd_ctl_elem_info_is_tlv_commandable.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:370
snd_ctl_elem_info_is_owner = _lib.snd_ctl_elem_info_is_owner
snd_ctl_elem_info_is_owner.restype = c_int
snd_ctl_elem_info_is_owner.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:371
snd_ctl_elem_info_is_user = _lib.snd_ctl_elem_info_is_user
snd_ctl_elem_info_is_user.restype = c_int
snd_ctl_elem_info_is_user.argtypes = [POINTER(snd_ctl_elem_info_t)]

__pid_t = c_int 	# /usr/include/gentoo-multilib/amd64/bits/types.h:146
pid_t = __pid_t 	# /usr/include/gentoo-multilib/amd64/unistd.h:229
# /usr/include/alsa/control.h:372
snd_ctl_elem_info_get_owner = _lib.snd_ctl_elem_info_get_owner
snd_ctl_elem_info_get_owner.restype = pid_t
snd_ctl_elem_info_get_owner.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:373
snd_ctl_elem_info_get_count = _lib.snd_ctl_elem_info_get_count
snd_ctl_elem_info_get_count.restype = c_uint
snd_ctl_elem_info_get_count.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:374
snd_ctl_elem_info_get_min = _lib.snd_ctl_elem_info_get_min
snd_ctl_elem_info_get_min.restype = c_long
snd_ctl_elem_info_get_min.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:375
snd_ctl_elem_info_get_max = _lib.snd_ctl_elem_info_get_max
snd_ctl_elem_info_get_max.restype = c_long
snd_ctl_elem_info_get_max.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:376
snd_ctl_elem_info_get_step = _lib.snd_ctl_elem_info_get_step
snd_ctl_elem_info_get_step.restype = c_long
snd_ctl_elem_info_get_step.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:377
snd_ctl_elem_info_get_min64 = _lib.snd_ctl_elem_info_get_min64
snd_ctl_elem_info_get_min64.restype = c_longlong
snd_ctl_elem_info_get_min64.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:378
snd_ctl_elem_info_get_max64 = _lib.snd_ctl_elem_info_get_max64
snd_ctl_elem_info_get_max64.restype = c_longlong
snd_ctl_elem_info_get_max64.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:379
snd_ctl_elem_info_get_step64 = _lib.snd_ctl_elem_info_get_step64
snd_ctl_elem_info_get_step64.restype = c_longlong
snd_ctl_elem_info_get_step64.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:380
snd_ctl_elem_info_get_items = _lib.snd_ctl_elem_info_get_items
snd_ctl_elem_info_get_items.restype = c_uint
snd_ctl_elem_info_get_items.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:381
snd_ctl_elem_info_set_item = _lib.snd_ctl_elem_info_set_item
snd_ctl_elem_info_set_item.restype = None
snd_ctl_elem_info_set_item.argtypes = [POINTER(snd_ctl_elem_info_t), c_uint]

# /usr/include/alsa/control.h:382
snd_ctl_elem_info_get_item_name = _lib.snd_ctl_elem_info_get_item_name
snd_ctl_elem_info_get_item_name.restype = c_char_p
snd_ctl_elem_info_get_item_name.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:383
snd_ctl_elem_info_get_dimensions = _lib.snd_ctl_elem_info_get_dimensions
snd_ctl_elem_info_get_dimensions.restype = c_int
snd_ctl_elem_info_get_dimensions.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:384
snd_ctl_elem_info_get_dimension = _lib.snd_ctl_elem_info_get_dimension
snd_ctl_elem_info_get_dimension.restype = c_int
snd_ctl_elem_info_get_dimension.argtypes = [POINTER(snd_ctl_elem_info_t), c_uint]

# /usr/include/alsa/control.h:385
snd_ctl_elem_info_get_id = _lib.snd_ctl_elem_info_get_id
snd_ctl_elem_info_get_id.restype = None
snd_ctl_elem_info_get_id.argtypes = [POINTER(snd_ctl_elem_info_t), POINTER(snd_ctl_elem_id_t)]

# /usr/include/alsa/control.h:386
snd_ctl_elem_info_get_numid = _lib.snd_ctl_elem_info_get_numid
snd_ctl_elem_info_get_numid.restype = c_uint
snd_ctl_elem_info_get_numid.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:387
snd_ctl_elem_info_get_interface = _lib.snd_ctl_elem_info_get_interface
snd_ctl_elem_info_get_interface.restype = snd_ctl_elem_iface_t
snd_ctl_elem_info_get_interface.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:388
snd_ctl_elem_info_get_device = _lib.snd_ctl_elem_info_get_device
snd_ctl_elem_info_get_device.restype = c_uint
snd_ctl_elem_info_get_device.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:389
snd_ctl_elem_info_get_subdevice = _lib.snd_ctl_elem_info_get_subdevice
snd_ctl_elem_info_get_subdevice.restype = c_uint
snd_ctl_elem_info_get_subdevice.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:390
snd_ctl_elem_info_get_name = _lib.snd_ctl_elem_info_get_name
snd_ctl_elem_info_get_name.restype = c_char_p
snd_ctl_elem_info_get_name.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:391
snd_ctl_elem_info_get_index = _lib.snd_ctl_elem_info_get_index
snd_ctl_elem_info_get_index.restype = c_uint
snd_ctl_elem_info_get_index.argtypes = [POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:392
snd_ctl_elem_info_set_id = _lib.snd_ctl_elem_info_set_id
snd_ctl_elem_info_set_id.restype = None
snd_ctl_elem_info_set_id.argtypes = [POINTER(snd_ctl_elem_info_t), POINTER(snd_ctl_elem_id_t)]

# /usr/include/alsa/control.h:393
snd_ctl_elem_info_set_numid = _lib.snd_ctl_elem_info_set_numid
snd_ctl_elem_info_set_numid.restype = None
snd_ctl_elem_info_set_numid.argtypes = [POINTER(snd_ctl_elem_info_t), c_uint]

# /usr/include/alsa/control.h:394
snd_ctl_elem_info_set_interface = _lib.snd_ctl_elem_info_set_interface
snd_ctl_elem_info_set_interface.restype = None
snd_ctl_elem_info_set_interface.argtypes = [POINTER(snd_ctl_elem_info_t), snd_ctl_elem_iface_t]

# /usr/include/alsa/control.h:395
snd_ctl_elem_info_set_device = _lib.snd_ctl_elem_info_set_device
snd_ctl_elem_info_set_device.restype = None
snd_ctl_elem_info_set_device.argtypes = [POINTER(snd_ctl_elem_info_t), c_uint]

# /usr/include/alsa/control.h:396
snd_ctl_elem_info_set_subdevice = _lib.snd_ctl_elem_info_set_subdevice
snd_ctl_elem_info_set_subdevice.restype = None
snd_ctl_elem_info_set_subdevice.argtypes = [POINTER(snd_ctl_elem_info_t), c_uint]

# /usr/include/alsa/control.h:397
snd_ctl_elem_info_set_name = _lib.snd_ctl_elem_info_set_name
snd_ctl_elem_info_set_name.restype = None
snd_ctl_elem_info_set_name.argtypes = [POINTER(snd_ctl_elem_info_t), c_char_p]

# /usr/include/alsa/control.h:398
snd_ctl_elem_info_set_index = _lib.snd_ctl_elem_info_set_index
snd_ctl_elem_info_set_index.restype = None
snd_ctl_elem_info_set_index.argtypes = [POINTER(snd_ctl_elem_info_t), c_uint]

# /usr/include/alsa/control.h:400
snd_ctl_elem_add_integer = _lib.snd_ctl_elem_add_integer
snd_ctl_elem_add_integer.restype = c_int
snd_ctl_elem_add_integer.argtypes = [POINTER(snd_ctl_t), POINTER(snd_ctl_elem_id_t), c_uint, c_long, c_long, c_long]

# /usr/include/alsa/control.h:401
snd_ctl_elem_add_integer64 = _lib.snd_ctl_elem_add_integer64
snd_ctl_elem_add_integer64.restype = c_int
snd_ctl_elem_add_integer64.argtypes = [POINTER(snd_ctl_t), POINTER(snd_ctl_elem_id_t), c_uint, c_longlong, c_longlong, c_longlong]

# /usr/include/alsa/control.h:402
snd_ctl_elem_add_boolean = _lib.snd_ctl_elem_add_boolean
snd_ctl_elem_add_boolean.restype = c_int
snd_ctl_elem_add_boolean.argtypes = [POINTER(snd_ctl_t), POINTER(snd_ctl_elem_id_t), c_uint]

# /usr/include/alsa/control.h:403
snd_ctl_elem_add_iec958 = _lib.snd_ctl_elem_add_iec958
snd_ctl_elem_add_iec958.restype = c_int
snd_ctl_elem_add_iec958.argtypes = [POINTER(snd_ctl_t), POINTER(snd_ctl_elem_id_t)]

# /usr/include/alsa/control.h:404
snd_ctl_elem_remove = _lib.snd_ctl_elem_remove
snd_ctl_elem_remove.restype = c_int
snd_ctl_elem_remove.argtypes = [POINTER(snd_ctl_t), POINTER(snd_ctl_elem_id_t)]

# /usr/include/alsa/control.h:406
snd_ctl_elem_value_sizeof = _lib.snd_ctl_elem_value_sizeof
snd_ctl_elem_value_sizeof.restype = c_size_t
snd_ctl_elem_value_sizeof.argtypes = []

# /usr/include/alsa/control.h:412
snd_ctl_elem_value_malloc = _lib.snd_ctl_elem_value_malloc
snd_ctl_elem_value_malloc.restype = c_int
snd_ctl_elem_value_malloc.argtypes = [POINTER(POINTER(snd_ctl_elem_value_t))]

# /usr/include/alsa/control.h:413
snd_ctl_elem_value_free = _lib.snd_ctl_elem_value_free
snd_ctl_elem_value_free.restype = None
snd_ctl_elem_value_free.argtypes = [POINTER(snd_ctl_elem_value_t)]

# /usr/include/alsa/control.h:414
snd_ctl_elem_value_clear = _lib.snd_ctl_elem_value_clear
snd_ctl_elem_value_clear.restype = None
snd_ctl_elem_value_clear.argtypes = [POINTER(snd_ctl_elem_value_t)]

# /usr/include/alsa/control.h:415
snd_ctl_elem_value_copy = _lib.snd_ctl_elem_value_copy
snd_ctl_elem_value_copy.restype = None
snd_ctl_elem_value_copy.argtypes = [POINTER(snd_ctl_elem_value_t), POINTER(snd_ctl_elem_value_t)]

# /usr/include/alsa/control.h:416
snd_ctl_elem_value_get_id = _lib.snd_ctl_elem_value_get_id
snd_ctl_elem_value_get_id.restype = None
snd_ctl_elem_value_get_id.argtypes = [POINTER(snd_ctl_elem_value_t), POINTER(snd_ctl_elem_id_t)]

# /usr/include/alsa/control.h:417
snd_ctl_elem_value_get_numid = _lib.snd_ctl_elem_value_get_numid
snd_ctl_elem_value_get_numid.restype = c_uint
snd_ctl_elem_value_get_numid.argtypes = [POINTER(snd_ctl_elem_value_t)]

# /usr/include/alsa/control.h:418
snd_ctl_elem_value_get_interface = _lib.snd_ctl_elem_value_get_interface
snd_ctl_elem_value_get_interface.restype = snd_ctl_elem_iface_t
snd_ctl_elem_value_get_interface.argtypes = [POINTER(snd_ctl_elem_value_t)]

# /usr/include/alsa/control.h:419
snd_ctl_elem_value_get_device = _lib.snd_ctl_elem_value_get_device
snd_ctl_elem_value_get_device.restype = c_uint
snd_ctl_elem_value_get_device.argtypes = [POINTER(snd_ctl_elem_value_t)]

# /usr/include/alsa/control.h:420
snd_ctl_elem_value_get_subdevice = _lib.snd_ctl_elem_value_get_subdevice
snd_ctl_elem_value_get_subdevice.restype = c_uint
snd_ctl_elem_value_get_subdevice.argtypes = [POINTER(snd_ctl_elem_value_t)]

# /usr/include/alsa/control.h:421
snd_ctl_elem_value_get_name = _lib.snd_ctl_elem_value_get_name
snd_ctl_elem_value_get_name.restype = c_char_p
snd_ctl_elem_value_get_name.argtypes = [POINTER(snd_ctl_elem_value_t)]

# /usr/include/alsa/control.h:422
snd_ctl_elem_value_get_index = _lib.snd_ctl_elem_value_get_index
snd_ctl_elem_value_get_index.restype = c_uint
snd_ctl_elem_value_get_index.argtypes = [POINTER(snd_ctl_elem_value_t)]

# /usr/include/alsa/control.h:423
snd_ctl_elem_value_set_id = _lib.snd_ctl_elem_value_set_id
snd_ctl_elem_value_set_id.restype = None
snd_ctl_elem_value_set_id.argtypes = [POINTER(snd_ctl_elem_value_t), POINTER(snd_ctl_elem_id_t)]

# /usr/include/alsa/control.h:424
snd_ctl_elem_value_set_numid = _lib.snd_ctl_elem_value_set_numid
snd_ctl_elem_value_set_numid.restype = None
snd_ctl_elem_value_set_numid.argtypes = [POINTER(snd_ctl_elem_value_t), c_uint]

# /usr/include/alsa/control.h:425
snd_ctl_elem_value_set_interface = _lib.snd_ctl_elem_value_set_interface
snd_ctl_elem_value_set_interface.restype = None
snd_ctl_elem_value_set_interface.argtypes = [POINTER(snd_ctl_elem_value_t), snd_ctl_elem_iface_t]

# /usr/include/alsa/control.h:426
snd_ctl_elem_value_set_device = _lib.snd_ctl_elem_value_set_device
snd_ctl_elem_value_set_device.restype = None
snd_ctl_elem_value_set_device.argtypes = [POINTER(snd_ctl_elem_value_t), c_uint]

# /usr/include/alsa/control.h:427
snd_ctl_elem_value_set_subdevice = _lib.snd_ctl_elem_value_set_subdevice
snd_ctl_elem_value_set_subdevice.restype = None
snd_ctl_elem_value_set_subdevice.argtypes = [POINTER(snd_ctl_elem_value_t), c_uint]

# /usr/include/alsa/control.h:428
snd_ctl_elem_value_set_name = _lib.snd_ctl_elem_value_set_name
snd_ctl_elem_value_set_name.restype = None
snd_ctl_elem_value_set_name.argtypes = [POINTER(snd_ctl_elem_value_t), c_char_p]

# /usr/include/alsa/control.h:429
snd_ctl_elem_value_set_index = _lib.snd_ctl_elem_value_set_index
snd_ctl_elem_value_set_index.restype = None
snd_ctl_elem_value_set_index.argtypes = [POINTER(snd_ctl_elem_value_t), c_uint]

# /usr/include/alsa/control.h:430
snd_ctl_elem_value_get_boolean = _lib.snd_ctl_elem_value_get_boolean
snd_ctl_elem_value_get_boolean.restype = c_int
snd_ctl_elem_value_get_boolean.argtypes = [POINTER(snd_ctl_elem_value_t), c_uint]

# /usr/include/alsa/control.h:431
snd_ctl_elem_value_get_integer = _lib.snd_ctl_elem_value_get_integer
snd_ctl_elem_value_get_integer.restype = c_long
snd_ctl_elem_value_get_integer.argtypes = [POINTER(snd_ctl_elem_value_t), c_uint]

# /usr/include/alsa/control.h:432
snd_ctl_elem_value_get_integer64 = _lib.snd_ctl_elem_value_get_integer64
snd_ctl_elem_value_get_integer64.restype = c_longlong
snd_ctl_elem_value_get_integer64.argtypes = [POINTER(snd_ctl_elem_value_t), c_uint]

# /usr/include/alsa/control.h:433
snd_ctl_elem_value_get_enumerated = _lib.snd_ctl_elem_value_get_enumerated
snd_ctl_elem_value_get_enumerated.restype = c_uint
snd_ctl_elem_value_get_enumerated.argtypes = [POINTER(snd_ctl_elem_value_t), c_uint]

# /usr/include/alsa/control.h:434
snd_ctl_elem_value_get_byte = _lib.snd_ctl_elem_value_get_byte
snd_ctl_elem_value_get_byte.restype = c_ubyte
snd_ctl_elem_value_get_byte.argtypes = [POINTER(snd_ctl_elem_value_t), c_uint]

# /usr/include/alsa/control.h:435
snd_ctl_elem_value_set_boolean = _lib.snd_ctl_elem_value_set_boolean
snd_ctl_elem_value_set_boolean.restype = None
snd_ctl_elem_value_set_boolean.argtypes = [POINTER(snd_ctl_elem_value_t), c_uint, c_long]

# /usr/include/alsa/control.h:436
snd_ctl_elem_value_set_integer = _lib.snd_ctl_elem_value_set_integer
snd_ctl_elem_value_set_integer.restype = None
snd_ctl_elem_value_set_integer.argtypes = [POINTER(snd_ctl_elem_value_t), c_uint, c_long]

# /usr/include/alsa/control.h:437
snd_ctl_elem_value_set_integer64 = _lib.snd_ctl_elem_value_set_integer64
snd_ctl_elem_value_set_integer64.restype = None
snd_ctl_elem_value_set_integer64.argtypes = [POINTER(snd_ctl_elem_value_t), c_uint, c_longlong]

# /usr/include/alsa/control.h:438
snd_ctl_elem_value_set_enumerated = _lib.snd_ctl_elem_value_set_enumerated
snd_ctl_elem_value_set_enumerated.restype = None
snd_ctl_elem_value_set_enumerated.argtypes = [POINTER(snd_ctl_elem_value_t), c_uint, c_uint]

# /usr/include/alsa/control.h:439
snd_ctl_elem_value_set_byte = _lib.snd_ctl_elem_value_set_byte
snd_ctl_elem_value_set_byte.restype = None
snd_ctl_elem_value_set_byte.argtypes = [POINTER(snd_ctl_elem_value_t), c_uint, c_ubyte]

# /usr/include/alsa/control.h:440
snd_ctl_elem_set_bytes = _lib.snd_ctl_elem_set_bytes
snd_ctl_elem_set_bytes.restype = None
snd_ctl_elem_set_bytes.argtypes = [POINTER(snd_ctl_elem_value_t), POINTER(None), c_size_t]

# /usr/include/alsa/control.h:441
snd_ctl_elem_value_get_bytes = _lib.snd_ctl_elem_value_get_bytes
snd_ctl_elem_value_get_bytes.restype = POINTER(c_void)
snd_ctl_elem_value_get_bytes.argtypes = [POINTER(snd_ctl_elem_value_t)]

# /usr/include/alsa/control.h:442
snd_ctl_elem_value_get_iec958 = _lib.snd_ctl_elem_value_get_iec958
snd_ctl_elem_value_get_iec958.restype = None
snd_ctl_elem_value_get_iec958.argtypes = [POINTER(snd_ctl_elem_value_t), POINTER(snd_aes_iec958_t)]

# /usr/include/alsa/control.h:443
snd_ctl_elem_value_set_iec958 = _lib.snd_ctl_elem_value_set_iec958
snd_ctl_elem_value_set_iec958.restype = None
snd_ctl_elem_value_set_iec958.argtypes = [POINTER(snd_ctl_elem_value_t), POINTER(snd_aes_iec958_t)]

class struct__snd_hctl_elem(Structure):
    __slots__ = [
    ]
struct__snd_hctl_elem._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_hctl_elem(Structure):
    __slots__ = [
    ]
struct__snd_hctl_elem._fields_ = [
    ('_opaque_struct', c_int)
]

snd_hctl_elem_t = struct__snd_hctl_elem 	# /usr/include/alsa/control.h:454
class struct__snd_hctl(Structure):
    __slots__ = [
    ]
struct__snd_hctl._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_hctl(Structure):
    __slots__ = [
    ]
struct__snd_hctl._fields_ = [
    ('_opaque_struct', c_int)
]

snd_hctl_t = struct__snd_hctl 	# /usr/include/alsa/control.h:457
snd_hctl_compare_t = CFUNCTYPE(c_int, POINTER(snd_hctl_elem_t), POINTER(snd_hctl_elem_t)) 	# /usr/include/alsa/control.h:465
# /usr/include/alsa/control.h:467
snd_hctl_compare_fast = _lib.snd_hctl_compare_fast
snd_hctl_compare_fast.restype = c_int
snd_hctl_compare_fast.argtypes = [POINTER(snd_hctl_elem_t), POINTER(snd_hctl_elem_t)]

snd_hctl_callback_t = CFUNCTYPE(c_int, POINTER(snd_hctl_t), c_uint, POINTER(snd_hctl_elem_t)) 	# /usr/include/alsa/control.h:476
snd_hctl_elem_callback_t = CFUNCTYPE(c_int, POINTER(snd_hctl_elem_t), c_uint) 	# /usr/include/alsa/control.h:485
# /usr/include/alsa/control.h:488
snd_hctl_open = _lib.snd_hctl_open
snd_hctl_open.restype = c_int
snd_hctl_open.argtypes = [POINTER(POINTER(snd_hctl_t)), c_char_p, c_int]

# /usr/include/alsa/control.h:489
snd_hctl_open_ctl = _lib.snd_hctl_open_ctl
snd_hctl_open_ctl.restype = c_int
snd_hctl_open_ctl.argtypes = [POINTER(POINTER(snd_hctl_t)), POINTER(snd_ctl_t)]

# /usr/include/alsa/control.h:490
snd_hctl_close = _lib.snd_hctl_close
snd_hctl_close.restype = c_int
snd_hctl_close.argtypes = [POINTER(snd_hctl_t)]

# /usr/include/alsa/control.h:491
snd_hctl_nonblock = _lib.snd_hctl_nonblock
snd_hctl_nonblock.restype = c_int
snd_hctl_nonblock.argtypes = [POINTER(snd_hctl_t), c_int]

# /usr/include/alsa/control.h:492
snd_hctl_poll_descriptors_count = _lib.snd_hctl_poll_descriptors_count
snd_hctl_poll_descriptors_count.restype = c_int
snd_hctl_poll_descriptors_count.argtypes = [POINTER(snd_hctl_t)]

class struct_pollfd(Structure):
    __slots__ = [
    ]
struct_pollfd._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/control.h:493
snd_hctl_poll_descriptors = _lib.snd_hctl_poll_descriptors
snd_hctl_poll_descriptors.restype = c_int
snd_hctl_poll_descriptors.argtypes = [POINTER(snd_hctl_t), POINTER(struct_pollfd), c_uint]

class struct_pollfd(Structure):
    __slots__ = [
    ]
struct_pollfd._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/control.h:494
snd_hctl_poll_descriptors_revents = _lib.snd_hctl_poll_descriptors_revents
snd_hctl_poll_descriptors_revents.restype = c_int
snd_hctl_poll_descriptors_revents.argtypes = [POINTER(snd_hctl_t), POINTER(struct_pollfd), c_uint, POINTER(c_ushort)]

# /usr/include/alsa/control.h:495
snd_hctl_get_count = _lib.snd_hctl_get_count
snd_hctl_get_count.restype = c_uint
snd_hctl_get_count.argtypes = [POINTER(snd_hctl_t)]

# /usr/include/alsa/control.h:496
snd_hctl_set_compare = _lib.snd_hctl_set_compare
snd_hctl_set_compare.restype = c_int
snd_hctl_set_compare.argtypes = [POINTER(snd_hctl_t), snd_hctl_compare_t]

# /usr/include/alsa/control.h:497
snd_hctl_first_elem = _lib.snd_hctl_first_elem
snd_hctl_first_elem.restype = POINTER(snd_hctl_elem_t)
snd_hctl_first_elem.argtypes = [POINTER(snd_hctl_t)]

# /usr/include/alsa/control.h:498
snd_hctl_last_elem = _lib.snd_hctl_last_elem
snd_hctl_last_elem.restype = POINTER(snd_hctl_elem_t)
snd_hctl_last_elem.argtypes = [POINTER(snd_hctl_t)]

# /usr/include/alsa/control.h:499
snd_hctl_find_elem = _lib.snd_hctl_find_elem
snd_hctl_find_elem.restype = POINTER(snd_hctl_elem_t)
snd_hctl_find_elem.argtypes = [POINTER(snd_hctl_t), POINTER(snd_ctl_elem_id_t)]

# /usr/include/alsa/control.h:500
snd_hctl_set_callback = _lib.snd_hctl_set_callback
snd_hctl_set_callback.restype = None
snd_hctl_set_callback.argtypes = [POINTER(snd_hctl_t), snd_hctl_callback_t]

# /usr/include/alsa/control.h:501
snd_hctl_set_callback_private = _lib.snd_hctl_set_callback_private
snd_hctl_set_callback_private.restype = None
snd_hctl_set_callback_private.argtypes = [POINTER(snd_hctl_t), POINTER(None)]

# /usr/include/alsa/control.h:502
snd_hctl_get_callback_private = _lib.snd_hctl_get_callback_private
snd_hctl_get_callback_private.restype = POINTER(c_void)
snd_hctl_get_callback_private.argtypes = [POINTER(snd_hctl_t)]

# /usr/include/alsa/control.h:503
snd_hctl_load = _lib.snd_hctl_load
snd_hctl_load.restype = c_int
snd_hctl_load.argtypes = [POINTER(snd_hctl_t)]

# /usr/include/alsa/control.h:504
snd_hctl_free = _lib.snd_hctl_free
snd_hctl_free.restype = c_int
snd_hctl_free.argtypes = [POINTER(snd_hctl_t)]

# /usr/include/alsa/control.h:505
snd_hctl_handle_events = _lib.snd_hctl_handle_events
snd_hctl_handle_events.restype = c_int
snd_hctl_handle_events.argtypes = [POINTER(snd_hctl_t)]

# /usr/include/alsa/control.h:506
snd_hctl_name = _lib.snd_hctl_name
snd_hctl_name.restype = c_char_p
snd_hctl_name.argtypes = [POINTER(snd_hctl_t)]

# /usr/include/alsa/control.h:507
snd_hctl_wait = _lib.snd_hctl_wait
snd_hctl_wait.restype = c_int
snd_hctl_wait.argtypes = [POINTER(snd_hctl_t), c_int]

# /usr/include/alsa/control.h:508
snd_hctl_ctl = _lib.snd_hctl_ctl
snd_hctl_ctl.restype = POINTER(snd_ctl_t)
snd_hctl_ctl.argtypes = [POINTER(snd_hctl_t)]

# /usr/include/alsa/control.h:510
snd_hctl_elem_next = _lib.snd_hctl_elem_next
snd_hctl_elem_next.restype = POINTER(snd_hctl_elem_t)
snd_hctl_elem_next.argtypes = [POINTER(snd_hctl_elem_t)]

# /usr/include/alsa/control.h:511
snd_hctl_elem_prev = _lib.snd_hctl_elem_prev
snd_hctl_elem_prev.restype = POINTER(snd_hctl_elem_t)
snd_hctl_elem_prev.argtypes = [POINTER(snd_hctl_elem_t)]

# /usr/include/alsa/control.h:512
snd_hctl_elem_info = _lib.snd_hctl_elem_info
snd_hctl_elem_info.restype = c_int
snd_hctl_elem_info.argtypes = [POINTER(snd_hctl_elem_t), POINTER(snd_ctl_elem_info_t)]

# /usr/include/alsa/control.h:513
snd_hctl_elem_read = _lib.snd_hctl_elem_read
snd_hctl_elem_read.restype = c_int
snd_hctl_elem_read.argtypes = [POINTER(snd_hctl_elem_t), POINTER(snd_ctl_elem_value_t)]

# /usr/include/alsa/control.h:514
snd_hctl_elem_write = _lib.snd_hctl_elem_write
snd_hctl_elem_write.restype = c_int
snd_hctl_elem_write.argtypes = [POINTER(snd_hctl_elem_t), POINTER(snd_ctl_elem_value_t)]

# /usr/include/alsa/control.h:515
snd_hctl_elem_tlv_read = _lib.snd_hctl_elem_tlv_read
snd_hctl_elem_tlv_read.restype = c_int
snd_hctl_elem_tlv_read.argtypes = [POINTER(snd_hctl_elem_t), POINTER(c_uint), c_uint]

# /usr/include/alsa/control.h:516
snd_hctl_elem_tlv_write = _lib.snd_hctl_elem_tlv_write
snd_hctl_elem_tlv_write.restype = c_int
snd_hctl_elem_tlv_write.argtypes = [POINTER(snd_hctl_elem_t), POINTER(c_uint)]

# /usr/include/alsa/control.h:517
snd_hctl_elem_tlv_command = _lib.snd_hctl_elem_tlv_command
snd_hctl_elem_tlv_command.restype = c_int
snd_hctl_elem_tlv_command.argtypes = [POINTER(snd_hctl_elem_t), POINTER(c_uint)]

# /usr/include/alsa/control.h:519
snd_hctl_elem_get_hctl = _lib.snd_hctl_elem_get_hctl
snd_hctl_elem_get_hctl.restype = POINTER(snd_hctl_t)
snd_hctl_elem_get_hctl.argtypes = [POINTER(snd_hctl_elem_t)]

# /usr/include/alsa/control.h:521
snd_hctl_elem_get_id = _lib.snd_hctl_elem_get_id
snd_hctl_elem_get_id.restype = None
snd_hctl_elem_get_id.argtypes = [POINTER(snd_hctl_elem_t), POINTER(snd_ctl_elem_id_t)]

# /usr/include/alsa/control.h:522
snd_hctl_elem_get_numid = _lib.snd_hctl_elem_get_numid
snd_hctl_elem_get_numid.restype = c_uint
snd_hctl_elem_get_numid.argtypes = [POINTER(snd_hctl_elem_t)]

# /usr/include/alsa/control.h:523
snd_hctl_elem_get_interface = _lib.snd_hctl_elem_get_interface
snd_hctl_elem_get_interface.restype = snd_ctl_elem_iface_t
snd_hctl_elem_get_interface.argtypes = [POINTER(snd_hctl_elem_t)]

# /usr/include/alsa/control.h:524
snd_hctl_elem_get_device = _lib.snd_hctl_elem_get_device
snd_hctl_elem_get_device.restype = c_uint
snd_hctl_elem_get_device.argtypes = [POINTER(snd_hctl_elem_t)]

# /usr/include/alsa/control.h:525
snd_hctl_elem_get_subdevice = _lib.snd_hctl_elem_get_subdevice
snd_hctl_elem_get_subdevice.restype = c_uint
snd_hctl_elem_get_subdevice.argtypes = [POINTER(snd_hctl_elem_t)]

# /usr/include/alsa/control.h:526
snd_hctl_elem_get_name = _lib.snd_hctl_elem_get_name
snd_hctl_elem_get_name.restype = c_char_p
snd_hctl_elem_get_name.argtypes = [POINTER(snd_hctl_elem_t)]

# /usr/include/alsa/control.h:527
snd_hctl_elem_get_index = _lib.snd_hctl_elem_get_index
snd_hctl_elem_get_index.restype = c_uint
snd_hctl_elem_get_index.argtypes = [POINTER(snd_hctl_elem_t)]

# /usr/include/alsa/control.h:528
snd_hctl_elem_set_callback = _lib.snd_hctl_elem_set_callback
snd_hctl_elem_set_callback.restype = None
snd_hctl_elem_set_callback.argtypes = [POINTER(snd_hctl_elem_t), snd_hctl_elem_callback_t]

# /usr/include/alsa/control.h:529
snd_hctl_elem_get_callback_private = _lib.snd_hctl_elem_get_callback_private
snd_hctl_elem_get_callback_private.restype = POINTER(c_void)
snd_hctl_elem_get_callback_private.argtypes = [POINTER(snd_hctl_elem_t)]

# /usr/include/alsa/control.h:530
snd_hctl_elem_set_callback_private = _lib.snd_hctl_elem_set_callback_private
snd_hctl_elem_set_callback_private.restype = None
snd_hctl_elem_set_callback_private.argtypes = [POINTER(snd_hctl_elem_t), POINTER(None)]

# /usr/include/alsa/control.h:543
snd_sctl_build = _lib.snd_sctl_build
snd_sctl_build.restype = c_int
snd_sctl_build.argtypes = [POINTER(POINTER(snd_sctl_t)), POINTER(snd_ctl_t), POINTER(snd_config_t), POINTER(snd_config_t), c_int]

# /usr/include/alsa/control.h:545
snd_sctl_free = _lib.snd_sctl_free
snd_sctl_free.restype = c_int
snd_sctl_free.argtypes = [POINTER(snd_sctl_t)]

# /usr/include/alsa/control.h:546
snd_sctl_install = _lib.snd_sctl_install
snd_sctl_install.restype = c_int
snd_sctl_install.argtypes = [POINTER(snd_sctl_t)]

# /usr/include/alsa/control.h:547
snd_sctl_remove = _lib.snd_sctl_remove
snd_sctl_remove.restype = c_int
snd_sctl_remove.argtypes = [POINTER(snd_sctl_t)]

class struct__snd_mixer(Structure):
    __slots__ = [
    ]
struct__snd_mixer._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_mixer(Structure):
    __slots__ = [
    ]
struct__snd_mixer._fields_ = [
    ('_opaque_struct', c_int)
]

snd_mixer_t = struct__snd_mixer 	# /usr/include/alsa/mixer.h:42
class struct__snd_mixer_class(Structure):
    __slots__ = [
    ]
struct__snd_mixer_class._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_mixer_class(Structure):
    __slots__ = [
    ]
struct__snd_mixer_class._fields_ = [
    ('_opaque_struct', c_int)
]

snd_mixer_class_t = struct__snd_mixer_class 	# /usr/include/alsa/mixer.h:44
class struct__snd_mixer_elem(Structure):
    __slots__ = [
    ]
struct__snd_mixer_elem._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_mixer_elem(Structure):
    __slots__ = [
    ]
struct__snd_mixer_elem._fields_ = [
    ('_opaque_struct', c_int)
]

snd_mixer_elem_t = struct__snd_mixer_elem 	# /usr/include/alsa/mixer.h:46
snd_mixer_callback_t = CFUNCTYPE(c_int, POINTER(snd_mixer_t), c_uint, POINTER(snd_mixer_elem_t)) 	# /usr/include/alsa/mixer.h:55
snd_mixer_elem_callback_t = CFUNCTYPE(c_int, POINTER(snd_mixer_elem_t), c_uint) 	# /usr/include/alsa/mixer.h:65
snd_mixer_compare_t = CFUNCTYPE(c_int, POINTER(snd_mixer_elem_t), POINTER(snd_mixer_elem_t)) 	# /usr/include/alsa/mixer.h:74
snd_mixer_event_t = CFUNCTYPE(c_int, POINTER(snd_mixer_class_t), c_uint, POINTER(snd_hctl_elem_t), POINTER(snd_mixer_elem_t)) 	# /usr/include/alsa/mixer.h:85
enum__snd_mixer_elem_type = c_int
SND_MIXER_ELEM_SIMPLE = 1
SND_MIXER_ELEM_LAST = 0
snd_mixer_elem_type_t = enum__snd_mixer_elem_type 	# /usr/include/alsa/mixer.h:94
# /usr/include/alsa/mixer.h:96
snd_mixer_open = _lib.snd_mixer_open
snd_mixer_open.restype = c_int
snd_mixer_open.argtypes = [POINTER(POINTER(snd_mixer_t)), c_int]

# /usr/include/alsa/mixer.h:97
snd_mixer_close = _lib.snd_mixer_close
snd_mixer_close.restype = c_int
snd_mixer_close.argtypes = [POINTER(snd_mixer_t)]

# /usr/include/alsa/mixer.h:98
snd_mixer_first_elem = _lib.snd_mixer_first_elem
snd_mixer_first_elem.restype = POINTER(snd_mixer_elem_t)
snd_mixer_first_elem.argtypes = [POINTER(snd_mixer_t)]

# /usr/include/alsa/mixer.h:99
snd_mixer_last_elem = _lib.snd_mixer_last_elem
snd_mixer_last_elem.restype = POINTER(snd_mixer_elem_t)
snd_mixer_last_elem.argtypes = [POINTER(snd_mixer_t)]

# /usr/include/alsa/mixer.h:100
snd_mixer_handle_events = _lib.snd_mixer_handle_events
snd_mixer_handle_events.restype = c_int
snd_mixer_handle_events.argtypes = [POINTER(snd_mixer_t)]

# /usr/include/alsa/mixer.h:101
snd_mixer_attach = _lib.snd_mixer_attach
snd_mixer_attach.restype = c_int
snd_mixer_attach.argtypes = [POINTER(snd_mixer_t), c_char_p]

# /usr/include/alsa/mixer.h:102
snd_mixer_attach_hctl = _lib.snd_mixer_attach_hctl
snd_mixer_attach_hctl.restype = c_int
snd_mixer_attach_hctl.argtypes = [POINTER(snd_mixer_t), POINTER(snd_hctl_t)]

# /usr/include/alsa/mixer.h:103
snd_mixer_detach = _lib.snd_mixer_detach
snd_mixer_detach.restype = c_int
snd_mixer_detach.argtypes = [POINTER(snd_mixer_t), c_char_p]

# /usr/include/alsa/mixer.h:104
snd_mixer_detach_hctl = _lib.snd_mixer_detach_hctl
snd_mixer_detach_hctl.restype = c_int
snd_mixer_detach_hctl.argtypes = [POINTER(snd_mixer_t), POINTER(snd_hctl_t)]

# /usr/include/alsa/mixer.h:105
snd_mixer_get_hctl = _lib.snd_mixer_get_hctl
snd_mixer_get_hctl.restype = c_int
snd_mixer_get_hctl.argtypes = [POINTER(snd_mixer_t), c_char_p, POINTER(POINTER(snd_hctl_t))]

# /usr/include/alsa/mixer.h:106
snd_mixer_poll_descriptors_count = _lib.snd_mixer_poll_descriptors_count
snd_mixer_poll_descriptors_count.restype = c_int
snd_mixer_poll_descriptors_count.argtypes = [POINTER(snd_mixer_t)]

class struct_pollfd(Structure):
    __slots__ = [
    ]
struct_pollfd._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/mixer.h:107
snd_mixer_poll_descriptors = _lib.snd_mixer_poll_descriptors
snd_mixer_poll_descriptors.restype = c_int
snd_mixer_poll_descriptors.argtypes = [POINTER(snd_mixer_t), POINTER(struct_pollfd), c_uint]

class struct_pollfd(Structure):
    __slots__ = [
    ]
struct_pollfd._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/mixer.h:108
snd_mixer_poll_descriptors_revents = _lib.snd_mixer_poll_descriptors_revents
snd_mixer_poll_descriptors_revents.restype = c_int
snd_mixer_poll_descriptors_revents.argtypes = [POINTER(snd_mixer_t), POINTER(struct_pollfd), c_uint, POINTER(c_ushort)]

# /usr/include/alsa/mixer.h:109
snd_mixer_load = _lib.snd_mixer_load
snd_mixer_load.restype = c_int
snd_mixer_load.argtypes = [POINTER(snd_mixer_t)]

# /usr/include/alsa/mixer.h:110
snd_mixer_free = _lib.snd_mixer_free
snd_mixer_free.restype = None
snd_mixer_free.argtypes = [POINTER(snd_mixer_t)]

# /usr/include/alsa/mixer.h:111
snd_mixer_wait = _lib.snd_mixer_wait
snd_mixer_wait.restype = c_int
snd_mixer_wait.argtypes = [POINTER(snd_mixer_t), c_int]

# /usr/include/alsa/mixer.h:112
snd_mixer_set_compare = _lib.snd_mixer_set_compare
snd_mixer_set_compare.restype = c_int
snd_mixer_set_compare.argtypes = [POINTER(snd_mixer_t), snd_mixer_compare_t]

# /usr/include/alsa/mixer.h:113
snd_mixer_set_callback = _lib.snd_mixer_set_callback
snd_mixer_set_callback.restype = None
snd_mixer_set_callback.argtypes = [POINTER(snd_mixer_t), snd_mixer_callback_t]

# /usr/include/alsa/mixer.h:114
snd_mixer_get_callback_private = _lib.snd_mixer_get_callback_private
snd_mixer_get_callback_private.restype = POINTER(c_void)
snd_mixer_get_callback_private.argtypes = [POINTER(snd_mixer_t)]

# /usr/include/alsa/mixer.h:115
snd_mixer_set_callback_private = _lib.snd_mixer_set_callback_private
snd_mixer_set_callback_private.restype = None
snd_mixer_set_callback_private.argtypes = [POINTER(snd_mixer_t), POINTER(None)]

# /usr/include/alsa/mixer.h:116
snd_mixer_get_count = _lib.snd_mixer_get_count
snd_mixer_get_count.restype = c_uint
snd_mixer_get_count.argtypes = [POINTER(snd_mixer_t)]

# /usr/include/alsa/mixer.h:117
snd_mixer_class_unregister = _lib.snd_mixer_class_unregister
snd_mixer_class_unregister.restype = c_int
snd_mixer_class_unregister.argtypes = [POINTER(snd_mixer_class_t)]

# /usr/include/alsa/mixer.h:119
snd_mixer_elem_next = _lib.snd_mixer_elem_next
snd_mixer_elem_next.restype = POINTER(snd_mixer_elem_t)
snd_mixer_elem_next.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:120
snd_mixer_elem_prev = _lib.snd_mixer_elem_prev
snd_mixer_elem_prev.restype = POINTER(snd_mixer_elem_t)
snd_mixer_elem_prev.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:121
snd_mixer_elem_set_callback = _lib.snd_mixer_elem_set_callback
snd_mixer_elem_set_callback.restype = None
snd_mixer_elem_set_callback.argtypes = [POINTER(snd_mixer_elem_t), snd_mixer_elem_callback_t]

# /usr/include/alsa/mixer.h:122
snd_mixer_elem_get_callback_private = _lib.snd_mixer_elem_get_callback_private
snd_mixer_elem_get_callback_private.restype = POINTER(c_void)
snd_mixer_elem_get_callback_private.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:123
snd_mixer_elem_set_callback_private = _lib.snd_mixer_elem_set_callback_private
snd_mixer_elem_set_callback_private.restype = None
snd_mixer_elem_set_callback_private.argtypes = [POINTER(snd_mixer_elem_t), POINTER(None)]

# /usr/include/alsa/mixer.h:124
snd_mixer_elem_get_type = _lib.snd_mixer_elem_get_type
snd_mixer_elem_get_type.restype = snd_mixer_elem_type_t
snd_mixer_elem_get_type.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:126
snd_mixer_class_register = _lib.snd_mixer_class_register
snd_mixer_class_register.restype = c_int
snd_mixer_class_register.argtypes = [POINTER(snd_mixer_class_t), POINTER(snd_mixer_t)]

'''
# XXX these two functions don't exist in my libasound.so
# /usr/include/alsa/mixer.h:127
snd_mixer_add_elem = _lib.snd_mixer_add_elem
snd_mixer_add_elem.restype = c_int
snd_mixer_add_elem.argtypes = [POINTER(snd_mixer_t), POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:128
snd_mixer_remove_elem = _lib.snd_mixer_remove_elem
snd_mixer_remove_elem.restype = c_int
snd_mixer_remove_elem.argtypes = [POINTER(snd_mixer_t), POINTER(snd_mixer_elem_t)]
'''

# /usr/include/alsa/mixer.h:129
snd_mixer_elem_new = _lib.snd_mixer_elem_new
snd_mixer_elem_new.restype = c_int
snd_mixer_elem_new.argtypes = [POINTER(POINTER(snd_mixer_elem_t)), snd_mixer_elem_type_t, c_int, POINTER(None), CFUNCTYPE(None, POINTER(snd_mixer_elem_t))]

# /usr/include/alsa/mixer.h:134
snd_mixer_elem_add = _lib.snd_mixer_elem_add
snd_mixer_elem_add.restype = c_int
snd_mixer_elem_add.argtypes = [POINTER(snd_mixer_elem_t), POINTER(snd_mixer_class_t)]

# /usr/include/alsa/mixer.h:135
snd_mixer_elem_remove = _lib.snd_mixer_elem_remove
snd_mixer_elem_remove.restype = c_int
snd_mixer_elem_remove.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:136
snd_mixer_elem_free = _lib.snd_mixer_elem_free
snd_mixer_elem_free.restype = None
snd_mixer_elem_free.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:137
snd_mixer_elem_info = _lib.snd_mixer_elem_info
snd_mixer_elem_info.restype = c_int
snd_mixer_elem_info.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:138
snd_mixer_elem_value = _lib.snd_mixer_elem_value
snd_mixer_elem_value.restype = c_int
snd_mixer_elem_value.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:139
snd_mixer_elem_attach = _lib.snd_mixer_elem_attach
snd_mixer_elem_attach.restype = c_int
snd_mixer_elem_attach.argtypes = [POINTER(snd_mixer_elem_t), POINTER(snd_hctl_elem_t)]

# /usr/include/alsa/mixer.h:140
snd_mixer_elem_detach = _lib.snd_mixer_elem_detach
snd_mixer_elem_detach.restype = c_int
snd_mixer_elem_detach.argtypes = [POINTER(snd_mixer_elem_t), POINTER(snd_hctl_elem_t)]

# /usr/include/alsa/mixer.h:141
snd_mixer_elem_empty = _lib.snd_mixer_elem_empty
snd_mixer_elem_empty.restype = c_int
snd_mixer_elem_empty.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:142
snd_mixer_elem_get_private = _lib.snd_mixer_elem_get_private
snd_mixer_elem_get_private.restype = POINTER(c_void)
snd_mixer_elem_get_private.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:144
snd_mixer_class_sizeof = _lib.snd_mixer_class_sizeof
snd_mixer_class_sizeof.restype = c_size_t
snd_mixer_class_sizeof.argtypes = []

# /usr/include/alsa/mixer.h:150
snd_mixer_class_malloc = _lib.snd_mixer_class_malloc
snd_mixer_class_malloc.restype = c_int
snd_mixer_class_malloc.argtypes = [POINTER(POINTER(snd_mixer_class_t))]

# /usr/include/alsa/mixer.h:151
snd_mixer_class_free = _lib.snd_mixer_class_free
snd_mixer_class_free.restype = None
snd_mixer_class_free.argtypes = [POINTER(snd_mixer_class_t)]

# /usr/include/alsa/mixer.h:152
snd_mixer_class_copy = _lib.snd_mixer_class_copy
snd_mixer_class_copy.restype = None
snd_mixer_class_copy.argtypes = [POINTER(snd_mixer_class_t), POINTER(snd_mixer_class_t)]

# /usr/include/alsa/mixer.h:153
snd_mixer_class_get_mixer = _lib.snd_mixer_class_get_mixer
snd_mixer_class_get_mixer.restype = POINTER(snd_mixer_t)
snd_mixer_class_get_mixer.argtypes = [POINTER(snd_mixer_class_t)]

# /usr/include/alsa/mixer.h:154
snd_mixer_class_get_event = _lib.snd_mixer_class_get_event
snd_mixer_class_get_event.restype = snd_mixer_event_t
snd_mixer_class_get_event.argtypes = [POINTER(snd_mixer_class_t)]

# /usr/include/alsa/mixer.h:155
snd_mixer_class_get_private = _lib.snd_mixer_class_get_private
snd_mixer_class_get_private.restype = POINTER(c_void)
snd_mixer_class_get_private.argtypes = [POINTER(snd_mixer_class_t)]

# /usr/include/alsa/mixer.h:156
snd_mixer_class_get_compare = _lib.snd_mixer_class_get_compare
snd_mixer_class_get_compare.restype = snd_mixer_compare_t
snd_mixer_class_get_compare.argtypes = [POINTER(snd_mixer_class_t)]

# /usr/include/alsa/mixer.h:157
snd_mixer_class_set_event = _lib.snd_mixer_class_set_event
snd_mixer_class_set_event.restype = c_int
snd_mixer_class_set_event.argtypes = [POINTER(snd_mixer_class_t), snd_mixer_event_t]

# /usr/include/alsa/mixer.h:158
snd_mixer_class_set_private = _lib.snd_mixer_class_set_private
snd_mixer_class_set_private.restype = c_int
snd_mixer_class_set_private.argtypes = [POINTER(snd_mixer_class_t), POINTER(None)]

# /usr/include/alsa/mixer.h:159
snd_mixer_class_set_private_free = _lib.snd_mixer_class_set_private_free
snd_mixer_class_set_private_free.restype = c_int
snd_mixer_class_set_private_free.argtypes = [POINTER(snd_mixer_class_t), CFUNCTYPE(None, POINTER(snd_mixer_class_t))]

# /usr/include/alsa/mixer.h:160
snd_mixer_class_set_compare = _lib.snd_mixer_class_set_compare
snd_mixer_class_set_compare.restype = c_int
snd_mixer_class_set_compare.argtypes = [POINTER(snd_mixer_class_t), snd_mixer_compare_t]

enum__snd_mixer_selem_channel_id = c_int
SND_MIXER_SCHN_UNKNOWN = 1
SND_MIXER_SCHN_FRONT_LEFT = 0
SND_MIXER_SCHN_FRONT_RIGHT = 1
SND_MIXER_SCHN_REAR_LEFT = 2
SND_MIXER_SCHN_REAR_RIGHT = 3
SND_MIXER_SCHN_FRONT_CENTER = 4
SND_MIXER_SCHN_WOOFER = 5
SND_MIXER_SCHN_SIDE_LEFT = 6
SND_MIXER_SCHN_SIDE_RIGHT = 7
SND_MIXER_SCHN_REAR_CENTER = 8
SND_MIXER_SCHN_LAST = 31
SND_MIXER_SCHN_MONO = 0
snd_mixer_selem_channel_id_t = enum__snd_mixer_selem_channel_id 	# /usr/include/alsa/mixer.h:196
class struct__snd_mixer_selem_id(Structure):
    __slots__ = [
    ]
struct__snd_mixer_selem_id._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_mixer_selem_id(Structure):
    __slots__ = [
    ]
struct__snd_mixer_selem_id._fields_ = [
    ('_opaque_struct', c_int)
]

snd_mixer_selem_id_t = struct__snd_mixer_selem_id 	# /usr/include/alsa/mixer.h:221
# /usr/include/alsa/mixer.h:223
snd_mixer_selem_channel_name = _lib.snd_mixer_selem_channel_name
snd_mixer_selem_channel_name.restype = c_char_p
snd_mixer_selem_channel_name.argtypes = [snd_mixer_selem_channel_id_t]

class struct_snd_mixer_selem_regopt(Structure):
    __slots__ = [
    ]
struct_snd_mixer_selem_regopt._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/mixer.h:225
snd_mixer_selem_register = _lib.snd_mixer_selem_register
snd_mixer_selem_register.restype = c_int
snd_mixer_selem_register.argtypes = [POINTER(snd_mixer_t), POINTER(struct_snd_mixer_selem_regopt), POINTER(POINTER(snd_mixer_class_t))]

# /usr/include/alsa/mixer.h:228
snd_mixer_selem_get_id = _lib.snd_mixer_selem_get_id
snd_mixer_selem_get_id.restype = None
snd_mixer_selem_get_id.argtypes = [POINTER(snd_mixer_elem_t), POINTER(snd_mixer_selem_id_t)]

# /usr/include/alsa/mixer.h:230
snd_mixer_selem_get_name = _lib.snd_mixer_selem_get_name
snd_mixer_selem_get_name.restype = c_char_p
snd_mixer_selem_get_name.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:231
snd_mixer_selem_get_index = _lib.snd_mixer_selem_get_index
snd_mixer_selem_get_index.restype = c_uint
snd_mixer_selem_get_index.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:232
snd_mixer_find_selem = _lib.snd_mixer_find_selem
snd_mixer_find_selem.restype = POINTER(snd_mixer_elem_t)
snd_mixer_find_selem.argtypes = [POINTER(snd_mixer_t), POINTER(snd_mixer_selem_id_t)]

# /usr/include/alsa/mixer.h:235
snd_mixer_selem_is_active = _lib.snd_mixer_selem_is_active
snd_mixer_selem_is_active.restype = c_int
snd_mixer_selem_is_active.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:236
snd_mixer_selem_is_playback_mono = _lib.snd_mixer_selem_is_playback_mono
snd_mixer_selem_is_playback_mono.restype = c_int
snd_mixer_selem_is_playback_mono.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:237
snd_mixer_selem_has_playback_channel = _lib.snd_mixer_selem_has_playback_channel
snd_mixer_selem_has_playback_channel.restype = c_int
snd_mixer_selem_has_playback_channel.argtypes = [POINTER(snd_mixer_elem_t), snd_mixer_selem_channel_id_t]

# /usr/include/alsa/mixer.h:238
snd_mixer_selem_is_capture_mono = _lib.snd_mixer_selem_is_capture_mono
snd_mixer_selem_is_capture_mono.restype = c_int
snd_mixer_selem_is_capture_mono.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:239
snd_mixer_selem_has_capture_channel = _lib.snd_mixer_selem_has_capture_channel
snd_mixer_selem_has_capture_channel.restype = c_int
snd_mixer_selem_has_capture_channel.argtypes = [POINTER(snd_mixer_elem_t), snd_mixer_selem_channel_id_t]

# /usr/include/alsa/mixer.h:240
snd_mixer_selem_get_capture_group = _lib.snd_mixer_selem_get_capture_group
snd_mixer_selem_get_capture_group.restype = c_int
snd_mixer_selem_get_capture_group.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:241
snd_mixer_selem_has_common_volume = _lib.snd_mixer_selem_has_common_volume
snd_mixer_selem_has_common_volume.restype = c_int
snd_mixer_selem_has_common_volume.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:242
snd_mixer_selem_has_playback_volume = _lib.snd_mixer_selem_has_playback_volume
snd_mixer_selem_has_playback_volume.restype = c_int
snd_mixer_selem_has_playback_volume.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:243
snd_mixer_selem_has_playback_volume_joined = _lib.snd_mixer_selem_has_playback_volume_joined
snd_mixer_selem_has_playback_volume_joined.restype = c_int
snd_mixer_selem_has_playback_volume_joined.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:244
snd_mixer_selem_has_capture_volume = _lib.snd_mixer_selem_has_capture_volume
snd_mixer_selem_has_capture_volume.restype = c_int
snd_mixer_selem_has_capture_volume.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:245
snd_mixer_selem_has_capture_volume_joined = _lib.snd_mixer_selem_has_capture_volume_joined
snd_mixer_selem_has_capture_volume_joined.restype = c_int
snd_mixer_selem_has_capture_volume_joined.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:246
snd_mixer_selem_has_common_switch = _lib.snd_mixer_selem_has_common_switch
snd_mixer_selem_has_common_switch.restype = c_int
snd_mixer_selem_has_common_switch.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:247
snd_mixer_selem_has_playback_switch = _lib.snd_mixer_selem_has_playback_switch
snd_mixer_selem_has_playback_switch.restype = c_int
snd_mixer_selem_has_playback_switch.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:248
snd_mixer_selem_has_playback_switch_joined = _lib.snd_mixer_selem_has_playback_switch_joined
snd_mixer_selem_has_playback_switch_joined.restype = c_int
snd_mixer_selem_has_playback_switch_joined.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:249
snd_mixer_selem_has_capture_switch = _lib.snd_mixer_selem_has_capture_switch
snd_mixer_selem_has_capture_switch.restype = c_int
snd_mixer_selem_has_capture_switch.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:250
snd_mixer_selem_has_capture_switch_joined = _lib.snd_mixer_selem_has_capture_switch_joined
snd_mixer_selem_has_capture_switch_joined.restype = c_int
snd_mixer_selem_has_capture_switch_joined.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:251
snd_mixer_selem_has_capture_switch_exclusive = _lib.snd_mixer_selem_has_capture_switch_exclusive
snd_mixer_selem_has_capture_switch_exclusive.restype = c_int
snd_mixer_selem_has_capture_switch_exclusive.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:253
snd_mixer_selem_get_playback_volume = _lib.snd_mixer_selem_get_playback_volume
snd_mixer_selem_get_playback_volume.restype = c_int
snd_mixer_selem_get_playback_volume.argtypes = [POINTER(snd_mixer_elem_t), snd_mixer_selem_channel_id_t, POINTER(c_long)]

# /usr/include/alsa/mixer.h:254
snd_mixer_selem_get_capture_volume = _lib.snd_mixer_selem_get_capture_volume
snd_mixer_selem_get_capture_volume.restype = c_int
snd_mixer_selem_get_capture_volume.argtypes = [POINTER(snd_mixer_elem_t), snd_mixer_selem_channel_id_t, POINTER(c_long)]

# /usr/include/alsa/mixer.h:255
snd_mixer_selem_get_playback_dB = _lib.snd_mixer_selem_get_playback_dB
snd_mixer_selem_get_playback_dB.restype = c_int
snd_mixer_selem_get_playback_dB.argtypes = [POINTER(snd_mixer_elem_t), snd_mixer_selem_channel_id_t, POINTER(c_long)]

# /usr/include/alsa/mixer.h:256
snd_mixer_selem_get_capture_dB = _lib.snd_mixer_selem_get_capture_dB
snd_mixer_selem_get_capture_dB.restype = c_int
snd_mixer_selem_get_capture_dB.argtypes = [POINTER(snd_mixer_elem_t), snd_mixer_selem_channel_id_t, POINTER(c_long)]

# /usr/include/alsa/mixer.h:257
snd_mixer_selem_get_playback_switch = _lib.snd_mixer_selem_get_playback_switch
snd_mixer_selem_get_playback_switch.restype = c_int
snd_mixer_selem_get_playback_switch.argtypes = [POINTER(snd_mixer_elem_t), snd_mixer_selem_channel_id_t, POINTER(c_int)]

# /usr/include/alsa/mixer.h:258
snd_mixer_selem_get_capture_switch = _lib.snd_mixer_selem_get_capture_switch
snd_mixer_selem_get_capture_switch.restype = c_int
snd_mixer_selem_get_capture_switch.argtypes = [POINTER(snd_mixer_elem_t), snd_mixer_selem_channel_id_t, POINTER(c_int)]

# /usr/include/alsa/mixer.h:259
snd_mixer_selem_set_playback_volume = _lib.snd_mixer_selem_set_playback_volume
snd_mixer_selem_set_playback_volume.restype = c_int
snd_mixer_selem_set_playback_volume.argtypes = [POINTER(snd_mixer_elem_t), snd_mixer_selem_channel_id_t, c_long]

# /usr/include/alsa/mixer.h:260
snd_mixer_selem_set_capture_volume = _lib.snd_mixer_selem_set_capture_volume
snd_mixer_selem_set_capture_volume.restype = c_int
snd_mixer_selem_set_capture_volume.argtypes = [POINTER(snd_mixer_elem_t), snd_mixer_selem_channel_id_t, c_long]

# /usr/include/alsa/mixer.h:261
snd_mixer_selem_set_playback_dB = _lib.snd_mixer_selem_set_playback_dB
snd_mixer_selem_set_playback_dB.restype = c_int
snd_mixer_selem_set_playback_dB.argtypes = [POINTER(snd_mixer_elem_t), snd_mixer_selem_channel_id_t, c_long, c_int]

# /usr/include/alsa/mixer.h:262
snd_mixer_selem_set_capture_dB = _lib.snd_mixer_selem_set_capture_dB
snd_mixer_selem_set_capture_dB.restype = c_int
snd_mixer_selem_set_capture_dB.argtypes = [POINTER(snd_mixer_elem_t), snd_mixer_selem_channel_id_t, c_long, c_int]

# /usr/include/alsa/mixer.h:263
snd_mixer_selem_set_playback_volume_all = _lib.snd_mixer_selem_set_playback_volume_all
snd_mixer_selem_set_playback_volume_all.restype = c_int
snd_mixer_selem_set_playback_volume_all.argtypes = [POINTER(snd_mixer_elem_t), c_long]

# /usr/include/alsa/mixer.h:264
snd_mixer_selem_set_capture_volume_all = _lib.snd_mixer_selem_set_capture_volume_all
snd_mixer_selem_set_capture_volume_all.restype = c_int
snd_mixer_selem_set_capture_volume_all.argtypes = [POINTER(snd_mixer_elem_t), c_long]

# /usr/include/alsa/mixer.h:265
snd_mixer_selem_set_playback_dB_all = _lib.snd_mixer_selem_set_playback_dB_all
snd_mixer_selem_set_playback_dB_all.restype = c_int
snd_mixer_selem_set_playback_dB_all.argtypes = [POINTER(snd_mixer_elem_t), c_long, c_int]

# /usr/include/alsa/mixer.h:266
snd_mixer_selem_set_capture_dB_all = _lib.snd_mixer_selem_set_capture_dB_all
snd_mixer_selem_set_capture_dB_all.restype = c_int
snd_mixer_selem_set_capture_dB_all.argtypes = [POINTER(snd_mixer_elem_t), c_long, c_int]

# /usr/include/alsa/mixer.h:267
snd_mixer_selem_set_playback_switch = _lib.snd_mixer_selem_set_playback_switch
snd_mixer_selem_set_playback_switch.restype = c_int
snd_mixer_selem_set_playback_switch.argtypes = [POINTER(snd_mixer_elem_t), snd_mixer_selem_channel_id_t, c_int]

# /usr/include/alsa/mixer.h:268
snd_mixer_selem_set_capture_switch = _lib.snd_mixer_selem_set_capture_switch
snd_mixer_selem_set_capture_switch.restype = c_int
snd_mixer_selem_set_capture_switch.argtypes = [POINTER(snd_mixer_elem_t), snd_mixer_selem_channel_id_t, c_int]

# /usr/include/alsa/mixer.h:269
snd_mixer_selem_set_playback_switch_all = _lib.snd_mixer_selem_set_playback_switch_all
snd_mixer_selem_set_playback_switch_all.restype = c_int
snd_mixer_selem_set_playback_switch_all.argtypes = [POINTER(snd_mixer_elem_t), c_int]

# /usr/include/alsa/mixer.h:270
snd_mixer_selem_set_capture_switch_all = _lib.snd_mixer_selem_set_capture_switch_all
snd_mixer_selem_set_capture_switch_all.restype = c_int
snd_mixer_selem_set_capture_switch_all.argtypes = [POINTER(snd_mixer_elem_t), c_int]

# /usr/include/alsa/mixer.h:271
snd_mixer_selem_get_playback_volume_range = _lib.snd_mixer_selem_get_playback_volume_range
snd_mixer_selem_get_playback_volume_range.restype = c_int
snd_mixer_selem_get_playback_volume_range.argtypes = [POINTER(snd_mixer_elem_t), POINTER(c_long), POINTER(c_long)]

# /usr/include/alsa/mixer.h:273
snd_mixer_selem_get_playback_dB_range = _lib.snd_mixer_selem_get_playback_dB_range
snd_mixer_selem_get_playback_dB_range.restype = c_int
snd_mixer_selem_get_playback_dB_range.argtypes = [POINTER(snd_mixer_elem_t), POINTER(c_long), POINTER(c_long)]

# /usr/include/alsa/mixer.h:275
snd_mixer_selem_set_playback_volume_range = _lib.snd_mixer_selem_set_playback_volume_range
snd_mixer_selem_set_playback_volume_range.restype = c_int
snd_mixer_selem_set_playback_volume_range.argtypes = [POINTER(snd_mixer_elem_t), c_long, c_long]

# /usr/include/alsa/mixer.h:277
snd_mixer_selem_get_capture_volume_range = _lib.snd_mixer_selem_get_capture_volume_range
snd_mixer_selem_get_capture_volume_range.restype = c_int
snd_mixer_selem_get_capture_volume_range.argtypes = [POINTER(snd_mixer_elem_t), POINTER(c_long), POINTER(c_long)]

# /usr/include/alsa/mixer.h:279
snd_mixer_selem_get_capture_dB_range = _lib.snd_mixer_selem_get_capture_dB_range
snd_mixer_selem_get_capture_dB_range.restype = c_int
snd_mixer_selem_get_capture_dB_range.argtypes = [POINTER(snd_mixer_elem_t), POINTER(c_long), POINTER(c_long)]

# /usr/include/alsa/mixer.h:281
snd_mixer_selem_set_capture_volume_range = _lib.snd_mixer_selem_set_capture_volume_range
snd_mixer_selem_set_capture_volume_range.restype = c_int
snd_mixer_selem_set_capture_volume_range.argtypes = [POINTER(snd_mixer_elem_t), c_long, c_long]

# /usr/include/alsa/mixer.h:284
snd_mixer_selem_is_enumerated = _lib.snd_mixer_selem_is_enumerated
snd_mixer_selem_is_enumerated.restype = c_int
snd_mixer_selem_is_enumerated.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:285
snd_mixer_selem_is_enum_playback = _lib.snd_mixer_selem_is_enum_playback
snd_mixer_selem_is_enum_playback.restype = c_int
snd_mixer_selem_is_enum_playback.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:286
snd_mixer_selem_is_enum_capture = _lib.snd_mixer_selem_is_enum_capture
snd_mixer_selem_is_enum_capture.restype = c_int
snd_mixer_selem_is_enum_capture.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:287
snd_mixer_selem_get_enum_items = _lib.snd_mixer_selem_get_enum_items
snd_mixer_selem_get_enum_items.restype = c_int
snd_mixer_selem_get_enum_items.argtypes = [POINTER(snd_mixer_elem_t)]

# /usr/include/alsa/mixer.h:288
snd_mixer_selem_get_enum_item_name = _lib.snd_mixer_selem_get_enum_item_name
snd_mixer_selem_get_enum_item_name.restype = c_int
snd_mixer_selem_get_enum_item_name.argtypes = [POINTER(snd_mixer_elem_t), c_uint, c_size_t, c_char_p]

# /usr/include/alsa/mixer.h:289
snd_mixer_selem_get_enum_item = _lib.snd_mixer_selem_get_enum_item
snd_mixer_selem_get_enum_item.restype = c_int
snd_mixer_selem_get_enum_item.argtypes = [POINTER(snd_mixer_elem_t), snd_mixer_selem_channel_id_t, POINTER(c_uint)]

# /usr/include/alsa/mixer.h:290
snd_mixer_selem_set_enum_item = _lib.snd_mixer_selem_set_enum_item
snd_mixer_selem_set_enum_item.restype = c_int
snd_mixer_selem_set_enum_item.argtypes = [POINTER(snd_mixer_elem_t), snd_mixer_selem_channel_id_t, c_uint]

# /usr/include/alsa/mixer.h:292
snd_mixer_selem_id_sizeof = _lib.snd_mixer_selem_id_sizeof
snd_mixer_selem_id_sizeof.restype = c_size_t
snd_mixer_selem_id_sizeof.argtypes = []

# /usr/include/alsa/mixer.h:298
snd_mixer_selem_id_malloc = _lib.snd_mixer_selem_id_malloc
snd_mixer_selem_id_malloc.restype = c_int
snd_mixer_selem_id_malloc.argtypes = [POINTER(POINTER(snd_mixer_selem_id_t))]

# /usr/include/alsa/mixer.h:299
snd_mixer_selem_id_free = _lib.snd_mixer_selem_id_free
snd_mixer_selem_id_free.restype = None
snd_mixer_selem_id_free.argtypes = [POINTER(snd_mixer_selem_id_t)]

# /usr/include/alsa/mixer.h:300
snd_mixer_selem_id_copy = _lib.snd_mixer_selem_id_copy
snd_mixer_selem_id_copy.restype = None
snd_mixer_selem_id_copy.argtypes = [POINTER(snd_mixer_selem_id_t), POINTER(snd_mixer_selem_id_t)]

# /usr/include/alsa/mixer.h:301
snd_mixer_selem_id_get_name = _lib.snd_mixer_selem_id_get_name
snd_mixer_selem_id_get_name.restype = c_char_p
snd_mixer_selem_id_get_name.argtypes = [POINTER(snd_mixer_selem_id_t)]

# /usr/include/alsa/mixer.h:302
snd_mixer_selem_id_get_index = _lib.snd_mixer_selem_id_get_index
snd_mixer_selem_id_get_index.restype = c_uint
snd_mixer_selem_id_get_index.argtypes = [POINTER(snd_mixer_selem_id_t)]

# /usr/include/alsa/mixer.h:303
snd_mixer_selem_id_set_name = _lib.snd_mixer_selem_id_set_name
snd_mixer_selem_id_set_name.restype = None
snd_mixer_selem_id_set_name.argtypes = [POINTER(snd_mixer_selem_id_t), c_char_p]

# /usr/include/alsa/mixer.h:304
snd_mixer_selem_id_set_index = _lib.snd_mixer_selem_id_set_index
snd_mixer_selem_id_set_index.restype = None
snd_mixer_selem_id_set_index.argtypes = [POINTER(snd_mixer_selem_id_t), c_uint]

snd_seq_event_type_t = c_ubyte 	# /usr/include/alsa/seq_event.h:41
class struct_snd_seq_addr(Structure):
    __slots__ = [
        'client',
        'port',
    ]
struct_snd_seq_addr._fields_ = [
    ('client', c_ubyte),
    ('port', c_ubyte),
]

snd_seq_addr_t = struct_snd_seq_addr 	# /usr/include/alsa/seq_event.h:239
class struct_snd_seq_connect(Structure):
    __slots__ = [
        'sender',
        'dest',
    ]
struct_snd_seq_connect._fields_ = [
    ('sender', snd_seq_addr_t),
    ('dest', snd_seq_addr_t),
]

snd_seq_connect_t = struct_snd_seq_connect 	# /usr/include/alsa/seq_event.h:245
class struct_snd_seq_real_time(Structure):
    __slots__ = [
        'tv_sec',
        'tv_nsec',
    ]
struct_snd_seq_real_time._fields_ = [
    ('tv_sec', c_uint),
    ('tv_nsec', c_uint),
]

snd_seq_real_time_t = struct_snd_seq_real_time 	# /usr/include/alsa/seq_event.h:252
snd_seq_tick_time_t = c_uint 	# /usr/include/alsa/seq_event.h:255
class struct_snd_seq_timestamp(Union):
    __slots__ = [
        'tick',
        'time',
    ]
struct_snd_seq_timestamp._fields_ = [
    ('tick', snd_seq_tick_time_t),
    ('time', struct_snd_seq_real_time),
]

snd_seq_timestamp_t = struct_snd_seq_timestamp 	# /usr/include/alsa/seq_event.h:261
SND_SEQ_TIME_STAMP_TICK = 0 	# /usr/include/alsa/seq_event.h:269
SND_SEQ_TIME_STAMP_REAL = 1 	# /usr/include/alsa/seq_event.h:270
SND_SEQ_TIME_STAMP_MASK = 1 	# /usr/include/alsa/seq_event.h:271
SND_SEQ_TIME_MODE_ABS = 0 	# /usr/include/alsa/seq_event.h:273
SND_SEQ_TIME_MODE_REL = 2 	# /usr/include/alsa/seq_event.h:274
SND_SEQ_TIME_MODE_MASK = 2 	# /usr/include/alsa/seq_event.h:275
SND_SEQ_EVENT_LENGTH_FIXED = 0 	# /usr/include/alsa/seq_event.h:277
SND_SEQ_EVENT_LENGTH_VARIABLE = 4 	# /usr/include/alsa/seq_event.h:278
SND_SEQ_EVENT_LENGTH_VARUSR = 8 	# /usr/include/alsa/seq_event.h:279
SND_SEQ_EVENT_LENGTH_MASK = 12 	# /usr/include/alsa/seq_event.h:280
SND_SEQ_PRIORITY_NORMAL = 0 	# /usr/include/alsa/seq_event.h:282
SND_SEQ_PRIORITY_HIGH = 16 	# /usr/include/alsa/seq_event.h:283
SND_SEQ_PRIORITY_MASK = 16 	# /usr/include/alsa/seq_event.h:284
class struct_snd_seq_ev_note(Structure):
    __slots__ = [
        'channel',
        'note',
        'velocity',
        'off_velocity',
        'duration',
    ]
struct_snd_seq_ev_note._fields_ = [
    ('channel', c_ubyte),
    ('note', c_ubyte),
    ('velocity', c_ubyte),
    ('off_velocity', c_ubyte),
    ('duration', c_uint),
]

snd_seq_ev_note_t = struct_snd_seq_ev_note 	# /usr/include/alsa/seq_event.h:294
class struct_snd_seq_ev_ctrl(Structure):
    __slots__ = [
        'channel',
        'unused',
        'param',
        'value',
    ]
struct_snd_seq_ev_ctrl._fields_ = [
    ('channel', c_ubyte),
    ('unused', c_ubyte * 3),
    ('param', c_uint),
    ('value', c_int),
]

snd_seq_ev_ctrl_t = struct_snd_seq_ev_ctrl 	# /usr/include/alsa/seq_event.h:302
class struct_snd_seq_ev_raw8(Structure):
    __slots__ = [
        'd',
    ]
struct_snd_seq_ev_raw8._fields_ = [
    ('d', c_ubyte * 12),
]

snd_seq_ev_raw8_t = struct_snd_seq_ev_raw8 	# /usr/include/alsa/seq_event.h:307
class struct_snd_seq_ev_raw32(Structure):
    __slots__ = [
        'd',
    ]
struct_snd_seq_ev_raw32._fields_ = [
    ('d', c_uint * 3),
]

snd_seq_ev_raw32_t = struct_snd_seq_ev_raw32 	# /usr/include/alsa/seq_event.h:312
class struct_snd_seq_ev_ext(Structure):
    __slots__ = [
        'len',
        'ptr',
    ]
struct_snd_seq_ev_ext._fields_ = [
    ('len', c_uint),
    ('ptr', POINTER(None)),
]

snd_seq_ev_ext_t = struct_snd_seq_ev_ext 	# /usr/include/alsa/seq_event.h:318
snd_seq_instr_cluster_t = c_uint 	# /usr/include/alsa/seq_event.h:321
class struct_snd_seq_instr(Structure):
    __slots__ = [
        'cluster',
        'std',
        'bank',
        'prg',
    ]
struct_snd_seq_instr._fields_ = [
    ('cluster', snd_seq_instr_cluster_t),
    ('std', c_uint),
    ('bank', c_ushort),
    ('prg', c_ushort),
]

snd_seq_instr_t = struct_snd_seq_instr 	# /usr/include/alsa/seq_event.h:329
class struct_snd_seq_ev_sample(Structure):
    __slots__ = [
        'std',
        'bank',
        'prg',
    ]
struct_snd_seq_ev_sample._fields_ = [
    ('std', c_uint),
    ('bank', c_ushort),
    ('prg', c_ushort),
]

snd_seq_ev_sample_t = struct_snd_seq_ev_sample 	# /usr/include/alsa/seq_event.h:336
class struct_snd_seq_ev_cluster(Structure):
    __slots__ = [
        'cluster',
    ]
struct_snd_seq_ev_cluster._fields_ = [
    ('cluster', snd_seq_instr_cluster_t),
]

snd_seq_ev_cluster_t = struct_snd_seq_ev_cluster 	# /usr/include/alsa/seq_event.h:341
snd_seq_position_t = c_uint 	# /usr/include/alsa/seq_event.h:344
enum_snd_seq_stop_mode = c_int
SND_SEQ_SAMPLE_STOP_IMMEDIATELY = 0
SND_SEQ_SAMPLE_STOP_VENVELOPE = 1
SND_SEQ_SAMPLE_STOP_LOOP = 2
snd_seq_stop_mode_t = enum_snd_seq_stop_mode 	# /usr/include/alsa/seq_event.h:351
snd_seq_frequency_t = c_int 	# /usr/include/alsa/seq_event.h:354
class struct_snd_seq_ev_volume(Structure):
    __slots__ = [
        'volume',
        'lr',
        'fr',
        'du',
    ]
struct_snd_seq_ev_volume._fields_ = [
    ('volume', c_short),
    ('lr', c_short),
    ('fr', c_short),
    ('du', c_short),
]

snd_seq_ev_volume_t = struct_snd_seq_ev_volume 	# /usr/include/alsa/seq_event.h:362
class struct_snd_seq_ev_loop(Structure):
    __slots__ = [
        'start',
        'end',
    ]
struct_snd_seq_ev_loop._fields_ = [
    ('start', c_uint),
    ('end', c_uint),
]

snd_seq_ev_loop_t = struct_snd_seq_ev_loop 	# /usr/include/alsa/seq_event.h:368
class struct_snd_seq_ev_sample_control(Structure):
    __slots__ = [
        'channel',
        'unused',
        'param',
    ]
class struct_anon_27(Union):
    __slots__ = [
        'sample',
        'cluster',
        'position',
        'stop_mode',
        'frequency',
        'volume',
        'loop',
        'raw8',
    ]
struct_anon_27._fields_ = [
    ('sample', snd_seq_ev_sample_t),
    ('cluster', snd_seq_ev_cluster_t),
    ('position', snd_seq_position_t),
    ('stop_mode', snd_seq_stop_mode_t),
    ('frequency', snd_seq_frequency_t),
    ('volume', snd_seq_ev_volume_t),
    ('loop', snd_seq_ev_loop_t),
    ('raw8', c_ubyte * 8),
]

struct_snd_seq_ev_sample_control._fields_ = [
    ('channel', c_ubyte),
    ('unused', c_ubyte * 3),
    ('param', struct_anon_27),
]

snd_seq_ev_sample_control_t = struct_snd_seq_ev_sample_control 	# /usr/include/alsa/seq_event.h:384
class struct_snd_seq_ev_instr_begin(Structure):
    __slots__ = [
        'timeout',
    ]
struct_snd_seq_ev_instr_begin._fields_ = [
    ('timeout', c_int),
]

snd_seq_ev_instr_begin_t = struct_snd_seq_ev_instr_begin 	# /usr/include/alsa/seq_event.h:391
class struct_snd_seq_result(Structure):
    __slots__ = [
        'event',
        'result',
    ]
struct_snd_seq_result._fields_ = [
    ('event', c_int),
    ('result', c_int),
]

snd_seq_result_t = struct_snd_seq_result 	# /usr/include/alsa/seq_event.h:397
class struct_snd_seq_queue_skew(Structure):
    __slots__ = [
        'value',
        'base',
    ]
struct_snd_seq_queue_skew._fields_ = [
    ('value', c_uint),
    ('base', c_uint),
]

snd_seq_queue_skew_t = struct_snd_seq_queue_skew 	# /usr/include/alsa/seq_event.h:403
class struct_snd_seq_ev_queue_control(Structure):
    __slots__ = [
        'queue',
        'unused',
        'param',
    ]
class struct_anon_28(Union):
    __slots__ = [
        'value',
        'time',
        'position',
        'skew',
        'd32',
        'd8',
    ]
struct_anon_28._fields_ = [
    ('value', c_int),
    ('time', snd_seq_timestamp_t),
    ('position', c_uint),
    ('skew', snd_seq_queue_skew_t),
    ('d32', c_uint * 2),
    ('d8', c_ubyte * 8),
]

struct_snd_seq_ev_queue_control._fields_ = [
    ('queue', c_ubyte),
    ('unused', c_ubyte * 3),
    ('param', struct_anon_28),
]

snd_seq_ev_queue_control_t = struct_snd_seq_ev_queue_control 	# /usr/include/alsa/seq_event.h:417
class struct_snd_seq_event(Structure):
    __slots__ = [
        'type',
        'flags',
        'tag',
        'queue',
        'time',
        'source',
        'dest',
        'data',
    ]
class struct_anon_29(Union):
    __slots__ = [
        'note',
        'control',
        'raw8',
        'raw32',
        'ext',
        'queue',
        'time',
        'addr',
        'connect',
        'result',
        'instr_begin',
        'sample',
    ]
struct_anon_29._fields_ = [
    ('note', snd_seq_ev_note_t),
    ('control', snd_seq_ev_ctrl_t),
    ('raw8', snd_seq_ev_raw8_t),
    ('raw32', snd_seq_ev_raw32_t),
    ('ext', snd_seq_ev_ext_t),
    ('queue', snd_seq_ev_queue_control_t),
    ('time', snd_seq_timestamp_t),
    ('addr', snd_seq_addr_t),
    ('connect', snd_seq_connect_t),
    ('result', snd_seq_result_t),
    ('instr_begin', snd_seq_ev_instr_begin_t),
    ('sample', snd_seq_ev_sample_control_t),
]

struct_snd_seq_event._fields_ = [
    ('type', snd_seq_event_type_t),
    ('flags', c_ubyte),
    ('tag', c_ubyte),
    ('queue', c_ubyte),
    ('time', snd_seq_timestamp_t),
    ('source', snd_seq_addr_t),
    ('dest', snd_seq_addr_t),
    ('data', struct_anon_29),
]

snd_seq_event_t = struct_snd_seq_event 	# /usr/include/alsa/seq_event.h:446
SND_SEQ_DLSYM_VERSION = 0 	# /usr/include/alsa/seq.h:44
class struct__snd_seq(Structure):
    __slots__ = [
    ]
struct__snd_seq._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_seq(Structure):
    __slots__ = [
    ]
struct__snd_seq._fields_ = [
    ('_opaque_struct', c_int)
]

snd_seq_t = struct__snd_seq 	# /usr/include/alsa/seq.h:47
SND_SEQ_OPEN_OUTPUT = 1 	# /usr/include/alsa/seq.h:61
SND_SEQ_OPEN_INPUT = 2 	# /usr/include/alsa/seq.h:62
SND_SEQ_OPEN_DUPLEX = 3 	# /usr/include/alsa/seq.h:63
SND_SEQ_NONBLOCK = 1 	# /usr/include/alsa/seq.h:68
enum__snd_seq_type = c_int
SND_SEQ_TYPE_HW = 1
SND_SEQ_TYPE_SHM = 2
SND_SEQ_TYPE_INET = 3
snd_seq_type_t = enum__snd_seq_type 	# /usr/include/alsa/seq.h:75
SND_SEQ_ADDRESS_UNKNOWN = 253 	# /usr/include/alsa/seq.h:78
SND_SEQ_ADDRESS_SUBSCRIBERS = 254 	# /usr/include/alsa/seq.h:79
SND_SEQ_ADDRESS_BROADCAST = 255 	# /usr/include/alsa/seq.h:80
SND_SEQ_CLIENT_SYSTEM = 0 	# /usr/include/alsa/seq.h:83
# /usr/include/alsa/seq.h:87
snd_seq_open = _lib.snd_seq_open
snd_seq_open.restype = c_int
snd_seq_open.argtypes = [POINTER(POINTER(snd_seq_t)), c_char_p, c_int, c_int]

# /usr/include/alsa/seq.h:88
snd_seq_open_lconf = _lib.snd_seq_open_lconf
snd_seq_open_lconf.restype = c_int
snd_seq_open_lconf.argtypes = [POINTER(POINTER(snd_seq_t)), c_char_p, c_int, c_int, POINTER(snd_config_t)]

# /usr/include/alsa/seq.h:89
snd_seq_name = _lib.snd_seq_name
snd_seq_name.restype = c_char_p
snd_seq_name.argtypes = [POINTER(snd_seq_t)]

# /usr/include/alsa/seq.h:90
snd_seq_type = _lib.snd_seq_type
snd_seq_type.restype = snd_seq_type_t
snd_seq_type.argtypes = [POINTER(snd_seq_t)]

# /usr/include/alsa/seq.h:91
snd_seq_close = _lib.snd_seq_close
snd_seq_close.restype = c_int
snd_seq_close.argtypes = [POINTER(snd_seq_t)]

# /usr/include/alsa/seq.h:92
snd_seq_poll_descriptors_count = _lib.snd_seq_poll_descriptors_count
snd_seq_poll_descriptors_count.restype = c_int
snd_seq_poll_descriptors_count.argtypes = [POINTER(snd_seq_t), c_short]

class struct_pollfd(Structure):
    __slots__ = [
    ]
struct_pollfd._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/seq.h:93
snd_seq_poll_descriptors = _lib.snd_seq_poll_descriptors
snd_seq_poll_descriptors.restype = c_int
snd_seq_poll_descriptors.argtypes = [POINTER(snd_seq_t), POINTER(struct_pollfd), c_uint, c_short]

class struct_pollfd(Structure):
    __slots__ = [
    ]
struct_pollfd._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/alsa/seq.h:94
snd_seq_poll_descriptors_revents = _lib.snd_seq_poll_descriptors_revents
snd_seq_poll_descriptors_revents.restype = c_int
snd_seq_poll_descriptors_revents.argtypes = [POINTER(snd_seq_t), POINTER(struct_pollfd), c_uint, POINTER(c_ushort)]

# /usr/include/alsa/seq.h:95
snd_seq_nonblock = _lib.snd_seq_nonblock
snd_seq_nonblock.restype = c_int
snd_seq_nonblock.argtypes = [POINTER(snd_seq_t), c_int]

# /usr/include/alsa/seq.h:96
snd_seq_client_id = _lib.snd_seq_client_id
snd_seq_client_id.restype = c_int
snd_seq_client_id.argtypes = [POINTER(snd_seq_t)]

# /usr/include/alsa/seq.h:98
snd_seq_get_output_buffer_size = _lib.snd_seq_get_output_buffer_size
snd_seq_get_output_buffer_size.restype = c_size_t
snd_seq_get_output_buffer_size.argtypes = [POINTER(snd_seq_t)]

# /usr/include/alsa/seq.h:99
snd_seq_get_input_buffer_size = _lib.snd_seq_get_input_buffer_size
snd_seq_get_input_buffer_size.restype = c_size_t
snd_seq_get_input_buffer_size.argtypes = [POINTER(snd_seq_t)]

# /usr/include/alsa/seq.h:100
snd_seq_set_output_buffer_size = _lib.snd_seq_set_output_buffer_size
snd_seq_set_output_buffer_size.restype = c_int
snd_seq_set_output_buffer_size.argtypes = [POINTER(snd_seq_t), c_size_t]

# /usr/include/alsa/seq.h:101
snd_seq_set_input_buffer_size = _lib.snd_seq_set_input_buffer_size
snd_seq_set_input_buffer_size.restype = c_int
snd_seq_set_input_buffer_size.argtypes = [POINTER(snd_seq_t), c_size_t]

class struct__snd_seq_system_info(Structure):
    __slots__ = [
    ]
struct__snd_seq_system_info._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_seq_system_info(Structure):
    __slots__ = [
    ]
struct__snd_seq_system_info._fields_ = [
    ('_opaque_struct', c_int)
]

snd_seq_system_info_t = struct__snd_seq_system_info 	# /usr/include/alsa/seq.h:104
# /usr/include/alsa/seq.h:106
snd_seq_system_info_sizeof = _lib.snd_seq_system_info_sizeof
snd_seq_system_info_sizeof.restype = c_size_t
snd_seq_system_info_sizeof.argtypes = []

# /usr/include/alsa/seq.h:110
snd_seq_system_info_malloc = _lib.snd_seq_system_info_malloc
snd_seq_system_info_malloc.restype = c_int
snd_seq_system_info_malloc.argtypes = [POINTER(POINTER(snd_seq_system_info_t))]

# /usr/include/alsa/seq.h:111
snd_seq_system_info_free = _lib.snd_seq_system_info_free
snd_seq_system_info_free.restype = None
snd_seq_system_info_free.argtypes = [POINTER(snd_seq_system_info_t)]

# /usr/include/alsa/seq.h:112
snd_seq_system_info_copy = _lib.snd_seq_system_info_copy
snd_seq_system_info_copy.restype = None
snd_seq_system_info_copy.argtypes = [POINTER(snd_seq_system_info_t), POINTER(snd_seq_system_info_t)]

# /usr/include/alsa/seq.h:114
snd_seq_system_info_get_queues = _lib.snd_seq_system_info_get_queues
snd_seq_system_info_get_queues.restype = c_int
snd_seq_system_info_get_queues.argtypes = [POINTER(snd_seq_system_info_t)]

# /usr/include/alsa/seq.h:115
snd_seq_system_info_get_clients = _lib.snd_seq_system_info_get_clients
snd_seq_system_info_get_clients.restype = c_int
snd_seq_system_info_get_clients.argtypes = [POINTER(snd_seq_system_info_t)]

# /usr/include/alsa/seq.h:116
snd_seq_system_info_get_ports = _lib.snd_seq_system_info_get_ports
snd_seq_system_info_get_ports.restype = c_int
snd_seq_system_info_get_ports.argtypes = [POINTER(snd_seq_system_info_t)]

# /usr/include/alsa/seq.h:117
snd_seq_system_info_get_channels = _lib.snd_seq_system_info_get_channels
snd_seq_system_info_get_channels.restype = c_int
snd_seq_system_info_get_channels.argtypes = [POINTER(snd_seq_system_info_t)]

# /usr/include/alsa/seq.h:118
snd_seq_system_info_get_cur_clients = _lib.snd_seq_system_info_get_cur_clients
snd_seq_system_info_get_cur_clients.restype = c_int
snd_seq_system_info_get_cur_clients.argtypes = [POINTER(snd_seq_system_info_t)]

# /usr/include/alsa/seq.h:119
snd_seq_system_info_get_cur_queues = _lib.snd_seq_system_info_get_cur_queues
snd_seq_system_info_get_cur_queues.restype = c_int
snd_seq_system_info_get_cur_queues.argtypes = [POINTER(snd_seq_system_info_t)]

# /usr/include/alsa/seq.h:121
snd_seq_system_info = _lib.snd_seq_system_info
snd_seq_system_info.restype = c_int
snd_seq_system_info.argtypes = [POINTER(snd_seq_t), POINTER(snd_seq_system_info_t)]

class struct__snd_seq_client_info(Structure):
    __slots__ = [
    ]
struct__snd_seq_client_info._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_seq_client_info(Structure):
    __slots__ = [
    ]
struct__snd_seq_client_info._fields_ = [
    ('_opaque_struct', c_int)
]

snd_seq_client_info_t = struct__snd_seq_client_info 	# /usr/include/alsa/seq.h:134
enum_snd_seq_client_type = c_int
SND_SEQ_USER_CLIENT = 1
SND_SEQ_KERNEL_CLIENT = 2
snd_seq_client_type_t = enum_snd_seq_client_type 	# /usr/include/alsa/seq.h:140
# /usr/include/alsa/seq.h:142
snd_seq_client_info_sizeof = _lib.snd_seq_client_info_sizeof
snd_seq_client_info_sizeof.restype = c_size_t
snd_seq_client_info_sizeof.argtypes = []

# /usr/include/alsa/seq.h:146
snd_seq_client_info_malloc = _lib.snd_seq_client_info_malloc
snd_seq_client_info_malloc.restype = c_int
snd_seq_client_info_malloc.argtypes = [POINTER(POINTER(snd_seq_client_info_t))]

# /usr/include/alsa/seq.h:147
snd_seq_client_info_free = _lib.snd_seq_client_info_free
snd_seq_client_info_free.restype = None
snd_seq_client_info_free.argtypes = [POINTER(snd_seq_client_info_t)]

# /usr/include/alsa/seq.h:148
snd_seq_client_info_copy = _lib.snd_seq_client_info_copy
snd_seq_client_info_copy.restype = None
snd_seq_client_info_copy.argtypes = [POINTER(snd_seq_client_info_t), POINTER(snd_seq_client_info_t)]

# /usr/include/alsa/seq.h:150
snd_seq_client_info_get_client = _lib.snd_seq_client_info_get_client
snd_seq_client_info_get_client.restype = c_int
snd_seq_client_info_get_client.argtypes = [POINTER(snd_seq_client_info_t)]

# /usr/include/alsa/seq.h:151
snd_seq_client_info_get_type = _lib.snd_seq_client_info_get_type
snd_seq_client_info_get_type.restype = snd_seq_client_type_t
snd_seq_client_info_get_type.argtypes = [POINTER(snd_seq_client_info_t)]

# /usr/include/alsa/seq.h:152
snd_seq_client_info_get_name = _lib.snd_seq_client_info_get_name
snd_seq_client_info_get_name.restype = c_char_p
snd_seq_client_info_get_name.argtypes = [POINTER(snd_seq_client_info_t)]

# /usr/include/alsa/seq.h:153
snd_seq_client_info_get_broadcast_filter = _lib.snd_seq_client_info_get_broadcast_filter
snd_seq_client_info_get_broadcast_filter.restype = c_int
snd_seq_client_info_get_broadcast_filter.argtypes = [POINTER(snd_seq_client_info_t)]

# /usr/include/alsa/seq.h:154
snd_seq_client_info_get_error_bounce = _lib.snd_seq_client_info_get_error_bounce
snd_seq_client_info_get_error_bounce.restype = c_int
snd_seq_client_info_get_error_bounce.argtypes = [POINTER(snd_seq_client_info_t)]

# /usr/include/alsa/seq.h:155
snd_seq_client_info_get_event_filter = _lib.snd_seq_client_info_get_event_filter
snd_seq_client_info_get_event_filter.restype = POINTER(c_ubyte)
snd_seq_client_info_get_event_filter.argtypes = [POINTER(snd_seq_client_info_t)]

# /usr/include/alsa/seq.h:156
snd_seq_client_info_get_num_ports = _lib.snd_seq_client_info_get_num_ports
snd_seq_client_info_get_num_ports.restype = c_int
snd_seq_client_info_get_num_ports.argtypes = [POINTER(snd_seq_client_info_t)]

# /usr/include/alsa/seq.h:157
snd_seq_client_info_get_event_lost = _lib.snd_seq_client_info_get_event_lost
snd_seq_client_info_get_event_lost.restype = c_int
snd_seq_client_info_get_event_lost.argtypes = [POINTER(snd_seq_client_info_t)]

# /usr/include/alsa/seq.h:159
snd_seq_client_info_set_client = _lib.snd_seq_client_info_set_client
snd_seq_client_info_set_client.restype = None
snd_seq_client_info_set_client.argtypes = [POINTER(snd_seq_client_info_t), c_int]

# /usr/include/alsa/seq.h:160
snd_seq_client_info_set_name = _lib.snd_seq_client_info_set_name
snd_seq_client_info_set_name.restype = None
snd_seq_client_info_set_name.argtypes = [POINTER(snd_seq_client_info_t), c_char_p]

# /usr/include/alsa/seq.h:161
snd_seq_client_info_set_broadcast_filter = _lib.snd_seq_client_info_set_broadcast_filter
snd_seq_client_info_set_broadcast_filter.restype = None
snd_seq_client_info_set_broadcast_filter.argtypes = [POINTER(snd_seq_client_info_t), c_int]

# /usr/include/alsa/seq.h:162
snd_seq_client_info_set_error_bounce = _lib.snd_seq_client_info_set_error_bounce
snd_seq_client_info_set_error_bounce.restype = None
snd_seq_client_info_set_error_bounce.argtypes = [POINTER(snd_seq_client_info_t), c_int]

# /usr/include/alsa/seq.h:163
snd_seq_client_info_set_event_filter = _lib.snd_seq_client_info_set_event_filter
snd_seq_client_info_set_event_filter.restype = None
snd_seq_client_info_set_event_filter.argtypes = [POINTER(snd_seq_client_info_t), POINTER(c_ubyte)]

# /usr/include/alsa/seq.h:165
snd_seq_get_client_info = _lib.snd_seq_get_client_info
snd_seq_get_client_info.restype = c_int
snd_seq_get_client_info.argtypes = [POINTER(snd_seq_t), POINTER(snd_seq_client_info_t)]

# /usr/include/alsa/seq.h:166
snd_seq_get_any_client_info = _lib.snd_seq_get_any_client_info
snd_seq_get_any_client_info.restype = c_int
snd_seq_get_any_client_info.argtypes = [POINTER(snd_seq_t), c_int, POINTER(snd_seq_client_info_t)]

# /usr/include/alsa/seq.h:167
snd_seq_set_client_info = _lib.snd_seq_set_client_info
snd_seq_set_client_info.restype = c_int
snd_seq_set_client_info.argtypes = [POINTER(snd_seq_t), POINTER(snd_seq_client_info_t)]

# /usr/include/alsa/seq.h:168
snd_seq_query_next_client = _lib.snd_seq_query_next_client
snd_seq_query_next_client.restype = c_int
snd_seq_query_next_client.argtypes = [POINTER(snd_seq_t), POINTER(snd_seq_client_info_t)]

class struct__snd_seq_client_pool(Structure):
    __slots__ = [
    ]
struct__snd_seq_client_pool._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_seq_client_pool(Structure):
    __slots__ = [
    ]
struct__snd_seq_client_pool._fields_ = [
    ('_opaque_struct', c_int)
]

snd_seq_client_pool_t = struct__snd_seq_client_pool 	# /usr/include/alsa/seq.h:174
# /usr/include/alsa/seq.h:176
snd_seq_client_pool_sizeof = _lib.snd_seq_client_pool_sizeof
snd_seq_client_pool_sizeof.restype = c_size_t
snd_seq_client_pool_sizeof.argtypes = []

# /usr/include/alsa/seq.h:180
snd_seq_client_pool_malloc = _lib.snd_seq_client_pool_malloc
snd_seq_client_pool_malloc.restype = c_int
snd_seq_client_pool_malloc.argtypes = [POINTER(POINTER(snd_seq_client_pool_t))]

# /usr/include/alsa/seq.h:181
snd_seq_client_pool_free = _lib.snd_seq_client_pool_free
snd_seq_client_pool_free.restype = None
snd_seq_client_pool_free.argtypes = [POINTER(snd_seq_client_pool_t)]

# /usr/include/alsa/seq.h:182
snd_seq_client_pool_copy = _lib.snd_seq_client_pool_copy
snd_seq_client_pool_copy.restype = None
snd_seq_client_pool_copy.argtypes = [POINTER(snd_seq_client_pool_t), POINTER(snd_seq_client_pool_t)]

# /usr/include/alsa/seq.h:184
snd_seq_client_pool_get_client = _lib.snd_seq_client_pool_get_client
snd_seq_client_pool_get_client.restype = c_int
snd_seq_client_pool_get_client.argtypes = [POINTER(snd_seq_client_pool_t)]

# /usr/include/alsa/seq.h:185
snd_seq_client_pool_get_output_pool = _lib.snd_seq_client_pool_get_output_pool
snd_seq_client_pool_get_output_pool.restype = c_size_t
snd_seq_client_pool_get_output_pool.argtypes = [POINTER(snd_seq_client_pool_t)]

# /usr/include/alsa/seq.h:186
snd_seq_client_pool_get_input_pool = _lib.snd_seq_client_pool_get_input_pool
snd_seq_client_pool_get_input_pool.restype = c_size_t
snd_seq_client_pool_get_input_pool.argtypes = [POINTER(snd_seq_client_pool_t)]

# /usr/include/alsa/seq.h:187
snd_seq_client_pool_get_output_room = _lib.snd_seq_client_pool_get_output_room
snd_seq_client_pool_get_output_room.restype = c_size_t
snd_seq_client_pool_get_output_room.argtypes = [POINTER(snd_seq_client_pool_t)]

# /usr/include/alsa/seq.h:188
snd_seq_client_pool_get_output_free = _lib.snd_seq_client_pool_get_output_free
snd_seq_client_pool_get_output_free.restype = c_size_t
snd_seq_client_pool_get_output_free.argtypes = [POINTER(snd_seq_client_pool_t)]

# /usr/include/alsa/seq.h:189
snd_seq_client_pool_get_input_free = _lib.snd_seq_client_pool_get_input_free
snd_seq_client_pool_get_input_free.restype = c_size_t
snd_seq_client_pool_get_input_free.argtypes = [POINTER(snd_seq_client_pool_t)]

# /usr/include/alsa/seq.h:190
snd_seq_client_pool_set_output_pool = _lib.snd_seq_client_pool_set_output_pool
snd_seq_client_pool_set_output_pool.restype = None
snd_seq_client_pool_set_output_pool.argtypes = [POINTER(snd_seq_client_pool_t), c_size_t]

# /usr/include/alsa/seq.h:191
snd_seq_client_pool_set_input_pool = _lib.snd_seq_client_pool_set_input_pool
snd_seq_client_pool_set_input_pool.restype = None
snd_seq_client_pool_set_input_pool.argtypes = [POINTER(snd_seq_client_pool_t), c_size_t]

# /usr/include/alsa/seq.h:192
snd_seq_client_pool_set_output_room = _lib.snd_seq_client_pool_set_output_room
snd_seq_client_pool_set_output_room.restype = None
snd_seq_client_pool_set_output_room.argtypes = [POINTER(snd_seq_client_pool_t), c_size_t]

# /usr/include/alsa/seq.h:194
snd_seq_get_client_pool = _lib.snd_seq_get_client_pool
snd_seq_get_client_pool.restype = c_int
snd_seq_get_client_pool.argtypes = [POINTER(snd_seq_t), POINTER(snd_seq_client_pool_t)]

# /usr/include/alsa/seq.h:195
snd_seq_set_client_pool = _lib.snd_seq_set_client_pool
snd_seq_set_client_pool.restype = c_int
snd_seq_set_client_pool.argtypes = [POINTER(snd_seq_t), POINTER(snd_seq_client_pool_t)]

class struct__snd_seq_port_info(Structure):
    __slots__ = [
    ]
struct__snd_seq_port_info._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_seq_port_info(Structure):
    __slots__ = [
    ]
struct__snd_seq_port_info._fields_ = [
    ('_opaque_struct', c_int)
]

snd_seq_port_info_t = struct__snd_seq_port_info 	# /usr/include/alsa/seq.h:209
SND_SEQ_PORT_SYSTEM_TIMER = 0 	# /usr/include/alsa/seq.h:212
SND_SEQ_PORT_SYSTEM_ANNOUNCE = 1 	# /usr/include/alsa/seq.h:213
SND_SEQ_PORT_CAP_READ = 1 	# /usr/include/alsa/seq.h:216
SND_SEQ_PORT_CAP_WRITE = 2 	# /usr/include/alsa/seq.h:217
SND_SEQ_PORT_CAP_SYNC_READ = 4 	# /usr/include/alsa/seq.h:219
SND_SEQ_PORT_CAP_SYNC_WRITE = 8 	# /usr/include/alsa/seq.h:220
SND_SEQ_PORT_CAP_DUPLEX = 16 	# /usr/include/alsa/seq.h:222
SND_SEQ_PORT_CAP_SUBS_READ = 32 	# /usr/include/alsa/seq.h:224
SND_SEQ_PORT_CAP_SUBS_WRITE = 64 	# /usr/include/alsa/seq.h:225
SND_SEQ_PORT_CAP_NO_EXPORT = 128 	# /usr/include/alsa/seq.h:226
SND_SEQ_PORT_TYPE_SPECIFIC = 1 	# /usr/include/alsa/seq.h:230
SND_SEQ_PORT_TYPE_MIDI_GENERIC = 2 	# /usr/include/alsa/seq.h:232
SND_SEQ_PORT_TYPE_MIDI_GM = 4 	# /usr/include/alsa/seq.h:234
SND_SEQ_PORT_TYPE_MIDI_GS = 8 	# /usr/include/alsa/seq.h:236
SND_SEQ_PORT_TYPE_MIDI_XG = 16 	# /usr/include/alsa/seq.h:238
SND_SEQ_PORT_TYPE_MIDI_MT32 = 32 	# /usr/include/alsa/seq.h:240
SND_SEQ_PORT_TYPE_MIDI_GM2 = 64 	# /usr/include/alsa/seq.h:242
SND_SEQ_PORT_TYPE_SYNTH = 1024 	# /usr/include/alsa/seq.h:245
SND_SEQ_PORT_TYPE_DIRECT_SAMPLE = 2048 	# /usr/include/alsa/seq.h:248
SND_SEQ_PORT_TYPE_SAMPLE = 4096 	# /usr/include/alsa/seq.h:251
SND_SEQ_PORT_TYPE_HARDWARE = 65536 	# /usr/include/alsa/seq.h:253
SND_SEQ_PORT_TYPE_SOFTWARE = 131072 	# /usr/include/alsa/seq.h:255
SND_SEQ_PORT_TYPE_SYNTHESIZER = 262144 	# /usr/include/alsa/seq.h:257
SND_SEQ_PORT_TYPE_PORT = 524288 	# /usr/include/alsa/seq.h:260
SND_SEQ_PORT_TYPE_APPLICATION = 1048576 	# /usr/include/alsa/seq.h:262
# /usr/include/alsa/seq.h:265
snd_seq_port_info_sizeof = _lib.snd_seq_port_info_sizeof
snd_seq_port_info_sizeof.restype = c_size_t
snd_seq_port_info_sizeof.argtypes = []

# /usr/include/alsa/seq.h:269
snd_seq_port_info_malloc = _lib.snd_seq_port_info_malloc
snd_seq_port_info_malloc.restype = c_int
snd_seq_port_info_malloc.argtypes = [POINTER(POINTER(snd_seq_port_info_t))]

# /usr/include/alsa/seq.h:270
snd_seq_port_info_free = _lib.snd_seq_port_info_free
snd_seq_port_info_free.restype = None
snd_seq_port_info_free.argtypes = [POINTER(snd_seq_port_info_t)]

# /usr/include/alsa/seq.h:271
snd_seq_port_info_copy = _lib.snd_seq_port_info_copy
snd_seq_port_info_copy.restype = None
snd_seq_port_info_copy.argtypes = [POINTER(snd_seq_port_info_t), POINTER(snd_seq_port_info_t)]

# /usr/include/alsa/seq.h:273
snd_seq_port_info_get_client = _lib.snd_seq_port_info_get_client
snd_seq_port_info_get_client.restype = c_int
snd_seq_port_info_get_client.argtypes = [POINTER(snd_seq_port_info_t)]

# /usr/include/alsa/seq.h:274
snd_seq_port_info_get_port = _lib.snd_seq_port_info_get_port
snd_seq_port_info_get_port.restype = c_int
snd_seq_port_info_get_port.argtypes = [POINTER(snd_seq_port_info_t)]

# /usr/include/alsa/seq.h:275
snd_seq_port_info_get_addr = _lib.snd_seq_port_info_get_addr
snd_seq_port_info_get_addr.restype = POINTER(snd_seq_addr_t)
snd_seq_port_info_get_addr.argtypes = [POINTER(snd_seq_port_info_t)]

# /usr/include/alsa/seq.h:276
snd_seq_port_info_get_name = _lib.snd_seq_port_info_get_name
snd_seq_port_info_get_name.restype = c_char_p
snd_seq_port_info_get_name.argtypes = [POINTER(snd_seq_port_info_t)]

# /usr/include/alsa/seq.h:277
snd_seq_port_info_get_capability = _lib.snd_seq_port_info_get_capability
snd_seq_port_info_get_capability.restype = c_uint
snd_seq_port_info_get_capability.argtypes = [POINTER(snd_seq_port_info_t)]

# /usr/include/alsa/seq.h:278
snd_seq_port_info_get_type = _lib.snd_seq_port_info_get_type
snd_seq_port_info_get_type.restype = c_uint
snd_seq_port_info_get_type.argtypes = [POINTER(snd_seq_port_info_t)]

# /usr/include/alsa/seq.h:279
snd_seq_port_info_get_midi_channels = _lib.snd_seq_port_info_get_midi_channels
snd_seq_port_info_get_midi_channels.restype = c_int
snd_seq_port_info_get_midi_channels.argtypes = [POINTER(snd_seq_port_info_t)]

# /usr/include/alsa/seq.h:280
snd_seq_port_info_get_midi_voices = _lib.snd_seq_port_info_get_midi_voices
snd_seq_port_info_get_midi_voices.restype = c_int
snd_seq_port_info_get_midi_voices.argtypes = [POINTER(snd_seq_port_info_t)]

# /usr/include/alsa/seq.h:281
snd_seq_port_info_get_synth_voices = _lib.snd_seq_port_info_get_synth_voices
snd_seq_port_info_get_synth_voices.restype = c_int
snd_seq_port_info_get_synth_voices.argtypes = [POINTER(snd_seq_port_info_t)]

# /usr/include/alsa/seq.h:282
snd_seq_port_info_get_read_use = _lib.snd_seq_port_info_get_read_use
snd_seq_port_info_get_read_use.restype = c_int
snd_seq_port_info_get_read_use.argtypes = [POINTER(snd_seq_port_info_t)]

# /usr/include/alsa/seq.h:283
snd_seq_port_info_get_write_use = _lib.snd_seq_port_info_get_write_use
snd_seq_port_info_get_write_use.restype = c_int
snd_seq_port_info_get_write_use.argtypes = [POINTER(snd_seq_port_info_t)]

# /usr/include/alsa/seq.h:284
snd_seq_port_info_get_port_specified = _lib.snd_seq_port_info_get_port_specified
snd_seq_port_info_get_port_specified.restype = c_int
snd_seq_port_info_get_port_specified.argtypes = [POINTER(snd_seq_port_info_t)]

# /usr/include/alsa/seq.h:285
snd_seq_port_info_get_timestamping = _lib.snd_seq_port_info_get_timestamping
snd_seq_port_info_get_timestamping.restype = c_int
snd_seq_port_info_get_timestamping.argtypes = [POINTER(snd_seq_port_info_t)]

# /usr/include/alsa/seq.h:286
snd_seq_port_info_get_timestamp_real = _lib.snd_seq_port_info_get_timestamp_real
snd_seq_port_info_get_timestamp_real.restype = c_int
snd_seq_port_info_get_timestamp_real.argtypes = [POINTER(snd_seq_port_info_t)]

# /usr/include/alsa/seq.h:287
snd_seq_port_info_get_timestamp_queue = _lib.snd_seq_port_info_get_timestamp_queue
snd_seq_port_info_get_timestamp_queue.restype = c_int
snd_seq_port_info_get_timestamp_queue.argtypes = [POINTER(snd_seq_port_info_t)]

# /usr/include/alsa/seq.h:289
snd_seq_port_info_set_client = _lib.snd_seq_port_info_set_client
snd_seq_port_info_set_client.restype = None
snd_seq_port_info_set_client.argtypes = [POINTER(snd_seq_port_info_t), c_int]

# /usr/include/alsa/seq.h:290
snd_seq_port_info_set_port = _lib.snd_seq_port_info_set_port
snd_seq_port_info_set_port.restype = None
snd_seq_port_info_set_port.argtypes = [POINTER(snd_seq_port_info_t), c_int]

# /usr/include/alsa/seq.h:291
snd_seq_port_info_set_addr = _lib.snd_seq_port_info_set_addr
snd_seq_port_info_set_addr.restype = None
snd_seq_port_info_set_addr.argtypes = [POINTER(snd_seq_port_info_t), POINTER(snd_seq_addr_t)]

# /usr/include/alsa/seq.h:292
snd_seq_port_info_set_name = _lib.snd_seq_port_info_set_name
snd_seq_port_info_set_name.restype = None
snd_seq_port_info_set_name.argtypes = [POINTER(snd_seq_port_info_t), c_char_p]

# /usr/include/alsa/seq.h:293
snd_seq_port_info_set_capability = _lib.snd_seq_port_info_set_capability
snd_seq_port_info_set_capability.restype = None
snd_seq_port_info_set_capability.argtypes = [POINTER(snd_seq_port_info_t), c_uint]

# /usr/include/alsa/seq.h:294
snd_seq_port_info_set_type = _lib.snd_seq_port_info_set_type
snd_seq_port_info_set_type.restype = None
snd_seq_port_info_set_type.argtypes = [POINTER(snd_seq_port_info_t), c_uint]

# /usr/include/alsa/seq.h:295
snd_seq_port_info_set_midi_channels = _lib.snd_seq_port_info_set_midi_channels
snd_seq_port_info_set_midi_channels.restype = None
snd_seq_port_info_set_midi_channels.argtypes = [POINTER(snd_seq_port_info_t), c_int]

# /usr/include/alsa/seq.h:296
snd_seq_port_info_set_midi_voices = _lib.snd_seq_port_info_set_midi_voices
snd_seq_port_info_set_midi_voices.restype = None
snd_seq_port_info_set_midi_voices.argtypes = [POINTER(snd_seq_port_info_t), c_int]

# /usr/include/alsa/seq.h:297
snd_seq_port_info_set_synth_voices = _lib.snd_seq_port_info_set_synth_voices
snd_seq_port_info_set_synth_voices.restype = None
snd_seq_port_info_set_synth_voices.argtypes = [POINTER(snd_seq_port_info_t), c_int]

# /usr/include/alsa/seq.h:298
snd_seq_port_info_set_port_specified = _lib.snd_seq_port_info_set_port_specified
snd_seq_port_info_set_port_specified.restype = None
snd_seq_port_info_set_port_specified.argtypes = [POINTER(snd_seq_port_info_t), c_int]

# /usr/include/alsa/seq.h:299
snd_seq_port_info_set_timestamping = _lib.snd_seq_port_info_set_timestamping
snd_seq_port_info_set_timestamping.restype = None
snd_seq_port_info_set_timestamping.argtypes = [POINTER(snd_seq_port_info_t), c_int]

# /usr/include/alsa/seq.h:300
snd_seq_port_info_set_timestamp_real = _lib.snd_seq_port_info_set_timestamp_real
snd_seq_port_info_set_timestamp_real.restype = None
snd_seq_port_info_set_timestamp_real.argtypes = [POINTER(snd_seq_port_info_t), c_int]

# /usr/include/alsa/seq.h:301
snd_seq_port_info_set_timestamp_queue = _lib.snd_seq_port_info_set_timestamp_queue
snd_seq_port_info_set_timestamp_queue.restype = None
snd_seq_port_info_set_timestamp_queue.argtypes = [POINTER(snd_seq_port_info_t), c_int]

# /usr/include/alsa/seq.h:303
snd_seq_create_port = _lib.snd_seq_create_port
snd_seq_create_port.restype = c_int
snd_seq_create_port.argtypes = [POINTER(snd_seq_t), POINTER(snd_seq_port_info_t)]

# /usr/include/alsa/seq.h:304
snd_seq_delete_port = _lib.snd_seq_delete_port
snd_seq_delete_port.restype = c_int
snd_seq_delete_port.argtypes = [POINTER(snd_seq_t), c_int]

# /usr/include/alsa/seq.h:305
snd_seq_get_port_info = _lib.snd_seq_get_port_info
snd_seq_get_port_info.restype = c_int
snd_seq_get_port_info.argtypes = [POINTER(snd_seq_t), c_int, POINTER(snd_seq_port_info_t)]

# /usr/include/alsa/seq.h:306
snd_seq_get_any_port_info = _lib.snd_seq_get_any_port_info
snd_seq_get_any_port_info.restype = c_int
snd_seq_get_any_port_info.argtypes = [POINTER(snd_seq_t), c_int, c_int, POINTER(snd_seq_port_info_t)]

# /usr/include/alsa/seq.h:307
snd_seq_set_port_info = _lib.snd_seq_set_port_info
snd_seq_set_port_info.restype = c_int
snd_seq_set_port_info.argtypes = [POINTER(snd_seq_t), c_int, POINTER(snd_seq_port_info_t)]

# /usr/include/alsa/seq.h:308
snd_seq_query_next_port = _lib.snd_seq_query_next_port
snd_seq_query_next_port.restype = c_int
snd_seq_query_next_port.argtypes = [POINTER(snd_seq_t), POINTER(snd_seq_port_info_t)]

class struct__snd_seq_port_subscribe(Structure):
    __slots__ = [
    ]
struct__snd_seq_port_subscribe._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_seq_port_subscribe(Structure):
    __slots__ = [
    ]
struct__snd_seq_port_subscribe._fields_ = [
    ('_opaque_struct', c_int)
]

snd_seq_port_subscribe_t = struct__snd_seq_port_subscribe 	# /usr/include/alsa/seq.h:321
# /usr/include/alsa/seq.h:323
snd_seq_port_subscribe_sizeof = _lib.snd_seq_port_subscribe_sizeof
snd_seq_port_subscribe_sizeof.restype = c_size_t
snd_seq_port_subscribe_sizeof.argtypes = []

# /usr/include/alsa/seq.h:327
snd_seq_port_subscribe_malloc = _lib.snd_seq_port_subscribe_malloc
snd_seq_port_subscribe_malloc.restype = c_int
snd_seq_port_subscribe_malloc.argtypes = [POINTER(POINTER(snd_seq_port_subscribe_t))]

# /usr/include/alsa/seq.h:328
snd_seq_port_subscribe_free = _lib.snd_seq_port_subscribe_free
snd_seq_port_subscribe_free.restype = None
snd_seq_port_subscribe_free.argtypes = [POINTER(snd_seq_port_subscribe_t)]

# /usr/include/alsa/seq.h:329
snd_seq_port_subscribe_copy = _lib.snd_seq_port_subscribe_copy
snd_seq_port_subscribe_copy.restype = None
snd_seq_port_subscribe_copy.argtypes = [POINTER(snd_seq_port_subscribe_t), POINTER(snd_seq_port_subscribe_t)]

# /usr/include/alsa/seq.h:331
snd_seq_port_subscribe_get_sender = _lib.snd_seq_port_subscribe_get_sender
snd_seq_port_subscribe_get_sender.restype = POINTER(snd_seq_addr_t)
snd_seq_port_subscribe_get_sender.argtypes = [POINTER(snd_seq_port_subscribe_t)]

# /usr/include/alsa/seq.h:332
snd_seq_port_subscribe_get_dest = _lib.snd_seq_port_subscribe_get_dest
snd_seq_port_subscribe_get_dest.restype = POINTER(snd_seq_addr_t)
snd_seq_port_subscribe_get_dest.argtypes = [POINTER(snd_seq_port_subscribe_t)]

# /usr/include/alsa/seq.h:333
snd_seq_port_subscribe_get_queue = _lib.snd_seq_port_subscribe_get_queue
snd_seq_port_subscribe_get_queue.restype = c_int
snd_seq_port_subscribe_get_queue.argtypes = [POINTER(snd_seq_port_subscribe_t)]

# /usr/include/alsa/seq.h:334
snd_seq_port_subscribe_get_exclusive = _lib.snd_seq_port_subscribe_get_exclusive
snd_seq_port_subscribe_get_exclusive.restype = c_int
snd_seq_port_subscribe_get_exclusive.argtypes = [POINTER(snd_seq_port_subscribe_t)]

# /usr/include/alsa/seq.h:335
snd_seq_port_subscribe_get_time_update = _lib.snd_seq_port_subscribe_get_time_update
snd_seq_port_subscribe_get_time_update.restype = c_int
snd_seq_port_subscribe_get_time_update.argtypes = [POINTER(snd_seq_port_subscribe_t)]

# /usr/include/alsa/seq.h:336
snd_seq_port_subscribe_get_time_real = _lib.snd_seq_port_subscribe_get_time_real
snd_seq_port_subscribe_get_time_real.restype = c_int
snd_seq_port_subscribe_get_time_real.argtypes = [POINTER(snd_seq_port_subscribe_t)]

# /usr/include/alsa/seq.h:338
snd_seq_port_subscribe_set_sender = _lib.snd_seq_port_subscribe_set_sender
snd_seq_port_subscribe_set_sender.restype = None
snd_seq_port_subscribe_set_sender.argtypes = [POINTER(snd_seq_port_subscribe_t), POINTER(snd_seq_addr_t)]

# /usr/include/alsa/seq.h:339
snd_seq_port_subscribe_set_dest = _lib.snd_seq_port_subscribe_set_dest
snd_seq_port_subscribe_set_dest.restype = None
snd_seq_port_subscribe_set_dest.argtypes = [POINTER(snd_seq_port_subscribe_t), POINTER(snd_seq_addr_t)]

# /usr/include/alsa/seq.h:340
snd_seq_port_subscribe_set_queue = _lib.snd_seq_port_subscribe_set_queue
snd_seq_port_subscribe_set_queue.restype = None
snd_seq_port_subscribe_set_queue.argtypes = [POINTER(snd_seq_port_subscribe_t), c_int]

# /usr/include/alsa/seq.h:341
snd_seq_port_subscribe_set_exclusive = _lib.snd_seq_port_subscribe_set_exclusive
snd_seq_port_subscribe_set_exclusive.restype = None
snd_seq_port_subscribe_set_exclusive.argtypes = [POINTER(snd_seq_port_subscribe_t), c_int]

# /usr/include/alsa/seq.h:342
snd_seq_port_subscribe_set_time_update = _lib.snd_seq_port_subscribe_set_time_update
snd_seq_port_subscribe_set_time_update.restype = None
snd_seq_port_subscribe_set_time_update.argtypes = [POINTER(snd_seq_port_subscribe_t), c_int]

# /usr/include/alsa/seq.h:343
snd_seq_port_subscribe_set_time_real = _lib.snd_seq_port_subscribe_set_time_real
snd_seq_port_subscribe_set_time_real.restype = None
snd_seq_port_subscribe_set_time_real.argtypes = [POINTER(snd_seq_port_subscribe_t), c_int]

# /usr/include/alsa/seq.h:345
snd_seq_get_port_subscription = _lib.snd_seq_get_port_subscription
snd_seq_get_port_subscription.restype = c_int
snd_seq_get_port_subscription.argtypes = [POINTER(snd_seq_t), POINTER(snd_seq_port_subscribe_t)]

# /usr/include/alsa/seq.h:346
snd_seq_subscribe_port = _lib.snd_seq_subscribe_port
snd_seq_subscribe_port.restype = c_int
snd_seq_subscribe_port.argtypes = [POINTER(snd_seq_t), POINTER(snd_seq_port_subscribe_t)]

# /usr/include/alsa/seq.h:347
snd_seq_unsubscribe_port = _lib.snd_seq_unsubscribe_port
snd_seq_unsubscribe_port.restype = c_int
snd_seq_unsubscribe_port.argtypes = [POINTER(snd_seq_t), POINTER(snd_seq_port_subscribe_t)]

class struct__snd_seq_query_subscribe(Structure):
    __slots__ = [
    ]
struct__snd_seq_query_subscribe._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_seq_query_subscribe(Structure):
    __slots__ = [
    ]
struct__snd_seq_query_subscribe._fields_ = [
    ('_opaque_struct', c_int)
]

snd_seq_query_subscribe_t = struct__snd_seq_query_subscribe 	# /usr/include/alsa/seq.h:353
enum_anon_30 = c_int
SND_SEQ_QUERY_SUBS_READ = 1
SND_SEQ_QUERY_SUBS_WRITE = 2
snd_seq_query_subs_type_t = enum_anon_30 	# /usr/include/alsa/seq.h:359
# /usr/include/alsa/seq.h:361
snd_seq_query_subscribe_sizeof = _lib.snd_seq_query_subscribe_sizeof
snd_seq_query_subscribe_sizeof.restype = c_size_t
snd_seq_query_subscribe_sizeof.argtypes = []

# /usr/include/alsa/seq.h:365
snd_seq_query_subscribe_malloc = _lib.snd_seq_query_subscribe_malloc
snd_seq_query_subscribe_malloc.restype = c_int
snd_seq_query_subscribe_malloc.argtypes = [POINTER(POINTER(snd_seq_query_subscribe_t))]

# /usr/include/alsa/seq.h:366
snd_seq_query_subscribe_free = _lib.snd_seq_query_subscribe_free
snd_seq_query_subscribe_free.restype = None
snd_seq_query_subscribe_free.argtypes = [POINTER(snd_seq_query_subscribe_t)]

# /usr/include/alsa/seq.h:367
snd_seq_query_subscribe_copy = _lib.snd_seq_query_subscribe_copy
snd_seq_query_subscribe_copy.restype = None
snd_seq_query_subscribe_copy.argtypes = [POINTER(snd_seq_query_subscribe_t), POINTER(snd_seq_query_subscribe_t)]

# /usr/include/alsa/seq.h:369
snd_seq_query_subscribe_get_client = _lib.snd_seq_query_subscribe_get_client
snd_seq_query_subscribe_get_client.restype = c_int
snd_seq_query_subscribe_get_client.argtypes = [POINTER(snd_seq_query_subscribe_t)]

# /usr/include/alsa/seq.h:370
snd_seq_query_subscribe_get_port = _lib.snd_seq_query_subscribe_get_port
snd_seq_query_subscribe_get_port.restype = c_int
snd_seq_query_subscribe_get_port.argtypes = [POINTER(snd_seq_query_subscribe_t)]

# /usr/include/alsa/seq.h:371
snd_seq_query_subscribe_get_root = _lib.snd_seq_query_subscribe_get_root
snd_seq_query_subscribe_get_root.restype = POINTER(snd_seq_addr_t)
snd_seq_query_subscribe_get_root.argtypes = [POINTER(snd_seq_query_subscribe_t)]

# /usr/include/alsa/seq.h:372
snd_seq_query_subscribe_get_type = _lib.snd_seq_query_subscribe_get_type
snd_seq_query_subscribe_get_type.restype = snd_seq_query_subs_type_t
snd_seq_query_subscribe_get_type.argtypes = [POINTER(snd_seq_query_subscribe_t)]

# /usr/include/alsa/seq.h:373
snd_seq_query_subscribe_get_index = _lib.snd_seq_query_subscribe_get_index
snd_seq_query_subscribe_get_index.restype = c_int
snd_seq_query_subscribe_get_index.argtypes = [POINTER(snd_seq_query_subscribe_t)]

# /usr/include/alsa/seq.h:374
snd_seq_query_subscribe_get_num_subs = _lib.snd_seq_query_subscribe_get_num_subs
snd_seq_query_subscribe_get_num_subs.restype = c_int
snd_seq_query_subscribe_get_num_subs.argtypes = [POINTER(snd_seq_query_subscribe_t)]

# /usr/include/alsa/seq.h:375
snd_seq_query_subscribe_get_addr = _lib.snd_seq_query_subscribe_get_addr
snd_seq_query_subscribe_get_addr.restype = POINTER(snd_seq_addr_t)
snd_seq_query_subscribe_get_addr.argtypes = [POINTER(snd_seq_query_subscribe_t)]

# /usr/include/alsa/seq.h:376
snd_seq_query_subscribe_get_queue = _lib.snd_seq_query_subscribe_get_queue
snd_seq_query_subscribe_get_queue.restype = c_int
snd_seq_query_subscribe_get_queue.argtypes = [POINTER(snd_seq_query_subscribe_t)]

# /usr/include/alsa/seq.h:377
snd_seq_query_subscribe_get_exclusive = _lib.snd_seq_query_subscribe_get_exclusive
snd_seq_query_subscribe_get_exclusive.restype = c_int
snd_seq_query_subscribe_get_exclusive.argtypes = [POINTER(snd_seq_query_subscribe_t)]

# /usr/include/alsa/seq.h:378
snd_seq_query_subscribe_get_time_update = _lib.snd_seq_query_subscribe_get_time_update
snd_seq_query_subscribe_get_time_update.restype = c_int
snd_seq_query_subscribe_get_time_update.argtypes = [POINTER(snd_seq_query_subscribe_t)]

# /usr/include/alsa/seq.h:379
snd_seq_query_subscribe_get_time_real = _lib.snd_seq_query_subscribe_get_time_real
snd_seq_query_subscribe_get_time_real.restype = c_int
snd_seq_query_subscribe_get_time_real.argtypes = [POINTER(snd_seq_query_subscribe_t)]

# /usr/include/alsa/seq.h:381
snd_seq_query_subscribe_set_client = _lib.snd_seq_query_subscribe_set_client
snd_seq_query_subscribe_set_client.restype = None
snd_seq_query_subscribe_set_client.argtypes = [POINTER(snd_seq_query_subscribe_t), c_int]

# /usr/include/alsa/seq.h:382
snd_seq_query_subscribe_set_port = _lib.snd_seq_query_subscribe_set_port
snd_seq_query_subscribe_set_port.restype = None
snd_seq_query_subscribe_set_port.argtypes = [POINTER(snd_seq_query_subscribe_t), c_int]

# /usr/include/alsa/seq.h:383
snd_seq_query_subscribe_set_root = _lib.snd_seq_query_subscribe_set_root
snd_seq_query_subscribe_set_root.restype = None
snd_seq_query_subscribe_set_root.argtypes = [POINTER(snd_seq_query_subscribe_t), POINTER(snd_seq_addr_t)]

# /usr/include/alsa/seq.h:384
snd_seq_query_subscribe_set_type = _lib.snd_seq_query_subscribe_set_type
snd_seq_query_subscribe_set_type.restype = None
snd_seq_query_subscribe_set_type.argtypes = [POINTER(snd_seq_query_subscribe_t), snd_seq_query_subs_type_t]

# /usr/include/alsa/seq.h:385
snd_seq_query_subscribe_set_index = _lib.snd_seq_query_subscribe_set_index
snd_seq_query_subscribe_set_index.restype = None
snd_seq_query_subscribe_set_index.argtypes = [POINTER(snd_seq_query_subscribe_t), c_int]

# /usr/include/alsa/seq.h:387
snd_seq_query_port_subscribers = _lib.snd_seq_query_port_subscribers
snd_seq_query_port_subscribers.restype = c_int
snd_seq_query_port_subscribers.argtypes = [POINTER(snd_seq_t), POINTER(snd_seq_query_subscribe_t)]

class struct__snd_seq_queue_info(Structure):
    __slots__ = [
    ]
struct__snd_seq_queue_info._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_seq_queue_info(Structure):
    __slots__ = [
    ]
struct__snd_seq_queue_info._fields_ = [
    ('_opaque_struct', c_int)
]

snd_seq_queue_info_t = struct__snd_seq_queue_info 	# /usr/include/alsa/seq.h:400
class struct__snd_seq_queue_status(Structure):
    __slots__ = [
    ]
struct__snd_seq_queue_status._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_seq_queue_status(Structure):
    __slots__ = [
    ]
struct__snd_seq_queue_status._fields_ = [
    ('_opaque_struct', c_int)
]

snd_seq_queue_status_t = struct__snd_seq_queue_status 	# /usr/include/alsa/seq.h:402
class struct__snd_seq_queue_tempo(Structure):
    __slots__ = [
    ]
struct__snd_seq_queue_tempo._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_seq_queue_tempo(Structure):
    __slots__ = [
    ]
struct__snd_seq_queue_tempo._fields_ = [
    ('_opaque_struct', c_int)
]

snd_seq_queue_tempo_t = struct__snd_seq_queue_tempo 	# /usr/include/alsa/seq.h:404
class struct__snd_seq_queue_timer(Structure):
    __slots__ = [
    ]
struct__snd_seq_queue_timer._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_seq_queue_timer(Structure):
    __slots__ = [
    ]
struct__snd_seq_queue_timer._fields_ = [
    ('_opaque_struct', c_int)
]

snd_seq_queue_timer_t = struct__snd_seq_queue_timer 	# /usr/include/alsa/seq.h:406
SND_SEQ_QUEUE_DIRECT = 253 	# /usr/include/alsa/seq.h:409
# /usr/include/alsa/seq.h:411
snd_seq_queue_info_sizeof = _lib.snd_seq_queue_info_sizeof
snd_seq_queue_info_sizeof.restype = c_size_t
snd_seq_queue_info_sizeof.argtypes = []

# /usr/include/alsa/seq.h:415
snd_seq_queue_info_malloc = _lib.snd_seq_queue_info_malloc
snd_seq_queue_info_malloc.restype = c_int
snd_seq_queue_info_malloc.argtypes = [POINTER(POINTER(snd_seq_queue_info_t))]

# /usr/include/alsa/seq.h:416
snd_seq_queue_info_free = _lib.snd_seq_queue_info_free
snd_seq_queue_info_free.restype = None
snd_seq_queue_info_free.argtypes = [POINTER(snd_seq_queue_info_t)]

# /usr/include/alsa/seq.h:417
snd_seq_queue_info_copy = _lib.snd_seq_queue_info_copy
snd_seq_queue_info_copy.restype = None
snd_seq_queue_info_copy.argtypes = [POINTER(snd_seq_queue_info_t), POINTER(snd_seq_queue_info_t)]

# /usr/include/alsa/seq.h:419
snd_seq_queue_info_get_queue = _lib.snd_seq_queue_info_get_queue
snd_seq_queue_info_get_queue.restype = c_int
snd_seq_queue_info_get_queue.argtypes = [POINTER(snd_seq_queue_info_t)]

# /usr/include/alsa/seq.h:420
snd_seq_queue_info_get_name = _lib.snd_seq_queue_info_get_name
snd_seq_queue_info_get_name.restype = c_char_p
snd_seq_queue_info_get_name.argtypes = [POINTER(snd_seq_queue_info_t)]

# /usr/include/alsa/seq.h:421
snd_seq_queue_info_get_owner = _lib.snd_seq_queue_info_get_owner
snd_seq_queue_info_get_owner.restype = c_int
snd_seq_queue_info_get_owner.argtypes = [POINTER(snd_seq_queue_info_t)]

# /usr/include/alsa/seq.h:422
snd_seq_queue_info_get_locked = _lib.snd_seq_queue_info_get_locked
snd_seq_queue_info_get_locked.restype = c_int
snd_seq_queue_info_get_locked.argtypes = [POINTER(snd_seq_queue_info_t)]

# /usr/include/alsa/seq.h:423
snd_seq_queue_info_get_flags = _lib.snd_seq_queue_info_get_flags
snd_seq_queue_info_get_flags.restype = c_uint
snd_seq_queue_info_get_flags.argtypes = [POINTER(snd_seq_queue_info_t)]

# /usr/include/alsa/seq.h:425
snd_seq_queue_info_set_name = _lib.snd_seq_queue_info_set_name
snd_seq_queue_info_set_name.restype = None
snd_seq_queue_info_set_name.argtypes = [POINTER(snd_seq_queue_info_t), c_char_p]

# /usr/include/alsa/seq.h:426
snd_seq_queue_info_set_owner = _lib.snd_seq_queue_info_set_owner
snd_seq_queue_info_set_owner.restype = None
snd_seq_queue_info_set_owner.argtypes = [POINTER(snd_seq_queue_info_t), c_int]

# /usr/include/alsa/seq.h:427
snd_seq_queue_info_set_locked = _lib.snd_seq_queue_info_set_locked
snd_seq_queue_info_set_locked.restype = None
snd_seq_queue_info_set_locked.argtypes = [POINTER(snd_seq_queue_info_t), c_int]

# /usr/include/alsa/seq.h:428
snd_seq_queue_info_set_flags = _lib.snd_seq_queue_info_set_flags
snd_seq_queue_info_set_flags.restype = None
snd_seq_queue_info_set_flags.argtypes = [POINTER(snd_seq_queue_info_t), c_uint]

# /usr/include/alsa/seq.h:430
snd_seq_create_queue = _lib.snd_seq_create_queue
snd_seq_create_queue.restype = c_int
snd_seq_create_queue.argtypes = [POINTER(snd_seq_t), POINTER(snd_seq_queue_info_t)]

# /usr/include/alsa/seq.h:431
snd_seq_alloc_named_queue = _lib.snd_seq_alloc_named_queue
snd_seq_alloc_named_queue.restype = c_int
snd_seq_alloc_named_queue.argtypes = [POINTER(snd_seq_t), c_char_p]

# /usr/include/alsa/seq.h:432
snd_seq_alloc_queue = _lib.snd_seq_alloc_queue
snd_seq_alloc_queue.restype = c_int
snd_seq_alloc_queue.argtypes = [POINTER(snd_seq_t)]

# /usr/include/alsa/seq.h:433
snd_seq_free_queue = _lib.snd_seq_free_queue
snd_seq_free_queue.restype = c_int
snd_seq_free_queue.argtypes = [POINTER(snd_seq_t), c_int]

# /usr/include/alsa/seq.h:434
snd_seq_get_queue_info = _lib.snd_seq_get_queue_info
snd_seq_get_queue_info.restype = c_int
snd_seq_get_queue_info.argtypes = [POINTER(snd_seq_t), c_int, POINTER(snd_seq_queue_info_t)]

# /usr/include/alsa/seq.h:435
snd_seq_set_queue_info = _lib.snd_seq_set_queue_info
snd_seq_set_queue_info.restype = c_int
snd_seq_set_queue_info.argtypes = [POINTER(snd_seq_t), c_int, POINTER(snd_seq_queue_info_t)]

# /usr/include/alsa/seq.h:436
snd_seq_query_named_queue = _lib.snd_seq_query_named_queue
snd_seq_query_named_queue.restype = c_int
snd_seq_query_named_queue.argtypes = [POINTER(snd_seq_t), c_char_p]

# /usr/include/alsa/seq.h:438
snd_seq_get_queue_usage = _lib.snd_seq_get_queue_usage
snd_seq_get_queue_usage.restype = c_int
snd_seq_get_queue_usage.argtypes = [POINTER(snd_seq_t), c_int]

# /usr/include/alsa/seq.h:439
snd_seq_set_queue_usage = _lib.snd_seq_set_queue_usage
snd_seq_set_queue_usage.restype = c_int
snd_seq_set_queue_usage.argtypes = [POINTER(snd_seq_t), c_int, c_int]

# /usr/include/alsa/seq.h:443
snd_seq_queue_status_sizeof = _lib.snd_seq_queue_status_sizeof
snd_seq_queue_status_sizeof.restype = c_size_t
snd_seq_queue_status_sizeof.argtypes = []

# /usr/include/alsa/seq.h:447
snd_seq_queue_status_malloc = _lib.snd_seq_queue_status_malloc
snd_seq_queue_status_malloc.restype = c_int
snd_seq_queue_status_malloc.argtypes = [POINTER(POINTER(snd_seq_queue_status_t))]

# /usr/include/alsa/seq.h:448
snd_seq_queue_status_free = _lib.snd_seq_queue_status_free
snd_seq_queue_status_free.restype = None
snd_seq_queue_status_free.argtypes = [POINTER(snd_seq_queue_status_t)]

# /usr/include/alsa/seq.h:449
snd_seq_queue_status_copy = _lib.snd_seq_queue_status_copy
snd_seq_queue_status_copy.restype = None
snd_seq_queue_status_copy.argtypes = [POINTER(snd_seq_queue_status_t), POINTER(snd_seq_queue_status_t)]

# /usr/include/alsa/seq.h:451
snd_seq_queue_status_get_queue = _lib.snd_seq_queue_status_get_queue
snd_seq_queue_status_get_queue.restype = c_int
snd_seq_queue_status_get_queue.argtypes = [POINTER(snd_seq_queue_status_t)]

# /usr/include/alsa/seq.h:452
snd_seq_queue_status_get_events = _lib.snd_seq_queue_status_get_events
snd_seq_queue_status_get_events.restype = c_int
snd_seq_queue_status_get_events.argtypes = [POINTER(snd_seq_queue_status_t)]

# /usr/include/alsa/seq.h:453
snd_seq_queue_status_get_tick_time = _lib.snd_seq_queue_status_get_tick_time
snd_seq_queue_status_get_tick_time.restype = snd_seq_tick_time_t
snd_seq_queue_status_get_tick_time.argtypes = [POINTER(snd_seq_queue_status_t)]

# /usr/include/alsa/seq.h:454
snd_seq_queue_status_get_real_time = _lib.snd_seq_queue_status_get_real_time
snd_seq_queue_status_get_real_time.restype = POINTER(snd_seq_real_time_t)
snd_seq_queue_status_get_real_time.argtypes = [POINTER(snd_seq_queue_status_t)]

# /usr/include/alsa/seq.h:455
snd_seq_queue_status_get_status = _lib.snd_seq_queue_status_get_status
snd_seq_queue_status_get_status.restype = c_uint
snd_seq_queue_status_get_status.argtypes = [POINTER(snd_seq_queue_status_t)]

# /usr/include/alsa/seq.h:457
snd_seq_get_queue_status = _lib.snd_seq_get_queue_status
snd_seq_get_queue_status.restype = c_int
snd_seq_get_queue_status.argtypes = [POINTER(snd_seq_t), c_int, POINTER(snd_seq_queue_status_t)]

# /usr/include/alsa/seq.h:461
snd_seq_queue_tempo_sizeof = _lib.snd_seq_queue_tempo_sizeof
snd_seq_queue_tempo_sizeof.restype = c_size_t
snd_seq_queue_tempo_sizeof.argtypes = []

# /usr/include/alsa/seq.h:465
snd_seq_queue_tempo_malloc = _lib.snd_seq_queue_tempo_malloc
snd_seq_queue_tempo_malloc.restype = c_int
snd_seq_queue_tempo_malloc.argtypes = [POINTER(POINTER(snd_seq_queue_tempo_t))]

# /usr/include/alsa/seq.h:466
snd_seq_queue_tempo_free = _lib.snd_seq_queue_tempo_free
snd_seq_queue_tempo_free.restype = None
snd_seq_queue_tempo_free.argtypes = [POINTER(snd_seq_queue_tempo_t)]

# /usr/include/alsa/seq.h:467
snd_seq_queue_tempo_copy = _lib.snd_seq_queue_tempo_copy
snd_seq_queue_tempo_copy.restype = None
snd_seq_queue_tempo_copy.argtypes = [POINTER(snd_seq_queue_tempo_t), POINTER(snd_seq_queue_tempo_t)]

# /usr/include/alsa/seq.h:469
snd_seq_queue_tempo_get_queue = _lib.snd_seq_queue_tempo_get_queue
snd_seq_queue_tempo_get_queue.restype = c_int
snd_seq_queue_tempo_get_queue.argtypes = [POINTER(snd_seq_queue_tempo_t)]

# /usr/include/alsa/seq.h:470
snd_seq_queue_tempo_get_tempo = _lib.snd_seq_queue_tempo_get_tempo
snd_seq_queue_tempo_get_tempo.restype = c_uint
snd_seq_queue_tempo_get_tempo.argtypes = [POINTER(snd_seq_queue_tempo_t)]

# /usr/include/alsa/seq.h:471
snd_seq_queue_tempo_get_ppq = _lib.snd_seq_queue_tempo_get_ppq
snd_seq_queue_tempo_get_ppq.restype = c_int
snd_seq_queue_tempo_get_ppq.argtypes = [POINTER(snd_seq_queue_tempo_t)]

# /usr/include/alsa/seq.h:472
snd_seq_queue_tempo_get_skew = _lib.snd_seq_queue_tempo_get_skew
snd_seq_queue_tempo_get_skew.restype = c_uint
snd_seq_queue_tempo_get_skew.argtypes = [POINTER(snd_seq_queue_tempo_t)]

# /usr/include/alsa/seq.h:473
snd_seq_queue_tempo_get_skew_base = _lib.snd_seq_queue_tempo_get_skew_base
snd_seq_queue_tempo_get_skew_base.restype = c_uint
snd_seq_queue_tempo_get_skew_base.argtypes = [POINTER(snd_seq_queue_tempo_t)]

# /usr/include/alsa/seq.h:474
snd_seq_queue_tempo_set_tempo = _lib.snd_seq_queue_tempo_set_tempo
snd_seq_queue_tempo_set_tempo.restype = None
snd_seq_queue_tempo_set_tempo.argtypes = [POINTER(snd_seq_queue_tempo_t), c_uint]

# /usr/include/alsa/seq.h:475
snd_seq_queue_tempo_set_ppq = _lib.snd_seq_queue_tempo_set_ppq
snd_seq_queue_tempo_set_ppq.restype = None
snd_seq_queue_tempo_set_ppq.argtypes = [POINTER(snd_seq_queue_tempo_t), c_int]

# /usr/include/alsa/seq.h:476
snd_seq_queue_tempo_set_skew = _lib.snd_seq_queue_tempo_set_skew
snd_seq_queue_tempo_set_skew.restype = None
snd_seq_queue_tempo_set_skew.argtypes = [POINTER(snd_seq_queue_tempo_t), c_uint]

# /usr/include/alsa/seq.h:477
snd_seq_queue_tempo_set_skew_base = _lib.snd_seq_queue_tempo_set_skew_base
snd_seq_queue_tempo_set_skew_base.restype = None
snd_seq_queue_tempo_set_skew_base.argtypes = [POINTER(snd_seq_queue_tempo_t), c_uint]

# /usr/include/alsa/seq.h:479
snd_seq_get_queue_tempo = _lib.snd_seq_get_queue_tempo
snd_seq_get_queue_tempo.restype = c_int
snd_seq_get_queue_tempo.argtypes = [POINTER(snd_seq_t), c_int, POINTER(snd_seq_queue_tempo_t)]

# /usr/include/alsa/seq.h:480
snd_seq_set_queue_tempo = _lib.snd_seq_set_queue_tempo
snd_seq_set_queue_tempo.restype = c_int
snd_seq_set_queue_tempo.argtypes = [POINTER(snd_seq_t), c_int, POINTER(snd_seq_queue_tempo_t)]

enum_anon_31 = c_int
SND_SEQ_TIMER_ALSA = 0
SND_SEQ_TIMER_MIDI_CLOCK = 1
SND_SEQ_TIMER_MIDI_TICK = 2
snd_seq_queue_timer_type_t = enum_anon_31 	# /usr/include/alsa/seq.h:490
# /usr/include/alsa/seq.h:492
snd_seq_queue_timer_sizeof = _lib.snd_seq_queue_timer_sizeof
snd_seq_queue_timer_sizeof.restype = c_size_t
snd_seq_queue_timer_sizeof.argtypes = []

# /usr/include/alsa/seq.h:496
snd_seq_queue_timer_malloc = _lib.snd_seq_queue_timer_malloc
snd_seq_queue_timer_malloc.restype = c_int
snd_seq_queue_timer_malloc.argtypes = [POINTER(POINTER(snd_seq_queue_timer_t))]

# /usr/include/alsa/seq.h:497
snd_seq_queue_timer_free = _lib.snd_seq_queue_timer_free
snd_seq_queue_timer_free.restype = None
snd_seq_queue_timer_free.argtypes = [POINTER(snd_seq_queue_timer_t)]

# /usr/include/alsa/seq.h:498
snd_seq_queue_timer_copy = _lib.snd_seq_queue_timer_copy
snd_seq_queue_timer_copy.restype = None
snd_seq_queue_timer_copy.argtypes = [POINTER(snd_seq_queue_timer_t), POINTER(snd_seq_queue_timer_t)]

# /usr/include/alsa/seq.h:500
snd_seq_queue_timer_get_queue = _lib.snd_seq_queue_timer_get_queue
snd_seq_queue_timer_get_queue.restype = c_int
snd_seq_queue_timer_get_queue.argtypes = [POINTER(snd_seq_queue_timer_t)]

# /usr/include/alsa/seq.h:501
snd_seq_queue_timer_get_type = _lib.snd_seq_queue_timer_get_type
snd_seq_queue_timer_get_type.restype = snd_seq_queue_timer_type_t
snd_seq_queue_timer_get_type.argtypes = [POINTER(snd_seq_queue_timer_t)]

# /usr/include/alsa/seq.h:502
snd_seq_queue_timer_get_id = _lib.snd_seq_queue_timer_get_id
snd_seq_queue_timer_get_id.restype = POINTER(snd_timer_id_t)
snd_seq_queue_timer_get_id.argtypes = [POINTER(snd_seq_queue_timer_t)]

# /usr/include/alsa/seq.h:503
snd_seq_queue_timer_get_resolution = _lib.snd_seq_queue_timer_get_resolution
snd_seq_queue_timer_get_resolution.restype = c_uint
snd_seq_queue_timer_get_resolution.argtypes = [POINTER(snd_seq_queue_timer_t)]

# /usr/include/alsa/seq.h:505
snd_seq_queue_timer_set_type = _lib.snd_seq_queue_timer_set_type
snd_seq_queue_timer_set_type.restype = None
snd_seq_queue_timer_set_type.argtypes = [POINTER(snd_seq_queue_timer_t), snd_seq_queue_timer_type_t]

# /usr/include/alsa/seq.h:506
snd_seq_queue_timer_set_id = _lib.snd_seq_queue_timer_set_id
snd_seq_queue_timer_set_id.restype = None
snd_seq_queue_timer_set_id.argtypes = [POINTER(snd_seq_queue_timer_t), POINTER(snd_timer_id_t)]

# /usr/include/alsa/seq.h:507
snd_seq_queue_timer_set_resolution = _lib.snd_seq_queue_timer_set_resolution
snd_seq_queue_timer_set_resolution.restype = None
snd_seq_queue_timer_set_resolution.argtypes = [POINTER(snd_seq_queue_timer_t), c_uint]

# /usr/include/alsa/seq.h:509
snd_seq_get_queue_timer = _lib.snd_seq_get_queue_timer
snd_seq_get_queue_timer.restype = c_int
snd_seq_get_queue_timer.argtypes = [POINTER(snd_seq_t), c_int, POINTER(snd_seq_queue_timer_t)]

# /usr/include/alsa/seq.h:510
snd_seq_set_queue_timer = _lib.snd_seq_set_queue_timer
snd_seq_set_queue_timer.restype = c_int
snd_seq_set_queue_timer.argtypes = [POINTER(snd_seq_t), c_int, POINTER(snd_seq_queue_timer_t)]

# /usr/include/alsa/seq.h:521
snd_seq_free_event = _lib.snd_seq_free_event
snd_seq_free_event.restype = c_int
snd_seq_free_event.argtypes = [POINTER(snd_seq_event_t)]

# /usr/include/alsa/seq.h:522
snd_seq_event_length = _lib.snd_seq_event_length
snd_seq_event_length.restype = ssize_t
snd_seq_event_length.argtypes = [POINTER(snd_seq_event_t)]

# /usr/include/alsa/seq.h:523
snd_seq_event_output = _lib.snd_seq_event_output
snd_seq_event_output.restype = c_int
snd_seq_event_output.argtypes = [POINTER(snd_seq_t), POINTER(snd_seq_event_t)]

# /usr/include/alsa/seq.h:524
snd_seq_event_output_buffer = _lib.snd_seq_event_output_buffer
snd_seq_event_output_buffer.restype = c_int
snd_seq_event_output_buffer.argtypes = [POINTER(snd_seq_t), POINTER(snd_seq_event_t)]

# /usr/include/alsa/seq.h:525
snd_seq_event_output_direct = _lib.snd_seq_event_output_direct
snd_seq_event_output_direct.restype = c_int
snd_seq_event_output_direct.argtypes = [POINTER(snd_seq_t), POINTER(snd_seq_event_t)]

# /usr/include/alsa/seq.h:526
snd_seq_event_input = _lib.snd_seq_event_input
snd_seq_event_input.restype = c_int
snd_seq_event_input.argtypes = [POINTER(snd_seq_t), POINTER(POINTER(snd_seq_event_t))]

# /usr/include/alsa/seq.h:527
snd_seq_event_input_pending = _lib.snd_seq_event_input_pending
snd_seq_event_input_pending.restype = c_int
snd_seq_event_input_pending.argtypes = [POINTER(snd_seq_t), c_int]

# /usr/include/alsa/seq.h:528
snd_seq_drain_output = _lib.snd_seq_drain_output
snd_seq_drain_output.restype = c_int
snd_seq_drain_output.argtypes = [POINTER(snd_seq_t)]

# /usr/include/alsa/seq.h:529
snd_seq_event_output_pending = _lib.snd_seq_event_output_pending
snd_seq_event_output_pending.restype = c_int
snd_seq_event_output_pending.argtypes = [POINTER(snd_seq_t)]

# /usr/include/alsa/seq.h:530
snd_seq_extract_output = _lib.snd_seq_extract_output
snd_seq_extract_output.restype = c_int
snd_seq_extract_output.argtypes = [POINTER(snd_seq_t), POINTER(POINTER(snd_seq_event_t))]

# /usr/include/alsa/seq.h:531
snd_seq_drop_output = _lib.snd_seq_drop_output
snd_seq_drop_output.restype = c_int
snd_seq_drop_output.argtypes = [POINTER(snd_seq_t)]

# /usr/include/alsa/seq.h:532
snd_seq_drop_output_buffer = _lib.snd_seq_drop_output_buffer
snd_seq_drop_output_buffer.restype = c_int
snd_seq_drop_output_buffer.argtypes = [POINTER(snd_seq_t)]

# /usr/include/alsa/seq.h:533
snd_seq_drop_input = _lib.snd_seq_drop_input
snd_seq_drop_input.restype = c_int
snd_seq_drop_input.argtypes = [POINTER(snd_seq_t)]

# /usr/include/alsa/seq.h:534
snd_seq_drop_input_buffer = _lib.snd_seq_drop_input_buffer
snd_seq_drop_input_buffer.restype = c_int
snd_seq_drop_input_buffer.argtypes = [POINTER(snd_seq_t)]

class struct__snd_seq_remove_events(Structure):
    __slots__ = [
    ]
struct__snd_seq_remove_events._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_seq_remove_events(Structure):
    __slots__ = [
    ]
struct__snd_seq_remove_events._fields_ = [
    ('_opaque_struct', c_int)
]

snd_seq_remove_events_t = struct__snd_seq_remove_events 	# /usr/include/alsa/seq.h:537
SND_SEQ_REMOVE_INPUT = 1 	# /usr/include/alsa/seq.h:540
SND_SEQ_REMOVE_OUTPUT = 2 	# /usr/include/alsa/seq.h:541
SND_SEQ_REMOVE_DEST = 4 	# /usr/include/alsa/seq.h:542
SND_SEQ_REMOVE_DEST_CHANNEL = 8 	# /usr/include/alsa/seq.h:543
SND_SEQ_REMOVE_TIME_BEFORE = 16 	# /usr/include/alsa/seq.h:544
SND_SEQ_REMOVE_TIME_AFTER = 32 	# /usr/include/alsa/seq.h:545
SND_SEQ_REMOVE_TIME_TICK = 64 	# /usr/include/alsa/seq.h:546
SND_SEQ_REMOVE_EVENT_TYPE = 128 	# /usr/include/alsa/seq.h:547
SND_SEQ_REMOVE_IGNORE_OFF = 256 	# /usr/include/alsa/seq.h:548
SND_SEQ_REMOVE_TAG_MATCH = 512 	# /usr/include/alsa/seq.h:549
# /usr/include/alsa/seq.h:551
snd_seq_remove_events_sizeof = _lib.snd_seq_remove_events_sizeof
snd_seq_remove_events_sizeof.restype = c_size_t
snd_seq_remove_events_sizeof.argtypes = []

# /usr/include/alsa/seq.h:555
snd_seq_remove_events_malloc = _lib.snd_seq_remove_events_malloc
snd_seq_remove_events_malloc.restype = c_int
snd_seq_remove_events_malloc.argtypes = [POINTER(POINTER(snd_seq_remove_events_t))]

# /usr/include/alsa/seq.h:556
snd_seq_remove_events_free = _lib.snd_seq_remove_events_free
snd_seq_remove_events_free.restype = None
snd_seq_remove_events_free.argtypes = [POINTER(snd_seq_remove_events_t)]

# /usr/include/alsa/seq.h:557
snd_seq_remove_events_copy = _lib.snd_seq_remove_events_copy
snd_seq_remove_events_copy.restype = None
snd_seq_remove_events_copy.argtypes = [POINTER(snd_seq_remove_events_t), POINTER(snd_seq_remove_events_t)]

# /usr/include/alsa/seq.h:559
snd_seq_remove_events_get_condition = _lib.snd_seq_remove_events_get_condition
snd_seq_remove_events_get_condition.restype = c_uint
snd_seq_remove_events_get_condition.argtypes = [POINTER(snd_seq_remove_events_t)]

# /usr/include/alsa/seq.h:560
snd_seq_remove_events_get_queue = _lib.snd_seq_remove_events_get_queue
snd_seq_remove_events_get_queue.restype = c_int
snd_seq_remove_events_get_queue.argtypes = [POINTER(snd_seq_remove_events_t)]

# /usr/include/alsa/seq.h:561
snd_seq_remove_events_get_time = _lib.snd_seq_remove_events_get_time
snd_seq_remove_events_get_time.restype = POINTER(snd_seq_timestamp_t)
snd_seq_remove_events_get_time.argtypes = [POINTER(snd_seq_remove_events_t)]

# /usr/include/alsa/seq.h:562
snd_seq_remove_events_get_dest = _lib.snd_seq_remove_events_get_dest
snd_seq_remove_events_get_dest.restype = POINTER(snd_seq_addr_t)
snd_seq_remove_events_get_dest.argtypes = [POINTER(snd_seq_remove_events_t)]

# /usr/include/alsa/seq.h:563
snd_seq_remove_events_get_channel = _lib.snd_seq_remove_events_get_channel
snd_seq_remove_events_get_channel.restype = c_int
snd_seq_remove_events_get_channel.argtypes = [POINTER(snd_seq_remove_events_t)]

# /usr/include/alsa/seq.h:564
snd_seq_remove_events_get_event_type = _lib.snd_seq_remove_events_get_event_type
snd_seq_remove_events_get_event_type.restype = c_int
snd_seq_remove_events_get_event_type.argtypes = [POINTER(snd_seq_remove_events_t)]

# /usr/include/alsa/seq.h:565
snd_seq_remove_events_get_tag = _lib.snd_seq_remove_events_get_tag
snd_seq_remove_events_get_tag.restype = c_int
snd_seq_remove_events_get_tag.argtypes = [POINTER(snd_seq_remove_events_t)]

# /usr/include/alsa/seq.h:567
snd_seq_remove_events_set_condition = _lib.snd_seq_remove_events_set_condition
snd_seq_remove_events_set_condition.restype = None
snd_seq_remove_events_set_condition.argtypes = [POINTER(snd_seq_remove_events_t), c_uint]

# /usr/include/alsa/seq.h:568
snd_seq_remove_events_set_queue = _lib.snd_seq_remove_events_set_queue
snd_seq_remove_events_set_queue.restype = None
snd_seq_remove_events_set_queue.argtypes = [POINTER(snd_seq_remove_events_t), c_int]

# /usr/include/alsa/seq.h:569
snd_seq_remove_events_set_time = _lib.snd_seq_remove_events_set_time
snd_seq_remove_events_set_time.restype = None
snd_seq_remove_events_set_time.argtypes = [POINTER(snd_seq_remove_events_t), POINTER(snd_seq_timestamp_t)]

# /usr/include/alsa/seq.h:570
snd_seq_remove_events_set_dest = _lib.snd_seq_remove_events_set_dest
snd_seq_remove_events_set_dest.restype = None
snd_seq_remove_events_set_dest.argtypes = [POINTER(snd_seq_remove_events_t), POINTER(snd_seq_addr_t)]

# /usr/include/alsa/seq.h:571
snd_seq_remove_events_set_channel = _lib.snd_seq_remove_events_set_channel
snd_seq_remove_events_set_channel.restype = None
snd_seq_remove_events_set_channel.argtypes = [POINTER(snd_seq_remove_events_t), c_int]

# /usr/include/alsa/seq.h:572
snd_seq_remove_events_set_event_type = _lib.snd_seq_remove_events_set_event_type
snd_seq_remove_events_set_event_type.restype = None
snd_seq_remove_events_set_event_type.argtypes = [POINTER(snd_seq_remove_events_t), c_int]

# /usr/include/alsa/seq.h:573
snd_seq_remove_events_set_tag = _lib.snd_seq_remove_events_set_tag
snd_seq_remove_events_set_tag.restype = None
snd_seq_remove_events_set_tag.argtypes = [POINTER(snd_seq_remove_events_t), c_int]

# /usr/include/alsa/seq.h:575
snd_seq_remove_events = _lib.snd_seq_remove_events
snd_seq_remove_events.restype = c_int
snd_seq_remove_events.argtypes = [POINTER(snd_seq_t), POINTER(snd_seq_remove_events_t)]

# /usr/include/alsa/seq.h:586
snd_seq_set_bit = _lib.snd_seq_set_bit
snd_seq_set_bit.restype = None
snd_seq_set_bit.argtypes = [c_int, POINTER(None)]

# /usr/include/alsa/seq.h:587
snd_seq_change_bit = _lib.snd_seq_change_bit
snd_seq_change_bit.restype = c_int
snd_seq_change_bit.argtypes = [c_int, POINTER(None)]

# /usr/include/alsa/seq.h:588
snd_seq_get_bit = _lib.snd_seq_get_bit
snd_seq_get_bit.restype = c_int
snd_seq_get_bit.argtypes = [c_int, POINTER(None)]

# /usr/include/alsa/seqmid.h:288
snd_seq_control_queue = _lib.snd_seq_control_queue
snd_seq_control_queue.restype = c_int
snd_seq_control_queue.argtypes = [POINTER(snd_seq_t), c_int, c_int, c_int, POINTER(snd_seq_event_t)]

# /usr/include/alsa/seqmid.h:328
snd_seq_create_simple_port = _lib.snd_seq_create_simple_port
snd_seq_create_simple_port.restype = c_int
snd_seq_create_simple_port.argtypes = [POINTER(snd_seq_t), c_char_p, c_uint, c_uint]

# /usr/include/alsa/seqmid.h:331
snd_seq_delete_simple_port = _lib.snd_seq_delete_simple_port
snd_seq_delete_simple_port.restype = c_int
snd_seq_delete_simple_port.argtypes = [POINTER(snd_seq_t), c_int]

# /usr/include/alsa/seqmid.h:336
snd_seq_connect_from = _lib.snd_seq_connect_from
snd_seq_connect_from.restype = c_int
snd_seq_connect_from.argtypes = [POINTER(snd_seq_t), c_int, c_int, c_int]

# /usr/include/alsa/seqmid.h:337
snd_seq_connect_to = _lib.snd_seq_connect_to
snd_seq_connect_to.restype = c_int
snd_seq_connect_to.argtypes = [POINTER(snd_seq_t), c_int, c_int, c_int]

# /usr/include/alsa/seqmid.h:338
snd_seq_disconnect_from = _lib.snd_seq_disconnect_from
snd_seq_disconnect_from.restype = c_int
snd_seq_disconnect_from.argtypes = [POINTER(snd_seq_t), c_int, c_int, c_int]

# /usr/include/alsa/seqmid.h:339
snd_seq_disconnect_to = _lib.snd_seq_disconnect_to
snd_seq_disconnect_to.restype = c_int
snd_seq_disconnect_to.argtypes = [POINTER(snd_seq_t), c_int, c_int, c_int]

# /usr/include/alsa/seqmid.h:344
snd_seq_set_client_name = _lib.snd_seq_set_client_name
snd_seq_set_client_name.restype = c_int
snd_seq_set_client_name.argtypes = [POINTER(snd_seq_t), c_char_p]

# /usr/include/alsa/seqmid.h:345
snd_seq_set_client_event_filter = _lib.snd_seq_set_client_event_filter
snd_seq_set_client_event_filter.restype = c_int
snd_seq_set_client_event_filter.argtypes = [POINTER(snd_seq_t), c_int]

# /usr/include/alsa/seqmid.h:346
snd_seq_set_client_pool_output = _lib.snd_seq_set_client_pool_output
snd_seq_set_client_pool_output.restype = c_int
snd_seq_set_client_pool_output.argtypes = [POINTER(snd_seq_t), c_size_t]

# /usr/include/alsa/seqmid.h:347
snd_seq_set_client_pool_output_room = _lib.snd_seq_set_client_pool_output_room
snd_seq_set_client_pool_output_room.restype = c_int
snd_seq_set_client_pool_output_room.argtypes = [POINTER(snd_seq_t), c_size_t]

# /usr/include/alsa/seqmid.h:348
snd_seq_set_client_pool_input = _lib.snd_seq_set_client_pool_input
snd_seq_set_client_pool_input.restype = c_int
snd_seq_set_client_pool_input.argtypes = [POINTER(snd_seq_t), c_size_t]

# /usr/include/alsa/seqmid.h:350
snd_seq_sync_output_queue = _lib.snd_seq_sync_output_queue
snd_seq_sync_output_queue.restype = c_int
snd_seq_sync_output_queue.argtypes = [POINTER(snd_seq_t)]

# /usr/include/alsa/seqmid.h:355
snd_seq_parse_address = _lib.snd_seq_parse_address
snd_seq_parse_address.restype = c_int
snd_seq_parse_address.argtypes = [POINTER(snd_seq_t), POINTER(snd_seq_addr_t), c_char_p]

# /usr/include/alsa/seqmid.h:360
snd_seq_reset_pool_output = _lib.snd_seq_reset_pool_output
snd_seq_reset_pool_output.restype = c_int
snd_seq_reset_pool_output.argtypes = [POINTER(snd_seq_t)]

# /usr/include/alsa/seqmid.h:361
snd_seq_reset_pool_input = _lib.snd_seq_reset_pool_input
snd_seq_reset_pool_input.restype = c_int
snd_seq_reset_pool_input.argtypes = [POINTER(snd_seq_t)]

class struct_snd_midi_event(Structure):
    __slots__ = [
    ]
struct_snd_midi_event._fields_ = [
    ('_opaque_struct', c_int)
]

class struct_snd_midi_event(Structure):
    __slots__ = [
    ]
struct_snd_midi_event._fields_ = [
    ('_opaque_struct', c_int)
]

snd_midi_event_t = struct_snd_midi_event 	# /usr/include/alsa/seq_midi_event.h:43
# /usr/include/alsa/seq_midi_event.h:45
snd_midi_event_new = _lib.snd_midi_event_new
snd_midi_event_new.restype = c_int
snd_midi_event_new.argtypes = [c_size_t, POINTER(POINTER(snd_midi_event_t))]

# /usr/include/alsa/seq_midi_event.h:46
snd_midi_event_resize_buffer = _lib.snd_midi_event_resize_buffer
snd_midi_event_resize_buffer.restype = c_int
snd_midi_event_resize_buffer.argtypes = [POINTER(snd_midi_event_t), c_size_t]

# /usr/include/alsa/seq_midi_event.h:47
snd_midi_event_free = _lib.snd_midi_event_free
snd_midi_event_free.restype = None
snd_midi_event_free.argtypes = [POINTER(snd_midi_event_t)]

# /usr/include/alsa/seq_midi_event.h:48
snd_midi_event_init = _lib.snd_midi_event_init
snd_midi_event_init.restype = None
snd_midi_event_init.argtypes = [POINTER(snd_midi_event_t)]

# /usr/include/alsa/seq_midi_event.h:49
snd_midi_event_reset_encode = _lib.snd_midi_event_reset_encode
snd_midi_event_reset_encode.restype = None
snd_midi_event_reset_encode.argtypes = [POINTER(snd_midi_event_t)]

# /usr/include/alsa/seq_midi_event.h:50
snd_midi_event_reset_decode = _lib.snd_midi_event_reset_decode
snd_midi_event_reset_decode.restype = None
snd_midi_event_reset_decode.argtypes = [POINTER(snd_midi_event_t)]

# /usr/include/alsa/seq_midi_event.h:51
snd_midi_event_no_status = _lib.snd_midi_event_no_status
snd_midi_event_no_status.restype = None
snd_midi_event_no_status.argtypes = [POINTER(snd_midi_event_t), c_int]

# /usr/include/alsa/seq_midi_event.h:53
snd_midi_event_encode = _lib.snd_midi_event_encode
snd_midi_event_encode.restype = c_long
snd_midi_event_encode.argtypes = [POINTER(snd_midi_event_t), POINTER(c_ubyte), c_long, POINTER(snd_seq_event_t)]

# /usr/include/alsa/seq_midi_event.h:54
snd_midi_event_encode_byte = _lib.snd_midi_event_encode_byte
snd_midi_event_encode_byte.restype = c_int
snd_midi_event_encode_byte.argtypes = [POINTER(snd_midi_event_t), c_int, POINTER(snd_seq_event_t)]

# /usr/include/alsa/seq_midi_event.h:56
snd_midi_event_decode = _lib.snd_midi_event_decode
snd_midi_event_decode.restype = c_long
snd_midi_event_decode.argtypes = [POINTER(snd_midi_event_t), POINTER(c_ubyte), c_long, POINTER(snd_seq_event_t)]

class struct__snd_instr_header(Structure):
    __slots__ = [
    ]
struct__snd_instr_header._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_instr_header(Structure):
    __slots__ = [
    ]
struct__snd_instr_header._fields_ = [
    ('_opaque_struct', c_int)
]

snd_instr_header_t = struct__snd_instr_header 	# /usr/include/alsa/instr.h:44
# /usr/include/alsa/instr.h:46
snd_instr_header_sizeof = _lib.snd_instr_header_sizeof
snd_instr_header_sizeof.restype = c_size_t
snd_instr_header_sizeof.argtypes = []

# /usr/include/alsa/instr.h:53
snd_instr_header_malloc = _lib.snd_instr_header_malloc
snd_instr_header_malloc.restype = c_int
snd_instr_header_malloc.argtypes = [POINTER(POINTER(snd_instr_header_t)), c_size_t]

# /usr/include/alsa/instr.h:54
snd_instr_header_free = _lib.snd_instr_header_free
snd_instr_header_free.restype = None
snd_instr_header_free.argtypes = [POINTER(snd_instr_header_t)]

# /usr/include/alsa/instr.h:55
snd_instr_header_copy = _lib.snd_instr_header_copy
snd_instr_header_copy.restype = None
snd_instr_header_copy.argtypes = [POINTER(snd_instr_header_t), POINTER(snd_instr_header_t)]

# /usr/include/alsa/instr.h:57
snd_instr_header_get_id = _lib.snd_instr_header_get_id
snd_instr_header_get_id.restype = POINTER(snd_seq_instr_t)
snd_instr_header_get_id.argtypes = [POINTER(snd_instr_header_t)]

# /usr/include/alsa/instr.h:58
snd_instr_header_get_cluster = _lib.snd_instr_header_get_cluster
snd_instr_header_get_cluster.restype = snd_seq_instr_cluster_t
snd_instr_header_get_cluster.argtypes = [POINTER(snd_instr_header_t)]

# /usr/include/alsa/instr.h:59
snd_instr_header_get_cmd = _lib.snd_instr_header_get_cmd
snd_instr_header_get_cmd.restype = c_uint
snd_instr_header_get_cmd.argtypes = [POINTER(snd_instr_header_t)]

# /usr/include/alsa/instr.h:60
snd_instr_header_get_len = _lib.snd_instr_header_get_len
snd_instr_header_get_len.restype = c_size_t
snd_instr_header_get_len.argtypes = [POINTER(snd_instr_header_t)]

# /usr/include/alsa/instr.h:61
snd_instr_header_get_name = _lib.snd_instr_header_get_name
snd_instr_header_get_name.restype = c_char_p
snd_instr_header_get_name.argtypes = [POINTER(snd_instr_header_t)]

# /usr/include/alsa/instr.h:62
snd_instr_header_get_type = _lib.snd_instr_header_get_type
snd_instr_header_get_type.restype = c_int
snd_instr_header_get_type.argtypes = [POINTER(snd_instr_header_t)]

# /usr/include/alsa/instr.h:63
snd_instr_header_get_format = _lib.snd_instr_header_get_format
snd_instr_header_get_format.restype = c_char_p
snd_instr_header_get_format.argtypes = [POINTER(snd_instr_header_t)]

# /usr/include/alsa/instr.h:64
snd_instr_header_get_alias = _lib.snd_instr_header_get_alias
snd_instr_header_get_alias.restype = POINTER(snd_seq_instr_t)
snd_instr_header_get_alias.argtypes = [POINTER(snd_instr_header_t)]

# /usr/include/alsa/instr.h:65
snd_instr_header_get_data = _lib.snd_instr_header_get_data
snd_instr_header_get_data.restype = POINTER(c_void)
snd_instr_header_get_data.argtypes = [POINTER(snd_instr_header_t)]

# /usr/include/alsa/instr.h:66
snd_instr_header_get_follow_alias = _lib.snd_instr_header_get_follow_alias
snd_instr_header_get_follow_alias.restype = c_int
snd_instr_header_get_follow_alias.argtypes = [POINTER(snd_instr_header_t)]

# /usr/include/alsa/instr.h:68
snd_instr_header_set_id = _lib.snd_instr_header_set_id
snd_instr_header_set_id.restype = None
snd_instr_header_set_id.argtypes = [POINTER(snd_instr_header_t), POINTER(snd_seq_instr_t)]

# /usr/include/alsa/instr.h:69
snd_instr_header_set_cluster = _lib.snd_instr_header_set_cluster
snd_instr_header_set_cluster.restype = None
snd_instr_header_set_cluster.argtypes = [POINTER(snd_instr_header_t), snd_seq_instr_cluster_t]

# /usr/include/alsa/instr.h:70
snd_instr_header_set_cmd = _lib.snd_instr_header_set_cmd
snd_instr_header_set_cmd.restype = None
snd_instr_header_set_cmd.argtypes = [POINTER(snd_instr_header_t), c_uint]

# /usr/include/alsa/instr.h:71
snd_instr_header_set_len = _lib.snd_instr_header_set_len
snd_instr_header_set_len.restype = None
snd_instr_header_set_len.argtypes = [POINTER(snd_instr_header_t), c_size_t]

# /usr/include/alsa/instr.h:72
snd_instr_header_set_name = _lib.snd_instr_header_set_name
snd_instr_header_set_name.restype = None
snd_instr_header_set_name.argtypes = [POINTER(snd_instr_header_t), c_char_p]

# /usr/include/alsa/instr.h:73
snd_instr_header_set_type = _lib.snd_instr_header_set_type
snd_instr_header_set_type.restype = None
snd_instr_header_set_type.argtypes = [POINTER(snd_instr_header_t), c_int]

# /usr/include/alsa/instr.h:74
snd_instr_header_set_format = _lib.snd_instr_header_set_format
snd_instr_header_set_format.restype = None
snd_instr_header_set_format.argtypes = [POINTER(snd_instr_header_t), c_char_p]

# /usr/include/alsa/instr.h:75
snd_instr_header_set_alias = _lib.snd_instr_header_set_alias
snd_instr_header_set_alias.restype = None
snd_instr_header_set_alias.argtypes = [POINTER(snd_instr_header_t), POINTER(snd_seq_instr_t)]

# /usr/include/alsa/instr.h:76
snd_instr_header_set_follow_alias = _lib.snd_instr_header_set_follow_alias
snd_instr_header_set_follow_alias.restype = None
snd_instr_header_set_follow_alias.argtypes = [POINTER(snd_instr_header_t), c_int]

SND_SEQ_INSTR_ATYPE_DATA = 0 	# /usr/include/alsa/instr.h:84
SND_SEQ_INSTR_ATYPE_ALIAS = 1 	# /usr/include/alsa/instr.h:85
SND_SEQ_INSTR_TYPE0_DLS1 = 1 	# /usr/include/alsa/instr.h:98
SND_SEQ_INSTR_TYPE0_DLS2 = 2 	# /usr/include/alsa/instr.h:99
SND_SEQ_INSTR_TYPE1_SIMPLE = 1 	# /usr/include/alsa/instr.h:100
SND_SEQ_INSTR_TYPE1_SOUNDFONT = 2 	# /usr/include/alsa/instr.h:101
SND_SEQ_INSTR_TYPE1_GUS_PATCH = 4 	# /usr/include/alsa/instr.h:102
SND_SEQ_INSTR_TYPE1_INTERWAVE = 8 	# /usr/include/alsa/instr.h:103
SND_SEQ_INSTR_TYPE2_OPL2_3 = 1 	# /usr/include/alsa/instr.h:104
SND_SEQ_INSTR_TYPE2_OPL4 = 2 	# /usr/include/alsa/instr.h:105
SND_SEQ_INSTR_PUT_CMD_CREATE = 0 	# /usr/include/alsa/instr.h:108
SND_SEQ_INSTR_PUT_CMD_REPLACE = 1 	# /usr/include/alsa/instr.h:109
SND_SEQ_INSTR_PUT_CMD_MODIFY = 2 	# /usr/include/alsa/instr.h:110
SND_SEQ_INSTR_PUT_CMD_ADD = 3 	# /usr/include/alsa/instr.h:111
SND_SEQ_INSTR_PUT_CMD_REMOVE = 4 	# /usr/include/alsa/instr.h:112
SND_SEQ_INSTR_GET_CMD_FULL = 0 	# /usr/include/alsa/instr.h:115
SND_SEQ_INSTR_GET_CMD_PARTIAL = 1 	# /usr/include/alsa/instr.h:116
SND_SEQ_INSTR_QUERY_FOLLOW_ALIAS = 1 	# /usr/include/alsa/instr.h:119
SND_SEQ_INSTR_FREE_CMD_ALL = 0 	# /usr/include/alsa/instr.h:122
SND_SEQ_INSTR_FREE_CMD_PRIVATE = 1 	# /usr/include/alsa/instr.h:123
SND_SEQ_INSTR_FREE_CMD_CLUSTER = 2 	# /usr/include/alsa/instr.h:124
SND_SEQ_INSTR_FREE_CMD_SINGLE = 3 	# /usr/include/alsa/instr.h:125
snd_instr_fm_t = None 	# /usr/include/alsa/instr.h:133
# /usr/include/alsa/instr.h:135
snd_instr_fm_convert_to_stream = _lib.snd_instr_fm_convert_to_stream
snd_instr_fm_convert_to_stream.restype = c_int
snd_instr_fm_convert_to_stream.argtypes = [POINTER(snd_instr_fm_t), c_char_p, POINTER(POINTER(snd_instr_header_t)), POINTER(c_size_t)]

# /usr/include/alsa/instr.h:136
snd_instr_fm_convert_from_stream = _lib.snd_instr_fm_convert_from_stream
snd_instr_fm_convert_from_stream.restype = c_int
snd_instr_fm_convert_from_stream.argtypes = [POINTER(snd_instr_header_t), c_size_t, POINTER(POINTER(snd_instr_fm_t))]

# /usr/include/alsa/instr.h:137
snd_instr_fm_free = _lib.snd_instr_fm_free
snd_instr_fm_free.restype = c_int
snd_instr_fm_free.argtypes = [POINTER(snd_instr_fm_t)]

snd_instr_simple_t = None 	# /usr/include/alsa/instr.h:145
# /usr/include/alsa/instr.h:147
snd_instr_simple_convert_to_stream = _lib.snd_instr_simple_convert_to_stream
snd_instr_simple_convert_to_stream.restype = c_int
snd_instr_simple_convert_to_stream.argtypes = [POINTER(snd_instr_simple_t), c_char_p, POINTER(POINTER(snd_instr_header_t)), POINTER(c_size_t)]

# /usr/include/alsa/instr.h:148
snd_instr_simple_convert_from_stream = _lib.snd_instr_simple_convert_from_stream
snd_instr_simple_convert_from_stream.restype = c_int
snd_instr_simple_convert_from_stream.argtypes = [POINTER(snd_instr_header_t), c_size_t, POINTER(POINTER(snd_instr_simple_t))]

# /usr/include/alsa/instr.h:149
snd_instr_simple_free = _lib.snd_instr_simple_free
snd_instr_simple_free.restype = c_int
snd_instr_simple_free.argtypes = [POINTER(snd_instr_simple_t)]

snd_instr_iwffff_t = None 	# /usr/include/alsa/instr.h:157
class struct__snd_iwffff_handle(Structure):
    __slots__ = [
    ]
struct__snd_iwffff_handle._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__snd_iwffff_handle(Structure):
    __slots__ = [
    ]
struct__snd_iwffff_handle._fields_ = [
    ('_opaque_struct', c_int)
]

snd_iwffff_handle_t = struct__snd_iwffff_handle 	# /usr/include/alsa/instr.h:159
# /usr/include/alsa/instr.h:161
snd_instr_iwffff_open = _lib.snd_instr_iwffff_open
snd_instr_iwffff_open.restype = c_int
snd_instr_iwffff_open.argtypes = [POINTER(POINTER(snd_iwffff_handle_t)), c_char_p, c_char_p]

# /usr/include/alsa/instr.h:162
snd_instr_iwffff_open_rom = _lib.snd_instr_iwffff_open_rom
snd_instr_iwffff_open_rom.restype = c_int
snd_instr_iwffff_open_rom.argtypes = [POINTER(POINTER(snd_iwffff_handle_t)), c_int, c_int, c_int]

# /usr/include/alsa/instr.h:163
snd_instr_iwffff_open_rom_file = _lib.snd_instr_iwffff_open_rom_file
snd_instr_iwffff_open_rom_file.restype = c_int
snd_instr_iwffff_open_rom_file.argtypes = [POINTER(POINTER(snd_iwffff_handle_t)), c_char_p, c_int, c_int]

# /usr/include/alsa/instr.h:164
snd_instr_iwffff_close = _lib.snd_instr_iwffff_close
snd_instr_iwffff_close.restype = c_int
snd_instr_iwffff_close.argtypes = [POINTER(snd_iwffff_handle_t)]

# /usr/include/alsa/instr.h:165
snd_instr_iwffff_load = _lib.snd_instr_iwffff_load
snd_instr_iwffff_load.restype = c_int
snd_instr_iwffff_load.argtypes = [POINTER(snd_iwffff_handle_t), c_int, c_int, POINTER(POINTER(snd_instr_iwffff_t))]

# /usr/include/alsa/instr.h:166
snd_instr_iwffff_convert_to_stream = _lib.snd_instr_iwffff_convert_to_stream
snd_instr_iwffff_convert_to_stream.restype = c_int
snd_instr_iwffff_convert_to_stream.argtypes = [POINTER(snd_instr_iwffff_t), c_char_p, POINTER(POINTER(snd_instr_header_t)), POINTER(c_size_t)]

# /usr/include/alsa/instr.h:167
snd_instr_iwffff_convert_from_stream = _lib.snd_instr_iwffff_convert_from_stream
snd_instr_iwffff_convert_from_stream.restype = c_int
snd_instr_iwffff_convert_from_stream.argtypes = [POINTER(snd_instr_header_t), c_size_t, POINTER(POINTER(snd_instr_iwffff_t))]

# /usr/include/alsa/instr.h:168
snd_instr_iwffff_free = _lib.snd_instr_iwffff_free
snd_instr_iwffff_free.restype = c_int
snd_instr_iwffff_free.argtypes = [POINTER(snd_instr_iwffff_t)]


__all__ = ['alsa_lisp_default_cfg_free', 'alsa_lisp', 'alsa_lisp_free',
'alsa_lisp_result_free', 'alsa_lisp_seq_first', 'alsa_lisp_seq_next',
'alsa_lisp_seq_count', 'alsa_lisp_seq_integer', 'alsa_lisp_seq_pointer',
'IEC958_AES0_PROFESSIONAL', 'IEC958_AES0_NONAUDIO',
'IEC958_AES0_PRO_EMPHASIS', 'IEC958_AES0_PRO_EMPHASIS_NOTID',
'IEC958_AES0_PRO_EMPHASIS_NONE', 'IEC958_AES0_PRO_EMPHASIS_5015',
'IEC958_AES0_PRO_EMPHASIS_CCITT', 'IEC958_AES0_PRO_FREQ_UNLOCKED',
'IEC958_AES0_PRO_FS', 'IEC958_AES0_PRO_FS_NOTID', 'IEC958_AES0_PRO_FS_44100',
'IEC958_AES0_PRO_FS_48000', 'IEC958_AES0_PRO_FS_32000',
'IEC958_AES0_CON_NOT_COPYRIGHT', 'IEC958_AES0_CON_EMPHASIS',
'IEC958_AES0_CON_EMPHASIS_NONE', 'IEC958_AES0_CON_EMPHASIS_5015',
'IEC958_AES0_CON_MODE', 'IEC958_AES1_PRO_MODE', 'IEC958_AES1_PRO_MODE_NOTID',
'IEC958_AES1_PRO_MODE_STEREOPHONIC', 'IEC958_AES1_PRO_MODE_SINGLE',
'IEC958_AES1_PRO_MODE_TWO', 'IEC958_AES1_PRO_MODE_PRIMARY',
'IEC958_AES1_PRO_MODE_BYTE3', 'IEC958_AES1_PRO_USERBITS',
'IEC958_AES1_PRO_USERBITS_NOTID', 'IEC958_AES1_PRO_USERBITS_192',
'IEC958_AES1_PRO_USERBITS_UDEF', 'IEC958_AES1_CON_CATEGORY',
'IEC958_AES1_CON_GENERAL', 'IEC958_AES1_CON_EXPERIMENTAL',
'IEC958_AES1_CON_SOLIDMEM_MASK', 'IEC958_AES1_CON_SOLIDMEM_ID',
'IEC958_AES1_CON_BROADCAST1_MASK', 'IEC958_AES1_CON_BROADCAST1_ID',
'IEC958_AES1_CON_DIGDIGCONV_MASK', 'IEC958_AES1_CON_DIGDIGCONV_ID',
'IEC958_AES1_CON_ADC_COPYRIGHT_MASK', 'IEC958_AES1_CON_ADC_COPYRIGHT_ID',
'IEC958_AES1_CON_ADC_MASK', 'IEC958_AES1_CON_ADC_ID',
'IEC958_AES1_CON_BROADCAST2_MASK', 'IEC958_AES1_CON_BROADCAST2_ID',
'IEC958_AES1_CON_LASEROPT_MASK', 'IEC958_AES1_CON_LASEROPT_ID',
'IEC958_AES1_CON_MUSICAL_MASK', 'IEC958_AES1_CON_MUSICAL_ID',
'IEC958_AES1_CON_MAGNETIC_MASK', 'IEC958_AES1_CON_MAGNETIC_ID',
'IEC958_AES1_CON_IEC908_CD', 'IEC958_AES1_CON_NON_IEC908_CD',
'IEC958_AES1_CON_PCM_CODER', 'IEC958_AES1_CON_SAMPLER',
'IEC958_AES1_CON_MIXER', 'IEC958_AES1_CON_RATE_CONVERTER',
'IEC958_AES1_CON_SYNTHESIZER', 'IEC958_AES1_CON_MICROPHONE',
'IEC958_AES1_CON_DAT', 'IEC958_AES1_CON_VCR', 'IEC958_AES1_CON_ORIGINAL',
'IEC958_AES2_PRO_SBITS', 'IEC958_AES2_PRO_SBITS_20',
'IEC958_AES2_PRO_SBITS_24', 'IEC958_AES2_PRO_SBITS_UDEF',
'IEC958_AES2_PRO_WORDLEN', 'IEC958_AES2_PRO_WORDLEN_NOTID',
'IEC958_AES2_PRO_WORDLEN_22_18', 'IEC958_AES2_PRO_WORDLEN_23_19',
'IEC958_AES2_PRO_WORDLEN_24_20', 'IEC958_AES2_PRO_WORDLEN_20_16',
'IEC958_AES2_CON_SOURCE', 'IEC958_AES2_CON_SOURCE_UNSPEC',
'IEC958_AES2_CON_CHANNEL', 'IEC958_AES2_CON_CHANNEL_UNSPEC',
'IEC958_AES3_CON_FS', 'IEC958_AES3_CON_FS_44100', 'IEC958_AES3_CON_FS_48000',
'IEC958_AES3_CON_FS_32000', 'IEC958_AES3_CON_CLOCK',
'IEC958_AES3_CON_CLOCK_1000PPM', 'IEC958_AES3_CON_CLOCK_50PPM',
'IEC958_AES3_CON_CLOCK_VARIABLE', 'MIDI_CHANNELS', 'MIDI_GM_DRUM_CHANNEL',
'MIDI_CMD_NOTE_OFF', 'MIDI_CMD_NOTE_ON', 'MIDI_CMD_NOTE_PRESSURE',
'MIDI_CMD_CONTROL', 'MIDI_CMD_PGM_CHANGE', 'MIDI_CMD_CHANNEL_PRESSURE',
'MIDI_CMD_BENDER', 'MIDI_CMD_COMMON_SYSEX', 'MIDI_CMD_COMMON_MTC_QUARTER',
'MIDI_CMD_COMMON_SONG_POS', 'MIDI_CMD_COMMON_SONG_SELECT',
'MIDI_CMD_COMMON_TUNE_REQUEST', 'MIDI_CMD_COMMON_SYSEX_END',
'MIDI_CMD_COMMON_CLOCK', 'MIDI_CMD_COMMON_START', 'MIDI_CMD_COMMON_CONTINUE',
'MIDI_CMD_COMMON_STOP', 'MIDI_CMD_COMMON_SENSING', 'MIDI_CMD_COMMON_RESET',
'MIDI_CTL_MSB_BANK', 'MIDI_CTL_MSB_MODWHEEL', 'MIDI_CTL_MSB_BREATH',
'MIDI_CTL_MSB_FOOT', 'MIDI_CTL_MSB_PORTAMENTO_TIME',
'MIDI_CTL_MSB_DATA_ENTRY', 'MIDI_CTL_MSB_MAIN_VOLUME', 'MIDI_CTL_MSB_BALANCE',
'MIDI_CTL_MSB_PAN', 'MIDI_CTL_MSB_EXPRESSION', 'MIDI_CTL_MSB_EFFECT1',
'MIDI_CTL_MSB_EFFECT2', 'MIDI_CTL_MSB_GENERAL_PURPOSE1',
'MIDI_CTL_MSB_GENERAL_PURPOSE2', 'MIDI_CTL_MSB_GENERAL_PURPOSE3',
'MIDI_CTL_MSB_GENERAL_PURPOSE4', 'MIDI_CTL_LSB_BANK', 'MIDI_CTL_LSB_MODWHEEL',
'MIDI_CTL_LSB_BREATH', 'MIDI_CTL_LSB_FOOT', 'MIDI_CTL_LSB_PORTAMENTO_TIME',
'MIDI_CTL_LSB_DATA_ENTRY', 'MIDI_CTL_LSB_MAIN_VOLUME', 'MIDI_CTL_LSB_BALANCE',
'MIDI_CTL_LSB_PAN', 'MIDI_CTL_LSB_EXPRESSION', 'MIDI_CTL_LSB_EFFECT1',
'MIDI_CTL_LSB_EFFECT2', 'MIDI_CTL_LSB_GENERAL_PURPOSE1',
'MIDI_CTL_LSB_GENERAL_PURPOSE2', 'MIDI_CTL_LSB_GENERAL_PURPOSE3',
'MIDI_CTL_LSB_GENERAL_PURPOSE4', 'MIDI_CTL_SUSTAIN', 'MIDI_CTL_PORTAMENTO',
'MIDI_CTL_SOSTENUTO', 'MIDI_CTL_SUSTENUTO', 'MIDI_CTL_SOFT_PEDAL',
'MIDI_CTL_LEGATO_FOOTSWITCH', 'MIDI_CTL_HOLD2',
'MIDI_CTL_SC1_SOUND_VARIATION', 'MIDI_CTL_SC2_TIMBRE',
'MIDI_CTL_SC3_RELEASE_TIME', 'MIDI_CTL_SC4_ATTACK_TIME',
'MIDI_CTL_SC5_BRIGHTNESS', 'MIDI_CTL_SC6', 'MIDI_CTL_SC7', 'MIDI_CTL_SC8',
'MIDI_CTL_SC9', 'MIDI_CTL_SC10', 'MIDI_CTL_GENERAL_PURPOSE5',
'MIDI_CTL_GENERAL_PURPOSE6', 'MIDI_CTL_GENERAL_PURPOSE7',
'MIDI_CTL_GENERAL_PURPOSE8', 'MIDI_CTL_PORTAMENTO_CONTROL',
'MIDI_CTL_E1_REVERB_DEPTH', 'MIDI_CTL_E2_TREMOLO_DEPTH',
'MIDI_CTL_E3_CHORUS_DEPTH', 'MIDI_CTL_E4_DETUNE_DEPTH',
'MIDI_CTL_E5_PHASER_DEPTH', 'MIDI_CTL_DATA_INCREMENT',
'MIDI_CTL_DATA_DECREMENT', 'MIDI_CTL_NONREG_PARM_NUM_LSB',
'MIDI_CTL_NONREG_PARM_NUM_MSB', 'MIDI_CTL_REGIST_PARM_NUM_LSB',
'MIDI_CTL_REGIST_PARM_NUM_MSB', 'MIDI_CTL_ALL_SOUNDS_OFF',
'MIDI_CTL_RESET_CONTROLLERS', 'MIDI_CTL_LOCAL_CONTROL_SWITCH',
'MIDI_CTL_ALL_NOTES_OFF', 'MIDI_CTL_OMNI_OFF', 'MIDI_CTL_OMNI_ON',
'MIDI_CTL_MONO1', 'MIDI_CTL_MONO2', 'IEC958_AES0_PROFESSIONAL',
'IEC958_AES0_NONAUDIO', 'IEC958_AES0_PRO_EMPHASIS',
'IEC958_AES0_PRO_EMPHASIS_NOTID', 'IEC958_AES0_PRO_EMPHASIS_NONE',
'IEC958_AES0_PRO_EMPHASIS_5015', 'IEC958_AES0_PRO_EMPHASIS_CCITT',
'IEC958_AES0_PRO_FREQ_UNLOCKED', 'IEC958_AES0_PRO_FS',
'IEC958_AES0_PRO_FS_NOTID', 'IEC958_AES0_PRO_FS_44100',
'IEC958_AES0_PRO_FS_48000', 'IEC958_AES0_PRO_FS_32000',
'IEC958_AES0_CON_NOT_COPYRIGHT', 'IEC958_AES0_CON_EMPHASIS',
'IEC958_AES0_CON_EMPHASIS_NONE', 'IEC958_AES0_CON_EMPHASIS_5015',
'IEC958_AES0_CON_MODE', 'IEC958_AES1_PRO_MODE', 'IEC958_AES1_PRO_MODE_NOTID',
'IEC958_AES1_PRO_MODE_STEREOPHONIC', 'IEC958_AES1_PRO_MODE_SINGLE',
'IEC958_AES1_PRO_MODE_TWO', 'IEC958_AES1_PRO_MODE_PRIMARY',
'IEC958_AES1_PRO_MODE_BYTE3', 'IEC958_AES1_PRO_USERBITS',
'IEC958_AES1_PRO_USERBITS_NOTID', 'IEC958_AES1_PRO_USERBITS_192',
'IEC958_AES1_PRO_USERBITS_UDEF', 'IEC958_AES1_CON_CATEGORY',
'IEC958_AES1_CON_GENERAL', 'IEC958_AES1_CON_EXPERIMENTAL',
'IEC958_AES1_CON_SOLIDMEM_MASK', 'IEC958_AES1_CON_SOLIDMEM_ID',
'IEC958_AES1_CON_BROADCAST1_MASK', 'IEC958_AES1_CON_BROADCAST1_ID',
'IEC958_AES1_CON_DIGDIGCONV_MASK', 'IEC958_AES1_CON_DIGDIGCONV_ID',
'IEC958_AES1_CON_ADC_COPYRIGHT_MASK', 'IEC958_AES1_CON_ADC_COPYRIGHT_ID',
'IEC958_AES1_CON_ADC_MASK', 'IEC958_AES1_CON_ADC_ID',
'IEC958_AES1_CON_BROADCAST2_MASK', 'IEC958_AES1_CON_BROADCAST2_ID',
'IEC958_AES1_CON_LASEROPT_MASK', 'IEC958_AES1_CON_LASEROPT_ID',
'IEC958_AES1_CON_MUSICAL_MASK', 'IEC958_AES1_CON_MUSICAL_ID',
'IEC958_AES1_CON_MAGNETIC_MASK', 'IEC958_AES1_CON_MAGNETIC_ID',
'IEC958_AES1_CON_IEC908_CD', 'IEC958_AES1_CON_NON_IEC908_CD',
'IEC958_AES1_CON_PCM_CODER', 'IEC958_AES1_CON_SAMPLER',
'IEC958_AES1_CON_MIXER', 'IEC958_AES1_CON_RATE_CONVERTER',
'IEC958_AES1_CON_SYNTHESIZER', 'IEC958_AES1_CON_MICROPHONE',
'IEC958_AES1_CON_DAT', 'IEC958_AES1_CON_VCR', 'IEC958_AES1_CON_ORIGINAL',
'IEC958_AES2_PRO_SBITS', 'IEC958_AES2_PRO_SBITS_20',
'IEC958_AES2_PRO_SBITS_24', 'IEC958_AES2_PRO_SBITS_UDEF',
'IEC958_AES2_PRO_WORDLEN', 'IEC958_AES2_PRO_WORDLEN_NOTID',
'IEC958_AES2_PRO_WORDLEN_22_18', 'IEC958_AES2_PRO_WORDLEN_23_19',
'IEC958_AES2_PRO_WORDLEN_24_20', 'IEC958_AES2_PRO_WORDLEN_20_16',
'IEC958_AES2_CON_SOURCE', 'IEC958_AES2_CON_SOURCE_UNSPEC',
'IEC958_AES2_CON_CHANNEL', 'IEC958_AES2_CON_CHANNEL_UNSPEC',
'IEC958_AES3_CON_FS', 'IEC958_AES3_CON_FS_44100', 'IEC958_AES3_CON_FS_48000',
'IEC958_AES3_CON_FS_32000', 'IEC958_AES3_CON_CLOCK',
'IEC958_AES3_CON_CLOCK_1000PPM', 'IEC958_AES3_CON_CLOCK_50PPM',
'IEC958_AES3_CON_CLOCK_VARIABLE', 'MIDI_CHANNELS', 'MIDI_GM_DRUM_CHANNEL',
'MIDI_CMD_NOTE_OFF', 'MIDI_CMD_NOTE_ON', 'MIDI_CMD_NOTE_PRESSURE',
'MIDI_CMD_CONTROL', 'MIDI_CMD_PGM_CHANGE', 'MIDI_CMD_CHANNEL_PRESSURE',
'MIDI_CMD_BENDER', 'MIDI_CMD_COMMON_SYSEX', 'MIDI_CMD_COMMON_MTC_QUARTER',
'MIDI_CMD_COMMON_SONG_POS', 'MIDI_CMD_COMMON_SONG_SELECT',
'MIDI_CMD_COMMON_TUNE_REQUEST', 'MIDI_CMD_COMMON_SYSEX_END',
'MIDI_CMD_COMMON_CLOCK', 'MIDI_CMD_COMMON_START', 'MIDI_CMD_COMMON_CONTINUE',
'MIDI_CMD_COMMON_STOP', 'MIDI_CMD_COMMON_SENSING', 'MIDI_CMD_COMMON_RESET',
'MIDI_CTL_MSB_BANK', 'MIDI_CTL_MSB_MODWHEEL', 'MIDI_CTL_MSB_BREATH',
'MIDI_CTL_MSB_FOOT', 'MIDI_CTL_MSB_PORTAMENTO_TIME',
'MIDI_CTL_MSB_DATA_ENTRY', 'MIDI_CTL_MSB_MAIN_VOLUME', 'MIDI_CTL_MSB_BALANCE',
'MIDI_CTL_MSB_PAN', 'MIDI_CTL_MSB_EXPRESSION', 'MIDI_CTL_MSB_EFFECT1',
'MIDI_CTL_MSB_EFFECT2', 'MIDI_CTL_MSB_GENERAL_PURPOSE1',
'MIDI_CTL_MSB_GENERAL_PURPOSE2', 'MIDI_CTL_MSB_GENERAL_PURPOSE3',
'MIDI_CTL_MSB_GENERAL_PURPOSE4', 'MIDI_CTL_LSB_BANK', 'MIDI_CTL_LSB_MODWHEEL',
'MIDI_CTL_LSB_BREATH', 'MIDI_CTL_LSB_FOOT', 'MIDI_CTL_LSB_PORTAMENTO_TIME',
'MIDI_CTL_LSB_DATA_ENTRY', 'MIDI_CTL_LSB_MAIN_VOLUME', 'MIDI_CTL_LSB_BALANCE',
'MIDI_CTL_LSB_PAN', 'MIDI_CTL_LSB_EXPRESSION', 'MIDI_CTL_LSB_EFFECT1',
'MIDI_CTL_LSB_EFFECT2', 'MIDI_CTL_LSB_GENERAL_PURPOSE1',
'MIDI_CTL_LSB_GENERAL_PURPOSE2', 'MIDI_CTL_LSB_GENERAL_PURPOSE3',
'MIDI_CTL_LSB_GENERAL_PURPOSE4', 'MIDI_CTL_SUSTAIN', 'MIDI_CTL_PORTAMENTO',
'MIDI_CTL_SOSTENUTO', 'MIDI_CTL_SUSTENUTO', 'MIDI_CTL_SOFT_PEDAL',
'MIDI_CTL_LEGATO_FOOTSWITCH', 'MIDI_CTL_HOLD2',
'MIDI_CTL_SC1_SOUND_VARIATION', 'MIDI_CTL_SC2_TIMBRE',
'MIDI_CTL_SC3_RELEASE_TIME', 'MIDI_CTL_SC4_ATTACK_TIME',
'MIDI_CTL_SC5_BRIGHTNESS', 'MIDI_CTL_SC6', 'MIDI_CTL_SC7', 'MIDI_CTL_SC8',
'MIDI_CTL_SC9', 'MIDI_CTL_SC10', 'MIDI_CTL_GENERAL_PURPOSE5',
'MIDI_CTL_GENERAL_PURPOSE6', 'MIDI_CTL_GENERAL_PURPOSE7',
'MIDI_CTL_GENERAL_PURPOSE8', 'MIDI_CTL_PORTAMENTO_CONTROL',
'MIDI_CTL_E1_REVERB_DEPTH', 'MIDI_CTL_E2_TREMOLO_DEPTH',
'MIDI_CTL_E3_CHORUS_DEPTH', 'MIDI_CTL_E4_DETUNE_DEPTH',
'MIDI_CTL_E5_PHASER_DEPTH', 'MIDI_CTL_DATA_INCREMENT',
'MIDI_CTL_DATA_DECREMENT', 'MIDI_CTL_NONREG_PARM_NUM_LSB',
'MIDI_CTL_NONREG_PARM_NUM_MSB', 'MIDI_CTL_REGIST_PARM_NUM_LSB',
'MIDI_CTL_REGIST_PARM_NUM_MSB', 'MIDI_CTL_ALL_SOUNDS_OFF',
'MIDI_CTL_RESET_CONTROLLERS', 'MIDI_CTL_LOCAL_CONTROL_SWITCH',
'MIDI_CTL_ALL_NOTES_OFF', 'MIDI_CTL_OMNI_OFF', 'MIDI_CTL_OMNI_ON',
'MIDI_CTL_MONO1', 'MIDI_CTL_MONO2', 'SND_LIB_MAJOR', 'SND_LIB_MINOR',
'SND_LIB_SUBMINOR', 'SND_LIB_EXTRAVER', 'SND_LIB_VERSION',
'snd_asoundlib_version', 'snd_dlopen', 'snd_dlsym', 'snd_dlclose',
'snd_async_handler_t', 'snd_async_callback_t', 'snd_async_add_handler',
'snd_async_del_handler', 'snd_async_handler_get_fd',
'snd_async_handler_get_signo', 'snd_async_handler_get_callback_private',
'snd_shm_area_create', 'snd_shm_area_share', 'snd_shm_area_destroy',
'snd_user_file', 'snd_timestamp_t', 'snd_htimestamp_t', 'snd_input_t',
'snd_input_type_t', 'SND_INPUT_STDIO', 'SND_INPUT_BUFFER',
'snd_input_stdio_open', 'snd_input_stdio_attach', 'snd_input_buffer_open',
'snd_input_close', 'snd_input_scanf', 'snd_input_gets', 'snd_input_getc',
'snd_input_ungetc', 'SND_ERROR_BEGIN', 'SND_ERROR_INCOMPATIBLE_VERSION',
'SND_ERROR_ALISP_NIL', 'snd_strerror', 'snd_lib_error_handler_t',
'snd_lib_error_set_handler', 'SND_CONFIG_DLSYM_VERSION_EVALUATE',
'SND_CONFIG_DLSYM_VERSION_HOOK', 'snd_config_type_t',
'SND_CONFIG_TYPE_INTEGER', 'SND_CONFIG_TYPE_INTEGER64',
'SND_CONFIG_TYPE_REAL', 'SND_CONFIG_TYPE_STRING', 'SND_CONFIG_TYPE_POINTER',
'SND_CONFIG_TYPE_COMPOUND', 'snd_config_t', 'snd_config_iterator_t',
'snd_config_update_t', 'snd_config_top', 'snd_config_load',
'snd_config_load_override', 'snd_config_save', 'snd_config_update',
'snd_config_update_r', 'snd_config_update_free',
'snd_config_update_free_global', 'snd_config_search', 'snd_config_searchv',
'snd_config_search_definition', 'snd_config_expand', 'snd_config_evaluate',
'snd_config_add', 'snd_config_delete', 'snd_config_delete_compound_members',
'snd_config_copy', 'snd_config_make', 'snd_config_make_integer',
'snd_config_make_integer64', 'snd_config_make_real', 'snd_config_make_string',
'snd_config_make_pointer', 'snd_config_make_compound',
'snd_config_imake_integer', 'snd_config_imake_integer64',
'snd_config_imake_real', 'snd_config_imake_string',
'snd_config_imake_pointer', 'snd_config_get_type', 'snd_config_set_id',
'snd_config_set_integer', 'snd_config_set_integer64', 'snd_config_set_real',
'snd_config_set_string', 'snd_config_set_ascii', 'snd_config_set_pointer',
'snd_config_get_id', 'snd_config_get_integer', 'snd_config_get_integer64',
'snd_config_get_real', 'snd_config_get_ireal', 'snd_config_get_string',
'snd_config_get_ascii', 'snd_config_get_pointer', 'snd_config_test_id',
'snd_config_iterator_first', 'snd_config_iterator_next',
'snd_config_iterator_end', 'snd_config_iterator_entry',
'snd_config_get_bool_ascii', 'snd_config_get_bool',
'snd_config_get_ctl_iface_ascii', 'snd_config_get_ctl_iface', 'snd_devname_t',
'snd_names_list', 'snd_names_list_free', 'SND_PCM_DLSYM_VERSION',
'snd_pcm_info_t', 'snd_pcm_hw_params_t', 'snd_pcm_sw_params_t',
'snd_pcm_status_t', 'snd_pcm_access_mask_t', 'snd_pcm_format_mask_t',
'snd_pcm_subformat_mask_t', 'snd_pcm_class_t', 'SND_PCM_CLASS_GENERIC',
'SND_PCM_CLASS_MULTI', 'SND_PCM_CLASS_MODEM', 'SND_PCM_CLASS_DIGITIZER',
'SND_PCM_CLASS_LAST', 'snd_pcm_subclass_t', 'SND_PCM_SUBCLASS_GENERIC_MIX',
'SND_PCM_SUBCLASS_MULTI_MIX', 'SND_PCM_SUBCLASS_LAST', 'snd_pcm_stream_t',
'SND_PCM_STREAM_PLAYBACK', 'SND_PCM_STREAM_CAPTURE', 'SND_PCM_STREAM_LAST',
'snd_pcm_access_t', 'SND_PCM_ACCESS_MMAP_INTERLEAVED',
'SND_PCM_ACCESS_MMAP_NONINTERLEAVED', 'SND_PCM_ACCESS_MMAP_COMPLEX',
'SND_PCM_ACCESS_RW_INTERLEAVED', 'SND_PCM_ACCESS_RW_NONINTERLEAVED',
'SND_PCM_ACCESS_LAST', 'snd_pcm_format_t', 'SND_PCM_FORMAT_UNKNOWN',
'SND_PCM_FORMAT_S8', 'SND_PCM_FORMAT_U8', 'SND_PCM_FORMAT_S16_LE',
'SND_PCM_FORMAT_S16_BE', 'SND_PCM_FORMAT_U16_LE', 'SND_PCM_FORMAT_U16_BE',
'SND_PCM_FORMAT_S24_LE', 'SND_PCM_FORMAT_S24_BE', 'SND_PCM_FORMAT_U24_LE',
'SND_PCM_FORMAT_U24_BE', 'SND_PCM_FORMAT_S32_LE', 'SND_PCM_FORMAT_S32_BE',
'SND_PCM_FORMAT_U32_LE', 'SND_PCM_FORMAT_U32_BE', 'SND_PCM_FORMAT_FLOAT_LE',
'SND_PCM_FORMAT_FLOAT_BE', 'SND_PCM_FORMAT_FLOAT64_LE',
'SND_PCM_FORMAT_FLOAT64_BE', 'SND_PCM_FORMAT_IEC958_SUBFRAME_LE',
'SND_PCM_FORMAT_IEC958_SUBFRAME_BE', 'SND_PCM_FORMAT_MU_LAW',
'SND_PCM_FORMAT_A_LAW', 'SND_PCM_FORMAT_IMA_ADPCM', 'SND_PCM_FORMAT_MPEG',
'SND_PCM_FORMAT_GSM', 'SND_PCM_FORMAT_SPECIAL', 'SND_PCM_FORMAT_S24_3LE',
'SND_PCM_FORMAT_S24_3BE', 'SND_PCM_FORMAT_U24_3LE', 'SND_PCM_FORMAT_U24_3BE',
'SND_PCM_FORMAT_S20_3LE', 'SND_PCM_FORMAT_S20_3BE', 'SND_PCM_FORMAT_U20_3LE',
'SND_PCM_FORMAT_U20_3BE', 'SND_PCM_FORMAT_S18_3LE', 'SND_PCM_FORMAT_S18_3BE',
'SND_PCM_FORMAT_U18_3LE', 'SND_PCM_FORMAT_U18_3BE', 'SND_PCM_FORMAT_LAST',
'SND_PCM_FORMAT_S16', 'SND_PCM_FORMAT_U16', 'SND_PCM_FORMAT_S24',
'SND_PCM_FORMAT_U24', 'SND_PCM_FORMAT_S32', 'SND_PCM_FORMAT_U32',
'SND_PCM_FORMAT_FLOAT', 'SND_PCM_FORMAT_FLOAT64',
'SND_PCM_FORMAT_IEC958_SUBFRAME', 'snd_pcm_subformat_t',
'SND_PCM_SUBFORMAT_STD', 'SND_PCM_SUBFORMAT_LAST', 'snd_pcm_state_t',
'SND_PCM_STATE_OPEN', 'SND_PCM_STATE_SETUP', 'SND_PCM_STATE_PREPARED',
'SND_PCM_STATE_RUNNING', 'SND_PCM_STATE_XRUN', 'SND_PCM_STATE_DRAINING',
'SND_PCM_STATE_PAUSED', 'SND_PCM_STATE_SUSPENDED',
'SND_PCM_STATE_DISCONNECTED', 'SND_PCM_STATE_LAST', 'snd_pcm_start_t',
'SND_PCM_START_DATA', 'SND_PCM_START_EXPLICIT', 'SND_PCM_START_LAST',
'snd_pcm_xrun_t', 'SND_PCM_XRUN_NONE', 'SND_PCM_XRUN_STOP',
'SND_PCM_XRUN_LAST', 'snd_pcm_tstamp_t', 'SND_PCM_TSTAMP_NONE',
'SND_PCM_TSTAMP_MMAP', 'SND_PCM_TSTAMP_LAST', 'snd_pcm_uframes_t',
'snd_pcm_sframes_t', 'SND_PCM_NONBLOCK', 'SND_PCM_ASYNC', 'snd_pcm_t',
'snd_pcm_type_t', 'snd_pcm_channel_area_t', 'snd_pcm_sync_id_t',
'snd_pcm_scope_t', 'snd_pcm_open', 'snd_pcm_open_lconf', 'snd_pcm_close',
'snd_pcm_name', 'snd_pcm_type', 'snd_pcm_stream',
'snd_pcm_poll_descriptors_count', 'snd_pcm_poll_descriptors',
'snd_pcm_poll_descriptors_revents', 'snd_pcm_nonblock',
'snd_async_add_pcm_handler', 'snd_async_handler_get_pcm', 'snd_pcm_info',
'snd_pcm_hw_params_current', 'snd_pcm_hw_params', 'snd_pcm_hw_free',
'snd_pcm_sw_params_current', 'snd_pcm_sw_params', 'snd_pcm_prepare',
'snd_pcm_reset', 'snd_pcm_status', 'snd_pcm_start', 'snd_pcm_drop',
'snd_pcm_drain', 'snd_pcm_pause', 'snd_pcm_state', 'snd_pcm_hwsync',
'snd_pcm_delay', 'snd_pcm_resume', 'snd_pcm_avail_update', 'snd_pcm_rewind',
'snd_pcm_forward', 'snd_pcm_writei', 'snd_pcm_readi', 'snd_pcm_writen',
'snd_pcm_readn', 'snd_pcm_wait', 'snd_pcm_link', 'snd_pcm_unlink',
'snd_pcm_recover', 'snd_pcm_set_params', 'snd_pcm_get_params',
'snd_pcm_info_sizeof', 'snd_pcm_info_malloc', 'snd_pcm_info_free',
'snd_pcm_info_copy', 'snd_pcm_info_get_device', 'snd_pcm_info_get_subdevice',
'snd_pcm_info_get_stream', 'snd_pcm_info_get_card', 'snd_pcm_info_get_id',
'snd_pcm_info_get_name', 'snd_pcm_info_get_subdevice_name',
'snd_pcm_info_get_class', 'snd_pcm_info_get_subclass',
'snd_pcm_info_get_subdevices_count', 'snd_pcm_info_get_subdevices_avail',
'snd_pcm_info_get_sync', 'snd_pcm_info_set_device',
'snd_pcm_info_set_subdevice', 'snd_pcm_info_set_stream',
'snd_pcm_hw_params_any', 'snd_pcm_hw_params_can_mmap_sample_resolution',
'snd_pcm_hw_params_is_double', 'snd_pcm_hw_params_is_batch',
'snd_pcm_hw_params_is_block_transfer', 'snd_pcm_hw_params_can_overrange',
'snd_pcm_hw_params_can_pause', 'snd_pcm_hw_params_can_resume',
'snd_pcm_hw_params_is_half_duplex', 'snd_pcm_hw_params_is_joint_duplex',
'snd_pcm_hw_params_can_sync_start', 'snd_pcm_hw_params_get_rate_numden',
'snd_pcm_hw_params_get_sbits', 'snd_pcm_hw_params_get_fifo_size',
'snd_pcm_hw_params_sizeof', 'snd_pcm_hw_params_malloc',
'snd_pcm_hw_params_free', 'snd_pcm_hw_params_copy',
'snd_pcm_hw_params_get_access', 'snd_pcm_hw_params_test_access',
'snd_pcm_hw_params_set_access', 'snd_pcm_hw_params_set_access_first',
'snd_pcm_hw_params_set_access_last', 'snd_pcm_hw_params_set_access_mask',
'snd_pcm_hw_params_get_access_mask', 'snd_pcm_hw_params_get_format',
'snd_pcm_hw_params_test_format', 'snd_pcm_hw_params_set_format',
'snd_pcm_hw_params_set_format_first', 'snd_pcm_hw_params_set_format_last',
'snd_pcm_hw_params_set_format_mask', 'snd_pcm_hw_params_get_format_mask',
'snd_pcm_hw_params_get_subformat', 'snd_pcm_hw_params_test_subformat',
'snd_pcm_hw_params_set_subformat', 'snd_pcm_hw_params_set_subformat_first',
'snd_pcm_hw_params_set_subformat_last',
'snd_pcm_hw_params_set_subformat_mask',
'snd_pcm_hw_params_get_subformat_mask', 'snd_pcm_hw_params_get_channels',
'snd_pcm_hw_params_get_channels_min', 'snd_pcm_hw_params_get_channels_max',
'snd_pcm_hw_params_test_channels', 'snd_pcm_hw_params_set_channels',
'snd_pcm_hw_params_set_channels_min', 'snd_pcm_hw_params_set_channels_max',
'snd_pcm_hw_params_set_channels_minmax',
'snd_pcm_hw_params_set_channels_near', 'snd_pcm_hw_params_set_channels_first',
'snd_pcm_hw_params_set_channels_last', 'snd_pcm_hw_params_get_rate',
'snd_pcm_hw_params_get_rate_min', 'snd_pcm_hw_params_get_rate_max',
'snd_pcm_hw_params_test_rate', 'snd_pcm_hw_params_set_rate',
'snd_pcm_hw_params_set_rate_min', 'snd_pcm_hw_params_set_rate_max',
'snd_pcm_hw_params_set_rate_minmax', 'snd_pcm_hw_params_set_rate_near',
'snd_pcm_hw_params_set_rate_first', 'snd_pcm_hw_params_set_rate_last',
'snd_pcm_hw_params_set_rate_resample', 'snd_pcm_hw_params_get_rate_resample',
'snd_pcm_hw_params_set_export_buffer', 'snd_pcm_hw_params_get_export_buffer',
'snd_pcm_hw_params_get_period_time', 'snd_pcm_hw_params_get_period_time_min',
'snd_pcm_hw_params_get_period_time_max', 'snd_pcm_hw_params_test_period_time',
'snd_pcm_hw_params_set_period_time', 'snd_pcm_hw_params_set_period_time_min',
'snd_pcm_hw_params_set_period_time_max',
'snd_pcm_hw_params_set_period_time_minmax',
'snd_pcm_hw_params_set_period_time_near',
'snd_pcm_hw_params_set_period_time_first',
'snd_pcm_hw_params_set_period_time_last', 'snd_pcm_hw_params_get_period_size',
'snd_pcm_hw_params_get_period_size_min',
'snd_pcm_hw_params_get_period_size_max', 'snd_pcm_hw_params_test_period_size',
'snd_pcm_hw_params_set_period_size', 'snd_pcm_hw_params_set_period_size_min',
'snd_pcm_hw_params_set_period_size_max',
'snd_pcm_hw_params_set_period_size_minmax',
'snd_pcm_hw_params_set_period_size_near',
'snd_pcm_hw_params_set_period_size_first',
'snd_pcm_hw_params_set_period_size_last',
'snd_pcm_hw_params_set_period_size_integer', 'snd_pcm_hw_params_get_periods',
'snd_pcm_hw_params_get_periods_min', 'snd_pcm_hw_params_get_periods_max',
'snd_pcm_hw_params_test_periods', 'snd_pcm_hw_params_set_periods',
'snd_pcm_hw_params_set_periods_min', 'snd_pcm_hw_params_set_periods_max',
'snd_pcm_hw_params_set_periods_minmax', 'snd_pcm_hw_params_set_periods_near',
'snd_pcm_hw_params_set_periods_first', 'snd_pcm_hw_params_set_periods_last',
'snd_pcm_hw_params_set_periods_integer', 'snd_pcm_hw_params_get_buffer_time',
'snd_pcm_hw_params_get_buffer_time_min',
'snd_pcm_hw_params_get_buffer_time_max', 'snd_pcm_hw_params_test_buffer_time',
'snd_pcm_hw_params_set_buffer_time', 'snd_pcm_hw_params_set_buffer_time_min',
'snd_pcm_hw_params_set_buffer_time_max',
'snd_pcm_hw_params_set_buffer_time_minmax',
'snd_pcm_hw_params_set_buffer_time_near',
'snd_pcm_hw_params_set_buffer_time_first',
'snd_pcm_hw_params_set_buffer_time_last', 'snd_pcm_hw_params_get_buffer_size',
'snd_pcm_hw_params_get_buffer_size_min',
'snd_pcm_hw_params_get_buffer_size_max', 'snd_pcm_hw_params_test_buffer_size',
'snd_pcm_hw_params_set_buffer_size', 'snd_pcm_hw_params_set_buffer_size_min',
'snd_pcm_hw_params_set_buffer_size_max',
'snd_pcm_hw_params_set_buffer_size_minmax',
'snd_pcm_hw_params_set_buffer_size_near',
'snd_pcm_hw_params_set_buffer_size_first',
'snd_pcm_hw_params_set_buffer_size_last', 'snd_pcm_hw_params_get_tick_time',
'snd_pcm_hw_params_get_tick_time_min', 'snd_pcm_hw_params_get_tick_time_max',
'snd_pcm_hw_params_test_tick_time', 'snd_pcm_hw_params_set_tick_time',
'snd_pcm_hw_params_set_tick_time_min', 'snd_pcm_hw_params_set_tick_time_max',
'snd_pcm_hw_params_set_tick_time_minmax',
'snd_pcm_hw_params_set_tick_time_near',
'snd_pcm_hw_params_set_tick_time_first',
'snd_pcm_hw_params_set_tick_time_last', 'snd_pcm_hw_params_get_min_align',
'snd_pcm_sw_params_sizeof', 'snd_pcm_sw_params_malloc',
'snd_pcm_sw_params_free', 'snd_pcm_sw_params_copy',
'snd_pcm_sw_params_get_boundary', 'snd_pcm_sw_params_set_tstamp_mode',
'snd_pcm_sw_params_get_tstamp_mode', 'snd_pcm_sw_params_set_sleep_min',
'snd_pcm_sw_params_get_sleep_min', 'snd_pcm_sw_params_set_avail_min',
'snd_pcm_sw_params_get_avail_min', 'snd_pcm_sw_params_set_xfer_align',
'snd_pcm_sw_params_get_xfer_align', 'snd_pcm_sw_params_set_start_threshold',
'snd_pcm_sw_params_get_start_threshold',
'snd_pcm_sw_params_set_stop_threshold',
'snd_pcm_sw_params_get_stop_threshold',
'snd_pcm_sw_params_set_silence_threshold',
'snd_pcm_sw_params_get_silence_threshold',
'snd_pcm_sw_params_set_silence_size', 'snd_pcm_sw_params_get_silence_size',
'snd_pcm_access_mask_sizeof', 'snd_pcm_access_mask_malloc',
'snd_pcm_access_mask_free', 'snd_pcm_access_mask_copy',
'snd_pcm_access_mask_none', 'snd_pcm_access_mask_any',
'snd_pcm_access_mask_test', 'snd_pcm_access_mask_empty',
'snd_pcm_access_mask_set', 'snd_pcm_access_mask_reset',
'snd_pcm_format_mask_sizeof', 'snd_pcm_format_mask_malloc',
'snd_pcm_format_mask_free', 'snd_pcm_format_mask_copy',
'snd_pcm_format_mask_none', 'snd_pcm_format_mask_any',
'snd_pcm_format_mask_test', 'snd_pcm_format_mask_empty',
'snd_pcm_format_mask_set', 'snd_pcm_format_mask_reset',
'snd_pcm_subformat_mask_sizeof', 'snd_pcm_subformat_mask_malloc',
'snd_pcm_subformat_mask_free', 'snd_pcm_subformat_mask_copy',
'snd_pcm_subformat_mask_none', 'snd_pcm_subformat_mask_any',
'snd_pcm_subformat_mask_test', 'snd_pcm_subformat_mask_empty',
'snd_pcm_subformat_mask_set', 'snd_pcm_subformat_mask_reset',
'snd_pcm_status_sizeof', 'snd_pcm_status_malloc', 'snd_pcm_status_free',
'snd_pcm_status_copy', 'snd_pcm_status_get_state',
'snd_pcm_status_get_trigger_tstamp', 'snd_pcm_status_get_trigger_htstamp',
'snd_pcm_status_get_tstamp', 'snd_pcm_status_get_htstamp',
'snd_pcm_status_get_delay', 'snd_pcm_status_get_avail',
'snd_pcm_status_get_avail_max', 'snd_pcm_status_get_overrange',
'snd_pcm_type_name', 'snd_pcm_stream_name', 'snd_pcm_access_name',
'snd_pcm_format_name', 'snd_pcm_format_description', 'snd_pcm_subformat_name',
'snd_pcm_subformat_description', 'snd_pcm_format_value',
'snd_pcm_tstamp_mode_name', 'snd_pcm_state_name', 'snd_pcm_dump',
'snd_pcm_dump_hw_setup', 'snd_pcm_dump_sw_setup', 'snd_pcm_dump_setup',
'snd_pcm_hw_params_dump', 'snd_pcm_sw_params_dump', 'snd_pcm_status_dump',
'snd_pcm_mmap_begin', 'snd_pcm_mmap_commit', 'snd_pcm_mmap_writei',
'snd_pcm_mmap_readi', 'snd_pcm_mmap_writen', 'snd_pcm_mmap_readn',
'snd_pcm_format_signed', 'snd_pcm_format_unsigned', 'snd_pcm_format_linear',
'snd_pcm_format_float', 'snd_pcm_format_little_endian',
'snd_pcm_format_big_endian', 'snd_pcm_format_cpu_endian',
'snd_pcm_format_width', 'snd_pcm_format_physical_width',
'snd_pcm_build_linear_format', 'snd_pcm_format_size',
'snd_pcm_format_silence', 'snd_pcm_format_silence_16',
'snd_pcm_format_silence_32', 'snd_pcm_format_silence_64',
'snd_pcm_format_set_silence', 'snd_pcm_bytes_to_frames',
'snd_pcm_frames_to_bytes', 'snd_pcm_bytes_to_samples',
'snd_pcm_samples_to_bytes', 'snd_pcm_area_silence', 'snd_pcm_areas_silence',
'snd_pcm_area_copy', 'snd_pcm_areas_copy', 'snd_pcm_hook_type_t',
'SND_PCM_HOOK_TYPE_HW_PARAMS', 'SND_PCM_HOOK_TYPE_HW_FREE',
'SND_PCM_HOOK_TYPE_CLOSE', 'SND_PCM_HOOK_TYPE_LAST', 'snd_pcm_hook_t',
'snd_pcm_hook_func_t', 'snd_pcm_hook_get_pcm', 'snd_pcm_hook_get_private',
'snd_pcm_hook_set_private', 'snd_pcm_hook_add', 'snd_pcm_hook_remove',
'snd_pcm_scope_ops_t', 'snd_pcm_meter_get_bufsize',
'snd_pcm_meter_get_channels', 'snd_pcm_meter_get_rate',
'snd_pcm_meter_get_now', 'snd_pcm_meter_get_boundary',
'snd_pcm_meter_add_scope', 'snd_pcm_meter_search_scope',
'snd_pcm_scope_malloc', 'snd_pcm_scope_set_ops', 'snd_pcm_scope_set_name',
'snd_pcm_scope_get_name', 'snd_pcm_scope_get_callback_private',
'snd_pcm_scope_set_callback_private', 'snd_pcm_scope_s16_open',
'snd_pcm_scope_s16_get_channel_buffer', 'snd_spcm_latency_t',
'SND_SPCM_LATENCY_STANDARD', 'SND_SPCM_LATENCY_MEDIUM',
'SND_SPCM_LATENCY_REALTIME', 'snd_spcm_xrun_type_t', 'SND_SPCM_XRUN_IGNORE',
'SND_SPCM_XRUN_STOP', 'snd_spcm_duplex_type_t', 'SND_SPCM_DUPLEX_LIBERAL',
'SND_SPCM_DUPLEX_PEDANTIC', 'snd_spcm_init', 'snd_spcm_init_duplex',
'snd_spcm_init_get_params', 'snd_pcm_start_mode_name',
'snd_pcm_xrun_mode_name', 'snd_pcm_sw_params_set_start_mode',
'snd_pcm_sw_params_get_start_mode', 'snd_pcm_sw_params_set_xrun_mode',
'snd_pcm_sw_params_get_xrun_mode', 'SND_RAWMIDI_DLSYM_VERSION',
'snd_rawmidi_info_t', 'snd_rawmidi_params_t', 'snd_rawmidi_status_t',
'snd_rawmidi_stream_t', 'SND_RAWMIDI_STREAM_OUTPUT',
'SND_RAWMIDI_STREAM_INPUT', 'SND_RAWMIDI_STREAM_LAST', 'SND_RAWMIDI_APPEND',
'SND_RAWMIDI_NONBLOCK', 'SND_RAWMIDI_SYNC', 'snd_rawmidi_t',
'snd_rawmidi_type_t', 'SND_RAWMIDI_TYPE_HW', 'SND_RAWMIDI_TYPE_SHM',
'SND_RAWMIDI_TYPE_INET', 'SND_RAWMIDI_TYPE_VIRTUAL', 'snd_rawmidi_open',
'snd_rawmidi_open_lconf', 'snd_rawmidi_close',
'snd_rawmidi_poll_descriptors_count', 'snd_rawmidi_poll_descriptors',
'snd_rawmidi_poll_descriptors_revents', 'snd_rawmidi_nonblock',
'snd_rawmidi_info_sizeof', 'snd_rawmidi_info_malloc', 'snd_rawmidi_info_free',
'snd_rawmidi_info_copy', 'snd_rawmidi_info_get_device',
'snd_rawmidi_info_get_subdevice', 'snd_rawmidi_info_get_stream',
'snd_rawmidi_info_get_card', 'snd_rawmidi_info_get_flags',
'snd_rawmidi_info_get_id', 'snd_rawmidi_info_get_name',
'snd_rawmidi_info_get_subdevice_name',
'snd_rawmidi_info_get_subdevices_count',
'snd_rawmidi_info_get_subdevices_avail', 'snd_rawmidi_info_set_device',
'snd_rawmidi_info_set_subdevice', 'snd_rawmidi_info_set_stream',
'snd_rawmidi_info', 'snd_rawmidi_params_sizeof', 'snd_rawmidi_params_malloc',
'snd_rawmidi_params_free', 'snd_rawmidi_params_copy',
'snd_rawmidi_params_set_buffer_size', 'snd_rawmidi_params_get_buffer_size',
'snd_rawmidi_params_set_avail_min', 'snd_rawmidi_params_get_avail_min',
'snd_rawmidi_params_set_no_active_sensing',
'snd_rawmidi_params_get_no_active_sensing', 'snd_rawmidi_params',
'snd_rawmidi_params_current', 'snd_rawmidi_status_sizeof',
'snd_rawmidi_status_malloc', 'snd_rawmidi_status_free',
'snd_rawmidi_status_copy', 'snd_rawmidi_status_get_tstamp',
'snd_rawmidi_status_get_avail', 'snd_rawmidi_status_get_xruns',
'snd_rawmidi_status', 'snd_rawmidi_drain', 'snd_rawmidi_drop',
'snd_rawmidi_write', 'snd_rawmidi_read', 'snd_rawmidi_name',
'snd_rawmidi_type', 'snd_rawmidi_stream', 'SND_TIMER_DLSYM_VERSION',
'SND_TIMER_QUERY_DLSYM_VERSION', 'snd_timer_id_t', 'snd_timer_ginfo_t',
'snd_timer_gparams_t', 'snd_timer_gstatus_t', 'snd_timer_info_t',
'snd_timer_params_t', 'snd_timer_status_t', 'snd_timer_class_t',
'SND_TIMER_CLASS_NONE', 'SND_TIMER_CLASS_SLAVE', 'SND_TIMER_CLASS_GLOBAL',
'SND_TIMER_CLASS_CARD', 'SND_TIMER_CLASS_PCM', 'SND_TIMER_CLASS_LAST',
'snd_timer_slave_class_t', 'SND_TIMER_SCLASS_NONE',
'SND_TIMER_SCLASS_APPLICATION', 'SND_TIMER_SCLASS_SEQUENCER',
'SND_TIMER_SCLASS_OSS_SEQUENCER', 'SND_TIMER_SCLASS_LAST',
'snd_timer_event_t', 'SND_TIMER_EVENT_RESOLUTION', 'SND_TIMER_EVENT_TICK',
'SND_TIMER_EVENT_START', 'SND_TIMER_EVENT_STOP', 'SND_TIMER_EVENT_CONTINUE',
'SND_TIMER_EVENT_PAUSE', 'SND_TIMER_EVENT_EARLY', 'SND_TIMER_EVENT_SUSPEND',
'SND_TIMER_EVENT_RESUME', 'SND_TIMER_EVENT_MSTART', 'SND_TIMER_EVENT_MSTOP',
'SND_TIMER_EVENT_MCONTINUE', 'SND_TIMER_EVENT_MPAUSE',
'SND_TIMER_EVENT_MSUSPEND', 'SND_TIMER_EVENT_MRESUME', 'snd_timer_read_t',
'snd_timer_tread_t', 'SND_TIMER_GLOBAL_SYSTEM', 'SND_TIMER_GLOBAL_RTC',
'SND_TIMER_GLOBAL_HPET', 'SND_TIMER_OPEN_NONBLOCK', 'SND_TIMER_OPEN_TREAD',
'snd_timer_type_t', 'SND_TIMER_TYPE_HW', 'SND_TIMER_TYPE_SHM',
'SND_TIMER_TYPE_INET', 'snd_timer_query_t', 'snd_timer_t',
'snd_timer_query_open', 'snd_timer_query_open_lconf', 'snd_timer_query_close',
'snd_timer_query_next_device', 'snd_timer_query_info',
'snd_timer_query_params', 'snd_timer_query_status', 'snd_timer_open',
'snd_timer_open_lconf', 'snd_timer_close', 'snd_async_add_timer_handler',
'snd_async_handler_get_timer', 'snd_timer_poll_descriptors_count',
'snd_timer_poll_descriptors', 'snd_timer_poll_descriptors_revents',
'snd_timer_info', 'snd_timer_params', 'snd_timer_status', 'snd_timer_start',
'snd_timer_stop', 'snd_timer_continue', 'snd_timer_read',
'snd_timer_id_sizeof', 'snd_timer_id_malloc', 'snd_timer_id_free',
'snd_timer_id_copy', 'snd_timer_id_set_class', 'snd_timer_id_get_class',
'snd_timer_id_set_sclass', 'snd_timer_id_get_sclass', 'snd_timer_id_set_card',
'snd_timer_id_get_card', 'snd_timer_id_set_device', 'snd_timer_id_get_device',
'snd_timer_id_set_subdevice', 'snd_timer_id_get_subdevice',
'snd_timer_ginfo_sizeof', 'snd_timer_ginfo_malloc', 'snd_timer_ginfo_free',
'snd_timer_ginfo_copy', 'snd_timer_ginfo_set_tid', 'snd_timer_ginfo_get_tid',
'snd_timer_ginfo_get_flags', 'snd_timer_ginfo_get_card',
'snd_timer_ginfo_get_id', 'snd_timer_ginfo_get_name',
'snd_timer_ginfo_get_resolution', 'snd_timer_ginfo_get_resolution_min',
'snd_timer_ginfo_get_resolution_max', 'snd_timer_ginfo_get_clients',
'snd_timer_info_sizeof', 'snd_timer_info_malloc', 'snd_timer_info_free',
'snd_timer_info_copy', 'snd_timer_info_is_slave', 'snd_timer_info_get_card',
'snd_timer_info_get_id', 'snd_timer_info_get_name',
'snd_timer_info_get_resolution', 'snd_timer_params_sizeof',
'snd_timer_params_malloc', 'snd_timer_params_free', 'snd_timer_params_copy',
'snd_timer_params_set_auto_start', 'snd_timer_params_get_auto_start',
'snd_timer_params_set_exclusive', 'snd_timer_params_get_exclusive',
'snd_timer_params_set_early_event', 'snd_timer_params_get_early_event',
'snd_timer_params_set_ticks', 'snd_timer_params_get_ticks',
'snd_timer_params_set_queue_size', 'snd_timer_params_get_queue_size',
'snd_timer_params_set_filter', 'snd_timer_params_get_filter',
'snd_timer_status_sizeof', 'snd_timer_status_malloc', 'snd_timer_status_free',
'snd_timer_status_copy', 'snd_timer_status_get_timestamp',
'snd_timer_status_get_resolution', 'snd_timer_status_get_lost',
'snd_timer_status_get_overrun', 'snd_timer_status_get_queue',
'snd_timer_info_get_ticks', 'SND_HWDEP_DLSYM_VERSION', 'snd_hwdep_info_t',
'snd_hwdep_dsp_status_t', 'snd_hwdep_dsp_image_t', 'snd_hwdep_iface_t',
'SND_HWDEP_IFACE_OPL2', 'SND_HWDEP_IFACE_OPL3', 'SND_HWDEP_IFACE_OPL4',
'SND_HWDEP_IFACE_SB16CSP', 'SND_HWDEP_IFACE_EMU10K1',
'SND_HWDEP_IFACE_YSS225', 'SND_HWDEP_IFACE_ICS2115', 'SND_HWDEP_IFACE_SSCAPE',
'SND_HWDEP_IFACE_VX', 'SND_HWDEP_IFACE_MIXART', 'SND_HWDEP_IFACE_USX2Y',
'SND_HWDEP_IFACE_EMUX_WAVETABLE', 'SND_HWDEP_IFACE_BLUETOOTH',
'SND_HWDEP_IFACE_USX2Y_PCM', 'SND_HWDEP_IFACE_PCXHR', 'SND_HWDEP_IFACE_SB_RC',
'SND_HWDEP_IFACE_LAST', 'SND_HWDEP_OPEN_READ', 'SND_HWDEP_OPEN_WRITE',
'SND_HWDEP_OPEN_DUPLEX', 'SND_HWDEP_OPEN_NONBLOCK', 'snd_hwdep_type_t',
'SND_HWDEP_TYPE_HW', 'SND_HWDEP_TYPE_SHM', 'SND_HWDEP_TYPE_INET',
'snd_hwdep_t', 'snd_hwdep_open', 'snd_hwdep_close',
'snd_hwdep_poll_descriptors', 'snd_hwdep_poll_descriptors_revents',
'snd_hwdep_nonblock', 'snd_hwdep_info', 'snd_hwdep_dsp_status',
'snd_hwdep_dsp_load', 'snd_hwdep_ioctl', 'snd_hwdep_write', 'snd_hwdep_read',
'snd_hwdep_info_sizeof', 'snd_hwdep_info_malloc', 'snd_hwdep_info_free',
'snd_hwdep_info_copy', 'snd_hwdep_info_get_device', 'snd_hwdep_info_get_card',
'snd_hwdep_info_get_id', 'snd_hwdep_info_get_name',
'snd_hwdep_info_get_iface', 'snd_hwdep_info_set_device',
'snd_hwdep_dsp_status_sizeof', 'snd_hwdep_dsp_status_malloc',
'snd_hwdep_dsp_status_free', 'snd_hwdep_dsp_status_copy',
'snd_hwdep_dsp_status_get_version', 'snd_hwdep_dsp_status_get_id',
'snd_hwdep_dsp_status_get_num_dsps', 'snd_hwdep_dsp_status_get_dsp_loaded',
'snd_hwdep_dsp_status_get_chip_ready', 'snd_hwdep_dsp_image_sizeof',
'snd_hwdep_dsp_image_malloc', 'snd_hwdep_dsp_image_free',
'snd_hwdep_dsp_image_copy', 'snd_hwdep_dsp_image_get_index',
'snd_hwdep_dsp_image_get_name', 'snd_hwdep_dsp_image_get_image',
'snd_hwdep_dsp_image_get_length', 'snd_hwdep_dsp_image_set_index',
'snd_hwdep_dsp_image_set_name', 'snd_hwdep_dsp_image_set_image',
'snd_hwdep_dsp_image_set_length', 'SND_CONTROL_DLSYM_VERSION',
'snd_aes_iec958_t', 'snd_ctl_card_info_t', 'snd_ctl_elem_id_t',
'snd_ctl_elem_list_t', 'snd_ctl_elem_info_t', 'snd_ctl_elem_value_t',
'snd_ctl_event_t', 'snd_ctl_elem_type_t', 'SND_CTL_ELEM_TYPE_NONE',
'SND_CTL_ELEM_TYPE_BOOLEAN', 'SND_CTL_ELEM_TYPE_INTEGER',
'SND_CTL_ELEM_TYPE_ENUMERATED', 'SND_CTL_ELEM_TYPE_BYTES',
'SND_CTL_ELEM_TYPE_IEC958', 'SND_CTL_ELEM_TYPE_INTEGER64',
'SND_CTL_ELEM_TYPE_LAST', 'snd_ctl_elem_iface_t', 'SND_CTL_ELEM_IFACE_CARD',
'SND_CTL_ELEM_IFACE_HWDEP', 'SND_CTL_ELEM_IFACE_MIXER',
'SND_CTL_ELEM_IFACE_PCM', 'SND_CTL_ELEM_IFACE_RAWMIDI',
'SND_CTL_ELEM_IFACE_TIMER', 'SND_CTL_ELEM_IFACE_SEQUENCER',
'SND_CTL_ELEM_IFACE_LAST', 'snd_ctl_event_type_t', 'SND_CTL_EVENT_ELEM',
'SND_CTL_EVENT_LAST', 'SND_CTL_EVENT_MASK_REMOVE', 'SND_CTL_EVENT_MASK_VALUE',
'SND_CTL_EVENT_MASK_INFO', 'SND_CTL_EVENT_MASK_ADD', 'SND_CTL_EVENT_MASK_TLV',
'SND_CTL_POWER_MASK', 'SND_CTL_POWER_D0', 'SND_CTL_POWER_D1',
'SND_CTL_POWER_D2', 'SND_CTL_POWER_D3', 'SND_CTL_POWER_D3hot',
'SND_CTL_POWER_D3cold', 'SND_CTL_TLVT_CONTAINER', 'SND_CTL_TLVT_DB_SCALE',
'SND_CTL_TLVT_DB_LINEAR', 'SND_CTL_TLVT_DB_RANGE', 'SND_CTL_TLV_DB_GAIN_MUTE',
'snd_ctl_type_t', 'SND_CTL_TYPE_HW', 'SND_CTL_TYPE_SHM', 'SND_CTL_TYPE_INET',
'SND_CTL_TYPE_EXT', 'SND_CTL_NONBLOCK', 'SND_CTL_ASYNC', 'SND_CTL_READONLY',
'snd_ctl_t', 'SND_SCTL_NOFREE', 'snd_sctl_t', 'snd_card_load',
'snd_card_next', 'snd_card_get_index', 'snd_card_get_name',
'snd_card_get_longname', 'snd_device_name_hint', 'snd_device_name_free_hint',
'snd_device_name_get_hint', 'snd_ctl_open', 'snd_ctl_open_lconf',
'snd_ctl_close', 'snd_ctl_nonblock', 'snd_async_add_ctl_handler',
'snd_async_handler_get_ctl', 'snd_ctl_poll_descriptors_count',
'snd_ctl_poll_descriptors', 'snd_ctl_poll_descriptors_revents',
'snd_ctl_subscribe_events', 'snd_ctl_card_info', 'snd_ctl_elem_list',
'snd_ctl_elem_info', 'snd_ctl_elem_read', 'snd_ctl_elem_write',
'snd_ctl_elem_lock', 'snd_ctl_elem_unlock', 'snd_ctl_elem_tlv_read',
'snd_ctl_elem_tlv_write', 'snd_ctl_elem_tlv_command',
'snd_ctl_hwdep_next_device', 'snd_ctl_hwdep_info', 'snd_ctl_pcm_next_device',
'snd_ctl_pcm_info', 'snd_ctl_pcm_prefer_subdevice',
'snd_ctl_rawmidi_next_device', 'snd_ctl_rawmidi_info',
'snd_ctl_rawmidi_prefer_subdevice', 'snd_ctl_set_power_state',
'snd_ctl_get_power_state', 'snd_ctl_read', 'snd_ctl_wait', 'snd_ctl_name',
'snd_ctl_type', 'snd_ctl_elem_type_name', 'snd_ctl_elem_iface_name',
'snd_ctl_event_type_name', 'snd_ctl_event_elem_get_mask',
'snd_ctl_event_elem_get_numid', 'snd_ctl_event_elem_get_id',
'snd_ctl_event_elem_get_interface', 'snd_ctl_event_elem_get_device',
'snd_ctl_event_elem_get_subdevice', 'snd_ctl_event_elem_get_name',
'snd_ctl_event_elem_get_index', 'snd_ctl_elem_list_alloc_space',
'snd_ctl_elem_list_free_space', 'snd_ctl_elem_id_sizeof',
'snd_ctl_elem_id_malloc', 'snd_ctl_elem_id_free', 'snd_ctl_elem_id_clear',
'snd_ctl_elem_id_copy', 'snd_ctl_elem_id_get_numid',
'snd_ctl_elem_id_get_interface', 'snd_ctl_elem_id_get_device',
'snd_ctl_elem_id_get_subdevice', 'snd_ctl_elem_id_get_name',
'snd_ctl_elem_id_get_index', 'snd_ctl_elem_id_set_numid',
'snd_ctl_elem_id_set_interface', 'snd_ctl_elem_id_set_device',
'snd_ctl_elem_id_set_subdevice', 'snd_ctl_elem_id_set_name',
'snd_ctl_elem_id_set_index', 'snd_ctl_card_info_sizeof',
'snd_ctl_card_info_malloc', 'snd_ctl_card_info_free',
'snd_ctl_card_info_clear', 'snd_ctl_card_info_copy',
'snd_ctl_card_info_get_card', 'snd_ctl_card_info_get_id',
'snd_ctl_card_info_get_driver', 'snd_ctl_card_info_get_name',
'snd_ctl_card_info_get_longname', 'snd_ctl_card_info_get_mixername',
'snd_ctl_card_info_get_components', 'snd_ctl_event_sizeof',
'snd_ctl_event_malloc', 'snd_ctl_event_free', 'snd_ctl_event_clear',
'snd_ctl_event_copy', 'snd_ctl_event_get_type', 'snd_ctl_elem_list_sizeof',
'snd_ctl_elem_list_malloc', 'snd_ctl_elem_list_free',
'snd_ctl_elem_list_clear', 'snd_ctl_elem_list_copy',
'snd_ctl_elem_list_set_offset', 'snd_ctl_elem_list_get_used',
'snd_ctl_elem_list_get_count', 'snd_ctl_elem_list_get_id',
'snd_ctl_elem_list_get_numid', 'snd_ctl_elem_list_get_interface',
'snd_ctl_elem_list_get_device', 'snd_ctl_elem_list_get_subdevice',
'snd_ctl_elem_list_get_name', 'snd_ctl_elem_list_get_index',
'snd_ctl_elem_info_sizeof', 'snd_ctl_elem_info_malloc',
'snd_ctl_elem_info_free', 'snd_ctl_elem_info_clear', 'snd_ctl_elem_info_copy',
'snd_ctl_elem_info_get_type', 'snd_ctl_elem_info_is_readable',
'snd_ctl_elem_info_is_writable', 'snd_ctl_elem_info_is_volatile',
'snd_ctl_elem_info_is_inactive', 'snd_ctl_elem_info_is_locked',
'snd_ctl_elem_info_is_tlv_readable', 'snd_ctl_elem_info_is_tlv_writable',
'snd_ctl_elem_info_is_tlv_commandable', 'snd_ctl_elem_info_is_owner',
'snd_ctl_elem_info_is_user', 'snd_ctl_elem_info_get_owner',
'snd_ctl_elem_info_get_count', 'snd_ctl_elem_info_get_min',
'snd_ctl_elem_info_get_max', 'snd_ctl_elem_info_get_step',
'snd_ctl_elem_info_get_min64', 'snd_ctl_elem_info_get_max64',
'snd_ctl_elem_info_get_step64', 'snd_ctl_elem_info_get_items',
'snd_ctl_elem_info_set_item', 'snd_ctl_elem_info_get_item_name',
'snd_ctl_elem_info_get_dimensions', 'snd_ctl_elem_info_get_dimension',
'snd_ctl_elem_info_get_id', 'snd_ctl_elem_info_get_numid',
'snd_ctl_elem_info_get_interface', 'snd_ctl_elem_info_get_device',
'snd_ctl_elem_info_get_subdevice', 'snd_ctl_elem_info_get_name',
'snd_ctl_elem_info_get_index', 'snd_ctl_elem_info_set_id',
'snd_ctl_elem_info_set_numid', 'snd_ctl_elem_info_set_interface',
'snd_ctl_elem_info_set_device', 'snd_ctl_elem_info_set_subdevice',
'snd_ctl_elem_info_set_name', 'snd_ctl_elem_info_set_index',
'snd_ctl_elem_add_integer', 'snd_ctl_elem_add_integer64',
'snd_ctl_elem_add_boolean', 'snd_ctl_elem_add_iec958', 'snd_ctl_elem_remove',
'snd_ctl_elem_value_sizeof', 'snd_ctl_elem_value_malloc',
'snd_ctl_elem_value_free', 'snd_ctl_elem_value_clear',
'snd_ctl_elem_value_copy', 'snd_ctl_elem_value_get_id',
'snd_ctl_elem_value_get_numid', 'snd_ctl_elem_value_get_interface',
'snd_ctl_elem_value_get_device', 'snd_ctl_elem_value_get_subdevice',
'snd_ctl_elem_value_get_name', 'snd_ctl_elem_value_get_index',
'snd_ctl_elem_value_set_id', 'snd_ctl_elem_value_set_numid',
'snd_ctl_elem_value_set_interface', 'snd_ctl_elem_value_set_device',
'snd_ctl_elem_value_set_subdevice', 'snd_ctl_elem_value_set_name',
'snd_ctl_elem_value_set_index', 'snd_ctl_elem_value_get_boolean',
'snd_ctl_elem_value_get_integer', 'snd_ctl_elem_value_get_integer64',
'snd_ctl_elem_value_get_enumerated', 'snd_ctl_elem_value_get_byte',
'snd_ctl_elem_value_set_boolean', 'snd_ctl_elem_value_set_integer',
'snd_ctl_elem_value_set_integer64', 'snd_ctl_elem_value_set_enumerated',
'snd_ctl_elem_value_set_byte', 'snd_ctl_elem_set_bytes',
'snd_ctl_elem_value_get_bytes', 'snd_ctl_elem_value_get_iec958',
'snd_ctl_elem_value_set_iec958', 'snd_hctl_elem_t', 'snd_hctl_t',
'snd_hctl_compare_t', 'snd_hctl_compare_fast', 'snd_hctl_callback_t',
'snd_hctl_elem_callback_t', 'snd_hctl_open', 'snd_hctl_open_ctl',
'snd_hctl_close', 'snd_hctl_nonblock', 'snd_hctl_poll_descriptors_count',
'snd_hctl_poll_descriptors', 'snd_hctl_poll_descriptors_revents',
'snd_hctl_get_count', 'snd_hctl_set_compare', 'snd_hctl_first_elem',
'snd_hctl_last_elem', 'snd_hctl_find_elem', 'snd_hctl_set_callback',
'snd_hctl_set_callback_private', 'snd_hctl_get_callback_private',
'snd_hctl_load', 'snd_hctl_free', 'snd_hctl_handle_events', 'snd_hctl_name',
'snd_hctl_wait', 'snd_hctl_ctl', 'snd_hctl_elem_next', 'snd_hctl_elem_prev',
'snd_hctl_elem_info', 'snd_hctl_elem_read', 'snd_hctl_elem_write',
'snd_hctl_elem_tlv_read', 'snd_hctl_elem_tlv_write',
'snd_hctl_elem_tlv_command', 'snd_hctl_elem_get_hctl', 'snd_hctl_elem_get_id',
'snd_hctl_elem_get_numid', 'snd_hctl_elem_get_interface',
'snd_hctl_elem_get_device', 'snd_hctl_elem_get_subdevice',
'snd_hctl_elem_get_name', 'snd_hctl_elem_get_index',
'snd_hctl_elem_set_callback', 'snd_hctl_elem_get_callback_private',
'snd_hctl_elem_set_callback_private', 'snd_sctl_build', 'snd_sctl_free',
'snd_sctl_install', 'snd_sctl_remove', 'snd_mixer_t', 'snd_mixer_class_t',
'snd_mixer_elem_t', 'snd_mixer_callback_t', 'snd_mixer_elem_callback_t',
'snd_mixer_compare_t', 'snd_mixer_event_t', 'snd_mixer_elem_type_t',
'SND_MIXER_ELEM_SIMPLE', 'SND_MIXER_ELEM_LAST', 'snd_mixer_open',
'snd_mixer_close', 'snd_mixer_first_elem', 'snd_mixer_last_elem',
'snd_mixer_handle_events', 'snd_mixer_attach', 'snd_mixer_attach_hctl',
'snd_mixer_detach', 'snd_mixer_detach_hctl', 'snd_mixer_get_hctl',
'snd_mixer_poll_descriptors_count', 'snd_mixer_poll_descriptors',
'snd_mixer_poll_descriptors_revents', 'snd_mixer_load', 'snd_mixer_free',
'snd_mixer_wait', 'snd_mixer_set_compare', 'snd_mixer_set_callback',
'snd_mixer_get_callback_private', 'snd_mixer_set_callback_private',
'snd_mixer_get_count', 'snd_mixer_class_unregister', 'snd_mixer_elem_next',
'snd_mixer_elem_prev', 'snd_mixer_elem_set_callback',
'snd_mixer_elem_get_callback_private', 'snd_mixer_elem_set_callback_private',
'snd_mixer_elem_get_type', 'snd_mixer_class_register', 'snd_mixer_add_elem',
'snd_mixer_remove_elem', 'snd_mixer_elem_new', 'snd_mixer_elem_add',
'snd_mixer_elem_remove', 'snd_mixer_elem_free', 'snd_mixer_elem_info',
'snd_mixer_elem_value', 'snd_mixer_elem_attach', 'snd_mixer_elem_detach',
'snd_mixer_elem_empty', 'snd_mixer_elem_get_private',
'snd_mixer_class_sizeof', 'snd_mixer_class_malloc', 'snd_mixer_class_free',
'snd_mixer_class_copy', 'snd_mixer_class_get_mixer',
'snd_mixer_class_get_event', 'snd_mixer_class_get_private',
'snd_mixer_class_get_compare', 'snd_mixer_class_set_event',
'snd_mixer_class_set_private', 'snd_mixer_class_set_private_free',
'snd_mixer_class_set_compare', 'snd_mixer_selem_channel_id_t',
'SND_MIXER_SCHN_UNKNOWN', 'SND_MIXER_SCHN_FRONT_LEFT',
'SND_MIXER_SCHN_FRONT_RIGHT', 'SND_MIXER_SCHN_REAR_LEFT',
'SND_MIXER_SCHN_REAR_RIGHT', 'SND_MIXER_SCHN_FRONT_CENTER',
'SND_MIXER_SCHN_WOOFER', 'SND_MIXER_SCHN_SIDE_LEFT',
'SND_MIXER_SCHN_SIDE_RIGHT', 'SND_MIXER_SCHN_REAR_CENTER',
'SND_MIXER_SCHN_LAST', 'SND_MIXER_SCHN_MONO', 'snd_mixer_selem_id_t',
'snd_mixer_selem_channel_name', 'snd_mixer_selem_register',
'snd_mixer_selem_get_id', 'snd_mixer_selem_get_name',
'snd_mixer_selem_get_index', 'snd_mixer_find_selem',
'snd_mixer_selem_is_active', 'snd_mixer_selem_is_playback_mono',
'snd_mixer_selem_has_playback_channel', 'snd_mixer_selem_is_capture_mono',
'snd_mixer_selem_has_capture_channel', 'snd_mixer_selem_get_capture_group',
'snd_mixer_selem_has_common_volume', 'snd_mixer_selem_has_playback_volume',
'snd_mixer_selem_has_playback_volume_joined',
'snd_mixer_selem_has_capture_volume',
'snd_mixer_selem_has_capture_volume_joined',
'snd_mixer_selem_has_common_switch', 'snd_mixer_selem_has_playback_switch',
'snd_mixer_selem_has_playback_switch_joined',
'snd_mixer_selem_has_capture_switch',
'snd_mixer_selem_has_capture_switch_joined',
'snd_mixer_selem_has_capture_switch_exclusive',
'snd_mixer_selem_get_playback_volume', 'snd_mixer_selem_get_capture_volume',
'snd_mixer_selem_get_playback_dB', 'snd_mixer_selem_get_capture_dB',
'snd_mixer_selem_get_playback_switch', 'snd_mixer_selem_get_capture_switch',
'snd_mixer_selem_set_playback_volume', 'snd_mixer_selem_set_capture_volume',
'snd_mixer_selem_set_playback_dB', 'snd_mixer_selem_set_capture_dB',
'snd_mixer_selem_set_playback_volume_all',
'snd_mixer_selem_set_capture_volume_all',
'snd_mixer_selem_set_playback_dB_all', 'snd_mixer_selem_set_capture_dB_all',
'snd_mixer_selem_set_playback_switch', 'snd_mixer_selem_set_capture_switch',
'snd_mixer_selem_set_playback_switch_all',
'snd_mixer_selem_set_capture_switch_all',
'snd_mixer_selem_get_playback_volume_range',
'snd_mixer_selem_get_playback_dB_range',
'snd_mixer_selem_set_playback_volume_range',
'snd_mixer_selem_get_capture_volume_range',
'snd_mixer_selem_get_capture_dB_range',
'snd_mixer_selem_set_capture_volume_range', 'snd_mixer_selem_is_enumerated',
'snd_mixer_selem_is_enum_playback', 'snd_mixer_selem_is_enum_capture',
'snd_mixer_selem_get_enum_items', 'snd_mixer_selem_get_enum_item_name',
'snd_mixer_selem_get_enum_item', 'snd_mixer_selem_set_enum_item',
'snd_mixer_selem_id_sizeof', 'snd_mixer_selem_id_malloc',
'snd_mixer_selem_id_free', 'snd_mixer_selem_id_copy',
'snd_mixer_selem_id_get_name', 'snd_mixer_selem_id_get_index',
'snd_mixer_selem_id_set_name', 'snd_mixer_selem_id_set_index',
'snd_seq_event_type_t', 'snd_seq_addr_t', 'snd_seq_connect_t',
'snd_seq_real_time_t', 'snd_seq_tick_time_t', 'snd_seq_timestamp_t',
'SND_SEQ_TIME_STAMP_TICK', 'SND_SEQ_TIME_STAMP_REAL',
'SND_SEQ_TIME_STAMP_MASK', 'SND_SEQ_TIME_MODE_ABS', 'SND_SEQ_TIME_MODE_REL',
'SND_SEQ_TIME_MODE_MASK', 'SND_SEQ_EVENT_LENGTH_FIXED',
'SND_SEQ_EVENT_LENGTH_VARIABLE', 'SND_SEQ_EVENT_LENGTH_VARUSR',
'SND_SEQ_EVENT_LENGTH_MASK', 'SND_SEQ_PRIORITY_NORMAL',
'SND_SEQ_PRIORITY_HIGH', 'SND_SEQ_PRIORITY_MASK', 'snd_seq_ev_note_t',
'snd_seq_ev_ctrl_t', 'snd_seq_ev_raw8_t', 'snd_seq_ev_raw32_t',
'snd_seq_ev_ext_t', 'snd_seq_instr_cluster_t', 'snd_seq_instr_t',
'snd_seq_ev_sample_t', 'snd_seq_ev_cluster_t', 'snd_seq_position_t',
'snd_seq_stop_mode_t', 'SND_SEQ_SAMPLE_STOP_IMMEDIATELY',
'SND_SEQ_SAMPLE_STOP_VENVELOPE', 'SND_SEQ_SAMPLE_STOP_LOOP',
'snd_seq_frequency_t', 'snd_seq_ev_volume_t', 'snd_seq_ev_loop_t',
'snd_seq_ev_sample_control_t', 'snd_seq_ev_instr_begin_t', 'snd_seq_result_t',
'snd_seq_queue_skew_t', 'snd_seq_ev_queue_control_t', 'snd_seq_event_t',
'SND_SEQ_DLSYM_VERSION', 'snd_seq_t', 'SND_SEQ_OPEN_OUTPUT',
'SND_SEQ_OPEN_INPUT', 'SND_SEQ_OPEN_DUPLEX', 'SND_SEQ_NONBLOCK',
'snd_seq_type_t', 'SND_SEQ_TYPE_HW', 'SND_SEQ_TYPE_SHM', 'SND_SEQ_TYPE_INET',
'SND_SEQ_ADDRESS_UNKNOWN', 'SND_SEQ_ADDRESS_SUBSCRIBERS',
'SND_SEQ_ADDRESS_BROADCAST', 'SND_SEQ_CLIENT_SYSTEM', 'snd_seq_open',
'snd_seq_open_lconf', 'snd_seq_name', 'snd_seq_type', 'snd_seq_close',
'snd_seq_poll_descriptors_count', 'snd_seq_poll_descriptors',
'snd_seq_poll_descriptors_revents', 'snd_seq_nonblock', 'snd_seq_client_id',
'snd_seq_get_output_buffer_size', 'snd_seq_get_input_buffer_size',
'snd_seq_set_output_buffer_size', 'snd_seq_set_input_buffer_size',
'snd_seq_system_info_t', 'snd_seq_system_info_sizeof',
'snd_seq_system_info_malloc', 'snd_seq_system_info_free',
'snd_seq_system_info_copy', 'snd_seq_system_info_get_queues',
'snd_seq_system_info_get_clients', 'snd_seq_system_info_get_ports',
'snd_seq_system_info_get_channels', 'snd_seq_system_info_get_cur_clients',
'snd_seq_system_info_get_cur_queues', 'snd_seq_system_info',
'snd_seq_client_info_t', 'snd_seq_client_type_t', 'SND_SEQ_USER_CLIENT',
'SND_SEQ_KERNEL_CLIENT', 'snd_seq_client_info_sizeof',
'snd_seq_client_info_malloc', 'snd_seq_client_info_free',
'snd_seq_client_info_copy', 'snd_seq_client_info_get_client',
'snd_seq_client_info_get_type', 'snd_seq_client_info_get_name',
'snd_seq_client_info_get_broadcast_filter',
'snd_seq_client_info_get_error_bounce',
'snd_seq_client_info_get_event_filter', 'snd_seq_client_info_get_num_ports',
'snd_seq_client_info_get_event_lost', 'snd_seq_client_info_set_client',
'snd_seq_client_info_set_name', 'snd_seq_client_info_set_broadcast_filter',
'snd_seq_client_info_set_error_bounce',
'snd_seq_client_info_set_event_filter', 'snd_seq_get_client_info',
'snd_seq_get_any_client_info', 'snd_seq_set_client_info',
'snd_seq_query_next_client', 'snd_seq_client_pool_t',
'snd_seq_client_pool_sizeof', 'snd_seq_client_pool_malloc',
'snd_seq_client_pool_free', 'snd_seq_client_pool_copy',
'snd_seq_client_pool_get_client', 'snd_seq_client_pool_get_output_pool',
'snd_seq_client_pool_get_input_pool', 'snd_seq_client_pool_get_output_room',
'snd_seq_client_pool_get_output_free', 'snd_seq_client_pool_get_input_free',
'snd_seq_client_pool_set_output_pool', 'snd_seq_client_pool_set_input_pool',
'snd_seq_client_pool_set_output_room', 'snd_seq_get_client_pool',
'snd_seq_set_client_pool', 'snd_seq_port_info_t', 'SND_SEQ_PORT_SYSTEM_TIMER',
'SND_SEQ_PORT_SYSTEM_ANNOUNCE', 'SND_SEQ_PORT_CAP_READ',
'SND_SEQ_PORT_CAP_WRITE', 'SND_SEQ_PORT_CAP_SYNC_READ',
'SND_SEQ_PORT_CAP_SYNC_WRITE', 'SND_SEQ_PORT_CAP_DUPLEX',
'SND_SEQ_PORT_CAP_SUBS_READ', 'SND_SEQ_PORT_CAP_SUBS_WRITE',
'SND_SEQ_PORT_CAP_NO_EXPORT', 'SND_SEQ_PORT_TYPE_SPECIFIC',
'SND_SEQ_PORT_TYPE_MIDI_GENERIC', 'SND_SEQ_PORT_TYPE_MIDI_GM',
'SND_SEQ_PORT_TYPE_MIDI_GS', 'SND_SEQ_PORT_TYPE_MIDI_XG',
'SND_SEQ_PORT_TYPE_MIDI_MT32', 'SND_SEQ_PORT_TYPE_MIDI_GM2',
'SND_SEQ_PORT_TYPE_SYNTH', 'SND_SEQ_PORT_TYPE_DIRECT_SAMPLE',
'SND_SEQ_PORT_TYPE_SAMPLE', 'SND_SEQ_PORT_TYPE_HARDWARE',
'SND_SEQ_PORT_TYPE_SOFTWARE', 'SND_SEQ_PORT_TYPE_SYNTHESIZER',
'SND_SEQ_PORT_TYPE_PORT', 'SND_SEQ_PORT_TYPE_APPLICATION',
'snd_seq_port_info_sizeof', 'snd_seq_port_info_malloc',
'snd_seq_port_info_free', 'snd_seq_port_info_copy',
'snd_seq_port_info_get_client', 'snd_seq_port_info_get_port',
'snd_seq_port_info_get_addr', 'snd_seq_port_info_get_name',
'snd_seq_port_info_get_capability', 'snd_seq_port_info_get_type',
'snd_seq_port_info_get_midi_channels', 'snd_seq_port_info_get_midi_voices',
'snd_seq_port_info_get_synth_voices', 'snd_seq_port_info_get_read_use',
'snd_seq_port_info_get_write_use', 'snd_seq_port_info_get_port_specified',
'snd_seq_port_info_get_timestamping', 'snd_seq_port_info_get_timestamp_real',
'snd_seq_port_info_get_timestamp_queue', 'snd_seq_port_info_set_client',
'snd_seq_port_info_set_port', 'snd_seq_port_info_set_addr',
'snd_seq_port_info_set_name', 'snd_seq_port_info_set_capability',
'snd_seq_port_info_set_type', 'snd_seq_port_info_set_midi_channels',
'snd_seq_port_info_set_midi_voices', 'snd_seq_port_info_set_synth_voices',
'snd_seq_port_info_set_port_specified', 'snd_seq_port_info_set_timestamping',
'snd_seq_port_info_set_timestamp_real',
'snd_seq_port_info_set_timestamp_queue', 'snd_seq_create_port',
'snd_seq_delete_port', 'snd_seq_get_port_info', 'snd_seq_get_any_port_info',
'snd_seq_set_port_info', 'snd_seq_query_next_port',
'snd_seq_port_subscribe_t', 'snd_seq_port_subscribe_sizeof',
'snd_seq_port_subscribe_malloc', 'snd_seq_port_subscribe_free',
'snd_seq_port_subscribe_copy', 'snd_seq_port_subscribe_get_sender',
'snd_seq_port_subscribe_get_dest', 'snd_seq_port_subscribe_get_queue',
'snd_seq_port_subscribe_get_exclusive',
'snd_seq_port_subscribe_get_time_update',
'snd_seq_port_subscribe_get_time_real', 'snd_seq_port_subscribe_set_sender',
'snd_seq_port_subscribe_set_dest', 'snd_seq_port_subscribe_set_queue',
'snd_seq_port_subscribe_set_exclusive',
'snd_seq_port_subscribe_set_time_update',
'snd_seq_port_subscribe_set_time_real', 'snd_seq_get_port_subscription',
'snd_seq_subscribe_port', 'snd_seq_unsubscribe_port',
'snd_seq_query_subscribe_t', 'snd_seq_query_subs_type_t',
'SND_SEQ_QUERY_SUBS_READ', 'SND_SEQ_QUERY_SUBS_WRITE',
'snd_seq_query_subscribe_sizeof', 'snd_seq_query_subscribe_malloc',
'snd_seq_query_subscribe_free', 'snd_seq_query_subscribe_copy',
'snd_seq_query_subscribe_get_client', 'snd_seq_query_subscribe_get_port',
'snd_seq_query_subscribe_get_root', 'snd_seq_query_subscribe_get_type',
'snd_seq_query_subscribe_get_index', 'snd_seq_query_subscribe_get_num_subs',
'snd_seq_query_subscribe_get_addr', 'snd_seq_query_subscribe_get_queue',
'snd_seq_query_subscribe_get_exclusive',
'snd_seq_query_subscribe_get_time_update',
'snd_seq_query_subscribe_get_time_real', 'snd_seq_query_subscribe_set_client',
'snd_seq_query_subscribe_set_port', 'snd_seq_query_subscribe_set_root',
'snd_seq_query_subscribe_set_type', 'snd_seq_query_subscribe_set_index',
'snd_seq_query_port_subscribers', 'snd_seq_queue_info_t',
'snd_seq_queue_status_t', 'snd_seq_queue_tempo_t', 'snd_seq_queue_timer_t',
'SND_SEQ_QUEUE_DIRECT', 'snd_seq_queue_info_sizeof',
'snd_seq_queue_info_malloc', 'snd_seq_queue_info_free',
'snd_seq_queue_info_copy', 'snd_seq_queue_info_get_queue',
'snd_seq_queue_info_get_name', 'snd_seq_queue_info_get_owner',
'snd_seq_queue_info_get_locked', 'snd_seq_queue_info_get_flags',
'snd_seq_queue_info_set_name', 'snd_seq_queue_info_set_owner',
'snd_seq_queue_info_set_locked', 'snd_seq_queue_info_set_flags',
'snd_seq_create_queue', 'snd_seq_alloc_named_queue', 'snd_seq_alloc_queue',
'snd_seq_free_queue', 'snd_seq_get_queue_info', 'snd_seq_set_queue_info',
'snd_seq_query_named_queue', 'snd_seq_get_queue_usage',
'snd_seq_set_queue_usage', 'snd_seq_queue_status_sizeof',
'snd_seq_queue_status_malloc', 'snd_seq_queue_status_free',
'snd_seq_queue_status_copy', 'snd_seq_queue_status_get_queue',
'snd_seq_queue_status_get_events', 'snd_seq_queue_status_get_tick_time',
'snd_seq_queue_status_get_real_time', 'snd_seq_queue_status_get_status',
'snd_seq_get_queue_status', 'snd_seq_queue_tempo_sizeof',
'snd_seq_queue_tempo_malloc', 'snd_seq_queue_tempo_free',
'snd_seq_queue_tempo_copy', 'snd_seq_queue_tempo_get_queue',
'snd_seq_queue_tempo_get_tempo', 'snd_seq_queue_tempo_get_ppq',
'snd_seq_queue_tempo_get_skew', 'snd_seq_queue_tempo_get_skew_base',
'snd_seq_queue_tempo_set_tempo', 'snd_seq_queue_tempo_set_ppq',
'snd_seq_queue_tempo_set_skew', 'snd_seq_queue_tempo_set_skew_base',
'snd_seq_get_queue_tempo', 'snd_seq_set_queue_tempo',
'snd_seq_queue_timer_type_t', 'SND_SEQ_TIMER_ALSA',
'SND_SEQ_TIMER_MIDI_CLOCK', 'SND_SEQ_TIMER_MIDI_TICK',
'snd_seq_queue_timer_sizeof', 'snd_seq_queue_timer_malloc',
'snd_seq_queue_timer_free', 'snd_seq_queue_timer_copy',
'snd_seq_queue_timer_get_queue', 'snd_seq_queue_timer_get_type',
'snd_seq_queue_timer_get_id', 'snd_seq_queue_timer_get_resolution',
'snd_seq_queue_timer_set_type', 'snd_seq_queue_timer_set_id',
'snd_seq_queue_timer_set_resolution', 'snd_seq_get_queue_timer',
'snd_seq_set_queue_timer', 'snd_seq_free_event', 'snd_seq_event_length',
'snd_seq_event_output', 'snd_seq_event_output_buffer',
'snd_seq_event_output_direct', 'snd_seq_event_input',
'snd_seq_event_input_pending', 'snd_seq_drain_output',
'snd_seq_event_output_pending', 'snd_seq_extract_output',
'snd_seq_drop_output', 'snd_seq_drop_output_buffer', 'snd_seq_drop_input',
'snd_seq_drop_input_buffer', 'snd_seq_remove_events_t',
'SND_SEQ_REMOVE_INPUT', 'SND_SEQ_REMOVE_OUTPUT', 'SND_SEQ_REMOVE_DEST',
'SND_SEQ_REMOVE_DEST_CHANNEL', 'SND_SEQ_REMOVE_TIME_BEFORE',
'SND_SEQ_REMOVE_TIME_AFTER', 'SND_SEQ_REMOVE_TIME_TICK',
'SND_SEQ_REMOVE_EVENT_TYPE', 'SND_SEQ_REMOVE_IGNORE_OFF',
'SND_SEQ_REMOVE_TAG_MATCH', 'snd_seq_remove_events_sizeof',
'snd_seq_remove_events_malloc', 'snd_seq_remove_events_free',
'snd_seq_remove_events_copy', 'snd_seq_remove_events_get_condition',
'snd_seq_remove_events_get_queue', 'snd_seq_remove_events_get_time',
'snd_seq_remove_events_get_dest', 'snd_seq_remove_events_get_channel',
'snd_seq_remove_events_get_event_type', 'snd_seq_remove_events_get_tag',
'snd_seq_remove_events_set_condition', 'snd_seq_remove_events_set_queue',
'snd_seq_remove_events_set_time', 'snd_seq_remove_events_set_dest',
'snd_seq_remove_events_set_channel', 'snd_seq_remove_events_set_event_type',
'snd_seq_remove_events_set_tag', 'snd_seq_remove_events', 'snd_seq_set_bit',
'snd_seq_change_bit', 'snd_seq_get_bit', 'snd_seq_control_queue',
'snd_seq_create_simple_port', 'snd_seq_delete_simple_port',
'snd_seq_connect_from', 'snd_seq_connect_to', 'snd_seq_disconnect_from',
'snd_seq_disconnect_to', 'snd_seq_set_client_name',
'snd_seq_set_client_event_filter', 'snd_seq_set_client_pool_output',
'snd_seq_set_client_pool_output_room', 'snd_seq_set_client_pool_input',
'snd_seq_sync_output_queue', 'snd_seq_parse_address',
'snd_seq_reset_pool_output', 'snd_seq_reset_pool_input', 'snd_midi_event_t',
'snd_midi_event_new', 'snd_midi_event_resize_buffer', 'snd_midi_event_free',
'snd_midi_event_init', 'snd_midi_event_reset_encode',
'snd_midi_event_reset_decode', 'snd_midi_event_no_status',
'snd_midi_event_encode', 'snd_midi_event_encode_byte',
'snd_midi_event_decode', 'snd_instr_header_t', 'snd_instr_header_sizeof',
'snd_instr_header_malloc', 'snd_instr_header_free', 'snd_instr_header_copy',
'snd_instr_header_get_id', 'snd_instr_header_get_cluster',
'snd_instr_header_get_cmd', 'snd_instr_header_get_len',
'snd_instr_header_get_name', 'snd_instr_header_get_type',
'snd_instr_header_get_format', 'snd_instr_header_get_alias',
'snd_instr_header_get_data', 'snd_instr_header_get_follow_alias',
'snd_instr_header_set_id', 'snd_instr_header_set_cluster',
'snd_instr_header_set_cmd', 'snd_instr_header_set_len',
'snd_instr_header_set_name', 'snd_instr_header_set_type',
'snd_instr_header_set_format', 'snd_instr_header_set_alias',
'snd_instr_header_set_follow_alias', 'SND_SEQ_INSTR_ATYPE_DATA',
'SND_SEQ_INSTR_ATYPE_ALIAS', 'SND_SEQ_INSTR_TYPE0_DLS1',
'SND_SEQ_INSTR_TYPE0_DLS2', 'SND_SEQ_INSTR_TYPE1_SIMPLE',
'SND_SEQ_INSTR_TYPE1_SOUNDFONT', 'SND_SEQ_INSTR_TYPE1_GUS_PATCH',
'SND_SEQ_INSTR_TYPE1_INTERWAVE', 'SND_SEQ_INSTR_TYPE2_OPL2_3',
'SND_SEQ_INSTR_TYPE2_OPL4', 'SND_SEQ_INSTR_PUT_CMD_CREATE',
'SND_SEQ_INSTR_PUT_CMD_REPLACE', 'SND_SEQ_INSTR_PUT_CMD_MODIFY',
'SND_SEQ_INSTR_PUT_CMD_ADD', 'SND_SEQ_INSTR_PUT_CMD_REMOVE',
'SND_SEQ_INSTR_GET_CMD_FULL', 'SND_SEQ_INSTR_GET_CMD_PARTIAL',
'SND_SEQ_INSTR_QUERY_FOLLOW_ALIAS', 'SND_SEQ_INSTR_FREE_CMD_ALL',
'SND_SEQ_INSTR_FREE_CMD_PRIVATE', 'SND_SEQ_INSTR_FREE_CMD_CLUSTER',
'SND_SEQ_INSTR_FREE_CMD_SINGLE', 'snd_instr_fm_t',
'snd_instr_fm_convert_to_stream', 'snd_instr_fm_convert_from_stream',
'snd_instr_fm_free', 'snd_instr_simple_t',
'snd_instr_simple_convert_to_stream', 'snd_instr_simple_convert_from_stream',
'snd_instr_simple_free', 'snd_instr_iwffff_t', 'snd_iwffff_handle_t',
'snd_instr_iwffff_open', 'snd_instr_iwffff_open_rom',
'snd_instr_iwffff_open_rom_file', 'snd_instr_iwffff_close',
'snd_instr_iwffff_load', 'snd_instr_iwffff_convert_to_stream',
'snd_instr_iwffff_convert_from_stream', 'snd_instr_iwffff_free']
