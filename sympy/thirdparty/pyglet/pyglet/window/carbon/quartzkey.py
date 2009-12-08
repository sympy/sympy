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

'''
'''

__docformat__ = 'restructuredtext'
__version__ = ''

from pyglet.window import key

# From SDL: src/video/quartz/SDL_QuartzKeys.h
# These are the Macintosh key scancode constants -- from Inside Macintosh

QZ_ESCAPE = 0x35
QZ_F1 = 0x7A
QZ_F2 = 0x78
QZ_F3 = 0x63
QZ_F4 = 0x76
QZ_F5 = 0x60
QZ_F6 = 0x61
QZ_F7 = 0x62
QZ_F8 = 0x64
QZ_F9 = 0x65
QZ_F10 = 0x6D
QZ_F11 = 0x67
QZ_F12 = 0x6F
QZ_PRINT = 0x69
QZ_SCROLLOCK = 0x6B
QZ_PAUSE = 0x71
QZ_POWER = 0x7F
QZ_BACKQUOTE = 0x32
QZ_1 = 0x12
QZ_2 = 0x13
QZ_3 = 0x14
QZ_4 = 0x15
QZ_5 = 0x17
QZ_6 = 0x16
QZ_7 = 0x1A
QZ_8 = 0x1C
QZ_9 = 0x19
QZ_0 = 0x1D
QZ_MINUS = 0x1B
QZ_EQUALS = 0x18
QZ_BACKSPACE = 0x33
QZ_INSERT = 0x72
QZ_HOME = 0x73
QZ_PAGEUP = 0x74
QZ_NUMLOCK = 0x47
QZ_KP_EQUALS = 0x51
QZ_KP_DIVIDE = 0x4B
QZ_KP_MULTIPLY = 0x43
QZ_TAB = 0x30
QZ_q = 0x0C
QZ_w = 0x0D
QZ_e = 0x0E
QZ_r = 0x0F
QZ_t = 0x11
QZ_y = 0x10
QZ_u = 0x20
QZ_i = 0x22
QZ_o = 0x1F
QZ_p = 0x23
QZ_LEFTBRACKET = 0x21
QZ_RIGHTBRACKET = 0x1E
QZ_BACKSLASH = 0x2A
QZ_DELETE = 0x75
QZ_END = 0x77
QZ_PAGEDOWN = 0x79
QZ_KP7 = 0x59
QZ_KP8 = 0x5B
QZ_KP9 = 0x5C
QZ_KP_MINUS = 0x4E
QZ_CAPSLOCK = 0x39
QZ_a = 0x00
QZ_s = 0x01
QZ_d = 0x02
QZ_f = 0x03
QZ_g = 0x05
QZ_h = 0x04
QZ_j = 0x26
QZ_k = 0x28
QZ_l = 0x25
QZ_SEMICOLON = 0x29
QZ_QUOTE = 0x27
QZ_RETURN = 0x24
QZ_KP4 = 0x56
QZ_KP5 = 0x57
QZ_KP6 = 0x58
QZ_KP_PLUS = 0x45
QZ_LSHIFT = 0x38
QZ_z = 0x06
QZ_x = 0x07
QZ_c = 0x08
QZ_v = 0x09
QZ_b = 0x0B
QZ_n = 0x2D
QZ_m = 0x2E
QZ_COMMA = 0x2B
QZ_PERIOD = 0x2F
QZ_SLASH = 0x2C
QZ_RSHIFT = 0x3C
QZ_UP = 0x7E
QZ_KP1 = 0x53
QZ_KP2 = 0x54
QZ_KP3 = 0x55
QZ_KP_ENTER = 0x4C
QZ_LCTRL = 0x3B
QZ_LALT = 0x3A
QZ_LMETA = 0x37
QZ_SPACE = 0x31
QZ_RMETA = 0x36
QZ_RALT = 0x3D
QZ_RCTRL = 0x3E
QZ_LEFT = 0x7B
QZ_DOWN = 0x7D
QZ_RIGHT = 0x7C
QZ_KP0 = 0x52
QZ_KP_PERIOD = 0x41
QZ_IBOOK_ENTER = 0x34
QZ_IBOOK_LEFT = 0x3B
QZ_IBOOK_RIGHT = 0x3C
QZ_IBOOK_DOWN = 0x3D
QZ_IBOOK_UP = 0x3E

keymap = {
    QZ_ESCAPE: key.ESCAPE,
    QZ_F1: key.F1,
    QZ_F2: key.F2,
    QZ_F3: key.F3,
    QZ_F4: key.F4,
    QZ_F5: key.F5,
    QZ_F6: key.F6,
    QZ_F7: key.F7,
    QZ_F8: key.F8,
    QZ_F9: key.F9,
    QZ_F10: key.F10,
    QZ_F11: key.F11,
    QZ_F12: key.F12,
    QZ_PRINT: key.PRINT,
    QZ_SCROLLOCK: key.SCROLLLOCK,
    QZ_PAUSE: key.PAUSE,
    #QZ_POWER: key.POWER,
    QZ_BACKQUOTE: key.QUOTELEFT,
    QZ_1: key._1,
    QZ_2: key._2,
    QZ_3: key._3,
    QZ_4: key._4,
    QZ_5: key._5,
    QZ_6: key._6,
    QZ_7: key._7,
    QZ_8: key._8,
    QZ_9: key._9,
    QZ_0: key._0,
    QZ_MINUS: key.MINUS,
    QZ_EQUALS: key.EQUAL,
    QZ_BACKSPACE: key.BACKSPACE,
    QZ_INSERT: key.INSERT,
    QZ_HOME: key.HOME,
    QZ_PAGEUP: key.PAGEUP,
    QZ_NUMLOCK: key.NUMLOCK,
    QZ_KP_EQUALS: key.NUM_EQUAL,
    QZ_KP_DIVIDE: key.NUM_DIVIDE,
    QZ_KP_MULTIPLY: key.NUM_MULTIPLY,
    QZ_TAB: key.TAB,
    QZ_q: key.Q,
    QZ_w: key.W,
    QZ_e: key.E,
    QZ_r: key.R,
    QZ_t: key.T,
    QZ_y: key.Y,
    QZ_u: key.U,
    QZ_i: key.I,
    QZ_o: key.O,
    QZ_p: key.P,
    QZ_LEFTBRACKET: key.BRACKETLEFT,
    QZ_RIGHTBRACKET: key.BRACKETRIGHT,
    QZ_BACKSLASH: key.BACKSLASH,
    QZ_DELETE: key.DELETE,
    QZ_END: key.END,
    QZ_PAGEDOWN: key.PAGEDOWN,
    QZ_KP7: key.NUM_7,
    QZ_KP8: key.NUM_8,
    QZ_KP9: key.NUM_9,
    QZ_KP_MINUS: key.NUM_SUBTRACT,
    QZ_CAPSLOCK: key.CAPSLOCK,
    QZ_a: key.A,
    QZ_s: key.S,
    QZ_d: key.D,
    QZ_f: key.F,
    QZ_g: key.G,
    QZ_h: key.H,
    QZ_j: key.J,
    QZ_k: key.K,
    QZ_l: key.L,
    QZ_SEMICOLON: key.SEMICOLON,
    QZ_QUOTE: key.APOSTROPHE,
    QZ_RETURN: key.RETURN,
    QZ_KP4: key.NUM_4,
    QZ_KP5: key.NUM_5,
    QZ_KP6: key.NUM_6,
    QZ_KP_PLUS: key.NUM_ADD,
    QZ_LSHIFT: key.LSHIFT,
    QZ_z: key.Z,
    QZ_x: key.X,
    QZ_c: key.C,
    QZ_v: key.V,
    QZ_b: key.B,
    QZ_n: key.N,
    QZ_m: key.M,
    QZ_COMMA: key.COMMA,
    QZ_PERIOD: key.PERIOD,
    QZ_SLASH: key.SLASH,
    QZ_RSHIFT: key.RSHIFT,
    QZ_UP: key.UP,
    QZ_KP1: key.NUM_1,
    QZ_KP2: key.NUM_2,
    QZ_KP3: key.NUM_3,
    QZ_KP_ENTER: key.NUM_ENTER,
    QZ_LCTRL: key.LCTRL,
    QZ_LALT: key.LALT,
    QZ_LMETA: key.LMETA,
    QZ_SPACE: key.SPACE,
    QZ_RMETA: key.RMETA,
    QZ_RALT: key.RALT,
    QZ_RCTRL: key.RCTRL,
    QZ_LEFT: key.LEFT,
    QZ_DOWN: key.DOWN,
    QZ_RIGHT: key.RIGHT,
    QZ_KP0: key.NUM_0,
    QZ_KP_PERIOD: key.NUM_DECIMAL,
    QZ_IBOOK_ENTER: key.ENTER,
    QZ_IBOOK_LEFT: key.LEFT,
    QZ_IBOOK_RIGHT: key.RIGHT,
    QZ_IBOOK_DOWN: key.DOWN,
    QZ_IBOOK_UP: key.UP,
}


