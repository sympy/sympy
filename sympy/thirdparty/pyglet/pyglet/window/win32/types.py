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
__version__ = '$Id: $'

from ctypes import *

BOOL = c_int
DWORD = c_uint32
BYTE = c_char
LONG = c_long
WORD = c_short

HANDLE = c_void_p
HWND = HANDLE
HMONITOR = HANDLE
HGLOBAL = HANDLE
HDC = HANDLE
HBITMAP = HANDLE
LPARAM = c_long

WNDPROC = WINFUNCTYPE(c_long, c_int, c_uint, c_int, c_int)
TIMERPROC = WINFUNCTYPE(None, HWND, c_uint, POINTER(c_uint), c_uint)

class RECT(Structure):
    _fields_ = [
        ('left', c_long),
        ('top', c_long),
        ('right', c_long),
        ('bottom', c_long)
    ]

class WNDCLASS(Structure):
    _fields_ = [
        ('style', c_uint),
        ('lpfnWndProc', WNDPROC),
        ('cbClsExtra', c_int),
        ('cbWndExtra', c_int),
        ('hInstance', c_int),
        ('hIcon', c_int),
        ('hCursor', c_int),
        ('hbrBackground', c_int),
        ('lpszMenuName', c_char_p),
        ('lpszClassName', c_wchar_p)
    ]

class POINT(Structure):
    _fields_ = [
        ('x', c_long),
        ('y', c_long)
    ]

class MSG(Structure):
    _fields_ = [
        ('hwnd', c_int),
        ('message', c_uint),
        ('wParam', c_int),
        ('lParam', c_int),
        ('time', c_int),
        ('pt', POINT)
    ]

class PIXELFORMATDESCRIPTOR(Structure):
    _fields_ = [
        ('nSize', c_ushort),
        ('nVersion', c_ushort),
        ('dwFlags', c_ulong),
        ('iPixelType', c_ubyte),
        ('cColorBits', c_ubyte),
        ('cRedBits', c_ubyte),
        ('cRedShift', c_ubyte),
        ('cGreenBits', c_ubyte),
        ('cGreenShift', c_ubyte),
        ('cBlueBits', c_ubyte),
        ('cBlueShift', c_ubyte),
        ('cAlphaBits', c_ubyte),
        ('cAlphaShift', c_ubyte),
        ('cAccumBits', c_ubyte),
        ('cAccumRedBits', c_ubyte),
        ('cAccumGreenBits', c_ubyte),
        ('cAccumBlueBits', c_ubyte),
        ('cAccumAlphaBits', c_ubyte),
        ('cDepthBits', c_ubyte),
        ('cStencilBits', c_ubyte),
        ('cAuxBuffers', c_ubyte),
        ('iLayerType', c_ubyte),
        ('bReserved', c_ubyte),
        ('dwLayerMask', c_ulong),
        ('dwVisibleMask', c_ulong),
        ('dwDamageMask', c_ulong)
    ]

class TRACKMOUSEEVENT(Structure):
    _fields_ = [
        ('cbSize', DWORD),
        ('dwFlags', DWORD),
        ('hwndTrack', HWND),
        ('dwHoverTime', DWORD)
    ]

class MINMAXINFO(Structure):
    _fields_ = [
        ('ptReserved', POINT),
        ('ptMaxSize', POINT),
        ('ptMaxPosition', POINT),
        ('ptMinTrackSize', POINT),
        ('ptMaxTrackSize', POINT)
    ]

class RGBQUAD(Structure):
    _fields_ = [
        ('rgbBlue', BYTE),
        ('rgbGreen', BYTE),
        ('rgbRed', BYTE),
        ('rgbReserved', BYTE),
    ]

class BITMAPINFOHEADER(Structure):
    _fields_ = [
        ('biSize', DWORD),
        ('biWidth', LONG),
        ('biHeight', LONG),
        ('biPlanes', WORD),
        ('biBitCount', WORD),
        ('biCompression', DWORD),
        ('biSizeImage', DWORD),
        ('biXPelsPerMeter', LONG),
        ('biYPelsPerMeter', LONG),
        ('biClrUsed', DWORD),
        ('biClrImportant', DWORD),
    ]

class BITMAPINFO(Structure):
    _fields_ = [
        ('bmiHeader', BITMAPINFOHEADER),
        ('bmiColors', RGBQUAD * 1)
    ]

class CIEXYZ(Structure):
    _fields_ = [
        ('ciexyzX', DWORD),
        ('ciexyzY', DWORD),
        ('ciexyzZ', DWORD),
    ]

class CIEXYZTRIPLE(Structure):
    _fields_ = [
        ('ciexyzRed', CIEXYZ),
        ('ciexyzBlue', CIEXYZ),
        ('ciexyzGreen', CIEXYZ),
    ]

class BITMAPV5HEADER(Structure):
    _fields_ = [
        ('bV5Size', DWORD),
        ('bV5Width', LONG),
        ('bV5Height', LONG),
        ('bV5Planes', WORD),
        ('bV5BitCount', WORD),
        ('bV5Compression', DWORD),
        ('bV5SizeImage', DWORD),
        ('bV5XPelsPerMeter', LONG),
        ('bV5YPelsPerMeter', LONG),
        ('bV5ClrUsed', DWORD),
        ('bV5ClrImportant', DWORD),
        ('bV5RedMask', DWORD),
        ('bV5GreenMask', DWORD),
        ('bV5BlueMask', DWORD),
        ('bV5AlphaMask', DWORD),
        ('bV5CSType', DWORD),
        ('bV5Endpoints', CIEXYZTRIPLE),
        ('bV5GammaRed', DWORD),
        ('bV5GammaGreen', DWORD),
        ('bV5GammaBlue', DWORD),
        ('bV5Intent', DWORD),
        ('bV5ProfileData', DWORD),
        ('bV5ProfileSize', DWORD),
        ('bV5Reserved', DWORD),
    ]

class ICONINFO(Structure):
    _fields_ = [
        ('fIcon', BOOL),
        ('xHotspot', DWORD),
        ('yHotspot', DWORD),
        ('hbmMask', HBITMAP),
        ('hbmColor', HBITMAP)
    ]
