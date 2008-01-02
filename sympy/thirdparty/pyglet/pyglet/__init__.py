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

'''pyglet is a cross-platform games and multimedia package.

Detailed documentation is available at http://www.pyglet.org
'''

__docformat__ = 'restructuredtext'
__version__ = '$Id: __init__.py 1558 2007-12-29 00:54:05Z Alex.Holkner $'

import os
import sys

#: The release version of this pyglet installation.  
#:
#: Valid only if pyglet was installed from a source or binary distribution
#: (i.e. not in a checked-out copy from SVN).
#: 
#: Use setuptools if you need to check for a specific release version, e.g.::
#:
#:    >>> import pyglet
#:    >>> from pkg_resources import parse_version
#:    >>> parse_version(pyglet.version) >= parse_version('1.0')
#:    False
#:
version = '1.0beta3'

def _require_ctypes_version(version):
    # Check ctypes version
    import ctypes
    req = [int(i) for i in version.split('.')]
    have = [int(i) for i in ctypes.__version__.split('.')]
    if not tuple(have) >= tuple(req):
        raise ImportError('pyglet requires ctypes %s or later.' % version)
_require_ctypes_version('1.0.0')

_enable_optimisations = not __debug__
if getattr(sys, 'frozen', None):
    _enable_optimisations = True

#: Global dict of pyglet options.  To change an option from its default, you
#: must import ``pyglet`` before any sub-packages.  For example::
#:
#:      import pyglet
#:      pyglet.options['debug_gl'] = False
#:
#: The default options can be overridden from the OS environment.  The
#: corresponding environment variable for each option key is prefaced by
#: ``PYGLET_``.  For example, in Bash you can set the ``debug_gl`` option with::
#:
#:      PYGLET_DEBUG_GL=True; export PYGLET_DEBUG_GL
#: 
#: For options requiring a tuple of values, separate each value with a comma.
#:
#: The non-development options are:
#:
#: audio
#:     A sequence of the names of audio modules to attempt to load, in
#:     order of preference.  Valid driver names are:
#:
#:     * directsound, the Windows DirectSound audio module (Windows only)
#:     * alsa, the ALSA audio module (Linux only) 
#:     * openal, the OpenAL audio module
#:     * silent, no audio
#: debug_gl
#:     If True, all calls to OpenGL functions are checked afterwards for
#:     errors using ``glGetError``.  This will severely impact performance,
#:     but provides useful exceptions at the point of failure.  By default,
#:     this option is enabled if ``__debug__`` is (i.e., if Python was not run
#:     with the -O option).  It is disabled by default when pyglet is "frozen"
#:     within a py2exe or py2app library archive.
#: vsync
#:     If set, the `pyglet.window.Window.vsync` property is ignored, and
#:     this option overrides it (to either force vsync on or off).  If unset,
#:     or set to None, the `pyglet.window.Window.vsync` property behaves
#:     as documented.
#:
options = {
    'audio': ('directsound', 'openal', 'alsa', 'silent'),
    'font': ('gdiplus', 'win32'), # ignored outside win32; win32 is deprecated
    'debug_font': False,
    'debug_gl': not _enable_optimisations,
    'debug_media': False,
    'debug_win32': False,
    'vsync': None,
}

_option_types = {
    'audio': tuple,
    'font': tuple,
    'debug_font': bool,
    'debug_gl': bool,
    'debug_media': bool,
    'debug_win32': bool,
    'vsync': bool,
}

def _read_environment():
    '''Read defaults for options from environment'''
    for key in options:
        env = 'PYGLET_%s' % key.upper()
        try:
            value = os.environ['PYGLET_%s' % key.upper()]
            if _option_types[key] is tuple:
                options[key] = value.split(',')
            elif _option_types[key] is bool:
                options[key] = value in ('true', 'TRUE', 'True', '1')
        except KeyError:
            pass
_read_environment()

if sys.platform == 'cygwin':
    # This hack pretends that the posix-like ctypes provides windows
    # functionality.  COM does not work with this hack, so there is no
    # DirectSound support.
    import ctypes
    ctypes.windll = ctypes.cdll
    ctypes.oledll = ctypes.cdll
    ctypes.WINFUNCTYPE = ctypes.CFUNCTYPE
    ctypes.HRESULT = ctypes.c_long
