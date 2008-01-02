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

'''Additional hand-coded GL extensions.

These are hand-wrapped extension tokens and functions that are in
the OpenGL Extension Registry but have not yet been added to either 
the registry's glext.h or nVidia's glext.h.  Remove wraps from here
when the headers are updated (and glext_arb.py or glext_nv.py are
regenerated).

When adding an extension here, include the name and URL, and any tokens and
functions appearing under "New Tokens" and "New Procedures" headings.  Don't
forget to add the GL_/gl prefix.

Unnumbered extensions in the registry are not included.
'''

__docformat__ = 'restructuredtext'
__version__ = '$Id: glext_missing.py 1322 2007-10-23 12:58:03Z Alex.Holkner $'

from ctypes import *
from pyglet.gl.lib import link_GL as _link_function
from pyglet.gl.lib import c_ptrdiff_t

# At time of writing, ABI glext.h was last updated 2005/06/20, so numbered
# non-ARB extensions from 312 on must be included here.

# GL_EXT_packed_depth_stencil
# http://oss.sgi.com/projects/ogl-sample/registry/EXT/packed_depth_stencil.txt

GL_DEPTH_STENCIL_EXT                 = 0x84F9
GL_UNSIGNED_INT_24_8_EXT             = 0x84FA
GL_DEPTH24_STENCIL8_EXT              = 0x88F0
GL_TEXTURE_STENCIL_SIZE_EXT          = 0x88F1

# GL_EXT_texture_sRGB
# http://oss.sgi.com/projects/ogl-sample/registry/EXT/texture_sRGB.txt

GL_SRGB_EXT                                 = 0x8C40
GL_SRGB8_EXT                                = 0x8C41
GL_SRGB_ALPHA_EXT                           = 0x8C42
GL_SRGB8_ALPHA8_EXT                         = 0x8C43
GL_SLUMINANCE_ALPHA_EXT                     = 0x8C44
GL_SLUMINANCE8_ALPHA8_EXT                   = 0x8C45
GL_SLUMINANCE_EXT                           = 0x8C46
GL_SLUMINANCE8_EXT                          = 0x8C47
GL_COMPRESSED_SRGB_EXT                      = 0x8C48
GL_COMPRESSED_SRGB_ALPHA_EXT                = 0x8C49
GL_COMPRESSED_SLUMINANCE_EXT                = 0x8C4A
GL_COMPRESSED_SLUMINANCE_ALPHA_EXT          = 0x8C4B
GL_COMPRESSED_SRGB_S3TC_DXT1_EXT            = 0x8C4C
GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT1_EXT      = 0x8C4D
GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT3_EXT      = 0x8C4E
GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT5_EXT      = 0x8C4F

# GL_EXT_stencil_clear_tag
# http://oss.sgi.com/projects/ogl-sample/registry/EXT/stencil_clear_tag.txt

GLuint = c_uint 	# /usr/include/GL/gl.h:62
GLsizei = c_int 	# /usr/include/GL/gl.h:59
glStencilClearTagEXT = _link_function(
    'glStencilClearTagEXT', None, [GLsizei, GLuint])

GL_STENCIL_TAG_BITS_EXT              = 0x88F2
GL_STENCIL_CLEAR_TAG_VALUE_EXT       = 0x88F3

# GL_EXT_framebuffer_blit
# http://oss.sgi.com/projects/ogl-sample/registry/EXT/framebuffer_blit.txt

GLenum = c_uint 	# /usr/include/GL/gl.h:53
GLint = c_int 	# /usr/include/GL/gl.h:58
glBlitFramebufferEXT = _link_function(
    'glBlitFramebufferEXT', None, [GLint, GLint, GLint, GLint,
                                   GLint, GLint, GLint, GLint,
                                   GLuint, GLenum])

GL_READ_FRAMEBUFFER_EXT              = 0x8CA8
GL_DRAW_FRAMEBUFFER_EXT              = 0x8CA9
GL_DRAW_FRAMEBUFFER_BINDING_EXT      = 0x8CA6
GL_READ_FRAMEBUFFER_BINDING_EXT      = 0x8CAA

# GL_EXT_framebuffer_multisample
# http://oss.sgi.com/projects/ogl-sample/registry/EXT/framebuffer_multisample.txt

GL_RENDERBUFFER_SAMPLES_EXT          = 0x8CAB

# GL_MESAX_texture_stack
# http://oss.sgi.com/projects/ogl-sample/registry/MESAX/texture_stack.txt

GL_TEXTURE_1D_STACK_MESAX            = 0x8759
GL_TEXTURE_2D_STACK_MESAX            = 0x875A
GL_PROXY_TEXTURE_1D_STACK_MESAX      = 0x875B
GL_PROXY_TEXTURE_2D_STACK_MESAX      = 0x875C
GL_TEXTURE_1D_STACK_BINDING_MESAX    = 0x875D
GL_TEXTURE_2D_STACK_BINDING_MESAX    = 0x875E
