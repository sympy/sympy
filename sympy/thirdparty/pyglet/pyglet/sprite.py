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

'''Display positioned, scaled and rotated images.

A sprite is an instance of an image displayed on-screen.  Multiple sprites can
display the same image at different positions on the screen.  Sprites can also
be scaled larger or smaller, rotated at any angle and drawn at a fractional
opacity.

The following complete example loads a ``"ball.png"`` image and creates a
sprite for that image.  The sprite is then drawn in the window's
draw event handler::

    import pyglet

    ball_image = pyglet.image.load('ball.png')
    ball = pyglet.sprite.Sprite(ball_image, x=50, y=50)

    window = pyglet.window.Window()

    @window.event
    def on_draw():
        ball.draw()

    pyglet.app.run()

The sprite can be moved by modifying the ``x`` and ``y`` properties.  Other
properties determine the sprite's rotation, scale and opacity.

The sprite's positioning, rotation and scaling all honor the original
image's anchor (anchor_x, anchor_y).


Drawing multiple sprites
========================

Sprites can be "batched" together and drawn at once more quickly than if each
of their ``draw`` methods were called individually.  The following example
creates one hundred ball sprites and adds each of them to a `Batch`.  The
entire batch of sprites is then drawn in one call::

    batch = pyglet.graphics.Batch()

    ball_sprites = []
    for i in range(100):
        x, y = i * 10, 50
        ball_sprites.append(pyglet.sprite.Sprite(ball_image, x, y, batch=batch)

    @window.event
    def on_draw():
        batch.draw()

Sprites can be freely modified in any way even after being added to a batch,
however a sprite can belong to at most one batch.  See the documentation for
`pyglet.graphics` for more details on batched rendering, and grouping of
sprites within batches.

:since: pyglet 1.1
'''

__docformat__ = 'restructuredtext'
__version__ = '$Id: sprite.py 2120 2008-06-17 12:08:05Z Alex.Holkner $'

import math
import sys

from pyglet.gl import *
from pyglet import clock
from pyglet import event
from pyglet import graphics
from pyglet import image

_is_epydoc = hasattr(sys, 'is_epydoc') and sys.is_epydoc

class SpriteGroup(graphics.Group):
    '''Shared sprite rendering group.

    The group is automatically coallesced with other sprite groups sharing the
    same parent group, texture and blend parameters.
    '''
    def __init__(self, texture, blend_src, blend_dest, parent=None):
        '''Create a sprite group.

        The group is created internally within `Sprite`; applications usually
        do not need to explicitly create it.

        :Parameters:
            `texture` : `Texture`
                The (top-level) texture containing the sprite image.
            `blend_src` : int
                OpenGL blend source mode; for example,
                ``GL_SRC_ALPHA``.
            `blend_dest` : int
                OpenGL blend destination mode; for example,
                ``GL_ONE_MINUS_SRC_ALPHA``.
            `parent` : `Group`
                Optional parent group.

        '''
        super(SpriteGroup, self).__init__(parent)
        self.texture = texture
        self.blend_src = blend_src
        self.blend_dest = blend_dest

    def set_state(self):
        glEnable(self.texture.target)
        glBindTexture(self.texture.target, self.texture.id)

        glPushAttrib(GL_COLOR_BUFFER_BIT)
        glEnable(GL_BLEND)
        glBlendFunc(self.blend_src, self.blend_dest)

    def unset_state(self):
        glPopAttrib()
        glDisable(self.texture.target)

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, self.texture)

    def __eq__(self, other):
        return (other.__class__ is self.__class__ and
                self.parent is other.parent and
                self.texture.target == other.texture.target and
                self.texture.id == other.texture.id and
                self.blend_src == other.blend_src and
                self.blend_dest == other.blend_dest)

    def __hash__(self):
        return hash((id(self.parent),
                     self.texture.id, self.texture.target,
                     self.blend_src, self.blend_dest))

class Sprite(event.EventDispatcher):
    '''Instance of an on-screen image.

    See the module documentation for usage.
    '''
    _batch = None
    _animation = None
    _rotation = 0
    _opacity = 255
    _rgb = (255, 255, 255)
    _scale = 1.0
    _visible = True
    _vertex_list = None

    def __init__(self,
                 img, x=0, y=0,
                 blend_src=GL_SRC_ALPHA,
                 blend_dest=GL_ONE_MINUS_SRC_ALPHA,
                 batch=None,
                 group=None,
                 usage='dynamic'):
        '''Create a sprite.

        :Parameters:
            `img` : `AbstractImage` or `Animation`
                Image or animation to display.
            `x` : int
                X coordinate of the sprite.
            `y` : int
                Y coordinate of the sprite.
            `blend_src` : int
                OpenGL blend source mode.  The default is suitable for
                compositing sprites drawn from back-to-front.
            `blend_dest` : int
                OpenGL blend destination mode.  The default is suitable for
                compositing sprites drawn from back-to-front.
            `batch` : `Batch`
                Optional batch to add the sprite to.
            `group` : `Group`
                Optional parent group of the sprite.
            `usage` : str
                Vertex buffer object usage hint, one of ``"none"`` (default),
                ``"stream"``, ``"dynamic"`` or ``"static"``.  Applies
                only to vertex data.

        '''
        if batch is not None:
            self._batch = batch

        self._x = x
        self._y = y

        if isinstance(img, image.Animation):
            self._animation = img
            self._frame_index = 0
            self._texture = img.frames[0].image.get_texture()
            self._next_dt = img.frames[0].duration
            if self._next_dt:
                clock.schedule_once(self._animate, self._next_dt)
        else:
            self._texture = img.get_texture()

        self._group = SpriteGroup(self._texture, blend_src, blend_dest, group)
        self._usage = usage
        self._create_vertex_list()

    def __del__(self):
        try:
            if self._vertex_list is not None:
                self._vertex_list.delete()
        except:
            pass

    def delete(self):
        '''Force immediate removal of the sprite from video memory.

        This is often necessary when using batches, as the Python garbage
        collector will not necessarily call the finalizer as soon as the
        sprite is garbage.
        '''
        if self._animation:
            clock.unschedule(self._animate)
        self._vertex_list.delete()
        self._vertex_list = None
        self._texture = None

        # Easy way to break circular reference, speeds up GC
        self._group = None

    def _animate(self, dt):
        self._frame_index += 1
        if self._frame_index >= len(self._animation.frames):
            self._frame_index = 0
            self.dispatch_event('on_animation_end')
            if self._vertex_list is None:
                return # Deleted in event handler.

        frame = self._animation.frames[self._frame_index]
        self._set_texture(frame.image.get_texture())

        if frame.duration is not None:
            duration = frame.duration - (self._next_dt - dt)
            duration = min(max(0, duration), frame.duration)
            clock.schedule_once(self._animate, duration)
            self._next_dt = duration
        else:
            self.dispatch_event('on_animation_end')

    def _set_batch(self, batch):
        if self._batch == batch:
            return

        if batch is not None and self._batch is not None:
            self._batch.migrate(self._vertex_list, GL_QUADS, self._group, batch)
            self._batch = batch
        else:
            self._vertex_list.delete()
            self._batch = batch
            self._create_vertex_list()

    def _get_batch(self):
        return self._batch

    batch = property(_get_batch, _set_batch,
                     doc='''Graphics batch.

    The sprite can be migrated from one batch to another, or removed from its
    batch (for individual drawing).  Note that this can be an expensive
    operation.

    :type: `Batch`
    ''')

    def _set_group(self, group):
        if self._group.parent == group:
            return

        self._group = SpriteGroup(self._texture,
                                  self._group.blend_src,
                                  self._group.blend_dest,
                                  group)

        if self._batch is not None:
           self._batch.migrate(self._vertex_list, GL_QUADS, self._group,
                               self._batch)

    def _get_group(self):
        return self._group.parent

    group = property(_get_group, _set_group,
                     doc='''Parent graphics group.

    The sprite can change its rendering group, however this can be an
    expensive operation.

    :type: `Group`
    ''')

    def _get_image(self):
        if self._animation:
            return self._animation
        return self._texture

    def _set_image(self, img):
        if self._animation is not None:
            clock.unschedule(self._animate)
            self._animation = None

        if isinstance(img, image.Animation):
            self._animation = img
            self._frame_index = 0
            self._set_texture(img.frames[0].image.get_texture())
            self._next_dt = img.frames[0].duration
            clock.schedule_once(self._animate, self._next_dt)
        else:
            self._set_texture(img.get_texture())
        self._update_position()

    image = property(_get_image, _set_image,
                     doc='''Image or animation to display.

    :type: `AbstractImage` or `Animation`
    ''')

    def _set_texture(self, texture):
        if texture.id is not self._texture.id:
            self._group = SpriteGroup(texture,
                                      self._group.blend_src,
                                      self._group.blend_dest,
                                      self._group.parent)
            if self._batch is None:
                self._vertex_list.tex_coords[:] = texture.tex_coords
            else:
                self._vertex_list.delete()
                self._texture = texture
                self._create_vertex_list()
        else:
            self._vertex_list.tex_coords[:] = texture.tex_coords
        self._texture = texture

    def _create_vertex_list(self):
        if self._batch is None:
            self._vertex_list = graphics.vertex_list(4,
                'v2i/%s' % self._usage, 
                'c4B', ('t3f', self._texture.tex_coords))
        else:
            self._vertex_list = self._batch.add(4, GL_QUADS, self._group,
                'v2i/%s' % self._usage, 
                'c4B', ('t3f', self._texture.tex_coords))
        self._update_position()
        self._update_color()

    def _update_position(self):
        img = self._texture
        if not self._visible:
            self._vertex_list.vertices[:] = [0, 0, 0, 0, 0, 0, 0, 0]
        elif self._rotation:
            x1 = -img.anchor_x * self._scale
            y1 = -img.anchor_y * self._scale
            x2 = x1 + img.width * self._scale
            y2 = y1 + img.height * self._scale
            x = self._x
            y = self._y

            r = -math.radians(self._rotation)
            cr = math.cos(r)
            sr = math.sin(r)
            ax = int(x1 * cr - y1 * sr + x)
            ay = int(x1 * sr + y1 * cr + y)
            bx = int(x2 * cr - y1 * sr + x)
            by = int(x2 * sr + y1 * cr + y)
            cx = int(x2 * cr - y2 * sr + x)
            cy = int(x2 * sr + y2 * cr + y)
            dx = int(x1 * cr - y2 * sr + x)
            dy = int(x1 * sr + y2 * cr + y)

            self._vertex_list.vertices[:] = [ax, ay, bx, by, cx, cy, dx, dy]
        elif self._scale != 1.0:
            x1 = int(self._x - img.anchor_x * self._scale)
            y1 = int(self._y - img.anchor_y * self._scale)
            x2 = int(x1 + img.width * self._scale)
            y2 = int(y1 + img.height * self._scale)
            self._vertex_list.vertices[:] = [x1, y1, x2, y1, x2, y2, x1, y2]
        else:
            x1 = int(self._x - img.anchor_x)
            y1 = int(self._y - img.anchor_y)
            x2 = x1 + img.width
            y2 = y1 + img.height
            self._vertex_list.vertices[:] = [x1, y1, x2, y1, x2, y2, x1, y2]

    def _update_color(self):
        r, g, b = self._rgb
        self._vertex_list.colors[:] = [r, g, b, int(self._opacity)] * 4

    def set_position(self, x, y):
        '''Set the X and Y coordinates of the sprite simultaneously.

        :Parameters:
            `x` : int
                X coordinate of the sprite.
            `y` : int
                Y coordinate of the sprite.

        '''
        self._x = x
        self._y = y
        self._update_position()

    position = property(lambda self: (self._x, self._y),
                        lambda self, t: self.set_position(*t),
                        doc='''The (x, y) coordinates of the sprite.

    :type: (int, int)
    ''')

    def _set_x(self, x):
        self._x = x
        self._update_position()

    x = property(lambda self: self._x, _set_x,
                 doc='''X coordinate of the sprite.

    :type: int
    ''')

    def _set_y(self, y):
        self._y = y
        self._update_position()

    y = property(lambda self: self._y, _set_y,
                 doc='''Y coordinate of the sprite.

    :type: int
    ''')

    def _set_rotation(self, rotation):
        self._rotation = rotation
        self._update_position()

    rotation = property(lambda self: self._rotation, _set_rotation,
                        doc='''Clockwise rotation of the sprite, in degrees.

    The sprite image will be rotated about its image's (anchor_x, anchor_y)
    position.

    :type: float
    ''')

    def _set_scale(self, scale):
        self._scale = scale
        self._update_position()

    scale = property(lambda self: self._scale, _set_scale,
                     doc='''Scaling factor.

    A scaling factor of 1 (the default) has no effect.  A scale of 2 will draw
    the sprite at twice the native size of its image.

    :type: float
    ''')

    width = property(lambda self: int(self._texture.width * self._scale),
                     doc='''Scaled width of the sprite.

    Read-only.  Invariant under rotation.

    :type: int
    ''')

    height = property(lambda self: int(self._texture.height * self._scale),
                      doc='''Scaled height of the sprite.

    Read-only.  Invariant under rotation.

    :type: int
    ''')

    def _set_opacity(self, opacity):
        self._opacity = opacity
        self._update_color()

    opacity = property(lambda self: self._opacity, _set_opacity,
                       doc='''Blend opacity.

    This property sets the alpha component of the colour of the sprite's
    vertices.  With the default blend mode (see the constructor), this
    allows the sprite to be drawn with fractional opacity, blending with the
    background.

    An opacity of 255 (the default) has no effect.  An opacity of 128 will
    make the sprite appear translucent.

    :type: int
    ''')

    def _set_color(self, rgb):
        self._rgb = map(int, rgb)
        self._update_color()

    color = property(lambda self: self._rgb, _set_color,
                       doc='''Blend color.

    This property sets the color of the sprite's vertices. This allows the
    sprite to be drawn with a color tint.
    
    The color is specified as an RGB tuple of integers ``(red, green, blue)``.
    Each color component must be in the range 0 (dark) to 255 (saturated).
    
    :type: (int, int, int)
    ''')

    def _set_visible(self, visible):
        self._visible = visible
        self._update_position()

    visible = property(lambda self: self._visible, _set_visible,
                       '''True if the sprite will be drawn.

    :type: bool
    ''')


    def draw(self):
        '''Draw the sprite at its current position.

        See the module documentation for hints on drawing multiple sprites
        efficiently.
        '''
        self._group.set_state_recursive()
        self._vertex_list.draw(GL_QUADS)
        self._group.unset_state_recursive()

    if _is_epydoc:
        def on_animation_end():
            '''The sprite animation reached the final frame.

            The event is triggered only if the sprite has an animation, not an
            image.  For looping animations, the event is triggered each time
            the animation loops.

            :event:
            '''

Sprite.register_event_type('on_animation_end')
