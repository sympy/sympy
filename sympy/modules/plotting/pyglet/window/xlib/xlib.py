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
'''Wrapper for X11

Generated with:
tools/genwrappers.py

Do not modify this file.
'''

__docformat__ =  'restructuredtext'
__version__ = '$Id: xlib.py 1033 2007-07-13 03:38:16Z Alex.Holkner $'

import ctypes
from ctypes import *

import pyglet.lib

_lib = pyglet.lib.load_library('X11')

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


import pyglet.gl.glx

XlibSpecificationRelease = 6 	# /usr/include/X11/Xlib.h:39
X_PROTOCOL = 11 	# /usr/include/X11/X.h:58
X_PROTOCOL_REVISION = 0 	# /usr/include/X11/X.h:59
XID = pyglet.gl.glx.XID
Mask = c_ulong 	# /usr/include/X11/X.h:75
Atom = c_ulong 	# /usr/include/X11/X.h:79
VisualID = pyglet.gl.glx.VisualID
Time = c_ulong 	# /usr/include/X11/X.h:82
Window = pyglet.gl.glx.Window
Drawable = XID 	# /usr/include/X11/X.h:102
Font = pyglet.gl.glx.Font
Pixmap = pyglet.gl.glx.Pixmap
Cursor = XID 	# /usr/include/X11/X.h:108
Colormap = XID 	# /usr/include/X11/X.h:109
GContext = XID 	# /usr/include/X11/X.h:110
KeySym = XID 	# /usr/include/X11/X.h:111
KeyCode = c_ubyte 	# /usr/include/X11/X.h:113
None_ = 0 	# /usr/include/X11/X.h:120
ParentRelative = 1 	# /usr/include/X11/X.h:123
CopyFromParent = 0 	# /usr/include/X11/X.h:126
PointerWindow = 0 	# /usr/include/X11/X.h:131
InputFocus = 1 	# /usr/include/X11/X.h:132
PointerRoot = 1 	# /usr/include/X11/X.h:134
AnyPropertyType = 0 	# /usr/include/X11/X.h:136
AnyKey = 0 	# /usr/include/X11/X.h:138
AnyButton = 0 	# /usr/include/X11/X.h:140
AllTemporary = 0 	# /usr/include/X11/X.h:142
CurrentTime = 0 	# /usr/include/X11/X.h:144
NoSymbol = 0 	# /usr/include/X11/X.h:146
NoEventMask = 0 	# /usr/include/X11/X.h:155
KeyPressMask = 1 	# /usr/include/X11/X.h:156
KeyReleaseMask = 2 	# /usr/include/X11/X.h:157
ButtonPressMask = 4 	# /usr/include/X11/X.h:158
ButtonReleaseMask = 8 	# /usr/include/X11/X.h:159
EnterWindowMask = 16 	# /usr/include/X11/X.h:160
LeaveWindowMask = 32 	# /usr/include/X11/X.h:161
PointerMotionMask = 64 	# /usr/include/X11/X.h:162
PointerMotionHintMask = 128 	# /usr/include/X11/X.h:163
Button1MotionMask = 256 	# /usr/include/X11/X.h:164
Button2MotionMask = 512 	# /usr/include/X11/X.h:165
Button3MotionMask = 1024 	# /usr/include/X11/X.h:166
Button4MotionMask = 2048 	# /usr/include/X11/X.h:167
Button5MotionMask = 4096 	# /usr/include/X11/X.h:168
ButtonMotionMask = 8192 	# /usr/include/X11/X.h:169
KeymapStateMask = 16384 	# /usr/include/X11/X.h:170
ExposureMask = 32768 	# /usr/include/X11/X.h:171
VisibilityChangeMask = 65536 	# /usr/include/X11/X.h:172
StructureNotifyMask = 131072 	# /usr/include/X11/X.h:173
ResizeRedirectMask = 262144 	# /usr/include/X11/X.h:174
SubstructureNotifyMask = 524288 	# /usr/include/X11/X.h:175
SubstructureRedirectMask = 1048576 	# /usr/include/X11/X.h:176
FocusChangeMask = 2097152 	# /usr/include/X11/X.h:177
PropertyChangeMask = 4194304 	# /usr/include/X11/X.h:178
ColormapChangeMask = 8388608 	# /usr/include/X11/X.h:179
OwnerGrabButtonMask = 16777216 	# /usr/include/X11/X.h:180
KeyPress = 2 	# /usr/include/X11/X.h:186
KeyRelease = 3 	# /usr/include/X11/X.h:187
ButtonPress = 4 	# /usr/include/X11/X.h:188
ButtonRelease = 5 	# /usr/include/X11/X.h:189
MotionNotify = 6 	# /usr/include/X11/X.h:190
EnterNotify = 7 	# /usr/include/X11/X.h:191
LeaveNotify = 8 	# /usr/include/X11/X.h:192
FocusIn = 9 	# /usr/include/X11/X.h:193
FocusOut = 10 	# /usr/include/X11/X.h:194
KeymapNotify = 11 	# /usr/include/X11/X.h:195
Expose = 12 	# /usr/include/X11/X.h:196
GraphicsExpose = 13 	# /usr/include/X11/X.h:197
NoExpose = 14 	# /usr/include/X11/X.h:198
VisibilityNotify = 15 	# /usr/include/X11/X.h:199
CreateNotify = 16 	# /usr/include/X11/X.h:200
DestroyNotify = 17 	# /usr/include/X11/X.h:201
UnmapNotify = 18 	# /usr/include/X11/X.h:202
MapNotify = 19 	# /usr/include/X11/X.h:203
MapRequest = 20 	# /usr/include/X11/X.h:204
ReparentNotify = 21 	# /usr/include/X11/X.h:205
ConfigureNotify = 22 	# /usr/include/X11/X.h:206
ConfigureRequest = 23 	# /usr/include/X11/X.h:207
GravityNotify = 24 	# /usr/include/X11/X.h:208
ResizeRequest = 25 	# /usr/include/X11/X.h:209
CirculateNotify = 26 	# /usr/include/X11/X.h:210
CirculateRequest = 27 	# /usr/include/X11/X.h:211
PropertyNotify = 28 	# /usr/include/X11/X.h:212
SelectionClear = 29 	# /usr/include/X11/X.h:213
SelectionRequest = 30 	# /usr/include/X11/X.h:214
SelectionNotify = 31 	# /usr/include/X11/X.h:215
ColormapNotify = 32 	# /usr/include/X11/X.h:216
ClientMessage = 33 	# /usr/include/X11/X.h:217
MappingNotify = 34 	# /usr/include/X11/X.h:218
LASTEvent = 35 	# /usr/include/X11/X.h:219
ShiftMask = 1 	# /usr/include/X11/X.h:225
LockMask = 2 	# /usr/include/X11/X.h:226
ControlMask = 4 	# /usr/include/X11/X.h:227
Mod1Mask = 8 	# /usr/include/X11/X.h:228
Mod2Mask = 16 	# /usr/include/X11/X.h:229
Mod3Mask = 32 	# /usr/include/X11/X.h:230
Mod4Mask = 64 	# /usr/include/X11/X.h:231
Mod5Mask = 128 	# /usr/include/X11/X.h:232
ShiftMapIndex = 0 	# /usr/include/X11/X.h:237
LockMapIndex = 1 	# /usr/include/X11/X.h:238
ControlMapIndex = 2 	# /usr/include/X11/X.h:239
Mod1MapIndex = 3 	# /usr/include/X11/X.h:240
Mod2MapIndex = 4 	# /usr/include/X11/X.h:241
Mod3MapIndex = 5 	# /usr/include/X11/X.h:242
Mod4MapIndex = 6 	# /usr/include/X11/X.h:243
Mod5MapIndex = 7 	# /usr/include/X11/X.h:244
Button1Mask = 256 	# /usr/include/X11/X.h:250
Button2Mask = 512 	# /usr/include/X11/X.h:251
Button3Mask = 1024 	# /usr/include/X11/X.h:252
Button4Mask = 2048 	# /usr/include/X11/X.h:253
Button5Mask = 4096 	# /usr/include/X11/X.h:254
AnyModifier = 32768 	# /usr/include/X11/X.h:256
Button1 = 1 	# /usr/include/X11/X.h:263
Button2 = 2 	# /usr/include/X11/X.h:264
Button3 = 3 	# /usr/include/X11/X.h:265
Button4 = 4 	# /usr/include/X11/X.h:266
Button5 = 5 	# /usr/include/X11/X.h:267
NotifyNormal = 0 	# /usr/include/X11/X.h:271
NotifyGrab = 1 	# /usr/include/X11/X.h:272
NotifyUngrab = 2 	# /usr/include/X11/X.h:273
NotifyWhileGrabbed = 3 	# /usr/include/X11/X.h:274
NotifyHint = 1 	# /usr/include/X11/X.h:276
NotifyAncestor = 0 	# /usr/include/X11/X.h:280
NotifyVirtual = 1 	# /usr/include/X11/X.h:281
NotifyInferior = 2 	# /usr/include/X11/X.h:282
NotifyNonlinear = 3 	# /usr/include/X11/X.h:283
NotifyNonlinearVirtual = 4 	# /usr/include/X11/X.h:284
NotifyPointer = 5 	# /usr/include/X11/X.h:285
NotifyPointerRoot = 6 	# /usr/include/X11/X.h:286
NotifyDetailNone = 7 	# /usr/include/X11/X.h:287
VisibilityUnobscured = 0 	# /usr/include/X11/X.h:291
VisibilityPartiallyObscured = 1 	# /usr/include/X11/X.h:292
VisibilityFullyObscured = 2 	# /usr/include/X11/X.h:293
PlaceOnTop = 0 	# /usr/include/X11/X.h:297
PlaceOnBottom = 1 	# /usr/include/X11/X.h:298
FamilyInternet = 0 	# /usr/include/X11/X.h:302
FamilyDECnet = 1 	# /usr/include/X11/X.h:303
FamilyChaos = 2 	# /usr/include/X11/X.h:304
FamilyInternet6 = 6 	# /usr/include/X11/X.h:305
FamilyServerInterpreted = 5 	# /usr/include/X11/X.h:308
PropertyNewValue = 0 	# /usr/include/X11/X.h:312
PropertyDelete = 1 	# /usr/include/X11/X.h:313
ColormapUninstalled = 0 	# /usr/include/X11/X.h:317
ColormapInstalled = 1 	# /usr/include/X11/X.h:318
GrabModeSync = 0 	# /usr/include/X11/X.h:322
GrabModeAsync = 1 	# /usr/include/X11/X.h:323
GrabSuccess = 0 	# /usr/include/X11/X.h:327
AlreadyGrabbed = 1 	# /usr/include/X11/X.h:328
GrabInvalidTime = 2 	# /usr/include/X11/X.h:329
GrabNotViewable = 3 	# /usr/include/X11/X.h:330
GrabFrozen = 4 	# /usr/include/X11/X.h:331
AsyncPointer = 0 	# /usr/include/X11/X.h:335
SyncPointer = 1 	# /usr/include/X11/X.h:336
ReplayPointer = 2 	# /usr/include/X11/X.h:337
AsyncKeyboard = 3 	# /usr/include/X11/X.h:338
SyncKeyboard = 4 	# /usr/include/X11/X.h:339
ReplayKeyboard = 5 	# /usr/include/X11/X.h:340
AsyncBoth = 6 	# /usr/include/X11/X.h:341
SyncBoth = 7 	# /usr/include/X11/X.h:342
RevertToParent = 2 	# /usr/include/X11/X.h:348
Success = 0 	# /usr/include/X11/X.h:354
BadRequest = 1 	# /usr/include/X11/X.h:355
BadValue = 2 	# /usr/include/X11/X.h:356
BadWindow = 3 	# /usr/include/X11/X.h:357
BadPixmap = 4 	# /usr/include/X11/X.h:358
BadAtom = 5 	# /usr/include/X11/X.h:359
BadCursor = 6 	# /usr/include/X11/X.h:360
BadFont = 7 	# /usr/include/X11/X.h:361
BadMatch = 8 	# /usr/include/X11/X.h:362
BadDrawable = 9 	# /usr/include/X11/X.h:363
BadAccess = 10 	# /usr/include/X11/X.h:364
BadAlloc = 11 	# /usr/include/X11/X.h:373
BadColor = 12 	# /usr/include/X11/X.h:374
BadGC = 13 	# /usr/include/X11/X.h:375
BadIDChoice = 14 	# /usr/include/X11/X.h:376
BadName = 15 	# /usr/include/X11/X.h:377
BadLength = 16 	# /usr/include/X11/X.h:378
BadImplementation = 17 	# /usr/include/X11/X.h:379
FirstExtensionError = 128 	# /usr/include/X11/X.h:381
LastExtensionError = 255 	# /usr/include/X11/X.h:382
InputOutput = 1 	# /usr/include/X11/X.h:391
InputOnly = 2 	# /usr/include/X11/X.h:392
CWBackPixmap = 1 	# /usr/include/X11/X.h:396
CWBackPixel = 2 	# /usr/include/X11/X.h:397
CWBorderPixmap = 4 	# /usr/include/X11/X.h:398
CWBorderPixel = 8 	# /usr/include/X11/X.h:399
CWBitGravity = 16 	# /usr/include/X11/X.h:400
CWWinGravity = 32 	# /usr/include/X11/X.h:401
CWBackingStore = 64 	# /usr/include/X11/X.h:402
CWBackingPlanes = 128 	# /usr/include/X11/X.h:403
CWBackingPixel = 256 	# /usr/include/X11/X.h:404
CWOverrideRedirect = 512 	# /usr/include/X11/X.h:405
CWSaveUnder = 1024 	# /usr/include/X11/X.h:406
CWEventMask = 2048 	# /usr/include/X11/X.h:407
CWDontPropagate = 4096 	# /usr/include/X11/X.h:408
CWColormap = 8192 	# /usr/include/X11/X.h:409
CWCursor = 16384 	# /usr/include/X11/X.h:410
CWX = 1 	# /usr/include/X11/X.h:414
CWY = 2 	# /usr/include/X11/X.h:415
CWWidth = 4 	# /usr/include/X11/X.h:416
CWHeight = 8 	# /usr/include/X11/X.h:417
CWBorderWidth = 16 	# /usr/include/X11/X.h:418
CWSibling = 32 	# /usr/include/X11/X.h:419
CWStackMode = 64 	# /usr/include/X11/X.h:420
ForgetGravity = 0 	# /usr/include/X11/X.h:425
NorthWestGravity = 1 	# /usr/include/X11/X.h:426
NorthGravity = 2 	# /usr/include/X11/X.h:427
NorthEastGravity = 3 	# /usr/include/X11/X.h:428
WestGravity = 4 	# /usr/include/X11/X.h:429
CenterGravity = 5 	# /usr/include/X11/X.h:430
EastGravity = 6 	# /usr/include/X11/X.h:431
SouthWestGravity = 7 	# /usr/include/X11/X.h:432
SouthGravity = 8 	# /usr/include/X11/X.h:433
SouthEastGravity = 9 	# /usr/include/X11/X.h:434
StaticGravity = 10 	# /usr/include/X11/X.h:435
UnmapGravity = 0 	# /usr/include/X11/X.h:439
NotUseful = 0 	# /usr/include/X11/X.h:443
WhenMapped = 1 	# /usr/include/X11/X.h:444
Always = 2 	# /usr/include/X11/X.h:445
IsUnmapped = 0 	# /usr/include/X11/X.h:449
IsUnviewable = 1 	# /usr/include/X11/X.h:450
IsViewable = 2 	# /usr/include/X11/X.h:451
SetModeInsert = 0 	# /usr/include/X11/X.h:455
SetModeDelete = 1 	# /usr/include/X11/X.h:456
DestroyAll = 0 	# /usr/include/X11/X.h:460
RetainPermanent = 1 	# /usr/include/X11/X.h:461
RetainTemporary = 2 	# /usr/include/X11/X.h:462
Above = 0 	# /usr/include/X11/X.h:466
Below = 1 	# /usr/include/X11/X.h:467
TopIf = 2 	# /usr/include/X11/X.h:468
BottomIf = 3 	# /usr/include/X11/X.h:469
Opposite = 4 	# /usr/include/X11/X.h:470
RaiseLowest = 0 	# /usr/include/X11/X.h:474
LowerHighest = 1 	# /usr/include/X11/X.h:475
PropModeReplace = 0 	# /usr/include/X11/X.h:479
PropModePrepend = 1 	# /usr/include/X11/X.h:480
PropModeAppend = 2 	# /usr/include/X11/X.h:481
GXclear = 0 	# /usr/include/X11/X.h:489
GXand = 1 	# /usr/include/X11/X.h:490
GXandReverse = 2 	# /usr/include/X11/X.h:491
GXcopy = 3 	# /usr/include/X11/X.h:492
GXandInverted = 4 	# /usr/include/X11/X.h:493
GXnoop = 5 	# /usr/include/X11/X.h:494
GXxor = 6 	# /usr/include/X11/X.h:495
GXor = 7 	# /usr/include/X11/X.h:496
GXnor = 8 	# /usr/include/X11/X.h:497
GXequiv = 9 	# /usr/include/X11/X.h:498
GXinvert = 10 	# /usr/include/X11/X.h:499
GXorReverse = 11 	# /usr/include/X11/X.h:500
GXcopyInverted = 12 	# /usr/include/X11/X.h:501
GXorInverted = 13 	# /usr/include/X11/X.h:502
GXnand = 14 	# /usr/include/X11/X.h:503
GXset = 0.0 	# /usr/include/X11/X.h:504
LineSolid = 0 	# /usr/include/X11/X.h:508
LineOnOffDash = 1 	# /usr/include/X11/X.h:509
LineDoubleDash = 2 	# /usr/include/X11/X.h:510
CapNotLast = 0 	# /usr/include/X11/X.h:514
CapButt = 1 	# /usr/include/X11/X.h:515
CapRound = 2 	# /usr/include/X11/X.h:516
CapProjecting = 3 	# /usr/include/X11/X.h:517
JoinMiter = 0 	# /usr/include/X11/X.h:521
JoinRound = 1 	# /usr/include/X11/X.h:522
JoinBevel = 2 	# /usr/include/X11/X.h:523
FillSolid = 0 	# /usr/include/X11/X.h:527
FillTiled = 1 	# /usr/include/X11/X.h:528
FillStippled = 2 	# /usr/include/X11/X.h:529
FillOpaqueStippled = 3 	# /usr/include/X11/X.h:530
EvenOddRule = 0 	# /usr/include/X11/X.h:534
WindingRule = 1 	# /usr/include/X11/X.h:535
ClipByChildren = 0 	# /usr/include/X11/X.h:539
IncludeInferiors = 1 	# /usr/include/X11/X.h:540
Unsorted = 0 	# /usr/include/X11/X.h:544
YSorted = 1 	# /usr/include/X11/X.h:545
YXSorted = 2 	# /usr/include/X11/X.h:546
YXBanded = 3 	# /usr/include/X11/X.h:547
CoordModeOrigin = 0 	# /usr/include/X11/X.h:551
CoordModePrevious = 1 	# /usr/include/X11/X.h:552
Complex = 0 	# /usr/include/X11/X.h:556
Nonconvex = 1 	# /usr/include/X11/X.h:557
Convex = 2 	# /usr/include/X11/X.h:558
ArcChord = 0 	# /usr/include/X11/X.h:562
ArcPieSlice = 1 	# /usr/include/X11/X.h:563
GCFunction = 1 	# /usr/include/X11/X.h:568
GCPlaneMask = 2 	# /usr/include/X11/X.h:569
GCForeground = 4 	# /usr/include/X11/X.h:570
GCBackground = 8 	# /usr/include/X11/X.h:571
GCLineWidth = 16 	# /usr/include/X11/X.h:572
GCLineStyle = 32 	# /usr/include/X11/X.h:573
GCCapStyle = 64 	# /usr/include/X11/X.h:574
GCJoinStyle = 128 	# /usr/include/X11/X.h:575
GCFillStyle = 256 	# /usr/include/X11/X.h:576
GCFillRule = 512 	# /usr/include/X11/X.h:577
GCTile = 1024 	# /usr/include/X11/X.h:578
GCStipple = 2048 	# /usr/include/X11/X.h:579
GCTileStipXOrigin = 4096 	# /usr/include/X11/X.h:580
GCTileStipYOrigin = 8192 	# /usr/include/X11/X.h:581
GCFont = 16384 	# /usr/include/X11/X.h:582
GCSubwindowMode = 32768 	# /usr/include/X11/X.h:583
GCGraphicsExposures = 65536 	# /usr/include/X11/X.h:584
GCClipXOrigin = 131072 	# /usr/include/X11/X.h:585
GCClipYOrigin = 262144 	# /usr/include/X11/X.h:586
GCClipMask = 524288 	# /usr/include/X11/X.h:587
GCDashOffset = 1048576 	# /usr/include/X11/X.h:588
GCDashList = 2097152 	# /usr/include/X11/X.h:589
GCArcMode = 4194304 	# /usr/include/X11/X.h:590
GCLastBit = 22 	# /usr/include/X11/X.h:592
FontLeftToRight = 0 	# /usr/include/X11/X.h:599
FontRightToLeft = 1 	# /usr/include/X11/X.h:600
FontChange = 255 	# /usr/include/X11/X.h:602
XYBitmap = 0 	# /usr/include/X11/X.h:610
XYPixmap = 1 	# /usr/include/X11/X.h:611
ZPixmap = 2 	# /usr/include/X11/X.h:612
AllocNone = 0 	# /usr/include/X11/X.h:620
AllocAll = 1 	# /usr/include/X11/X.h:621
DoRed = 1 	# /usr/include/X11/X.h:626
DoGreen = 2 	# /usr/include/X11/X.h:627
DoBlue = 4 	# /usr/include/X11/X.h:628
CursorShape = 0 	# /usr/include/X11/X.h:636
TileShape = 1 	# /usr/include/X11/X.h:637
StippleShape = 2 	# /usr/include/X11/X.h:638
AutoRepeatModeOff = 0 	# /usr/include/X11/X.h:644
AutoRepeatModeOn = 1 	# /usr/include/X11/X.h:645
AutoRepeatModeDefault = 2 	# /usr/include/X11/X.h:646
LedModeOff = 0 	# /usr/include/X11/X.h:648
LedModeOn = 1 	# /usr/include/X11/X.h:649
KBKeyClickPercent = 1 	# /usr/include/X11/X.h:653
KBBellPercent = 2 	# /usr/include/X11/X.h:654
KBBellPitch = 4 	# /usr/include/X11/X.h:655
KBBellDuration = 8 	# /usr/include/X11/X.h:656
KBLed = 16 	# /usr/include/X11/X.h:657
KBLedMode = 32 	# /usr/include/X11/X.h:658
KBKey = 64 	# /usr/include/X11/X.h:659
KBAutoRepeatMode = 128 	# /usr/include/X11/X.h:660
MappingSuccess = 0 	# /usr/include/X11/X.h:662
MappingBusy = 1 	# /usr/include/X11/X.h:663
MappingFailed = 2 	# /usr/include/X11/X.h:664
MappingModifier = 0 	# /usr/include/X11/X.h:666
MappingKeyboard = 1 	# /usr/include/X11/X.h:667
MappingPointer = 2 	# /usr/include/X11/X.h:668
DontPreferBlanking = 0 	# /usr/include/X11/X.h:674
PreferBlanking = 1 	# /usr/include/X11/X.h:675
DefaultBlanking = 2 	# /usr/include/X11/X.h:676
DisableScreenSaver = 0 	# /usr/include/X11/X.h:678
DisableScreenInterval = 0 	# /usr/include/X11/X.h:679
DontAllowExposures = 0 	# /usr/include/X11/X.h:681
AllowExposures = 1 	# /usr/include/X11/X.h:682
DefaultExposures = 2 	# /usr/include/X11/X.h:683
ScreenSaverReset = 0 	# /usr/include/X11/X.h:687
ScreenSaverActive = 1 	# /usr/include/X11/X.h:688
HostInsert = 0 	# /usr/include/X11/X.h:696
HostDelete = 1 	# /usr/include/X11/X.h:697
EnableAccess = 1 	# /usr/include/X11/X.h:701
DisableAccess = 0 	# /usr/include/X11/X.h:702
StaticGray = 0 	# /usr/include/X11/X.h:708
GrayScale = 1 	# /usr/include/X11/X.h:709
StaticColor = 2 	# /usr/include/X11/X.h:710
PseudoColor = 3 	# /usr/include/X11/X.h:711
TrueColor = 4 	# /usr/include/X11/X.h:712
DirectColor = 5 	# /usr/include/X11/X.h:713
LSBFirst = 0 	# /usr/include/X11/X.h:718
MSBFirst = 1 	# /usr/include/X11/X.h:719
# /usr/include/X11/Xlib.h:94
_Xmblen = _lib._Xmblen
_Xmblen.restype = c_int
_Xmblen.argtypes = [c_char_p, c_int]

X_HAVE_UTF8_STRING = 1 	# /usr/include/X11/Xlib.h:106
XPointer = pyglet.gl.glx.XPointer
Bool = 0 	# /usr/include/X11/Xlib.h:110
Status = 0 	# /usr/include/X11/Xlib.h:111
True = 1 	# /usr/include/X11/Xlib.h:112
False = 0 	# /usr/include/X11/Xlib.h:113
QueuedAlready = 0 	# /usr/include/X11/Xlib.h:115
QueuedAfterReading = 1 	# /usr/include/X11/Xlib.h:116
QueuedAfterFlush = 2 	# /usr/include/X11/Xlib.h:117
XExtData = pyglet.gl.glx.XExtData
class struct_anon_8(Structure):
    __slots__ = [
        'extension',
        'major_opcode',
        'first_event',
        'first_error',
    ]
struct_anon_8._fields_ = [
    ('extension', c_int),
    ('major_opcode', c_int),
    ('first_event', c_int),
    ('first_error', c_int),
]

XExtCodes = struct_anon_8 	# /usr/include/X11/Xlib.h:197
class struct_anon_9(Structure):
    __slots__ = [
        'depth',
        'bits_per_pixel',
        'scanline_pad',
    ]
struct_anon_9._fields_ = [
    ('depth', c_int),
    ('bits_per_pixel', c_int),
    ('scanline_pad', c_int),
]

XPixmapFormatValues = struct_anon_9 	# /usr/include/X11/Xlib.h:207
class struct_anon_10(Structure):
    __slots__ = [
        'function',
        'plane_mask',
        'foreground',
        'background',
        'line_width',
        'line_style',
        'cap_style',
        'join_style',
        'fill_style',
        'fill_rule',
        'arc_mode',
        'tile',
        'stipple',
        'ts_x_origin',
        'ts_y_origin',
        'font',
        'subwindow_mode',
        'graphics_exposures',
        'clip_x_origin',
        'clip_y_origin',
        'clip_mask',
        'dash_offset',
        'dashes',
    ]
struct_anon_10._fields_ = [
    ('function', c_int),
    ('plane_mask', c_ulong),
    ('foreground', c_ulong),
    ('background', c_ulong),
    ('line_width', c_int),
    ('line_style', c_int),
    ('cap_style', c_int),
    ('join_style', c_int),
    ('fill_style', c_int),
    ('fill_rule', c_int),
    ('arc_mode', c_int),
    ('tile', Pixmap),
    ('stipple', Pixmap),
    ('ts_x_origin', c_int),
    ('ts_y_origin', c_int),
    ('font', Font),
    ('subwindow_mode', c_int),
    ('graphics_exposures', c_int),
    ('clip_x_origin', c_int),
    ('clip_y_origin', c_int),
    ('clip_mask', Pixmap),
    ('dash_offset', c_int),
    ('dashes', c_char),
]

XGCValues = struct_anon_10 	# /usr/include/X11/Xlib.h:239
class struct__XGC(Structure):
    __slots__ = [
    ]
struct__XGC._fields_ = [
    ('_opaque_struct', c_int)
]

GC = POINTER(struct__XGC) 	# /usr/include/X11/Xlib.h:254
Visual = pyglet.gl.glx.Visual
class struct_anon_12(Structure):
    __slots__ = [
        'depth',
        'nvisuals',
        'visuals',
    ]
struct_anon_12._fields_ = [
    ('depth', c_int),
    ('nvisuals', c_int),
    ('visuals', POINTER(Visual)),
]

Depth = struct_anon_12 	# /usr/include/X11/Xlib.h:279
class struct_anon_13(Structure):
    __slots__ = [
        'ext_data',
        'display',
        'root',
        'width',
        'height',
        'mwidth',
        'mheight',
        'ndepths',
        'depths',
        'root_depth',
        'root_visual',
        'default_gc',
        'cmap',
        'white_pixel',
        'black_pixel',
        'max_maps',
        'min_maps',
        'backing_store',
        'save_unders',
        'root_input_mask',
    ]
class struct__XDisplay(Structure):
    __slots__ = [
    ]
struct__XDisplay._fields_ = [
    ('_opaque_struct', c_int)
]

struct_anon_13._fields_ = [
    ('ext_data', POINTER(XExtData)),
    ('display', POINTER(struct__XDisplay)),
    ('root', Window),
    ('width', c_int),
    ('height', c_int),
    ('mwidth', c_int),
    ('mheight', c_int),
    ('ndepths', c_int),
    ('depths', POINTER(Depth)),
    ('root_depth', c_int),
    ('root_visual', POINTER(Visual)),
    ('default_gc', GC),
    ('cmap', Colormap),
    ('white_pixel', c_ulong),
    ('black_pixel', c_ulong),
    ('max_maps', c_int),
    ('min_maps', c_int),
    ('backing_store', c_int),
    ('save_unders', c_int),
    ('root_input_mask', c_long),
]

Screen = struct_anon_13 	# /usr/include/X11/Xlib.h:307
class struct_anon_14(Structure):
    __slots__ = [
        'ext_data',
        'depth',
        'bits_per_pixel',
        'scanline_pad',
    ]
struct_anon_14._fields_ = [
    ('ext_data', POINTER(XExtData)),
    ('depth', c_int),
    ('bits_per_pixel', c_int),
    ('scanline_pad', c_int),
]

ScreenFormat = struct_anon_14 	# /usr/include/X11/Xlib.h:317
class struct_anon_15(Structure):
    __slots__ = [
        'background_pixmap',
        'background_pixel',
        'border_pixmap',
        'border_pixel',
        'bit_gravity',
        'win_gravity',
        'backing_store',
        'backing_planes',
        'backing_pixel',
        'save_under',
        'event_mask',
        'do_not_propagate_mask',
        'override_redirect',
        'colormap',
        'cursor',
    ]
struct_anon_15._fields_ = [
    ('background_pixmap', Pixmap),
    ('background_pixel', c_ulong),
    ('border_pixmap', Pixmap),
    ('border_pixel', c_ulong),
    ('bit_gravity', c_int),
    ('win_gravity', c_int),
    ('backing_store', c_int),
    ('backing_planes', c_ulong),
    ('backing_pixel', c_ulong),
    ('save_under', c_int),
    ('event_mask', c_long),
    ('do_not_propagate_mask', c_long),
    ('override_redirect', c_int),
    ('colormap', Colormap),
    ('cursor', Cursor),
]

XSetWindowAttributes = struct_anon_15 	# /usr/include/X11/Xlib.h:338
class struct_anon_16(Structure):
    __slots__ = [
        'x',
        'y',
        'width',
        'height',
        'border_width',
        'depth',
        'visual',
        'root',
        'class',
        'bit_gravity',
        'win_gravity',
        'backing_store',
        'backing_planes',
        'backing_pixel',
        'save_under',
        'colormap',
        'map_installed',
        'map_state',
        'all_event_masks',
        'your_event_mask',
        'do_not_propagate_mask',
        'override_redirect',
        'screen',
    ]
struct_anon_16._fields_ = [
    ('x', c_int),
    ('y', c_int),
    ('width', c_int),
    ('height', c_int),
    ('border_width', c_int),
    ('depth', c_int),
    ('visual', POINTER(Visual)),
    ('root', Window),
    ('class', c_int),
    ('bit_gravity', c_int),
    ('win_gravity', c_int),
    ('backing_store', c_int),
    ('backing_planes', c_ulong),
    ('backing_pixel', c_ulong),
    ('save_under', c_int),
    ('colormap', Colormap),
    ('map_installed', c_int),
    ('map_state', c_int),
    ('all_event_masks', c_long),
    ('your_event_mask', c_long),
    ('do_not_propagate_mask', c_long),
    ('override_redirect', c_int),
    ('screen', POINTER(Screen)),
]

XWindowAttributes = struct_anon_16 	# /usr/include/X11/Xlib.h:366
class struct_anon_17(Structure):
    __slots__ = [
        'family',
        'length',
        'address',
    ]
struct_anon_17._fields_ = [
    ('family', c_int),
    ('length', c_int),
    ('address', c_char_p),
]

XHostAddress = struct_anon_17 	# /usr/include/X11/Xlib.h:377
class struct_anon_18(Structure):
    __slots__ = [
        'typelength',
        'valuelength',
        'type',
        'value',
    ]
struct_anon_18._fields_ = [
    ('typelength', c_int),
    ('valuelength', c_int),
    ('type', c_char_p),
    ('value', c_char_p),
]

XServerInterpretedAddress = struct_anon_18 	# /usr/include/X11/Xlib.h:387
class struct__XImage(Structure):
    __slots__ = [
        'width',
        'height',
        'xoffset',
        'format',
        'data',
        'byte_order',
        'bitmap_unit',
        'bitmap_bit_order',
        'bitmap_pad',
        'depth',
        'bytes_per_line',
        'bits_per_pixel',
        'red_mask',
        'green_mask',
        'blue_mask',
        'obdata',
        'f',
    ]
class struct_funcs(Structure):
    __slots__ = [
        'create_image',
        'destroy_image',
        'get_pixel',
        'put_pixel',
        'sub_image',
        'add_pixel',
    ]
struct_funcs._fields_ = [
    ('create_image', POINTER(CFUNCTYPE(POINTER(struct__XImage), POINTER(struct__XDisplay), POINTER(Visual), c_uint, c_int, c_int, c_char_p, c_uint, c_uint, c_int, c_int))),
    ('destroy_image', POINTER(CFUNCTYPE(c_int, POINTER(struct__XImage)))),
    ('get_pixel', POINTER(CFUNCTYPE(c_ulong, POINTER(struct__XImage), c_int, c_int))),
    ('put_pixel', POINTER(CFUNCTYPE(c_int, POINTER(struct__XImage), c_int, c_int, c_ulong))),
    ('sub_image', POINTER(CFUNCTYPE(POINTER(struct__XImage), POINTER(struct__XImage), c_int, c_int, c_uint, c_uint))),
    ('add_pixel', POINTER(CFUNCTYPE(c_int, POINTER(struct__XImage), c_long))),
]

struct__XImage._fields_ = [
    ('width', c_int),
    ('height', c_int),
    ('xoffset', c_int),
    ('format', c_int),
    ('data', c_char_p),
    ('byte_order', c_int),
    ('bitmap_unit', c_int),
    ('bitmap_bit_order', c_int),
    ('bitmap_pad', c_int),
    ('depth', c_int),
    ('bytes_per_line', c_int),
    ('bits_per_pixel', c_int),
    ('red_mask', c_ulong),
    ('green_mask', c_ulong),
    ('blue_mask', c_ulong),
    ('obdata', XPointer),
    ('f', struct_funcs),
]

XImage = struct__XImage 	# /usr/include/X11/Xlib.h:426
class struct_anon_19(Structure):
    __slots__ = [
        'x',
        'y',
        'width',
        'height',
        'border_width',
        'sibling',
        'stack_mode',
    ]
struct_anon_19._fields_ = [
    ('x', c_int),
    ('y', c_int),
    ('width', c_int),
    ('height', c_int),
    ('border_width', c_int),
    ('sibling', Window),
    ('stack_mode', c_int),
]

XWindowChanges = struct_anon_19 	# /usr/include/X11/Xlib.h:437
class struct_anon_20(Structure):
    __slots__ = [
        'pixel',
        'red',
        'green',
        'blue',
        'flags',
        'pad',
    ]
struct_anon_20._fields_ = [
    ('pixel', c_ulong),
    ('red', c_ushort),
    ('green', c_ushort),
    ('blue', c_ushort),
    ('flags', c_char),
    ('pad', c_char),
]

XColor = struct_anon_20 	# /usr/include/X11/Xlib.h:447
class struct_anon_21(Structure):
    __slots__ = [
        'x1',
        'y1',
        'x2',
        'y2',
    ]
struct_anon_21._fields_ = [
    ('x1', c_short),
    ('y1', c_short),
    ('x2', c_short),
    ('y2', c_short),
]

XSegment = struct_anon_21 	# /usr/include/X11/Xlib.h:456
class struct_anon_22(Structure):
    __slots__ = [
        'x',
        'y',
    ]
struct_anon_22._fields_ = [
    ('x', c_short),
    ('y', c_short),
]

XPoint = struct_anon_22 	# /usr/include/X11/Xlib.h:460
class struct_anon_23(Structure):
    __slots__ = [
        'x',
        'y',
        'width',
        'height',
    ]
struct_anon_23._fields_ = [
    ('x', c_short),
    ('y', c_short),
    ('width', c_ushort),
    ('height', c_ushort),
]

XRectangle = struct_anon_23 	# /usr/include/X11/Xlib.h:465
class struct_anon_24(Structure):
    __slots__ = [
        'x',
        'y',
        'width',
        'height',
        'angle1',
        'angle2',
    ]
struct_anon_24._fields_ = [
    ('x', c_short),
    ('y', c_short),
    ('width', c_ushort),
    ('height', c_ushort),
    ('angle1', c_short),
    ('angle2', c_short),
]

XArc = struct_anon_24 	# /usr/include/X11/Xlib.h:471
class struct_anon_25(Structure):
    __slots__ = [
        'key_click_percent',
        'bell_percent',
        'bell_pitch',
        'bell_duration',
        'led',
        'led_mode',
        'key',
        'auto_repeat_mode',
    ]
struct_anon_25._fields_ = [
    ('key_click_percent', c_int),
    ('bell_percent', c_int),
    ('bell_pitch', c_int),
    ('bell_duration', c_int),
    ('led', c_int),
    ('led_mode', c_int),
    ('key', c_int),
    ('auto_repeat_mode', c_int),
]

XKeyboardControl = struct_anon_25 	# /usr/include/X11/Xlib.h:485
class struct_anon_26(Structure):
    __slots__ = [
        'key_click_percent',
        'bell_percent',
        'bell_pitch',
        'bell_duration',
        'led_mask',
        'global_auto_repeat',
        'auto_repeats',
    ]
struct_anon_26._fields_ = [
    ('key_click_percent', c_int),
    ('bell_percent', c_int),
    ('bell_pitch', c_uint),
    ('bell_duration', c_uint),
    ('led_mask', c_ulong),
    ('global_auto_repeat', c_int),
    ('auto_repeats', c_char * 32),
]

XKeyboardState = struct_anon_26 	# /usr/include/X11/Xlib.h:496
class struct_anon_27(Structure):
    __slots__ = [
        'time',
        'x',
        'y',
    ]
struct_anon_27._fields_ = [
    ('time', Time),
    ('x', c_short),
    ('y', c_short),
]

XTimeCoord = struct_anon_27 	# /usr/include/X11/Xlib.h:503
class struct_anon_28(Structure):
    __slots__ = [
        'max_keypermod',
        'modifiermap',
    ]
struct_anon_28._fields_ = [
    ('max_keypermod', c_int),
    ('modifiermap', POINTER(KeyCode)),
]

XModifierKeymap = struct_anon_28 	# /usr/include/X11/Xlib.h:510
Display = pyglet.gl.glx.Display
class struct_anon_29(Structure):
    __slots__ = [
        'ext_data',
        'private1',
        'fd',
        'private2',
        'proto_major_version',
        'proto_minor_version',
        'vendor',
        'private3',
        'private4',
        'private5',
        'private6',
        'resource_alloc',
        'byte_order',
        'bitmap_unit',
        'bitmap_pad',
        'bitmap_bit_order',
        'nformats',
        'pixmap_format',
        'private8',
        'release',
        'private9',
        'private10',
        'qlen',
        'last_request_read',
        'request',
        'private11',
        'private12',
        'private13',
        'private14',
        'max_request_size',
        'db',
        'private15',
        'display_name',
        'default_screen',
        'nscreens',
        'screens',
        'motion_buffer',
        'private16',
        'min_keycode',
        'max_keycode',
        'private17',
        'private18',
        'private19',
        'xdefaults',
    ]
class struct__XPrivate(Structure):
    __slots__ = [
    ]
struct__XPrivate._fields_ = [
    ('_opaque_struct', c_int)
]

class struct__XrmHashBucketRec(Structure):
    __slots__ = [
    ]
struct__XrmHashBucketRec._fields_ = [
    ('_opaque_struct', c_int)
]

struct_anon_29._fields_ = [
    ('ext_data', POINTER(XExtData)),
    ('private1', POINTER(struct__XPrivate)),
    ('fd', c_int),
    ('private2', c_int),
    ('proto_major_version', c_int),
    ('proto_minor_version', c_int),
    ('vendor', c_char_p),
    ('private3', XID),
    ('private4', XID),
    ('private5', XID),
    ('private6', c_int),
    ('resource_alloc', POINTER(CFUNCTYPE(XID, POINTER(struct__XDisplay)))),
    ('byte_order', c_int),
    ('bitmap_unit', c_int),
    ('bitmap_pad', c_int),
    ('bitmap_bit_order', c_int),
    ('nformats', c_int),
    ('pixmap_format', POINTER(ScreenFormat)),
    ('private8', c_int),
    ('release', c_int),
    ('private9', POINTER(struct__XPrivate)),
    ('private10', POINTER(struct__XPrivate)),
    ('qlen', c_int),
    ('last_request_read', c_ulong),
    ('request', c_ulong),
    ('private11', XPointer),
    ('private12', XPointer),
    ('private13', XPointer),
    ('private14', XPointer),
    ('max_request_size', c_uint),
    ('db', POINTER(struct__XrmHashBucketRec)),
    ('private15', POINTER(CFUNCTYPE(c_int, POINTER(struct__XDisplay)))),
    ('display_name', c_char_p),
    ('default_screen', c_int),
    ('nscreens', c_int),
    ('screens', POINTER(Screen)),
    ('motion_buffer', c_ulong),
    ('private16', c_ulong),
    ('min_keycode', c_int),
    ('max_keycode', c_int),
    ('private17', XPointer),
    ('private18', XPointer),
    ('private19', c_int),
    ('xdefaults', c_char_p),
]

_XPrivDisplay = POINTER(struct_anon_29) 	# /usr/include/X11/Xlib.h:582
class struct_anon_30(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'window',
        'root',
        'subwindow',
        'time',
        'x',
        'y',
        'x_root',
        'y_root',
        'state',
        'keycode',
        'same_screen',
    ]
struct_anon_30._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('window', Window),
    ('root', Window),
    ('subwindow', Window),
    ('time', Time),
    ('x', c_int),
    ('y', c_int),
    ('x_root', c_int),
    ('y_root', c_int),
    ('state', c_uint),
    ('keycode', c_uint),
    ('same_screen', c_int),
]

XKeyEvent = struct_anon_30 	# /usr/include/X11/Xlib.h:603
XKeyPressedEvent = XKeyEvent 	# /usr/include/X11/Xlib.h:604
XKeyReleasedEvent = XKeyEvent 	# /usr/include/X11/Xlib.h:605
class struct_anon_31(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'window',
        'root',
        'subwindow',
        'time',
        'x',
        'y',
        'x_root',
        'y_root',
        'state',
        'button',
        'same_screen',
    ]
struct_anon_31._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('window', Window),
    ('root', Window),
    ('subwindow', Window),
    ('time', Time),
    ('x', c_int),
    ('y', c_int),
    ('x_root', c_int),
    ('y_root', c_int),
    ('state', c_uint),
    ('button', c_uint),
    ('same_screen', c_int),
]

XButtonEvent = struct_anon_31 	# /usr/include/X11/Xlib.h:621
XButtonPressedEvent = XButtonEvent 	# /usr/include/X11/Xlib.h:622
XButtonReleasedEvent = XButtonEvent 	# /usr/include/X11/Xlib.h:623
class struct_anon_32(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'window',
        'root',
        'subwindow',
        'time',
        'x',
        'y',
        'x_root',
        'y_root',
        'state',
        'is_hint',
        'same_screen',
    ]
struct_anon_32._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('window', Window),
    ('root', Window),
    ('subwindow', Window),
    ('time', Time),
    ('x', c_int),
    ('y', c_int),
    ('x_root', c_int),
    ('y_root', c_int),
    ('state', c_uint),
    ('is_hint', c_char),
    ('same_screen', c_int),
]

XMotionEvent = struct_anon_32 	# /usr/include/X11/Xlib.h:639
XPointerMovedEvent = XMotionEvent 	# /usr/include/X11/Xlib.h:640
class struct_anon_33(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'window',
        'root',
        'subwindow',
        'time',
        'x',
        'y',
        'x_root',
        'y_root',
        'mode',
        'detail',
        'same_screen',
        'focus',
        'state',
    ]
struct_anon_33._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('window', Window),
    ('root', Window),
    ('subwindow', Window),
    ('time', Time),
    ('x', c_int),
    ('y', c_int),
    ('x_root', c_int),
    ('y_root', c_int),
    ('mode', c_int),
    ('detail', c_int),
    ('same_screen', c_int),
    ('focus', c_int),
    ('state', c_uint),
]

XCrossingEvent = struct_anon_33 	# /usr/include/X11/Xlib.h:662
XEnterWindowEvent = XCrossingEvent 	# /usr/include/X11/Xlib.h:663
XLeaveWindowEvent = XCrossingEvent 	# /usr/include/X11/Xlib.h:664
class struct_anon_34(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'window',
        'mode',
        'detail',
    ]
struct_anon_34._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('window', Window),
    ('mode', c_int),
    ('detail', c_int),
]

XFocusChangeEvent = struct_anon_34 	# /usr/include/X11/Xlib.h:680
XFocusInEvent = XFocusChangeEvent 	# /usr/include/X11/Xlib.h:681
XFocusOutEvent = XFocusChangeEvent 	# /usr/include/X11/Xlib.h:682
class struct_anon_35(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'window',
        'key_vector',
    ]
struct_anon_35._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('window', Window),
    ('key_vector', c_char * 32),
]

XKeymapEvent = struct_anon_35 	# /usr/include/X11/Xlib.h:692
class struct_anon_36(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'window',
        'x',
        'y',
        'width',
        'height',
        'count',
    ]
struct_anon_36._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('window', Window),
    ('x', c_int),
    ('y', c_int),
    ('width', c_int),
    ('height', c_int),
    ('count', c_int),
]

XExposeEvent = struct_anon_36 	# /usr/include/X11/Xlib.h:703
class struct_anon_37(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'drawable',
        'x',
        'y',
        'width',
        'height',
        'count',
        'major_code',
        'minor_code',
    ]
struct_anon_37._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('drawable', Drawable),
    ('x', c_int),
    ('y', c_int),
    ('width', c_int),
    ('height', c_int),
    ('count', c_int),
    ('major_code', c_int),
    ('minor_code', c_int),
]

XGraphicsExposeEvent = struct_anon_37 	# /usr/include/X11/Xlib.h:716
class struct_anon_38(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'drawable',
        'major_code',
        'minor_code',
    ]
struct_anon_38._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('drawable', Drawable),
    ('major_code', c_int),
    ('minor_code', c_int),
]

XNoExposeEvent = struct_anon_38 	# /usr/include/X11/Xlib.h:726
class struct_anon_39(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'window',
        'state',
    ]
struct_anon_39._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('window', Window),
    ('state', c_int),
]

XVisibilityEvent = struct_anon_39 	# /usr/include/X11/Xlib.h:735
class struct_anon_40(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'parent',
        'window',
        'x',
        'y',
        'width',
        'height',
        'border_width',
        'override_redirect',
    ]
struct_anon_40._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('parent', Window),
    ('window', Window),
    ('x', c_int),
    ('y', c_int),
    ('width', c_int),
    ('height', c_int),
    ('border_width', c_int),
    ('override_redirect', c_int),
]

XCreateWindowEvent = struct_anon_40 	# /usr/include/X11/Xlib.h:748
class struct_anon_41(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'event',
        'window',
    ]
struct_anon_41._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('event', Window),
    ('window', Window),
]

XDestroyWindowEvent = struct_anon_41 	# /usr/include/X11/Xlib.h:757
class struct_anon_42(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'event',
        'window',
        'from_configure',
    ]
struct_anon_42._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('event', Window),
    ('window', Window),
    ('from_configure', c_int),
]

XUnmapEvent = struct_anon_42 	# /usr/include/X11/Xlib.h:767
class struct_anon_43(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'event',
        'window',
        'override_redirect',
    ]
struct_anon_43._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('event', Window),
    ('window', Window),
    ('override_redirect', c_int),
]

XMapEvent = struct_anon_43 	# /usr/include/X11/Xlib.h:777
class struct_anon_44(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'parent',
        'window',
    ]
struct_anon_44._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('parent', Window),
    ('window', Window),
]

XMapRequestEvent = struct_anon_44 	# /usr/include/X11/Xlib.h:786
class struct_anon_45(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'event',
        'window',
        'parent',
        'x',
        'y',
        'override_redirect',
    ]
struct_anon_45._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('event', Window),
    ('window', Window),
    ('parent', Window),
    ('x', c_int),
    ('y', c_int),
    ('override_redirect', c_int),
]

XReparentEvent = struct_anon_45 	# /usr/include/X11/Xlib.h:798
class struct_anon_46(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'event',
        'window',
        'x',
        'y',
        'width',
        'height',
        'border_width',
        'above',
        'override_redirect',
    ]
struct_anon_46._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('event', Window),
    ('window', Window),
    ('x', c_int),
    ('y', c_int),
    ('width', c_int),
    ('height', c_int),
    ('border_width', c_int),
    ('above', Window),
    ('override_redirect', c_int),
]

XConfigureEvent = struct_anon_46 	# /usr/include/X11/Xlib.h:812
class struct_anon_47(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'event',
        'window',
        'x',
        'y',
    ]
struct_anon_47._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('event', Window),
    ('window', Window),
    ('x', c_int),
    ('y', c_int),
]

XGravityEvent = struct_anon_47 	# /usr/include/X11/Xlib.h:822
class struct_anon_48(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'window',
        'width',
        'height',
    ]
struct_anon_48._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('window', Window),
    ('width', c_int),
    ('height', c_int),
]

XResizeRequestEvent = struct_anon_48 	# /usr/include/X11/Xlib.h:831
class struct_anon_49(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'parent',
        'window',
        'x',
        'y',
        'width',
        'height',
        'border_width',
        'above',
        'detail',
        'value_mask',
    ]
struct_anon_49._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('parent', Window),
    ('window', Window),
    ('x', c_int),
    ('y', c_int),
    ('width', c_int),
    ('height', c_int),
    ('border_width', c_int),
    ('above', Window),
    ('detail', c_int),
    ('value_mask', c_ulong),
]

XConfigureRequestEvent = struct_anon_49 	# /usr/include/X11/Xlib.h:846
class struct_anon_50(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'event',
        'window',
        'place',
    ]
struct_anon_50._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('event', Window),
    ('window', Window),
    ('place', c_int),
]

XCirculateEvent = struct_anon_50 	# /usr/include/X11/Xlib.h:856
class struct_anon_51(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'parent',
        'window',
        'place',
    ]
struct_anon_51._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('parent', Window),
    ('window', Window),
    ('place', c_int),
]

XCirculateRequestEvent = struct_anon_51 	# /usr/include/X11/Xlib.h:866
class struct_anon_52(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'window',
        'atom',
        'time',
        'state',
    ]
struct_anon_52._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('window', Window),
    ('atom', Atom),
    ('time', Time),
    ('state', c_int),
]

XPropertyEvent = struct_anon_52 	# /usr/include/X11/Xlib.h:877
class struct_anon_53(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'window',
        'selection',
        'time',
    ]
struct_anon_53._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('window', Window),
    ('selection', Atom),
    ('time', Time),
]

XSelectionClearEvent = struct_anon_53 	# /usr/include/X11/Xlib.h:887
class struct_anon_54(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'owner',
        'requestor',
        'selection',
        'target',
        'property',
        'time',
    ]
struct_anon_54._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('owner', Window),
    ('requestor', Window),
    ('selection', Atom),
    ('target', Atom),
    ('property', Atom),
    ('time', Time),
]

XSelectionRequestEvent = struct_anon_54 	# /usr/include/X11/Xlib.h:900
class struct_anon_55(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'requestor',
        'selection',
        'target',
        'property',
        'time',
    ]
struct_anon_55._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('requestor', Window),
    ('selection', Atom),
    ('target', Atom),
    ('property', Atom),
    ('time', Time),
]

XSelectionEvent = struct_anon_55 	# /usr/include/X11/Xlib.h:912
class struct_anon_56(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'window',
        'colormap',
        'new',
        'state',
    ]
struct_anon_56._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('window', Window),
    ('colormap', Colormap),
    ('new', c_int),
    ('state', c_int),
]

XColormapEvent = struct_anon_56 	# /usr/include/X11/Xlib.h:927
class struct_anon_57(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'window',
        'message_type',
        'format',
        'data',
    ]
class struct_anon_58(Union):
    __slots__ = [
        'b',
        's',
        'l',
    ]
struct_anon_58._fields_ = [
    ('b', c_char * 20),
    ('s', c_short * 10),
    ('l', c_long * 5),
]

struct_anon_57._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('window', Window),
    ('message_type', Atom),
    ('format', c_int),
    ('data', struct_anon_58),
]

XClientMessageEvent = struct_anon_57 	# /usr/include/X11/Xlib.h:942
class struct_anon_59(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'window',
        'request',
        'first_keycode',
        'count',
    ]
struct_anon_59._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('window', Window),
    ('request', c_int),
    ('first_keycode', c_int),
    ('count', c_int),
]

XMappingEvent = struct_anon_59 	# /usr/include/X11/Xlib.h:954
class struct_anon_60(Structure):
    __slots__ = [
        'type',
        'display',
        'resourceid',
        'serial',
        'error_code',
        'request_code',
        'minor_code',
    ]
struct_anon_60._fields_ = [
    ('type', c_int),
    ('display', POINTER(Display)),
    ('resourceid', XID),
    ('serial', c_ulong),
    ('error_code', c_ubyte),
    ('request_code', c_ubyte),
    ('minor_code', c_ubyte),
]

XErrorEvent = struct_anon_60 	# /usr/include/X11/Xlib.h:964
class struct_anon_61(Structure):
    __slots__ = [
        'type',
        'serial',
        'send_event',
        'display',
        'window',
    ]
struct_anon_61._fields_ = [
    ('type', c_int),
    ('serial', c_ulong),
    ('send_event', c_int),
    ('display', POINTER(Display)),
    ('window', Window),
]

XAnyEvent = struct_anon_61 	# /usr/include/X11/Xlib.h:972
class struct__XEvent(Union):
    __slots__ = [
        'type',
        'xany',
        'xkey',
        'xbutton',
        'xmotion',
        'xcrossing',
        'xfocus',
        'xexpose',
        'xgraphicsexpose',
        'xnoexpose',
        'xvisibility',
        'xcreatewindow',
        'xdestroywindow',
        'xunmap',
        'xmap',
        'xmaprequest',
        'xreparent',
        'xconfigure',
        'xgravity',
        'xresizerequest',
        'xconfigurerequest',
        'xcirculate',
        'xcirculaterequest',
        'xproperty',
        'xselectionclear',
        'xselectionrequest',
        'xselection',
        'xcolormap',
        'xclient',
        'xmapping',
        'xerror',
        'xkeymap',
        'pad',
    ]
struct__XEvent._fields_ = [
    ('type', c_int),
    ('xany', XAnyEvent),
    ('xkey', XKeyEvent),
    ('xbutton', XButtonEvent),
    ('xmotion', XMotionEvent),
    ('xcrossing', XCrossingEvent),
    ('xfocus', XFocusChangeEvent),
    ('xexpose', XExposeEvent),
    ('xgraphicsexpose', XGraphicsExposeEvent),
    ('xnoexpose', XNoExposeEvent),
    ('xvisibility', XVisibilityEvent),
    ('xcreatewindow', XCreateWindowEvent),
    ('xdestroywindow', XDestroyWindowEvent),
    ('xunmap', XUnmapEvent),
    ('xmap', XMapEvent),
    ('xmaprequest', XMapRequestEvent),
    ('xreparent', XReparentEvent),
    ('xconfigure', XConfigureEvent),
    ('xgravity', XGravityEvent),
    ('xresizerequest', XResizeRequestEvent),
    ('xconfigurerequest', XConfigureRequestEvent),
    ('xcirculate', XCirculateEvent),
    ('xcirculaterequest', XCirculateRequestEvent),
    ('xproperty', XPropertyEvent),
    ('xselectionclear', XSelectionClearEvent),
    ('xselectionrequest', XSelectionRequestEvent),
    ('xselection', XSelectionEvent),
    ('xcolormap', XColormapEvent),
    ('xclient', XClientMessageEvent),
    ('xmapping', XMappingEvent),
    ('xerror', XErrorEvent),
    ('xkeymap', XKeymapEvent),
    ('pad', c_long * 24),
]

XEvent = struct__XEvent 	# /usr/include/X11/Xlib.h:1012
class struct_anon_62(Structure):
    __slots__ = [
        'lbearing',
        'rbearing',
        'width',
        'ascent',
        'descent',
        'attributes',
    ]
struct_anon_62._fields_ = [
    ('lbearing', c_short),
    ('rbearing', c_short),
    ('width', c_short),
    ('ascent', c_short),
    ('descent', c_short),
    ('attributes', c_ushort),
]

XCharStruct = struct_anon_62 	# /usr/include/X11/Xlib.h:1027
class struct_anon_63(Structure):
    __slots__ = [
        'name',
        'card32',
    ]
struct_anon_63._fields_ = [
    ('name', Atom),
    ('card32', c_ulong),
]

XFontProp = struct_anon_63 	# /usr/include/X11/Xlib.h:1036
class struct_anon_64(Structure):
    __slots__ = [
        'ext_data',
        'fid',
        'direction',
        'min_char_or_byte2',
        'max_char_or_byte2',
        'min_byte1',
        'max_byte1',
        'all_chars_exist',
        'default_char',
        'n_properties',
        'properties',
        'min_bounds',
        'max_bounds',
        'per_char',
        'ascent',
        'descent',
    ]
struct_anon_64._fields_ = [
    ('ext_data', POINTER(XExtData)),
    ('fid', Font),
    ('direction', c_uint),
    ('min_char_or_byte2', c_uint),
    ('max_char_or_byte2', c_uint),
    ('min_byte1', c_uint),
    ('max_byte1', c_uint),
    ('all_chars_exist', c_int),
    ('default_char', c_uint),
    ('n_properties', c_int),
    ('properties', POINTER(XFontProp)),
    ('min_bounds', XCharStruct),
    ('max_bounds', XCharStruct),
    ('per_char', POINTER(XCharStruct)),
    ('ascent', c_int),
    ('descent', c_int),
]

XFontStruct = struct_anon_64 	# /usr/include/X11/Xlib.h:1055
class struct_anon_65(Structure):
    __slots__ = [
        'chars',
        'nchars',
        'delta',
        'font',
    ]
struct_anon_65._fields_ = [
    ('chars', c_char_p),
    ('nchars', c_int),
    ('delta', c_int),
    ('font', Font),
]

XTextItem = struct_anon_65 	# /usr/include/X11/Xlib.h:1065
class struct_anon_66(Structure):
    __slots__ = [
        'byte1',
        'byte2',
    ]
struct_anon_66._fields_ = [
    ('byte1', c_ubyte),
    ('byte2', c_ubyte),
]

XChar2b = struct_anon_66 	# /usr/include/X11/Xlib.h:1070
class struct_anon_67(Structure):
    __slots__ = [
        'chars',
        'nchars',
        'delta',
        'font',
    ]
struct_anon_67._fields_ = [
    ('chars', POINTER(XChar2b)),
    ('nchars', c_int),
    ('delta', c_int),
    ('font', Font),
]

XTextItem16 = struct_anon_67 	# /usr/include/X11/Xlib.h:1077
class struct_anon_68(Union):
    __slots__ = [
        'display',
        'gc',
        'visual',
        'screen',
        'pixmap_format',
        'font',
    ]
struct_anon_68._fields_ = [
    ('display', POINTER(Display)),
    ('gc', GC),
    ('visual', POINTER(Visual)),
    ('screen', POINTER(Screen)),
    ('pixmap_format', POINTER(ScreenFormat)),
    ('font', POINTER(XFontStruct)),
]

XEDataObject = struct_anon_68 	# /usr/include/X11/Xlib.h:1085
class struct_anon_69(Structure):
    __slots__ = [
        'max_ink_extent',
        'max_logical_extent',
    ]
struct_anon_69._fields_ = [
    ('max_ink_extent', XRectangle),
    ('max_logical_extent', XRectangle),
]

XFontSetExtents = struct_anon_69 	# /usr/include/X11/Xlib.h:1090
class struct__XOM(Structure):
    __slots__ = [
    ]
struct__XOM._fields_ = [
    ('_opaque_struct', c_int)
]

XOM = POINTER(struct__XOM) 	# /usr/include/X11/Xlib.h:1096
class struct__XOC(Structure):
    __slots__ = [
    ]
struct__XOC._fields_ = [
    ('_opaque_struct', c_int)
]

XOC = POINTER(struct__XOC) 	# /usr/include/X11/Xlib.h:1097
XFontSet = POINTER(struct__XOC) 	# /usr/include/X11/Xlib.h:1097
class struct_anon_70(Structure):
    __slots__ = [
        'chars',
        'nchars',
        'delta',
        'font_set',
    ]
struct_anon_70._fields_ = [
    ('chars', c_char_p),
    ('nchars', c_int),
    ('delta', c_int),
    ('font_set', XFontSet),
]

XmbTextItem = struct_anon_70 	# /usr/include/X11/Xlib.h:1104
class struct_anon_71(Structure):
    __slots__ = [
        'chars',
        'nchars',
        'delta',
        'font_set',
    ]
struct_anon_71._fields_ = [
    ('chars', c_wchar_p),
    ('nchars', c_int),
    ('delta', c_int),
    ('font_set', XFontSet),
]

XwcTextItem = struct_anon_71 	# /usr/include/X11/Xlib.h:1111
class struct_anon_72(Structure):
    __slots__ = [
        'charset_count',
        'charset_list',
    ]
struct_anon_72._fields_ = [
    ('charset_count', c_int),
    ('charset_list', POINTER(c_char_p)),
]

XOMCharSetList = struct_anon_72 	# /usr/include/X11/Xlib.h:1127
enum_anon_73 = c_int
XOMOrientation_LTR_TTB = 1
XOMOrientation_RTL_TTB = 2
XOMOrientation_TTB_LTR = 3
XOMOrientation_TTB_RTL = 4
XOMOrientation_Context = 5
XOrientation = enum_anon_73 	# /usr/include/X11/Xlib.h:1135
class struct_anon_74(Structure):
    __slots__ = [
        'num_orientation',
        'orientation',
    ]
struct_anon_74._fields_ = [
    ('num_orientation', c_int),
    ('orientation', POINTER(XOrientation)),
]

XOMOrientation = struct_anon_74 	# /usr/include/X11/Xlib.h:1140
class struct_anon_75(Structure):
    __slots__ = [
        'num_font',
        'font_struct_list',
        'font_name_list',
    ]
struct_anon_75._fields_ = [
    ('num_font', c_int),
    ('font_struct_list', POINTER(POINTER(XFontStruct))),
    ('font_name_list', POINTER(c_char_p)),
]

XOMFontInfo = struct_anon_75 	# /usr/include/X11/Xlib.h:1146
class struct__XIM(Structure):
    __slots__ = [
    ]
struct__XIM._fields_ = [
    ('_opaque_struct', c_int)
]

XIM = POINTER(struct__XIM) 	# /usr/include/X11/Xlib.h:1148
class struct__XIC(Structure):
    __slots__ = [
    ]
struct__XIC._fields_ = [
    ('_opaque_struct', c_int)
]

XIC = POINTER(struct__XIC) 	# /usr/include/X11/Xlib.h:1149
XIMProc = CFUNCTYPE(None, XIM, XPointer, XPointer) 	# /usr/include/X11/Xlib.h:1151
XICProc = CFUNCTYPE(c_int, XIC, XPointer, XPointer) 	# /usr/include/X11/Xlib.h:1157
XIDProc = CFUNCTYPE(None, POINTER(Display), XPointer, XPointer) 	# /usr/include/X11/Xlib.h:1163
XIMStyle = c_ulong 	# /usr/include/X11/Xlib.h:1169
class struct_anon_76(Structure):
    __slots__ = [
        'count_styles',
        'supported_styles',
    ]
struct_anon_76._fields_ = [
    ('count_styles', c_ushort),
    ('supported_styles', POINTER(XIMStyle)),
]

XIMStyles = struct_anon_76 	# /usr/include/X11/Xlib.h:1174
XIMPreeditArea = 1 	# /usr/include/X11/Xlib.h:1176
XIMPreeditCallbacks = 2 	# /usr/include/X11/Xlib.h:1177
XIMPreeditPosition = 4 	# /usr/include/X11/Xlib.h:1178
XIMPreeditNothing = 8 	# /usr/include/X11/Xlib.h:1179
XIMPreeditNone = 16 	# /usr/include/X11/Xlib.h:1180
XIMStatusArea = 256 	# /usr/include/X11/Xlib.h:1181
XIMStatusCallbacks = 512 	# /usr/include/X11/Xlib.h:1182
XIMStatusNothing = 1024 	# /usr/include/X11/Xlib.h:1183
XIMStatusNone = 2048 	# /usr/include/X11/Xlib.h:1184
XBufferOverflow = -1 	# /usr/include/X11/Xlib.h:1230
XLookupNone = 1 	# /usr/include/X11/Xlib.h:1231
XLookupChars = 2 	# /usr/include/X11/Xlib.h:1232
XLookupKeySym = 3 	# /usr/include/X11/Xlib.h:1233
XLookupBoth = 4 	# /usr/include/X11/Xlib.h:1234
XVaNestedList = POINTER(None) 	# /usr/include/X11/Xlib.h:1236
class struct_anon_77(Structure):
    __slots__ = [
        'client_data',
        'callback',
    ]
struct_anon_77._fields_ = [
    ('client_data', XPointer),
    ('callback', XIMProc),
]

XIMCallback = struct_anon_77 	# /usr/include/X11/Xlib.h:1241
class struct_anon_78(Structure):
    __slots__ = [
        'client_data',
        'callback',
    ]
struct_anon_78._fields_ = [
    ('client_data', XPointer),
    ('callback', XICProc),
]

XICCallback = struct_anon_78 	# /usr/include/X11/Xlib.h:1246
XIMFeedback = c_ulong 	# /usr/include/X11/Xlib.h:1248
XIMReverse = 1 	# /usr/include/X11/Xlib.h:1250
XIMUnderline = 2 	# /usr/include/X11/Xlib.h:1251
XIMHighlight = 4 	# /usr/include/X11/Xlib.h:1252
XIMPrimary = 32 	# /usr/include/X11/Xlib.h:1253
XIMSecondary = 64 	# /usr/include/X11/Xlib.h:1254
XIMTertiary = 128 	# /usr/include/X11/Xlib.h:1255
XIMVisibleToForward = 256 	# /usr/include/X11/Xlib.h:1256
XIMVisibleToBackword = 512 	# /usr/include/X11/Xlib.h:1257
XIMVisibleToCenter = 1024 	# /usr/include/X11/Xlib.h:1258
class struct__XIMText(Structure):
    __slots__ = [
        'length',
        'feedback',
        'encoding_is_wchar',
        'string',
    ]
class struct_anon_79(Union):
    __slots__ = [
        'multi_byte',
        'wide_char',
    ]
struct_anon_79._fields_ = [
    ('multi_byte', c_char_p),
    ('wide_char', c_wchar_p),
]

struct__XIMText._fields_ = [
    ('length', c_ushort),
    ('feedback', POINTER(XIMFeedback)),
    ('encoding_is_wchar', c_int),
    ('string', struct_anon_79),
]

XIMText = struct__XIMText 	# /usr/include/X11/Xlib.h:1268
XIMPreeditState = c_ulong 	# /usr/include/X11/Xlib.h:1270
XIMPreeditUnKnown = 0 	# /usr/include/X11/Xlib.h:1272
XIMPreeditEnable = 1 	# /usr/include/X11/Xlib.h:1273
XIMPreeditDisable = 2 	# /usr/include/X11/Xlib.h:1274
class struct__XIMPreeditStateNotifyCallbackStruct(Structure):
    __slots__ = [
        'state',
    ]
struct__XIMPreeditStateNotifyCallbackStruct._fields_ = [
    ('state', XIMPreeditState),
]

XIMPreeditStateNotifyCallbackStruct = struct__XIMPreeditStateNotifyCallbackStruct 	# /usr/include/X11/Xlib.h:1278
XIMResetState = c_ulong 	# /usr/include/X11/Xlib.h:1280
XIMInitialState = 1 	# /usr/include/X11/Xlib.h:1282
XIMPreserveState = 2 	# /usr/include/X11/Xlib.h:1283
XIMStringConversionFeedback = c_ulong 	# /usr/include/X11/Xlib.h:1285
XIMStringConversionLeftEdge = 1 	# /usr/include/X11/Xlib.h:1287
XIMStringConversionRightEdge = 2 	# /usr/include/X11/Xlib.h:1288
XIMStringConversionTopEdge = 4 	# /usr/include/X11/Xlib.h:1289
XIMStringConversionBottomEdge = 8 	# /usr/include/X11/Xlib.h:1290
XIMStringConversionConcealed = 16 	# /usr/include/X11/Xlib.h:1291
XIMStringConversionWrapped = 32 	# /usr/include/X11/Xlib.h:1292
class struct__XIMStringConversionText(Structure):
    __slots__ = [
        'length',
        'feedback',
        'encoding_is_wchar',
        'string',
    ]
class struct_anon_80(Union):
    __slots__ = [
        'mbs',
        'wcs',
    ]
struct_anon_80._fields_ = [
    ('mbs', c_char_p),
    ('wcs', c_wchar_p),
]

struct__XIMStringConversionText._fields_ = [
    ('length', c_ushort),
    ('feedback', POINTER(XIMStringConversionFeedback)),
    ('encoding_is_wchar', c_int),
    ('string', struct_anon_80),
]

XIMStringConversionText = struct__XIMStringConversionText 	# /usr/include/X11/Xlib.h:1302
XIMStringConversionPosition = c_ushort 	# /usr/include/X11/Xlib.h:1304
XIMStringConversionType = c_ushort 	# /usr/include/X11/Xlib.h:1306
XIMStringConversionBuffer = 1 	# /usr/include/X11/Xlib.h:1308
XIMStringConversionLine = 2 	# /usr/include/X11/Xlib.h:1309
XIMStringConversionWord = 3 	# /usr/include/X11/Xlib.h:1310
XIMStringConversionChar = 4 	# /usr/include/X11/Xlib.h:1311
XIMStringConversionOperation = c_ushort 	# /usr/include/X11/Xlib.h:1313
XIMStringConversionSubstitution = 1 	# /usr/include/X11/Xlib.h:1315
XIMStringConversionRetrieval = 2 	# /usr/include/X11/Xlib.h:1316
enum_anon_81 = c_int
XIMForwardChar = 1
XIMBackwardChar = 2
XIMForwardWord = 3
XIMBackwardWord = 4
XIMCaretUp = 5
XIMCaretDown = 6
XIMNextLine = 7
XIMPreviousLine = 8
XIMLineStart = 9
XIMLineEnd = 10
XIMAbsolutePosition = 11
XIMDontChange = 12
XIMCaretDirection = enum_anon_81 	# /usr/include/X11/Xlib.h:1326
class struct__XIMStringConversionCallbackStruct(Structure):
    __slots__ = [
        'position',
        'direction',
        'operation',
        'factor',
        'text',
    ]
struct__XIMStringConversionCallbackStruct._fields_ = [
    ('position', XIMStringConversionPosition),
    ('direction', XIMCaretDirection),
    ('operation', XIMStringConversionOperation),
    ('factor', c_ushort),
    ('text', POINTER(XIMStringConversionText)),
]

XIMStringConversionCallbackStruct = struct__XIMStringConversionCallbackStruct 	# /usr/include/X11/Xlib.h:1334
class struct__XIMPreeditDrawCallbackStruct(Structure):
    __slots__ = [
        'caret',
        'chg_first',
        'chg_length',
        'text',
    ]
struct__XIMPreeditDrawCallbackStruct._fields_ = [
    ('caret', c_int),
    ('chg_first', c_int),
    ('chg_length', c_int),
    ('text', POINTER(XIMText)),
]

XIMPreeditDrawCallbackStruct = struct__XIMPreeditDrawCallbackStruct 	# /usr/include/X11/Xlib.h:1341
enum_anon_82 = c_int
XIMIsInvisible = 1
XIMIsPrimary = 2
XIMIsSecondary = 3
XIMCaretStyle = enum_anon_82 	# /usr/include/X11/Xlib.h:1347
class struct__XIMPreeditCaretCallbackStruct(Structure):
    __slots__ = [
        'position',
        'direction',
        'style',
    ]
struct__XIMPreeditCaretCallbackStruct._fields_ = [
    ('position', c_int),
    ('direction', XIMCaretDirection),
    ('style', XIMCaretStyle),
]

XIMPreeditCaretCallbackStruct = struct__XIMPreeditCaretCallbackStruct 	# /usr/include/X11/Xlib.h:1353
enum_anon_83 = c_int
XIMTextType = 1
XIMBitmapType = 2
XIMStatusDataType = enum_anon_83 	# /usr/include/X11/Xlib.h:1358
class struct__XIMStatusDrawCallbackStruct(Structure):
    __slots__ = [
        'type',
        'data',
    ]
class struct_anon_84(Union):
    __slots__ = [
        'text',
        'bitmap',
    ]
struct_anon_84._fields_ = [
    ('text', POINTER(XIMText)),
    ('bitmap', Pixmap),
]

struct__XIMStatusDrawCallbackStruct._fields_ = [
    ('type', XIMStatusDataType),
    ('data', struct_anon_84),
]

XIMStatusDrawCallbackStruct = struct__XIMStatusDrawCallbackStruct 	# /usr/include/X11/Xlib.h:1366
class struct__XIMHotKeyTrigger(Structure):
    __slots__ = [
        'keysym',
        'modifier',
        'modifier_mask',
    ]
struct__XIMHotKeyTrigger._fields_ = [
    ('keysym', KeySym),
    ('modifier', c_int),
    ('modifier_mask', c_int),
]

XIMHotKeyTrigger = struct__XIMHotKeyTrigger 	# /usr/include/X11/Xlib.h:1372
class struct__XIMHotKeyTriggers(Structure):
    __slots__ = [
        'num_hot_key',
        'key',
    ]
struct__XIMHotKeyTriggers._fields_ = [
    ('num_hot_key', c_int),
    ('key', POINTER(XIMHotKeyTrigger)),
]

XIMHotKeyTriggers = struct__XIMHotKeyTriggers 	# /usr/include/X11/Xlib.h:1377
XIMHotKeyState = c_ulong 	# /usr/include/X11/Xlib.h:1379
XIMHotKeyStateON = 1 	# /usr/include/X11/Xlib.h:1381
XIMHotKeyStateOFF = 2 	# /usr/include/X11/Xlib.h:1382
class struct_anon_85(Structure):
    __slots__ = [
        'count_values',
        'supported_values',
    ]
struct_anon_85._fields_ = [
    ('count_values', c_ushort),
    ('supported_values', POINTER(c_char_p)),
]

XIMValuesList = struct_anon_85 	# /usr/include/X11/Xlib.h:1387
# /usr/include/X11/Xlib.h:1397
XLoadQueryFont = _lib.XLoadQueryFont
XLoadQueryFont.restype = POINTER(XFontStruct)
XLoadQueryFont.argtypes = [POINTER(Display), c_char_p]

# /usr/include/X11/Xlib.h:1402
XQueryFont = _lib.XQueryFont
XQueryFont.restype = POINTER(XFontStruct)
XQueryFont.argtypes = [POINTER(Display), XID]

# /usr/include/X11/Xlib.h:1408
XGetMotionEvents = _lib.XGetMotionEvents
XGetMotionEvents.restype = POINTER(XTimeCoord)
XGetMotionEvents.argtypes = [POINTER(Display), Window, Time, Time, POINTER(c_int)]

# /usr/include/X11/Xlib.h:1416
XDeleteModifiermapEntry = _lib.XDeleteModifiermapEntry
XDeleteModifiermapEntry.restype = POINTER(XModifierKeymap)
XDeleteModifiermapEntry.argtypes = [POINTER(XModifierKeymap), KeyCode, c_int]

# /usr/include/X11/Xlib.h:1426
XGetModifierMapping = _lib.XGetModifierMapping
XGetModifierMapping.restype = POINTER(XModifierKeymap)
XGetModifierMapping.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1430
XInsertModifiermapEntry = _lib.XInsertModifiermapEntry
XInsertModifiermapEntry.restype = POINTER(XModifierKeymap)
XInsertModifiermapEntry.argtypes = [POINTER(XModifierKeymap), KeyCode, c_int]

# /usr/include/X11/Xlib.h:1440
XNewModifiermap = _lib.XNewModifiermap
XNewModifiermap.restype = POINTER(XModifierKeymap)
XNewModifiermap.argtypes = [c_int]

# /usr/include/X11/Xlib.h:1444
XCreateImage = _lib.XCreateImage
XCreateImage.restype = POINTER(XImage)
XCreateImage.argtypes = [POINTER(Display), POINTER(Visual), c_uint, c_int, c_int, c_char_p, c_uint, c_uint, c_int, c_int]

# /usr/include/X11/Xlib.h:1456
XInitImage = _lib.XInitImage
XInitImage.restype = c_int
XInitImage.argtypes = [POINTER(XImage)]

# /usr/include/X11/Xlib.h:1459
XGetImage = _lib.XGetImage
XGetImage.restype = POINTER(XImage)
XGetImage.argtypes = [POINTER(Display), Drawable, c_int, c_int, c_uint, c_uint, c_ulong, c_int]

# /usr/include/X11/Xlib.h:1469
XGetSubImage = _lib.XGetSubImage
XGetSubImage.restype = POINTER(XImage)
XGetSubImage.argtypes = [POINTER(Display), Drawable, c_int, c_int, c_uint, c_uint, c_ulong, c_int, POINTER(XImage), c_int, c_int]

# /usr/include/X11/Xlib.h:1486
XOpenDisplay = _lib.XOpenDisplay
XOpenDisplay.restype = POINTER(Display)
XOpenDisplay.argtypes = [c_char_p]

# /usr/include/X11/Xlib.h:1490
XrmInitialize = _lib.XrmInitialize
XrmInitialize.restype = None
XrmInitialize.argtypes = []

# /usr/include/X11/Xlib.h:1494
XFetchBytes = _lib.XFetchBytes
XFetchBytes.restype = c_char_p
XFetchBytes.argtypes = [POINTER(Display), POINTER(c_int)]

# /usr/include/X11/Xlib.h:1498
XFetchBuffer = _lib.XFetchBuffer
XFetchBuffer.restype = c_char_p
XFetchBuffer.argtypes = [POINTER(Display), POINTER(c_int), c_int]

# /usr/include/X11/Xlib.h:1503
XGetAtomName = _lib.XGetAtomName
XGetAtomName.restype = c_char_p
XGetAtomName.argtypes = [POINTER(Display), Atom]

# /usr/include/X11/Xlib.h:1507
XGetAtomNames = _lib.XGetAtomNames
XGetAtomNames.restype = c_int
XGetAtomNames.argtypes = [POINTER(Display), POINTER(Atom), c_int, POINTER(c_char_p)]

# /usr/include/X11/Xlib.h:1513
XGetDefault = _lib.XGetDefault
XGetDefault.restype = c_char_p
XGetDefault.argtypes = [POINTER(Display), c_char_p, c_char_p]

# /usr/include/X11/Xlib.h:1518
XDisplayName = _lib.XDisplayName
XDisplayName.restype = c_char_p
XDisplayName.argtypes = [c_char_p]

# /usr/include/X11/Xlib.h:1521
XKeysymToString = _lib.XKeysymToString
XKeysymToString.restype = c_char_p
XKeysymToString.argtypes = [KeySym]

# /usr/include/X11/Xlib.h:1525
XSynchronize = _lib.XSynchronize
XSynchronize.restype = POINTER(CFUNCTYPE(c_int, POINTER(Display)))
XSynchronize.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:1531
XSetAfterFunction = _lib.XSetAfterFunction
XSetAfterFunction.restype = POINTER(CFUNCTYPE(c_int, POINTER(Display)))
XSetAfterFunction.argtypes = [POINTER(Display), CFUNCTYPE(c_int, POINTER(Display))]

# /usr/include/X11/Xlib.h:1539
XInternAtom = _lib.XInternAtom
XInternAtom.restype = Atom
XInternAtom.argtypes = [POINTER(Display), c_char_p, c_int]

# /usr/include/X11/Xlib.h:1544
XInternAtoms = _lib.XInternAtoms
XInternAtoms.restype = c_int
XInternAtoms.argtypes = [POINTER(Display), POINTER(c_char_p), c_int, c_int, POINTER(Atom)]

# /usr/include/X11/Xlib.h:1551
XCopyColormapAndFree = _lib.XCopyColormapAndFree
XCopyColormapAndFree.restype = Colormap
XCopyColormapAndFree.argtypes = [POINTER(Display), Colormap]

# /usr/include/X11/Xlib.h:1555
XCreateColormap = _lib.XCreateColormap
XCreateColormap.restype = Colormap
XCreateColormap.argtypes = [POINTER(Display), Window, POINTER(Visual), c_int]

# /usr/include/X11/Xlib.h:1561
XCreatePixmapCursor = _lib.XCreatePixmapCursor
XCreatePixmapCursor.restype = Cursor
XCreatePixmapCursor.argtypes = [POINTER(Display), Pixmap, Pixmap, POINTER(XColor), POINTER(XColor), c_uint, c_uint]

# /usr/include/X11/Xlib.h:1570
XCreateGlyphCursor = _lib.XCreateGlyphCursor
XCreateGlyphCursor.restype = Cursor
XCreateGlyphCursor.argtypes = [POINTER(Display), Font, Font, c_uint, c_uint, POINTER(XColor), POINTER(XColor)]

# /usr/include/X11/Xlib.h:1579
XCreateFontCursor = _lib.XCreateFontCursor
XCreateFontCursor.restype = Cursor
XCreateFontCursor.argtypes = [POINTER(Display), c_uint]

# /usr/include/X11/Xlib.h:1583
XLoadFont = _lib.XLoadFont
XLoadFont.restype = Font
XLoadFont.argtypes = [POINTER(Display), c_char_p]

# /usr/include/X11/Xlib.h:1587
XCreateGC = _lib.XCreateGC
XCreateGC.restype = GC
XCreateGC.argtypes = [POINTER(Display), Drawable, c_ulong, POINTER(XGCValues)]

# /usr/include/X11/Xlib.h:1593
XGContextFromGC = _lib.XGContextFromGC
XGContextFromGC.restype = GContext
XGContextFromGC.argtypes = [GC]

# /usr/include/X11/Xlib.h:1596
XFlushGC = _lib.XFlushGC
XFlushGC.restype = None
XFlushGC.argtypes = [POINTER(Display), GC]

# /usr/include/X11/Xlib.h:1600
XCreatePixmap = _lib.XCreatePixmap
XCreatePixmap.restype = Pixmap
XCreatePixmap.argtypes = [POINTER(Display), Drawable, c_uint, c_uint, c_uint]

# /usr/include/X11/Xlib.h:1607
XCreateBitmapFromData = _lib.XCreateBitmapFromData
XCreateBitmapFromData.restype = Pixmap
XCreateBitmapFromData.argtypes = [POINTER(Display), Drawable, c_char_p, c_uint, c_uint]

# /usr/include/X11/Xlib.h:1614
XCreatePixmapFromBitmapData = _lib.XCreatePixmapFromBitmapData
XCreatePixmapFromBitmapData.restype = Pixmap
XCreatePixmapFromBitmapData.argtypes = [POINTER(Display), Drawable, c_char_p, c_uint, c_uint, c_ulong, c_ulong, c_uint]

# /usr/include/X11/Xlib.h:1624
XCreateSimpleWindow = _lib.XCreateSimpleWindow
XCreateSimpleWindow.restype = Window
XCreateSimpleWindow.argtypes = [POINTER(Display), Window, c_int, c_int, c_uint, c_uint, c_uint, c_ulong, c_ulong]

# /usr/include/X11/Xlib.h:1635
XGetSelectionOwner = _lib.XGetSelectionOwner
XGetSelectionOwner.restype = Window
XGetSelectionOwner.argtypes = [POINTER(Display), Atom]

# /usr/include/X11/Xlib.h:1639
XCreateWindow = _lib.XCreateWindow
XCreateWindow.restype = Window
XCreateWindow.argtypes = [POINTER(Display), Window, c_int, c_int, c_uint, c_uint, c_uint, c_int, c_uint, POINTER(Visual), c_ulong, POINTER(XSetWindowAttributes)]

# /usr/include/X11/Xlib.h:1653
XListInstalledColormaps = _lib.XListInstalledColormaps
XListInstalledColormaps.restype = POINTER(Colormap)
XListInstalledColormaps.argtypes = [POINTER(Display), Window, POINTER(c_int)]

# /usr/include/X11/Xlib.h:1658
XListFonts = _lib.XListFonts
XListFonts.restype = POINTER(c_char_p)
XListFonts.argtypes = [POINTER(Display), c_char_p, c_int, POINTER(c_int)]

# /usr/include/X11/Xlib.h:1664
XListFontsWithInfo = _lib.XListFontsWithInfo
XListFontsWithInfo.restype = POINTER(c_char_p)
XListFontsWithInfo.argtypes = [POINTER(Display), c_char_p, c_int, POINTER(c_int), POINTER(POINTER(XFontStruct))]

# /usr/include/X11/Xlib.h:1671
XGetFontPath = _lib.XGetFontPath
XGetFontPath.restype = POINTER(c_char_p)
XGetFontPath.argtypes = [POINTER(Display), POINTER(c_int)]

# /usr/include/X11/Xlib.h:1675
XListExtensions = _lib.XListExtensions
XListExtensions.restype = POINTER(c_char_p)
XListExtensions.argtypes = [POINTER(Display), POINTER(c_int)]

# /usr/include/X11/Xlib.h:1679
XListProperties = _lib.XListProperties
XListProperties.restype = POINTER(Atom)
XListProperties.argtypes = [POINTER(Display), Window, POINTER(c_int)]

# /usr/include/X11/Xlib.h:1684
XListHosts = _lib.XListHosts
XListHosts.restype = POINTER(XHostAddress)
XListHosts.argtypes = [POINTER(Display), POINTER(c_int), POINTER(c_int)]

# /usr/include/X11/Xlib.h:1689
XKeycodeToKeysym = _lib.XKeycodeToKeysym
XKeycodeToKeysym.restype = KeySym
XKeycodeToKeysym.argtypes = [POINTER(Display), KeyCode, c_int]

# /usr/include/X11/Xlib.h:1698
XLookupKeysym = _lib.XLookupKeysym
XLookupKeysym.restype = KeySym
XLookupKeysym.argtypes = [POINTER(XKeyEvent), c_int]

# /usr/include/X11/Xlib.h:1702
XGetKeyboardMapping = _lib.XGetKeyboardMapping
XGetKeyboardMapping.restype = POINTER(KeySym)
XGetKeyboardMapping.argtypes = [POINTER(Display), KeyCode, c_int, POINTER(c_int)]

# /usr/include/X11/Xlib.h:1712
XStringToKeysym = _lib.XStringToKeysym
XStringToKeysym.restype = KeySym
XStringToKeysym.argtypes = [c_char_p]

# /usr/include/X11/Xlib.h:1715
XMaxRequestSize = _lib.XMaxRequestSize
XMaxRequestSize.restype = c_long
XMaxRequestSize.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1718
XExtendedMaxRequestSize = _lib.XExtendedMaxRequestSize
XExtendedMaxRequestSize.restype = c_long
XExtendedMaxRequestSize.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1721
XResourceManagerString = _lib.XResourceManagerString
XResourceManagerString.restype = c_char_p
XResourceManagerString.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1724
XScreenResourceString = _lib.XScreenResourceString
XScreenResourceString.restype = c_char_p
XScreenResourceString.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:1727
XDisplayMotionBufferSize = _lib.XDisplayMotionBufferSize
XDisplayMotionBufferSize.restype = c_ulong
XDisplayMotionBufferSize.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1730
XVisualIDFromVisual = _lib.XVisualIDFromVisual
XVisualIDFromVisual.restype = VisualID
XVisualIDFromVisual.argtypes = [POINTER(Visual)]

# /usr/include/X11/Xlib.h:1736
XInitThreads = _lib.XInitThreads
XInitThreads.restype = c_int
XInitThreads.argtypes = []

# /usr/include/X11/Xlib.h:1740
XLockDisplay = _lib.XLockDisplay
XLockDisplay.restype = None
XLockDisplay.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1744
XUnlockDisplay = _lib.XUnlockDisplay
XUnlockDisplay.restype = None
XUnlockDisplay.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1750
XInitExtension = _lib.XInitExtension
XInitExtension.restype = POINTER(XExtCodes)
XInitExtension.argtypes = [POINTER(Display), c_char_p]

# /usr/include/X11/Xlib.h:1755
XAddExtension = _lib.XAddExtension
XAddExtension.restype = POINTER(XExtCodes)
XAddExtension.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1758
XFindOnExtensionList = _lib.XFindOnExtensionList
XFindOnExtensionList.restype = POINTER(XExtData)
XFindOnExtensionList.argtypes = [POINTER(POINTER(XExtData)), c_int]

# /usr/include/X11/Xlib.h:1762
XEHeadOfExtensionList = _lib.XEHeadOfExtensionList
XEHeadOfExtensionList.restype = POINTER(POINTER(XExtData))
XEHeadOfExtensionList.argtypes = [XEDataObject]

# /usr/include/X11/Xlib.h:1767
XRootWindow = _lib.XRootWindow
XRootWindow.restype = Window
XRootWindow.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:1771
XDefaultRootWindow = _lib.XDefaultRootWindow
XDefaultRootWindow.restype = Window
XDefaultRootWindow.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1774
XRootWindowOfScreen = _lib.XRootWindowOfScreen
XRootWindowOfScreen.restype = Window
XRootWindowOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:1777
XDefaultVisual = _lib.XDefaultVisual
XDefaultVisual.restype = POINTER(Visual)
XDefaultVisual.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:1781
XDefaultVisualOfScreen = _lib.XDefaultVisualOfScreen
XDefaultVisualOfScreen.restype = POINTER(Visual)
XDefaultVisualOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:1784
XDefaultGC = _lib.XDefaultGC
XDefaultGC.restype = GC
XDefaultGC.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:1788
XDefaultGCOfScreen = _lib.XDefaultGCOfScreen
XDefaultGCOfScreen.restype = GC
XDefaultGCOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:1791
XBlackPixel = _lib.XBlackPixel
XBlackPixel.restype = c_ulong
XBlackPixel.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:1795
XWhitePixel = _lib.XWhitePixel
XWhitePixel.restype = c_ulong
XWhitePixel.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:1799
XAllPlanes = _lib.XAllPlanes
XAllPlanes.restype = c_ulong
XAllPlanes.argtypes = []

# /usr/include/X11/Xlib.h:1802
XBlackPixelOfScreen = _lib.XBlackPixelOfScreen
XBlackPixelOfScreen.restype = c_ulong
XBlackPixelOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:1805
XWhitePixelOfScreen = _lib.XWhitePixelOfScreen
XWhitePixelOfScreen.restype = c_ulong
XWhitePixelOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:1808
XNextRequest = _lib.XNextRequest
XNextRequest.restype = c_ulong
XNextRequest.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1811
XLastKnownRequestProcessed = _lib.XLastKnownRequestProcessed
XLastKnownRequestProcessed.restype = c_ulong
XLastKnownRequestProcessed.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1814
XServerVendor = _lib.XServerVendor
XServerVendor.restype = c_char_p
XServerVendor.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1817
XDisplayString = _lib.XDisplayString
XDisplayString.restype = c_char_p
XDisplayString.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1820
XDefaultColormap = _lib.XDefaultColormap
XDefaultColormap.restype = Colormap
XDefaultColormap.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:1824
XDefaultColormapOfScreen = _lib.XDefaultColormapOfScreen
XDefaultColormapOfScreen.restype = Colormap
XDefaultColormapOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:1827
XDisplayOfScreen = _lib.XDisplayOfScreen
XDisplayOfScreen.restype = POINTER(Display)
XDisplayOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:1830
XScreenOfDisplay = _lib.XScreenOfDisplay
XScreenOfDisplay.restype = POINTER(Screen)
XScreenOfDisplay.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:1834
XDefaultScreenOfDisplay = _lib.XDefaultScreenOfDisplay
XDefaultScreenOfDisplay.restype = POINTER(Screen)
XDefaultScreenOfDisplay.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1837
XEventMaskOfScreen = _lib.XEventMaskOfScreen
XEventMaskOfScreen.restype = c_long
XEventMaskOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:1841
XScreenNumberOfScreen = _lib.XScreenNumberOfScreen
XScreenNumberOfScreen.restype = c_int
XScreenNumberOfScreen.argtypes = [POINTER(Screen)]

XErrorHandler = CFUNCTYPE(c_int, POINTER(Display), POINTER(XErrorEvent)) 	# /usr/include/X11/Xlib.h:1845
# /usr/include/X11/Xlib.h:1850
XSetErrorHandler = _lib.XSetErrorHandler
XSetErrorHandler.restype = XErrorHandler
XSetErrorHandler.argtypes = [XErrorHandler]

XIOErrorHandler = CFUNCTYPE(c_int, POINTER(Display)) 	# /usr/include/X11/Xlib.h:1855
# /usr/include/X11/Xlib.h:1859
XSetIOErrorHandler = _lib.XSetIOErrorHandler
XSetIOErrorHandler.restype = XIOErrorHandler
XSetIOErrorHandler.argtypes = [XIOErrorHandler]

# /usr/include/X11/Xlib.h:1864
XListPixmapFormats = _lib.XListPixmapFormats
XListPixmapFormats.restype = POINTER(XPixmapFormatValues)
XListPixmapFormats.argtypes = [POINTER(Display), POINTER(c_int)]

# /usr/include/X11/Xlib.h:1868
XListDepths = _lib.XListDepths
XListDepths.restype = POINTER(c_int)
XListDepths.argtypes = [POINTER(Display), c_int, POINTER(c_int)]

# /usr/include/X11/Xlib.h:1876
XReconfigureWMWindow = _lib.XReconfigureWMWindow
XReconfigureWMWindow.restype = c_int
XReconfigureWMWindow.argtypes = [POINTER(Display), Window, c_int, c_uint, POINTER(XWindowChanges)]

# /usr/include/X11/Xlib.h:1884
XGetWMProtocols = _lib.XGetWMProtocols
XGetWMProtocols.restype = c_int
XGetWMProtocols.argtypes = [POINTER(Display), Window, POINTER(POINTER(Atom)), POINTER(c_int)]

# /usr/include/X11/Xlib.h:1890
XSetWMProtocols = _lib.XSetWMProtocols
XSetWMProtocols.restype = c_int
XSetWMProtocols.argtypes = [POINTER(Display), Window, POINTER(Atom), c_int]

# /usr/include/X11/Xlib.h:1896
XIconifyWindow = _lib.XIconifyWindow
XIconifyWindow.restype = c_int
XIconifyWindow.argtypes = [POINTER(Display), Window, c_int]

# /usr/include/X11/Xlib.h:1901
XWithdrawWindow = _lib.XWithdrawWindow
XWithdrawWindow.restype = c_int
XWithdrawWindow.argtypes = [POINTER(Display), Window, c_int]

# /usr/include/X11/Xlib.h:1906
XGetCommand = _lib.XGetCommand
XGetCommand.restype = c_int
XGetCommand.argtypes = [POINTER(Display), Window, POINTER(POINTER(c_char_p)), POINTER(c_int)]

# /usr/include/X11/Xlib.h:1912
XGetWMColormapWindows = _lib.XGetWMColormapWindows
XGetWMColormapWindows.restype = c_int
XGetWMColormapWindows.argtypes = [POINTER(Display), Window, POINTER(POINTER(Window)), POINTER(c_int)]

# /usr/include/X11/Xlib.h:1918
XSetWMColormapWindows = _lib.XSetWMColormapWindows
XSetWMColormapWindows.restype = c_int
XSetWMColormapWindows.argtypes = [POINTER(Display), Window, POINTER(Window), c_int]

# /usr/include/X11/Xlib.h:1924
XFreeStringList = _lib.XFreeStringList
XFreeStringList.restype = None
XFreeStringList.argtypes = [POINTER(c_char_p)]

# /usr/include/X11/Xlib.h:1927
XSetTransientForHint = _lib.XSetTransientForHint
XSetTransientForHint.restype = c_int
XSetTransientForHint.argtypes = [POINTER(Display), Window, Window]

# /usr/include/X11/Xlib.h:1935
XActivateScreenSaver = _lib.XActivateScreenSaver
XActivateScreenSaver.restype = c_int
XActivateScreenSaver.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1939
XAddHost = _lib.XAddHost
XAddHost.restype = c_int
XAddHost.argtypes = [POINTER(Display), POINTER(XHostAddress)]

# /usr/include/X11/Xlib.h:1944
XAddHosts = _lib.XAddHosts
XAddHosts.restype = c_int
XAddHosts.argtypes = [POINTER(Display), POINTER(XHostAddress), c_int]

class struct__XExtData(Structure):
    __slots__ = [
    ]
struct__XExtData._fields_ = [
    ('_opaque_struct', c_int)
]

# /usr/include/X11/Xlib.h:1950
XAddToExtensionList = _lib.XAddToExtensionList
XAddToExtensionList.restype = c_int
XAddToExtensionList.argtypes = [POINTER(POINTER(struct__XExtData)), POINTER(XExtData)]

# /usr/include/X11/Xlib.h:1955
XAddToSaveSet = _lib.XAddToSaveSet
XAddToSaveSet.restype = c_int
XAddToSaveSet.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:1960
XAllocColor = _lib.XAllocColor
XAllocColor.restype = c_int
XAllocColor.argtypes = [POINTER(Display), Colormap, POINTER(XColor)]

# /usr/include/X11/Xlib.h:1966
XAllocColorCells = _lib.XAllocColorCells
XAllocColorCells.restype = c_int
XAllocColorCells.argtypes = [POINTER(Display), Colormap, c_int, POINTER(c_ulong), c_uint, POINTER(c_ulong), c_uint]

# /usr/include/X11/Xlib.h:1976
XAllocColorPlanes = _lib.XAllocColorPlanes
XAllocColorPlanes.restype = c_int
XAllocColorPlanes.argtypes = [POINTER(Display), Colormap, c_int, POINTER(c_ulong), c_int, c_int, c_int, c_int, POINTER(c_ulong), POINTER(c_ulong), POINTER(c_ulong)]

# /usr/include/X11/Xlib.h:1990
XAllocNamedColor = _lib.XAllocNamedColor
XAllocNamedColor.restype = c_int
XAllocNamedColor.argtypes = [POINTER(Display), Colormap, c_char_p, POINTER(XColor), POINTER(XColor)]

# /usr/include/X11/Xlib.h:1998
XAllowEvents = _lib.XAllowEvents
XAllowEvents.restype = c_int
XAllowEvents.argtypes = [POINTER(Display), c_int, Time]

# /usr/include/X11/Xlib.h:2004
XAutoRepeatOff = _lib.XAutoRepeatOff
XAutoRepeatOff.restype = c_int
XAutoRepeatOff.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2008
XAutoRepeatOn = _lib.XAutoRepeatOn
XAutoRepeatOn.restype = c_int
XAutoRepeatOn.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2012
XBell = _lib.XBell
XBell.restype = c_int
XBell.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:2017
XBitmapBitOrder = _lib.XBitmapBitOrder
XBitmapBitOrder.restype = c_int
XBitmapBitOrder.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2021
XBitmapPad = _lib.XBitmapPad
XBitmapPad.restype = c_int
XBitmapPad.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2025
XBitmapUnit = _lib.XBitmapUnit
XBitmapUnit.restype = c_int
XBitmapUnit.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2029
XCellsOfScreen = _lib.XCellsOfScreen
XCellsOfScreen.restype = c_int
XCellsOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:2033
XChangeActivePointerGrab = _lib.XChangeActivePointerGrab
XChangeActivePointerGrab.restype = c_int
XChangeActivePointerGrab.argtypes = [POINTER(Display), c_uint, Cursor, Time]

# /usr/include/X11/Xlib.h:2040
XChangeGC = _lib.XChangeGC
XChangeGC.restype = c_int
XChangeGC.argtypes = [POINTER(Display), GC, c_ulong, POINTER(XGCValues)]

# /usr/include/X11/Xlib.h:2047
XChangeKeyboardControl = _lib.XChangeKeyboardControl
XChangeKeyboardControl.restype = c_int
XChangeKeyboardControl.argtypes = [POINTER(Display), c_ulong, POINTER(XKeyboardControl)]

# /usr/include/X11/Xlib.h:2053
XChangeKeyboardMapping = _lib.XChangeKeyboardMapping
XChangeKeyboardMapping.restype = c_int
XChangeKeyboardMapping.argtypes = [POINTER(Display), c_int, c_int, POINTER(KeySym), c_int]

# /usr/include/X11/Xlib.h:2061
XChangePointerControl = _lib.XChangePointerControl
XChangePointerControl.restype = c_int
XChangePointerControl.argtypes = [POINTER(Display), c_int, c_int, c_int, c_int, c_int]

# /usr/include/X11/Xlib.h:2070
XChangeProperty = _lib.XChangeProperty
XChangeProperty.restype = c_int
XChangeProperty.argtypes = [POINTER(Display), Window, Atom, Atom, c_int, c_int, POINTER(c_ubyte), c_int]

# /usr/include/X11/Xlib.h:2081
XChangeSaveSet = _lib.XChangeSaveSet
XChangeSaveSet.restype = c_int
XChangeSaveSet.argtypes = [POINTER(Display), Window, c_int]

# /usr/include/X11/Xlib.h:2087
XChangeWindowAttributes = _lib.XChangeWindowAttributes
XChangeWindowAttributes.restype = c_int
XChangeWindowAttributes.argtypes = [POINTER(Display), Window, c_ulong, POINTER(XSetWindowAttributes)]

# /usr/include/X11/Xlib.h:2094
XCheckIfEvent = _lib.XCheckIfEvent
XCheckIfEvent.restype = c_int
XCheckIfEvent.argtypes = [POINTER(Display), POINTER(XEvent), CFUNCTYPE(c_int, POINTER(Display), POINTER(XEvent), XPointer), XPointer]

# /usr/include/X11/Xlib.h:2105
XCheckMaskEvent = _lib.XCheckMaskEvent
XCheckMaskEvent.restype = c_int
XCheckMaskEvent.argtypes = [POINTER(Display), c_long, POINTER(XEvent)]

# /usr/include/X11/Xlib.h:2111
XCheckTypedEvent = _lib.XCheckTypedEvent
XCheckTypedEvent.restype = c_int
XCheckTypedEvent.argtypes = [POINTER(Display), c_int, POINTER(XEvent)]

# /usr/include/X11/Xlib.h:2117
XCheckTypedWindowEvent = _lib.XCheckTypedWindowEvent
XCheckTypedWindowEvent.restype = c_int
XCheckTypedWindowEvent.argtypes = [POINTER(Display), Window, c_int, POINTER(XEvent)]

# /usr/include/X11/Xlib.h:2124
XCheckWindowEvent = _lib.XCheckWindowEvent
XCheckWindowEvent.restype = c_int
XCheckWindowEvent.argtypes = [POINTER(Display), Window, c_long, POINTER(XEvent)]

# /usr/include/X11/Xlib.h:2131
XCirculateSubwindows = _lib.XCirculateSubwindows
XCirculateSubwindows.restype = c_int
XCirculateSubwindows.argtypes = [POINTER(Display), Window, c_int]

# /usr/include/X11/Xlib.h:2137
XCirculateSubwindowsDown = _lib.XCirculateSubwindowsDown
XCirculateSubwindowsDown.restype = c_int
XCirculateSubwindowsDown.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:2142
XCirculateSubwindowsUp = _lib.XCirculateSubwindowsUp
XCirculateSubwindowsUp.restype = c_int
XCirculateSubwindowsUp.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:2147
XClearArea = _lib.XClearArea
XClearArea.restype = c_int
XClearArea.argtypes = [POINTER(Display), Window, c_int, c_int, c_uint, c_uint, c_int]

# /usr/include/X11/Xlib.h:2157
XClearWindow = _lib.XClearWindow
XClearWindow.restype = c_int
XClearWindow.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:2162
XCloseDisplay = _lib.XCloseDisplay
XCloseDisplay.restype = c_int
XCloseDisplay.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2166
XConfigureWindow = _lib.XConfigureWindow
XConfigureWindow.restype = c_int
XConfigureWindow.argtypes = [POINTER(Display), Window, c_uint, POINTER(XWindowChanges)]

# /usr/include/X11/Xlib.h:2173
XConnectionNumber = _lib.XConnectionNumber
XConnectionNumber.restype = c_int
XConnectionNumber.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2177
XConvertSelection = _lib.XConvertSelection
XConvertSelection.restype = c_int
XConvertSelection.argtypes = [POINTER(Display), Atom, Atom, Atom, Window, Time]

# /usr/include/X11/Xlib.h:2186
XCopyArea = _lib.XCopyArea
XCopyArea.restype = c_int
XCopyArea.argtypes = [POINTER(Display), Drawable, Drawable, GC, c_int, c_int, c_uint, c_uint, c_int, c_int]

# /usr/include/X11/Xlib.h:2199
XCopyGC = _lib.XCopyGC
XCopyGC.restype = c_int
XCopyGC.argtypes = [POINTER(Display), GC, c_ulong, GC]

# /usr/include/X11/Xlib.h:2206
XCopyPlane = _lib.XCopyPlane
XCopyPlane.restype = c_int
XCopyPlane.argtypes = [POINTER(Display), Drawable, Drawable, GC, c_int, c_int, c_uint, c_uint, c_int, c_int, c_ulong]

# /usr/include/X11/Xlib.h:2220
XDefaultDepth = _lib.XDefaultDepth
XDefaultDepth.restype = c_int
XDefaultDepth.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:2225
XDefaultDepthOfScreen = _lib.XDefaultDepthOfScreen
XDefaultDepthOfScreen.restype = c_int
XDefaultDepthOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:2229
XDefaultScreen = _lib.XDefaultScreen
XDefaultScreen.restype = c_int
XDefaultScreen.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2233
XDefineCursor = _lib.XDefineCursor
XDefineCursor.restype = c_int
XDefineCursor.argtypes = [POINTER(Display), Window, Cursor]

# /usr/include/X11/Xlib.h:2239
XDeleteProperty = _lib.XDeleteProperty
XDeleteProperty.restype = c_int
XDeleteProperty.argtypes = [POINTER(Display), Window, Atom]

# /usr/include/X11/Xlib.h:2245
XDestroyWindow = _lib.XDestroyWindow
XDestroyWindow.restype = c_int
XDestroyWindow.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:2250
XDestroySubwindows = _lib.XDestroySubwindows
XDestroySubwindows.restype = c_int
XDestroySubwindows.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:2255
XDoesBackingStore = _lib.XDoesBackingStore
XDoesBackingStore.restype = c_int
XDoesBackingStore.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:2259
XDoesSaveUnders = _lib.XDoesSaveUnders
XDoesSaveUnders.restype = c_int
XDoesSaveUnders.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:2263
XDisableAccessControl = _lib.XDisableAccessControl
XDisableAccessControl.restype = c_int
XDisableAccessControl.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2268
XDisplayCells = _lib.XDisplayCells
XDisplayCells.restype = c_int
XDisplayCells.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:2273
XDisplayHeight = _lib.XDisplayHeight
XDisplayHeight.restype = c_int
XDisplayHeight.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:2278
XDisplayHeightMM = _lib.XDisplayHeightMM
XDisplayHeightMM.restype = c_int
XDisplayHeightMM.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:2283
XDisplayKeycodes = _lib.XDisplayKeycodes
XDisplayKeycodes.restype = c_int
XDisplayKeycodes.argtypes = [POINTER(Display), POINTER(c_int), POINTER(c_int)]

# /usr/include/X11/Xlib.h:2289
XDisplayPlanes = _lib.XDisplayPlanes
XDisplayPlanes.restype = c_int
XDisplayPlanes.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:2294
XDisplayWidth = _lib.XDisplayWidth
XDisplayWidth.restype = c_int
XDisplayWidth.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:2299
XDisplayWidthMM = _lib.XDisplayWidthMM
XDisplayWidthMM.restype = c_int
XDisplayWidthMM.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:2304
XDrawArc = _lib.XDrawArc
XDrawArc.restype = c_int
XDrawArc.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, c_uint, c_uint, c_int, c_int]

# /usr/include/X11/Xlib.h:2316
XDrawArcs = _lib.XDrawArcs
XDrawArcs.restype = c_int
XDrawArcs.argtypes = [POINTER(Display), Drawable, GC, POINTER(XArc), c_int]

# /usr/include/X11/Xlib.h:2324
XDrawImageString = _lib.XDrawImageString
XDrawImageString.restype = c_int
XDrawImageString.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, c_char_p, c_int]

# /usr/include/X11/Xlib.h:2334
XDrawImageString16 = _lib.XDrawImageString16
XDrawImageString16.restype = c_int
XDrawImageString16.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, POINTER(XChar2b), c_int]

# /usr/include/X11/Xlib.h:2344
XDrawLine = _lib.XDrawLine
XDrawLine.restype = c_int
XDrawLine.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, c_int, c_int]

# /usr/include/X11/Xlib.h:2354
XDrawLines = _lib.XDrawLines
XDrawLines.restype = c_int
XDrawLines.argtypes = [POINTER(Display), Drawable, GC, POINTER(XPoint), c_int, c_int]

# /usr/include/X11/Xlib.h:2363
XDrawPoint = _lib.XDrawPoint
XDrawPoint.restype = c_int
XDrawPoint.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int]

# /usr/include/X11/Xlib.h:2371
XDrawPoints = _lib.XDrawPoints
XDrawPoints.restype = c_int
XDrawPoints.argtypes = [POINTER(Display), Drawable, GC, POINTER(XPoint), c_int, c_int]

# /usr/include/X11/Xlib.h:2380
XDrawRectangle = _lib.XDrawRectangle
XDrawRectangle.restype = c_int
XDrawRectangle.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, c_uint, c_uint]

# /usr/include/X11/Xlib.h:2390
XDrawRectangles = _lib.XDrawRectangles
XDrawRectangles.restype = c_int
XDrawRectangles.argtypes = [POINTER(Display), Drawable, GC, POINTER(XRectangle), c_int]

# /usr/include/X11/Xlib.h:2398
XDrawSegments = _lib.XDrawSegments
XDrawSegments.restype = c_int
XDrawSegments.argtypes = [POINTER(Display), Drawable, GC, POINTER(XSegment), c_int]

# /usr/include/X11/Xlib.h:2406
XDrawString = _lib.XDrawString
XDrawString.restype = c_int
XDrawString.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, c_char_p, c_int]

# /usr/include/X11/Xlib.h:2416
XDrawString16 = _lib.XDrawString16
XDrawString16.restype = c_int
XDrawString16.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, POINTER(XChar2b), c_int]

# /usr/include/X11/Xlib.h:2426
XDrawText = _lib.XDrawText
XDrawText.restype = c_int
XDrawText.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, POINTER(XTextItem), c_int]

# /usr/include/X11/Xlib.h:2436
XDrawText16 = _lib.XDrawText16
XDrawText16.restype = c_int
XDrawText16.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, POINTER(XTextItem16), c_int]

# /usr/include/X11/Xlib.h:2446
XEnableAccessControl = _lib.XEnableAccessControl
XEnableAccessControl.restype = c_int
XEnableAccessControl.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2450
XEventsQueued = _lib.XEventsQueued
XEventsQueued.restype = c_int
XEventsQueued.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:2455
XFetchName = _lib.XFetchName
XFetchName.restype = c_int
XFetchName.argtypes = [POINTER(Display), Window, POINTER(c_char_p)]

# /usr/include/X11/Xlib.h:2461
XFillArc = _lib.XFillArc
XFillArc.restype = c_int
XFillArc.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, c_uint, c_uint, c_int, c_int]

# /usr/include/X11/Xlib.h:2473
XFillArcs = _lib.XFillArcs
XFillArcs.restype = c_int
XFillArcs.argtypes = [POINTER(Display), Drawable, GC, POINTER(XArc), c_int]

# /usr/include/X11/Xlib.h:2481
XFillPolygon = _lib.XFillPolygon
XFillPolygon.restype = c_int
XFillPolygon.argtypes = [POINTER(Display), Drawable, GC, POINTER(XPoint), c_int, c_int, c_int]

# /usr/include/X11/Xlib.h:2491
XFillRectangle = _lib.XFillRectangle
XFillRectangle.restype = c_int
XFillRectangle.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, c_uint, c_uint]

# /usr/include/X11/Xlib.h:2501
XFillRectangles = _lib.XFillRectangles
XFillRectangles.restype = c_int
XFillRectangles.argtypes = [POINTER(Display), Drawable, GC, POINTER(XRectangle), c_int]

# /usr/include/X11/Xlib.h:2509
XFlush = _lib.XFlush
XFlush.restype = c_int
XFlush.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2513
XForceScreenSaver = _lib.XForceScreenSaver
XForceScreenSaver.restype = c_int
XForceScreenSaver.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:2518
XFree = _lib.XFree
XFree.restype = c_int
XFree.argtypes = [POINTER(None)]

# /usr/include/X11/Xlib.h:2522
XFreeColormap = _lib.XFreeColormap
XFreeColormap.restype = c_int
XFreeColormap.argtypes = [POINTER(Display), Colormap]

# /usr/include/X11/Xlib.h:2527
XFreeColors = _lib.XFreeColors
XFreeColors.restype = c_int
XFreeColors.argtypes = [POINTER(Display), Colormap, POINTER(c_ulong), c_int, c_ulong]

# /usr/include/X11/Xlib.h:2535
XFreeCursor = _lib.XFreeCursor
XFreeCursor.restype = c_int
XFreeCursor.argtypes = [POINTER(Display), Cursor]

# /usr/include/X11/Xlib.h:2540
XFreeExtensionList = _lib.XFreeExtensionList
XFreeExtensionList.restype = c_int
XFreeExtensionList.argtypes = [POINTER(c_char_p)]

# /usr/include/X11/Xlib.h:2544
XFreeFont = _lib.XFreeFont
XFreeFont.restype = c_int
XFreeFont.argtypes = [POINTER(Display), POINTER(XFontStruct)]

# /usr/include/X11/Xlib.h:2549
XFreeFontInfo = _lib.XFreeFontInfo
XFreeFontInfo.restype = c_int
XFreeFontInfo.argtypes = [POINTER(c_char_p), POINTER(XFontStruct), c_int]

# /usr/include/X11/Xlib.h:2555
XFreeFontNames = _lib.XFreeFontNames
XFreeFontNames.restype = c_int
XFreeFontNames.argtypes = [POINTER(c_char_p)]

# /usr/include/X11/Xlib.h:2559
XFreeFontPath = _lib.XFreeFontPath
XFreeFontPath.restype = c_int
XFreeFontPath.argtypes = [POINTER(c_char_p)]

# /usr/include/X11/Xlib.h:2563
XFreeGC = _lib.XFreeGC
XFreeGC.restype = c_int
XFreeGC.argtypes = [POINTER(Display), GC]

# /usr/include/X11/Xlib.h:2568
XFreeModifiermap = _lib.XFreeModifiermap
XFreeModifiermap.restype = c_int
XFreeModifiermap.argtypes = [POINTER(XModifierKeymap)]

# /usr/include/X11/Xlib.h:2572
XFreePixmap = _lib.XFreePixmap
XFreePixmap.restype = c_int
XFreePixmap.argtypes = [POINTER(Display), Pixmap]

# /usr/include/X11/Xlib.h:2577
XGeometry = _lib.XGeometry
XGeometry.restype = c_int
XGeometry.argtypes = [POINTER(Display), c_int, c_char_p, c_char_p, c_uint, c_uint, c_uint, c_int, c_int, POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int)]

# /usr/include/X11/Xlib.h:2593
XGetErrorDatabaseText = _lib.XGetErrorDatabaseText
XGetErrorDatabaseText.restype = c_int
XGetErrorDatabaseText.argtypes = [POINTER(Display), c_char_p, c_char_p, c_char_p, c_char_p, c_int]

# /usr/include/X11/Xlib.h:2602
XGetErrorText = _lib.XGetErrorText
XGetErrorText.restype = c_int
XGetErrorText.argtypes = [POINTER(Display), c_int, c_char_p, c_int]

# /usr/include/X11/Xlib.h:2609
XGetFontProperty = _lib.XGetFontProperty
XGetFontProperty.restype = c_int
XGetFontProperty.argtypes = [POINTER(XFontStruct), Atom, POINTER(c_ulong)]

# /usr/include/X11/Xlib.h:2615
XGetGCValues = _lib.XGetGCValues
XGetGCValues.restype = c_int
XGetGCValues.argtypes = [POINTER(Display), GC, c_ulong, POINTER(XGCValues)]

# /usr/include/X11/Xlib.h:2622
XGetGeometry = _lib.XGetGeometry
XGetGeometry.restype = c_int
XGetGeometry.argtypes = [POINTER(Display), Drawable, POINTER(Window), POINTER(c_int), POINTER(c_int), POINTER(c_uint), POINTER(c_uint), POINTER(c_uint), POINTER(c_uint)]

# /usr/include/X11/Xlib.h:2634
XGetIconName = _lib.XGetIconName
XGetIconName.restype = c_int
XGetIconName.argtypes = [POINTER(Display), Window, POINTER(c_char_p)]

# /usr/include/X11/Xlib.h:2640
XGetInputFocus = _lib.XGetInputFocus
XGetInputFocus.restype = c_int
XGetInputFocus.argtypes = [POINTER(Display), POINTER(Window), POINTER(c_int)]

# /usr/include/X11/Xlib.h:2646
XGetKeyboardControl = _lib.XGetKeyboardControl
XGetKeyboardControl.restype = c_int
XGetKeyboardControl.argtypes = [POINTER(Display), POINTER(XKeyboardState)]

# /usr/include/X11/Xlib.h:2651
XGetPointerControl = _lib.XGetPointerControl
XGetPointerControl.restype = c_int
XGetPointerControl.argtypes = [POINTER(Display), POINTER(c_int), POINTER(c_int), POINTER(c_int)]

# /usr/include/X11/Xlib.h:2658
XGetPointerMapping = _lib.XGetPointerMapping
XGetPointerMapping.restype = c_int
XGetPointerMapping.argtypes = [POINTER(Display), POINTER(c_ubyte), c_int]

# /usr/include/X11/Xlib.h:2664
XGetScreenSaver = _lib.XGetScreenSaver
XGetScreenSaver.restype = c_int
XGetScreenSaver.argtypes = [POINTER(Display), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int)]

# /usr/include/X11/Xlib.h:2672
XGetTransientForHint = _lib.XGetTransientForHint
XGetTransientForHint.restype = c_int
XGetTransientForHint.argtypes = [POINTER(Display), Window, POINTER(Window)]

# /usr/include/X11/Xlib.h:2678
XGetWindowProperty = _lib.XGetWindowProperty
XGetWindowProperty.restype = c_int
XGetWindowProperty.argtypes = [POINTER(Display), Window, Atom, c_long, c_long, c_int, Atom, POINTER(Atom), POINTER(c_int), POINTER(c_ulong), POINTER(c_ulong), POINTER(POINTER(c_ubyte))]

# /usr/include/X11/Xlib.h:2693
XGetWindowAttributes = _lib.XGetWindowAttributes
XGetWindowAttributes.restype = c_int
XGetWindowAttributes.argtypes = [POINTER(Display), Window, POINTER(XWindowAttributes)]

# /usr/include/X11/Xlib.h:2699
XGrabButton = _lib.XGrabButton
XGrabButton.restype = c_int
XGrabButton.argtypes = [POINTER(Display), c_uint, c_uint, Window, c_int, c_uint, c_int, c_int, Window, Cursor]

# /usr/include/X11/Xlib.h:2712
XGrabKey = _lib.XGrabKey
XGrabKey.restype = c_int
XGrabKey.argtypes = [POINTER(Display), c_int, c_uint, Window, c_int, c_int, c_int]

# /usr/include/X11/Xlib.h:2722
XGrabKeyboard = _lib.XGrabKeyboard
XGrabKeyboard.restype = c_int
XGrabKeyboard.argtypes = [POINTER(Display), Window, c_int, c_int, c_int, Time]

# /usr/include/X11/Xlib.h:2731
XGrabPointer = _lib.XGrabPointer
XGrabPointer.restype = c_int
XGrabPointer.argtypes = [POINTER(Display), Window, c_int, c_uint, c_int, c_int, Window, Cursor, Time]

# /usr/include/X11/Xlib.h:2743
XGrabServer = _lib.XGrabServer
XGrabServer.restype = c_int
XGrabServer.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2747
XHeightMMOfScreen = _lib.XHeightMMOfScreen
XHeightMMOfScreen.restype = c_int
XHeightMMOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:2751
XHeightOfScreen = _lib.XHeightOfScreen
XHeightOfScreen.restype = c_int
XHeightOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:2755
XIfEvent = _lib.XIfEvent
XIfEvent.restype = c_int
XIfEvent.argtypes = [POINTER(Display), POINTER(XEvent), CFUNCTYPE(c_int, POINTER(Display), POINTER(XEvent), XPointer), XPointer]

# /usr/include/X11/Xlib.h:2766
XImageByteOrder = _lib.XImageByteOrder
XImageByteOrder.restype = c_int
XImageByteOrder.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2770
XInstallColormap = _lib.XInstallColormap
XInstallColormap.restype = c_int
XInstallColormap.argtypes = [POINTER(Display), Colormap]

# /usr/include/X11/Xlib.h:2775
XKeysymToKeycode = _lib.XKeysymToKeycode
XKeysymToKeycode.restype = KeyCode
XKeysymToKeycode.argtypes = [POINTER(Display), KeySym]

# /usr/include/X11/Xlib.h:2780
XKillClient = _lib.XKillClient
XKillClient.restype = c_int
XKillClient.argtypes = [POINTER(Display), XID]

# /usr/include/X11/Xlib.h:2785
XLookupColor = _lib.XLookupColor
XLookupColor.restype = c_int
XLookupColor.argtypes = [POINTER(Display), Colormap, c_char_p, POINTER(XColor), POINTER(XColor)]

# /usr/include/X11/Xlib.h:2793
XLowerWindow = _lib.XLowerWindow
XLowerWindow.restype = c_int
XLowerWindow.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:2798
XMapRaised = _lib.XMapRaised
XMapRaised.restype = c_int
XMapRaised.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:2803
XMapSubwindows = _lib.XMapSubwindows
XMapSubwindows.restype = c_int
XMapSubwindows.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:2808
XMapWindow = _lib.XMapWindow
XMapWindow.restype = c_int
XMapWindow.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:2813
XMaskEvent = _lib.XMaskEvent
XMaskEvent.restype = c_int
XMaskEvent.argtypes = [POINTER(Display), c_long, POINTER(XEvent)]

# /usr/include/X11/Xlib.h:2819
XMaxCmapsOfScreen = _lib.XMaxCmapsOfScreen
XMaxCmapsOfScreen.restype = c_int
XMaxCmapsOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:2823
XMinCmapsOfScreen = _lib.XMinCmapsOfScreen
XMinCmapsOfScreen.restype = c_int
XMinCmapsOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:2827
XMoveResizeWindow = _lib.XMoveResizeWindow
XMoveResizeWindow.restype = c_int
XMoveResizeWindow.argtypes = [POINTER(Display), Window, c_int, c_int, c_uint, c_uint]

# /usr/include/X11/Xlib.h:2836
XMoveWindow = _lib.XMoveWindow
XMoveWindow.restype = c_int
XMoveWindow.argtypes = [POINTER(Display), Window, c_int, c_int]

# /usr/include/X11/Xlib.h:2843
XNextEvent = _lib.XNextEvent
XNextEvent.restype = c_int
XNextEvent.argtypes = [POINTER(Display), POINTER(XEvent)]

# /usr/include/X11/Xlib.h:2848
XNoOp = _lib.XNoOp
XNoOp.restype = c_int
XNoOp.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2852
XParseColor = _lib.XParseColor
XParseColor.restype = c_int
XParseColor.argtypes = [POINTER(Display), Colormap, c_char_p, POINTER(XColor)]

# /usr/include/X11/Xlib.h:2859
XParseGeometry = _lib.XParseGeometry
XParseGeometry.restype = c_int
XParseGeometry.argtypes = [c_char_p, POINTER(c_int), POINTER(c_int), POINTER(c_uint), POINTER(c_uint)]

# /usr/include/X11/Xlib.h:2867
XPeekEvent = _lib.XPeekEvent
XPeekEvent.restype = c_int
XPeekEvent.argtypes = [POINTER(Display), POINTER(XEvent)]

# /usr/include/X11/Xlib.h:2872
XPeekIfEvent = _lib.XPeekIfEvent
XPeekIfEvent.restype = c_int
XPeekIfEvent.argtypes = [POINTER(Display), POINTER(XEvent), CFUNCTYPE(c_int, POINTER(Display), POINTER(XEvent), XPointer), XPointer]

# /usr/include/X11/Xlib.h:2883
XPending = _lib.XPending
XPending.restype = c_int
XPending.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2887
XPlanesOfScreen = _lib.XPlanesOfScreen
XPlanesOfScreen.restype = c_int
XPlanesOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:2891
XProtocolRevision = _lib.XProtocolRevision
XProtocolRevision.restype = c_int
XProtocolRevision.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2895
XProtocolVersion = _lib.XProtocolVersion
XProtocolVersion.restype = c_int
XProtocolVersion.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2900
XPutBackEvent = _lib.XPutBackEvent
XPutBackEvent.restype = c_int
XPutBackEvent.argtypes = [POINTER(Display), POINTER(XEvent)]

# /usr/include/X11/Xlib.h:2905
XPutImage = _lib.XPutImage
XPutImage.restype = c_int
XPutImage.argtypes = [POINTER(Display), Drawable, GC, POINTER(XImage), c_int, c_int, c_int, c_int, c_uint, c_uint]

# /usr/include/X11/Xlib.h:2918
XQLength = _lib.XQLength
XQLength.restype = c_int
XQLength.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2922
XQueryBestCursor = _lib.XQueryBestCursor
XQueryBestCursor.restype = c_int
XQueryBestCursor.argtypes = [POINTER(Display), Drawable, c_uint, c_uint, POINTER(c_uint), POINTER(c_uint)]

# /usr/include/X11/Xlib.h:2931
XQueryBestSize = _lib.XQueryBestSize
XQueryBestSize.restype = c_int
XQueryBestSize.argtypes = [POINTER(Display), c_int, Drawable, c_uint, c_uint, POINTER(c_uint), POINTER(c_uint)]

# /usr/include/X11/Xlib.h:2941
XQueryBestStipple = _lib.XQueryBestStipple
XQueryBestStipple.restype = c_int
XQueryBestStipple.argtypes = [POINTER(Display), Drawable, c_uint, c_uint, POINTER(c_uint), POINTER(c_uint)]

# /usr/include/X11/Xlib.h:2950
XQueryBestTile = _lib.XQueryBestTile
XQueryBestTile.restype = c_int
XQueryBestTile.argtypes = [POINTER(Display), Drawable, c_uint, c_uint, POINTER(c_uint), POINTER(c_uint)]

# /usr/include/X11/Xlib.h:2959
XQueryColor = _lib.XQueryColor
XQueryColor.restype = c_int
XQueryColor.argtypes = [POINTER(Display), Colormap, POINTER(XColor)]

# /usr/include/X11/Xlib.h:2965
XQueryColors = _lib.XQueryColors
XQueryColors.restype = c_int
XQueryColors.argtypes = [POINTER(Display), Colormap, POINTER(XColor), c_int]

# /usr/include/X11/Xlib.h:2972
XQueryExtension = _lib.XQueryExtension
XQueryExtension.restype = c_int
XQueryExtension.argtypes = [POINTER(Display), c_char_p, POINTER(c_int), POINTER(c_int), POINTER(c_int)]

# /usr/include/X11/Xlib.h:2980
XQueryKeymap = _lib.XQueryKeymap
XQueryKeymap.restype = c_int
XQueryKeymap.argtypes = [POINTER(Display), c_char * 32]

# /usr/include/X11/Xlib.h:2985
XQueryPointer = _lib.XQueryPointer
XQueryPointer.restype = c_int
XQueryPointer.argtypes = [POINTER(Display), Window, POINTER(Window), POINTER(Window), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_uint)]

# /usr/include/X11/Xlib.h:2997
XQueryTextExtents = _lib.XQueryTextExtents
XQueryTextExtents.restype = c_int
XQueryTextExtents.argtypes = [POINTER(Display), XID, c_char_p, c_int, POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(XCharStruct)]

# /usr/include/X11/Xlib.h:3008
XQueryTextExtents16 = _lib.XQueryTextExtents16
XQueryTextExtents16.restype = c_int
XQueryTextExtents16.argtypes = [POINTER(Display), XID, POINTER(XChar2b), c_int, POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(XCharStruct)]

# /usr/include/X11/Xlib.h:3019
XQueryTree = _lib.XQueryTree
XQueryTree.restype = c_int
XQueryTree.argtypes = [POINTER(Display), Window, POINTER(Window), POINTER(Window), POINTER(POINTER(Window)), POINTER(c_uint)]

# /usr/include/X11/Xlib.h:3028
XRaiseWindow = _lib.XRaiseWindow
XRaiseWindow.restype = c_int
XRaiseWindow.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:3033
XReadBitmapFile = _lib.XReadBitmapFile
XReadBitmapFile.restype = c_int
XReadBitmapFile.argtypes = [POINTER(Display), Drawable, c_char_p, POINTER(c_uint), POINTER(c_uint), POINTER(Pixmap), POINTER(c_int), POINTER(c_int)]

# /usr/include/X11/Xlib.h:3044
XReadBitmapFileData = _lib.XReadBitmapFileData
XReadBitmapFileData.restype = c_int
XReadBitmapFileData.argtypes = [c_char_p, POINTER(c_uint), POINTER(c_uint), POINTER(POINTER(c_ubyte)), POINTER(c_int), POINTER(c_int)]

# /usr/include/X11/Xlib.h:3053
XRebindKeysym = _lib.XRebindKeysym
XRebindKeysym.restype = c_int
XRebindKeysym.argtypes = [POINTER(Display), KeySym, POINTER(KeySym), c_int, POINTER(c_ubyte), c_int]

# /usr/include/X11/Xlib.h:3062
XRecolorCursor = _lib.XRecolorCursor
XRecolorCursor.restype = c_int
XRecolorCursor.argtypes = [POINTER(Display), Cursor, POINTER(XColor), POINTER(XColor)]

# /usr/include/X11/Xlib.h:3069
XRefreshKeyboardMapping = _lib.XRefreshKeyboardMapping
XRefreshKeyboardMapping.restype = c_int
XRefreshKeyboardMapping.argtypes = [POINTER(XMappingEvent)]

# /usr/include/X11/Xlib.h:3073
XRemoveFromSaveSet = _lib.XRemoveFromSaveSet
XRemoveFromSaveSet.restype = c_int
XRemoveFromSaveSet.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:3078
XRemoveHost = _lib.XRemoveHost
XRemoveHost.restype = c_int
XRemoveHost.argtypes = [POINTER(Display), POINTER(XHostAddress)]

# /usr/include/X11/Xlib.h:3083
XRemoveHosts = _lib.XRemoveHosts
XRemoveHosts.restype = c_int
XRemoveHosts.argtypes = [POINTER(Display), POINTER(XHostAddress), c_int]

# /usr/include/X11/Xlib.h:3089
XReparentWindow = _lib.XReparentWindow
XReparentWindow.restype = c_int
XReparentWindow.argtypes = [POINTER(Display), Window, Window, c_int, c_int]

# /usr/include/X11/Xlib.h:3097
XResetScreenSaver = _lib.XResetScreenSaver
XResetScreenSaver.restype = c_int
XResetScreenSaver.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:3101
XResizeWindow = _lib.XResizeWindow
XResizeWindow.restype = c_int
XResizeWindow.argtypes = [POINTER(Display), Window, c_uint, c_uint]

# /usr/include/X11/Xlib.h:3108
XRestackWindows = _lib.XRestackWindows
XRestackWindows.restype = c_int
XRestackWindows.argtypes = [POINTER(Display), POINTER(Window), c_int]

# /usr/include/X11/Xlib.h:3114
XRotateBuffers = _lib.XRotateBuffers
XRotateBuffers.restype = c_int
XRotateBuffers.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:3119
XRotateWindowProperties = _lib.XRotateWindowProperties
XRotateWindowProperties.restype = c_int
XRotateWindowProperties.argtypes = [POINTER(Display), Window, POINTER(Atom), c_int, c_int]

# /usr/include/X11/Xlib.h:3127
XScreenCount = _lib.XScreenCount
XScreenCount.restype = c_int
XScreenCount.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:3131
XSelectInput = _lib.XSelectInput
XSelectInput.restype = c_int
XSelectInput.argtypes = [POINTER(Display), Window, c_long]

# /usr/include/X11/Xlib.h:3137
XSendEvent = _lib.XSendEvent
XSendEvent.restype = c_int
XSendEvent.argtypes = [POINTER(Display), Window, c_int, c_long, POINTER(XEvent)]

# /usr/include/X11/Xlib.h:3145
XSetAccessControl = _lib.XSetAccessControl
XSetAccessControl.restype = c_int
XSetAccessControl.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:3150
XSetArcMode = _lib.XSetArcMode
XSetArcMode.restype = c_int
XSetArcMode.argtypes = [POINTER(Display), GC, c_int]

# /usr/include/X11/Xlib.h:3156
XSetBackground = _lib.XSetBackground
XSetBackground.restype = c_int
XSetBackground.argtypes = [POINTER(Display), GC, c_ulong]

# /usr/include/X11/Xlib.h:3162
XSetClipMask = _lib.XSetClipMask
XSetClipMask.restype = c_int
XSetClipMask.argtypes = [POINTER(Display), GC, Pixmap]

# /usr/include/X11/Xlib.h:3168
XSetClipOrigin = _lib.XSetClipOrigin
XSetClipOrigin.restype = c_int
XSetClipOrigin.argtypes = [POINTER(Display), GC, c_int, c_int]

# /usr/include/X11/Xlib.h:3175
XSetClipRectangles = _lib.XSetClipRectangles
XSetClipRectangles.restype = c_int
XSetClipRectangles.argtypes = [POINTER(Display), GC, c_int, c_int, POINTER(XRectangle), c_int, c_int]

# /usr/include/X11/Xlib.h:3185
XSetCloseDownMode = _lib.XSetCloseDownMode
XSetCloseDownMode.restype = c_int
XSetCloseDownMode.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:3190
XSetCommand = _lib.XSetCommand
XSetCommand.restype = c_int
XSetCommand.argtypes = [POINTER(Display), Window, POINTER(c_char_p), c_int]

# /usr/include/X11/Xlib.h:3197
XSetDashes = _lib.XSetDashes
XSetDashes.restype = c_int
XSetDashes.argtypes = [POINTER(Display), GC, c_int, c_char_p, c_int]

# /usr/include/X11/Xlib.h:3205
XSetFillRule = _lib.XSetFillRule
XSetFillRule.restype = c_int
XSetFillRule.argtypes = [POINTER(Display), GC, c_int]

# /usr/include/X11/Xlib.h:3211
XSetFillStyle = _lib.XSetFillStyle
XSetFillStyle.restype = c_int
XSetFillStyle.argtypes = [POINTER(Display), GC, c_int]

# /usr/include/X11/Xlib.h:3217
XSetFont = _lib.XSetFont
XSetFont.restype = c_int
XSetFont.argtypes = [POINTER(Display), GC, Font]

# /usr/include/X11/Xlib.h:3223
XSetFontPath = _lib.XSetFontPath
XSetFontPath.restype = c_int
XSetFontPath.argtypes = [POINTER(Display), POINTER(c_char_p), c_int]

# /usr/include/X11/Xlib.h:3229
XSetForeground = _lib.XSetForeground
XSetForeground.restype = c_int
XSetForeground.argtypes = [POINTER(Display), GC, c_ulong]

# /usr/include/X11/Xlib.h:3235
XSetFunction = _lib.XSetFunction
XSetFunction.restype = c_int
XSetFunction.argtypes = [POINTER(Display), GC, c_int]

# /usr/include/X11/Xlib.h:3241
XSetGraphicsExposures = _lib.XSetGraphicsExposures
XSetGraphicsExposures.restype = c_int
XSetGraphicsExposures.argtypes = [POINTER(Display), GC, c_int]

# /usr/include/X11/Xlib.h:3247
XSetIconName = _lib.XSetIconName
XSetIconName.restype = c_int
XSetIconName.argtypes = [POINTER(Display), Window, c_char_p]

# /usr/include/X11/Xlib.h:3253
XSetInputFocus = _lib.XSetInputFocus
XSetInputFocus.restype = c_int
XSetInputFocus.argtypes = [POINTER(Display), Window, c_int, Time]

# /usr/include/X11/Xlib.h:3260
XSetLineAttributes = _lib.XSetLineAttributes
XSetLineAttributes.restype = c_int
XSetLineAttributes.argtypes = [POINTER(Display), GC, c_uint, c_int, c_int, c_int]

# /usr/include/X11/Xlib.h:3269
XSetModifierMapping = _lib.XSetModifierMapping
XSetModifierMapping.restype = c_int
XSetModifierMapping.argtypes = [POINTER(Display), POINTER(XModifierKeymap)]

# /usr/include/X11/Xlib.h:3274
XSetPlaneMask = _lib.XSetPlaneMask
XSetPlaneMask.restype = c_int
XSetPlaneMask.argtypes = [POINTER(Display), GC, c_ulong]

# /usr/include/X11/Xlib.h:3280
XSetPointerMapping = _lib.XSetPointerMapping
XSetPointerMapping.restype = c_int
XSetPointerMapping.argtypes = [POINTER(Display), POINTER(c_ubyte), c_int]

# /usr/include/X11/Xlib.h:3286
XSetScreenSaver = _lib.XSetScreenSaver
XSetScreenSaver.restype = c_int
XSetScreenSaver.argtypes = [POINTER(Display), c_int, c_int, c_int, c_int]

# /usr/include/X11/Xlib.h:3294
XSetSelectionOwner = _lib.XSetSelectionOwner
XSetSelectionOwner.restype = c_int
XSetSelectionOwner.argtypes = [POINTER(Display), Atom, Window, Time]

# /usr/include/X11/Xlib.h:3301
XSetState = _lib.XSetState
XSetState.restype = c_int
XSetState.argtypes = [POINTER(Display), GC, c_ulong, c_ulong, c_int, c_ulong]

# /usr/include/X11/Xlib.h:3310
XSetStipple = _lib.XSetStipple
XSetStipple.restype = c_int
XSetStipple.argtypes = [POINTER(Display), GC, Pixmap]

# /usr/include/X11/Xlib.h:3316
XSetSubwindowMode = _lib.XSetSubwindowMode
XSetSubwindowMode.restype = c_int
XSetSubwindowMode.argtypes = [POINTER(Display), GC, c_int]

# /usr/include/X11/Xlib.h:3322
XSetTSOrigin = _lib.XSetTSOrigin
XSetTSOrigin.restype = c_int
XSetTSOrigin.argtypes = [POINTER(Display), GC, c_int, c_int]

# /usr/include/X11/Xlib.h:3329
XSetTile = _lib.XSetTile
XSetTile.restype = c_int
XSetTile.argtypes = [POINTER(Display), GC, Pixmap]

# /usr/include/X11/Xlib.h:3335
XSetWindowBackground = _lib.XSetWindowBackground
XSetWindowBackground.restype = c_int
XSetWindowBackground.argtypes = [POINTER(Display), Window, c_ulong]

# /usr/include/X11/Xlib.h:3341
XSetWindowBackgroundPixmap = _lib.XSetWindowBackgroundPixmap
XSetWindowBackgroundPixmap.restype = c_int
XSetWindowBackgroundPixmap.argtypes = [POINTER(Display), Window, Pixmap]

# /usr/include/X11/Xlib.h:3347
XSetWindowBorder = _lib.XSetWindowBorder
XSetWindowBorder.restype = c_int
XSetWindowBorder.argtypes = [POINTER(Display), Window, c_ulong]

# /usr/include/X11/Xlib.h:3353
XSetWindowBorderPixmap = _lib.XSetWindowBorderPixmap
XSetWindowBorderPixmap.restype = c_int
XSetWindowBorderPixmap.argtypes = [POINTER(Display), Window, Pixmap]

# /usr/include/X11/Xlib.h:3359
XSetWindowBorderWidth = _lib.XSetWindowBorderWidth
XSetWindowBorderWidth.restype = c_int
XSetWindowBorderWidth.argtypes = [POINTER(Display), Window, c_uint]

# /usr/include/X11/Xlib.h:3365
XSetWindowColormap = _lib.XSetWindowColormap
XSetWindowColormap.restype = c_int
XSetWindowColormap.argtypes = [POINTER(Display), Window, Colormap]

# /usr/include/X11/Xlib.h:3371
XStoreBuffer = _lib.XStoreBuffer
XStoreBuffer.restype = c_int
XStoreBuffer.argtypes = [POINTER(Display), c_char_p, c_int, c_int]

# /usr/include/X11/Xlib.h:3378
XStoreBytes = _lib.XStoreBytes
XStoreBytes.restype = c_int
XStoreBytes.argtypes = [POINTER(Display), c_char_p, c_int]

# /usr/include/X11/Xlib.h:3384
XStoreColor = _lib.XStoreColor
XStoreColor.restype = c_int
XStoreColor.argtypes = [POINTER(Display), Colormap, POINTER(XColor)]

# /usr/include/X11/Xlib.h:3390
XStoreColors = _lib.XStoreColors
XStoreColors.restype = c_int
XStoreColors.argtypes = [POINTER(Display), Colormap, POINTER(XColor), c_int]

# /usr/include/X11/Xlib.h:3397
XStoreName = _lib.XStoreName
XStoreName.restype = c_int
XStoreName.argtypes = [POINTER(Display), Window, c_char_p]

# /usr/include/X11/Xlib.h:3403
XStoreNamedColor = _lib.XStoreNamedColor
XStoreNamedColor.restype = c_int
XStoreNamedColor.argtypes = [POINTER(Display), Colormap, c_char_p, c_ulong, c_int]

# /usr/include/X11/Xlib.h:3411
XSync = _lib.XSync
XSync.restype = c_int
XSync.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:3416
XTextExtents = _lib.XTextExtents
XTextExtents.restype = c_int
XTextExtents.argtypes = [POINTER(XFontStruct), c_char_p, c_int, POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(XCharStruct)]

# /usr/include/X11/Xlib.h:3426
XTextExtents16 = _lib.XTextExtents16
XTextExtents16.restype = c_int
XTextExtents16.argtypes = [POINTER(XFontStruct), POINTER(XChar2b), c_int, POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(XCharStruct)]

# /usr/include/X11/Xlib.h:3436
XTextWidth = _lib.XTextWidth
XTextWidth.restype = c_int
XTextWidth.argtypes = [POINTER(XFontStruct), c_char_p, c_int]

# /usr/include/X11/Xlib.h:3442
XTextWidth16 = _lib.XTextWidth16
XTextWidth16.restype = c_int
XTextWidth16.argtypes = [POINTER(XFontStruct), POINTER(XChar2b), c_int]

# /usr/include/X11/Xlib.h:3448
XTranslateCoordinates = _lib.XTranslateCoordinates
XTranslateCoordinates.restype = c_int
XTranslateCoordinates.argtypes = [POINTER(Display), Window, Window, c_int, c_int, POINTER(c_int), POINTER(c_int), POINTER(Window)]

# /usr/include/X11/Xlib.h:3459
XUndefineCursor = _lib.XUndefineCursor
XUndefineCursor.restype = c_int
XUndefineCursor.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:3464
XUngrabButton = _lib.XUngrabButton
XUngrabButton.restype = c_int
XUngrabButton.argtypes = [POINTER(Display), c_uint, c_uint, Window]

# /usr/include/X11/Xlib.h:3471
XUngrabKey = _lib.XUngrabKey
XUngrabKey.restype = c_int
XUngrabKey.argtypes = [POINTER(Display), c_int, c_uint, Window]

# /usr/include/X11/Xlib.h:3478
XUngrabKeyboard = _lib.XUngrabKeyboard
XUngrabKeyboard.restype = c_int
XUngrabKeyboard.argtypes = [POINTER(Display), Time]

# /usr/include/X11/Xlib.h:3483
XUngrabPointer = _lib.XUngrabPointer
XUngrabPointer.restype = c_int
XUngrabPointer.argtypes = [POINTER(Display), Time]

# /usr/include/X11/Xlib.h:3488
XUngrabServer = _lib.XUngrabServer
XUngrabServer.restype = c_int
XUngrabServer.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:3492
XUninstallColormap = _lib.XUninstallColormap
XUninstallColormap.restype = c_int
XUninstallColormap.argtypes = [POINTER(Display), Colormap]

# /usr/include/X11/Xlib.h:3497
XUnloadFont = _lib.XUnloadFont
XUnloadFont.restype = c_int
XUnloadFont.argtypes = [POINTER(Display), Font]

# /usr/include/X11/Xlib.h:3502
XUnmapSubwindows = _lib.XUnmapSubwindows
XUnmapSubwindows.restype = c_int
XUnmapSubwindows.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:3507
XUnmapWindow = _lib.XUnmapWindow
XUnmapWindow.restype = c_int
XUnmapWindow.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:3512
XVendorRelease = _lib.XVendorRelease
XVendorRelease.restype = c_int
XVendorRelease.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:3516
XWarpPointer = _lib.XWarpPointer
XWarpPointer.restype = c_int
XWarpPointer.argtypes = [POINTER(Display), Window, Window, c_int, c_int, c_uint, c_uint, c_int, c_int]

# /usr/include/X11/Xlib.h:3528
XWidthMMOfScreen = _lib.XWidthMMOfScreen
XWidthMMOfScreen.restype = c_int
XWidthMMOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:3532
XWidthOfScreen = _lib.XWidthOfScreen
XWidthOfScreen.restype = c_int
XWidthOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:3536
XWindowEvent = _lib.XWindowEvent
XWindowEvent.restype = c_int
XWindowEvent.argtypes = [POINTER(Display), Window, c_long, POINTER(XEvent)]

# /usr/include/X11/Xlib.h:3543
XWriteBitmapFile = _lib.XWriteBitmapFile
XWriteBitmapFile.restype = c_int
XWriteBitmapFile.argtypes = [POINTER(Display), c_char_p, Pixmap, c_uint, c_uint, c_int, c_int]

# /usr/include/X11/Xlib.h:3553
XSupportsLocale = _lib.XSupportsLocale
XSupportsLocale.restype = c_int
XSupportsLocale.argtypes = []

# /usr/include/X11/Xlib.h:3555
XSetLocaleModifiers = _lib.XSetLocaleModifiers
XSetLocaleModifiers.restype = c_char_p
XSetLocaleModifiers.argtypes = [c_char_p]

# /usr/include/X11/Xlib.h:3559
XOpenOM = _lib.XOpenOM
XOpenOM.restype = XOM
XOpenOM.argtypes = [POINTER(Display), POINTER(struct__XrmHashBucketRec), c_char_p, c_char_p]

# /usr/include/X11/Xlib.h:3566
XCloseOM = _lib.XCloseOM
XCloseOM.restype = c_int
XCloseOM.argtypes = [XOM]

# /usr/include/X11/Xlib.h:3570
XSetOMValues = _lib.XSetOMValues
XSetOMValues.restype = c_char_p
XSetOMValues.argtypes = [XOM]

# /usr/include/X11/Xlib.h:3575
XGetOMValues = _lib.XGetOMValues
XGetOMValues.restype = c_char_p
XGetOMValues.argtypes = [XOM]

# /usr/include/X11/Xlib.h:3580
XDisplayOfOM = _lib.XDisplayOfOM
XDisplayOfOM.restype = POINTER(Display)
XDisplayOfOM.argtypes = [XOM]

# /usr/include/X11/Xlib.h:3584
XLocaleOfOM = _lib.XLocaleOfOM
XLocaleOfOM.restype = c_char_p
XLocaleOfOM.argtypes = [XOM]

# /usr/include/X11/Xlib.h:3588
XCreateOC = _lib.XCreateOC
XCreateOC.restype = XOC
XCreateOC.argtypes = [XOM]

# /usr/include/X11/Xlib.h:3593
XDestroyOC = _lib.XDestroyOC
XDestroyOC.restype = None
XDestroyOC.argtypes = [XOC]

# /usr/include/X11/Xlib.h:3597
XOMOfOC = _lib.XOMOfOC
XOMOfOC.restype = XOM
XOMOfOC.argtypes = [XOC]

# /usr/include/X11/Xlib.h:3601
XSetOCValues = _lib.XSetOCValues
XSetOCValues.restype = c_char_p
XSetOCValues.argtypes = [XOC]

# /usr/include/X11/Xlib.h:3606
XGetOCValues = _lib.XGetOCValues
XGetOCValues.restype = c_char_p
XGetOCValues.argtypes = [XOC]

# /usr/include/X11/Xlib.h:3611
XCreateFontSet = _lib.XCreateFontSet
XCreateFontSet.restype = XFontSet
XCreateFontSet.argtypes = [POINTER(Display), c_char_p, POINTER(POINTER(c_char_p)), POINTER(c_int), POINTER(c_char_p)]

# /usr/include/X11/Xlib.h:3619
XFreeFontSet = _lib.XFreeFontSet
XFreeFontSet.restype = None
XFreeFontSet.argtypes = [POINTER(Display), XFontSet]

# /usr/include/X11/Xlib.h:3624
XFontsOfFontSet = _lib.XFontsOfFontSet
XFontsOfFontSet.restype = c_int
XFontsOfFontSet.argtypes = [XFontSet, POINTER(POINTER(POINTER(XFontStruct))), POINTER(POINTER(c_char_p))]

# /usr/include/X11/Xlib.h:3630
XBaseFontNameListOfFontSet = _lib.XBaseFontNameListOfFontSet
XBaseFontNameListOfFontSet.restype = c_char_p
XBaseFontNameListOfFontSet.argtypes = [XFontSet]

# /usr/include/X11/Xlib.h:3634
XLocaleOfFontSet = _lib.XLocaleOfFontSet
XLocaleOfFontSet.restype = c_char_p
XLocaleOfFontSet.argtypes = [XFontSet]

# /usr/include/X11/Xlib.h:3638
XContextDependentDrawing = _lib.XContextDependentDrawing
XContextDependentDrawing.restype = c_int
XContextDependentDrawing.argtypes = [XFontSet]

# /usr/include/X11/Xlib.h:3642
XDirectionalDependentDrawing = _lib.XDirectionalDependentDrawing
XDirectionalDependentDrawing.restype = c_int
XDirectionalDependentDrawing.argtypes = [XFontSet]

# /usr/include/X11/Xlib.h:3646
XContextualDrawing = _lib.XContextualDrawing
XContextualDrawing.restype = c_int
XContextualDrawing.argtypes = [XFontSet]

# /usr/include/X11/Xlib.h:3650
XExtentsOfFontSet = _lib.XExtentsOfFontSet
XExtentsOfFontSet.restype = POINTER(XFontSetExtents)
XExtentsOfFontSet.argtypes = [XFontSet]

# /usr/include/X11/Xlib.h:3654
XmbTextEscapement = _lib.XmbTextEscapement
XmbTextEscapement.restype = c_int
XmbTextEscapement.argtypes = [XFontSet, c_char_p, c_int]

# /usr/include/X11/Xlib.h:3660
XwcTextEscapement = _lib.XwcTextEscapement
XwcTextEscapement.restype = c_int
XwcTextEscapement.argtypes = [XFontSet, c_wchar_p, c_int]

# /usr/include/X11/Xlib.h:3666
Xutf8TextEscapement = _lib.Xutf8TextEscapement
Xutf8TextEscapement.restype = c_int
Xutf8TextEscapement.argtypes = [XFontSet, c_char_p, c_int]

# /usr/include/X11/Xlib.h:3672
XmbTextExtents = _lib.XmbTextExtents
XmbTextExtents.restype = c_int
XmbTextExtents.argtypes = [XFontSet, c_char_p, c_int, POINTER(XRectangle), POINTER(XRectangle)]

# /usr/include/X11/Xlib.h:3680
XwcTextExtents = _lib.XwcTextExtents
XwcTextExtents.restype = c_int
XwcTextExtents.argtypes = [XFontSet, c_wchar_p, c_int, POINTER(XRectangle), POINTER(XRectangle)]

# /usr/include/X11/Xlib.h:3688
Xutf8TextExtents = _lib.Xutf8TextExtents
Xutf8TextExtents.restype = c_int
Xutf8TextExtents.argtypes = [XFontSet, c_char_p, c_int, POINTER(XRectangle), POINTER(XRectangle)]

# /usr/include/X11/Xlib.h:3696
XmbTextPerCharExtents = _lib.XmbTextPerCharExtents
XmbTextPerCharExtents.restype = c_int
XmbTextPerCharExtents.argtypes = [XFontSet, c_char_p, c_int, POINTER(XRectangle), POINTER(XRectangle), c_int, POINTER(c_int), POINTER(XRectangle), POINTER(XRectangle)]

# /usr/include/X11/Xlib.h:3708
XwcTextPerCharExtents = _lib.XwcTextPerCharExtents
XwcTextPerCharExtents.restype = c_int
XwcTextPerCharExtents.argtypes = [XFontSet, c_wchar_p, c_int, POINTER(XRectangle), POINTER(XRectangle), c_int, POINTER(c_int), POINTER(XRectangle), POINTER(XRectangle)]

# /usr/include/X11/Xlib.h:3720
Xutf8TextPerCharExtents = _lib.Xutf8TextPerCharExtents
Xutf8TextPerCharExtents.restype = c_int
Xutf8TextPerCharExtents.argtypes = [XFontSet, c_char_p, c_int, POINTER(XRectangle), POINTER(XRectangle), c_int, POINTER(c_int), POINTER(XRectangle), POINTER(XRectangle)]

# /usr/include/X11/Xlib.h:3732
XmbDrawText = _lib.XmbDrawText
XmbDrawText.restype = None
XmbDrawText.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, POINTER(XmbTextItem), c_int]

# /usr/include/X11/Xlib.h:3742
XwcDrawText = _lib.XwcDrawText
XwcDrawText.restype = None
XwcDrawText.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, POINTER(XwcTextItem), c_int]

# /usr/include/X11/Xlib.h:3752
Xutf8DrawText = _lib.Xutf8DrawText
Xutf8DrawText.restype = None
Xutf8DrawText.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, POINTER(XmbTextItem), c_int]

# /usr/include/X11/Xlib.h:3762
XmbDrawString = _lib.XmbDrawString
XmbDrawString.restype = None
XmbDrawString.argtypes = [POINTER(Display), Drawable, XFontSet, GC, c_int, c_int, c_char_p, c_int]

# /usr/include/X11/Xlib.h:3773
XwcDrawString = _lib.XwcDrawString
XwcDrawString.restype = None
XwcDrawString.argtypes = [POINTER(Display), Drawable, XFontSet, GC, c_int, c_int, c_wchar_p, c_int]

# /usr/include/X11/Xlib.h:3784
Xutf8DrawString = _lib.Xutf8DrawString
Xutf8DrawString.restype = None
Xutf8DrawString.argtypes = [POINTER(Display), Drawable, XFontSet, GC, c_int, c_int, c_char_p, c_int]

# /usr/include/X11/Xlib.h:3795
XmbDrawImageString = _lib.XmbDrawImageString
XmbDrawImageString.restype = None
XmbDrawImageString.argtypes = [POINTER(Display), Drawable, XFontSet, GC, c_int, c_int, c_char_p, c_int]

# /usr/include/X11/Xlib.h:3806
XwcDrawImageString = _lib.XwcDrawImageString
XwcDrawImageString.restype = None
XwcDrawImageString.argtypes = [POINTER(Display), Drawable, XFontSet, GC, c_int, c_int, c_wchar_p, c_int]

# /usr/include/X11/Xlib.h:3817
Xutf8DrawImageString = _lib.Xutf8DrawImageString
Xutf8DrawImageString.restype = None
Xutf8DrawImageString.argtypes = [POINTER(Display), Drawable, XFontSet, GC, c_int, c_int, c_char_p, c_int]

# /usr/include/X11/Xlib.h:3828
XOpenIM = _lib.XOpenIM
XOpenIM.restype = XIM
XOpenIM.argtypes = [POINTER(Display), POINTER(struct__XrmHashBucketRec), c_char_p, c_char_p]

# /usr/include/X11/Xlib.h:3835
XCloseIM = _lib.XCloseIM
XCloseIM.restype = c_int
XCloseIM.argtypes = [XIM]

# /usr/include/X11/Xlib.h:3839
XGetIMValues = _lib.XGetIMValues
XGetIMValues.restype = c_char_p
XGetIMValues.argtypes = [XIM]

# /usr/include/X11/Xlib.h:3843
XSetIMValues = _lib.XSetIMValues
XSetIMValues.restype = c_char_p
XSetIMValues.argtypes = [XIM]

# /usr/include/X11/Xlib.h:3847
XDisplayOfIM = _lib.XDisplayOfIM
XDisplayOfIM.restype = POINTER(Display)
XDisplayOfIM.argtypes = [XIM]

# /usr/include/X11/Xlib.h:3851
XLocaleOfIM = _lib.XLocaleOfIM
XLocaleOfIM.restype = c_char_p
XLocaleOfIM.argtypes = [XIM]

# /usr/include/X11/Xlib.h:3855
XCreateIC = _lib.XCreateIC
XCreateIC.restype = XIC
XCreateIC.argtypes = [XIM]

# /usr/include/X11/Xlib.h:3859
XDestroyIC = _lib.XDestroyIC
XDestroyIC.restype = None
XDestroyIC.argtypes = [XIC]

# /usr/include/X11/Xlib.h:3863
XSetICFocus = _lib.XSetICFocus
XSetICFocus.restype = None
XSetICFocus.argtypes = [XIC]

# /usr/include/X11/Xlib.h:3867
XUnsetICFocus = _lib.XUnsetICFocus
XUnsetICFocus.restype = None
XUnsetICFocus.argtypes = [XIC]

# /usr/include/X11/Xlib.h:3871
XwcResetIC = _lib.XwcResetIC
XwcResetIC.restype = c_wchar_p
XwcResetIC.argtypes = [XIC]

# /usr/include/X11/Xlib.h:3875
XmbResetIC = _lib.XmbResetIC
XmbResetIC.restype = c_char_p
XmbResetIC.argtypes = [XIC]

# /usr/include/X11/Xlib.h:3879
Xutf8ResetIC = _lib.Xutf8ResetIC
Xutf8ResetIC.restype = c_char_p
Xutf8ResetIC.argtypes = [XIC]

# /usr/include/X11/Xlib.h:3883
XSetICValues = _lib.XSetICValues
XSetICValues.restype = c_char_p
XSetICValues.argtypes = [XIC]

# /usr/include/X11/Xlib.h:3887
XGetICValues = _lib.XGetICValues
XGetICValues.restype = c_char_p
XGetICValues.argtypes = [XIC]

# /usr/include/X11/Xlib.h:3891
XIMOfIC = _lib.XIMOfIC
XIMOfIC.restype = XIM
XIMOfIC.argtypes = [XIC]

# /usr/include/X11/Xlib.h:3895
XFilterEvent = _lib.XFilterEvent
XFilterEvent.restype = c_int
XFilterEvent.argtypes = [POINTER(XEvent), Window]

# /usr/include/X11/Xlib.h:3900
XmbLookupString = _lib.XmbLookupString
XmbLookupString.restype = c_int
XmbLookupString.argtypes = [XIC, POINTER(XKeyPressedEvent), c_char_p, c_int, POINTER(KeySym), POINTER(c_int)]

# /usr/include/X11/Xlib.h:3909
XwcLookupString = _lib.XwcLookupString
XwcLookupString.restype = c_int
XwcLookupString.argtypes = [XIC, POINTER(XKeyPressedEvent), c_wchar_p, c_int, POINTER(KeySym), POINTER(c_int)]

# /usr/include/X11/Xlib.h:3918
Xutf8LookupString = _lib.Xutf8LookupString
Xutf8LookupString.restype = c_int
Xutf8LookupString.argtypes = [XIC, POINTER(XKeyPressedEvent), c_char_p, c_int, POINTER(KeySym), POINTER(c_int)]

# /usr/include/X11/Xlib.h:3927
XVaCreateNestedList = _lib.XVaCreateNestedList
XVaCreateNestedList.restype = XVaNestedList
XVaCreateNestedList.argtypes = [c_int]

# /usr/include/X11/Xlib.h:3933
XRegisterIMInstantiateCallback = _lib.XRegisterIMInstantiateCallback
XRegisterIMInstantiateCallback.restype = c_int
XRegisterIMInstantiateCallback.argtypes = [POINTER(Display), POINTER(struct__XrmHashBucketRec), c_char_p, c_char_p, XIDProc, XPointer]

# /usr/include/X11/Xlib.h:3942
XUnregisterIMInstantiateCallback = _lib.XUnregisterIMInstantiateCallback
XUnregisterIMInstantiateCallback.restype = c_int
XUnregisterIMInstantiateCallback.argtypes = [POINTER(Display), POINTER(struct__XrmHashBucketRec), c_char_p, c_char_p, XIDProc, XPointer]

XConnectionWatchProc = CFUNCTYPE(None, POINTER(Display), XPointer, c_int, c_int, POINTER(XPointer)) 	# /usr/include/X11/Xlib.h:3951
# /usr/include/X11/Xlib.h:3960
XInternalConnectionNumbers = _lib.XInternalConnectionNumbers
XInternalConnectionNumbers.restype = c_int
XInternalConnectionNumbers.argtypes = [POINTER(Display), POINTER(POINTER(c_int)), POINTER(c_int)]

# /usr/include/X11/Xlib.h:3966
XProcessInternalConnection = _lib.XProcessInternalConnection
XProcessInternalConnection.restype = None
XProcessInternalConnection.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:3971
XAddConnectionWatch = _lib.XAddConnectionWatch
XAddConnectionWatch.restype = c_int
XAddConnectionWatch.argtypes = [POINTER(Display), XConnectionWatchProc, XPointer]

# /usr/include/X11/Xlib.h:3977
XRemoveConnectionWatch = _lib.XRemoveConnectionWatch
XRemoveConnectionWatch.restype = None
XRemoveConnectionWatch.argtypes = [POINTER(Display), XConnectionWatchProc, XPointer]

# /usr/include/X11/Xlib.h:3983
XSetAuthorization = _lib.XSetAuthorization
XSetAuthorization.restype = None
XSetAuthorization.argtypes = [c_char_p, c_int, c_char_p, c_int]

# /usr/include/X11/Xlib.h:3990
_Xmbtowc = _lib._Xmbtowc
_Xmbtowc.restype = c_int
_Xmbtowc.argtypes = [c_wchar_p, c_char_p, c_int]

# /usr/include/X11/Xlib.h:4001
_Xwctomb = _lib._Xwctomb
_Xwctomb.restype = c_int
_Xwctomb.argtypes = [c_char_p, c_wchar]

X_PROTOCOL = 11 	# /usr/include/X11/X.h:4066
X_PROTOCOL_REVISION = 0 	# /usr/include/X11/X.h:4067
None_ = 0 	# /usr/include/X11/X.h:4128
ParentRelative = 1 	# /usr/include/X11/X.h:4131
CopyFromParent = 0 	# /usr/include/X11/X.h:4134
PointerWindow = 0 	# /usr/include/X11/X.h:4139
InputFocus = 1 	# /usr/include/X11/X.h:4140
PointerRoot = 1 	# /usr/include/X11/X.h:4142
AnyPropertyType = 0 	# /usr/include/X11/X.h:4144
AnyKey = 0 	# /usr/include/X11/X.h:4146
AnyButton = 0 	# /usr/include/X11/X.h:4148
AllTemporary = 0 	# /usr/include/X11/X.h:4150
CurrentTime = 0 	# /usr/include/X11/X.h:4152
NoSymbol = 0 	# /usr/include/X11/X.h:4154
NoEventMask = 0 	# /usr/include/X11/X.h:4163
KeyPressMask = 1 	# /usr/include/X11/X.h:4164
KeyReleaseMask = 2 	# /usr/include/X11/X.h:4165
ButtonPressMask = 4 	# /usr/include/X11/X.h:4166
ButtonReleaseMask = 8 	# /usr/include/X11/X.h:4167
EnterWindowMask = 16 	# /usr/include/X11/X.h:4168
LeaveWindowMask = 32 	# /usr/include/X11/X.h:4169
PointerMotionMask = 64 	# /usr/include/X11/X.h:4170
PointerMotionHintMask = 128 	# /usr/include/X11/X.h:4171
Button1MotionMask = 256 	# /usr/include/X11/X.h:4172
Button2MotionMask = 512 	# /usr/include/X11/X.h:4173
Button3MotionMask = 1024 	# /usr/include/X11/X.h:4174
Button4MotionMask = 2048 	# /usr/include/X11/X.h:4175
Button5MotionMask = 4096 	# /usr/include/X11/X.h:4176
ButtonMotionMask = 8192 	# /usr/include/X11/X.h:4177
KeymapStateMask = 16384 	# /usr/include/X11/X.h:4178
ExposureMask = 32768 	# /usr/include/X11/X.h:4179
VisibilityChangeMask = 65536 	# /usr/include/X11/X.h:4180
StructureNotifyMask = 131072 	# /usr/include/X11/X.h:4181
ResizeRedirectMask = 262144 	# /usr/include/X11/X.h:4182
SubstructureNotifyMask = 524288 	# /usr/include/X11/X.h:4183
SubstructureRedirectMask = 1048576 	# /usr/include/X11/X.h:4184
FocusChangeMask = 2097152 	# /usr/include/X11/X.h:4185
PropertyChangeMask = 4194304 	# /usr/include/X11/X.h:4186
ColormapChangeMask = 8388608 	# /usr/include/X11/X.h:4187
OwnerGrabButtonMask = 16777216 	# /usr/include/X11/X.h:4188
KeyPress = 2 	# /usr/include/X11/X.h:4194
KeyRelease = 3 	# /usr/include/X11/X.h:4195
ButtonPress = 4 	# /usr/include/X11/X.h:4196
ButtonRelease = 5 	# /usr/include/X11/X.h:4197
MotionNotify = 6 	# /usr/include/X11/X.h:4198
EnterNotify = 7 	# /usr/include/X11/X.h:4199
LeaveNotify = 8 	# /usr/include/X11/X.h:4200
FocusIn = 9 	# /usr/include/X11/X.h:4201
FocusOut = 10 	# /usr/include/X11/X.h:4202
KeymapNotify = 11 	# /usr/include/X11/X.h:4203
Expose = 12 	# /usr/include/X11/X.h:4204
GraphicsExpose = 13 	# /usr/include/X11/X.h:4205
NoExpose = 14 	# /usr/include/X11/X.h:4206
VisibilityNotify = 15 	# /usr/include/X11/X.h:4207
CreateNotify = 16 	# /usr/include/X11/X.h:4208
DestroyNotify = 17 	# /usr/include/X11/X.h:4209
UnmapNotify = 18 	# /usr/include/X11/X.h:4210
MapNotify = 19 	# /usr/include/X11/X.h:4211
MapRequest = 20 	# /usr/include/X11/X.h:4212
ReparentNotify = 21 	# /usr/include/X11/X.h:4213
ConfigureNotify = 22 	# /usr/include/X11/X.h:4214
ConfigureRequest = 23 	# /usr/include/X11/X.h:4215
GravityNotify = 24 	# /usr/include/X11/X.h:4216
ResizeRequest = 25 	# /usr/include/X11/X.h:4217
CirculateNotify = 26 	# /usr/include/X11/X.h:4218
CirculateRequest = 27 	# /usr/include/X11/X.h:4219
PropertyNotify = 28 	# /usr/include/X11/X.h:4220
SelectionClear = 29 	# /usr/include/X11/X.h:4221
SelectionRequest = 30 	# /usr/include/X11/X.h:4222
SelectionNotify = 31 	# /usr/include/X11/X.h:4223
ColormapNotify = 32 	# /usr/include/X11/X.h:4224
ClientMessage = 33 	# /usr/include/X11/X.h:4225
MappingNotify = 34 	# /usr/include/X11/X.h:4226
LASTEvent = 35 	# /usr/include/X11/X.h:4227
ShiftMask = 1 	# /usr/include/X11/X.h:4233
LockMask = 2 	# /usr/include/X11/X.h:4234
ControlMask = 4 	# /usr/include/X11/X.h:4235
Mod1Mask = 8 	# /usr/include/X11/X.h:4236
Mod2Mask = 16 	# /usr/include/X11/X.h:4237
Mod3Mask = 32 	# /usr/include/X11/X.h:4238
Mod4Mask = 64 	# /usr/include/X11/X.h:4239
Mod5Mask = 128 	# /usr/include/X11/X.h:4240
ShiftMapIndex = 0 	# /usr/include/X11/X.h:4245
LockMapIndex = 1 	# /usr/include/X11/X.h:4246
ControlMapIndex = 2 	# /usr/include/X11/X.h:4247
Mod1MapIndex = 3 	# /usr/include/X11/X.h:4248
Mod2MapIndex = 4 	# /usr/include/X11/X.h:4249
Mod3MapIndex = 5 	# /usr/include/X11/X.h:4250
Mod4MapIndex = 6 	# /usr/include/X11/X.h:4251
Mod5MapIndex = 7 	# /usr/include/X11/X.h:4252
Button1Mask = 256 	# /usr/include/X11/X.h:4258
Button2Mask = 512 	# /usr/include/X11/X.h:4259
Button3Mask = 1024 	# /usr/include/X11/X.h:4260
Button4Mask = 2048 	# /usr/include/X11/X.h:4261
Button5Mask = 4096 	# /usr/include/X11/X.h:4262
AnyModifier = 32768 	# /usr/include/X11/X.h:4264
Button1 = 1 	# /usr/include/X11/X.h:4271
Button2 = 2 	# /usr/include/X11/X.h:4272
Button3 = 3 	# /usr/include/X11/X.h:4273
Button4 = 4 	# /usr/include/X11/X.h:4274
Button5 = 5 	# /usr/include/X11/X.h:4275
NotifyNormal = 0 	# /usr/include/X11/X.h:4279
NotifyGrab = 1 	# /usr/include/X11/X.h:4280
NotifyUngrab = 2 	# /usr/include/X11/X.h:4281
NotifyWhileGrabbed = 3 	# /usr/include/X11/X.h:4282
NotifyHint = 1 	# /usr/include/X11/X.h:4284
NotifyAncestor = 0 	# /usr/include/X11/X.h:4288
NotifyVirtual = 1 	# /usr/include/X11/X.h:4289
NotifyInferior = 2 	# /usr/include/X11/X.h:4290
NotifyNonlinear = 3 	# /usr/include/X11/X.h:4291
NotifyNonlinearVirtual = 4 	# /usr/include/X11/X.h:4292
NotifyPointer = 5 	# /usr/include/X11/X.h:4293
NotifyPointerRoot = 6 	# /usr/include/X11/X.h:4294
NotifyDetailNone = 7 	# /usr/include/X11/X.h:4295
VisibilityUnobscured = 0 	# /usr/include/X11/X.h:4299
VisibilityPartiallyObscured = 1 	# /usr/include/X11/X.h:4300
VisibilityFullyObscured = 2 	# /usr/include/X11/X.h:4301
PlaceOnTop = 0 	# /usr/include/X11/X.h:4305
PlaceOnBottom = 1 	# /usr/include/X11/X.h:4306
FamilyInternet = 0 	# /usr/include/X11/X.h:4310
FamilyDECnet = 1 	# /usr/include/X11/X.h:4311
FamilyChaos = 2 	# /usr/include/X11/X.h:4312
FamilyInternet6 = 6 	# /usr/include/X11/X.h:4313
FamilyServerInterpreted = 5 	# /usr/include/X11/X.h:4316
PropertyNewValue = 0 	# /usr/include/X11/X.h:4320
PropertyDelete = 1 	# /usr/include/X11/X.h:4321
ColormapUninstalled = 0 	# /usr/include/X11/X.h:4325
ColormapInstalled = 1 	# /usr/include/X11/X.h:4326
GrabModeSync = 0 	# /usr/include/X11/X.h:4330
GrabModeAsync = 1 	# /usr/include/X11/X.h:4331
GrabSuccess = 0 	# /usr/include/X11/X.h:4335
AlreadyGrabbed = 1 	# /usr/include/X11/X.h:4336
GrabInvalidTime = 2 	# /usr/include/X11/X.h:4337
GrabNotViewable = 3 	# /usr/include/X11/X.h:4338
GrabFrozen = 4 	# /usr/include/X11/X.h:4339
AsyncPointer = 0 	# /usr/include/X11/X.h:4343
SyncPointer = 1 	# /usr/include/X11/X.h:4344
ReplayPointer = 2 	# /usr/include/X11/X.h:4345
AsyncKeyboard = 3 	# /usr/include/X11/X.h:4346
SyncKeyboard = 4 	# /usr/include/X11/X.h:4347
ReplayKeyboard = 5 	# /usr/include/X11/X.h:4348
AsyncBoth = 6 	# /usr/include/X11/X.h:4349
SyncBoth = 7 	# /usr/include/X11/X.h:4350
RevertToParent = 2 	# /usr/include/X11/X.h:4356
Success = 0 	# /usr/include/X11/X.h:4362
BadRequest = 1 	# /usr/include/X11/X.h:4363
BadValue = 2 	# /usr/include/X11/X.h:4364
BadWindow = 3 	# /usr/include/X11/X.h:4365
BadPixmap = 4 	# /usr/include/X11/X.h:4366
BadAtom = 5 	# /usr/include/X11/X.h:4367
BadCursor = 6 	# /usr/include/X11/X.h:4368
BadFont = 7 	# /usr/include/X11/X.h:4369
BadMatch = 8 	# /usr/include/X11/X.h:4370
BadDrawable = 9 	# /usr/include/X11/X.h:4371
BadAccess = 10 	# /usr/include/X11/X.h:4372
BadAlloc = 11 	# /usr/include/X11/X.h:4381
BadColor = 12 	# /usr/include/X11/X.h:4382
BadGC = 13 	# /usr/include/X11/X.h:4383
BadIDChoice = 14 	# /usr/include/X11/X.h:4384
BadName = 15 	# /usr/include/X11/X.h:4385
BadLength = 16 	# /usr/include/X11/X.h:4386
BadImplementation = 17 	# /usr/include/X11/X.h:4387
FirstExtensionError = 128 	# /usr/include/X11/X.h:4389
LastExtensionError = 255 	# /usr/include/X11/X.h:4390
InputOutput = 1 	# /usr/include/X11/X.h:4399
InputOnly = 2 	# /usr/include/X11/X.h:4400
CWBackPixmap = 1 	# /usr/include/X11/X.h:4404
CWBackPixel = 2 	# /usr/include/X11/X.h:4405
CWBorderPixmap = 4 	# /usr/include/X11/X.h:4406
CWBorderPixel = 8 	# /usr/include/X11/X.h:4407
CWBitGravity = 16 	# /usr/include/X11/X.h:4408
CWWinGravity = 32 	# /usr/include/X11/X.h:4409
CWBackingStore = 64 	# /usr/include/X11/X.h:4410
CWBackingPlanes = 128 	# /usr/include/X11/X.h:4411
CWBackingPixel = 256 	# /usr/include/X11/X.h:4412
CWOverrideRedirect = 512 	# /usr/include/X11/X.h:4413
CWSaveUnder = 1024 	# /usr/include/X11/X.h:4414
CWEventMask = 2048 	# /usr/include/X11/X.h:4415
CWDontPropagate = 4096 	# /usr/include/X11/X.h:4416
CWColormap = 8192 	# /usr/include/X11/X.h:4417
CWCursor = 16384 	# /usr/include/X11/X.h:4418
CWX = 1 	# /usr/include/X11/X.h:4422
CWY = 2 	# /usr/include/X11/X.h:4423
CWWidth = 4 	# /usr/include/X11/X.h:4424
CWHeight = 8 	# /usr/include/X11/X.h:4425
CWBorderWidth = 16 	# /usr/include/X11/X.h:4426
CWSibling = 32 	# /usr/include/X11/X.h:4427
CWStackMode = 64 	# /usr/include/X11/X.h:4428
ForgetGravity = 0 	# /usr/include/X11/X.h:4433
NorthWestGravity = 1 	# /usr/include/X11/X.h:4434
NorthGravity = 2 	# /usr/include/X11/X.h:4435
NorthEastGravity = 3 	# /usr/include/X11/X.h:4436
WestGravity = 4 	# /usr/include/X11/X.h:4437
CenterGravity = 5 	# /usr/include/X11/X.h:4438
EastGravity = 6 	# /usr/include/X11/X.h:4439
SouthWestGravity = 7 	# /usr/include/X11/X.h:4440
SouthGravity = 8 	# /usr/include/X11/X.h:4441
SouthEastGravity = 9 	# /usr/include/X11/X.h:4442
StaticGravity = 10 	# /usr/include/X11/X.h:4443
UnmapGravity = 0 	# /usr/include/X11/X.h:4447
NotUseful = 0 	# /usr/include/X11/X.h:4451
WhenMapped = 1 	# /usr/include/X11/X.h:4452
Always = 2 	# /usr/include/X11/X.h:4453
IsUnmapped = 0 	# /usr/include/X11/X.h:4457
IsUnviewable = 1 	# /usr/include/X11/X.h:4458
IsViewable = 2 	# /usr/include/X11/X.h:4459
SetModeInsert = 0 	# /usr/include/X11/X.h:4463
SetModeDelete = 1 	# /usr/include/X11/X.h:4464
DestroyAll = 0 	# /usr/include/X11/X.h:4468
RetainPermanent = 1 	# /usr/include/X11/X.h:4469
RetainTemporary = 2 	# /usr/include/X11/X.h:4470
Above = 0 	# /usr/include/X11/X.h:4474
Below = 1 	# /usr/include/X11/X.h:4475
TopIf = 2 	# /usr/include/X11/X.h:4476
BottomIf = 3 	# /usr/include/X11/X.h:4477
Opposite = 4 	# /usr/include/X11/X.h:4478
RaiseLowest = 0 	# /usr/include/X11/X.h:4482
LowerHighest = 1 	# /usr/include/X11/X.h:4483
PropModeReplace = 0 	# /usr/include/X11/X.h:4487
PropModePrepend = 1 	# /usr/include/X11/X.h:4488
PropModeAppend = 2 	# /usr/include/X11/X.h:4489
GXclear = 0 	# /usr/include/X11/X.h:4497
GXand = 1 	# /usr/include/X11/X.h:4498
GXandReverse = 2 	# /usr/include/X11/X.h:4499
GXcopy = 3 	# /usr/include/X11/X.h:4500
GXandInverted = 4 	# /usr/include/X11/X.h:4501
GXnoop = 5 	# /usr/include/X11/X.h:4502
GXxor = 6 	# /usr/include/X11/X.h:4503
GXor = 7 	# /usr/include/X11/X.h:4504
GXnor = 8 	# /usr/include/X11/X.h:4505
GXequiv = 9 	# /usr/include/X11/X.h:4506
GXinvert = 10 	# /usr/include/X11/X.h:4507
GXorReverse = 11 	# /usr/include/X11/X.h:4508
GXcopyInverted = 12 	# /usr/include/X11/X.h:4509
GXorInverted = 13 	# /usr/include/X11/X.h:4510
GXnand = 14 	# /usr/include/X11/X.h:4511
GXset = 0.0 	# /usr/include/X11/X.h:4512
LineSolid = 0 	# /usr/include/X11/X.h:4516
LineOnOffDash = 1 	# /usr/include/X11/X.h:4517
LineDoubleDash = 2 	# /usr/include/X11/X.h:4518
CapNotLast = 0 	# /usr/include/X11/X.h:4522
CapButt = 1 	# /usr/include/X11/X.h:4523
CapRound = 2 	# /usr/include/X11/X.h:4524
CapProjecting = 3 	# /usr/include/X11/X.h:4525
JoinMiter = 0 	# /usr/include/X11/X.h:4529
JoinRound = 1 	# /usr/include/X11/X.h:4530
JoinBevel = 2 	# /usr/include/X11/X.h:4531
FillSolid = 0 	# /usr/include/X11/X.h:4535
FillTiled = 1 	# /usr/include/X11/X.h:4536
FillStippled = 2 	# /usr/include/X11/X.h:4537
FillOpaqueStippled = 3 	# /usr/include/X11/X.h:4538
EvenOddRule = 0 	# /usr/include/X11/X.h:4542
WindingRule = 1 	# /usr/include/X11/X.h:4543
ClipByChildren = 0 	# /usr/include/X11/X.h:4547
IncludeInferiors = 1 	# /usr/include/X11/X.h:4548
Unsorted = 0 	# /usr/include/X11/X.h:4552
YSorted = 1 	# /usr/include/X11/X.h:4553
YXSorted = 2 	# /usr/include/X11/X.h:4554
YXBanded = 3 	# /usr/include/X11/X.h:4555
CoordModeOrigin = 0 	# /usr/include/X11/X.h:4559
CoordModePrevious = 1 	# /usr/include/X11/X.h:4560
Complex = 0 	# /usr/include/X11/X.h:4564
Nonconvex = 1 	# /usr/include/X11/X.h:4565
Convex = 2 	# /usr/include/X11/X.h:4566
ArcChord = 0 	# /usr/include/X11/X.h:4570
ArcPieSlice = 1 	# /usr/include/X11/X.h:4571
GCFunction = 1 	# /usr/include/X11/X.h:4576
GCPlaneMask = 2 	# /usr/include/X11/X.h:4577
GCForeground = 4 	# /usr/include/X11/X.h:4578
GCBackground = 8 	# /usr/include/X11/X.h:4579
GCLineWidth = 16 	# /usr/include/X11/X.h:4580
GCLineStyle = 32 	# /usr/include/X11/X.h:4581
GCCapStyle = 64 	# /usr/include/X11/X.h:4582
GCJoinStyle = 128 	# /usr/include/X11/X.h:4583
GCFillStyle = 256 	# /usr/include/X11/X.h:4584
GCFillRule = 512 	# /usr/include/X11/X.h:4585
GCTile = 1024 	# /usr/include/X11/X.h:4586
GCStipple = 2048 	# /usr/include/X11/X.h:4587
GCTileStipXOrigin = 4096 	# /usr/include/X11/X.h:4588
GCTileStipYOrigin = 8192 	# /usr/include/X11/X.h:4589
GCFont = 16384 	# /usr/include/X11/X.h:4590
GCSubwindowMode = 32768 	# /usr/include/X11/X.h:4591
GCGraphicsExposures = 65536 	# /usr/include/X11/X.h:4592
GCClipXOrigin = 131072 	# /usr/include/X11/X.h:4593
GCClipYOrigin = 262144 	# /usr/include/X11/X.h:4594
GCClipMask = 524288 	# /usr/include/X11/X.h:4595
GCDashOffset = 1048576 	# /usr/include/X11/X.h:4596
GCDashList = 2097152 	# /usr/include/X11/X.h:4597
GCArcMode = 4194304 	# /usr/include/X11/X.h:4598
GCLastBit = 22 	# /usr/include/X11/X.h:4600
FontLeftToRight = 0 	# /usr/include/X11/X.h:4607
FontRightToLeft = 1 	# /usr/include/X11/X.h:4608
FontChange = 255 	# /usr/include/X11/X.h:4610
XYBitmap = 0 	# /usr/include/X11/X.h:4618
XYPixmap = 1 	# /usr/include/X11/X.h:4619
ZPixmap = 2 	# /usr/include/X11/X.h:4620
AllocNone = 0 	# /usr/include/X11/X.h:4628
AllocAll = 1 	# /usr/include/X11/X.h:4629
DoRed = 1 	# /usr/include/X11/X.h:4634
DoGreen = 2 	# /usr/include/X11/X.h:4635
DoBlue = 4 	# /usr/include/X11/X.h:4636
CursorShape = 0 	# /usr/include/X11/X.h:4644
TileShape = 1 	# /usr/include/X11/X.h:4645
StippleShape = 2 	# /usr/include/X11/X.h:4646
AutoRepeatModeOff = 0 	# /usr/include/X11/X.h:4652
AutoRepeatModeOn = 1 	# /usr/include/X11/X.h:4653
AutoRepeatModeDefault = 2 	# /usr/include/X11/X.h:4654
LedModeOff = 0 	# /usr/include/X11/X.h:4656
LedModeOn = 1 	# /usr/include/X11/X.h:4657
KBKeyClickPercent = 1 	# /usr/include/X11/X.h:4661
KBBellPercent = 2 	# /usr/include/X11/X.h:4662
KBBellPitch = 4 	# /usr/include/X11/X.h:4663
KBBellDuration = 8 	# /usr/include/X11/X.h:4664
KBLed = 16 	# /usr/include/X11/X.h:4665
KBLedMode = 32 	# /usr/include/X11/X.h:4666
KBKey = 64 	# /usr/include/X11/X.h:4667
KBAutoRepeatMode = 128 	# /usr/include/X11/X.h:4668
MappingSuccess = 0 	# /usr/include/X11/X.h:4670
MappingBusy = 1 	# /usr/include/X11/X.h:4671
MappingFailed = 2 	# /usr/include/X11/X.h:4672
MappingModifier = 0 	# /usr/include/X11/X.h:4674
MappingKeyboard = 1 	# /usr/include/X11/X.h:4675
MappingPointer = 2 	# /usr/include/X11/X.h:4676
DontPreferBlanking = 0 	# /usr/include/X11/X.h:4682
PreferBlanking = 1 	# /usr/include/X11/X.h:4683
DefaultBlanking = 2 	# /usr/include/X11/X.h:4684
DisableScreenSaver = 0 	# /usr/include/X11/X.h:4686
DisableScreenInterval = 0 	# /usr/include/X11/X.h:4687
DontAllowExposures = 0 	# /usr/include/X11/X.h:4689
AllowExposures = 1 	# /usr/include/X11/X.h:4690
DefaultExposures = 2 	# /usr/include/X11/X.h:4691
ScreenSaverReset = 0 	# /usr/include/X11/X.h:4695
ScreenSaverActive = 1 	# /usr/include/X11/X.h:4696
HostInsert = 0 	# /usr/include/X11/X.h:4704
HostDelete = 1 	# /usr/include/X11/X.h:4705
EnableAccess = 1 	# /usr/include/X11/X.h:4709
DisableAccess = 0 	# /usr/include/X11/X.h:4710
StaticGray = 0 	# /usr/include/X11/X.h:4716
GrayScale = 1 	# /usr/include/X11/X.h:4717
StaticColor = 2 	# /usr/include/X11/X.h:4718
PseudoColor = 3 	# /usr/include/X11/X.h:4719
TrueColor = 4 	# /usr/include/X11/X.h:4720
DirectColor = 5 	# /usr/include/X11/X.h:4721
LSBFirst = 0 	# /usr/include/X11/X.h:4726
MSBFirst = 1 	# /usr/include/X11/X.h:4727
XlibSpecificationRelease = 6 	# /usr/include/X11/Xlib.h:39
X_PROTOCOL = 11 	# /usr/include/X11/X.h:58
X_PROTOCOL_REVISION = 0 	# /usr/include/X11/X.h:59
None_ = 0 	# /usr/include/X11/X.h:120
ParentRelative = 1 	# /usr/include/X11/X.h:123
CopyFromParent = 0 	# /usr/include/X11/X.h:126
PointerWindow = 0 	# /usr/include/X11/X.h:131
InputFocus = 1 	# /usr/include/X11/X.h:132
PointerRoot = 1 	# /usr/include/X11/X.h:134
AnyPropertyType = 0 	# /usr/include/X11/X.h:136
AnyKey = 0 	# /usr/include/X11/X.h:138
AnyButton = 0 	# /usr/include/X11/X.h:140
AllTemporary = 0 	# /usr/include/X11/X.h:142
CurrentTime = 0 	# /usr/include/X11/X.h:144
NoSymbol = 0 	# /usr/include/X11/X.h:146
NoEventMask = 0 	# /usr/include/X11/X.h:155
KeyPressMask = 1 	# /usr/include/X11/X.h:156
KeyReleaseMask = 2 	# /usr/include/X11/X.h:157
ButtonPressMask = 4 	# /usr/include/X11/X.h:158
ButtonReleaseMask = 8 	# /usr/include/X11/X.h:159
EnterWindowMask = 16 	# /usr/include/X11/X.h:160
LeaveWindowMask = 32 	# /usr/include/X11/X.h:161
PointerMotionMask = 64 	# /usr/include/X11/X.h:162
PointerMotionHintMask = 128 	# /usr/include/X11/X.h:163
Button1MotionMask = 256 	# /usr/include/X11/X.h:164
Button2MotionMask = 512 	# /usr/include/X11/X.h:165
Button3MotionMask = 1024 	# /usr/include/X11/X.h:166
Button4MotionMask = 2048 	# /usr/include/X11/X.h:167
Button5MotionMask = 4096 	# /usr/include/X11/X.h:168
ButtonMotionMask = 8192 	# /usr/include/X11/X.h:169
KeymapStateMask = 16384 	# /usr/include/X11/X.h:170
ExposureMask = 32768 	# /usr/include/X11/X.h:171
VisibilityChangeMask = 65536 	# /usr/include/X11/X.h:172
StructureNotifyMask = 131072 	# /usr/include/X11/X.h:173
ResizeRedirectMask = 262144 	# /usr/include/X11/X.h:174
SubstructureNotifyMask = 524288 	# /usr/include/X11/X.h:175
SubstructureRedirectMask = 1048576 	# /usr/include/X11/X.h:176
FocusChangeMask = 2097152 	# /usr/include/X11/X.h:177
PropertyChangeMask = 4194304 	# /usr/include/X11/X.h:178
ColormapChangeMask = 8388608 	# /usr/include/X11/X.h:179
OwnerGrabButtonMask = 16777216 	# /usr/include/X11/X.h:180
KeyPress = 2 	# /usr/include/X11/X.h:186
KeyRelease = 3 	# /usr/include/X11/X.h:187
ButtonPress = 4 	# /usr/include/X11/X.h:188
ButtonRelease = 5 	# /usr/include/X11/X.h:189
MotionNotify = 6 	# /usr/include/X11/X.h:190
EnterNotify = 7 	# /usr/include/X11/X.h:191
LeaveNotify = 8 	# /usr/include/X11/X.h:192
FocusIn = 9 	# /usr/include/X11/X.h:193
FocusOut = 10 	# /usr/include/X11/X.h:194
KeymapNotify = 11 	# /usr/include/X11/X.h:195
Expose = 12 	# /usr/include/X11/X.h:196
GraphicsExpose = 13 	# /usr/include/X11/X.h:197
NoExpose = 14 	# /usr/include/X11/X.h:198
VisibilityNotify = 15 	# /usr/include/X11/X.h:199
CreateNotify = 16 	# /usr/include/X11/X.h:200
DestroyNotify = 17 	# /usr/include/X11/X.h:201
UnmapNotify = 18 	# /usr/include/X11/X.h:202
MapNotify = 19 	# /usr/include/X11/X.h:203
MapRequest = 20 	# /usr/include/X11/X.h:204
ReparentNotify = 21 	# /usr/include/X11/X.h:205
ConfigureNotify = 22 	# /usr/include/X11/X.h:206
ConfigureRequest = 23 	# /usr/include/X11/X.h:207
GravityNotify = 24 	# /usr/include/X11/X.h:208
ResizeRequest = 25 	# /usr/include/X11/X.h:209
CirculateNotify = 26 	# /usr/include/X11/X.h:210
CirculateRequest = 27 	# /usr/include/X11/X.h:211
PropertyNotify = 28 	# /usr/include/X11/X.h:212
SelectionClear = 29 	# /usr/include/X11/X.h:213
SelectionRequest = 30 	# /usr/include/X11/X.h:214
SelectionNotify = 31 	# /usr/include/X11/X.h:215
ColormapNotify = 32 	# /usr/include/X11/X.h:216
ClientMessage = 33 	# /usr/include/X11/X.h:217
MappingNotify = 34 	# /usr/include/X11/X.h:218
LASTEvent = 35 	# /usr/include/X11/X.h:219
ShiftMask = 1 	# /usr/include/X11/X.h:225
LockMask = 2 	# /usr/include/X11/X.h:226
ControlMask = 4 	# /usr/include/X11/X.h:227
Mod1Mask = 8 	# /usr/include/X11/X.h:228
Mod2Mask = 16 	# /usr/include/X11/X.h:229
Mod3Mask = 32 	# /usr/include/X11/X.h:230
Mod4Mask = 64 	# /usr/include/X11/X.h:231
Mod5Mask = 128 	# /usr/include/X11/X.h:232
ShiftMapIndex = 0 	# /usr/include/X11/X.h:237
LockMapIndex = 1 	# /usr/include/X11/X.h:238
ControlMapIndex = 2 	# /usr/include/X11/X.h:239
Mod1MapIndex = 3 	# /usr/include/X11/X.h:240
Mod2MapIndex = 4 	# /usr/include/X11/X.h:241
Mod3MapIndex = 5 	# /usr/include/X11/X.h:242
Mod4MapIndex = 6 	# /usr/include/X11/X.h:243
Mod5MapIndex = 7 	# /usr/include/X11/X.h:244
Button1Mask = 256 	# /usr/include/X11/X.h:250
Button2Mask = 512 	# /usr/include/X11/X.h:251
Button3Mask = 1024 	# /usr/include/X11/X.h:252
Button4Mask = 2048 	# /usr/include/X11/X.h:253
Button5Mask = 4096 	# /usr/include/X11/X.h:254
AnyModifier = 32768 	# /usr/include/X11/X.h:256
Button1 = 1 	# /usr/include/X11/X.h:263
Button2 = 2 	# /usr/include/X11/X.h:264
Button3 = 3 	# /usr/include/X11/X.h:265
Button4 = 4 	# /usr/include/X11/X.h:266
Button5 = 5 	# /usr/include/X11/X.h:267
NotifyNormal = 0 	# /usr/include/X11/X.h:271
NotifyGrab = 1 	# /usr/include/X11/X.h:272
NotifyUngrab = 2 	# /usr/include/X11/X.h:273
NotifyWhileGrabbed = 3 	# /usr/include/X11/X.h:274
NotifyHint = 1 	# /usr/include/X11/X.h:276
NotifyAncestor = 0 	# /usr/include/X11/X.h:280
NotifyVirtual = 1 	# /usr/include/X11/X.h:281
NotifyInferior = 2 	# /usr/include/X11/X.h:282
NotifyNonlinear = 3 	# /usr/include/X11/X.h:283
NotifyNonlinearVirtual = 4 	# /usr/include/X11/X.h:284
NotifyPointer = 5 	# /usr/include/X11/X.h:285
NotifyPointerRoot = 6 	# /usr/include/X11/X.h:286
NotifyDetailNone = 7 	# /usr/include/X11/X.h:287
VisibilityUnobscured = 0 	# /usr/include/X11/X.h:291
VisibilityPartiallyObscured = 1 	# /usr/include/X11/X.h:292
VisibilityFullyObscured = 2 	# /usr/include/X11/X.h:293
PlaceOnTop = 0 	# /usr/include/X11/X.h:297
PlaceOnBottom = 1 	# /usr/include/X11/X.h:298
FamilyInternet = 0 	# /usr/include/X11/X.h:302
FamilyDECnet = 1 	# /usr/include/X11/X.h:303
FamilyChaos = 2 	# /usr/include/X11/X.h:304
FamilyInternet6 = 6 	# /usr/include/X11/X.h:305
FamilyServerInterpreted = 5 	# /usr/include/X11/X.h:308
PropertyNewValue = 0 	# /usr/include/X11/X.h:312
PropertyDelete = 1 	# /usr/include/X11/X.h:313
ColormapUninstalled = 0 	# /usr/include/X11/X.h:317
ColormapInstalled = 1 	# /usr/include/X11/X.h:318
GrabModeSync = 0 	# /usr/include/X11/X.h:322
GrabModeAsync = 1 	# /usr/include/X11/X.h:323
GrabSuccess = 0 	# /usr/include/X11/X.h:327
AlreadyGrabbed = 1 	# /usr/include/X11/X.h:328
GrabInvalidTime = 2 	# /usr/include/X11/X.h:329
GrabNotViewable = 3 	# /usr/include/X11/X.h:330
GrabFrozen = 4 	# /usr/include/X11/X.h:331
AsyncPointer = 0 	# /usr/include/X11/X.h:335
SyncPointer = 1 	# /usr/include/X11/X.h:336
ReplayPointer = 2 	# /usr/include/X11/X.h:337
AsyncKeyboard = 3 	# /usr/include/X11/X.h:338
SyncKeyboard = 4 	# /usr/include/X11/X.h:339
ReplayKeyboard = 5 	# /usr/include/X11/X.h:340
AsyncBoth = 6 	# /usr/include/X11/X.h:341
SyncBoth = 7 	# /usr/include/X11/X.h:342
RevertToParent = 2 	# /usr/include/X11/X.h:348
Success = 0 	# /usr/include/X11/X.h:354
BadRequest = 1 	# /usr/include/X11/X.h:355
BadValue = 2 	# /usr/include/X11/X.h:356
BadWindow = 3 	# /usr/include/X11/X.h:357
BadPixmap = 4 	# /usr/include/X11/X.h:358
BadAtom = 5 	# /usr/include/X11/X.h:359
BadCursor = 6 	# /usr/include/X11/X.h:360
BadFont = 7 	# /usr/include/X11/X.h:361
BadMatch = 8 	# /usr/include/X11/X.h:362
BadDrawable = 9 	# /usr/include/X11/X.h:363
BadAccess = 10 	# /usr/include/X11/X.h:364
BadAlloc = 11 	# /usr/include/X11/X.h:373
BadColor = 12 	# /usr/include/X11/X.h:374
BadGC = 13 	# /usr/include/X11/X.h:375
BadIDChoice = 14 	# /usr/include/X11/X.h:376
BadName = 15 	# /usr/include/X11/X.h:377
BadLength = 16 	# /usr/include/X11/X.h:378
BadImplementation = 17 	# /usr/include/X11/X.h:379
FirstExtensionError = 128 	# /usr/include/X11/X.h:381
LastExtensionError = 255 	# /usr/include/X11/X.h:382
InputOutput = 1 	# /usr/include/X11/X.h:391
InputOnly = 2 	# /usr/include/X11/X.h:392
CWBackPixmap = 1 	# /usr/include/X11/X.h:396
CWBackPixel = 2 	# /usr/include/X11/X.h:397
CWBorderPixmap = 4 	# /usr/include/X11/X.h:398
CWBorderPixel = 8 	# /usr/include/X11/X.h:399
CWBitGravity = 16 	# /usr/include/X11/X.h:400
CWWinGravity = 32 	# /usr/include/X11/X.h:401
CWBackingStore = 64 	# /usr/include/X11/X.h:402
CWBackingPlanes = 128 	# /usr/include/X11/X.h:403
CWBackingPixel = 256 	# /usr/include/X11/X.h:404
CWOverrideRedirect = 512 	# /usr/include/X11/X.h:405
CWSaveUnder = 1024 	# /usr/include/X11/X.h:406
CWEventMask = 2048 	# /usr/include/X11/X.h:407
CWDontPropagate = 4096 	# /usr/include/X11/X.h:408
CWColormap = 8192 	# /usr/include/X11/X.h:409
CWCursor = 16384 	# /usr/include/X11/X.h:410
CWX = 1 	# /usr/include/X11/X.h:414
CWY = 2 	# /usr/include/X11/X.h:415
CWWidth = 4 	# /usr/include/X11/X.h:416
CWHeight = 8 	# /usr/include/X11/X.h:417
CWBorderWidth = 16 	# /usr/include/X11/X.h:418
CWSibling = 32 	# /usr/include/X11/X.h:419
CWStackMode = 64 	# /usr/include/X11/X.h:420
ForgetGravity = 0 	# /usr/include/X11/X.h:425
NorthWestGravity = 1 	# /usr/include/X11/X.h:426
NorthGravity = 2 	# /usr/include/X11/X.h:427
NorthEastGravity = 3 	# /usr/include/X11/X.h:428
WestGravity = 4 	# /usr/include/X11/X.h:429
CenterGravity = 5 	# /usr/include/X11/X.h:430
EastGravity = 6 	# /usr/include/X11/X.h:431
SouthWestGravity = 7 	# /usr/include/X11/X.h:432
SouthGravity = 8 	# /usr/include/X11/X.h:433
SouthEastGravity = 9 	# /usr/include/X11/X.h:434
StaticGravity = 10 	# /usr/include/X11/X.h:435
UnmapGravity = 0 	# /usr/include/X11/X.h:439
NotUseful = 0 	# /usr/include/X11/X.h:443
WhenMapped = 1 	# /usr/include/X11/X.h:444
Always = 2 	# /usr/include/X11/X.h:445
IsUnmapped = 0 	# /usr/include/X11/X.h:449
IsUnviewable = 1 	# /usr/include/X11/X.h:450
IsViewable = 2 	# /usr/include/X11/X.h:451
SetModeInsert = 0 	# /usr/include/X11/X.h:455
SetModeDelete = 1 	# /usr/include/X11/X.h:456
DestroyAll = 0 	# /usr/include/X11/X.h:460
RetainPermanent = 1 	# /usr/include/X11/X.h:461
RetainTemporary = 2 	# /usr/include/X11/X.h:462
Above = 0 	# /usr/include/X11/X.h:466
Below = 1 	# /usr/include/X11/X.h:467
TopIf = 2 	# /usr/include/X11/X.h:468
BottomIf = 3 	# /usr/include/X11/X.h:469
Opposite = 4 	# /usr/include/X11/X.h:470
RaiseLowest = 0 	# /usr/include/X11/X.h:474
LowerHighest = 1 	# /usr/include/X11/X.h:475
PropModeReplace = 0 	# /usr/include/X11/X.h:479
PropModePrepend = 1 	# /usr/include/X11/X.h:480
PropModeAppend = 2 	# /usr/include/X11/X.h:481
GXclear = 0 	# /usr/include/X11/X.h:489
GXand = 1 	# /usr/include/X11/X.h:490
GXandReverse = 2 	# /usr/include/X11/X.h:491
GXcopy = 3 	# /usr/include/X11/X.h:492
GXandInverted = 4 	# /usr/include/X11/X.h:493
GXnoop = 5 	# /usr/include/X11/X.h:494
GXxor = 6 	# /usr/include/X11/X.h:495
GXor = 7 	# /usr/include/X11/X.h:496
GXnor = 8 	# /usr/include/X11/X.h:497
GXequiv = 9 	# /usr/include/X11/X.h:498
GXinvert = 10 	# /usr/include/X11/X.h:499
GXorReverse = 11 	# /usr/include/X11/X.h:500
GXcopyInverted = 12 	# /usr/include/X11/X.h:501
GXorInverted = 13 	# /usr/include/X11/X.h:502
GXnand = 14 	# /usr/include/X11/X.h:503
GXset = 0.0 	# /usr/include/X11/X.h:504
LineSolid = 0 	# /usr/include/X11/X.h:508
LineOnOffDash = 1 	# /usr/include/X11/X.h:509
LineDoubleDash = 2 	# /usr/include/X11/X.h:510
CapNotLast = 0 	# /usr/include/X11/X.h:514
CapButt = 1 	# /usr/include/X11/X.h:515
CapRound = 2 	# /usr/include/X11/X.h:516
CapProjecting = 3 	# /usr/include/X11/X.h:517
JoinMiter = 0 	# /usr/include/X11/X.h:521
JoinRound = 1 	# /usr/include/X11/X.h:522
JoinBevel = 2 	# /usr/include/X11/X.h:523
FillSolid = 0 	# /usr/include/X11/X.h:527
FillTiled = 1 	# /usr/include/X11/X.h:528
FillStippled = 2 	# /usr/include/X11/X.h:529
FillOpaqueStippled = 3 	# /usr/include/X11/X.h:530
EvenOddRule = 0 	# /usr/include/X11/X.h:534
WindingRule = 1 	# /usr/include/X11/X.h:535
ClipByChildren = 0 	# /usr/include/X11/X.h:539
IncludeInferiors = 1 	# /usr/include/X11/X.h:540
Unsorted = 0 	# /usr/include/X11/X.h:544
YSorted = 1 	# /usr/include/X11/X.h:545
YXSorted = 2 	# /usr/include/X11/X.h:546
YXBanded = 3 	# /usr/include/X11/X.h:547
CoordModeOrigin = 0 	# /usr/include/X11/X.h:551
CoordModePrevious = 1 	# /usr/include/X11/X.h:552
Complex = 0 	# /usr/include/X11/X.h:556
Nonconvex = 1 	# /usr/include/X11/X.h:557
Convex = 2 	# /usr/include/X11/X.h:558
ArcChord = 0 	# /usr/include/X11/X.h:562
ArcPieSlice = 1 	# /usr/include/X11/X.h:563
GCFunction = 1 	# /usr/include/X11/X.h:568
GCPlaneMask = 2 	# /usr/include/X11/X.h:569
GCForeground = 4 	# /usr/include/X11/X.h:570
GCBackground = 8 	# /usr/include/X11/X.h:571
GCLineWidth = 16 	# /usr/include/X11/X.h:572
GCLineStyle = 32 	# /usr/include/X11/X.h:573
GCCapStyle = 64 	# /usr/include/X11/X.h:574
GCJoinStyle = 128 	# /usr/include/X11/X.h:575
GCFillStyle = 256 	# /usr/include/X11/X.h:576
GCFillRule = 512 	# /usr/include/X11/X.h:577
GCTile = 1024 	# /usr/include/X11/X.h:578
GCStipple = 2048 	# /usr/include/X11/X.h:579
GCTileStipXOrigin = 4096 	# /usr/include/X11/X.h:580
GCTileStipYOrigin = 8192 	# /usr/include/X11/X.h:581
GCFont = 16384 	# /usr/include/X11/X.h:582
GCSubwindowMode = 32768 	# /usr/include/X11/X.h:583
GCGraphicsExposures = 65536 	# /usr/include/X11/X.h:584
GCClipXOrigin = 131072 	# /usr/include/X11/X.h:585
GCClipYOrigin = 262144 	# /usr/include/X11/X.h:586
GCClipMask = 524288 	# /usr/include/X11/X.h:587
GCDashOffset = 1048576 	# /usr/include/X11/X.h:588
GCDashList = 2097152 	# /usr/include/X11/X.h:589
GCArcMode = 4194304 	# /usr/include/X11/X.h:590
GCLastBit = 22 	# /usr/include/X11/X.h:592
FontLeftToRight = 0 	# /usr/include/X11/X.h:599
FontRightToLeft = 1 	# /usr/include/X11/X.h:600
FontChange = 255 	# /usr/include/X11/X.h:602
XYBitmap = 0 	# /usr/include/X11/X.h:610
XYPixmap = 1 	# /usr/include/X11/X.h:611
ZPixmap = 2 	# /usr/include/X11/X.h:612
AllocNone = 0 	# /usr/include/X11/X.h:620
AllocAll = 1 	# /usr/include/X11/X.h:621
DoRed = 1 	# /usr/include/X11/X.h:626
DoGreen = 2 	# /usr/include/X11/X.h:627
DoBlue = 4 	# /usr/include/X11/X.h:628
CursorShape = 0 	# /usr/include/X11/X.h:636
TileShape = 1 	# /usr/include/X11/X.h:637
StippleShape = 2 	# /usr/include/X11/X.h:638
AutoRepeatModeOff = 0 	# /usr/include/X11/X.h:644
AutoRepeatModeOn = 1 	# /usr/include/X11/X.h:645
AutoRepeatModeDefault = 2 	# /usr/include/X11/X.h:646
LedModeOff = 0 	# /usr/include/X11/X.h:648
LedModeOn = 1 	# /usr/include/X11/X.h:649
KBKeyClickPercent = 1 	# /usr/include/X11/X.h:653
KBBellPercent = 2 	# /usr/include/X11/X.h:654
KBBellPitch = 4 	# /usr/include/X11/X.h:655
KBBellDuration = 8 	# /usr/include/X11/X.h:656
KBLed = 16 	# /usr/include/X11/X.h:657
KBLedMode = 32 	# /usr/include/X11/X.h:658
KBKey = 64 	# /usr/include/X11/X.h:659
KBAutoRepeatMode = 128 	# /usr/include/X11/X.h:660
MappingSuccess = 0 	# /usr/include/X11/X.h:662
MappingBusy = 1 	# /usr/include/X11/X.h:663
MappingFailed = 2 	# /usr/include/X11/X.h:664
MappingModifier = 0 	# /usr/include/X11/X.h:666
MappingKeyboard = 1 	# /usr/include/X11/X.h:667
MappingPointer = 2 	# /usr/include/X11/X.h:668
DontPreferBlanking = 0 	# /usr/include/X11/X.h:674
PreferBlanking = 1 	# /usr/include/X11/X.h:675
DefaultBlanking = 2 	# /usr/include/X11/X.h:676
DisableScreenSaver = 0 	# /usr/include/X11/X.h:678
DisableScreenInterval = 0 	# /usr/include/X11/X.h:679
DontAllowExposures = 0 	# /usr/include/X11/X.h:681
AllowExposures = 1 	# /usr/include/X11/X.h:682
DefaultExposures = 2 	# /usr/include/X11/X.h:683
ScreenSaverReset = 0 	# /usr/include/X11/X.h:687
ScreenSaverActive = 1 	# /usr/include/X11/X.h:688
HostInsert = 0 	# /usr/include/X11/X.h:696
HostDelete = 1 	# /usr/include/X11/X.h:697
EnableAccess = 1 	# /usr/include/X11/X.h:701
DisableAccess = 0 	# /usr/include/X11/X.h:702
StaticGray = 0 	# /usr/include/X11/X.h:708
GrayScale = 1 	# /usr/include/X11/X.h:709
StaticColor = 2 	# /usr/include/X11/X.h:710
PseudoColor = 3 	# /usr/include/X11/X.h:711
TrueColor = 4 	# /usr/include/X11/X.h:712
DirectColor = 5 	# /usr/include/X11/X.h:713
LSBFirst = 0 	# /usr/include/X11/X.h:718
MSBFirst = 1 	# /usr/include/X11/X.h:719
# /usr/include/X11/Xlib.h:94
_Xmblen = _lib._Xmblen
_Xmblen.restype = c_int
_Xmblen.argtypes = [c_char_p, c_int]

X_HAVE_UTF8_STRING = 1 	# /usr/include/X11/Xlib.h:106
Bool = 0 	# /usr/include/X11/Xlib.h:110
Status = 0 	# /usr/include/X11/Xlib.h:111
True = 1 	# /usr/include/X11/Xlib.h:112
False = 0 	# /usr/include/X11/Xlib.h:113
QueuedAlready = 0 	# /usr/include/X11/Xlib.h:115
QueuedAfterReading = 1 	# /usr/include/X11/Xlib.h:116
QueuedAfterFlush = 2 	# /usr/include/X11/Xlib.h:117
XIMPreeditArea = 1 	# /usr/include/X11/Xlib.h:1176
XIMPreeditCallbacks = 2 	# /usr/include/X11/Xlib.h:1177
XIMPreeditPosition = 4 	# /usr/include/X11/Xlib.h:1178
XIMPreeditNothing = 8 	# /usr/include/X11/Xlib.h:1179
XIMPreeditNone = 16 	# /usr/include/X11/Xlib.h:1180
XIMStatusArea = 256 	# /usr/include/X11/Xlib.h:1181
XIMStatusCallbacks = 512 	# /usr/include/X11/Xlib.h:1182
XIMStatusNothing = 1024 	# /usr/include/X11/Xlib.h:1183
XIMStatusNone = 2048 	# /usr/include/X11/Xlib.h:1184
XBufferOverflow = -1 	# /usr/include/X11/Xlib.h:1230
XLookupNone = 1 	# /usr/include/X11/Xlib.h:1231
XLookupChars = 2 	# /usr/include/X11/Xlib.h:1232
XLookupKeySym = 3 	# /usr/include/X11/Xlib.h:1233
XLookupBoth = 4 	# /usr/include/X11/Xlib.h:1234
XIMReverse = 1 	# /usr/include/X11/Xlib.h:1250
XIMUnderline = 2 	# /usr/include/X11/Xlib.h:1251
XIMHighlight = 4 	# /usr/include/X11/Xlib.h:1252
XIMPrimary = 32 	# /usr/include/X11/Xlib.h:1253
XIMSecondary = 64 	# /usr/include/X11/Xlib.h:1254
XIMTertiary = 128 	# /usr/include/X11/Xlib.h:1255
XIMVisibleToForward = 256 	# /usr/include/X11/Xlib.h:1256
XIMVisibleToBackword = 512 	# /usr/include/X11/Xlib.h:1257
XIMVisibleToCenter = 1024 	# /usr/include/X11/Xlib.h:1258
XIMPreeditUnKnown = 0 	# /usr/include/X11/Xlib.h:1272
XIMPreeditEnable = 1 	# /usr/include/X11/Xlib.h:1273
XIMPreeditDisable = 2 	# /usr/include/X11/Xlib.h:1274
XIMInitialState = 1 	# /usr/include/X11/Xlib.h:1282
XIMPreserveState = 2 	# /usr/include/X11/Xlib.h:1283
XIMStringConversionLeftEdge = 1 	# /usr/include/X11/Xlib.h:1287
XIMStringConversionRightEdge = 2 	# /usr/include/X11/Xlib.h:1288
XIMStringConversionTopEdge = 4 	# /usr/include/X11/Xlib.h:1289
XIMStringConversionBottomEdge = 8 	# /usr/include/X11/Xlib.h:1290
XIMStringConversionConcealed = 16 	# /usr/include/X11/Xlib.h:1291
XIMStringConversionWrapped = 32 	# /usr/include/X11/Xlib.h:1292
XIMStringConversionBuffer = 1 	# /usr/include/X11/Xlib.h:1308
XIMStringConversionLine = 2 	# /usr/include/X11/Xlib.h:1309
XIMStringConversionWord = 3 	# /usr/include/X11/Xlib.h:1310
XIMStringConversionChar = 4 	# /usr/include/X11/Xlib.h:1311
XIMStringConversionSubstitution = 1 	# /usr/include/X11/Xlib.h:1315
XIMStringConversionRetrieval = 2 	# /usr/include/X11/Xlib.h:1316
XIMHotKeyStateON = 1 	# /usr/include/X11/Xlib.h:1381
XIMHotKeyStateOFF = 2 	# /usr/include/X11/Xlib.h:1382
# /usr/include/X11/Xlib.h:1397
XLoadQueryFont = _lib.XLoadQueryFont
XLoadQueryFont.restype = POINTER(XFontStruct)
XLoadQueryFont.argtypes = [POINTER(Display), c_char_p]

# /usr/include/X11/Xlib.h:1402
XQueryFont = _lib.XQueryFont
XQueryFont.restype = POINTER(XFontStruct)
XQueryFont.argtypes = [POINTER(Display), XID]

# /usr/include/X11/Xlib.h:1408
XGetMotionEvents = _lib.XGetMotionEvents
XGetMotionEvents.restype = POINTER(XTimeCoord)
XGetMotionEvents.argtypes = [POINTER(Display), Window, Time, Time, POINTER(c_int)]

# /usr/include/X11/Xlib.h:1416
XDeleteModifiermapEntry = _lib.XDeleteModifiermapEntry
XDeleteModifiermapEntry.restype = POINTER(XModifierKeymap)
XDeleteModifiermapEntry.argtypes = [POINTER(XModifierKeymap), KeyCode, c_int]

# /usr/include/X11/Xlib.h:1426
XGetModifierMapping = _lib.XGetModifierMapping
XGetModifierMapping.restype = POINTER(XModifierKeymap)
XGetModifierMapping.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1430
XInsertModifiermapEntry = _lib.XInsertModifiermapEntry
XInsertModifiermapEntry.restype = POINTER(XModifierKeymap)
XInsertModifiermapEntry.argtypes = [POINTER(XModifierKeymap), KeyCode, c_int]

# /usr/include/X11/Xlib.h:1440
XNewModifiermap = _lib.XNewModifiermap
XNewModifiermap.restype = POINTER(XModifierKeymap)
XNewModifiermap.argtypes = [c_int]

# /usr/include/X11/Xlib.h:1444
XCreateImage = _lib.XCreateImage
XCreateImage.restype = POINTER(XImage)
XCreateImage.argtypes = [POINTER(Display), POINTER(Visual), c_uint, c_int, c_int, c_char_p, c_uint, c_uint, c_int, c_int]

# /usr/include/X11/Xlib.h:1456
XInitImage = _lib.XInitImage
XInitImage.restype = c_int
XInitImage.argtypes = [POINTER(XImage)]

# /usr/include/X11/Xlib.h:1459
XGetImage = _lib.XGetImage
XGetImage.restype = POINTER(XImage)
XGetImage.argtypes = [POINTER(Display), Drawable, c_int, c_int, c_uint, c_uint, c_ulong, c_int]

# /usr/include/X11/Xlib.h:1469
XGetSubImage = _lib.XGetSubImage
XGetSubImage.restype = POINTER(XImage)
XGetSubImage.argtypes = [POINTER(Display), Drawable, c_int, c_int, c_uint, c_uint, c_ulong, c_int, POINTER(XImage), c_int, c_int]

# /usr/include/X11/Xlib.h:1486
XOpenDisplay = _lib.XOpenDisplay
XOpenDisplay.restype = POINTER(Display)
XOpenDisplay.argtypes = [c_char_p]

# /usr/include/X11/Xlib.h:1490
XrmInitialize = _lib.XrmInitialize
XrmInitialize.restype = None
XrmInitialize.argtypes = []

# /usr/include/X11/Xlib.h:1494
XFetchBytes = _lib.XFetchBytes
XFetchBytes.restype = c_char_p
XFetchBytes.argtypes = [POINTER(Display), POINTER(c_int)]

# /usr/include/X11/Xlib.h:1498
XFetchBuffer = _lib.XFetchBuffer
XFetchBuffer.restype = c_char_p
XFetchBuffer.argtypes = [POINTER(Display), POINTER(c_int), c_int]

# /usr/include/X11/Xlib.h:1503
XGetAtomName = _lib.XGetAtomName
XGetAtomName.restype = c_char_p
XGetAtomName.argtypes = [POINTER(Display), Atom]

# /usr/include/X11/Xlib.h:1507
XGetAtomNames = _lib.XGetAtomNames
XGetAtomNames.restype = c_int
XGetAtomNames.argtypes = [POINTER(Display), POINTER(Atom), c_int, POINTER(c_char_p)]

# /usr/include/X11/Xlib.h:1513
XGetDefault = _lib.XGetDefault
XGetDefault.restype = c_char_p
XGetDefault.argtypes = [POINTER(Display), c_char_p, c_char_p]

# /usr/include/X11/Xlib.h:1518
XDisplayName = _lib.XDisplayName
XDisplayName.restype = c_char_p
XDisplayName.argtypes = [c_char_p]

# /usr/include/X11/Xlib.h:1521
XKeysymToString = _lib.XKeysymToString
XKeysymToString.restype = c_char_p
XKeysymToString.argtypes = [KeySym]

# /usr/include/X11/Xlib.h:1525
XSynchronize = _lib.XSynchronize
XSynchronize.restype = POINTER(CFUNCTYPE(c_int, POINTER(Display)))
XSynchronize.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:1531
XSetAfterFunction = _lib.XSetAfterFunction
XSetAfterFunction.restype = POINTER(CFUNCTYPE(c_int, POINTER(Display)))
XSetAfterFunction.argtypes = [POINTER(Display), CFUNCTYPE(c_int, POINTER(Display))]

# /usr/include/X11/Xlib.h:1539
XInternAtom = _lib.XInternAtom
XInternAtom.restype = Atom
XInternAtom.argtypes = [POINTER(Display), c_char_p, c_int]

# /usr/include/X11/Xlib.h:1544
XInternAtoms = _lib.XInternAtoms
XInternAtoms.restype = c_int
XInternAtoms.argtypes = [POINTER(Display), POINTER(c_char_p), c_int, c_int, POINTER(Atom)]

# /usr/include/X11/Xlib.h:1551
XCopyColormapAndFree = _lib.XCopyColormapAndFree
XCopyColormapAndFree.restype = Colormap
XCopyColormapAndFree.argtypes = [POINTER(Display), Colormap]

# /usr/include/X11/Xlib.h:1555
XCreateColormap = _lib.XCreateColormap
XCreateColormap.restype = Colormap
XCreateColormap.argtypes = [POINTER(Display), Window, POINTER(Visual), c_int]

# /usr/include/X11/Xlib.h:1561
XCreatePixmapCursor = _lib.XCreatePixmapCursor
XCreatePixmapCursor.restype = Cursor
XCreatePixmapCursor.argtypes = [POINTER(Display), Pixmap, Pixmap, POINTER(XColor), POINTER(XColor), c_uint, c_uint]

# /usr/include/X11/Xlib.h:1570
XCreateGlyphCursor = _lib.XCreateGlyphCursor
XCreateGlyphCursor.restype = Cursor
XCreateGlyphCursor.argtypes = [POINTER(Display), Font, Font, c_uint, c_uint, POINTER(XColor), POINTER(XColor)]

# /usr/include/X11/Xlib.h:1579
XCreateFontCursor = _lib.XCreateFontCursor
XCreateFontCursor.restype = Cursor
XCreateFontCursor.argtypes = [POINTER(Display), c_uint]

# /usr/include/X11/Xlib.h:1583
XLoadFont = _lib.XLoadFont
XLoadFont.restype = Font
XLoadFont.argtypes = [POINTER(Display), c_char_p]

# /usr/include/X11/Xlib.h:1587
XCreateGC = _lib.XCreateGC
XCreateGC.restype = GC
XCreateGC.argtypes = [POINTER(Display), Drawable, c_ulong, POINTER(XGCValues)]

# /usr/include/X11/Xlib.h:1593
XGContextFromGC = _lib.XGContextFromGC
XGContextFromGC.restype = GContext
XGContextFromGC.argtypes = [GC]

# /usr/include/X11/Xlib.h:1596
XFlushGC = _lib.XFlushGC
XFlushGC.restype = None
XFlushGC.argtypes = [POINTER(Display), GC]

# /usr/include/X11/Xlib.h:1600
XCreatePixmap = _lib.XCreatePixmap
XCreatePixmap.restype = Pixmap
XCreatePixmap.argtypes = [POINTER(Display), Drawable, c_uint, c_uint, c_uint]

# /usr/include/X11/Xlib.h:1607
XCreateBitmapFromData = _lib.XCreateBitmapFromData
XCreateBitmapFromData.restype = Pixmap
XCreateBitmapFromData.argtypes = [POINTER(Display), Drawable, c_char_p, c_uint, c_uint]

# /usr/include/X11/Xlib.h:1614
XCreatePixmapFromBitmapData = _lib.XCreatePixmapFromBitmapData
XCreatePixmapFromBitmapData.restype = Pixmap
XCreatePixmapFromBitmapData.argtypes = [POINTER(Display), Drawable, c_char_p, c_uint, c_uint, c_ulong, c_ulong, c_uint]

# /usr/include/X11/Xlib.h:1624
XCreateSimpleWindow = _lib.XCreateSimpleWindow
XCreateSimpleWindow.restype = Window
XCreateSimpleWindow.argtypes = [POINTER(Display), Window, c_int, c_int, c_uint, c_uint, c_uint, c_ulong, c_ulong]

# /usr/include/X11/Xlib.h:1635
XGetSelectionOwner = _lib.XGetSelectionOwner
XGetSelectionOwner.restype = Window
XGetSelectionOwner.argtypes = [POINTER(Display), Atom]

# /usr/include/X11/Xlib.h:1639
XCreateWindow = _lib.XCreateWindow
XCreateWindow.restype = Window
XCreateWindow.argtypes = [POINTER(Display), Window, c_int, c_int, c_uint, c_uint, c_uint, c_int, c_uint, POINTER(Visual), c_ulong, POINTER(XSetWindowAttributes)]

# /usr/include/X11/Xlib.h:1653
XListInstalledColormaps = _lib.XListInstalledColormaps
XListInstalledColormaps.restype = POINTER(Colormap)
XListInstalledColormaps.argtypes = [POINTER(Display), Window, POINTER(c_int)]

# /usr/include/X11/Xlib.h:1658
XListFonts = _lib.XListFonts
XListFonts.restype = POINTER(c_char_p)
XListFonts.argtypes = [POINTER(Display), c_char_p, c_int, POINTER(c_int)]

# /usr/include/X11/Xlib.h:1664
XListFontsWithInfo = _lib.XListFontsWithInfo
XListFontsWithInfo.restype = POINTER(c_char_p)
XListFontsWithInfo.argtypes = [POINTER(Display), c_char_p, c_int, POINTER(c_int), POINTER(POINTER(XFontStruct))]

# /usr/include/X11/Xlib.h:1671
XGetFontPath = _lib.XGetFontPath
XGetFontPath.restype = POINTER(c_char_p)
XGetFontPath.argtypes = [POINTER(Display), POINTER(c_int)]

# /usr/include/X11/Xlib.h:1675
XListExtensions = _lib.XListExtensions
XListExtensions.restype = POINTER(c_char_p)
XListExtensions.argtypes = [POINTER(Display), POINTER(c_int)]

# /usr/include/X11/Xlib.h:1679
XListProperties = _lib.XListProperties
XListProperties.restype = POINTER(Atom)
XListProperties.argtypes = [POINTER(Display), Window, POINTER(c_int)]

# /usr/include/X11/Xlib.h:1684
XListHosts = _lib.XListHosts
XListHosts.restype = POINTER(XHostAddress)
XListHosts.argtypes = [POINTER(Display), POINTER(c_int), POINTER(c_int)]

# /usr/include/X11/Xlib.h:1689
XKeycodeToKeysym = _lib.XKeycodeToKeysym
XKeycodeToKeysym.restype = KeySym
XKeycodeToKeysym.argtypes = [POINTER(Display), KeyCode, c_int]

# /usr/include/X11/Xlib.h:1698
XLookupKeysym = _lib.XLookupKeysym
XLookupKeysym.restype = KeySym
XLookupKeysym.argtypes = [POINTER(XKeyEvent), c_int]

# /usr/include/X11/Xlib.h:1702
XGetKeyboardMapping = _lib.XGetKeyboardMapping
XGetKeyboardMapping.restype = POINTER(KeySym)
XGetKeyboardMapping.argtypes = [POINTER(Display), KeyCode, c_int, POINTER(c_int)]

# /usr/include/X11/Xlib.h:1712
XStringToKeysym = _lib.XStringToKeysym
XStringToKeysym.restype = KeySym
XStringToKeysym.argtypes = [c_char_p]

# /usr/include/X11/Xlib.h:1715
XMaxRequestSize = _lib.XMaxRequestSize
XMaxRequestSize.restype = c_long
XMaxRequestSize.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1718
XExtendedMaxRequestSize = _lib.XExtendedMaxRequestSize
XExtendedMaxRequestSize.restype = c_long
XExtendedMaxRequestSize.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1721
XResourceManagerString = _lib.XResourceManagerString
XResourceManagerString.restype = c_char_p
XResourceManagerString.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1724
XScreenResourceString = _lib.XScreenResourceString
XScreenResourceString.restype = c_char_p
XScreenResourceString.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:1727
XDisplayMotionBufferSize = _lib.XDisplayMotionBufferSize
XDisplayMotionBufferSize.restype = c_ulong
XDisplayMotionBufferSize.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1730
XVisualIDFromVisual = _lib.XVisualIDFromVisual
XVisualIDFromVisual.restype = VisualID
XVisualIDFromVisual.argtypes = [POINTER(Visual)]

# /usr/include/X11/Xlib.h:1736
XInitThreads = _lib.XInitThreads
XInitThreads.restype = c_int
XInitThreads.argtypes = []

# /usr/include/X11/Xlib.h:1740
XLockDisplay = _lib.XLockDisplay
XLockDisplay.restype = None
XLockDisplay.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1744
XUnlockDisplay = _lib.XUnlockDisplay
XUnlockDisplay.restype = None
XUnlockDisplay.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1750
XInitExtension = _lib.XInitExtension
XInitExtension.restype = POINTER(XExtCodes)
XInitExtension.argtypes = [POINTER(Display), c_char_p]

# /usr/include/X11/Xlib.h:1755
XAddExtension = _lib.XAddExtension
XAddExtension.restype = POINTER(XExtCodes)
XAddExtension.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1758
XFindOnExtensionList = _lib.XFindOnExtensionList
XFindOnExtensionList.restype = POINTER(XExtData)
XFindOnExtensionList.argtypes = [POINTER(POINTER(XExtData)), c_int]

# /usr/include/X11/Xlib.h:1762
XEHeadOfExtensionList = _lib.XEHeadOfExtensionList
XEHeadOfExtensionList.restype = POINTER(POINTER(XExtData))
XEHeadOfExtensionList.argtypes = [XEDataObject]

# /usr/include/X11/Xlib.h:1767
XRootWindow = _lib.XRootWindow
XRootWindow.restype = Window
XRootWindow.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:1771
XDefaultRootWindow = _lib.XDefaultRootWindow
XDefaultRootWindow.restype = Window
XDefaultRootWindow.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1774
XRootWindowOfScreen = _lib.XRootWindowOfScreen
XRootWindowOfScreen.restype = Window
XRootWindowOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:1777
XDefaultVisual = _lib.XDefaultVisual
XDefaultVisual.restype = POINTER(Visual)
XDefaultVisual.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:1781
XDefaultVisualOfScreen = _lib.XDefaultVisualOfScreen
XDefaultVisualOfScreen.restype = POINTER(Visual)
XDefaultVisualOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:1784
XDefaultGC = _lib.XDefaultGC
XDefaultGC.restype = GC
XDefaultGC.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:1788
XDefaultGCOfScreen = _lib.XDefaultGCOfScreen
XDefaultGCOfScreen.restype = GC
XDefaultGCOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:1791
XBlackPixel = _lib.XBlackPixel
XBlackPixel.restype = c_ulong
XBlackPixel.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:1795
XWhitePixel = _lib.XWhitePixel
XWhitePixel.restype = c_ulong
XWhitePixel.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:1799
XAllPlanes = _lib.XAllPlanes
XAllPlanes.restype = c_ulong
XAllPlanes.argtypes = []

# /usr/include/X11/Xlib.h:1802
XBlackPixelOfScreen = _lib.XBlackPixelOfScreen
XBlackPixelOfScreen.restype = c_ulong
XBlackPixelOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:1805
XWhitePixelOfScreen = _lib.XWhitePixelOfScreen
XWhitePixelOfScreen.restype = c_ulong
XWhitePixelOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:1808
XNextRequest = _lib.XNextRequest
XNextRequest.restype = c_ulong
XNextRequest.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1811
XLastKnownRequestProcessed = _lib.XLastKnownRequestProcessed
XLastKnownRequestProcessed.restype = c_ulong
XLastKnownRequestProcessed.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1814
XServerVendor = _lib.XServerVendor
XServerVendor.restype = c_char_p
XServerVendor.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1817
XDisplayString = _lib.XDisplayString
XDisplayString.restype = c_char_p
XDisplayString.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1820
XDefaultColormap = _lib.XDefaultColormap
XDefaultColormap.restype = Colormap
XDefaultColormap.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:1824
XDefaultColormapOfScreen = _lib.XDefaultColormapOfScreen
XDefaultColormapOfScreen.restype = Colormap
XDefaultColormapOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:1827
XDisplayOfScreen = _lib.XDisplayOfScreen
XDisplayOfScreen.restype = POINTER(Display)
XDisplayOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:1830
XScreenOfDisplay = _lib.XScreenOfDisplay
XScreenOfDisplay.restype = POINTER(Screen)
XScreenOfDisplay.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:1834
XDefaultScreenOfDisplay = _lib.XDefaultScreenOfDisplay
XDefaultScreenOfDisplay.restype = POINTER(Screen)
XDefaultScreenOfDisplay.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1837
XEventMaskOfScreen = _lib.XEventMaskOfScreen
XEventMaskOfScreen.restype = c_long
XEventMaskOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:1841
XScreenNumberOfScreen = _lib.XScreenNumberOfScreen
XScreenNumberOfScreen.restype = c_int
XScreenNumberOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:1850
XSetErrorHandler = _lib.XSetErrorHandler
XSetErrorHandler.restype = XErrorHandler
XSetErrorHandler.argtypes = [XErrorHandler]

# /usr/include/X11/Xlib.h:1859
XSetIOErrorHandler = _lib.XSetIOErrorHandler
XSetIOErrorHandler.restype = XIOErrorHandler
XSetIOErrorHandler.argtypes = [XIOErrorHandler]

# /usr/include/X11/Xlib.h:1864
XListPixmapFormats = _lib.XListPixmapFormats
XListPixmapFormats.restype = POINTER(XPixmapFormatValues)
XListPixmapFormats.argtypes = [POINTER(Display), POINTER(c_int)]

# /usr/include/X11/Xlib.h:1868
XListDepths = _lib.XListDepths
XListDepths.restype = POINTER(c_int)
XListDepths.argtypes = [POINTER(Display), c_int, POINTER(c_int)]

# /usr/include/X11/Xlib.h:1876
XReconfigureWMWindow = _lib.XReconfigureWMWindow
XReconfigureWMWindow.restype = c_int
XReconfigureWMWindow.argtypes = [POINTER(Display), Window, c_int, c_uint, POINTER(XWindowChanges)]

# /usr/include/X11/Xlib.h:1884
XGetWMProtocols = _lib.XGetWMProtocols
XGetWMProtocols.restype = c_int
XGetWMProtocols.argtypes = [POINTER(Display), Window, POINTER(POINTER(Atom)), POINTER(c_int)]

# /usr/include/X11/Xlib.h:1890
XSetWMProtocols = _lib.XSetWMProtocols
XSetWMProtocols.restype = c_int
XSetWMProtocols.argtypes = [POINTER(Display), Window, POINTER(Atom), c_int]

# /usr/include/X11/Xlib.h:1896
XIconifyWindow = _lib.XIconifyWindow
XIconifyWindow.restype = c_int
XIconifyWindow.argtypes = [POINTER(Display), Window, c_int]

# /usr/include/X11/Xlib.h:1901
XWithdrawWindow = _lib.XWithdrawWindow
XWithdrawWindow.restype = c_int
XWithdrawWindow.argtypes = [POINTER(Display), Window, c_int]

# /usr/include/X11/Xlib.h:1906
XGetCommand = _lib.XGetCommand
XGetCommand.restype = c_int
XGetCommand.argtypes = [POINTER(Display), Window, POINTER(POINTER(c_char_p)), POINTER(c_int)]

# /usr/include/X11/Xlib.h:1912
XGetWMColormapWindows = _lib.XGetWMColormapWindows
XGetWMColormapWindows.restype = c_int
XGetWMColormapWindows.argtypes = [POINTER(Display), Window, POINTER(POINTER(Window)), POINTER(c_int)]

# /usr/include/X11/Xlib.h:1918
XSetWMColormapWindows = _lib.XSetWMColormapWindows
XSetWMColormapWindows.restype = c_int
XSetWMColormapWindows.argtypes = [POINTER(Display), Window, POINTER(Window), c_int]

# /usr/include/X11/Xlib.h:1924
XFreeStringList = _lib.XFreeStringList
XFreeStringList.restype = None
XFreeStringList.argtypes = [POINTER(c_char_p)]

# /usr/include/X11/Xlib.h:1927
XSetTransientForHint = _lib.XSetTransientForHint
XSetTransientForHint.restype = c_int
XSetTransientForHint.argtypes = [POINTER(Display), Window, Window]

# /usr/include/X11/Xlib.h:1935
XActivateScreenSaver = _lib.XActivateScreenSaver
XActivateScreenSaver.restype = c_int
XActivateScreenSaver.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:1939
XAddHost = _lib.XAddHost
XAddHost.restype = c_int
XAddHost.argtypes = [POINTER(Display), POINTER(XHostAddress)]

# /usr/include/X11/Xlib.h:1944
XAddHosts = _lib.XAddHosts
XAddHosts.restype = c_int
XAddHosts.argtypes = [POINTER(Display), POINTER(XHostAddress), c_int]

# /usr/include/X11/Xlib.h:1950
XAddToExtensionList = _lib.XAddToExtensionList
XAddToExtensionList.restype = c_int
XAddToExtensionList.argtypes = [POINTER(POINTER(struct__XExtData)), POINTER(XExtData)]

# /usr/include/X11/Xlib.h:1955
XAddToSaveSet = _lib.XAddToSaveSet
XAddToSaveSet.restype = c_int
XAddToSaveSet.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:1960
XAllocColor = _lib.XAllocColor
XAllocColor.restype = c_int
XAllocColor.argtypes = [POINTER(Display), Colormap, POINTER(XColor)]

# /usr/include/X11/Xlib.h:1966
XAllocColorCells = _lib.XAllocColorCells
XAllocColorCells.restype = c_int
XAllocColorCells.argtypes = [POINTER(Display), Colormap, c_int, POINTER(c_ulong), c_uint, POINTER(c_ulong), c_uint]

# /usr/include/X11/Xlib.h:1976
XAllocColorPlanes = _lib.XAllocColorPlanes
XAllocColorPlanes.restype = c_int
XAllocColorPlanes.argtypes = [POINTER(Display), Colormap, c_int, POINTER(c_ulong), c_int, c_int, c_int, c_int, POINTER(c_ulong), POINTER(c_ulong), POINTER(c_ulong)]

# /usr/include/X11/Xlib.h:1990
XAllocNamedColor = _lib.XAllocNamedColor
XAllocNamedColor.restype = c_int
XAllocNamedColor.argtypes = [POINTER(Display), Colormap, c_char_p, POINTER(XColor), POINTER(XColor)]

# /usr/include/X11/Xlib.h:1998
XAllowEvents = _lib.XAllowEvents
XAllowEvents.restype = c_int
XAllowEvents.argtypes = [POINTER(Display), c_int, Time]

# /usr/include/X11/Xlib.h:2004
XAutoRepeatOff = _lib.XAutoRepeatOff
XAutoRepeatOff.restype = c_int
XAutoRepeatOff.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2008
XAutoRepeatOn = _lib.XAutoRepeatOn
XAutoRepeatOn.restype = c_int
XAutoRepeatOn.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2012
XBell = _lib.XBell
XBell.restype = c_int
XBell.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:2017
XBitmapBitOrder = _lib.XBitmapBitOrder
XBitmapBitOrder.restype = c_int
XBitmapBitOrder.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2021
XBitmapPad = _lib.XBitmapPad
XBitmapPad.restype = c_int
XBitmapPad.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2025
XBitmapUnit = _lib.XBitmapUnit
XBitmapUnit.restype = c_int
XBitmapUnit.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2029
XCellsOfScreen = _lib.XCellsOfScreen
XCellsOfScreen.restype = c_int
XCellsOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:2033
XChangeActivePointerGrab = _lib.XChangeActivePointerGrab
XChangeActivePointerGrab.restype = c_int
XChangeActivePointerGrab.argtypes = [POINTER(Display), c_uint, Cursor, Time]

# /usr/include/X11/Xlib.h:2040
XChangeGC = _lib.XChangeGC
XChangeGC.restype = c_int
XChangeGC.argtypes = [POINTER(Display), GC, c_ulong, POINTER(XGCValues)]

# /usr/include/X11/Xlib.h:2047
XChangeKeyboardControl = _lib.XChangeKeyboardControl
XChangeKeyboardControl.restype = c_int
XChangeKeyboardControl.argtypes = [POINTER(Display), c_ulong, POINTER(XKeyboardControl)]

# /usr/include/X11/Xlib.h:2053
XChangeKeyboardMapping = _lib.XChangeKeyboardMapping
XChangeKeyboardMapping.restype = c_int
XChangeKeyboardMapping.argtypes = [POINTER(Display), c_int, c_int, POINTER(KeySym), c_int]

# /usr/include/X11/Xlib.h:2061
XChangePointerControl = _lib.XChangePointerControl
XChangePointerControl.restype = c_int
XChangePointerControl.argtypes = [POINTER(Display), c_int, c_int, c_int, c_int, c_int]

# /usr/include/X11/Xlib.h:2070
XChangeProperty = _lib.XChangeProperty
XChangeProperty.restype = c_int
XChangeProperty.argtypes = [POINTER(Display), Window, Atom, Atom, c_int, c_int, POINTER(c_ubyte), c_int]

# /usr/include/X11/Xlib.h:2081
XChangeSaveSet = _lib.XChangeSaveSet
XChangeSaveSet.restype = c_int
XChangeSaveSet.argtypes = [POINTER(Display), Window, c_int]

# /usr/include/X11/Xlib.h:2087
XChangeWindowAttributes = _lib.XChangeWindowAttributes
XChangeWindowAttributes.restype = c_int
XChangeWindowAttributes.argtypes = [POINTER(Display), Window, c_ulong, POINTER(XSetWindowAttributes)]

# /usr/include/X11/Xlib.h:2094
XCheckIfEvent = _lib.XCheckIfEvent
XCheckIfEvent.restype = c_int
XCheckIfEvent.argtypes = [POINTER(Display), POINTER(XEvent), CFUNCTYPE(c_int, POINTER(Display), POINTER(XEvent), XPointer), XPointer]

# /usr/include/X11/Xlib.h:2105
XCheckMaskEvent = _lib.XCheckMaskEvent
XCheckMaskEvent.restype = c_int
XCheckMaskEvent.argtypes = [POINTER(Display), c_long, POINTER(XEvent)]

# /usr/include/X11/Xlib.h:2111
XCheckTypedEvent = _lib.XCheckTypedEvent
XCheckTypedEvent.restype = c_int
XCheckTypedEvent.argtypes = [POINTER(Display), c_int, POINTER(XEvent)]

# /usr/include/X11/Xlib.h:2117
XCheckTypedWindowEvent = _lib.XCheckTypedWindowEvent
XCheckTypedWindowEvent.restype = c_int
XCheckTypedWindowEvent.argtypes = [POINTER(Display), Window, c_int, POINTER(XEvent)]

# /usr/include/X11/Xlib.h:2124
XCheckWindowEvent = _lib.XCheckWindowEvent
XCheckWindowEvent.restype = c_int
XCheckWindowEvent.argtypes = [POINTER(Display), Window, c_long, POINTER(XEvent)]

# /usr/include/X11/Xlib.h:2131
XCirculateSubwindows = _lib.XCirculateSubwindows
XCirculateSubwindows.restype = c_int
XCirculateSubwindows.argtypes = [POINTER(Display), Window, c_int]

# /usr/include/X11/Xlib.h:2137
XCirculateSubwindowsDown = _lib.XCirculateSubwindowsDown
XCirculateSubwindowsDown.restype = c_int
XCirculateSubwindowsDown.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:2142
XCirculateSubwindowsUp = _lib.XCirculateSubwindowsUp
XCirculateSubwindowsUp.restype = c_int
XCirculateSubwindowsUp.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:2147
XClearArea = _lib.XClearArea
XClearArea.restype = c_int
XClearArea.argtypes = [POINTER(Display), Window, c_int, c_int, c_uint, c_uint, c_int]

# /usr/include/X11/Xlib.h:2157
XClearWindow = _lib.XClearWindow
XClearWindow.restype = c_int
XClearWindow.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:2162
XCloseDisplay = _lib.XCloseDisplay
XCloseDisplay.restype = c_int
XCloseDisplay.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2166
XConfigureWindow = _lib.XConfigureWindow
XConfigureWindow.restype = c_int
XConfigureWindow.argtypes = [POINTER(Display), Window, c_uint, POINTER(XWindowChanges)]

# /usr/include/X11/Xlib.h:2173
XConnectionNumber = _lib.XConnectionNumber
XConnectionNumber.restype = c_int
XConnectionNumber.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2177
XConvertSelection = _lib.XConvertSelection
XConvertSelection.restype = c_int
XConvertSelection.argtypes = [POINTER(Display), Atom, Atom, Atom, Window, Time]

# /usr/include/X11/Xlib.h:2186
XCopyArea = _lib.XCopyArea
XCopyArea.restype = c_int
XCopyArea.argtypes = [POINTER(Display), Drawable, Drawable, GC, c_int, c_int, c_uint, c_uint, c_int, c_int]

# /usr/include/X11/Xlib.h:2199
XCopyGC = _lib.XCopyGC
XCopyGC.restype = c_int
XCopyGC.argtypes = [POINTER(Display), GC, c_ulong, GC]

# /usr/include/X11/Xlib.h:2206
XCopyPlane = _lib.XCopyPlane
XCopyPlane.restype = c_int
XCopyPlane.argtypes = [POINTER(Display), Drawable, Drawable, GC, c_int, c_int, c_uint, c_uint, c_int, c_int, c_ulong]

# /usr/include/X11/Xlib.h:2220
XDefaultDepth = _lib.XDefaultDepth
XDefaultDepth.restype = c_int
XDefaultDepth.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:2225
XDefaultDepthOfScreen = _lib.XDefaultDepthOfScreen
XDefaultDepthOfScreen.restype = c_int
XDefaultDepthOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:2229
XDefaultScreen = _lib.XDefaultScreen
XDefaultScreen.restype = c_int
XDefaultScreen.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2233
XDefineCursor = _lib.XDefineCursor
XDefineCursor.restype = c_int
XDefineCursor.argtypes = [POINTER(Display), Window, Cursor]

# /usr/include/X11/Xlib.h:2239
XDeleteProperty = _lib.XDeleteProperty
XDeleteProperty.restype = c_int
XDeleteProperty.argtypes = [POINTER(Display), Window, Atom]

# /usr/include/X11/Xlib.h:2245
XDestroyWindow = _lib.XDestroyWindow
XDestroyWindow.restype = c_int
XDestroyWindow.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:2250
XDestroySubwindows = _lib.XDestroySubwindows
XDestroySubwindows.restype = c_int
XDestroySubwindows.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:2255
XDoesBackingStore = _lib.XDoesBackingStore
XDoesBackingStore.restype = c_int
XDoesBackingStore.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:2259
XDoesSaveUnders = _lib.XDoesSaveUnders
XDoesSaveUnders.restype = c_int
XDoesSaveUnders.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:2263
XDisableAccessControl = _lib.XDisableAccessControl
XDisableAccessControl.restype = c_int
XDisableAccessControl.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2268
XDisplayCells = _lib.XDisplayCells
XDisplayCells.restype = c_int
XDisplayCells.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:2273
XDisplayHeight = _lib.XDisplayHeight
XDisplayHeight.restype = c_int
XDisplayHeight.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:2278
XDisplayHeightMM = _lib.XDisplayHeightMM
XDisplayHeightMM.restype = c_int
XDisplayHeightMM.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:2283
XDisplayKeycodes = _lib.XDisplayKeycodes
XDisplayKeycodes.restype = c_int
XDisplayKeycodes.argtypes = [POINTER(Display), POINTER(c_int), POINTER(c_int)]

# /usr/include/X11/Xlib.h:2289
XDisplayPlanes = _lib.XDisplayPlanes
XDisplayPlanes.restype = c_int
XDisplayPlanes.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:2294
XDisplayWidth = _lib.XDisplayWidth
XDisplayWidth.restype = c_int
XDisplayWidth.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:2299
XDisplayWidthMM = _lib.XDisplayWidthMM
XDisplayWidthMM.restype = c_int
XDisplayWidthMM.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:2304
XDrawArc = _lib.XDrawArc
XDrawArc.restype = c_int
XDrawArc.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, c_uint, c_uint, c_int, c_int]

# /usr/include/X11/Xlib.h:2316
XDrawArcs = _lib.XDrawArcs
XDrawArcs.restype = c_int
XDrawArcs.argtypes = [POINTER(Display), Drawable, GC, POINTER(XArc), c_int]

# /usr/include/X11/Xlib.h:2324
XDrawImageString = _lib.XDrawImageString
XDrawImageString.restype = c_int
XDrawImageString.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, c_char_p, c_int]

# /usr/include/X11/Xlib.h:2334
XDrawImageString16 = _lib.XDrawImageString16
XDrawImageString16.restype = c_int
XDrawImageString16.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, POINTER(XChar2b), c_int]

# /usr/include/X11/Xlib.h:2344
XDrawLine = _lib.XDrawLine
XDrawLine.restype = c_int
XDrawLine.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, c_int, c_int]

# /usr/include/X11/Xlib.h:2354
XDrawLines = _lib.XDrawLines
XDrawLines.restype = c_int
XDrawLines.argtypes = [POINTER(Display), Drawable, GC, POINTER(XPoint), c_int, c_int]

# /usr/include/X11/Xlib.h:2363
XDrawPoint = _lib.XDrawPoint
XDrawPoint.restype = c_int
XDrawPoint.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int]

# /usr/include/X11/Xlib.h:2371
XDrawPoints = _lib.XDrawPoints
XDrawPoints.restype = c_int
XDrawPoints.argtypes = [POINTER(Display), Drawable, GC, POINTER(XPoint), c_int, c_int]

# /usr/include/X11/Xlib.h:2380
XDrawRectangle = _lib.XDrawRectangle
XDrawRectangle.restype = c_int
XDrawRectangle.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, c_uint, c_uint]

# /usr/include/X11/Xlib.h:2390
XDrawRectangles = _lib.XDrawRectangles
XDrawRectangles.restype = c_int
XDrawRectangles.argtypes = [POINTER(Display), Drawable, GC, POINTER(XRectangle), c_int]

# /usr/include/X11/Xlib.h:2398
XDrawSegments = _lib.XDrawSegments
XDrawSegments.restype = c_int
XDrawSegments.argtypes = [POINTER(Display), Drawable, GC, POINTER(XSegment), c_int]

# /usr/include/X11/Xlib.h:2406
XDrawString = _lib.XDrawString
XDrawString.restype = c_int
XDrawString.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, c_char_p, c_int]

# /usr/include/X11/Xlib.h:2416
XDrawString16 = _lib.XDrawString16
XDrawString16.restype = c_int
XDrawString16.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, POINTER(XChar2b), c_int]

# /usr/include/X11/Xlib.h:2426
XDrawText = _lib.XDrawText
XDrawText.restype = c_int
XDrawText.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, POINTER(XTextItem), c_int]

# /usr/include/X11/Xlib.h:2436
XDrawText16 = _lib.XDrawText16
XDrawText16.restype = c_int
XDrawText16.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, POINTER(XTextItem16), c_int]

# /usr/include/X11/Xlib.h:2446
XEnableAccessControl = _lib.XEnableAccessControl
XEnableAccessControl.restype = c_int
XEnableAccessControl.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2450
XEventsQueued = _lib.XEventsQueued
XEventsQueued.restype = c_int
XEventsQueued.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:2455
XFetchName = _lib.XFetchName
XFetchName.restype = c_int
XFetchName.argtypes = [POINTER(Display), Window, POINTER(c_char_p)]

# /usr/include/X11/Xlib.h:2461
XFillArc = _lib.XFillArc
XFillArc.restype = c_int
XFillArc.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, c_uint, c_uint, c_int, c_int]

# /usr/include/X11/Xlib.h:2473
XFillArcs = _lib.XFillArcs
XFillArcs.restype = c_int
XFillArcs.argtypes = [POINTER(Display), Drawable, GC, POINTER(XArc), c_int]

# /usr/include/X11/Xlib.h:2481
XFillPolygon = _lib.XFillPolygon
XFillPolygon.restype = c_int
XFillPolygon.argtypes = [POINTER(Display), Drawable, GC, POINTER(XPoint), c_int, c_int, c_int]

# /usr/include/X11/Xlib.h:2491
XFillRectangle = _lib.XFillRectangle
XFillRectangle.restype = c_int
XFillRectangle.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, c_uint, c_uint]

# /usr/include/X11/Xlib.h:2501
XFillRectangles = _lib.XFillRectangles
XFillRectangles.restype = c_int
XFillRectangles.argtypes = [POINTER(Display), Drawable, GC, POINTER(XRectangle), c_int]

# /usr/include/X11/Xlib.h:2509
XFlush = _lib.XFlush
XFlush.restype = c_int
XFlush.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2513
XForceScreenSaver = _lib.XForceScreenSaver
XForceScreenSaver.restype = c_int
XForceScreenSaver.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:2518
XFree = _lib.XFree
XFree.restype = c_int
XFree.argtypes = [POINTER(None)]

# /usr/include/X11/Xlib.h:2522
XFreeColormap = _lib.XFreeColormap
XFreeColormap.restype = c_int
XFreeColormap.argtypes = [POINTER(Display), Colormap]

# /usr/include/X11/Xlib.h:2527
XFreeColors = _lib.XFreeColors
XFreeColors.restype = c_int
XFreeColors.argtypes = [POINTER(Display), Colormap, POINTER(c_ulong), c_int, c_ulong]

# /usr/include/X11/Xlib.h:2535
XFreeCursor = _lib.XFreeCursor
XFreeCursor.restype = c_int
XFreeCursor.argtypes = [POINTER(Display), Cursor]

# /usr/include/X11/Xlib.h:2540
XFreeExtensionList = _lib.XFreeExtensionList
XFreeExtensionList.restype = c_int
XFreeExtensionList.argtypes = [POINTER(c_char_p)]

# /usr/include/X11/Xlib.h:2544
XFreeFont = _lib.XFreeFont
XFreeFont.restype = c_int
XFreeFont.argtypes = [POINTER(Display), POINTER(XFontStruct)]

# /usr/include/X11/Xlib.h:2549
XFreeFontInfo = _lib.XFreeFontInfo
XFreeFontInfo.restype = c_int
XFreeFontInfo.argtypes = [POINTER(c_char_p), POINTER(XFontStruct), c_int]

# /usr/include/X11/Xlib.h:2555
XFreeFontNames = _lib.XFreeFontNames
XFreeFontNames.restype = c_int
XFreeFontNames.argtypes = [POINTER(c_char_p)]

# /usr/include/X11/Xlib.h:2559
XFreeFontPath = _lib.XFreeFontPath
XFreeFontPath.restype = c_int
XFreeFontPath.argtypes = [POINTER(c_char_p)]

# /usr/include/X11/Xlib.h:2563
XFreeGC = _lib.XFreeGC
XFreeGC.restype = c_int
XFreeGC.argtypes = [POINTER(Display), GC]

# /usr/include/X11/Xlib.h:2568
XFreeModifiermap = _lib.XFreeModifiermap
XFreeModifiermap.restype = c_int
XFreeModifiermap.argtypes = [POINTER(XModifierKeymap)]

# /usr/include/X11/Xlib.h:2572
XFreePixmap = _lib.XFreePixmap
XFreePixmap.restype = c_int
XFreePixmap.argtypes = [POINTER(Display), Pixmap]

# /usr/include/X11/Xlib.h:2577
XGeometry = _lib.XGeometry
XGeometry.restype = c_int
XGeometry.argtypes = [POINTER(Display), c_int, c_char_p, c_char_p, c_uint, c_uint, c_uint, c_int, c_int, POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int)]

# /usr/include/X11/Xlib.h:2593
XGetErrorDatabaseText = _lib.XGetErrorDatabaseText
XGetErrorDatabaseText.restype = c_int
XGetErrorDatabaseText.argtypes = [POINTER(Display), c_char_p, c_char_p, c_char_p, c_char_p, c_int]

# /usr/include/X11/Xlib.h:2602
XGetErrorText = _lib.XGetErrorText
XGetErrorText.restype = c_int
XGetErrorText.argtypes = [POINTER(Display), c_int, c_char_p, c_int]

# /usr/include/X11/Xlib.h:2609
XGetFontProperty = _lib.XGetFontProperty
XGetFontProperty.restype = c_int
XGetFontProperty.argtypes = [POINTER(XFontStruct), Atom, POINTER(c_ulong)]

# /usr/include/X11/Xlib.h:2615
XGetGCValues = _lib.XGetGCValues
XGetGCValues.restype = c_int
XGetGCValues.argtypes = [POINTER(Display), GC, c_ulong, POINTER(XGCValues)]

# /usr/include/X11/Xlib.h:2622
XGetGeometry = _lib.XGetGeometry
XGetGeometry.restype = c_int
XGetGeometry.argtypes = [POINTER(Display), Drawable, POINTER(Window), POINTER(c_int), POINTER(c_int), POINTER(c_uint), POINTER(c_uint), POINTER(c_uint), POINTER(c_uint)]

# /usr/include/X11/Xlib.h:2634
XGetIconName = _lib.XGetIconName
XGetIconName.restype = c_int
XGetIconName.argtypes = [POINTER(Display), Window, POINTER(c_char_p)]

# /usr/include/X11/Xlib.h:2640
XGetInputFocus = _lib.XGetInputFocus
XGetInputFocus.restype = c_int
XGetInputFocus.argtypes = [POINTER(Display), POINTER(Window), POINTER(c_int)]

# /usr/include/X11/Xlib.h:2646
XGetKeyboardControl = _lib.XGetKeyboardControl
XGetKeyboardControl.restype = c_int
XGetKeyboardControl.argtypes = [POINTER(Display), POINTER(XKeyboardState)]

# /usr/include/X11/Xlib.h:2651
XGetPointerControl = _lib.XGetPointerControl
XGetPointerControl.restype = c_int
XGetPointerControl.argtypes = [POINTER(Display), POINTER(c_int), POINTER(c_int), POINTER(c_int)]

# /usr/include/X11/Xlib.h:2658
XGetPointerMapping = _lib.XGetPointerMapping
XGetPointerMapping.restype = c_int
XGetPointerMapping.argtypes = [POINTER(Display), POINTER(c_ubyte), c_int]

# /usr/include/X11/Xlib.h:2664
XGetScreenSaver = _lib.XGetScreenSaver
XGetScreenSaver.restype = c_int
XGetScreenSaver.argtypes = [POINTER(Display), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int)]

# /usr/include/X11/Xlib.h:2672
XGetTransientForHint = _lib.XGetTransientForHint
XGetTransientForHint.restype = c_int
XGetTransientForHint.argtypes = [POINTER(Display), Window, POINTER(Window)]

# /usr/include/X11/Xlib.h:2678
XGetWindowProperty = _lib.XGetWindowProperty
XGetWindowProperty.restype = c_int
XGetWindowProperty.argtypes = [POINTER(Display), Window, Atom, c_long, c_long, c_int, Atom, POINTER(Atom), POINTER(c_int), POINTER(c_ulong), POINTER(c_ulong), POINTER(POINTER(c_ubyte))]

# /usr/include/X11/Xlib.h:2693
XGetWindowAttributes = _lib.XGetWindowAttributes
XGetWindowAttributes.restype = c_int
XGetWindowAttributes.argtypes = [POINTER(Display), Window, POINTER(XWindowAttributes)]

# /usr/include/X11/Xlib.h:2699
XGrabButton = _lib.XGrabButton
XGrabButton.restype = c_int
XGrabButton.argtypes = [POINTER(Display), c_uint, c_uint, Window, c_int, c_uint, c_int, c_int, Window, Cursor]

# /usr/include/X11/Xlib.h:2712
XGrabKey = _lib.XGrabKey
XGrabKey.restype = c_int
XGrabKey.argtypes = [POINTER(Display), c_int, c_uint, Window, c_int, c_int, c_int]

# /usr/include/X11/Xlib.h:2722
XGrabKeyboard = _lib.XGrabKeyboard
XGrabKeyboard.restype = c_int
XGrabKeyboard.argtypes = [POINTER(Display), Window, c_int, c_int, c_int, Time]

# /usr/include/X11/Xlib.h:2731
XGrabPointer = _lib.XGrabPointer
XGrabPointer.restype = c_int
XGrabPointer.argtypes = [POINTER(Display), Window, c_int, c_uint, c_int, c_int, Window, Cursor, Time]

# /usr/include/X11/Xlib.h:2743
XGrabServer = _lib.XGrabServer
XGrabServer.restype = c_int
XGrabServer.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2747
XHeightMMOfScreen = _lib.XHeightMMOfScreen
XHeightMMOfScreen.restype = c_int
XHeightMMOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:2751
XHeightOfScreen = _lib.XHeightOfScreen
XHeightOfScreen.restype = c_int
XHeightOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:2755
XIfEvent = _lib.XIfEvent
XIfEvent.restype = c_int
XIfEvent.argtypes = [POINTER(Display), POINTER(XEvent), CFUNCTYPE(c_int, POINTER(Display), POINTER(XEvent), XPointer), XPointer]

# /usr/include/X11/Xlib.h:2766
XImageByteOrder = _lib.XImageByteOrder
XImageByteOrder.restype = c_int
XImageByteOrder.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2770
XInstallColormap = _lib.XInstallColormap
XInstallColormap.restype = c_int
XInstallColormap.argtypes = [POINTER(Display), Colormap]

# /usr/include/X11/Xlib.h:2775
XKeysymToKeycode = _lib.XKeysymToKeycode
XKeysymToKeycode.restype = KeyCode
XKeysymToKeycode.argtypes = [POINTER(Display), KeySym]

# /usr/include/X11/Xlib.h:2780
XKillClient = _lib.XKillClient
XKillClient.restype = c_int
XKillClient.argtypes = [POINTER(Display), XID]

# /usr/include/X11/Xlib.h:2785
XLookupColor = _lib.XLookupColor
XLookupColor.restype = c_int
XLookupColor.argtypes = [POINTER(Display), Colormap, c_char_p, POINTER(XColor), POINTER(XColor)]

# /usr/include/X11/Xlib.h:2793
XLowerWindow = _lib.XLowerWindow
XLowerWindow.restype = c_int
XLowerWindow.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:2798
XMapRaised = _lib.XMapRaised
XMapRaised.restype = c_int
XMapRaised.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:2803
XMapSubwindows = _lib.XMapSubwindows
XMapSubwindows.restype = c_int
XMapSubwindows.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:2808
XMapWindow = _lib.XMapWindow
XMapWindow.restype = c_int
XMapWindow.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:2813
XMaskEvent = _lib.XMaskEvent
XMaskEvent.restype = c_int
XMaskEvent.argtypes = [POINTER(Display), c_long, POINTER(XEvent)]

# /usr/include/X11/Xlib.h:2819
XMaxCmapsOfScreen = _lib.XMaxCmapsOfScreen
XMaxCmapsOfScreen.restype = c_int
XMaxCmapsOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:2823
XMinCmapsOfScreen = _lib.XMinCmapsOfScreen
XMinCmapsOfScreen.restype = c_int
XMinCmapsOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:2827
XMoveResizeWindow = _lib.XMoveResizeWindow
XMoveResizeWindow.restype = c_int
XMoveResizeWindow.argtypes = [POINTER(Display), Window, c_int, c_int, c_uint, c_uint]

# /usr/include/X11/Xlib.h:2836
XMoveWindow = _lib.XMoveWindow
XMoveWindow.restype = c_int
XMoveWindow.argtypes = [POINTER(Display), Window, c_int, c_int]

# /usr/include/X11/Xlib.h:2843
XNextEvent = _lib.XNextEvent
XNextEvent.restype = c_int
XNextEvent.argtypes = [POINTER(Display), POINTER(XEvent)]

# /usr/include/X11/Xlib.h:2848
XNoOp = _lib.XNoOp
XNoOp.restype = c_int
XNoOp.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2852
XParseColor = _lib.XParseColor
XParseColor.restype = c_int
XParseColor.argtypes = [POINTER(Display), Colormap, c_char_p, POINTER(XColor)]

# /usr/include/X11/Xlib.h:2859
XParseGeometry = _lib.XParseGeometry
XParseGeometry.restype = c_int
XParseGeometry.argtypes = [c_char_p, POINTER(c_int), POINTER(c_int), POINTER(c_uint), POINTER(c_uint)]

# /usr/include/X11/Xlib.h:2867
XPeekEvent = _lib.XPeekEvent
XPeekEvent.restype = c_int
XPeekEvent.argtypes = [POINTER(Display), POINTER(XEvent)]

# /usr/include/X11/Xlib.h:2872
XPeekIfEvent = _lib.XPeekIfEvent
XPeekIfEvent.restype = c_int
XPeekIfEvent.argtypes = [POINTER(Display), POINTER(XEvent), CFUNCTYPE(c_int, POINTER(Display), POINTER(XEvent), XPointer), XPointer]

# /usr/include/X11/Xlib.h:2883
XPending = _lib.XPending
XPending.restype = c_int
XPending.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2887
XPlanesOfScreen = _lib.XPlanesOfScreen
XPlanesOfScreen.restype = c_int
XPlanesOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:2891
XProtocolRevision = _lib.XProtocolRevision
XProtocolRevision.restype = c_int
XProtocolRevision.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2895
XProtocolVersion = _lib.XProtocolVersion
XProtocolVersion.restype = c_int
XProtocolVersion.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2900
XPutBackEvent = _lib.XPutBackEvent
XPutBackEvent.restype = c_int
XPutBackEvent.argtypes = [POINTER(Display), POINTER(XEvent)]

# /usr/include/X11/Xlib.h:2905
XPutImage = _lib.XPutImage
XPutImage.restype = c_int
XPutImage.argtypes = [POINTER(Display), Drawable, GC, POINTER(XImage), c_int, c_int, c_int, c_int, c_uint, c_uint]

# /usr/include/X11/Xlib.h:2918
XQLength = _lib.XQLength
XQLength.restype = c_int
XQLength.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:2922
XQueryBestCursor = _lib.XQueryBestCursor
XQueryBestCursor.restype = c_int
XQueryBestCursor.argtypes = [POINTER(Display), Drawable, c_uint, c_uint, POINTER(c_uint), POINTER(c_uint)]

# /usr/include/X11/Xlib.h:2931
XQueryBestSize = _lib.XQueryBestSize
XQueryBestSize.restype = c_int
XQueryBestSize.argtypes = [POINTER(Display), c_int, Drawable, c_uint, c_uint, POINTER(c_uint), POINTER(c_uint)]

# /usr/include/X11/Xlib.h:2941
XQueryBestStipple = _lib.XQueryBestStipple
XQueryBestStipple.restype = c_int
XQueryBestStipple.argtypes = [POINTER(Display), Drawable, c_uint, c_uint, POINTER(c_uint), POINTER(c_uint)]

# /usr/include/X11/Xlib.h:2950
XQueryBestTile = _lib.XQueryBestTile
XQueryBestTile.restype = c_int
XQueryBestTile.argtypes = [POINTER(Display), Drawable, c_uint, c_uint, POINTER(c_uint), POINTER(c_uint)]

# /usr/include/X11/Xlib.h:2959
XQueryColor = _lib.XQueryColor
XQueryColor.restype = c_int
XQueryColor.argtypes = [POINTER(Display), Colormap, POINTER(XColor)]

# /usr/include/X11/Xlib.h:2965
XQueryColors = _lib.XQueryColors
XQueryColors.restype = c_int
XQueryColors.argtypes = [POINTER(Display), Colormap, POINTER(XColor), c_int]

# /usr/include/X11/Xlib.h:2972
XQueryExtension = _lib.XQueryExtension
XQueryExtension.restype = c_int
XQueryExtension.argtypes = [POINTER(Display), c_char_p, POINTER(c_int), POINTER(c_int), POINTER(c_int)]

# /usr/include/X11/Xlib.h:2980
XQueryKeymap = _lib.XQueryKeymap
XQueryKeymap.restype = c_int
XQueryKeymap.argtypes = [POINTER(Display), c_char * 32]

# /usr/include/X11/Xlib.h:2985
XQueryPointer = _lib.XQueryPointer
XQueryPointer.restype = c_int
XQueryPointer.argtypes = [POINTER(Display), Window, POINTER(Window), POINTER(Window), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_uint)]

# /usr/include/X11/Xlib.h:2997
XQueryTextExtents = _lib.XQueryTextExtents
XQueryTextExtents.restype = c_int
XQueryTextExtents.argtypes = [POINTER(Display), XID, c_char_p, c_int, POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(XCharStruct)]

# /usr/include/X11/Xlib.h:3008
XQueryTextExtents16 = _lib.XQueryTextExtents16
XQueryTextExtents16.restype = c_int
XQueryTextExtents16.argtypes = [POINTER(Display), XID, POINTER(XChar2b), c_int, POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(XCharStruct)]

# /usr/include/X11/Xlib.h:3019
XQueryTree = _lib.XQueryTree
XQueryTree.restype = c_int
XQueryTree.argtypes = [POINTER(Display), Window, POINTER(Window), POINTER(Window), POINTER(POINTER(Window)), POINTER(c_uint)]

# /usr/include/X11/Xlib.h:3028
XRaiseWindow = _lib.XRaiseWindow
XRaiseWindow.restype = c_int
XRaiseWindow.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:3033
XReadBitmapFile = _lib.XReadBitmapFile
XReadBitmapFile.restype = c_int
XReadBitmapFile.argtypes = [POINTER(Display), Drawable, c_char_p, POINTER(c_uint), POINTER(c_uint), POINTER(Pixmap), POINTER(c_int), POINTER(c_int)]

# /usr/include/X11/Xlib.h:3044
XReadBitmapFileData = _lib.XReadBitmapFileData
XReadBitmapFileData.restype = c_int
XReadBitmapFileData.argtypes = [c_char_p, POINTER(c_uint), POINTER(c_uint), POINTER(POINTER(c_ubyte)), POINTER(c_int), POINTER(c_int)]

# /usr/include/X11/Xlib.h:3053
XRebindKeysym = _lib.XRebindKeysym
XRebindKeysym.restype = c_int
XRebindKeysym.argtypes = [POINTER(Display), KeySym, POINTER(KeySym), c_int, POINTER(c_ubyte), c_int]

# /usr/include/X11/Xlib.h:3062
XRecolorCursor = _lib.XRecolorCursor
XRecolorCursor.restype = c_int
XRecolorCursor.argtypes = [POINTER(Display), Cursor, POINTER(XColor), POINTER(XColor)]

# /usr/include/X11/Xlib.h:3069
XRefreshKeyboardMapping = _lib.XRefreshKeyboardMapping
XRefreshKeyboardMapping.restype = c_int
XRefreshKeyboardMapping.argtypes = [POINTER(XMappingEvent)]

# /usr/include/X11/Xlib.h:3073
XRemoveFromSaveSet = _lib.XRemoveFromSaveSet
XRemoveFromSaveSet.restype = c_int
XRemoveFromSaveSet.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:3078
XRemoveHost = _lib.XRemoveHost
XRemoveHost.restype = c_int
XRemoveHost.argtypes = [POINTER(Display), POINTER(XHostAddress)]

# /usr/include/X11/Xlib.h:3083
XRemoveHosts = _lib.XRemoveHosts
XRemoveHosts.restype = c_int
XRemoveHosts.argtypes = [POINTER(Display), POINTER(XHostAddress), c_int]

# /usr/include/X11/Xlib.h:3089
XReparentWindow = _lib.XReparentWindow
XReparentWindow.restype = c_int
XReparentWindow.argtypes = [POINTER(Display), Window, Window, c_int, c_int]

# /usr/include/X11/Xlib.h:3097
XResetScreenSaver = _lib.XResetScreenSaver
XResetScreenSaver.restype = c_int
XResetScreenSaver.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:3101
XResizeWindow = _lib.XResizeWindow
XResizeWindow.restype = c_int
XResizeWindow.argtypes = [POINTER(Display), Window, c_uint, c_uint]

# /usr/include/X11/Xlib.h:3108
XRestackWindows = _lib.XRestackWindows
XRestackWindows.restype = c_int
XRestackWindows.argtypes = [POINTER(Display), POINTER(Window), c_int]

# /usr/include/X11/Xlib.h:3114
XRotateBuffers = _lib.XRotateBuffers
XRotateBuffers.restype = c_int
XRotateBuffers.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:3119
XRotateWindowProperties = _lib.XRotateWindowProperties
XRotateWindowProperties.restype = c_int
XRotateWindowProperties.argtypes = [POINTER(Display), Window, POINTER(Atom), c_int, c_int]

# /usr/include/X11/Xlib.h:3127
XScreenCount = _lib.XScreenCount
XScreenCount.restype = c_int
XScreenCount.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:3131
XSelectInput = _lib.XSelectInput
XSelectInput.restype = c_int
XSelectInput.argtypes = [POINTER(Display), Window, c_long]

# /usr/include/X11/Xlib.h:3137
XSendEvent = _lib.XSendEvent
XSendEvent.restype = c_int
XSendEvent.argtypes = [POINTER(Display), Window, c_int, c_long, POINTER(XEvent)]

# /usr/include/X11/Xlib.h:3145
XSetAccessControl = _lib.XSetAccessControl
XSetAccessControl.restype = c_int
XSetAccessControl.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:3150
XSetArcMode = _lib.XSetArcMode
XSetArcMode.restype = c_int
XSetArcMode.argtypes = [POINTER(Display), GC, c_int]

# /usr/include/X11/Xlib.h:3156
XSetBackground = _lib.XSetBackground
XSetBackground.restype = c_int
XSetBackground.argtypes = [POINTER(Display), GC, c_ulong]

# /usr/include/X11/Xlib.h:3162
XSetClipMask = _lib.XSetClipMask
XSetClipMask.restype = c_int
XSetClipMask.argtypes = [POINTER(Display), GC, Pixmap]

# /usr/include/X11/Xlib.h:3168
XSetClipOrigin = _lib.XSetClipOrigin
XSetClipOrigin.restype = c_int
XSetClipOrigin.argtypes = [POINTER(Display), GC, c_int, c_int]

# /usr/include/X11/Xlib.h:3175
XSetClipRectangles = _lib.XSetClipRectangles
XSetClipRectangles.restype = c_int
XSetClipRectangles.argtypes = [POINTER(Display), GC, c_int, c_int, POINTER(XRectangle), c_int, c_int]

# /usr/include/X11/Xlib.h:3185
XSetCloseDownMode = _lib.XSetCloseDownMode
XSetCloseDownMode.restype = c_int
XSetCloseDownMode.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:3190
XSetCommand = _lib.XSetCommand
XSetCommand.restype = c_int
XSetCommand.argtypes = [POINTER(Display), Window, POINTER(c_char_p), c_int]

# /usr/include/X11/Xlib.h:3197
XSetDashes = _lib.XSetDashes
XSetDashes.restype = c_int
XSetDashes.argtypes = [POINTER(Display), GC, c_int, c_char_p, c_int]

# /usr/include/X11/Xlib.h:3205
XSetFillRule = _lib.XSetFillRule
XSetFillRule.restype = c_int
XSetFillRule.argtypes = [POINTER(Display), GC, c_int]

# /usr/include/X11/Xlib.h:3211
XSetFillStyle = _lib.XSetFillStyle
XSetFillStyle.restype = c_int
XSetFillStyle.argtypes = [POINTER(Display), GC, c_int]

# /usr/include/X11/Xlib.h:3217
XSetFont = _lib.XSetFont
XSetFont.restype = c_int
XSetFont.argtypes = [POINTER(Display), GC, Font]

# /usr/include/X11/Xlib.h:3223
XSetFontPath = _lib.XSetFontPath
XSetFontPath.restype = c_int
XSetFontPath.argtypes = [POINTER(Display), POINTER(c_char_p), c_int]

# /usr/include/X11/Xlib.h:3229
XSetForeground = _lib.XSetForeground
XSetForeground.restype = c_int
XSetForeground.argtypes = [POINTER(Display), GC, c_ulong]

# /usr/include/X11/Xlib.h:3235
XSetFunction = _lib.XSetFunction
XSetFunction.restype = c_int
XSetFunction.argtypes = [POINTER(Display), GC, c_int]

# /usr/include/X11/Xlib.h:3241
XSetGraphicsExposures = _lib.XSetGraphicsExposures
XSetGraphicsExposures.restype = c_int
XSetGraphicsExposures.argtypes = [POINTER(Display), GC, c_int]

# /usr/include/X11/Xlib.h:3247
XSetIconName = _lib.XSetIconName
XSetIconName.restype = c_int
XSetIconName.argtypes = [POINTER(Display), Window, c_char_p]

# /usr/include/X11/Xlib.h:3253
XSetInputFocus = _lib.XSetInputFocus
XSetInputFocus.restype = c_int
XSetInputFocus.argtypes = [POINTER(Display), Window, c_int, Time]

# /usr/include/X11/Xlib.h:3260
XSetLineAttributes = _lib.XSetLineAttributes
XSetLineAttributes.restype = c_int
XSetLineAttributes.argtypes = [POINTER(Display), GC, c_uint, c_int, c_int, c_int]

# /usr/include/X11/Xlib.h:3269
XSetModifierMapping = _lib.XSetModifierMapping
XSetModifierMapping.restype = c_int
XSetModifierMapping.argtypes = [POINTER(Display), POINTER(XModifierKeymap)]

# /usr/include/X11/Xlib.h:3274
XSetPlaneMask = _lib.XSetPlaneMask
XSetPlaneMask.restype = c_int
XSetPlaneMask.argtypes = [POINTER(Display), GC, c_ulong]

# /usr/include/X11/Xlib.h:3280
XSetPointerMapping = _lib.XSetPointerMapping
XSetPointerMapping.restype = c_int
XSetPointerMapping.argtypes = [POINTER(Display), POINTER(c_ubyte), c_int]

# /usr/include/X11/Xlib.h:3286
XSetScreenSaver = _lib.XSetScreenSaver
XSetScreenSaver.restype = c_int
XSetScreenSaver.argtypes = [POINTER(Display), c_int, c_int, c_int, c_int]

# /usr/include/X11/Xlib.h:3294
XSetSelectionOwner = _lib.XSetSelectionOwner
XSetSelectionOwner.restype = c_int
XSetSelectionOwner.argtypes = [POINTER(Display), Atom, Window, Time]

# /usr/include/X11/Xlib.h:3301
XSetState = _lib.XSetState
XSetState.restype = c_int
XSetState.argtypes = [POINTER(Display), GC, c_ulong, c_ulong, c_int, c_ulong]

# /usr/include/X11/Xlib.h:3310
XSetStipple = _lib.XSetStipple
XSetStipple.restype = c_int
XSetStipple.argtypes = [POINTER(Display), GC, Pixmap]

# /usr/include/X11/Xlib.h:3316
XSetSubwindowMode = _lib.XSetSubwindowMode
XSetSubwindowMode.restype = c_int
XSetSubwindowMode.argtypes = [POINTER(Display), GC, c_int]

# /usr/include/X11/Xlib.h:3322
XSetTSOrigin = _lib.XSetTSOrigin
XSetTSOrigin.restype = c_int
XSetTSOrigin.argtypes = [POINTER(Display), GC, c_int, c_int]

# /usr/include/X11/Xlib.h:3329
XSetTile = _lib.XSetTile
XSetTile.restype = c_int
XSetTile.argtypes = [POINTER(Display), GC, Pixmap]

# /usr/include/X11/Xlib.h:3335
XSetWindowBackground = _lib.XSetWindowBackground
XSetWindowBackground.restype = c_int
XSetWindowBackground.argtypes = [POINTER(Display), Window, c_ulong]

# /usr/include/X11/Xlib.h:3341
XSetWindowBackgroundPixmap = _lib.XSetWindowBackgroundPixmap
XSetWindowBackgroundPixmap.restype = c_int
XSetWindowBackgroundPixmap.argtypes = [POINTER(Display), Window, Pixmap]

# /usr/include/X11/Xlib.h:3347
XSetWindowBorder = _lib.XSetWindowBorder
XSetWindowBorder.restype = c_int
XSetWindowBorder.argtypes = [POINTER(Display), Window, c_ulong]

# /usr/include/X11/Xlib.h:3353
XSetWindowBorderPixmap = _lib.XSetWindowBorderPixmap
XSetWindowBorderPixmap.restype = c_int
XSetWindowBorderPixmap.argtypes = [POINTER(Display), Window, Pixmap]

# /usr/include/X11/Xlib.h:3359
XSetWindowBorderWidth = _lib.XSetWindowBorderWidth
XSetWindowBorderWidth.restype = c_int
XSetWindowBorderWidth.argtypes = [POINTER(Display), Window, c_uint]

# /usr/include/X11/Xlib.h:3365
XSetWindowColormap = _lib.XSetWindowColormap
XSetWindowColormap.restype = c_int
XSetWindowColormap.argtypes = [POINTER(Display), Window, Colormap]

# /usr/include/X11/Xlib.h:3371
XStoreBuffer = _lib.XStoreBuffer
XStoreBuffer.restype = c_int
XStoreBuffer.argtypes = [POINTER(Display), c_char_p, c_int, c_int]

# /usr/include/X11/Xlib.h:3378
XStoreBytes = _lib.XStoreBytes
XStoreBytes.restype = c_int
XStoreBytes.argtypes = [POINTER(Display), c_char_p, c_int]

# /usr/include/X11/Xlib.h:3384
XStoreColor = _lib.XStoreColor
XStoreColor.restype = c_int
XStoreColor.argtypes = [POINTER(Display), Colormap, POINTER(XColor)]

# /usr/include/X11/Xlib.h:3390
XStoreColors = _lib.XStoreColors
XStoreColors.restype = c_int
XStoreColors.argtypes = [POINTER(Display), Colormap, POINTER(XColor), c_int]

# /usr/include/X11/Xlib.h:3397
XStoreName = _lib.XStoreName
XStoreName.restype = c_int
XStoreName.argtypes = [POINTER(Display), Window, c_char_p]

# /usr/include/X11/Xlib.h:3403
XStoreNamedColor = _lib.XStoreNamedColor
XStoreNamedColor.restype = c_int
XStoreNamedColor.argtypes = [POINTER(Display), Colormap, c_char_p, c_ulong, c_int]

# /usr/include/X11/Xlib.h:3411
XSync = _lib.XSync
XSync.restype = c_int
XSync.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:3416
XTextExtents = _lib.XTextExtents
XTextExtents.restype = c_int
XTextExtents.argtypes = [POINTER(XFontStruct), c_char_p, c_int, POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(XCharStruct)]

# /usr/include/X11/Xlib.h:3426
XTextExtents16 = _lib.XTextExtents16
XTextExtents16.restype = c_int
XTextExtents16.argtypes = [POINTER(XFontStruct), POINTER(XChar2b), c_int, POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(XCharStruct)]

# /usr/include/X11/Xlib.h:3436
XTextWidth = _lib.XTextWidth
XTextWidth.restype = c_int
XTextWidth.argtypes = [POINTER(XFontStruct), c_char_p, c_int]

# /usr/include/X11/Xlib.h:3442
XTextWidth16 = _lib.XTextWidth16
XTextWidth16.restype = c_int
XTextWidth16.argtypes = [POINTER(XFontStruct), POINTER(XChar2b), c_int]

# /usr/include/X11/Xlib.h:3448
XTranslateCoordinates = _lib.XTranslateCoordinates
XTranslateCoordinates.restype = c_int
XTranslateCoordinates.argtypes = [POINTER(Display), Window, Window, c_int, c_int, POINTER(c_int), POINTER(c_int), POINTER(Window)]

# /usr/include/X11/Xlib.h:3459
XUndefineCursor = _lib.XUndefineCursor
XUndefineCursor.restype = c_int
XUndefineCursor.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:3464
XUngrabButton = _lib.XUngrabButton
XUngrabButton.restype = c_int
XUngrabButton.argtypes = [POINTER(Display), c_uint, c_uint, Window]

# /usr/include/X11/Xlib.h:3471
XUngrabKey = _lib.XUngrabKey
XUngrabKey.restype = c_int
XUngrabKey.argtypes = [POINTER(Display), c_int, c_uint, Window]

# /usr/include/X11/Xlib.h:3478
XUngrabKeyboard = _lib.XUngrabKeyboard
XUngrabKeyboard.restype = c_int
XUngrabKeyboard.argtypes = [POINTER(Display), Time]

# /usr/include/X11/Xlib.h:3483
XUngrabPointer = _lib.XUngrabPointer
XUngrabPointer.restype = c_int
XUngrabPointer.argtypes = [POINTER(Display), Time]

# /usr/include/X11/Xlib.h:3488
XUngrabServer = _lib.XUngrabServer
XUngrabServer.restype = c_int
XUngrabServer.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:3492
XUninstallColormap = _lib.XUninstallColormap
XUninstallColormap.restype = c_int
XUninstallColormap.argtypes = [POINTER(Display), Colormap]

# /usr/include/X11/Xlib.h:3497
XUnloadFont = _lib.XUnloadFont
XUnloadFont.restype = c_int
XUnloadFont.argtypes = [POINTER(Display), Font]

# /usr/include/X11/Xlib.h:3502
XUnmapSubwindows = _lib.XUnmapSubwindows
XUnmapSubwindows.restype = c_int
XUnmapSubwindows.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:3507
XUnmapWindow = _lib.XUnmapWindow
XUnmapWindow.restype = c_int
XUnmapWindow.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xlib.h:3512
XVendorRelease = _lib.XVendorRelease
XVendorRelease.restype = c_int
XVendorRelease.argtypes = [POINTER(Display)]

# /usr/include/X11/Xlib.h:3516
XWarpPointer = _lib.XWarpPointer
XWarpPointer.restype = c_int
XWarpPointer.argtypes = [POINTER(Display), Window, Window, c_int, c_int, c_uint, c_uint, c_int, c_int]

# /usr/include/X11/Xlib.h:3528
XWidthMMOfScreen = _lib.XWidthMMOfScreen
XWidthMMOfScreen.restype = c_int
XWidthMMOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:3532
XWidthOfScreen = _lib.XWidthOfScreen
XWidthOfScreen.restype = c_int
XWidthOfScreen.argtypes = [POINTER(Screen)]

# /usr/include/X11/Xlib.h:3536
XWindowEvent = _lib.XWindowEvent
XWindowEvent.restype = c_int
XWindowEvent.argtypes = [POINTER(Display), Window, c_long, POINTER(XEvent)]

# /usr/include/X11/Xlib.h:3543
XWriteBitmapFile = _lib.XWriteBitmapFile
XWriteBitmapFile.restype = c_int
XWriteBitmapFile.argtypes = [POINTER(Display), c_char_p, Pixmap, c_uint, c_uint, c_int, c_int]

# /usr/include/X11/Xlib.h:3553
XSupportsLocale = _lib.XSupportsLocale
XSupportsLocale.restype = c_int
XSupportsLocale.argtypes = []

# /usr/include/X11/Xlib.h:3555
XSetLocaleModifiers = _lib.XSetLocaleModifiers
XSetLocaleModifiers.restype = c_char_p
XSetLocaleModifiers.argtypes = [c_char_p]

# /usr/include/X11/Xlib.h:3559
XOpenOM = _lib.XOpenOM
XOpenOM.restype = XOM
XOpenOM.argtypes = [POINTER(Display), POINTER(struct__XrmHashBucketRec), c_char_p, c_char_p]

# /usr/include/X11/Xlib.h:3566
XCloseOM = _lib.XCloseOM
XCloseOM.restype = c_int
XCloseOM.argtypes = [XOM]

# /usr/include/X11/Xlib.h:3570
XSetOMValues = _lib.XSetOMValues
XSetOMValues.restype = c_char_p
XSetOMValues.argtypes = [XOM]

# /usr/include/X11/Xlib.h:3575
XGetOMValues = _lib.XGetOMValues
XGetOMValues.restype = c_char_p
XGetOMValues.argtypes = [XOM]

# /usr/include/X11/Xlib.h:3580
XDisplayOfOM = _lib.XDisplayOfOM
XDisplayOfOM.restype = POINTER(Display)
XDisplayOfOM.argtypes = [XOM]

# /usr/include/X11/Xlib.h:3584
XLocaleOfOM = _lib.XLocaleOfOM
XLocaleOfOM.restype = c_char_p
XLocaleOfOM.argtypes = [XOM]

# /usr/include/X11/Xlib.h:3588
XCreateOC = _lib.XCreateOC
XCreateOC.restype = XOC
XCreateOC.argtypes = [XOM]

# /usr/include/X11/Xlib.h:3593
XDestroyOC = _lib.XDestroyOC
XDestroyOC.restype = None
XDestroyOC.argtypes = [XOC]

# /usr/include/X11/Xlib.h:3597
XOMOfOC = _lib.XOMOfOC
XOMOfOC.restype = XOM
XOMOfOC.argtypes = [XOC]

# /usr/include/X11/Xlib.h:3601
XSetOCValues = _lib.XSetOCValues
XSetOCValues.restype = c_char_p
XSetOCValues.argtypes = [XOC]

# /usr/include/X11/Xlib.h:3606
XGetOCValues = _lib.XGetOCValues
XGetOCValues.restype = c_char_p
XGetOCValues.argtypes = [XOC]

# /usr/include/X11/Xlib.h:3611
XCreateFontSet = _lib.XCreateFontSet
XCreateFontSet.restype = XFontSet
XCreateFontSet.argtypes = [POINTER(Display), c_char_p, POINTER(POINTER(c_char_p)), POINTER(c_int), POINTER(c_char_p)]

# /usr/include/X11/Xlib.h:3619
XFreeFontSet = _lib.XFreeFontSet
XFreeFontSet.restype = None
XFreeFontSet.argtypes = [POINTER(Display), XFontSet]

# /usr/include/X11/Xlib.h:3624
XFontsOfFontSet = _lib.XFontsOfFontSet
XFontsOfFontSet.restype = c_int
XFontsOfFontSet.argtypes = [XFontSet, POINTER(POINTER(POINTER(XFontStruct))), POINTER(POINTER(c_char_p))]

# /usr/include/X11/Xlib.h:3630
XBaseFontNameListOfFontSet = _lib.XBaseFontNameListOfFontSet
XBaseFontNameListOfFontSet.restype = c_char_p
XBaseFontNameListOfFontSet.argtypes = [XFontSet]

# /usr/include/X11/Xlib.h:3634
XLocaleOfFontSet = _lib.XLocaleOfFontSet
XLocaleOfFontSet.restype = c_char_p
XLocaleOfFontSet.argtypes = [XFontSet]

# /usr/include/X11/Xlib.h:3638
XContextDependentDrawing = _lib.XContextDependentDrawing
XContextDependentDrawing.restype = c_int
XContextDependentDrawing.argtypes = [XFontSet]

# /usr/include/X11/Xlib.h:3642
XDirectionalDependentDrawing = _lib.XDirectionalDependentDrawing
XDirectionalDependentDrawing.restype = c_int
XDirectionalDependentDrawing.argtypes = [XFontSet]

# /usr/include/X11/Xlib.h:3646
XContextualDrawing = _lib.XContextualDrawing
XContextualDrawing.restype = c_int
XContextualDrawing.argtypes = [XFontSet]

# /usr/include/X11/Xlib.h:3650
XExtentsOfFontSet = _lib.XExtentsOfFontSet
XExtentsOfFontSet.restype = POINTER(XFontSetExtents)
XExtentsOfFontSet.argtypes = [XFontSet]

# /usr/include/X11/Xlib.h:3654
XmbTextEscapement = _lib.XmbTextEscapement
XmbTextEscapement.restype = c_int
XmbTextEscapement.argtypes = [XFontSet, c_char_p, c_int]

# /usr/include/X11/Xlib.h:3660
XwcTextEscapement = _lib.XwcTextEscapement
XwcTextEscapement.restype = c_int
XwcTextEscapement.argtypes = [XFontSet, c_wchar_p, c_int]

# /usr/include/X11/Xlib.h:3666
Xutf8TextEscapement = _lib.Xutf8TextEscapement
Xutf8TextEscapement.restype = c_int
Xutf8TextEscapement.argtypes = [XFontSet, c_char_p, c_int]

# /usr/include/X11/Xlib.h:3672
XmbTextExtents = _lib.XmbTextExtents
XmbTextExtents.restype = c_int
XmbTextExtents.argtypes = [XFontSet, c_char_p, c_int, POINTER(XRectangle), POINTER(XRectangle)]

# /usr/include/X11/Xlib.h:3680
XwcTextExtents = _lib.XwcTextExtents
XwcTextExtents.restype = c_int
XwcTextExtents.argtypes = [XFontSet, c_wchar_p, c_int, POINTER(XRectangle), POINTER(XRectangle)]

# /usr/include/X11/Xlib.h:3688
Xutf8TextExtents = _lib.Xutf8TextExtents
Xutf8TextExtents.restype = c_int
Xutf8TextExtents.argtypes = [XFontSet, c_char_p, c_int, POINTER(XRectangle), POINTER(XRectangle)]

# /usr/include/X11/Xlib.h:3696
XmbTextPerCharExtents = _lib.XmbTextPerCharExtents
XmbTextPerCharExtents.restype = c_int
XmbTextPerCharExtents.argtypes = [XFontSet, c_char_p, c_int, POINTER(XRectangle), POINTER(XRectangle), c_int, POINTER(c_int), POINTER(XRectangle), POINTER(XRectangle)]

# /usr/include/X11/Xlib.h:3708
XwcTextPerCharExtents = _lib.XwcTextPerCharExtents
XwcTextPerCharExtents.restype = c_int
XwcTextPerCharExtents.argtypes = [XFontSet, c_wchar_p, c_int, POINTER(XRectangle), POINTER(XRectangle), c_int, POINTER(c_int), POINTER(XRectangle), POINTER(XRectangle)]

# /usr/include/X11/Xlib.h:3720
Xutf8TextPerCharExtents = _lib.Xutf8TextPerCharExtents
Xutf8TextPerCharExtents.restype = c_int
Xutf8TextPerCharExtents.argtypes = [XFontSet, c_char_p, c_int, POINTER(XRectangle), POINTER(XRectangle), c_int, POINTER(c_int), POINTER(XRectangle), POINTER(XRectangle)]

# /usr/include/X11/Xlib.h:3732
XmbDrawText = _lib.XmbDrawText
XmbDrawText.restype = None
XmbDrawText.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, POINTER(XmbTextItem), c_int]

# /usr/include/X11/Xlib.h:3742
XwcDrawText = _lib.XwcDrawText
XwcDrawText.restype = None
XwcDrawText.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, POINTER(XwcTextItem), c_int]

# /usr/include/X11/Xlib.h:3752
Xutf8DrawText = _lib.Xutf8DrawText
Xutf8DrawText.restype = None
Xutf8DrawText.argtypes = [POINTER(Display), Drawable, GC, c_int, c_int, POINTER(XmbTextItem), c_int]

# /usr/include/X11/Xlib.h:3762
XmbDrawString = _lib.XmbDrawString
XmbDrawString.restype = None
XmbDrawString.argtypes = [POINTER(Display), Drawable, XFontSet, GC, c_int, c_int, c_char_p, c_int]

# /usr/include/X11/Xlib.h:3773
XwcDrawString = _lib.XwcDrawString
XwcDrawString.restype = None
XwcDrawString.argtypes = [POINTER(Display), Drawable, XFontSet, GC, c_int, c_int, c_wchar_p, c_int]

# /usr/include/X11/Xlib.h:3784
Xutf8DrawString = _lib.Xutf8DrawString
Xutf8DrawString.restype = None
Xutf8DrawString.argtypes = [POINTER(Display), Drawable, XFontSet, GC, c_int, c_int, c_char_p, c_int]

# /usr/include/X11/Xlib.h:3795
XmbDrawImageString = _lib.XmbDrawImageString
XmbDrawImageString.restype = None
XmbDrawImageString.argtypes = [POINTER(Display), Drawable, XFontSet, GC, c_int, c_int, c_char_p, c_int]

# /usr/include/X11/Xlib.h:3806
XwcDrawImageString = _lib.XwcDrawImageString
XwcDrawImageString.restype = None
XwcDrawImageString.argtypes = [POINTER(Display), Drawable, XFontSet, GC, c_int, c_int, c_wchar_p, c_int]

# /usr/include/X11/Xlib.h:3817
Xutf8DrawImageString = _lib.Xutf8DrawImageString
Xutf8DrawImageString.restype = None
Xutf8DrawImageString.argtypes = [POINTER(Display), Drawable, XFontSet, GC, c_int, c_int, c_char_p, c_int]

# /usr/include/X11/Xlib.h:3828
XOpenIM = _lib.XOpenIM
XOpenIM.restype = XIM
XOpenIM.argtypes = [POINTER(Display), POINTER(struct__XrmHashBucketRec), c_char_p, c_char_p]

# /usr/include/X11/Xlib.h:3835
XCloseIM = _lib.XCloseIM
XCloseIM.restype = c_int
XCloseIM.argtypes = [XIM]

# /usr/include/X11/Xlib.h:3839
XGetIMValues = _lib.XGetIMValues
XGetIMValues.restype = c_char_p
XGetIMValues.argtypes = [XIM]

# /usr/include/X11/Xlib.h:3843
XSetIMValues = _lib.XSetIMValues
XSetIMValues.restype = c_char_p
XSetIMValues.argtypes = [XIM]

# /usr/include/X11/Xlib.h:3847
XDisplayOfIM = _lib.XDisplayOfIM
XDisplayOfIM.restype = POINTER(Display)
XDisplayOfIM.argtypes = [XIM]

# /usr/include/X11/Xlib.h:3851
XLocaleOfIM = _lib.XLocaleOfIM
XLocaleOfIM.restype = c_char_p
XLocaleOfIM.argtypes = [XIM]

# /usr/include/X11/Xlib.h:3855
XCreateIC = _lib.XCreateIC
XCreateIC.restype = XIC
XCreateIC.argtypes = [XIM]

# /usr/include/X11/Xlib.h:3859
XDestroyIC = _lib.XDestroyIC
XDestroyIC.restype = None
XDestroyIC.argtypes = [XIC]

# /usr/include/X11/Xlib.h:3863
XSetICFocus = _lib.XSetICFocus
XSetICFocus.restype = None
XSetICFocus.argtypes = [XIC]

# /usr/include/X11/Xlib.h:3867
XUnsetICFocus = _lib.XUnsetICFocus
XUnsetICFocus.restype = None
XUnsetICFocus.argtypes = [XIC]

# /usr/include/X11/Xlib.h:3871
XwcResetIC = _lib.XwcResetIC
XwcResetIC.restype = c_wchar_p
XwcResetIC.argtypes = [XIC]

# /usr/include/X11/Xlib.h:3875
XmbResetIC = _lib.XmbResetIC
XmbResetIC.restype = c_char_p
XmbResetIC.argtypes = [XIC]

# /usr/include/X11/Xlib.h:3879
Xutf8ResetIC = _lib.Xutf8ResetIC
Xutf8ResetIC.restype = c_char_p
Xutf8ResetIC.argtypes = [XIC]

# /usr/include/X11/Xlib.h:3883
XSetICValues = _lib.XSetICValues
XSetICValues.restype = c_char_p
XSetICValues.argtypes = [XIC]

# /usr/include/X11/Xlib.h:3887
XGetICValues = _lib.XGetICValues
XGetICValues.restype = c_char_p
XGetICValues.argtypes = [XIC]

# /usr/include/X11/Xlib.h:3891
XIMOfIC = _lib.XIMOfIC
XIMOfIC.restype = XIM
XIMOfIC.argtypes = [XIC]

# /usr/include/X11/Xlib.h:3895
XFilterEvent = _lib.XFilterEvent
XFilterEvent.restype = c_int
XFilterEvent.argtypes = [POINTER(XEvent), Window]

# /usr/include/X11/Xlib.h:3900
XmbLookupString = _lib.XmbLookupString
XmbLookupString.restype = c_int
XmbLookupString.argtypes = [XIC, POINTER(XKeyPressedEvent), c_char_p, c_int, POINTER(KeySym), POINTER(c_int)]

# /usr/include/X11/Xlib.h:3909
XwcLookupString = _lib.XwcLookupString
XwcLookupString.restype = c_int
XwcLookupString.argtypes = [XIC, POINTER(XKeyPressedEvent), c_wchar_p, c_int, POINTER(KeySym), POINTER(c_int)]

# /usr/include/X11/Xlib.h:3918
Xutf8LookupString = _lib.Xutf8LookupString
Xutf8LookupString.restype = c_int
Xutf8LookupString.argtypes = [XIC, POINTER(XKeyPressedEvent), c_char_p, c_int, POINTER(KeySym), POINTER(c_int)]

# /usr/include/X11/Xlib.h:3927
XVaCreateNestedList = _lib.XVaCreateNestedList
XVaCreateNestedList.restype = XVaNestedList
XVaCreateNestedList.argtypes = [c_int]

# /usr/include/X11/Xlib.h:3933
XRegisterIMInstantiateCallback = _lib.XRegisterIMInstantiateCallback
XRegisterIMInstantiateCallback.restype = c_int
XRegisterIMInstantiateCallback.argtypes = [POINTER(Display), POINTER(struct__XrmHashBucketRec), c_char_p, c_char_p, XIDProc, XPointer]

# /usr/include/X11/Xlib.h:3942
XUnregisterIMInstantiateCallback = _lib.XUnregisterIMInstantiateCallback
XUnregisterIMInstantiateCallback.restype = c_int
XUnregisterIMInstantiateCallback.argtypes = [POINTER(Display), POINTER(struct__XrmHashBucketRec), c_char_p, c_char_p, XIDProc, XPointer]

# /usr/include/X11/Xlib.h:3960
XInternalConnectionNumbers = _lib.XInternalConnectionNumbers
XInternalConnectionNumbers.restype = c_int
XInternalConnectionNumbers.argtypes = [POINTER(Display), POINTER(POINTER(c_int)), POINTER(c_int)]

# /usr/include/X11/Xlib.h:3966
XProcessInternalConnection = _lib.XProcessInternalConnection
XProcessInternalConnection.restype = None
XProcessInternalConnection.argtypes = [POINTER(Display), c_int]

# /usr/include/X11/Xlib.h:3971
XAddConnectionWatch = _lib.XAddConnectionWatch
XAddConnectionWatch.restype = c_int
XAddConnectionWatch.argtypes = [POINTER(Display), XConnectionWatchProc, XPointer]

# /usr/include/X11/Xlib.h:3977
XRemoveConnectionWatch = _lib.XRemoveConnectionWatch
XRemoveConnectionWatch.restype = None
XRemoveConnectionWatch.argtypes = [POINTER(Display), XConnectionWatchProc, XPointer]

# /usr/include/X11/Xlib.h:3983
XSetAuthorization = _lib.XSetAuthorization
XSetAuthorization.restype = None
XSetAuthorization.argtypes = [c_char_p, c_int, c_char_p, c_int]

# /usr/include/X11/Xlib.h:3990
_Xmbtowc = _lib._Xmbtowc
_Xmbtowc.restype = c_int
_Xmbtowc.argtypes = [c_wchar_p, c_char_p, c_int]

# /usr/include/X11/Xlib.h:4001
_Xwctomb = _lib._Xwctomb
_Xwctomb.restype = c_int
_Xwctomb.argtypes = [c_char_p, c_wchar]

NoValue = 0 	# /usr/include/X11/Xutil.h:4791
XValue = 1 	# /usr/include/X11/Xutil.h:4792
YValue = 2 	# /usr/include/X11/Xutil.h:4793
WidthValue = 4 	# /usr/include/X11/Xutil.h:4794
HeightValue = 8 	# /usr/include/X11/Xutil.h:4795
AllValues = 0 	# /usr/include/X11/Xutil.h:4796
XNegative = 16 	# /usr/include/X11/Xutil.h:4797
YNegative = 32 	# /usr/include/X11/Xutil.h:4798
class struct_anon_86(Structure):
    __slots__ = [
        'flags',
        'x',
        'y',
        'width',
        'height',
        'min_width',
        'min_height',
        'max_width',
        'max_height',
        'width_inc',
        'height_inc',
        'min_aspect',
        'max_aspect',
        'base_width',
        'base_height',
        'win_gravity',
    ]
class struct_anon_87(Structure):
    __slots__ = [
        'x',
        'y',
    ]
struct_anon_87._fields_ = [
    ('x', c_int),
    ('y', c_int),
]

class struct_anon_88(Structure):
    __slots__ = [
        'x',
        'y',
    ]
struct_anon_88._fields_ = [
    ('x', c_int),
    ('y', c_int),
]

struct_anon_86._fields_ = [
    ('flags', c_long),
    ('x', c_int),
    ('y', c_int),
    ('width', c_int),
    ('height', c_int),
    ('min_width', c_int),
    ('min_height', c_int),
    ('max_width', c_int),
    ('max_height', c_int),
    ('width_inc', c_int),
    ('height_inc', c_int),
    ('min_aspect', struct_anon_87),
    ('max_aspect', struct_anon_88),
    ('base_width', c_int),
    ('base_height', c_int),
    ('win_gravity', c_int),
]

XSizeHints = struct_anon_86 	# /usr/include/X11/Xutil.h:4817
USPosition = 1 	# /usr/include/X11/Xutil.h:4825
USSize = 2 	# /usr/include/X11/Xutil.h:4826
PPosition = 4 	# /usr/include/X11/Xutil.h:4828
PSize = 8 	# /usr/include/X11/Xutil.h:4829
PMinSize = 16 	# /usr/include/X11/Xutil.h:4830
PMaxSize = 32 	# /usr/include/X11/Xutil.h:4831
PResizeInc = 64 	# /usr/include/X11/Xutil.h:4832
PAspect = 128 	# /usr/include/X11/Xutil.h:4833
PBaseSize = 256 	# /usr/include/X11/Xutil.h:4834
PWinGravity = 512 	# /usr/include/X11/Xutil.h:4835
PAllHints = 252 	# /usr/include/X11/Xutil.h:4838
class struct_anon_89(Structure):
    __slots__ = [
        'flags',
        'input',
        'initial_state',
        'icon_pixmap',
        'icon_window',
        'icon_x',
        'icon_y',
        'icon_mask',
        'window_group',
    ]
struct_anon_89._fields_ = [
    ('flags', c_long),
    ('input', c_int),
    ('initial_state', c_int),
    ('icon_pixmap', Pixmap),
    ('icon_window', Window),
    ('icon_x', c_int),
    ('icon_y', c_int),
    ('icon_mask', Pixmap),
    ('window_group', XID),
]

XWMHints = struct_anon_89 	# /usr/include/X11/Xutil.h:4853
InputHint = 1 	# /usr/include/X11/Xutil.h:4857
StateHint = 2 	# /usr/include/X11/Xutil.h:4858
IconPixmapHint = 4 	# /usr/include/X11/Xutil.h:4859
IconWindowHint = 8 	# /usr/include/X11/Xutil.h:4860
IconPositionHint = 16 	# /usr/include/X11/Xutil.h:4861
IconMaskHint = 32 	# /usr/include/X11/Xutil.h:4862
WindowGroupHint = 64 	# /usr/include/X11/Xutil.h:4863
AllHints = 127 	# /usr/include/X11/Xutil.h:4864
XUrgencyHint = 256 	# /usr/include/X11/Xutil.h:4866
WithdrawnState = 0 	# /usr/include/X11/Xutil.h:4869
NormalState = 1 	# /usr/include/X11/Xutil.h:4870
IconicState = 3 	# /usr/include/X11/Xutil.h:4871
DontCareState = 0 	# /usr/include/X11/Xutil.h:4876
ZoomState = 2 	# /usr/include/X11/Xutil.h:4877
InactiveState = 4 	# /usr/include/X11/Xutil.h:4878
class struct_anon_90(Structure):
    __slots__ = [
        'value',
        'encoding',
        'format',
        'nitems',
    ]
struct_anon_90._fields_ = [
    ('value', POINTER(c_ubyte)),
    ('encoding', Atom),
    ('format', c_int),
    ('nitems', c_ulong),
]

XTextProperty = struct_anon_90 	# /usr/include/X11/Xutil.h:4891
XNoMemory = -1 	# /usr/include/X11/Xutil.h:4893
XLocaleNotSupported = -2 	# /usr/include/X11/Xutil.h:4894
XConverterNotFound = -3 	# /usr/include/X11/Xutil.h:4895
enum_anon_91 = c_int
XStringStyle = 1
XCompoundTextStyle = 2
XTextStyle = 3
XStdICCTextStyle = 4
XUTF8StringStyle = 5
XICCEncodingStyle = enum_anon_91 	# /usr/include/X11/Xutil.h:4904
class struct_anon_92(Structure):
    __slots__ = [
        'min_width',
        'min_height',
        'max_width',
        'max_height',
        'width_inc',
        'height_inc',
    ]
struct_anon_92._fields_ = [
    ('min_width', c_int),
    ('min_height', c_int),
    ('max_width', c_int),
    ('max_height', c_int),
    ('width_inc', c_int),
    ('height_inc', c_int),
]

XIconSize = struct_anon_92 	# /usr/include/X11/Xutil.h:4910
class struct_anon_93(Structure):
    __slots__ = [
        'res_name',
        'res_class',
    ]
struct_anon_93._fields_ = [
    ('res_name', c_char_p),
    ('res_class', c_char_p),
]

XClassHint = struct_anon_93 	# /usr/include/X11/Xutil.h:4915
class struct__XComposeStatus(Structure):
    __slots__ = [
        'compose_ptr',
        'chars_matched',
    ]
struct__XComposeStatus._fields_ = [
    ('compose_ptr', XPointer),
    ('chars_matched', c_int),
]

XComposeStatus = struct__XComposeStatus 	# /usr/include/X11/Xutil.h:4957
class struct__XRegion(Structure):
    __slots__ = [
    ]
struct__XRegion._fields_ = [
    ('_opaque_struct', c_int)
]

Region = POINTER(struct__XRegion) 	# /usr/include/X11/Xutil.h:4996
RectangleOut = 0 	# /usr/include/X11/Xutil.h:5000
RectangleIn = 1 	# /usr/include/X11/Xutil.h:5001
RectanglePart = 2 	# /usr/include/X11/Xutil.h:5002
XVisualInfo = pyglet.gl.glx.XVisualInfo
VisualNoMask = 0 	# /usr/include/X11/Xutil.h:5027
VisualIDMask = 1 	# /usr/include/X11/Xutil.h:5028
VisualScreenMask = 2 	# /usr/include/X11/Xutil.h:5029
VisualDepthMask = 4 	# /usr/include/X11/Xutil.h:5030
VisualClassMask = 8 	# /usr/include/X11/Xutil.h:5031
VisualRedMaskMask = 16 	# /usr/include/X11/Xutil.h:5032
VisualGreenMaskMask = 32 	# /usr/include/X11/Xutil.h:5033
VisualBlueMaskMask = 64 	# /usr/include/X11/Xutil.h:5034
VisualColormapSizeMask = 128 	# /usr/include/X11/Xutil.h:5035
VisualBitsPerRGBMask = 256 	# /usr/include/X11/Xutil.h:5036
VisualAllMask = 1 	# /usr/include/X11/Xutil.h:5037
class struct_anon_95(Structure):
    __slots__ = [
        'colormap',
        'red_max',
        'red_mult',
        'green_max',
        'green_mult',
        'blue_max',
        'blue_mult',
        'base_pixel',
        'visualid',
        'killid',
    ]
struct_anon_95._fields_ = [
    ('colormap', Colormap),
    ('red_max', c_ulong),
    ('red_mult', c_ulong),
    ('green_max', c_ulong),
    ('green_mult', c_ulong),
    ('blue_max', c_ulong),
    ('blue_mult', c_ulong),
    ('base_pixel', c_ulong),
    ('visualid', VisualID),
    ('killid', XID),
]

XStandardColormap = struct_anon_95 	# /usr/include/X11/Xutil.h:5054
BitmapSuccess = 0 	# /usr/include/X11/Xutil.h:5062
BitmapOpenFailed = 1 	# /usr/include/X11/Xutil.h:5063
BitmapFileInvalid = 2 	# /usr/include/X11/Xutil.h:5064
BitmapNoMemory = 3 	# /usr/include/X11/Xutil.h:5065
XCSUCCESS = 0 	# /usr/include/X11/Xutil.h:5076
XCNOMEM = 1 	# /usr/include/X11/Xutil.h:5077
XCNOENT = 2 	# /usr/include/X11/Xutil.h:5078
XContext = c_int 	# /usr/include/X11/Xutil.h:5080
# /usr/include/X11/Xutil.h:5089
XAllocClassHint = _lib.XAllocClassHint
XAllocClassHint.restype = POINTER(XClassHint)
XAllocClassHint.argtypes = []

# /usr/include/X11/Xutil.h:5093
XAllocIconSize = _lib.XAllocIconSize
XAllocIconSize.restype = POINTER(XIconSize)
XAllocIconSize.argtypes = []

# /usr/include/X11/Xutil.h:5097
XAllocSizeHints = _lib.XAllocSizeHints
XAllocSizeHints.restype = POINTER(XSizeHints)
XAllocSizeHints.argtypes = []

# /usr/include/X11/Xutil.h:5101
XAllocStandardColormap = _lib.XAllocStandardColormap
XAllocStandardColormap.restype = POINTER(XStandardColormap)
XAllocStandardColormap.argtypes = []

# /usr/include/X11/Xutil.h:5105
XAllocWMHints = _lib.XAllocWMHints
XAllocWMHints.restype = POINTER(XWMHints)
XAllocWMHints.argtypes = []

# /usr/include/X11/Xutil.h:5109
XClipBox = _lib.XClipBox
XClipBox.restype = c_int
XClipBox.argtypes = [Region, POINTER(XRectangle)]

# /usr/include/X11/Xutil.h:5114
XCreateRegion = _lib.XCreateRegion
XCreateRegion.restype = Region
XCreateRegion.argtypes = []

# /usr/include/X11/Xutil.h:5118
XDefaultString = _lib.XDefaultString
XDefaultString.restype = c_char_p
XDefaultString.argtypes = []

# /usr/include/X11/Xutil.h:5120
XDeleteContext = _lib.XDeleteContext
XDeleteContext.restype = c_int
XDeleteContext.argtypes = [POINTER(Display), XID, XContext]

# /usr/include/X11/Xutil.h:5126
XDestroyRegion = _lib.XDestroyRegion
XDestroyRegion.restype = c_int
XDestroyRegion.argtypes = [Region]

# /usr/include/X11/Xutil.h:5130
XEmptyRegion = _lib.XEmptyRegion
XEmptyRegion.restype = c_int
XEmptyRegion.argtypes = [Region]

# /usr/include/X11/Xutil.h:5134
XEqualRegion = _lib.XEqualRegion
XEqualRegion.restype = c_int
XEqualRegion.argtypes = [Region, Region]

# /usr/include/X11/Xutil.h:5139
XFindContext = _lib.XFindContext
XFindContext.restype = c_int
XFindContext.argtypes = [POINTER(Display), XID, XContext, POINTER(XPointer)]

# /usr/include/X11/Xutil.h:5146
XGetClassHint = _lib.XGetClassHint
XGetClassHint.restype = c_int
XGetClassHint.argtypes = [POINTER(Display), Window, POINTER(XClassHint)]

# /usr/include/X11/Xutil.h:5152
XGetIconSizes = _lib.XGetIconSizes
XGetIconSizes.restype = c_int
XGetIconSizes.argtypes = [POINTER(Display), Window, POINTER(POINTER(XIconSize)), POINTER(c_int)]

# /usr/include/X11/Xutil.h:5159
XGetNormalHints = _lib.XGetNormalHints
XGetNormalHints.restype = c_int
XGetNormalHints.argtypes = [POINTER(Display), Window, POINTER(XSizeHints)]

# /usr/include/X11/Xutil.h:5165
XGetRGBColormaps = _lib.XGetRGBColormaps
XGetRGBColormaps.restype = c_int
XGetRGBColormaps.argtypes = [POINTER(Display), Window, POINTER(POINTER(XStandardColormap)), POINTER(c_int), Atom]

# /usr/include/X11/Xutil.h:5173
XGetSizeHints = _lib.XGetSizeHints
XGetSizeHints.restype = c_int
XGetSizeHints.argtypes = [POINTER(Display), Window, POINTER(XSizeHints), Atom]

# /usr/include/X11/Xutil.h:5180
XGetStandardColormap = _lib.XGetStandardColormap
XGetStandardColormap.restype = c_int
XGetStandardColormap.argtypes = [POINTER(Display), Window, POINTER(XStandardColormap), Atom]

# /usr/include/X11/Xutil.h:5187
XGetTextProperty = _lib.XGetTextProperty
XGetTextProperty.restype = c_int
XGetTextProperty.argtypes = [POINTER(Display), Window, POINTER(XTextProperty), Atom]

# /usr/include/X11/Xutil.h:5194
XGetVisualInfo = _lib.XGetVisualInfo
XGetVisualInfo.restype = POINTER(XVisualInfo)
XGetVisualInfo.argtypes = [POINTER(Display), c_long, POINTER(XVisualInfo), POINTER(c_int)]

# /usr/include/X11/Xutil.h:5201
XGetWMClientMachine = _lib.XGetWMClientMachine
XGetWMClientMachine.restype = c_int
XGetWMClientMachine.argtypes = [POINTER(Display), Window, POINTER(XTextProperty)]

# /usr/include/X11/Xutil.h:5207
XGetWMHints = _lib.XGetWMHints
XGetWMHints.restype = POINTER(XWMHints)
XGetWMHints.argtypes = [POINTER(Display), Window]

# /usr/include/X11/Xutil.h:5212
XGetWMIconName = _lib.XGetWMIconName
XGetWMIconName.restype = c_int
XGetWMIconName.argtypes = [POINTER(Display), Window, POINTER(XTextProperty)]

# /usr/include/X11/Xutil.h:5218
XGetWMName = _lib.XGetWMName
XGetWMName.restype = c_int
XGetWMName.argtypes = [POINTER(Display), Window, POINTER(XTextProperty)]

# /usr/include/X11/Xutil.h:5224
XGetWMNormalHints = _lib.XGetWMNormalHints
XGetWMNormalHints.restype = c_int
XGetWMNormalHints.argtypes = [POINTER(Display), Window, POINTER(XSizeHints), POINTER(c_long)]

# /usr/include/X11/Xutil.h:5231
XGetWMSizeHints = _lib.XGetWMSizeHints
XGetWMSizeHints.restype = c_int
XGetWMSizeHints.argtypes = [POINTER(Display), Window, POINTER(XSizeHints), POINTER(c_long), Atom]

# /usr/include/X11/Xutil.h:5239
XGetZoomHints = _lib.XGetZoomHints
XGetZoomHints.restype = c_int
XGetZoomHints.argtypes = [POINTER(Display), Window, POINTER(XSizeHints)]

# /usr/include/X11/Xutil.h:5245
XIntersectRegion = _lib.XIntersectRegion
XIntersectRegion.restype = c_int
XIntersectRegion.argtypes = [Region, Region, Region]

# /usr/include/X11/Xutil.h:5251
XConvertCase = _lib.XConvertCase
XConvertCase.restype = None
XConvertCase.argtypes = [KeySym, POINTER(KeySym), POINTER(KeySym)]

# /usr/include/X11/Xutil.h:5257
XLookupString = _lib.XLookupString
XLookupString.restype = c_int
XLookupString.argtypes = [POINTER(XKeyEvent), c_char_p, c_int, POINTER(KeySym), POINTER(XComposeStatus)]

# /usr/include/X11/Xutil.h:5265
XMatchVisualInfo = _lib.XMatchVisualInfo
XMatchVisualInfo.restype = c_int
XMatchVisualInfo.argtypes = [POINTER(Display), c_int, c_int, c_int, POINTER(XVisualInfo)]

# /usr/include/X11/Xutil.h:5273
XOffsetRegion = _lib.XOffsetRegion
XOffsetRegion.restype = c_int
XOffsetRegion.argtypes = [Region, c_int, c_int]

# /usr/include/X11/Xutil.h:5279
XPointInRegion = _lib.XPointInRegion
XPointInRegion.restype = c_int
XPointInRegion.argtypes = [Region, c_int, c_int]

# /usr/include/X11/Xutil.h:5285
XPolygonRegion = _lib.XPolygonRegion
XPolygonRegion.restype = Region
XPolygonRegion.argtypes = [POINTER(XPoint), c_int, c_int]

# /usr/include/X11/Xutil.h:5291
XRectInRegion = _lib.XRectInRegion
XRectInRegion.restype = c_int
XRectInRegion.argtypes = [Region, c_int, c_int, c_uint, c_uint]

# /usr/include/X11/Xutil.h:5299
XSaveContext = _lib.XSaveContext
XSaveContext.restype = c_int
XSaveContext.argtypes = [POINTER(Display), XID, XContext, c_char_p]

# /usr/include/X11/Xutil.h:5306
XSetClassHint = _lib.XSetClassHint
XSetClassHint.restype = c_int
XSetClassHint.argtypes = [POINTER(Display), Window, POINTER(XClassHint)]

# /usr/include/X11/Xutil.h:5312
XSetIconSizes = _lib.XSetIconSizes
XSetIconSizes.restype = c_int
XSetIconSizes.argtypes = [POINTER(Display), Window, POINTER(XIconSize), c_int]

# /usr/include/X11/Xutil.h:5319
XSetNormalHints = _lib.XSetNormalHints
XSetNormalHints.restype = c_int
XSetNormalHints.argtypes = [POINTER(Display), Window, POINTER(XSizeHints)]

# /usr/include/X11/Xutil.h:5325
XSetRGBColormaps = _lib.XSetRGBColormaps
XSetRGBColormaps.restype = None
XSetRGBColormaps.argtypes = [POINTER(Display), Window, POINTER(XStandardColormap), c_int, Atom]

# /usr/include/X11/Xutil.h:5333
XSetSizeHints = _lib.XSetSizeHints
XSetSizeHints.restype = c_int
XSetSizeHints.argtypes = [POINTER(Display), Window, POINTER(XSizeHints), Atom]

# /usr/include/X11/Xutil.h:5340
XSetStandardProperties = _lib.XSetStandardProperties
XSetStandardProperties.restype = c_int
XSetStandardProperties.argtypes = [POINTER(Display), Window, c_char_p, c_char_p, Pixmap, POINTER(c_char_p), c_int, POINTER(XSizeHints)]

# /usr/include/X11/Xutil.h:5351
XSetTextProperty = _lib.XSetTextProperty
XSetTextProperty.restype = None
XSetTextProperty.argtypes = [POINTER(Display), Window, POINTER(XTextProperty), Atom]

# /usr/include/X11/Xutil.h:5358
XSetWMClientMachine = _lib.XSetWMClientMachine
XSetWMClientMachine.restype = None
XSetWMClientMachine.argtypes = [POINTER(Display), Window, POINTER(XTextProperty)]

# /usr/include/X11/Xutil.h:5364
XSetWMHints = _lib.XSetWMHints
XSetWMHints.restype = c_int
XSetWMHints.argtypes = [POINTER(Display), Window, POINTER(XWMHints)]

# /usr/include/X11/Xutil.h:5370
XSetWMIconName = _lib.XSetWMIconName
XSetWMIconName.restype = None
XSetWMIconName.argtypes = [POINTER(Display), Window, POINTER(XTextProperty)]

# /usr/include/X11/Xutil.h:5376
XSetWMName = _lib.XSetWMName
XSetWMName.restype = None
XSetWMName.argtypes = [POINTER(Display), Window, POINTER(XTextProperty)]

# /usr/include/X11/Xutil.h:5382
XSetWMNormalHints = _lib.XSetWMNormalHints
XSetWMNormalHints.restype = None
XSetWMNormalHints.argtypes = [POINTER(Display), Window, POINTER(XSizeHints)]

# /usr/include/X11/Xutil.h:5388
XSetWMProperties = _lib.XSetWMProperties
XSetWMProperties.restype = None
XSetWMProperties.argtypes = [POINTER(Display), Window, POINTER(XTextProperty), POINTER(XTextProperty), POINTER(c_char_p), c_int, POINTER(XSizeHints), POINTER(XWMHints), POINTER(XClassHint)]

# /usr/include/X11/Xutil.h:5400
XmbSetWMProperties = _lib.XmbSetWMProperties
XmbSetWMProperties.restype = None
XmbSetWMProperties.argtypes = [POINTER(Display), Window, c_char_p, c_char_p, POINTER(c_char_p), c_int, POINTER(XSizeHints), POINTER(XWMHints), POINTER(XClassHint)]

# /usr/include/X11/Xutil.h:5412
Xutf8SetWMProperties = _lib.Xutf8SetWMProperties
Xutf8SetWMProperties.restype = None
Xutf8SetWMProperties.argtypes = [POINTER(Display), Window, c_char_p, c_char_p, POINTER(c_char_p), c_int, POINTER(XSizeHints), POINTER(XWMHints), POINTER(XClassHint)]

# /usr/include/X11/Xutil.h:5424
XSetWMSizeHints = _lib.XSetWMSizeHints
XSetWMSizeHints.restype = None
XSetWMSizeHints.argtypes = [POINTER(Display), Window, POINTER(XSizeHints), Atom]

# /usr/include/X11/Xutil.h:5431
XSetRegion = _lib.XSetRegion
XSetRegion.restype = c_int
XSetRegion.argtypes = [POINTER(Display), GC, Region]

# /usr/include/X11/Xutil.h:5437
XSetStandardColormap = _lib.XSetStandardColormap
XSetStandardColormap.restype = None
XSetStandardColormap.argtypes = [POINTER(Display), Window, POINTER(XStandardColormap), Atom]

# /usr/include/X11/Xutil.h:5444
XSetZoomHints = _lib.XSetZoomHints
XSetZoomHints.restype = c_int
XSetZoomHints.argtypes = [POINTER(Display), Window, POINTER(XSizeHints)]

# /usr/include/X11/Xutil.h:5450
XShrinkRegion = _lib.XShrinkRegion
XShrinkRegion.restype = c_int
XShrinkRegion.argtypes = [Region, c_int, c_int]

# /usr/include/X11/Xutil.h:5456
XStringListToTextProperty = _lib.XStringListToTextProperty
XStringListToTextProperty.restype = c_int
XStringListToTextProperty.argtypes = [POINTER(c_char_p), c_int, POINTER(XTextProperty)]

# /usr/include/X11/Xutil.h:5462
XSubtractRegion = _lib.XSubtractRegion
XSubtractRegion.restype = c_int
XSubtractRegion.argtypes = [Region, Region, Region]

# /usr/include/X11/Xutil.h:5468
XmbTextListToTextProperty = _lib.XmbTextListToTextProperty
XmbTextListToTextProperty.restype = c_int
XmbTextListToTextProperty.argtypes = [POINTER(Display), POINTER(c_char_p), c_int, XICCEncodingStyle, POINTER(XTextProperty)]

# /usr/include/X11/Xutil.h:5476
XwcTextListToTextProperty = _lib.XwcTextListToTextProperty
XwcTextListToTextProperty.restype = c_int
XwcTextListToTextProperty.argtypes = [POINTER(Display), POINTER(c_wchar_p), c_int, XICCEncodingStyle, POINTER(XTextProperty)]

# /usr/include/X11/Xutil.h:5484
Xutf8TextListToTextProperty = _lib.Xutf8TextListToTextProperty
Xutf8TextListToTextProperty.restype = c_int
Xutf8TextListToTextProperty.argtypes = [POINTER(Display), POINTER(c_char_p), c_int, XICCEncodingStyle, POINTER(XTextProperty)]

# /usr/include/X11/Xutil.h:5492
XwcFreeStringList = _lib.XwcFreeStringList
XwcFreeStringList.restype = None
XwcFreeStringList.argtypes = [POINTER(c_wchar_p)]

# /usr/include/X11/Xutil.h:5496
XTextPropertyToStringList = _lib.XTextPropertyToStringList
XTextPropertyToStringList.restype = c_int
XTextPropertyToStringList.argtypes = [POINTER(XTextProperty), POINTER(POINTER(c_char_p)), POINTER(c_int)]

# /usr/include/X11/Xutil.h:5502
XmbTextPropertyToTextList = _lib.XmbTextPropertyToTextList
XmbTextPropertyToTextList.restype = c_int
XmbTextPropertyToTextList.argtypes = [POINTER(Display), POINTER(XTextProperty), POINTER(POINTER(c_char_p)), POINTER(c_int)]

# /usr/include/X11/Xutil.h:5509
XwcTextPropertyToTextList = _lib.XwcTextPropertyToTextList
XwcTextPropertyToTextList.restype = c_int
XwcTextPropertyToTextList.argtypes = [POINTER(Display), POINTER(XTextProperty), POINTER(POINTER(c_wchar_p)), POINTER(c_int)]

# /usr/include/X11/Xutil.h:5516
Xutf8TextPropertyToTextList = _lib.Xutf8TextPropertyToTextList
Xutf8TextPropertyToTextList.restype = c_int
Xutf8TextPropertyToTextList.argtypes = [POINTER(Display), POINTER(XTextProperty), POINTER(POINTER(c_char_p)), POINTER(c_int)]

# /usr/include/X11/Xutil.h:5523
XUnionRectWithRegion = _lib.XUnionRectWithRegion
XUnionRectWithRegion.restype = c_int
XUnionRectWithRegion.argtypes = [POINTER(XRectangle), Region, Region]

# /usr/include/X11/Xutil.h:5529
XUnionRegion = _lib.XUnionRegion
XUnionRegion.restype = c_int
XUnionRegion.argtypes = [Region, Region, Region]

# /usr/include/X11/Xutil.h:5535
XWMGeometry = _lib.XWMGeometry
XWMGeometry.restype = c_int
XWMGeometry.argtypes = [POINTER(Display), c_int, c_char_p, c_char_p, c_uint, POINTER(XSizeHints), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int)]

# /usr/include/X11/Xutil.h:5549
XXorRegion = _lib.XXorRegion
XXorRegion.restype = c_int
XXorRegion.argtypes = [Region, Region, Region]


__all__ = ['XlibSpecificationRelease', 'X_PROTOCOL', 'X_PROTOCOL_REVISION',
'XID', 'Mask', 'Atom', 'VisualID', 'Time', 'Window', 'Drawable', 'Font',
'Pixmap', 'Cursor', 'Colormap', 'GContext', 'KeySym', 'KeyCode', 'None_',
'ParentRelative', 'CopyFromParent', 'PointerWindow', 'InputFocus',
'PointerRoot', 'AnyPropertyType', 'AnyKey', 'AnyButton', 'AllTemporary',
'CurrentTime', 'NoSymbol', 'NoEventMask', 'KeyPressMask', 'KeyReleaseMask',
'ButtonPressMask', 'ButtonReleaseMask', 'EnterWindowMask', 'LeaveWindowMask',
'PointerMotionMask', 'PointerMotionHintMask', 'Button1MotionMask',
'Button2MotionMask', 'Button3MotionMask', 'Button4MotionMask',
'Button5MotionMask', 'ButtonMotionMask', 'KeymapStateMask', 'ExposureMask',
'VisibilityChangeMask', 'StructureNotifyMask', 'ResizeRedirectMask',
'SubstructureNotifyMask', 'SubstructureRedirectMask', 'FocusChangeMask',
'PropertyChangeMask', 'ColormapChangeMask', 'OwnerGrabButtonMask', 'KeyPress',
'KeyRelease', 'ButtonPress', 'ButtonRelease', 'MotionNotify', 'EnterNotify',
'LeaveNotify', 'FocusIn', 'FocusOut', 'KeymapNotify', 'Expose',
'GraphicsExpose', 'NoExpose', 'VisibilityNotify', 'CreateNotify',
'DestroyNotify', 'UnmapNotify', 'MapNotify', 'MapRequest', 'ReparentNotify',
'ConfigureNotify', 'ConfigureRequest', 'GravityNotify', 'ResizeRequest',
'CirculateNotify', 'CirculateRequest', 'PropertyNotify', 'SelectionClear',
'SelectionRequest', 'SelectionNotify', 'ColormapNotify', 'ClientMessage',
'MappingNotify', 'LASTEvent', 'ShiftMask', 'LockMask', 'ControlMask',
'Mod1Mask', 'Mod2Mask', 'Mod3Mask', 'Mod4Mask', 'Mod5Mask', 'ShiftMapIndex',
'LockMapIndex', 'ControlMapIndex', 'Mod1MapIndex', 'Mod2MapIndex',
'Mod3MapIndex', 'Mod4MapIndex', 'Mod5MapIndex', 'Button1Mask', 'Button2Mask',
'Button3Mask', 'Button4Mask', 'Button5Mask', 'AnyModifier', 'Button1',
'Button2', 'Button3', 'Button4', 'Button5', 'NotifyNormal', 'NotifyGrab',
'NotifyUngrab', 'NotifyWhileGrabbed', 'NotifyHint', 'NotifyAncestor',
'NotifyVirtual', 'NotifyInferior', 'NotifyNonlinear',
'NotifyNonlinearVirtual', 'NotifyPointer', 'NotifyPointerRoot',
'NotifyDetailNone', 'VisibilityUnobscured', 'VisibilityPartiallyObscured',
'VisibilityFullyObscured', 'PlaceOnTop', 'PlaceOnBottom', 'FamilyInternet',
'FamilyDECnet', 'FamilyChaos', 'FamilyInternet6', 'FamilyServerInterpreted',
'PropertyNewValue', 'PropertyDelete', 'ColormapUninstalled',
'ColormapInstalled', 'GrabModeSync', 'GrabModeAsync', 'GrabSuccess',
'AlreadyGrabbed', 'GrabInvalidTime', 'GrabNotViewable', 'GrabFrozen',
'AsyncPointer', 'SyncPointer', 'ReplayPointer', 'AsyncKeyboard',
'SyncKeyboard', 'ReplayKeyboard', 'AsyncBoth', 'SyncBoth', 'RevertToParent',
'Success', 'BadRequest', 'BadValue', 'BadWindow', 'BadPixmap', 'BadAtom',
'BadCursor', 'BadFont', 'BadMatch', 'BadDrawable', 'BadAccess', 'BadAlloc',
'BadColor', 'BadGC', 'BadIDChoice', 'BadName', 'BadLength',
'BadImplementation', 'FirstExtensionError', 'LastExtensionError',
'InputOutput', 'InputOnly', 'CWBackPixmap', 'CWBackPixel', 'CWBorderPixmap',
'CWBorderPixel', 'CWBitGravity', 'CWWinGravity', 'CWBackingStore',
'CWBackingPlanes', 'CWBackingPixel', 'CWOverrideRedirect', 'CWSaveUnder',
'CWEventMask', 'CWDontPropagate', 'CWColormap', 'CWCursor', 'CWX', 'CWY',
'CWWidth', 'CWHeight', 'CWBorderWidth', 'CWSibling', 'CWStackMode',
'ForgetGravity', 'NorthWestGravity', 'NorthGravity', 'NorthEastGravity',
'WestGravity', 'CenterGravity', 'EastGravity', 'SouthWestGravity',
'SouthGravity', 'SouthEastGravity', 'StaticGravity', 'UnmapGravity',
'NotUseful', 'WhenMapped', 'Always', 'IsUnmapped', 'IsUnviewable',
'IsViewable', 'SetModeInsert', 'SetModeDelete', 'DestroyAll',
'RetainPermanent', 'RetainTemporary', 'Above', 'Below', 'TopIf', 'BottomIf',
'Opposite', 'RaiseLowest', 'LowerHighest', 'PropModeReplace',
'PropModePrepend', 'PropModeAppend', 'GXclear', 'GXand', 'GXandReverse',
'GXcopy', 'GXandInverted', 'GXnoop', 'GXxor', 'GXor', 'GXnor', 'GXequiv',
'GXinvert', 'GXorReverse', 'GXcopyInverted', 'GXorInverted', 'GXnand',
'GXset', 'LineSolid', 'LineOnOffDash', 'LineDoubleDash', 'CapNotLast',
'CapButt', 'CapRound', 'CapProjecting', 'JoinMiter', 'JoinRound', 'JoinBevel',
'FillSolid', 'FillTiled', 'FillStippled', 'FillOpaqueStippled', 'EvenOddRule',
'WindingRule', 'ClipByChildren', 'IncludeInferiors', 'Unsorted', 'YSorted',
'YXSorted', 'YXBanded', 'CoordModeOrigin', 'CoordModePrevious', 'Complex',
'Nonconvex', 'Convex', 'ArcChord', 'ArcPieSlice', 'GCFunction', 'GCPlaneMask',
'GCForeground', 'GCBackground', 'GCLineWidth', 'GCLineStyle', 'GCCapStyle',
'GCJoinStyle', 'GCFillStyle', 'GCFillRule', 'GCTile', 'GCStipple',
'GCTileStipXOrigin', 'GCTileStipYOrigin', 'GCFont', 'GCSubwindowMode',
'GCGraphicsExposures', 'GCClipXOrigin', 'GCClipYOrigin', 'GCClipMask',
'GCDashOffset', 'GCDashList', 'GCArcMode', 'GCLastBit', 'FontLeftToRight',
'FontRightToLeft', 'FontChange', 'XYBitmap', 'XYPixmap', 'ZPixmap',
'AllocNone', 'AllocAll', 'DoRed', 'DoGreen', 'DoBlue', 'CursorShape',
'TileShape', 'StippleShape', 'AutoRepeatModeOff', 'AutoRepeatModeOn',
'AutoRepeatModeDefault', 'LedModeOff', 'LedModeOn', 'KBKeyClickPercent',
'KBBellPercent', 'KBBellPitch', 'KBBellDuration', 'KBLed', 'KBLedMode',
'KBKey', 'KBAutoRepeatMode', 'MappingSuccess', 'MappingBusy', 'MappingFailed',
'MappingModifier', 'MappingKeyboard', 'MappingPointer', 'DontPreferBlanking',
'PreferBlanking', 'DefaultBlanking', 'DisableScreenSaver',
'DisableScreenInterval', 'DontAllowExposures', 'AllowExposures',
'DefaultExposures', 'ScreenSaverReset', 'ScreenSaverActive', 'HostInsert',
'HostDelete', 'EnableAccess', 'DisableAccess', 'StaticGray', 'GrayScale',
'StaticColor', 'PseudoColor', 'TrueColor', 'DirectColor', 'LSBFirst',
'MSBFirst', '_Xmblen', 'X_HAVE_UTF8_STRING', 'XPointer', 'Bool', 'Status',
'True', 'False', 'QueuedAlready', 'QueuedAfterReading', 'QueuedAfterFlush',
'XExtData', 'XExtCodes', 'XPixmapFormatValues', 'XGCValues', 'GC', 'Visual',
'Depth', 'Screen', 'ScreenFormat', 'XSetWindowAttributes',
'XWindowAttributes', 'XHostAddress', 'XServerInterpretedAddress', 'XImage',
'XWindowChanges', 'XColor', 'XSegment', 'XPoint', 'XRectangle', 'XArc',
'XKeyboardControl', 'XKeyboardState', 'XTimeCoord', 'XModifierKeymap',
'Display', '_XPrivDisplay', 'XKeyEvent', 'XKeyPressedEvent',
'XKeyReleasedEvent', 'XButtonEvent', 'XButtonPressedEvent',
'XButtonReleasedEvent', 'XMotionEvent', 'XPointerMovedEvent',
'XCrossingEvent', 'XEnterWindowEvent', 'XLeaveWindowEvent',
'XFocusChangeEvent', 'XFocusInEvent', 'XFocusOutEvent', 'XKeymapEvent',
'XExposeEvent', 'XGraphicsExposeEvent', 'XNoExposeEvent', 'XVisibilityEvent',
'XCreateWindowEvent', 'XDestroyWindowEvent', 'XUnmapEvent', 'XMapEvent',
'XMapRequestEvent', 'XReparentEvent', 'XConfigureEvent', 'XGravityEvent',
'XResizeRequestEvent', 'XConfigureRequestEvent', 'XCirculateEvent',
'XCirculateRequestEvent', 'XPropertyEvent', 'XSelectionClearEvent',
'XSelectionRequestEvent', 'XSelectionEvent', 'XColormapEvent',
'XClientMessageEvent', 'XMappingEvent', 'XErrorEvent', 'XAnyEvent', 'XEvent',
'XCharStruct', 'XFontProp', 'XFontStruct', 'XTextItem', 'XChar2b',
'XTextItem16', 'XEDataObject', 'XFontSetExtents', 'XOM', 'XOC', 'XFontSet',
'XmbTextItem', 'XwcTextItem', 'XOMCharSetList', 'XOrientation',
'XOMOrientation_LTR_TTB', 'XOMOrientation_RTL_TTB', 'XOMOrientation_TTB_LTR',
'XOMOrientation_TTB_RTL', 'XOMOrientation_Context', 'XOMOrientation',
'XOMFontInfo', 'XIM', 'XIC', 'XIMProc', 'XICProc', 'XIDProc', 'XIMStyle',
'XIMStyles', 'XIMPreeditArea', 'XIMPreeditCallbacks', 'XIMPreeditPosition',
'XIMPreeditNothing', 'XIMPreeditNone', 'XIMStatusArea', 'XIMStatusCallbacks',
'XIMStatusNothing', 'XIMStatusNone', 'XBufferOverflow', 'XLookupNone',
'XLookupChars', 'XLookupKeySym', 'XLookupBoth', 'XVaNestedList',
'XIMCallback', 'XICCallback', 'XIMFeedback', 'XIMReverse', 'XIMUnderline',
'XIMHighlight', 'XIMPrimary', 'XIMSecondary', 'XIMTertiary',
'XIMVisibleToForward', 'XIMVisibleToBackword', 'XIMVisibleToCenter',
'XIMText', 'XIMPreeditState', 'XIMPreeditUnKnown', 'XIMPreeditEnable',
'XIMPreeditDisable', 'XIMPreeditStateNotifyCallbackStruct', 'XIMResetState',
'XIMInitialState', 'XIMPreserveState', 'XIMStringConversionFeedback',
'XIMStringConversionLeftEdge', 'XIMStringConversionRightEdge',
'XIMStringConversionTopEdge', 'XIMStringConversionBottomEdge',
'XIMStringConversionConcealed', 'XIMStringConversionWrapped',
'XIMStringConversionText', 'XIMStringConversionPosition',
'XIMStringConversionType', 'XIMStringConversionBuffer',
'XIMStringConversionLine', 'XIMStringConversionWord',
'XIMStringConversionChar', 'XIMStringConversionOperation',
'XIMStringConversionSubstitution', 'XIMStringConversionRetrieval',
'XIMCaretDirection', 'XIMForwardChar', 'XIMBackwardChar', 'XIMForwardWord',
'XIMBackwardWord', 'XIMCaretUp', 'XIMCaretDown', 'XIMNextLine',
'XIMPreviousLine', 'XIMLineStart', 'XIMLineEnd', 'XIMAbsolutePosition',
'XIMDontChange', 'XIMStringConversionCallbackStruct',
'XIMPreeditDrawCallbackStruct', 'XIMCaretStyle', 'XIMIsInvisible',
'XIMIsPrimary', 'XIMIsSecondary', 'XIMPreeditCaretCallbackStruct',
'XIMStatusDataType', 'XIMTextType', 'XIMBitmapType',
'XIMStatusDrawCallbackStruct', 'XIMHotKeyTrigger', 'XIMHotKeyTriggers',
'XIMHotKeyState', 'XIMHotKeyStateON', 'XIMHotKeyStateOFF', 'XIMValuesList',
'XLoadQueryFont', 'XQueryFont', 'XGetMotionEvents', 'XDeleteModifiermapEntry',
'XGetModifierMapping', 'XInsertModifiermapEntry', 'XNewModifiermap',
'XCreateImage', 'XInitImage', 'XGetImage', 'XGetSubImage', 'XOpenDisplay',
'XrmInitialize', 'XFetchBytes', 'XFetchBuffer', 'XGetAtomName',
'XGetAtomNames', 'XGetDefault', 'XDisplayName', 'XKeysymToString',
'XSynchronize', 'XSetAfterFunction', 'XInternAtom', 'XInternAtoms',
'XCopyColormapAndFree', 'XCreateColormap', 'XCreatePixmapCursor',
'XCreateGlyphCursor', 'XCreateFontCursor', 'XLoadFont', 'XCreateGC',
'XGContextFromGC', 'XFlushGC', 'XCreatePixmap', 'XCreateBitmapFromData',
'XCreatePixmapFromBitmapData', 'XCreateSimpleWindow', 'XGetSelectionOwner',
'XCreateWindow', 'XListInstalledColormaps', 'XListFonts',
'XListFontsWithInfo', 'XGetFontPath', 'XListExtensions', 'XListProperties',
'XListHosts', 'XKeycodeToKeysym', 'XLookupKeysym', 'XGetKeyboardMapping',
'XStringToKeysym', 'XMaxRequestSize', 'XExtendedMaxRequestSize',
'XResourceManagerString', 'XScreenResourceString', 'XDisplayMotionBufferSize',
'XVisualIDFromVisual', 'XInitThreads', 'XLockDisplay', 'XUnlockDisplay',
'XInitExtension', 'XAddExtension', 'XFindOnExtensionList',
'XEHeadOfExtensionList', 'XRootWindow', 'XDefaultRootWindow',
'XRootWindowOfScreen', 'XDefaultVisual', 'XDefaultVisualOfScreen',
'XDefaultGC', 'XDefaultGCOfScreen', 'XBlackPixel', 'XWhitePixel',
'XAllPlanes', 'XBlackPixelOfScreen', 'XWhitePixelOfScreen', 'XNextRequest',
'XLastKnownRequestProcessed', 'XServerVendor', 'XDisplayString',
'XDefaultColormap', 'XDefaultColormapOfScreen', 'XDisplayOfScreen',
'XScreenOfDisplay', 'XDefaultScreenOfDisplay', 'XEventMaskOfScreen',
'XScreenNumberOfScreen', 'XErrorHandler', 'XSetErrorHandler',
'XIOErrorHandler', 'XSetIOErrorHandler', 'XListPixmapFormats', 'XListDepths',
'XReconfigureWMWindow', 'XGetWMProtocols', 'XSetWMProtocols',
'XIconifyWindow', 'XWithdrawWindow', 'XGetCommand', 'XGetWMColormapWindows',
'XSetWMColormapWindows', 'XFreeStringList', 'XSetTransientForHint',
'XActivateScreenSaver', 'XAddHost', 'XAddHosts', 'XAddToExtensionList',
'XAddToSaveSet', 'XAllocColor', 'XAllocColorCells', 'XAllocColorPlanes',
'XAllocNamedColor', 'XAllowEvents', 'XAutoRepeatOff', 'XAutoRepeatOn',
'XBell', 'XBitmapBitOrder', 'XBitmapPad', 'XBitmapUnit', 'XCellsOfScreen',
'XChangeActivePointerGrab', 'XChangeGC', 'XChangeKeyboardControl',
'XChangeKeyboardMapping', 'XChangePointerControl', 'XChangeProperty',
'XChangeSaveSet', 'XChangeWindowAttributes', 'XCheckIfEvent',
'XCheckMaskEvent', 'XCheckTypedEvent', 'XCheckTypedWindowEvent',
'XCheckWindowEvent', 'XCirculateSubwindows', 'XCirculateSubwindowsDown',
'XCirculateSubwindowsUp', 'XClearArea', 'XClearWindow', 'XCloseDisplay',
'XConfigureWindow', 'XConnectionNumber', 'XConvertSelection', 'XCopyArea',
'XCopyGC', 'XCopyPlane', 'XDefaultDepth', 'XDefaultDepthOfScreen',
'XDefaultScreen', 'XDefineCursor', 'XDeleteProperty', 'XDestroyWindow',
'XDestroySubwindows', 'XDoesBackingStore', 'XDoesSaveUnders',
'XDisableAccessControl', 'XDisplayCells', 'XDisplayHeight',
'XDisplayHeightMM', 'XDisplayKeycodes', 'XDisplayPlanes', 'XDisplayWidth',
'XDisplayWidthMM', 'XDrawArc', 'XDrawArcs', 'XDrawImageString',
'XDrawImageString16', 'XDrawLine', 'XDrawLines', 'XDrawPoint', 'XDrawPoints',
'XDrawRectangle', 'XDrawRectangles', 'XDrawSegments', 'XDrawString',
'XDrawString16', 'XDrawText', 'XDrawText16', 'XEnableAccessControl',
'XEventsQueued', 'XFetchName', 'XFillArc', 'XFillArcs', 'XFillPolygon',
'XFillRectangle', 'XFillRectangles', 'XFlush', 'XForceScreenSaver', 'XFree',
'XFreeColormap', 'XFreeColors', 'XFreeCursor', 'XFreeExtensionList',
'XFreeFont', 'XFreeFontInfo', 'XFreeFontNames', 'XFreeFontPath', 'XFreeGC',
'XFreeModifiermap', 'XFreePixmap', 'XGeometry', 'XGetErrorDatabaseText',
'XGetErrorText', 'XGetFontProperty', 'XGetGCValues', 'XGetGeometry',
'XGetIconName', 'XGetInputFocus', 'XGetKeyboardControl', 'XGetPointerControl',
'XGetPointerMapping', 'XGetScreenSaver', 'XGetTransientForHint',
'XGetWindowProperty', 'XGetWindowAttributes', 'XGrabButton', 'XGrabKey',
'XGrabKeyboard', 'XGrabPointer', 'XGrabServer', 'XHeightMMOfScreen',
'XHeightOfScreen', 'XIfEvent', 'XImageByteOrder', 'XInstallColormap',
'XKeysymToKeycode', 'XKillClient', 'XLookupColor', 'XLowerWindow',
'XMapRaised', 'XMapSubwindows', 'XMapWindow', 'XMaskEvent',
'XMaxCmapsOfScreen', 'XMinCmapsOfScreen', 'XMoveResizeWindow', 'XMoveWindow',
'XNextEvent', 'XNoOp', 'XParseColor', 'XParseGeometry', 'XPeekEvent',
'XPeekIfEvent', 'XPending', 'XPlanesOfScreen', 'XProtocolRevision',
'XProtocolVersion', 'XPutBackEvent', 'XPutImage', 'XQLength',
'XQueryBestCursor', 'XQueryBestSize', 'XQueryBestStipple', 'XQueryBestTile',
'XQueryColor', 'XQueryColors', 'XQueryExtension', 'XQueryKeymap',
'XQueryPointer', 'XQueryTextExtents', 'XQueryTextExtents16', 'XQueryTree',
'XRaiseWindow', 'XReadBitmapFile', 'XReadBitmapFileData', 'XRebindKeysym',
'XRecolorCursor', 'XRefreshKeyboardMapping', 'XRemoveFromSaveSet',
'XRemoveHost', 'XRemoveHosts', 'XReparentWindow', 'XResetScreenSaver',
'XResizeWindow', 'XRestackWindows', 'XRotateBuffers',
'XRotateWindowProperties', 'XScreenCount', 'XSelectInput', 'XSendEvent',
'XSetAccessControl', 'XSetArcMode', 'XSetBackground', 'XSetClipMask',
'XSetClipOrigin', 'XSetClipRectangles', 'XSetCloseDownMode', 'XSetCommand',
'XSetDashes', 'XSetFillRule', 'XSetFillStyle', 'XSetFont', 'XSetFontPath',
'XSetForeground', 'XSetFunction', 'XSetGraphicsExposures', 'XSetIconName',
'XSetInputFocus', 'XSetLineAttributes', 'XSetModifierMapping',
'XSetPlaneMask', 'XSetPointerMapping', 'XSetScreenSaver',
'XSetSelectionOwner', 'XSetState', 'XSetStipple', 'XSetSubwindowMode',
'XSetTSOrigin', 'XSetTile', 'XSetWindowBackground',
'XSetWindowBackgroundPixmap', 'XSetWindowBorder', 'XSetWindowBorderPixmap',
'XSetWindowBorderWidth', 'XSetWindowColormap', 'XStoreBuffer', 'XStoreBytes',
'XStoreColor', 'XStoreColors', 'XStoreName', 'XStoreNamedColor', 'XSync',
'XTextExtents', 'XTextExtents16', 'XTextWidth', 'XTextWidth16',
'XTranslateCoordinates', 'XUndefineCursor', 'XUngrabButton', 'XUngrabKey',
'XUngrabKeyboard', 'XUngrabPointer', 'XUngrabServer', 'XUninstallColormap',
'XUnloadFont', 'XUnmapSubwindows', 'XUnmapWindow', 'XVendorRelease',
'XWarpPointer', 'XWidthMMOfScreen', 'XWidthOfScreen', 'XWindowEvent',
'XWriteBitmapFile', 'XSupportsLocale', 'XSetLocaleModifiers', 'XOpenOM',
'XCloseOM', 'XSetOMValues', 'XGetOMValues', 'XDisplayOfOM', 'XLocaleOfOM',
'XCreateOC', 'XDestroyOC', 'XOMOfOC', 'XSetOCValues', 'XGetOCValues',
'XCreateFontSet', 'XFreeFontSet', 'XFontsOfFontSet',
'XBaseFontNameListOfFontSet', 'XLocaleOfFontSet', 'XContextDependentDrawing',
'XDirectionalDependentDrawing', 'XContextualDrawing', 'XExtentsOfFontSet',
'XmbTextEscapement', 'XwcTextEscapement', 'Xutf8TextEscapement',
'XmbTextExtents', 'XwcTextExtents', 'Xutf8TextExtents',
'XmbTextPerCharExtents', 'XwcTextPerCharExtents', 'Xutf8TextPerCharExtents',
'XmbDrawText', 'XwcDrawText', 'Xutf8DrawText', 'XmbDrawString',
'XwcDrawString', 'Xutf8DrawString', 'XmbDrawImageString',
'XwcDrawImageString', 'Xutf8DrawImageString', 'XOpenIM', 'XCloseIM',
'XGetIMValues', 'XSetIMValues', 'XDisplayOfIM', 'XLocaleOfIM', 'XCreateIC',
'XDestroyIC', 'XSetICFocus', 'XUnsetICFocus', 'XwcResetIC', 'XmbResetIC',
'Xutf8ResetIC', 'XSetICValues', 'XGetICValues', 'XIMOfIC', 'XFilterEvent',
'XmbLookupString', 'XwcLookupString', 'Xutf8LookupString',
'XVaCreateNestedList', 'XRegisterIMInstantiateCallback',
'XUnregisterIMInstantiateCallback', 'XConnectionWatchProc',
'XInternalConnectionNumbers', 'XProcessInternalConnection',
'XAddConnectionWatch', 'XRemoveConnectionWatch', 'XSetAuthorization',
'_Xmbtowc', '_Xwctomb', 'X_PROTOCOL', 'X_PROTOCOL_REVISION', 'None_',
'ParentRelative', 'CopyFromParent', 'PointerWindow', 'InputFocus',
'PointerRoot', 'AnyPropertyType', 'AnyKey', 'AnyButton', 'AllTemporary',
'CurrentTime', 'NoSymbol', 'NoEventMask', 'KeyPressMask', 'KeyReleaseMask',
'ButtonPressMask', 'ButtonReleaseMask', 'EnterWindowMask', 'LeaveWindowMask',
'PointerMotionMask', 'PointerMotionHintMask', 'Button1MotionMask',
'Button2MotionMask', 'Button3MotionMask', 'Button4MotionMask',
'Button5MotionMask', 'ButtonMotionMask', 'KeymapStateMask', 'ExposureMask',
'VisibilityChangeMask', 'StructureNotifyMask', 'ResizeRedirectMask',
'SubstructureNotifyMask', 'SubstructureRedirectMask', 'FocusChangeMask',
'PropertyChangeMask', 'ColormapChangeMask', 'OwnerGrabButtonMask', 'KeyPress',
'KeyRelease', 'ButtonPress', 'ButtonRelease', 'MotionNotify', 'EnterNotify',
'LeaveNotify', 'FocusIn', 'FocusOut', 'KeymapNotify', 'Expose',
'GraphicsExpose', 'NoExpose', 'VisibilityNotify', 'CreateNotify',
'DestroyNotify', 'UnmapNotify', 'MapNotify', 'MapRequest', 'ReparentNotify',
'ConfigureNotify', 'ConfigureRequest', 'GravityNotify', 'ResizeRequest',
'CirculateNotify', 'CirculateRequest', 'PropertyNotify', 'SelectionClear',
'SelectionRequest', 'SelectionNotify', 'ColormapNotify', 'ClientMessage',
'MappingNotify', 'LASTEvent', 'ShiftMask', 'LockMask', 'ControlMask',
'Mod1Mask', 'Mod2Mask', 'Mod3Mask', 'Mod4Mask', 'Mod5Mask', 'ShiftMapIndex',
'LockMapIndex', 'ControlMapIndex', 'Mod1MapIndex', 'Mod2MapIndex',
'Mod3MapIndex', 'Mod4MapIndex', 'Mod5MapIndex', 'Button1Mask', 'Button2Mask',
'Button3Mask', 'Button4Mask', 'Button5Mask', 'AnyModifier', 'Button1',
'Button2', 'Button3', 'Button4', 'Button5', 'NotifyNormal', 'NotifyGrab',
'NotifyUngrab', 'NotifyWhileGrabbed', 'NotifyHint', 'NotifyAncestor',
'NotifyVirtual', 'NotifyInferior', 'NotifyNonlinear',
'NotifyNonlinearVirtual', 'NotifyPointer', 'NotifyPointerRoot',
'NotifyDetailNone', 'VisibilityUnobscured', 'VisibilityPartiallyObscured',
'VisibilityFullyObscured', 'PlaceOnTop', 'PlaceOnBottom', 'FamilyInternet',
'FamilyDECnet', 'FamilyChaos', 'FamilyInternet6', 'FamilyServerInterpreted',
'PropertyNewValue', 'PropertyDelete', 'ColormapUninstalled',
'ColormapInstalled', 'GrabModeSync', 'GrabModeAsync', 'GrabSuccess',
'AlreadyGrabbed', 'GrabInvalidTime', 'GrabNotViewable', 'GrabFrozen',
'AsyncPointer', 'SyncPointer', 'ReplayPointer', 'AsyncKeyboard',
'SyncKeyboard', 'ReplayKeyboard', 'AsyncBoth', 'SyncBoth', 'RevertToParent',
'Success', 'BadRequest', 'BadValue', 'BadWindow', 'BadPixmap', 'BadAtom',
'BadCursor', 'BadFont', 'BadMatch', 'BadDrawable', 'BadAccess', 'BadAlloc',
'BadColor', 'BadGC', 'BadIDChoice', 'BadName', 'BadLength',
'BadImplementation', 'FirstExtensionError', 'LastExtensionError',
'InputOutput', 'InputOnly', 'CWBackPixmap', 'CWBackPixel', 'CWBorderPixmap',
'CWBorderPixel', 'CWBitGravity', 'CWWinGravity', 'CWBackingStore',
'CWBackingPlanes', 'CWBackingPixel', 'CWOverrideRedirect', 'CWSaveUnder',
'CWEventMask', 'CWDontPropagate', 'CWColormap', 'CWCursor', 'CWX', 'CWY',
'CWWidth', 'CWHeight', 'CWBorderWidth', 'CWSibling', 'CWStackMode',
'ForgetGravity', 'NorthWestGravity', 'NorthGravity', 'NorthEastGravity',
'WestGravity', 'CenterGravity', 'EastGravity', 'SouthWestGravity',
'SouthGravity', 'SouthEastGravity', 'StaticGravity', 'UnmapGravity',
'NotUseful', 'WhenMapped', 'Always', 'IsUnmapped', 'IsUnviewable',
'IsViewable', 'SetModeInsert', 'SetModeDelete', 'DestroyAll',
'RetainPermanent', 'RetainTemporary', 'Above', 'Below', 'TopIf', 'BottomIf',
'Opposite', 'RaiseLowest', 'LowerHighest', 'PropModeReplace',
'PropModePrepend', 'PropModeAppend', 'GXclear', 'GXand', 'GXandReverse',
'GXcopy', 'GXandInverted', 'GXnoop', 'GXxor', 'GXor', 'GXnor', 'GXequiv',
'GXinvert', 'GXorReverse', 'GXcopyInverted', 'GXorInverted', 'GXnand',
'GXset', 'LineSolid', 'LineOnOffDash', 'LineDoubleDash', 'CapNotLast',
'CapButt', 'CapRound', 'CapProjecting', 'JoinMiter', 'JoinRound', 'JoinBevel',
'FillSolid', 'FillTiled', 'FillStippled', 'FillOpaqueStippled', 'EvenOddRule',
'WindingRule', 'ClipByChildren', 'IncludeInferiors', 'Unsorted', 'YSorted',
'YXSorted', 'YXBanded', 'CoordModeOrigin', 'CoordModePrevious', 'Complex',
'Nonconvex', 'Convex', 'ArcChord', 'ArcPieSlice', 'GCFunction', 'GCPlaneMask',
'GCForeground', 'GCBackground', 'GCLineWidth', 'GCLineStyle', 'GCCapStyle',
'GCJoinStyle', 'GCFillStyle', 'GCFillRule', 'GCTile', 'GCStipple',
'GCTileStipXOrigin', 'GCTileStipYOrigin', 'GCFont', 'GCSubwindowMode',
'GCGraphicsExposures', 'GCClipXOrigin', 'GCClipYOrigin', 'GCClipMask',
'GCDashOffset', 'GCDashList', 'GCArcMode', 'GCLastBit', 'FontLeftToRight',
'FontRightToLeft', 'FontChange', 'XYBitmap', 'XYPixmap', 'ZPixmap',
'AllocNone', 'AllocAll', 'DoRed', 'DoGreen', 'DoBlue', 'CursorShape',
'TileShape', 'StippleShape', 'AutoRepeatModeOff', 'AutoRepeatModeOn',
'AutoRepeatModeDefault', 'LedModeOff', 'LedModeOn', 'KBKeyClickPercent',
'KBBellPercent', 'KBBellPitch', 'KBBellDuration', 'KBLed', 'KBLedMode',
'KBKey', 'KBAutoRepeatMode', 'MappingSuccess', 'MappingBusy', 'MappingFailed',
'MappingModifier', 'MappingKeyboard', 'MappingPointer', 'DontPreferBlanking',
'PreferBlanking', 'DefaultBlanking', 'DisableScreenSaver',
'DisableScreenInterval', 'DontAllowExposures', 'AllowExposures',
'DefaultExposures', 'ScreenSaverReset', 'ScreenSaverActive', 'HostInsert',
'HostDelete', 'EnableAccess', 'DisableAccess', 'StaticGray', 'GrayScale',
'StaticColor', 'PseudoColor', 'TrueColor', 'DirectColor', 'LSBFirst',
'MSBFirst', 'XlibSpecificationRelease', 'X_PROTOCOL', 'X_PROTOCOL_REVISION',
'None_', 'ParentRelative', 'CopyFromParent', 'PointerWindow', 'InputFocus',
'PointerRoot', 'AnyPropertyType', 'AnyKey', 'AnyButton', 'AllTemporary',
'CurrentTime', 'NoSymbol', 'NoEventMask', 'KeyPressMask', 'KeyReleaseMask',
'ButtonPressMask', 'ButtonReleaseMask', 'EnterWindowMask', 'LeaveWindowMask',
'PointerMotionMask', 'PointerMotionHintMask', 'Button1MotionMask',
'Button2MotionMask', 'Button3MotionMask', 'Button4MotionMask',
'Button5MotionMask', 'ButtonMotionMask', 'KeymapStateMask', 'ExposureMask',
'VisibilityChangeMask', 'StructureNotifyMask', 'ResizeRedirectMask',
'SubstructureNotifyMask', 'SubstructureRedirectMask', 'FocusChangeMask',
'PropertyChangeMask', 'ColormapChangeMask', 'OwnerGrabButtonMask', 'KeyPress',
'KeyRelease', 'ButtonPress', 'ButtonRelease', 'MotionNotify', 'EnterNotify',
'LeaveNotify', 'FocusIn', 'FocusOut', 'KeymapNotify', 'Expose',
'GraphicsExpose', 'NoExpose', 'VisibilityNotify', 'CreateNotify',
'DestroyNotify', 'UnmapNotify', 'MapNotify', 'MapRequest', 'ReparentNotify',
'ConfigureNotify', 'ConfigureRequest', 'GravityNotify', 'ResizeRequest',
'CirculateNotify', 'CirculateRequest', 'PropertyNotify', 'SelectionClear',
'SelectionRequest', 'SelectionNotify', 'ColormapNotify', 'ClientMessage',
'MappingNotify', 'LASTEvent', 'ShiftMask', 'LockMask', 'ControlMask',
'Mod1Mask', 'Mod2Mask', 'Mod3Mask', 'Mod4Mask', 'Mod5Mask', 'ShiftMapIndex',
'LockMapIndex', 'ControlMapIndex', 'Mod1MapIndex', 'Mod2MapIndex',
'Mod3MapIndex', 'Mod4MapIndex', 'Mod5MapIndex', 'Button1Mask', 'Button2Mask',
'Button3Mask', 'Button4Mask', 'Button5Mask', 'AnyModifier', 'Button1',
'Button2', 'Button3', 'Button4', 'Button5', 'NotifyNormal', 'NotifyGrab',
'NotifyUngrab', 'NotifyWhileGrabbed', 'NotifyHint', 'NotifyAncestor',
'NotifyVirtual', 'NotifyInferior', 'NotifyNonlinear',
'NotifyNonlinearVirtual', 'NotifyPointer', 'NotifyPointerRoot',
'NotifyDetailNone', 'VisibilityUnobscured', 'VisibilityPartiallyObscured',
'VisibilityFullyObscured', 'PlaceOnTop', 'PlaceOnBottom', 'FamilyInternet',
'FamilyDECnet', 'FamilyChaos', 'FamilyInternet6', 'FamilyServerInterpreted',
'PropertyNewValue', 'PropertyDelete', 'ColormapUninstalled',
'ColormapInstalled', 'GrabModeSync', 'GrabModeAsync', 'GrabSuccess',
'AlreadyGrabbed', 'GrabInvalidTime', 'GrabNotViewable', 'GrabFrozen',
'AsyncPointer', 'SyncPointer', 'ReplayPointer', 'AsyncKeyboard',
'SyncKeyboard', 'ReplayKeyboard', 'AsyncBoth', 'SyncBoth', 'RevertToParent',
'Success', 'BadRequest', 'BadValue', 'BadWindow', 'BadPixmap', 'BadAtom',
'BadCursor', 'BadFont', 'BadMatch', 'BadDrawable', 'BadAccess', 'BadAlloc',
'BadColor', 'BadGC', 'BadIDChoice', 'BadName', 'BadLength',
'BadImplementation', 'FirstExtensionError', 'LastExtensionError',
'InputOutput', 'InputOnly', 'CWBackPixmap', 'CWBackPixel', 'CWBorderPixmap',
'CWBorderPixel', 'CWBitGravity', 'CWWinGravity', 'CWBackingStore',
'CWBackingPlanes', 'CWBackingPixel', 'CWOverrideRedirect', 'CWSaveUnder',
'CWEventMask', 'CWDontPropagate', 'CWColormap', 'CWCursor', 'CWX', 'CWY',
'CWWidth', 'CWHeight', 'CWBorderWidth', 'CWSibling', 'CWStackMode',
'ForgetGravity', 'NorthWestGravity', 'NorthGravity', 'NorthEastGravity',
'WestGravity', 'CenterGravity', 'EastGravity', 'SouthWestGravity',
'SouthGravity', 'SouthEastGravity', 'StaticGravity', 'UnmapGravity',
'NotUseful', 'WhenMapped', 'Always', 'IsUnmapped', 'IsUnviewable',
'IsViewable', 'SetModeInsert', 'SetModeDelete', 'DestroyAll',
'RetainPermanent', 'RetainTemporary', 'Above', 'Below', 'TopIf', 'BottomIf',
'Opposite', 'RaiseLowest', 'LowerHighest', 'PropModeReplace',
'PropModePrepend', 'PropModeAppend', 'GXclear', 'GXand', 'GXandReverse',
'GXcopy', 'GXandInverted', 'GXnoop', 'GXxor', 'GXor', 'GXnor', 'GXequiv',
'GXinvert', 'GXorReverse', 'GXcopyInverted', 'GXorInverted', 'GXnand',
'GXset', 'LineSolid', 'LineOnOffDash', 'LineDoubleDash', 'CapNotLast',
'CapButt', 'CapRound', 'CapProjecting', 'JoinMiter', 'JoinRound', 'JoinBevel',
'FillSolid', 'FillTiled', 'FillStippled', 'FillOpaqueStippled', 'EvenOddRule',
'WindingRule', 'ClipByChildren', 'IncludeInferiors', 'Unsorted', 'YSorted',
'YXSorted', 'YXBanded', 'CoordModeOrigin', 'CoordModePrevious', 'Complex',
'Nonconvex', 'Convex', 'ArcChord', 'ArcPieSlice', 'GCFunction', 'GCPlaneMask',
'GCForeground', 'GCBackground', 'GCLineWidth', 'GCLineStyle', 'GCCapStyle',
'GCJoinStyle', 'GCFillStyle', 'GCFillRule', 'GCTile', 'GCStipple',
'GCTileStipXOrigin', 'GCTileStipYOrigin', 'GCFont', 'GCSubwindowMode',
'GCGraphicsExposures', 'GCClipXOrigin', 'GCClipYOrigin', 'GCClipMask',
'GCDashOffset', 'GCDashList', 'GCArcMode', 'GCLastBit', 'FontLeftToRight',
'FontRightToLeft', 'FontChange', 'XYBitmap', 'XYPixmap', 'ZPixmap',
'AllocNone', 'AllocAll', 'DoRed', 'DoGreen', 'DoBlue', 'CursorShape',
'TileShape', 'StippleShape', 'AutoRepeatModeOff', 'AutoRepeatModeOn',
'AutoRepeatModeDefault', 'LedModeOff', 'LedModeOn', 'KBKeyClickPercent',
'KBBellPercent', 'KBBellPitch', 'KBBellDuration', 'KBLed', 'KBLedMode',
'KBKey', 'KBAutoRepeatMode', 'MappingSuccess', 'MappingBusy', 'MappingFailed',
'MappingModifier', 'MappingKeyboard', 'MappingPointer', 'DontPreferBlanking',
'PreferBlanking', 'DefaultBlanking', 'DisableScreenSaver',
'DisableScreenInterval', 'DontAllowExposures', 'AllowExposures',
'DefaultExposures', 'ScreenSaverReset', 'ScreenSaverActive', 'HostInsert',
'HostDelete', 'EnableAccess', 'DisableAccess', 'StaticGray', 'GrayScale',
'StaticColor', 'PseudoColor', 'TrueColor', 'DirectColor', 'LSBFirst',
'MSBFirst', '_Xmblen', 'X_HAVE_UTF8_STRING', 'Bool', 'Status', 'True',
'False', 'QueuedAlready', 'QueuedAfterReading', 'QueuedAfterFlush',
'XIMPreeditArea', 'XIMPreeditCallbacks', 'XIMPreeditPosition',
'XIMPreeditNothing', 'XIMPreeditNone', 'XIMStatusArea', 'XIMStatusCallbacks',
'XIMStatusNothing', 'XIMStatusNone', 'XBufferOverflow', 'XLookupNone',
'XLookupChars', 'XLookupKeySym', 'XLookupBoth', 'XIMReverse', 'XIMUnderline',
'XIMHighlight', 'XIMPrimary', 'XIMSecondary', 'XIMTertiary',
'XIMVisibleToForward', 'XIMVisibleToBackword', 'XIMVisibleToCenter',
'XIMPreeditUnKnown', 'XIMPreeditEnable', 'XIMPreeditDisable',
'XIMInitialState', 'XIMPreserveState', 'XIMStringConversionLeftEdge',
'XIMStringConversionRightEdge', 'XIMStringConversionTopEdge',
'XIMStringConversionBottomEdge', 'XIMStringConversionConcealed',
'XIMStringConversionWrapped', 'XIMStringConversionBuffer',
'XIMStringConversionLine', 'XIMStringConversionWord',
'XIMStringConversionChar', 'XIMStringConversionSubstitution',
'XIMStringConversionRetrieval', 'XIMHotKeyStateON', 'XIMHotKeyStateOFF',
'XLoadQueryFont', 'XQueryFont', 'XGetMotionEvents', 'XDeleteModifiermapEntry',
'XGetModifierMapping', 'XInsertModifiermapEntry', 'XNewModifiermap',
'XCreateImage', 'XInitImage', 'XGetImage', 'XGetSubImage', 'XOpenDisplay',
'XrmInitialize', 'XFetchBytes', 'XFetchBuffer', 'XGetAtomName',
'XGetAtomNames', 'XGetDefault', 'XDisplayName', 'XKeysymToString',
'XSynchronize', 'XSetAfterFunction', 'XInternAtom', 'XInternAtoms',
'XCopyColormapAndFree', 'XCreateColormap', 'XCreatePixmapCursor',
'XCreateGlyphCursor', 'XCreateFontCursor', 'XLoadFont', 'XCreateGC',
'XGContextFromGC', 'XFlushGC', 'XCreatePixmap', 'XCreateBitmapFromData',
'XCreatePixmapFromBitmapData', 'XCreateSimpleWindow', 'XGetSelectionOwner',
'XCreateWindow', 'XListInstalledColormaps', 'XListFonts',
'XListFontsWithInfo', 'XGetFontPath', 'XListExtensions', 'XListProperties',
'XListHosts', 'XKeycodeToKeysym', 'XLookupKeysym', 'XGetKeyboardMapping',
'XStringToKeysym', 'XMaxRequestSize', 'XExtendedMaxRequestSize',
'XResourceManagerString', 'XScreenResourceString', 'XDisplayMotionBufferSize',
'XVisualIDFromVisual', 'XInitThreads', 'XLockDisplay', 'XUnlockDisplay',
'XInitExtension', 'XAddExtension', 'XFindOnExtensionList',
'XEHeadOfExtensionList', 'XRootWindow', 'XDefaultRootWindow',
'XRootWindowOfScreen', 'XDefaultVisual', 'XDefaultVisualOfScreen',
'XDefaultGC', 'XDefaultGCOfScreen', 'XBlackPixel', 'XWhitePixel',
'XAllPlanes', 'XBlackPixelOfScreen', 'XWhitePixelOfScreen', 'XNextRequest',
'XLastKnownRequestProcessed', 'XServerVendor', 'XDisplayString',
'XDefaultColormap', 'XDefaultColormapOfScreen', 'XDisplayOfScreen',
'XScreenOfDisplay', 'XDefaultScreenOfDisplay', 'XEventMaskOfScreen',
'XScreenNumberOfScreen', 'XSetErrorHandler', 'XSetIOErrorHandler',
'XListPixmapFormats', 'XListDepths', 'XReconfigureWMWindow',
'XGetWMProtocols', 'XSetWMProtocols', 'XIconifyWindow', 'XWithdrawWindow',
'XGetCommand', 'XGetWMColormapWindows', 'XSetWMColormapWindows',
'XFreeStringList', 'XSetTransientForHint', 'XActivateScreenSaver', 'XAddHost',
'XAddHosts', 'XAddToExtensionList', 'XAddToSaveSet', 'XAllocColor',
'XAllocColorCells', 'XAllocColorPlanes', 'XAllocNamedColor', 'XAllowEvents',
'XAutoRepeatOff', 'XAutoRepeatOn', 'XBell', 'XBitmapBitOrder', 'XBitmapPad',
'XBitmapUnit', 'XCellsOfScreen', 'XChangeActivePointerGrab', 'XChangeGC',
'XChangeKeyboardControl', 'XChangeKeyboardMapping', 'XChangePointerControl',
'XChangeProperty', 'XChangeSaveSet', 'XChangeWindowAttributes',
'XCheckIfEvent', 'XCheckMaskEvent', 'XCheckTypedEvent',
'XCheckTypedWindowEvent', 'XCheckWindowEvent', 'XCirculateSubwindows',
'XCirculateSubwindowsDown', 'XCirculateSubwindowsUp', 'XClearArea',
'XClearWindow', 'XCloseDisplay', 'XConfigureWindow', 'XConnectionNumber',
'XConvertSelection', 'XCopyArea', 'XCopyGC', 'XCopyPlane', 'XDefaultDepth',
'XDefaultDepthOfScreen', 'XDefaultScreen', 'XDefineCursor', 'XDeleteProperty',
'XDestroyWindow', 'XDestroySubwindows', 'XDoesBackingStore',
'XDoesSaveUnders', 'XDisableAccessControl', 'XDisplayCells', 'XDisplayHeight',
'XDisplayHeightMM', 'XDisplayKeycodes', 'XDisplayPlanes', 'XDisplayWidth',
'XDisplayWidthMM', 'XDrawArc', 'XDrawArcs', 'XDrawImageString',
'XDrawImageString16', 'XDrawLine', 'XDrawLines', 'XDrawPoint', 'XDrawPoints',
'XDrawRectangle', 'XDrawRectangles', 'XDrawSegments', 'XDrawString',
'XDrawString16', 'XDrawText', 'XDrawText16', 'XEnableAccessControl',
'XEventsQueued', 'XFetchName', 'XFillArc', 'XFillArcs', 'XFillPolygon',
'XFillRectangle', 'XFillRectangles', 'XFlush', 'XForceScreenSaver', 'XFree',
'XFreeColormap', 'XFreeColors', 'XFreeCursor', 'XFreeExtensionList',
'XFreeFont', 'XFreeFontInfo', 'XFreeFontNames', 'XFreeFontPath', 'XFreeGC',
'XFreeModifiermap', 'XFreePixmap', 'XGeometry', 'XGetErrorDatabaseText',
'XGetErrorText', 'XGetFontProperty', 'XGetGCValues', 'XGetGeometry',
'XGetIconName', 'XGetInputFocus', 'XGetKeyboardControl', 'XGetPointerControl',
'XGetPointerMapping', 'XGetScreenSaver', 'XGetTransientForHint',
'XGetWindowProperty', 'XGetWindowAttributes', 'XGrabButton', 'XGrabKey',
'XGrabKeyboard', 'XGrabPointer', 'XGrabServer', 'XHeightMMOfScreen',
'XHeightOfScreen', 'XIfEvent', 'XImageByteOrder', 'XInstallColormap',
'XKeysymToKeycode', 'XKillClient', 'XLookupColor', 'XLowerWindow',
'XMapRaised', 'XMapSubwindows', 'XMapWindow', 'XMaskEvent',
'XMaxCmapsOfScreen', 'XMinCmapsOfScreen', 'XMoveResizeWindow', 'XMoveWindow',
'XNextEvent', 'XNoOp', 'XParseColor', 'XParseGeometry', 'XPeekEvent',
'XPeekIfEvent', 'XPending', 'XPlanesOfScreen', 'XProtocolRevision',
'XProtocolVersion', 'XPutBackEvent', 'XPutImage', 'XQLength',
'XQueryBestCursor', 'XQueryBestSize', 'XQueryBestStipple', 'XQueryBestTile',
'XQueryColor', 'XQueryColors', 'XQueryExtension', 'XQueryKeymap',
'XQueryPointer', 'XQueryTextExtents', 'XQueryTextExtents16', 'XQueryTree',
'XRaiseWindow', 'XReadBitmapFile', 'XReadBitmapFileData', 'XRebindKeysym',
'XRecolorCursor', 'XRefreshKeyboardMapping', 'XRemoveFromSaveSet',
'XRemoveHost', 'XRemoveHosts', 'XReparentWindow', 'XResetScreenSaver',
'XResizeWindow', 'XRestackWindows', 'XRotateBuffers',
'XRotateWindowProperties', 'XScreenCount', 'XSelectInput', 'XSendEvent',
'XSetAccessControl', 'XSetArcMode', 'XSetBackground', 'XSetClipMask',
'XSetClipOrigin', 'XSetClipRectangles', 'XSetCloseDownMode', 'XSetCommand',
'XSetDashes', 'XSetFillRule', 'XSetFillStyle', 'XSetFont', 'XSetFontPath',
'XSetForeground', 'XSetFunction', 'XSetGraphicsExposures', 'XSetIconName',
'XSetInputFocus', 'XSetLineAttributes', 'XSetModifierMapping',
'XSetPlaneMask', 'XSetPointerMapping', 'XSetScreenSaver',
'XSetSelectionOwner', 'XSetState', 'XSetStipple', 'XSetSubwindowMode',
'XSetTSOrigin', 'XSetTile', 'XSetWindowBackground',
'XSetWindowBackgroundPixmap', 'XSetWindowBorder', 'XSetWindowBorderPixmap',
'XSetWindowBorderWidth', 'XSetWindowColormap', 'XStoreBuffer', 'XStoreBytes',
'XStoreColor', 'XStoreColors', 'XStoreName', 'XStoreNamedColor', 'XSync',
'XTextExtents', 'XTextExtents16', 'XTextWidth', 'XTextWidth16',
'XTranslateCoordinates', 'XUndefineCursor', 'XUngrabButton', 'XUngrabKey',
'XUngrabKeyboard', 'XUngrabPointer', 'XUngrabServer', 'XUninstallColormap',
'XUnloadFont', 'XUnmapSubwindows', 'XUnmapWindow', 'XVendorRelease',
'XWarpPointer', 'XWidthMMOfScreen', 'XWidthOfScreen', 'XWindowEvent',
'XWriteBitmapFile', 'XSupportsLocale', 'XSetLocaleModifiers', 'XOpenOM',
'XCloseOM', 'XSetOMValues', 'XGetOMValues', 'XDisplayOfOM', 'XLocaleOfOM',
'XCreateOC', 'XDestroyOC', 'XOMOfOC', 'XSetOCValues', 'XGetOCValues',
'XCreateFontSet', 'XFreeFontSet', 'XFontsOfFontSet',
'XBaseFontNameListOfFontSet', 'XLocaleOfFontSet', 'XContextDependentDrawing',
'XDirectionalDependentDrawing', 'XContextualDrawing', 'XExtentsOfFontSet',
'XmbTextEscapement', 'XwcTextEscapement', 'Xutf8TextEscapement',
'XmbTextExtents', 'XwcTextExtents', 'Xutf8TextExtents',
'XmbTextPerCharExtents', 'XwcTextPerCharExtents', 'Xutf8TextPerCharExtents',
'XmbDrawText', 'XwcDrawText', 'Xutf8DrawText', 'XmbDrawString',
'XwcDrawString', 'Xutf8DrawString', 'XmbDrawImageString',
'XwcDrawImageString', 'Xutf8DrawImageString', 'XOpenIM', 'XCloseIM',
'XGetIMValues', 'XSetIMValues', 'XDisplayOfIM', 'XLocaleOfIM', 'XCreateIC',
'XDestroyIC', 'XSetICFocus', 'XUnsetICFocus', 'XwcResetIC', 'XmbResetIC',
'Xutf8ResetIC', 'XSetICValues', 'XGetICValues', 'XIMOfIC', 'XFilterEvent',
'XmbLookupString', 'XwcLookupString', 'Xutf8LookupString',
'XVaCreateNestedList', 'XRegisterIMInstantiateCallback',
'XUnregisterIMInstantiateCallback', 'XInternalConnectionNumbers',
'XProcessInternalConnection', 'XAddConnectionWatch', 'XRemoveConnectionWatch',
'XSetAuthorization', '_Xmbtowc', '_Xwctomb', 'NoValue', 'XValue', 'YValue',
'WidthValue', 'HeightValue', 'AllValues', 'XNegative', 'YNegative',
'XSizeHints', 'USPosition', 'USSize', 'PPosition', 'PSize', 'PMinSize',
'PMaxSize', 'PResizeInc', 'PAspect', 'PBaseSize', 'PWinGravity', 'PAllHints',
'XWMHints', 'InputHint', 'StateHint', 'IconPixmapHint', 'IconWindowHint',
'IconPositionHint', 'IconMaskHint', 'WindowGroupHint', 'AllHints',
'XUrgencyHint', 'WithdrawnState', 'NormalState', 'IconicState',
'DontCareState', 'ZoomState', 'InactiveState', 'XTextProperty', 'XNoMemory',
'XLocaleNotSupported', 'XConverterNotFound', 'XICCEncodingStyle',
'XStringStyle', 'XCompoundTextStyle', 'XTextStyle', 'XStdICCTextStyle',
'XUTF8StringStyle', 'XIconSize', 'XClassHint', 'XComposeStatus', 'Region',
'RectangleOut', 'RectangleIn', 'RectanglePart', 'XVisualInfo', 'VisualNoMask',
'VisualIDMask', 'VisualScreenMask', 'VisualDepthMask', 'VisualClassMask',
'VisualRedMaskMask', 'VisualGreenMaskMask', 'VisualBlueMaskMask',
'VisualColormapSizeMask', 'VisualBitsPerRGBMask', 'VisualAllMask',
'XStandardColormap', 'BitmapSuccess', 'BitmapOpenFailed', 'BitmapFileInvalid',
'BitmapNoMemory', 'XCSUCCESS', 'XCNOMEM', 'XCNOENT', 'XContext',
'XAllocClassHint', 'XAllocIconSize', 'XAllocSizeHints',
'XAllocStandardColormap', 'XAllocWMHints', 'XClipBox', 'XCreateRegion',
'XDefaultString', 'XDeleteContext', 'XDestroyRegion', 'XEmptyRegion',
'XEqualRegion', 'XFindContext', 'XGetClassHint', 'XGetIconSizes',
'XGetNormalHints', 'XGetRGBColormaps', 'XGetSizeHints',
'XGetStandardColormap', 'XGetTextProperty', 'XGetVisualInfo',
'XGetWMClientMachine', 'XGetWMHints', 'XGetWMIconName', 'XGetWMName',
'XGetWMNormalHints', 'XGetWMSizeHints', 'XGetZoomHints', 'XIntersectRegion',
'XConvertCase', 'XLookupString', 'XMatchVisualInfo', 'XOffsetRegion',
'XPointInRegion', 'XPolygonRegion', 'XRectInRegion', 'XSaveContext',
'XSetClassHint', 'XSetIconSizes', 'XSetNormalHints', 'XSetRGBColormaps',
'XSetSizeHints', 'XSetStandardProperties', 'XSetTextProperty',
'XSetWMClientMachine', 'XSetWMHints', 'XSetWMIconName', 'XSetWMName',
'XSetWMNormalHints', 'XSetWMProperties', 'XmbSetWMProperties',
'Xutf8SetWMProperties', 'XSetWMSizeHints', 'XSetRegion',
'XSetStandardColormap', 'XSetZoomHints', 'XShrinkRegion',
'XStringListToTextProperty', 'XSubtractRegion', 'XmbTextListToTextProperty',
'XwcTextListToTextProperty', 'Xutf8TextListToTextProperty',
'XwcFreeStringList', 'XTextPropertyToStringList', 'XmbTextPropertyToTextList',
'XwcTextPropertyToTextList', 'Xutf8TextPropertyToTextList',
'XUnionRectWithRegion', 'XUnionRegion', 'XWMGeometry', 'XXorRegion']
