"""
Functions to convert between different LaTeX color spaces.

Main purpose is to be able to find one valid color space among those provided.

Will always return the original string if not possible to convert and not say
anything about it.
"""

# -*- coding: utf-8 -*-

from enum import Enum
from collections import namedtuple

ColorType = namedtuple("ColorType", ["value", "int_type", "to_str"])


class ColorTypes(Enum):
    @property
    def int_type(self):
        return self.value.int_type

    def to_str(self, s):
        return self.value.to_str(s)

    STRING = ColorType(1, None, lambda x: "unknown")
    rgb = ColorType(
        2, False, lambda x: "rgb {}".format(" ".join(["{:.5g}".format(f) for f in x]))
    )
    RGB = ColorType(
        3, True, lambda x: "RGB {}".format(" ".join(["{}".format(f) for f in x]))
    )
    HTML = ColorType(
        4, True, lambda x: "HTML {}".format("".join(["{:02X}".format(f) for f in x]))
    )
    cmyk = ColorType(
        5,
        False,
        lambda x: "rgb {}".format(
            " ".join(["{:.5g}".format(f) for f in rgb_to_cmyk(x)])
        ),
    )
    rgbhex = ColorType(
        6, True, lambda x: "#{}".format("".join(["{:02X}".format(f) for f in x]))
    )
    gray = ColorType(7, False, lambda x: "gray {:.5g}".format(rgb_to_cmyk(x)[3]))


def get_color_type(s):
    s = s.strip()
    if s.startswith("rgb"):
        p = [x for x in s.split(" ") if x]
        if len(p) == 4:
            return ColorTypes.rgb
    if s.startswith("RGB"):
        p = [x for x in s.split(" ") if x]
        if len(p) == 4:
            return ColorTypes.RGB
    if s.startswith("cmyk"):
        p = [x for x in s.split(" ") if x]
        if len(p) == 5:
            return ColorTypes.cmyk
    if s.startswith("HTML"):
        p = [x for x in s.split(" ") if x]
        if len(p) == 2 and len(p[1]) == 6:
            return ColorTypes.HTML
    if s.startswith("#"):
        if len(s) == 7:
            return ColorTypes.rgbhex
    if s.startswith("gray"):
        p = [x for x in s.split(" ") if x]
        if len(p) == 2:
            return ColorTypes.gray
    return ColorTypes.STRING


def convert_color_to(s, valid_types):
    t = get_color_type(s)
    if t == ColorTypes.STRING or t in valid_types:
        return s

    val, int_type = _get_color_value(s, t)
    int_types = [x for x in valid_types if x.int_type == True]
    float_types = [x for x in valid_types if x.int_type == False]
    if int_type == True:
        if int_types:
            return int_types[0].to_str(val)
        elif float_types:
            return float_types[0].to_str([x / 255 for x in val])
    elif int_type == False:
        if float_types:
            return float_types[0].to_str(val)
        elif int_types:
            return int_types[0].to_str([round(x * 255) for x in val])
    return s


def _get_color_value(s, t=None):
    if t is None:
        t = get_color_type(s)
        if t == ColorTypes.STRING:
            return (None, None)
    if t == ColorTypes.gray:
        p = [x for x in s.split(" ") if x]
        if len(p) == 2:
            try:
                f = float(p[1])
                return ([f, f, f], t.int_type)
            except ValueError:
                pass
    elif t == ColorTypes.rgb:
        p = [x for x in s.split(" ") if x]
        if len(p) == 4:
            try:
                r = float(p[1])
                g = float(p[2])
                b = float(p[3])
                return ([r, g, b], t.int_type)
            except ValueError:
                pass
    elif t == ColorTypes.RGB:
        p = [x for x in s.split(" ") if x]
        if len(p) == 4:
            try:
                R = int(p[1])
                G = int(p[2])
                B = int(p[3])
                return ([R, G, B], t.int_type)
            except ValueError:
                pass
    elif t == ColorTypes.cmyk:
        p = [x for x in s.split(" ") if x]
        if len(p) == 5:
            try:
                c = float(p[1])
                m = float(p[2])
                y = float(p[3])
                k = float(p[4])
                return (cmyk_to_rgb([c, m, y, k]), t.int_type)
            except ValueError:
                pass
    elif t == ColorTypes.HTML:
        p = [x for x in s.split(" ") if x]
        if len(p) == 2 and len(p[1]) == 6:
            try:
                R = int(p[1][0:2], 16)
                G = int(p[1][2:4], 16)
                B = int(p[1][4:], 16)
                return ([R, G, B], t.int_type)
            except ValueError:
                pass
    elif t == ColorTypes.rgbhex:
        p = s.strip()
        if len(p) == 7:
            try:
                R = int(p[1:3], 16)
                G = int(p[3:5], 16)
                B = int(p[5:], 16)
                return ([R, G, B], t.int_type)
            except ValueError:
                pass
    return (None, None)


def rgb_to_cmyk(rgb):
    """
    Convert from rgb to cmyk.

    Assume that rgb values are floating-point between 0 and 1.
    """
    r = rgb[0]
    g = rgb[1]
    b = rgb[2]
    k = 1 - max(r, g, b)
    if k != 1:
        c = (1 - r - k) / (1 - k)
        m = (1 - g - k) / (1 - k)
        y = (1 - b - k) / (1 - k)
    else:
        c = 0
        m = 0
        y = 0
    return [c, m, y, k]


def cmyk_to_rgb(cmyk):
    """
    Convert from cmyk to rgb.

    Assume that cmyk values are floating-point between 0 and 1.
    """
    c = cmyk[0]
    m = cmyk[1]
    y = cmyk[2]
    k = cmyk[3]
    r = (1 - c) * (1 - k)
    g = (1 - m) * (1 - k)
    b = (1 - y) * (1 - k)
    return [r, g, b]


_latex_colors = {
    "red": ("red", [1, 0, 0]),
    "green": ("green", [0, 1, 0]),
    "blue": ("blue", [0, 0, 1]),
    "brown": ("brown", [0.75, 0.5, 0.25]),
    "lime": ("lime", [0.75, 1, 0]),
    "orange": ("orange", [1, 0.5, 0]),
    "pink": ("pink", [1, 0.75, 0.75]),
    "purple": ("purple", [0.75, 0, 0.25]),
    "teal": ("teal", [0, 0.5, 0.5]),
    "violet": ("violet", [0.5, 0, 0.5]),
    "cyan": ("cyan", [0, 1, 1]),
    "magenta": ("magenta", [1, 0, 1]),
    "yellow": ("yellow", [1, 1, 0]),
    "olive": ("olive", [0.5, 0.5, 0]),
    "black": ("black", [0, 0, 0]),
    "darkgray": ("darkgray", [0.25, 0.25, 0.25]),
    "gray": ("gray", [0.5, 0.5, 0.5]),
    "lightgray": ("lightgray", [0.75, 0.75, 0.75]),
    "white": ("white", [1, 1, 1]),
}

_svg_colors = {
    "aliceblue": ("AliceBlue", [0.94, 0.972, 1]),
    "antiquewhite": ("AntiqueWhite", [0.98, 0.92, 0.844]),
    "aqua": ("Aqua", [0, 1, 1]),
    "aquamarine": ("Aquamarine", [0.498, 1, 0.83]),
    "azure": ("Azure", [0.94, 1, 1]),
    "beige": ("Beige", [0.96, 0.96, 0.864]),
    "bisque": ("Bisque", [1, 0.894, 0.77]),
    "black": ("Black", [0, 0, 0]),
    "blanchedalmond": ("BlanchedAlmond", [1, 0.92, 0.804]),
    "blue": ("Blue", [0, 0, 1]),
    "blueviolet": ("BlueViolet", [0.54, 0.17, 0.888]),
    "brown": ("Brown", [0.648, 0.165, 0.165]),
    "burlywood": ("BurlyWood", [0.87, 0.72, 0.53]),
    "cadetblue": ("CadetBlue", [0.372, 0.62, 0.628]),
    "chartreuse": ("Chartreuse", [0.498, 1, 0]),
    "chocolate": ("Chocolate", [0.824, 0.41, 0.116]),
    "coral": ("Coral", [1, 0.498, 0.312]),
    "cornflowerblue": ("CornflowerBlue", [0.392, 0.585, 0.93]),
    "cornsilk": ("Cornsilk", [1, 0.972, 0.864]),
    "crimson": ("Crimson", [0.864, 0.08, 0.235]),
    "cyan": ("Cyan", [0, 1, 1]),
    "darkblue": ("DarkBlue", [0, 0, 0.545]),
    "darkcyan": ("DarkCyan", [0, 0.545, 0.545]),
    "darkgoldenrod": ("DarkGoldenrod", [0.72, 0.525, 0.044]),
    "darkgray": ("DarkGray", [0.664, 0.664, 0.664]),
    "darkgreen": ("DarkGreen", [0, 0.392, 0]),
    "darkgrey": ("DarkGrey", [0.664, 0.664, 0.664]),
    "darkkhaki": ("DarkKhaki", [0.74, 0.716, 0.42]),
    "darkmagenta": ("DarkMagenta", [0.545, 0, 0.545]),
    "darkolivegreen": ("DarkOliveGreen", [0.332, 0.42, 0.185]),
    "darkorange": ("DarkOrange", [1, 0.55, 0]),
    "darkorchid": ("DarkOrchid", [0.6, 0.196, 0.8]),
    "darkred": ("DarkRed", [0.545, 0, 0]),
    "darksalmon": ("DarkSalmon", [0.912, 0.59, 0.48]),
    "darkseagreen": ("DarkSeaGreen", [0.56, 0.736, 0.56]),
    "darkslateblue": ("DarkSlateBlue", [0.284, 0.24, 0.545]),
    "darkslategray": ("DarkSlateGray", [0.185, 0.31, 0.31]),
    "darkslategrey": ("DarkSlateGrey", [0.185, 0.31, 0.31]),
    "darkturquoise": ("DarkTurquoise", [0, 0.808, 0.82]),
    "darkviolet": ("DarkViolet", [0.58, 0, 0.828]),
    "deeppink": ("DeepPink", [1, 0.08, 0.576]),
    "deepskyblue": ("DeepSkyBlue", [0, 0.75, 1]),
    "dimgray": ("DimGray", [0.41, 0.41, 0.41]),
    "dimgrey": ("DimGrey", [0.41, 0.41, 0.41]),
    "dodgerblue": ("DodgerBlue", [0.116, 0.565, 1]),
    "firebrick": ("FireBrick", [0.698, 0.132, 0.132]),
    "floralwhite": ("FloralWhite", [1, 0.98, 0.94]),
    "forestgreen": ("ForestGreen", [0.132, 0.545, 0.132]),
    "fuchsia": ("Fuchsia", [1, 0, 1]),
    "gainsboro": ("Gainsboro", [0.864, 0.864, 0.864]),
    "ghostwhite": ("GhostWhite", [0.972, 0.972, 1]),
    "gold": ("Gold", [1, 0.844, 0]),
    "goldenrod": ("Goldenrod", [0.855, 0.648, 0.125]),
    "gray": ("Gray", [0.5, 0.5, 0.5]),
    "green": ("Green", [0, 0.5, 0]),
    "greenyellow": ("GreenYellow", [0.68, 1, 0.185]),
    "grey": ("Grey", [0.5, 0.5, 0.5]),
    "honeydew": ("Honeydew", [0.94, 1, 0.94]),
    "hotpink": ("HotPink", [1, 0.41, 0.705]),
    "indianred": ("IndianRed", [0.804, 0.36, 0.36]),
    "indigo": ("Indigo", [0.294, 0, 0.51]),
    "ivory": ("Ivory", [1, 1, 0.94]),
    "khaki": ("Khaki", [0.94, 0.9, 0.55]),
    "lavender": ("Lavender", [0.9, 0.9, 0.98]),
    "lavenderblush": ("LavenderBlush", [1, 0.94, 0.96]),
    "lawngreen": ("LawnGreen", [0.488, 0.99, 0]),
    "lemonchiffon": ("LemonChiffon", [1, 0.98, 0.804]),
    "lightblue": ("LightBlue", [0.68, 0.848, 0.9]),
    "lightcoral": ("LightCoral", [0.94, 0.5, 0.5]),
    "lightcyan": ("LightCyan", [0.88, 1, 1]),
    "lightgoldenrod": ("LightGoldenrod", [0.933, 0.867, 0.51]),
    "lightgoldenrodyellow": ("LightGoldenrodYellow", [0.98, 0.98, 0.824]),
    "lightgray": ("LightGray", [0.828, 0.828, 0.828]),
    "lightgreen": ("LightGreen", [0.565, 0.932, 0.565]),
    "lightgrey": ("LightGrey", [0.828, 0.828, 0.828]),
    "lightpink": ("LightPink", [1, 0.712, 0.756]),
    "lightsalmon": ("LightSalmon", [1, 0.628, 0.48]),
    "lightseagreen": ("LightSeaGreen", [0.125, 0.698, 0.668]),
    "lightskyblue": ("LightSkyBlue", [0.53, 0.808, 0.98]),
    "lightslateblue": ("LightSlateBlue", [0.518, 0.44, 1]),
    "lightslategray": ("LightSlateGray", [0.468, 0.532, 0.6]),
    "lightslategrey": ("LightSlateGrey", [0.468, 0.532, 0.6]),
    "lightsteelblue": ("LightSteelBlue", [0.69, 0.77, 0.87]),
    "lightyellow": ("LightYellow", [1, 1, 0.88]),
    "lime": ("Lime", [0, 1, 0]),
    "limegreen": ("LimeGreen", [0.196, 0.804, 0.196]),
    "linen": ("Linen", [0.98, 0.94, 0.9]),
    "magenta": ("Magenta", [1, 0, 1]),
    "maroon": ("Maroon", [0.5, 0, 0]),
    "mediumaquamarine": ("MediumAquamarine", [0.4, 0.804, 0.668]),
    "mediumblue": ("MediumBlue", [0, 0, 0.804]),
    "mediumorchid": ("MediumOrchid", [0.73, 0.332, 0.828]),
    "mediumpurple": ("MediumPurple", [0.576, 0.44, 0.86]),
    "mediumseagreen": ("MediumSeaGreen", [0.235, 0.7, 0.444]),
    "mediumslateblue": ("MediumSlateBlue", [0.484, 0.408, 0.932]),
    "mediumspringgreen": ("MediumSpringGreen", [0, 0.98, 0.604]),
    "mediumturquoise": ("MediumTurquoise", [0.284, 0.82, 0.8]),
    "mediumvioletred": ("MediumVioletRed", [0.78, 0.084, 0.52]),
    "midnightblue": ("MidnightBlue", [0.098, 0.098, 0.44]),
    "mintcream": ("MintCream", [0.96, 1, 0.98]),
    "mistyrose": ("MistyRose", [1, 0.894, 0.884]),
    "moccasin": ("Moccasin", [1, 0.894, 0.71]),
    "navajowhite": ("NavajoWhite", [1, 0.87, 0.68]),
    "navy": ("Navy", [0, 0, 0.5]),
    "navyblue": ("NavyBlue", [0, 0, 0.5]),
    "oldlace": ("OldLace", [0.992, 0.96, 0.9]),
    "olive": ("Olive", [0.5, 0.5, 0]),
    "olivedrab": ("OliveDrab", [0.42, 0.556, 0.136]),
    "orange": ("Orange", [1, 0.648, 0]),
    "orangered": ("OrangeRed", [1, 0.27, 0]),
    "orchid": ("Orchid", [0.855, 0.44, 0.84]),
    "palegoldenrod": ("PaleGoldenrod", [0.932, 0.91, 0.668]),
    "palegreen": ("PaleGreen", [0.596, 0.985, 0.596]),
    "paleturquoise": ("PaleTurquoise", [0.688, 0.932, 0.932]),
    "palevioletred": ("PaleVioletRed", [0.86, 0.44, 0.576]),
    "papayawhip": ("PapayaWhip", [1, 0.936, 0.835]),
    "peachpuff": ("PeachPuff", [1, 0.855, 0.725]),
    "peru": ("Peru", [0.804, 0.52, 0.248]),
    "pink": ("Pink", [1, 0.752, 0.796]),
    "plum": ("Plum", [0.868, 0.628, 0.868]),
    "powderblue": ("PowderBlue", [0.69, 0.88, 0.9]),
    "purple": ("Purple", [0.5, 0, 0.5]),
    "red": ("Red", [1, 0, 0]),
    "rosybrown": ("RosyBrown", [0.736, 0.56, 0.56]),
    "royalblue": ("RoyalBlue", [0.255, 0.41, 0.884]),
    "saddlebrown": ("SaddleBrown", [0.545, 0.27, 0.075]),
    "salmon": ("Salmon", [0.98, 0.5, 0.448]),
    "sandybrown": ("SandyBrown", [0.956, 0.644, 0.376]),
    "seagreen": ("SeaGreen", [0.18, 0.545, 0.34]),
    "seashell": ("Seashell", [1, 0.96, 0.932]),
    "sienna": ("Sienna", [0.628, 0.32, 0.176]),
    "silver": ("Silver", [0.752, 0.752, 0.752]),
    "skyblue": ("SkyBlue", [0.53, 0.808, 0.92]),
    "slateblue": ("SlateBlue", [0.415, 0.352, 0.804]),
    "slategray": ("SlateGray", [0.44, 0.5, 0.565]),
    "slategrey": ("SlateGrey", [0.44, 0.5, 0.565]),
    "snow": ("Snow", [1, 0.98, 0.98]),
    "springgreen": ("SpringGreen", [0, 1, 0.498]),
    "steelblue": ("SteelBlue", [0.275, 0.51, 0.705]),
    "tan": ("Tan", [0.824, 0.705, 0.55]),
    "teal": ("Teal", [0, 0.5, 0.5]),
    "thistle": ("Thistle", [0.848, 0.75, 0.848]),
    "tomato": ("Tomato", [1, 0.39, 0.28]),
    "turquoise": ("Turquoise", [0.25, 0.88, 0.815]),
    "violet": ("Violet", [0.932, 0.51, 0.932]),
    "violetred": ("VioletRed", [0.816, 0.125, 0.565]),
    "wheat": ("Wheat", [0.96, 0.87, 0.7]),
    "white": ("White", [1, 1, 1]),
    "whitesmoke": ("WhiteSmoke", [0.96, 0.96, 0.96]),
    "yellow": ("Yellow", [1, 1, 0]),
    "yellowgreen": ("YellowGreen", [0.604, 0.804, 0.196]),
}

_x11_colors = {
    "antiquewhite1": ("AntiqueWhite1", (1, 0.936, 0.86)),
    "antiquewhite2": ("AntiqueWhite2", (0.932, 0.875, 0.8)),
    "antiquewhite3": ("AntiqueWhite3", (0.804, 0.752, 0.69)),
    "antiquewhite4": ("AntiqueWhite4", (0.545, 0.512, 0.47)),
    "aquamarine1": ("Aquamarine1", (0.498, 1, 0.83)),
    "aquamarine2": ("Aquamarine2", (0.464, 0.932, 0.776)),
    "aquamarine3": ("Aquamarine3", (0.4, 0.804, 0.668)),
    "aquamarine4": ("Aquamarine4", (0.27, 0.545, 0.455)),
    "azure1": ("Azure1", (0.94, 1, 1)),
    "azure2": ("Azure2", (0.88, 0.932, 0.932)),
    "azure3": ("Azure3", (0.756, 0.804, 0.804)),
    "azure4": ("Azure4", (0.512, 0.545, 0.545)),
    "bisque1": ("Bisque1", (1, 0.894, 0.77)),
    "bisque2": ("Bisque2", (0.932, 0.835, 0.716)),
    "bisque3": ("Bisque3", (0.804, 0.716, 0.62)),
    "bisque4": ("Bisque4", (0.545, 0.49, 0.42)),
    "blue1": ("Blue1", (0, 0, 1)),
    "blue2": ("Blue2", (0, 0, 0.932)),
    "blue3": ("Blue3", (0, 0, 0.804)),
    "blue4": ("Blue4", (0, 0, 0.545)),
    "brown1": ("Brown1", (1, 0.25, 0.25)),
    "brown2": ("Brown2", (0.932, 0.23, 0.23)),
    "brown3": ("Brown3", (0.804, 0.2, 0.2)),
    "brown4": ("Brown4", (0.545, 0.136, 0.136)),
    "burlywood1": ("Burlywood1", (1, 0.828, 0.608)),
    "burlywood2": ("Burlywood2", (0.932, 0.772, 0.57)),
    "burlywood3": ("Burlywood3", (0.804, 0.668, 0.49)),
    "burlywood4": ("Burlywood4", (0.545, 0.45, 0.332)),
    "cadetblue1": ("CadetBlue1", (0.596, 0.96, 1)),
    "cadetblue2": ("CadetBlue2", (0.556, 0.898, 0.932)),
    "cadetblue3": ("CadetBlue3", (0.48, 0.772, 0.804)),
    "cadetblue4": ("CadetBlue4", (0.325, 0.525, 0.545)),
    "chartreuse1": ("Chartreuse1", (0.498, 1, 0)),
    "chartreuse2": ("Chartreuse2", (0.464, 0.932, 0)),
    "chartreuse3": ("Chartreuse3", (0.4, 0.804, 0)),
    "chartreuse4": ("Chartreuse4", (0.27, 0.545, 0)),
    "chocolate1": ("Chocolate1", (1, 0.498, 0.14)),
    "chocolate2": ("Chocolate2", (0.932, 0.464, 0.13)),
    "chocolate3": ("Chocolate3", (0.804, 0.4, 0.112)),
    "chocolate4": ("Chocolate4", (0.545, 0.27, 0.075)),
    "coral1": ("Coral1", (1, 0.448, 0.336)),
    "coral2": ("Coral2", (0.932, 0.415, 0.312)),
    "coral3": ("Coral3", (0.804, 0.356, 0.27)),
    "coral4": ("Coral4", (0.545, 0.244, 0.185)),
    "cornsilk1": ("Cornsilk1", (1, 0.972, 0.864)),
    "cornsilk2": ("Cornsilk2", (0.932, 0.91, 0.804)),
    "cornsilk3": ("Cornsilk3", (0.804, 0.785, 0.694)),
    "cornsilk4": ("Cornsilk4", (0.545, 0.532, 0.47)),
    "cyan1": ("Cyan1", (0, 1, 1)),
    "cyan2": ("Cyan2", (0, 0.932, 0.932)),
    "cyan3": ("Cyan3", (0, 0.804, 0.804)),
    "cyan4": ("Cyan4", (0, 0.545, 0.545)),
    "darkgoldenrod1": ("DarkGoldenrod1", (1, 0.725, 0.06)),
    "darkgoldenrod2": ("DarkGoldenrod2", (0.932, 0.68, 0.055)),
    "darkgoldenrod3": ("DarkGoldenrod3", (0.804, 0.585, 0.048)),
    "darkgoldenrod4": ("DarkGoldenrod4", (0.545, 0.396, 0.03)),
    "darkolivegreen1": ("DarkOliveGreen1", (0.792, 1, 0.44)),
    "darkolivegreen2": ("DarkOliveGreen2", (0.736, 0.932, 0.408)),
    "darkolivegreen3": ("DarkOliveGreen3", (0.635, 0.804, 0.352)),
    "darkolivegreen4": ("DarkOliveGreen4", (0.43, 0.545, 0.24)),
    "darkorange1": ("DarkOrange1", (1, 0.498, 0)),
    "darkorange2": ("DarkOrange2", (0.932, 0.464, 0)),
    "darkorange3": ("DarkOrange3", (0.804, 0.4, 0)),
    "darkorange4": ("DarkOrange4", (0.545, 0.27, 0)),
    "darkorchid1": ("DarkOrchid1", (0.75, 0.244, 1)),
    "darkorchid2": ("DarkOrchid2", (0.698, 0.228, 0.932)),
    "darkorchid3": ("DarkOrchid3", (0.604, 0.196, 0.804)),
    "darkorchid4": ("DarkOrchid4", (0.408, 0.132, 0.545)),
    "darkseagreen1": ("DarkSeaGreen1", (0.756, 1, 0.756)),
    "darkseagreen2": ("DarkSeaGreen2", (0.705, 0.932, 0.705)),
    "darkseagreen3": ("DarkSeaGreen3", (0.608, 0.804, 0.608)),
    "darkseagreen4": ("DarkSeaGreen4", (0.41, 0.545, 0.41)),
    "darkslategray1": ("DarkSlateGray1", (0.592, 1, 1)),
    "darkslategray2": ("DarkSlateGray2", (0.552, 0.932, 0.932)),
    "darkslategray3": ("DarkSlateGray3", (0.475, 0.804, 0.804)),
    "darkslategray4": ("DarkSlateGray4", (0.32, 0.545, 0.545)),
    "deeppink1": ("DeepPink1", (1, 0.08, 0.576)),
    "deeppink2": ("DeepPink2", (0.932, 0.07, 0.536)),
    "deeppink3": ("DeepPink3", (0.804, 0.064, 0.464)),
    "deeppink4": ("DeepPink4", (0.545, 0.04, 0.312)),
    "deepskyblue1": ("DeepSkyBlue1", (0, 0.75, 1)),
    "deepskyblue2": ("DeepSkyBlue2", (0, 0.698, 0.932)),
    "deepskyblue3": ("DeepSkyBlue3", (0, 0.604, 0.804)),
    "deepskyblue4": ("DeepSkyBlue4", (0, 0.408, 0.545)),
    "dodgerblue1": ("DodgerBlue1", (0.116, 0.565, 1)),
    "dodgerblue2": ("DodgerBlue2", (0.11, 0.525, 0.932)),
    "dodgerblue3": ("DodgerBlue3", (0.094, 0.455, 0.804)),
    "dodgerblue4": ("DodgerBlue4", (0.064, 0.305, 0.545)),
    "firebrick1": ("Firebrick1", (1, 0.19, 0.19)),
    "firebrick2": ("Firebrick2", (0.932, 0.172, 0.172)),
    "firebrick3": ("Firebrick3", (0.804, 0.15, 0.15)),
    "firebrick4": ("Firebrick4", (0.545, 0.1, 0.1)),
    "gold1": ("Gold1", (1, 0.844, 0)),
    "gold2": ("Gold2", (0.932, 0.79, 0)),
    "gold3": ("Gold3", (0.804, 0.68, 0)),
    "gold4": ("Gold4", (0.545, 0.46, 0)),
    "goldenrod1": ("Goldenrod1", (1, 0.756, 0.145)),
    "goldenrod2": ("Goldenrod2", (0.932, 0.705, 0.132)),
    "goldenrod3": ("Goldenrod3", (0.804, 0.608, 0.112)),
    "goldenrod4": ("Goldenrod4", (0.545, 0.41, 0.08)),
    "green1": ("Green1", (0, 1, 0)),
    "green2": ("Green2", (0, 0.932, 0)),
    "green3": ("Green3", (0, 0.804, 0)),
    "green4": ("Green4", (0, 0.545, 0)),
    "honeydew1": ("Honeydew1", (0.94, 1, 0.94)),
    "honeydew2": ("Honeydew2", (0.88, 0.932, 0.88)),
    "honeydew3": ("Honeydew3", (0.756, 0.804, 0.756)),
    "honeydew4": ("Honeydew4", (0.512, 0.545, 0.512)),
    "hotpink1": ("HotPink1", (1, 0.43, 0.705)),
    "hotpink2": ("HotPink2", (0.932, 0.415, 0.655)),
    "hotpink3": ("HotPink3", (0.804, 0.376, 0.565)),
    "hotpink4": ("HotPink4", (0.545, 0.228, 0.385)),
    "indianred1": ("IndianRed1", (1, 0.415, 0.415)),
    "indianred2": ("IndianRed2", (0.932, 0.39, 0.39)),
    "indianred3": ("IndianRed3", (0.804, 0.332, 0.332)),
    "indianred4": ("IndianRed4", (0.545, 0.228, 0.228)),
    "ivory1": ("Ivory1", (1, 1, 0.94)),
    "ivory2": ("Ivory2", (0.932, 0.932, 0.88)),
    "ivory3": ("Ivory3", (0.804, 0.804, 0.756)),
    "ivory4": ("Ivory4", (0.545, 0.545, 0.512)),
    "khaki1": ("Khaki1", (1, 0.965, 0.56)),
    "khaki2": ("Khaki2", (0.932, 0.9, 0.52)),
    "khaki3": ("Khaki3", (0.804, 0.776, 0.45)),
    "khaki4": ("Khaki4", (0.545, 0.525, 0.305)),
    "lavenderblush1": ("LavenderBlush1", (1, 0.94, 0.96)),
    "lavenderblush2": ("LavenderBlush2", (0.932, 0.88, 0.898)),
    "lavenderblush3": ("LavenderBlush3", (0.804, 0.756, 0.772)),
    "lavenderblush4": ("LavenderBlush4", (0.545, 0.512, 0.525)),
    "lemonchiffon1": ("LemonChiffon1", (1, 0.98, 0.804)),
    "lemonchiffon2": ("LemonChiffon2", (0.932, 0.912, 0.75)),
    "lemonchiffon3": ("LemonChiffon3", (0.804, 0.79, 0.648)),
    "lemonchiffon4": ("LemonChiffon4", (0.545, 0.536, 0.44)),
    "lightblue1": ("LightBlue1", (0.75, 0.936, 1)),
    "lightblue2": ("LightBlue2", (0.698, 0.875, 0.932)),
    "lightblue3": ("LightBlue3", (0.604, 0.752, 0.804)),
    "lightblue4": ("LightBlue4", (0.408, 0.512, 0.545)),
    "lightcyan1": ("LightCyan1", (0.88, 1, 1)),
    "lightcyan2": ("LightCyan2", (0.82, 0.932, 0.932)),
    "lightcyan3": ("LightCyan3", (0.705, 0.804, 0.804)),
    "lightcyan4": ("LightCyan4", (0.48, 0.545, 0.545)),
    "lightgoldenrod1": ("LightGoldenrod1", (1, 0.925, 0.545)),
    "lightgoldenrod2": ("LightGoldenrod2", (0.932, 0.864, 0.51)),
    "lightgoldenrod3": ("LightGoldenrod3", (0.804, 0.745, 0.44)),
    "lightgoldenrod4": ("LightGoldenrod4", (0.545, 0.505, 0.298)),
    "lightpink1": ("LightPink1", (1, 0.684, 0.725)),
    "lightpink2": ("LightPink2", (0.932, 0.635, 0.68)),
    "lightpink3": ("LightPink3", (0.804, 0.55, 0.585)),
    "lightpink4": ("LightPink4", (0.545, 0.372, 0.396)),
    "lightsalmon1": ("LightSalmon1", (1, 0.628, 0.48)),
    "lightsalmon2": ("LightSalmon2", (0.932, 0.585, 0.448)),
    "lightsalmon3": ("LightSalmon3", (0.804, 0.505, 0.385)),
    "lightsalmon4": ("LightSalmon4", (0.545, 0.34, 0.26)),
    "lightskyblue1": ("LightSkyBlue1", (0.69, 0.888, 1)),
    "lightskyblue2": ("LightSkyBlue2", (0.644, 0.828, 0.932)),
    "lightskyblue3": ("LightSkyBlue3", (0.552, 0.712, 0.804)),
    "lightskyblue4": ("LightSkyBlue4", (0.376, 0.484, 0.545)),
    "lightsteelblue1": ("LightSteelBlue1", (0.792, 0.884, 1)),
    "lightsteelblue2": ("LightSteelBlue2", (0.736, 0.824, 0.932)),
    "lightsteelblue3": ("LightSteelBlue3", (0.635, 0.71, 0.804)),
    "lightsteelblue4": ("LightSteelBlue4", (0.43, 0.484, 0.545)),
    "lightyellow1": ("LightYellow1", (1, 1, 0.88)),
    "lightyellow2": ("LightYellow2", (0.932, 0.932, 0.82)),
    "lightyellow3": ("LightYellow3", (0.804, 0.804, 0.705)),
    "lightyellow4": ("LightYellow4", (0.545, 0.545, 0.48)),
    "magenta1": ("Magenta1", (1, 0, 1)),
    "magenta2": ("Magenta2", (0.932, 0, 0.932)),
    "magenta3": ("Magenta3", (0.804, 0, 0.804)),
    "magenta4": ("Magenta4", (0.545, 0, 0.545)),
    "maroon1": ("Maroon1", (1, 0.204, 0.7)),
    "maroon2": ("Maroon2", (0.932, 0.19, 0.655)),
    "maroon3": ("Maroon3", (0.804, 0.16, 0.565)),
    "maroon4": ("Maroon4", (0.545, 0.11, 0.385)),
    "mediumorchid1": ("MediumOrchid1", (0.88, 0.4, 1)),
    "mediumorchid2": ("MediumOrchid2", (0.82, 0.372, 0.932)),
    "mediumorchid3": ("MediumOrchid3", (0.705, 0.32, 0.804)),
    "mediumorchid4": ("MediumOrchid4", (0.48, 0.215, 0.545)),
    "mediumpurple1": ("MediumPurple1", (0.67, 0.51, 1)),
    "mediumpurple2": ("MediumPurple2", (0.624, 0.475, 0.932)),
    "mediumpurple3": ("MediumPurple3", (0.536, 0.408, 0.804)),
    "mediumpurple4": ("MediumPurple4", (0.365, 0.28, 0.545)),
    "mistyrose1": ("MistyRose1", (1, 0.894, 0.884)),
    "mistyrose2": ("MistyRose2", (0.932, 0.835, 0.824)),
    "mistyrose3": ("MistyRose3", (0.804, 0.716, 0.71)),
    "mistyrose4": ("MistyRose4", (0.545, 0.49, 0.484)),
    "navajowhite1": ("NavajoWhite1", (1, 0.87, 0.68)),
    "navajowhite2": ("NavajoWhite2", (0.932, 0.81, 0.63)),
    "navajowhite3": ("NavajoWhite3", (0.804, 0.7, 0.545)),
    "navajowhite4": ("NavajoWhite4", (0.545, 0.475, 0.37)),
    "olivedrab1": ("OliveDrab1", (0.752, 1, 0.244)),
    "olivedrab2": ("OliveDrab2", (0.7, 0.932, 0.228)),
    "olivedrab3": ("OliveDrab3", (0.604, 0.804, 0.196)),
    "olivedrab4": ("OliveDrab4", (0.41, 0.545, 0.132)),
    "orange1": ("Orange1", (1, 0.648, 0)),
    "orange2": ("Orange2", (0.932, 0.604, 0)),
    "orange3": ("Orange3", (0.804, 0.52, 0)),
    "orange4": ("Orange4", (0.545, 0.352, 0)),
    "orangered1": ("OrangeRed1", (1, 0.27, 0)),
    "orangered2": ("OrangeRed2", (0.932, 0.25, 0)),
    "orangered3": ("OrangeRed3", (0.804, 0.215, 0)),
    "orangered4": ("OrangeRed4", (0.545, 0.145, 0)),
    "orchid1": ("Orchid1", (1, 0.512, 0.98)),
    "orchid2": ("Orchid2", (0.932, 0.48, 0.912)),
    "orchid3": ("Orchid3", (0.804, 0.41, 0.79)),
    "orchid4": ("Orchid4", (0.545, 0.28, 0.536)),
    "palegreen1": ("PaleGreen1", (0.604, 1, 0.604)),
    "palegreen2": ("PaleGreen2", (0.565, 0.932, 0.565)),
    "palegreen3": ("PaleGreen3", (0.488, 0.804, 0.488)),
    "palegreen4": ("PaleGreen4", (0.33, 0.545, 0.33)),
    "paleturquoise1": ("PaleTurquoise1", (0.732, 1, 1)),
    "paleturquoise2": ("PaleTurquoise2", (0.684, 0.932, 0.932)),
    "paleturquoise3": ("PaleTurquoise3", (0.59, 0.804, 0.804)),
    "paleturquoise4": ("PaleTurquoise4", (0.4, 0.545, 0.545)),
    "palevioletred1": ("PaleVioletRed1", (1, 0.51, 0.67)),
    "palevioletred2": ("PaleVioletRed2", (0.932, 0.475, 0.624)),
    "palevioletred3": ("PaleVioletRed3", (0.804, 0.408, 0.536)),
    "palevioletred4": ("PaleVioletRed4", (0.545, 0.28, 0.365)),
    "peachpuff1": ("PeachPuff1", (1, 0.855, 0.725)),
    "peachpuff2": ("PeachPuff2", (0.932, 0.796, 0.68)),
    "peachpuff3": ("PeachPuff3", (0.804, 0.688, 0.585)),
    "peachpuff4": ("PeachPuff4", (0.545, 0.468, 0.396)),
    "pink1": ("Pink1", (1, 0.71, 0.772)),
    "pink2": ("Pink2", (0.932, 0.664, 0.72)),
    "pink3": ("Pink3", (0.804, 0.57, 0.62)),
    "pink4": ("Pink4", (0.545, 0.39, 0.424)),
    "plum1": ("Plum1", (1, 0.732, 1)),
    "plum2": ("Plum2", (0.932, 0.684, 0.932)),
    "plum3": ("Plum3", (0.804, 0.59, 0.804)),
    "plum4": ("Plum4", (0.545, 0.4, 0.545)),
    "purple1": ("Purple1", (0.608, 0.19, 1)),
    "purple2": ("Purple2", (0.57, 0.172, 0.932)),
    "purple3": ("Purple3", (0.49, 0.15, 0.804)),
    "purple4": ("Purple4", (0.332, 0.1, 0.545)),
    "red1": ("Red1", (1, 0, 0)),
    "red2": ("Red2", (0.932, 0, 0)),
    "red3": ("Red3", (0.804, 0, 0)),
    "red4": ("Red4", (0.545, 0, 0)),
    "rosybrown1": ("RosyBrown1", (1, 0.756, 0.756)),
    "rosybrown2": ("RosyBrown2", (0.932, 0.705, 0.705)),
    "rosybrown3": ("RosyBrown3", (0.804, 0.608, 0.608)),
    "rosybrown4": ("RosyBrown4", (0.545, 0.41, 0.41)),
    "royalblue1": ("RoyalBlue1", (0.284, 0.464, 1)),
    "royalblue2": ("RoyalBlue2", (0.264, 0.43, 0.932)),
    "royalblue3": ("RoyalBlue3", (0.228, 0.372, 0.804)),
    "royalblue4": ("RoyalBlue4", (0.152, 0.25, 0.545)),
    "salmon1": ("Salmon1", (1, 0.55, 0.41)),
    "salmon2": ("Salmon2", (0.932, 0.51, 0.385)),
    "salmon3": ("Salmon3", (0.804, 0.44, 0.33)),
    "salmon4": ("Salmon4", (0.545, 0.298, 0.224)),
    "seagreen1": ("SeaGreen1", (0.33, 1, 0.624)),
    "seagreen2": ("SeaGreen2", (0.305, 0.932, 0.58)),
    "seagreen3": ("SeaGreen3", (0.264, 0.804, 0.5)),
    "seagreen4": ("SeaGreen4", (0.18, 0.545, 0.34)),
    "seashell1": ("Seashell1", (1, 0.96, 0.932)),
    "seashell2": ("Seashell2", (0.932, 0.898, 0.87)),
    "seashell3": ("Seashell3", (0.804, 0.772, 0.75)),
    "seashell4": ("Seashell4", (0.545, 0.525, 0.51)),
    "sienna1": ("Sienna1", (1, 0.51, 0.28)),
    "sienna2": ("Sienna2", (0.932, 0.475, 0.26)),
    "sienna3": ("Sienna3", (0.804, 0.408, 0.224)),
    "sienna4": ("Sienna4", (0.545, 0.28, 0.15)),
    "skyblue1": ("SkyBlue1", (0.53, 0.808, 1)),
    "skyblue2": ("SkyBlue2", (0.494, 0.752, 0.932)),
    "skyblue3": ("SkyBlue3", (0.424, 0.65, 0.804)),
    "skyblue4": ("SkyBlue4", (0.29, 0.44, 0.545)),
    "slateblue1": ("SlateBlue1", (0.512, 0.435, 1)),
    "slateblue2": ("SlateBlue2", (0.48, 0.404, 0.932)),
    "slateblue3": ("SlateBlue3", (0.41, 0.35, 0.804)),
    "slateblue4": ("SlateBlue4", (0.28, 0.235, 0.545)),
    "slategray1": ("SlateGray1", (0.776, 0.888, 1)),
    "slategray2": ("SlateGray2", (0.725, 0.828, 0.932)),
    "slategray3": ("SlateGray3", (0.624, 0.712, 0.804)),
    "slategray4": ("SlateGray4", (0.424, 0.484, 0.545)),
    "snow1": ("Snow1", (1, 0.98, 0.98)),
    "snow2": ("Snow2", (0.932, 0.912, 0.912)),
    "snow3": ("Snow3", (0.804, 0.79, 0.79)),
    "snow4": ("Snow4", (0.545, 0.536, 0.536)),
    "springgreen1": ("SpringGreen1", (0, 1, 0.498)),
    "springgreen2": ("SpringGreen2", (0, 0.932, 0.464)),
    "springgreen3": ("SpringGreen3", (0, 0.804, 0.4)),
    "springgreen4": ("SpringGreen4", (0, 0.545, 0.27)),
    "steelblue1": ("SteelBlue1", (0.39, 0.72, 1)),
    "steelblue2": ("SteelBlue2", (0.36, 0.675, 0.932)),
    "steelblue3": ("SteelBlue3", (0.31, 0.58, 0.804)),
    "steelblue4": ("SteelBlue4", (0.21, 0.392, 0.545)),
    "tan1": ("Tan1", (1, 0.648, 0.31)),
    "tan2": ("Tan2", (0.932, 0.604, 0.288)),
    "tan3": ("Tan3", (0.804, 0.52, 0.248)),
    "tan4": ("Tan4", (0.545, 0.352, 0.17)),
    "thistle1": ("Thistle1", (1, 0.884, 1)),
    "thistle2": ("Thistle2", (0.932, 0.824, 0.932)),
    "thistle3": ("Thistle3", (0.804, 0.71, 0.804)),
    "thistle4": ("Thistle4", (0.545, 0.484, 0.545)),
    "tomato1": ("Tomato1", (1, 0.39, 0.28)),
    "tomato2": ("Tomato2", (0.932, 0.36, 0.26)),
    "tomato3": ("Tomato3", (0.804, 0.31, 0.224)),
    "tomato4": ("Tomato4", (0.545, 0.21, 0.15)),
    "turquoise1": ("Turquoise1", (0, 0.96, 1)),
    "turquoise2": ("Turquoise2", (0, 0.898, 0.932)),
    "turquoise3": ("Turquoise3", (0, 0.772, 0.804)),
    "turquoise4": ("Turquoise4", (0, 0.525, 0.545)),
    "violetred1": ("VioletRed1", (1, 0.244, 0.59)),
    "violetred2": ("VioletRed2", (0.932, 0.228, 0.55)),
    "violetred3": ("VioletRed3", (0.804, 0.196, 0.47)),
    "violetred4": ("VioletRed4", (0.545, 0.132, 0.32)),
    "wheat1": ("Wheat1", (1, 0.905, 0.73)),
    "wheat2": ("Wheat2", (0.932, 0.848, 0.684)),
    "wheat3": ("Wheat3", (0.804, 0.73, 0.59)),
    "wheat4": ("Wheat4", (0.545, 0.494, 0.4)),
    "yellow1": ("Yellow1", (1, 1, 0)),
    "yellow2": ("Yellow2", (0.932, 0.932, 0)),
    "yellow3": ("Yellow3", (0.804, 0.804, 0)),
    "yellow4": ("Yellow4", (0.545, 0.545, 0)),
    "gray0": ("Gray0", (0.745, 0.745, 0.745)),
    "green0": ("Green0", (0, 1, 0)),
    "grey0": ("Grey0", (0.745, 0.745, 0.745)),
    "maroon0": ("Maroon0", (0.69, 0.19, 0.376)),
    "purple0": ("Purple0", (0.628, 0.125, 0.94)),
}

_dvips_colors = {
    "greenyellow": ("GreenYellow", [0.85, 1, 0.31]),
    "yellow": ("Yellow", [1, 1, 0]),
    "goldenrod": ("Goldenrod", [1, 0.9, 0.16]),
    "dandelion": ("Dandelion", [1, 0.71, 0.16]),
    "apricot": ("Apricot", [1, 0.68, 0.48]),
    "peach": ("Peach", [1, 0.5, 0.3]),
    "melon": ("Melon", [1, 0.54, 0.5]),
    "yelloworange": ("YellowOrange", [1, 0.58, 0]),
    "orange": ("Orange", [1, 0.39, 0.13]),
    "burntorange": ("BurntOrange", [1, 0.49, 0]),
    "bittersweet": ("Bittersweet", [0.76, 0.19, 0.0]),
    "redorange": ("RedOrange", [1, 0.23, 0.13]),
    "mahogany": ("Mahogany", [0.65, 0.0975, 0.0845]),
    "maroon": ("Maroon", [0.68, 0.0884, 0.2176],),
    "brickred": ("BrickRed", [0.72, 0.0792, 0.0432]),
    "red": ("Red", [1, 0, 0]),
    "orangered": ("OrangeRed", [1, 0, 0.5]),
    "rubinered": ("RubineRed", [1, 0, 0.87]),
    "wildstrawberry": ("WildStrawberry", [1, 0.04, 0.61]),
    "salmon": ("Salmon", [1, 0.47, 0.62]),
    "carnationpink": ("CarnationPink", [1, 0.37, 1]),
    "magenta": ("Magenta", [1, 0, 1]),
    "violetred": ("VioletRed", [1, 0.19, 1]),
    "rhodamine": ("Rhodamine", [1, 0.18, 1]),
    "mulberry": ("Mulberry", [0.6468, 0.098, 0.98]),
    "redviolet": ("RedViolet", [0.6138, 0.066, 0.66],),
    "fuchsia": ("Fuchsia", [0.4876, 0.0828, 0.92]),
    "lavender": ("Lavender", [1, 0.52, 1]),
    "thistle": ("Thistle", [0.88, 0.41, 1]),
    "orchid": ("Orchid", [0.68, 0.36, 1]),
    "darkorchid": ("DarkOrchid", [0.6, 0.2, 0.8]),
    "purple": ("Purple", [0.55, 0.14, 1]),
    "plum": ("Plum", [0.5, 0, 1]),
    "violet": ("Violet", [0.21, 0.12, 1]),
    "royalpurple": ("RoyalPurple", [0.25, 0.1, 1]),
    "blueviolet": ("BlueViolet", [0.1344, 0.0864, 0.96]),
    "periwinkle": ("Periwinkle", [0.43, 0.45, 1]),
    "cadetblue": ("CadetBlue", [0.38, 0.43, 0.77]),
    "cornflowerblue": ("CornflowerBlue", [0.35, 0.87, 1]),
    "midnightblue": ("MidnightBlue", [0.0114, 0.4959, 0.57],),
    "navyblue": ("NavyBlue", [0.06, 0.46, 1]),
    "royalblue": ("RoyalBlue", [0, 0.5, 1]),
    "blue": ("Blue", [0, 0, 1]),
    "cerulean": ("Cerulean", [0.06, 0.89, 1]),
    "cyan": ("Cyan", [0, 1, 1]),
    "processblue": ("ProcessBlue", [0.04, 1, 1]),
    "skyblue": ("SkyBlue", [0.38, 1, 0.88]),
    "turquoise": ("Turquoise", [0.15, 1, 0.8]),
    "tealblue": ("TealBlue", [0.1372, 0.98, 0.6468]),
    "aquamarine": ("Aquamarine", [0.18, 1, 0.7]),
    "bluegreen": ("BlueGreen", [0.15, 1, 0.67]),
    "emerald": ("Emerald", [0, 1, 0.5]),
    "junglegreen": ("JungleGreen", [0.01, 1, 0.48]),
    "seagreen": ("SeaGreen", [0.31, 1, 0.5]),
    "green": ("Green", [0, 1, 0]),
    "forestgreen": ("ForestGreen", [0.0792, 0.88, 0.1056]),
    "pinegreen": ("PineGreen", [0.06, 0.75, 0.3075]),
    "limegreen": ("LimeGreen", [0.5, 1, 0]),
    "yellowgreen": ("YellowGreen", [0.56, 1, 0.26]),
    "springgreen": ("SpringGreen", [0.74, 1, 0.24]),
    "olivegreen": ("OliveGreen", [0.216, 0.6, 0.03]),
    "rawsienna": ("RawSienna", [0.55, 0.154, 0.0]),
    "sepia": ("Sepia", [0.3, 0.051, 0.0]),
    "brown": ("Brown", [0.4, 0.076, 0.0]),
    "tan": ("Tan", [0.86, 0.58, 0.44]),
    "gray": ("Gray", [0.5, 0.5, 0.5]),
    "black": ("Black", [0, 0, 0]),
    "white": ("White", [1, 1, 1]),
}
