"""
Pygments styles used for syntax highlighting.

These are based on the Sphinx style (see
https://github.com/sphinx-doc/sphinx/blob/master/sphinx/pygments_styles.py)
for light mode and the Friendly style for dark mode.

The styles here have been adjusted so that they are WCAG AA compatible. The
tool at https://github.com/mpchadwick/pygments-high-contrast-stylesheets was
used to identify colors that should be adjusted.

"""
from pygments.style import Style
from pygments.styles.friendly import FriendlyStyle
from pygments.styles.native import NativeStyle
from pygments.token import Comment, Generic, Literal, Name, Number, Text

class SphinxHighContrastStyle(Style):
    """
    Like Sphinx (which is like friendly, but a bit darker to enhance contrast
    on the green background) but with higher contrast colors.

    """

    @property
    def _pre_style(self):
        # This is used instead of the default 125% so that multiline Unicode
        # pprint output looks good
        return 'line-height: 120%;'

    background_color = '#eeffcc'
    default_style = ''

    styles = FriendlyStyle.styles
    styles.update({
        # These are part of the Sphinx modification to "friendly"
        Generic.Output: '#333',
        Number: '#208050',

        # These are adjusted from "friendly" (Comment is adjusted from
        # "sphinx") to have better color contrast against the background.
        Comment: 'italic #3c7a88',
        Comment.Hashbang: 'italic #3c7a88',
        Comment.Multiline: 'italic #3c7a88',
        Comment.PreprocFile: 'italic #3c7a88',
        Comment.Single: 'italic #3c7a88',
        Comment.Special: '#3a7784 bg:#fff0f0',
        Generic.Error: '#e60000',
        Generic.Inserted: '#008200',
        Generic.Prompt: 'bold #b75709',
        Name.Class: 'bold #0e7ba6',
        Name.Constant: '#2b79a1',
        Name.Entity: 'bold #c54629',
        Name.Namespace: 'bold #0e7ba6',
        Name.Variable: '#ab40cd',
        Text.Whitespace: '#707070',
        Literal.String.Interpol: 'italic #3973b7',
        Literal.String.Other: '#b75709',
        Name.Variable.Class: '#ab40cd',
        Name.Variable.Global: '#ab40cd',
        Name.Variable.Instance: '#ab40cd',
        Name.Variable.Magic: '#ab40cd',
    })



class NativeHighContrastStyle(NativeStyle):
    """
    Like native, but with higher contrast colors.
    """
    @property
    def _pre_style(self):
        # This is used instead of the default 125% so that multiline Unicode
        # pprint output looks good
        return 'line-height: 120%;'

    styles = NativeStyle.styles

    # These are adjusted to have better color contrast against the background
    styles.update({
        Comment.Preproc: 'bold #e15a5a',
        Comment.Special: 'bold #f75050 bg:#520000',
        Generic.Deleted: '#e75959',
        Generic.Error: '#e75959',
        Generic.Traceback: '#e75959',
        Literal.Number: '#438dc4',
        Name.Builtin: '#2594a1',
        # We also remove the underline here from the original style
        Name.Class: '#548bd3',
        Name.Function: '#548bd3',
        # We also remove the underline here from the original style
        Name.Namespace: '#548bd3',
        Text.Whitespace: '#878787',
        Literal.Number.Bin: '#438dc4',
        Literal.Number.Float: '#438dc4',
        Literal.Number.Hex: '#438dc4',
        Literal.Number.Integer: '#438dc4',
        Literal.Number.Oct: '#438dc4',
        Name.Builtin.Pseudo: '#2594a1',
        Name.Function.Magic: '#548bd3',
        Literal.Number.Integer.Long: '#438dc4',
    })
