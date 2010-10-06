"""Tools for manipulating integer sequences. """

import webbrowser

_oeis_url = 'http://www.research.att.com/~njas/sequences/?q='

def oeis(seq, tab=True):
    """The On-Line Encyclopedia of Integer Sequences. """
    url = _oeis_url + ','.join(map(str, seq))

    if tab:
        webbrowser.open_new_tab(url)
    else:
        webbrowser.open_new(url)
