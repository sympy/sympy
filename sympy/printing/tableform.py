from sympy.core.containers import Tuple

class TableForm(object):
    """
    Create a nice table representation of data.

    Example::

        >>> from sympy import TableForm
        >>> t = TableForm([[5, 7], [4, 2], [10, 3]])

    Then use "print t" to print the table. You can use the SymPy's printing
    system to produce tables in any format (ascii, latex, html, ...).

    """

    def __init__(self, data, **kwarg):
        """
        Creates a TableForm.

        Parameters:

            data ...        2D data to be put into the table

            headings ...    gives the labels for entries in each dimension:
                            None ... no labels in any dimension
                            "automatic" ... gives successive integer labels
                            [[l1, l2, ...], ...] gives labels for each entry in
                                each dimension (can be None for some dimension)
                            [default: (None, None)]

            alignment ...   alignment of the columns controlled by values;
                            "left" or "<"
                            "center" or "^"
                            "right" or ">"
                            [default: left]

            wipe_zeros ...  Wipe zeros, don't show them in the table
                            [default: True].

        Example:

        >>> from sympy import TableForm
        >>> t = TableForm([[5, 7], [4, 2], [10, 3]])

        """
        # We only support 2D data. Check the consistency:
        _w = len(data[0])
        _h = len(data)
        assert all(len(line) == _w for line in data)
        _lines = Tuple(*data)

        headings = kwarg.get("headings", [None, None])
        if headings == "automatic":
            _headings = [range(1, _h + 1), range(1, _w + 1)]
        else:
            h1, h2 = headings
            if h1 == "automatic":
                h1 = range(1, _h + 1)
            if h2 == "automatic":
                h2 = range(1, _w + 1)
            _headings = [h1, h2]

        alignment = kwarg.get("alignment", "left")
        _alignment = {'<': 'left', '>': 'right', '^': 'center'}.get(alignment, alignment)

        _column_formats = kwarg.get("column_formats", [None]*_w)

        _wipe_zeros = kwarg.get("wipe_zeros", True)

        self._w = _w
        self._h = _h
        self._lines = _lines
        self._headings = _headings
        self._alignment = _alignment
        self._column_formats = _column_formats
        self._wipe_zeros = _wipe_zeros

    def __repr__(self):
        from str import sstr
        return sstr(self, order=None)

    def __str__(self):
        from str import sstr
        return sstr(self, order=None)

    def as_str(self):
        # obsolated ?
        return str(self)

    def as_latex(self):
        from latex import latex
        return latex(self)

    def _sympystr(self, p):
        """
        Returns the string representation of 'self'.

        Example:

        >>> from sympy import TableForm
        >>> t = TableForm([[5, 7], [4, 2], [10, 3]])
        >>> s = t.as_str()

        """
        column_widths = [0] * self._w
        lines = []
        for line in self._lines:
            new_line = []
            for i in range(self._w):
                # Format the item somehow if needed:
                s = str(line[i])
                if self._wipe_zeros and (s == "0"):
                    s = " "
                w = len(s)
                if w > column_widths[i]:
                    column_widths[i] = w
                new_line.append(s)
            lines.append(new_line)

        # Check heading:
        if self._headings[1]:
            new_line = []
            for i in range(self._w):
                # Format the item somehow if needed:
                s = str(self._headings[1][i])
                w = len(s)
                if w > column_widths[i]:
                    column_widths[i] = w
                new_line.append(s)
            self._headings[1] = new_line

        format_str = []
        for w in column_widths:
            if self._alignment == "left":
                align = "-"
            elif self._alignment == "right":
                align = ""
            elif self._alignment == "center":
                align = ""
            else:
                raise NotImplementedError()
            format_str += ["%" + align + str(w) + "s"]
        format_str = ' '.join(format_str)
        format_str += "\n"

        if self._headings[0]:
            self._headings[0] = [str(x) for x in self._headings[0]]
            heading_width = max([len(x) for x in self._headings[0]])
            format_str = "%" + str(heading_width) + "s | " + format_str

        s = []
        if self._headings[1]:
            d = self._headings[1]
            if self._headings[0]:
                d = [""] + d
            first_line = format_str % tuple(d)
            s.append(first_line)
            s.append("-" * (len(first_line) - 2) + "\n")
        for i, line in enumerate(lines):
            d = [l if self._alignment != 'center' else
                 l.center(column_widths[j]) for j,l in enumerate(line)]
            if self._headings[0]:
                d = [self._headings[0][i]] + d
            s.append(format_str % tuple(d))
        return ''.join(s)

    def _latex(self, printer):
        """
        Returns the string representation of 'self'.
        """
        # Check heading:
        if self._headings[1]:
            new_line = []
            for i in range(self._w):
                # Format the item somehow if needed:
                s = str(self._headings[1][i])
                w = len(s)
                new_line.append(s)
            self._headings[1] = new_line

        if self._headings[0]:
            self._headings[0] = [str(x) for x in self._headings[0]]

        align_char = {'left': 'l', 'right':'r', 'center':'c'}.get(self._alignment)
        align_list = [align_char] * len(self._lines[0])
        if self._headings[0]:
            align_list = [align_char] + align_list
        s = r"\begin{tabular}{" + " ".join(align_list) + "}\n"

        if self._headings[1]:
            d = self._headings[1]
            if self._headings[0]:
                d = [""] + d
            first_line = " & ".join(d) + r" \\" + "\n"
            s += first_line
            s += r"\hline" + "\n"
        for i, line in enumerate(self._lines):
            d = []
            for j, x in enumerate(line):
                if self._column_formats[j]:
                    d.append(self._column_formats[j] % x)
                else:
                    v = printer._print(x)
                    if self._wipe_zeros and (v == "0"):
                        d.append(" ")
                    else:
                        d.append("$%s$" % v)
            if self._headings[0]:
                d = [self._headings[0][i]] + d
            s += " & ".join(d) + r" \\" + "\n"
        s += r"\end{tabular}"
        return s
