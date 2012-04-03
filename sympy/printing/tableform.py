from sympy.core.containers import Tuple

class TableForm(object):
    """
    Allows to create nice table representation of data.

    Example::

        >>> from sympy import TableForm
        >>> t = TableForm([[5, 7], [4, 2], [10, 3]])

    Then use "print t" to print the table. You can use the SymPy's printing
    system to produce tables in any format (ascii, latex, html, ...).

    """

    def __new__(cls, data, headings=None, alignment="left",
            column_formats=None):
        """
        Creates a TableForm.

        Parameters:

            data ... 2D data to be put into the table
            headings ... gives the labels for entries in each dimension:
                         None ... no labels in any dimension
                         "automatic" ... gives successive integer labels
                         [[l1, l2, ...], ...] gives labels for each entry in
                             each dimension (can be None for some dimension)
            alignment ... alignment of the columns controlled by values;
                          "left" or "<"
                          "center" or "^"
                          "right" or ">"

        Example:

        >>> from sympy import TableForm
        >>> t = TableForm([[5, 7], [4, 2], [10, 3]])

        """
        # We only support 2D data. Check the consistency:
        from sympy.core.numbers import Integer
        _w = Integer(len(data[0]))
        _h = Integer(len(data))
        for line in data:
            assert len(line) == _w
        _lines = Tuple(*data)

        if headings is None:
            _headings = [None, None]
        elif headings == "automatic":
            _headings = [range(1, _h + 1), range(1, _w + 1)]
        else:
            h1, h2 = headings
            if h1 == "automatic":
                h1 = range(1, _h + 1)
            if h2 == "automatic":
                h2 = range(1, _w + 1)
            _headings = [h1, h2]

        _alignment = {'<': 'left', '>': 'right', '^': 'center'}.get(alignment, alignment)
        if column_formats:
            _column_formats = column_formats
        else:
            _column_formats = [None]*_w

        obj = object.__new__(cls)
        obj._w = _w
        obj._h = _h
        obj._lines = _lines
        obj._headings = _headings
        obj._alignment = _alignment
        obj._column_formats = _column_formats
        return obj

    def __repr__(self):
        from sstr import sstr
        return sstr(self, order=None)

    def __str__(self):
        from sstr import sstr
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

    def _latex(self, p):
        """
        Returns the string representation of 'self'.

        Example:

        >>> from sympy import TableForm
        >>> t = TableForm([[5, 7], [4, 2], [10, 3]])
        >>> s = t.as_latex()

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

        s = r"\begin{tabular}{" + " ".join(["c" for x in self._lines[0]]) + "}\n"
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
                    d.append(str(x))
            if self._headings[0]:
                d = [self._headings[0][i]] + d
            s += " & ".join(d) + r" \\" + "\n"
        s += r"\end{tabular}"
        return s
