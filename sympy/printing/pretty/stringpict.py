"""
Prettyprinter by Jurjen Bos.

Classes used to render mathematical expressions as text-based pictures.
"""

import shutil

from .pretty_symbology import hobj, vobj, xsym, xobj, pretty_use_unicode, line_width, center
from sympy.utilities.exceptions import sympy_deprecation_warning

_GLOBAL_WRAP_LINE = None


class stringPict:
    """An ASCII picture.
    The pictures are represented as a list of equal length strings.
    """

    #special value for stringPict.below
    LINE = 'line'

    def __init__(self, s, baseline=0):
        """Initialize from string.
        Multiline strings are centered.
        """
        assert isinstance(s, str)
        assert isinstance(baseline, int)
        self.s = s
        #picture is a string that just can be printed
        self.picture = stringPict.equalLengths(s.splitlines())
        #baseline is the line number of the "base line"
        self.baseline = baseline
        self.binding = None

    @staticmethod
    def equalLengths(lines):
        # empty lines
        if not lines:
            return ['']

        width = max(line_width(line) for line in lines)
        return [center(line, width) for line in lines]

    def height(self):
        """The height of the picture in characters."""
        return len(self.picture)

    def width(self):
        """The width of the picture in characters."""
        return line_width(self.picture[0])

    @staticmethod
    def next(*args):
        """Put a string of stringPicts next to each other.
        Returns string, baseline arguments for stringPict.
        """
        #convert everything to stringPicts
        objects = []
        for arg in args:
            if isinstance(arg, str):
                arg = stringPict(arg)
            objects.append(arg)

        #make a list of pictures, with equal height and baseline
        newBaseline = max((obj.baseline for obj in objects))
        newHeightBelowBaseline = max(obj.height() - obj.baseline for obj in objects)
        newHeight = newBaseline + newHeightBelowBaseline

        pictures = []
        for obj in objects:
            oneEmptyLine = [' ' * obj.width()]
            basePadding = newBaseline - obj.baseline
            totalPadding = newHeight - obj.height()
            pictures.append(
                oneEmptyLine * basePadding
                + obj.picture
                + oneEmptyLine * (totalPadding - basePadding)
            )

        result = [''.join(lines) for lines in zip(*pictures)]
        return '\n'.join(result), newBaseline


    def v_align(self, *args, align="t"):
        """
        Align a list of elements, keeping
        the baseline relative to self.
        """
        objects = []
        for arg in args:
            if isinstance(arg, str):
                arg = self.__class__(arg)
            objects.append(arg)

        baseline = self.baseline
        height = self.height()
        max_height = max((height, max(arg.height() for obj in objects)))
        texts = (str(obj) for obj in objects)

        if align=="b":
            if max_height > height:
                baseline = baseline + max_height - height

            offsets = (max_height-obj.height() for obj in objects)
            texts = (offset*"\n"+ text if offset>0 else text for offset, obj, text in  zip(offsets, objects, texts))
        elif align=="c":
            if max_height > height:
                baseline = baseline + (max_height - height+1)//2
            offsets = ( (max_height-obj.height()+1)//2 for obj in objects)
            texts = (offset*"\n"+ text if offset>0 else text for offset, obj, text in  zip(offsets, objects, texts))
        elif align!="t":
            raise ValueError()

        result = [self.__class__(text, baseline) for text in texts]
        return result

    def right(self, *args, align=""):
        r"""Put pictures next to this one.
        Returns string, baseline arguments for stringPict.
        (Multiline) strings are allowed, and are given a baseline of 0.

        Examples
        ========

        >>> from sympy.printing.pretty.stringpict import stringPict
        >>> print(stringPict("10").right(" + ",stringPict("1\r-\r2",1)))
             1
        10 + -
             2

        """

        if align=="":
            return self.__class__(*stringPict.next(self, *args))

        return self.__class__(*stringPict.next(*self.v_align( *((self,) + args), align=align)))

    def left(self, *args, align=""):
        """
        Put pictures (left to right) at left.
        Returns string, baseline arguments for stringPict.

        The parameter `align` controls the relative alignments
        of the elements. By default (`align=""`) aligns the
        elements regarding the base line.

        The other possible values for the `align` parameter
        are

            - "t": tops are aligned
            - "c": centers are aligned
            - "b": bottoms are aligned

        In all the cases, the resulting baseline is localized
        on to the `self` baseline.

        Examples
        ========

        >>> from sympy.printing.pretty.stringpict import stringPict
        >>> stringPict('a\\n-\\nb', 1).left(
        ...     "numerator-> ",align="t").right(" <-denominator",
        ...      align="b").parenthesis("","}").right(" fraction")
        numerator-> a              \\
                    -              > fraction
                    b <-denominator/
        """

        if align=="":
            return self.__class__(*stringPict.next(*(args + (self,))))
        return self.__class__(
            *stringPict.next(
                *self.v_align( *(args + (self,)), align=align)
            )
        )

    @staticmethod
    def stack(*args, align="c"):
        """Put pictures on top of each other,
        from top to bottom.
        Returns string, baseline arguments for stringPict.
        The baseline is the baseline of the second picture.
        By default, everything is centered. If `align` is set to
        'l' ('r') the elements are aligned to the  left (right) margin.
        Baseline is the baseline of the second picture.
        Strings are allowed.
        The special value stringPict.LINE is a row of '-' extended to the width.
        """
        #convert everything to stringPicts; keep LINE
        objects = []
        for arg in args:
            if arg is not stringPict.LINE and isinstance(arg, str):
                arg = stringPict(arg)
            objects.append(arg)

        #compute new width
        newWidth = max(obj.width() for obj in objects if obj is not stringPict.LINE)

        lineObj = stringPict(hobj('-', newWidth))

        #replace LINE with proper lines
        for i, obj in enumerate(objects):
            if obj is stringPict.LINE:
                objects[i] = lineObj

        #stack the pictures, and center the result
        if align == "c":
            newPicture = [
                center(line, newWidth) for obj in objects for line in obj.picture
            ]
        elif align == "l":
            newPicture = [
                line + (newWidth - len(line)) * ' '
                for obj in objects
                for line in obj.picture
            ]
        elif align == "r":
            newPicture = [
                (newWidth - len(line)) * ' ' + line
                for obj in objects
                for line in obj.picture
            ]
        else:
            raise ValueError(
                f"the align parameter must be one of 'l'(left), 'r'(right) or 'c' (center). Got {align}."
            )

        newBaseline = objects[0].height() + objects[1].baseline
        return "\n".join(newPicture), newBaseline

    def below(self, *args, align="c"):
        """Put pictures under this picture.
        Returns string, baseline arguments for stringPict.
        Baseline is baseline of top picture

        Examples
        ========

        Let' s consider first how to draw a fraction::

            >>> from sympy.printing.pretty.stringpict import stringPict
            >>> print(stringPict("x+3").below(
            ...       stringPict.LINE, '3')) #doctest: +NORMALIZE_WHITESPACE
            x+3
            ---
             3

        The optional argument `align` controls the alignment.
        The default value is `'c'`. For `align='l'`,
        the text is left aligned::

            >>> print(stringPict("Humpty Dumpty").below(
            ...       stringPict("falls off the wall..."), align="l"))
            Humpty Dumpty
            falls off the wall...

        and for `align='r'`, the elements are right aligned::

            >>> print(stringPict("Humpty Dumpty").below(
            ...       stringPict("falls off the wall..."), align="r"))
                    Humpty Dumpty
            falls off the wall...

        """
        s, baseline = stringPict.stack(self, *args, align=align)
        return self.__class__(s, self.baseline)

    def above(self, *args, align="c"):
        """
        Put pictures above this picture.

        Returns string, baseline arguments for stringPict.
        Baseline is baseline of bottom picture.
        the horizontal alignment is controlled by the keyword
        parameter `align`. See stringPict.below() for details and
        examples.
        """
        string, baseline = stringPict.stack(*(args + (self,)), align=align)
        baseline = len(string.splitlines()) - self.height() + self.baseline
        return self.__class__(string, baseline)

    def parens(self, left="(", right=")", ifascii_nougly=False):
        """Put parentheses around self.
        Returns string, baseline arguments for stringPict.

        left or right can be None or empty string which means 'no paren from
        that side'.
        """

        # TODO: Add a deprecation warning in favor of
        # `self.parenthesis(...)`?

        h = self.height()
        b = self.baseline

        # XXX this is a hack -- ascii parens are ugly!
        if ifascii_nougly and not pretty_use_unicode():
            h = 1
            b = 0

        res = self

        if left:
            lparen = stringPict(vobj(left, h), baseline=b)
            res = lparen.right(self)
        if right:
            rparen = stringPict(vobj(right, h), baseline=b)
            res = res.right(rparen)

        return "\n".join(res.picture), res.baseline

    def parenthesis(self, left="(", right=")", ifascii_nougly=False):
        """Returns a new object of the same type than self,
        sorrounded by parenthesis of the type specified by
        the arguments `left` and `right`.

        `left` or `right` can be None or empty string which means 'no paren from
        that side'
        """
        h = self.height()
        b = self.baseline

        # XXX this is a hack -- ascii parens are ugly!
        if ifascii_nougly and not pretty_use_unicode():
            h = 1
            b = 0

        res = self

        if left:
            lparen = stringPict(vobj(left, h), baseline=b)
            res = lparen.right(self)
        if right:
            rparen = stringPict(vobj(right, h), baseline=b)
            res = res.right(rparen)

        return self.__class__("\n".join(res.picture), res.baseline)

    def leftslash(self):
        """Precede object by a slash of the proper size."""
        # XXX not used anywhere ?
        height = max(self.baseline, self.height() - 1 - self.baseline) * 2 + 1
        slash = '\n'.join(
            ' ' * (height - i - 1) + xobj('/', 1) + ' ' * i for i in range(height)
        )
        return self.left(stringPict(slash, height // 2))

    def root(self, n=None):
        """
        Produce a nice root symbol.
        Produces ugly results for big n inserts.

        Examples
        ========

        >>> from sympy.printing.pretty.stringpict import stringPict, prettyForm
        >>> print(stringPict("x+3").root().right(" + a"))
          _____
        \\/ x+3 + a

        >>> print(stringPict("x+3").root(stringPict("3")).right(" + a"))
        3 _____
        \\/ x+3 + a

        >>> print((prettyForm("x")**stringPict("a")).root().right(" + a"))
           ____
          /  a
        \\/  x + a

        >>> print((prettyForm("x")**stringPict("a")).root(stringPict("3")).right(" + a"))
           ____
        3 /  a
        \\/  x + a

        >>> print((prettyForm("x+3")/prettyForm("y")).root().right(" + a"))
            _____
           / x+3
          /  --- + a
        \\/   y

        >>> print((prettyForm("x+3")/prettyForm("y")).root(stringPict("3")).right(" + a"))
            _____
           / x+3
        3 /  --- + a
        \\/    y

        For indices with more than one line, use the Pow form:

        >>> print((prettyForm("x+3")/prettyForm("y")).root(
        ...       prettyForm("3")/prettyForm("5")).right(" + a"))
                 /3\\
             1 / |-|
                 \\5/
        /x+3\\
        |---|         + a
        \\ y /
        """
        # Decide if using a square root symbol or
        # an base - exponent form:
        if n is not None:
            if isinstance(n, str):
                n = n.ljust(2)
                n = stringPict(n)
            elif n.width()<2:
                n = stringPict(str(n).ljust(2))
            if n.height() > 1:
                exponent =  n.parenthesis().left(stringPict("1 / "), align="c")
                return self ** exponent

        # put line over expression
        result = self.above(hobj('_', 2 + self.width()))
        #construct right half of root symbol
        height = self.height()
        _zZ = xobj('/', 1)
        root_sign = prettyForm(xobj('\\', 1)+ _zZ)
        if n is not None:
            root_sign = root_sign.above(n, align="r")
        if height>1:
            slash = '\n'.join(' ' * (height - i - 2) + _zZ + ' ' * i for i in range(height-1))
            # TODO: To improve the use of the space, consider
            # using a vertical line instead '/', like
            #    -
            #   |x
            # 20|-
            #  \|2
            #
            # # remove the `.ljust.(2)` in `n` and
            # # replace the previous line by
            # _zZ = xobj('|', 1)
            # slash = "\n".join(height*[_zZ])
            #
            # but this requires to change many tests.
            slash = stringPict(" ").above(slash, align='l')
            root_sign = root_sign.right(slash, align="b")

        return result.left(root_sign, align="b")

    def render(self, *args, **kwargs):
        """Return the string form of self.

           Unless the argument line_break is set to False, it will
           break the expression in a form that can be printed
           on the terminal without being broken up.
         """
        if _GLOBAL_WRAP_LINE is not None:
            kwargs["wrap_line"] = _GLOBAL_WRAP_LINE

        if kwargs["wrap_line"] is False:
            return '\n'.join(self.picture)

        if kwargs["num_columns"] is not None:
            # Read the argument num_columns if it is not None
            ncols = kwargs["num_columns"]
        else:
            # Attempt to get a terminal width
            ncols = self.terminal_width()

        if ncols <= 0:
            ncols = 80

        # If smaller than the terminal width, no need to correct
        if self.width() <= ncols:
            return type(self.picture[0])(self)

        """
        Break long-lines in a visually pleasing format.
        without overflow indicators | with overflow indicators
        |   2  2        3     |     |   2  2        3    ↪|
        |6*x *y  + 4*x*y  +   |     |6*x *y  + 4*x*y  +  ↪|
        |                     |     |                     |
        |     3    4    4     |     |↪      3    4    4   |
        |4*y*x  + x  + y      |     |↪ 4*y*x  + x  + y    |
        |a*c*e + a*c*f + a*d  |     |a*c*e + a*c*f + a*d ↪|
        |*e + a*d*f + b*c*e   |     |                     |
        |+ b*c*f + b*d*e + b  |     |↪ *e + a*d*f + b*c* ↪|
        |*d*f                 |     |                     |
        |                     |     |↪ e + b*c*f + b*d*e ↪|
        |                     |     |                     |
        |                     |     |↪ + b*d*f            |
        """

        overflow_first = ""
        if kwargs["use_unicode"] or pretty_use_unicode():
            overflow_start = "\N{RIGHTWARDS ARROW WITH HOOK} "
            overflow_end = " \N{RIGHTWARDS ARROW WITH HOOK}"
        else:
            overflow_start = "> "
            overflow_end = " >"

        def chunks(line):
            """Yields consecutive chunks of line_width ncols"""
            prefix = overflow_first
            width, start = line_width(prefix + overflow_end), 0
            for i, x in enumerate(line):
                wx = line_width(x)
                # Only flush the screen when the current character overflows.
                # This way, combining marks can be appended even when width == ncols.
                if width + wx > ncols:
                    yield prefix + line[start:i] + overflow_end
                    prefix = overflow_start
                    width, start = line_width(prefix + overflow_end), i
                width += wx
            yield prefix + line[start:]

        # Concurrently assemble chunks of all lines into individual screens
        pictures = zip(*map(chunks, self.picture))

        # Join lines of each screen into sub-pictures
        pictures = ['\n'.join(picture) for picture in pictures]

        # Add spacers between sub-pictures
        return "\n\n".join(pictures)

    def subindex(self, sub_index):
        """Add a `subindex`.

        Examples
        ========

        >>> from sympy.printing.pretty.stringpict import stringPict
        >>> print( (stringPict("a").below("-").below(stringPict("b"))
        ...         ).parenthesis().subindex(stringPict("a=4")))
        /a\\
        |-|
        \\b/
            a=4
        """
        width_self = self.width()
        width_index = sub_index.width()
        extended_self = self.right(stringPict(width_index * ' '))
        sub_index = sub_index.left(stringPict(width_self * ' '))
        return extended_self.below(sub_index)

    def subsuperindices(self, sub_index, super_index):
        """
        Add  sub and super indices.

        Examples
        ========

        >>> from sympy.printing.pretty.stringpict import stringPict
        >>> print(stringPict("T").subsuperindices(
        ...       stringPict("a"),stringPict("b,c")))
         b,c
        T
         a
        """
        width_self = self.width()
        width_index = max(sub_index.width() ,super_index.width())
        extended_self = self.right(stringPict(width_index * ' '))
        super_index = super_index.left(stringPict(width_self * ' '))
        sub_index = sub_index.left(stringPict(width_self * ' '))
        return extended_self.above(super_index,
                                   align="l").below(
                                       sub_index,
                                       align="l")

    def superindex(self, super_index):
        """
        Add a `super index`.

        Examples
        ========

        >>> from sympy.printing.pretty.stringpict import stringPict
        >>> print(stringPict("e").superindex(stringPict("-s(x)")))
         -s(x)
        e

        """
        width_self = self.width()
        width_index = super_index.width()
        extended_self = self.right(stringPict(width_index * ' '))
        super_index = super_index.left(stringPict(width_self * ' '))
        return extended_self.above(super_index)

    def terminal_width(self):
        """Return the terminal width if possible, otherwise return 0."""
        size = shutil.get_terminal_size(fallback=(0, 0))
        return size.columns

    def __eq__(self, o):
        if isinstance(o, str):
            return '\n'.join(self.picture) == o
        elif isinstance(o, stringPict):
            return o.picture == self.picture
        return False

    def __hash__(self):
        return super().__hash__()

    def __str__(self):
        return '\n'.join(self.picture)

    def __repr__(self):
        return "stringPict(%r,%d)" % ('\n'.join(self.picture), self.baseline)

    def __getitem__(self, index):
        return self.picture[index]

    def __len__(self):
        return len(self.s)


class prettyForm(stringPict):
    """
    Extension of the stringPict class that knows about basic math applications,
    optimizing double minus signs.

    "Binding" is interpreted as follows::

        ATOM this is an atom: never needs to be parenthesized
        FUNC this is a function application: parenthesize if added (?)
        DIV  this is a division: make wider division if divided
        POW  this is a power: only parenthesize if exponent
        MUL  this is a multiplication: parenthesize if powered
        ADD  this is an addition: parenthesize if multiplied or powered
        NEG  this is a negative number: optimize if added, parenthesize if
             multiplied or powered
        OPEN this is an open object: parenthesize if added, multiplied, or
             powered (example: Piecewise)
    """

    ATOM, FUNC, DIV, POW, MUL, ADD, NEG, OPEN = range(8)

    def __init__(self, s, baseline=0, binding=0, unicode=None):
        """Initialize from stringPict and binding power."""
        assert isinstance(s, str)
        assert isinstance(baseline, int)
        stringPict.__init__(self, s, baseline)
        self.binding = binding
        if unicode is not None:
            sympy_deprecation_warning(
                """
                The unicode argument to prettyForm is deprecated. Only the s
                argument (the first positional argument) should be passed.
                """,
                deprecated_since_version="1.7",
                active_deprecations_target="deprecated-pretty-printing-functions")
        self._unicode = unicode or s

    @property
    def unicode(self):
        sympy_deprecation_warning(
            """
            The prettyForm.unicode attribute is deprecated. Use the
            prettyForm.s attribute instead.
            """,
            deprecated_since_version="1.7",
            active_deprecations_target="deprecated-pretty-printing-functions")
        return self._unicode

    # Note: code to handle subtraction is in _print_Add

    def __add__(self, *others):
        """Make a pretty addition.
        Addition of negative numbers is simplified.
        """
        arg = self
        if arg.binding > prettyForm.NEG:
            arg = arg.parenthesis()
        result = [arg]
        for arg in others:
            # add parentheses for weak binders
            if arg.binding > prettyForm.NEG:
                arg = arg.parenthesis()
            #use existing minus sign if available
            if arg.binding != prettyForm.NEG:
                result.append(' + ')
            result.append(arg)
        return prettyForm(binding=prettyForm.ADD, *stringPict.next(*result))

    def __truediv__(self, den, slashed=False):
        """Make a pretty division; stacked or slashed."""
        if slashed:
            raise NotImplementedError("Can't do slashed fraction yet")
        num = self
        if num.binding == prettyForm.DIV:
            num = num.parenthesis()
        if den.binding == prettyForm.DIV:
            den = den.parenthesis()

        if num.binding == prettyForm.NEG:
            num = num.right(' ')

        return prettyForm(
            binding=prettyForm.DIV, *stringPict.stack(num, stringPict.LINE, den)
        )

    def __mul__(self, *others):
        """Make a pretty multiplication.
        Parentheses are needed around +, - and neg.
        """
        quantity = {"degree": "\N{DEGREE SIGN}"}

        if len(others) == 0:
            return self  #We aren't actually multiplying... So nothing to do here.

        #add parens on args that need them
        arg = self
        if arg.binding > prettyForm.MUL and arg.binding != prettyForm.NEG:
            arg = arg.parenthesis()
        result = [arg]
        for arg in others:
            if arg.picture[0] not in quantity.values():
                result.append(xsym('*'))
            #add parentheses for weak binders
            if arg.binding > prettyForm.MUL and arg.binding != prettyForm.NEG:
                arg = arg.parenthesis()
            result.append(arg)

        len_res = len(result)
        for i in range(len_res):
            if i < len_res - 1 and result[i] == "-1" and result[i + 1] == xsym('*'):
                # substitute -1 by -, like in -1*x -> -x
                result.pop(i)
                result.pop(i)
                result.insert(i, '-')
        if result[0][0] == '-':
            # if there is a - sign in front of all
            # This test was failing to catch a prettyForm.__mul__(prettyForm("-1", 0, 6)) being negative
            bin = prettyForm.NEG
            if result[0] == '-':
                right = result[1]
                if right.picture[right.baseline][0] == '-':
                    result[0] = '- '
        else:
            bin = prettyForm.MUL
        return prettyForm(binding=bin, *stringPict.next(*result))

    def __repr__(self):
        return "prettyForm(%r,%d,%d)" % (
            '\n'.join(self.picture),
            self.baseline,
            self.binding,
        )

    def __pow__(self, b):
        """Make a pretty power."""
        a = self
        use_inline_func_form = False
        if b.binding == prettyForm.POW:
            b = b.parenthesis()
        if a.binding > prettyForm.FUNC:
            a = a.parenthesis()
        elif a.binding == prettyForm.FUNC:
            # heuristic for when to use inline power
            if b.height() > 1:
                a = a.parenthesis()
            else:
                use_inline_func_form = True

        if use_inline_func_form:
            #         2
            #  sin  +   + (x)
            func = a.prettyFunc.superindex(b)
            return func.right(a.prettyArgs)
        else:
            result = a.superindex(b)
            result.binding = prettyForm.POW
            return result

    simpleFunctions = ["sin", "cos", "tan"]

    @staticmethod
    def apply(function, *args):
        """Functions of one or more variables."""
        if function in prettyForm.simpleFunctions:
            #simple function: use only space if possible
            assert len(args) == 1, "Simple function %s must have 1 argument" % function
            arg = args[0].__pretty__()
            if arg.binding <= prettyForm.DIV:
                #optimization: no parentheses necessary
                result = arg.left(function + ' ')
                result.binding = prettyForm.FUNC
                return result
        argumentList = []
        for arg in args:
            argumentList.append(',')
            argumentList.append(arg.__pretty__())
        if argumentList:
            first_arg, *rest = argumentList
            argument_pict = first_arg.right(*rest)
        else:
            argument_pict = stringPict("")
        argument_pict = argument_pict.parenthesis()
        func_pict = argument_pict.left(function)
        func_pict.binding = prettyForm.ATOM
        return func_pict
