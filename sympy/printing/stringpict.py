"""Prettyprinter by Jurjen Bos.
(I hate spammers: mail me at pietjepuk314 at the reverse of ku.oc.oohay).
All objects have a method that create a "stringPict",
that can be used in the str method for pretty printing.

Updates by Jason gedge (email <my last name> at cs mun ca)
    - terminal_string() method
    - minor fixes and changes (mostly to prettyForm)

TODO:
    - Allow left/center/right alignment options for above/below and
      top/center/bottom alignment options for left/right
"""

class stringPict:
    """A ASCII picture.
    The pictures are represented as a list of equal length strings.
    """
    #special value for stringPict.below
    LINE = 'line'

    def __init__(self, s, baseline=0):
        """Initialize from string.
        Multiline strings are centered.
        """
        #picture is a string that just can be printed
        self.picture = stringPict.equalLengths(s.splitlines())
        #baseline is the line number of the "base line"
        self.baseline = baseline
        self.binding = None

    def __len__(self):
        return len(str(self))

    @staticmethod
    def equalLengths(lines):
        width = max(len(line) for line in lines)
        return [line.center(width) for line in lines]

    def height(self):
        return len(self.picture)

    def width(self):
        return len(self.picture[0])

    @staticmethod
    def next(*args):
        """Put a string of stringPicts next to each other.
        Returns string, baseline arguments for stringPict.
        """
        #convert everything to stringPicts
        objects = []
        for arg in args:
            if isinstance(arg, str): arg = stringPict(arg)
            objects.append(arg)

        #make a list of pictures, with equal height and baseline
        newBaseline = max(obj.baseline for obj in objects)
        newHeightBelowBaseline = max(
            obj.height()-obj.baseline
            for obj in objects)
        newHeight = newBaseline + newHeightBelowBaseline

        pictures = []
        for obj in objects:
            oneEmptyLine = [' '*obj.width()]
            basePadding = newBaseline-obj.baseline
            totalPadding = newHeight-obj.height()
            pictures.append(
                oneEmptyLine * basePadding +
                obj.picture +
                oneEmptyLine * (totalPadding-basePadding))

        result = [''.join(lines) for lines in zip(*pictures)]
        return '\n'.join(result), newBaseline

    def right(self, *args):
        """Put pictures next to this one.
        Returns string, baseline arguments for stringPict.
        (Multiline) strings are allowed, and are given a baseline of 0.
        >>> print stringPict("10").right(" + ",stringPict("1\r-\r2",1))[0]
             1
        10 + -
             2
        """
        return stringPict.next(self, *args)

    def left(self, *args):
        """Put pictures (left to right) at left.
        Returns string, baseline arguments for stringPict.
        """
        return stringPict.next(*(args+(self,)))

    @staticmethod
    def stack(*args):
        """Put pictures on top of each other,
        from top to bottom.
        Returns string, baseline arguments for stringPict.
        The baseline is the baseline of the second picture.
        Everything is centered.
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
        newWidth = max(
            obj.width()
            for obj in objects
            if obj is not stringPict.LINE)

        lineObj = stringPict('-'*newWidth)

        #replace LINE with proper lines
        for i, obj in enumerate(objects):
            if obj is stringPict.LINE:
                objects[i] = lineObj

        #stack the pictures, and center the result
        newPicture = []
        for obj in objects:
            newPicture.extend(obj.picture)
        newPicture = [line.center(newWidth) for line in newPicture]
        newBaseline = objects[0].height()+objects[1].baseline
        return '\n'.join(newPicture), newBaseline

    def below(self, *args):
         """Put pictures under this picture.
         Returns string, baseline arguments for stringPict.
         Baseline is baseline of top picture
         >>> print stringPict("x+3").below(stringPict.LINE, '3')[0] #doctest: +NORMALIZE_WHITESPACE
         x+3
         ---
          3
         """
         s, baseline = stringPict.stack(self, *args)
         return s, self.baseline

    def above(self, *args):
        """Put pictures above this picture.
        Returns string, baseline arguments for stringPict.
        Baseline is baseline of bottom picture.
        """
        string, baseline = stringPict.stack(*(args+(self,)))
        baseline = len(string.splitlines())-self.height()+self.baseline
        return string, baseline

    def parens(self):
        """Put parentheses around self.
        Returns string, baseline arguments for stringPict.
        """
        height = self.height()
        if height==1:
            return stringPict('(').right(self, ')')
        else:
            verticalBar = '\n' + '|\n' * (self.height()-2)
            lparen = stringPict('/'+verticalBar+'\\',self.baseline)
            rparen = stringPict('\\'+verticalBar+'/',self.baseline)
            return lparen.right(self, rparen)

    def leftslash(self):
        """Precede object by a slash of the proper size.
        """
        height = max(
            self.baseline,
            self.height()-1-self.baseline)*2 + 1
        slash = '\n'.join(
            ' '*(height-i-1)+'/'+' '*i
            for i in range(height)
            )
        return self.left(stringPict(slash, height//2))

    def root(self, n=None):
        """Produce a nice root symbol.
        Produces ugly results for big n inserts.
        """
        #put line over expression
        result = self.above('_'*self.width())
        #construct right half of root symbol
        height = self.height()
        slash = '\n'.join(
            ' ' * (height-i-1) + '/' + ' ' * i
            for i in range(height)
            )
        slash = stringPict(slash, height-1)
        #left half of root symbol
        if height > 2:
            downline = stringPict('\\ \n \\',1)
        else:
            downline = stringPict('\\')
        #put n on top, as low as possible
        if n is not None and n.width()>downline.width():
            downline = downline.left(' '*(n.width()-downline.width()))
            downline = downline.above(n)
        #build root symbol
        root = downline.right(slash)
        #glue it on at the proper height
        #normally, the root symbel is as high as self
        #which is one less than result
        #this moves the root symbol one down
        #if the root became higher, the baseline has to grow too
        root.baseline = result.baseline-result.height()+root.height()
        return result.left(root)

    def terminal_string(self):
        """Return the string form of self such that it can be printed
           on the terminal without being broken up.
        """

        # Attempt to get a terminal width
        ncols = 0
        try:
            import curses
            curses.setupterm()
            ncols = curses.tigetnum('cols')
        except (ImportError, TypeError):
            pass
        
        ncols -= 2
        if ncols <= 0:
            ncols = 78

        # If smaller than the terminal width, no need to correct
        if self.width() <= ncols:
            return type(self.picture[0])(self)

        i = 0
        svals = []
        while i < self.width():
            svals.extend([ sval[i:i+ncols] for sval in self.picture ])
            svals.append("") # a vertical spacer
            i += ncols

        del svals[-1] #  Get rid of the last spacer
        _str = type(self.picture[0])
        return _str.join(_str("\n"), svals)

    def __eq__(self, o):
        if isinstance(o, str):
            return '\n'.join(self.picture) == o
        elif isinstace(o, stringPict):
            return o.picture == self.picture
        return False

    def __str__(self):
        return str.join('\n', self.picture)

    def __unicode__(self):
        return unicode.join(u'\n', self.picture)

    def __repr__(self):
        return "stringPict(%r,%d)"%('\n'.join(self.picture), self.baseline)

    def __getitem__(self, index):
        return self.picture[index]

    def __len__(self):
        return len(self.s)

class prettyForm(stringPict):
    """Extension of the stringPict class that knows about
    basic math applications, optimizing double minus signs.
    "Binding" is interpreted as follows:
    ATOM this is an atom: never needs to be parenthesised
    DIV  this is a division: make wider division if divided 
    POW  this is a power: only parenthesise if exponent
    MUL  this is a multiplication: parenthesise if powered
    ADD  this is an addition: parenthesise if multiplied or powered
    NEG  this is a negative number: optimise if added, parenthesise if multiplied or powered
    FUNC this is a function application: parenthesise if multiplied, powered, or added
    """
    ATOM, DIV, POW, MUL, ADD, NEG, FUNC = range(7)

    def __init__(self, s, baseline=0, binding=0, unicode=None):
        """Initialize from stringPict and binding power."""
        stringPict.__init__(self, s, baseline)
        self.binding = binding
        self.unicode = unicode or s

    def __add__(self, *others):
        """Make a pretty addition.
        Addition of negative numbers is simplified.
        """
        arg = self
        if arg.binding > prettyForm.NEG: arg = stringPict(*arg.parens())
        result = [arg]
        for arg in others:
            #add parentheses for weak binders
            if arg.binding > prettyForm.NEG: 
                arg = stringPict(*arg.parens())
            #use existing minus sign if available
            if arg.binding != prettyForm.NEG:
                result.append(' + ')
            result.append(arg)
        return prettyForm(binding=prettyForm.ADD, *stringPict.next(*result))

    def __div__(self, den, slashed = False):
        """Make a pretty division; stacked or slashed.
        """
        if slashed: raise NotImplentedError, "Can't do slashed fraction yet"
        num = self
        if num.binding==prettyForm.DIV: num = stringPict(*num.parens())
        if den.binding==prettyForm.DIV: den = stringPict(*den.parens())

        return prettyForm(binding=prettyForm.DIV, *stringPict.stack(
            num,
            stringPict.LINE,
            den))

    def __mul__(self, *others):
        """Make a pretty multiplication.
        Parentheses are needed around +, - and neg.
        """
        args = self
        if args.binding > prettyForm.MUL: arg = stringPict(*args.parens())
        result = [args]
        for arg in others:
            result.append('*')
            #add parentheses for weak binders
            if arg.binding > prettyForm.MUL: arg = stringPict(*arg.parens())
            result.append(arg)
        len_res = len(result)
        for i in xrange(len_res):
            if i < len_res-1 and result[i] == '-1' and result[i+1] == "*":
                # substitute -1 by -, like in -1*x -> -x
                result.pop(i)
                result.pop(i)
                result.insert(i, '-')
        if result[0][0] == '-':
            # if there is a - sign in front of all
            bin = prettyForm.NEG
        else:
            bin = prettyForm.MUL
        return prettyForm(binding=bin, *stringPict.next(*result))

    def __repr__(self):
        return "prettyForm(%r,%d,%d)"%(
            '\n'.join(self.picture),
            self.baseline,
            self.binding)

    def __pow__(self, b):
        """Make a pretty power.
        """
        a = self
        if a.binding > prettyForm.ATOM: a = stringPict(*a.parens())
        if b.binding == prettyForm.POW: b = stringPict(*b.parens())
        top = stringPict(*b.left(' '*a.width()))
        bottom = stringPict(*a.right(' '*b.width()))
        return prettyForm(binding=prettyForm.POW, *bottom.above(top))

    simpleFunctions = ["sin", "cos", "tan"]
    @staticmethod
    def apply(function, *args):
        """Functions of one or more variables.
        """
        if function in prettyForm.simpleFunctions:
            #simple function: use only space if possible
            assert len(args)==1, "Simple function %s must have 1 argument"%function
            arg = args[0].__pretty__()
            if arg.binding <= prettyForm.DIV:
                #optimization: no parentheses necessary
                return prettyForm(binding=prettyForm.FUNC, *arg.left(function+' '))
        argumentList = []
        for arg in args:
            argumentList.append(',')
            argumentList.append(arg.__pretty__())
        argumentList = stringPict(*stringPict.next(*argumentList[1:]))
        argumentList = stringPict(*argumentList.parens())
        return prettyForm(binding=prettyForm.ATOM, *argumentList.left(function))
