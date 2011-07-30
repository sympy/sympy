class DataMatrix(object):

    def _format_str(self, strfunc, rowsep='\n'): 
        # Handle zero dimensions:
        if self.rows == 0 or self.cols == 0:
            return '[]'
        # Build table of string representations of the elements
        res = []
        # Track per-column max lengths for pretty alignment
        maxlen = [0] * self.cols
        for i in range(self.rows):
            res.append([])
            for j in range(self.cols):
                string = strfunc(self[i,j])
                res[-1].append(string)
                maxlen[j] = max(len(string), maxlen[j])
        # Patch strings together
        for i, row in enumerate(res):
            for j, elem in enumerate(row):
                # Pad each element up to maxlen so the columns line up
                row[j] = elem.rjust(maxlen[j])
            res[i] = "[" + ", ".join(row) + "]"
        return rowsep.join(res)

    def __str__(self):
        from sympy.printing.str import sstr
        return sstr(self)

    __repr__ = __str__

    def __len__(self):
        """
        Return the number of elements of self.

        Implemented mainly so bool(Matrix()) == False.
        """
        return self.rows * self.cols

    def __array__(self):
        return matrix2numpy(self)

    def hash(self): # Remove this ?
        """Compute a hash every time, because the matrix elements
        could change."""
        return hash(self.__str__() )

    @property
    def shape(self):
        return (self.rows, self.cols)

    def __mathml__(self):
        mml = ""
        for i in range(self.rows):
            mml += "<matrixrow>"
            for j in range(self.cols):
                mml += self[i,j].__mathml__()
            mml += "</matrixrow>"
        return "<matrix>" + mml + "</matrix>"
