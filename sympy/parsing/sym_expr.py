from sympy.parsing.fortran.fortran_parser import src_to_sympy
from sympy.printing import pycode, ccode, fcode
from sympy.core.basic import Basic

class SymPyExpression(Basic):
    """Class to store and handle SymPy expressions

    This class will hold SymPy Expressions and handle the API for the conversion to and from different languages.

    It works with the C and the Fortran Parser to generate SymPy expressions which are stored here and which can be converted to multiple language's source code.

    Notes
    =====

    The module and its API are currently under development and experimental and can be changed during development.

    Examples
    ========
    >>> src = '''\
    ... integer function f(a,b)
    ... integer, intent(in) :: a, b
    ... integer :: r
    ... end function
    ... '''
    >>> a = SymPyExpression(src, 'f')
    >>> a.convert_to_python()
    ['def f(a, b):\n    f = 0\n    r = 0\n    return f']


    """
    def __init__(self, source_code = None, mode = None):
        """Constructor for SymPyExpression class"""
        super(SymPyExpressions, self).__init__()
        if mode == 'f' or mode == 'F':
            if source_code:
                self._expr = src_to_sympy(source_code)
            else:
                self._expr = []
        elif mode == 'c'or mode == 'C':
            if source_code:
                self._expr = src_to_c(source_code)
            else:
                self._expr = []
        else:
            self._expr = []

    def convert_to_expr(self, src_code, mode):
        """Converts the given source code to sympy Expressions

        Attributes
        ==========

        src_code : String
            the source code or filename of the source code that is to be converted
        mode: String
            the mode to determine which parser is to be used according to the language of the source code

        Examples
        ========

        >>> src = '''\
        ... integer function f(a,b)
        ... integer, intent(in) :: a, b
        ... integer :: r
        ... end function
        ... '''
        >>> a = SymPyExpression()
        >>> a.convert_to_expr(src,'f')
        >>> a.convert_to_python()
        ['def f(a, b):\n    f = 0\n    r = 0\n    return f']

        """
        if mode.lower() == 'f':
            self._expr = src_to_sympy(src_code)
        elif mode.lower() == 'c':
            self._expr = src_to_c(src_code)
        else:
            raise NotImplementedError("The langauge parser has not been implemented. Invalid Input!")

    def convert_to_python(self):
        """Returns a list with python code for the sympy expressions"""
        self._pycode = []
        for iter in self._expr:
            self._pycode.append(pycode(iter))
        return self._pycode

    def convert_to_c(self):
        """Returns a list with the c source code for the sympy expressions"""
        self._ccode = []
        for iter in self._expr:
            self._ccode.append(ccode(iter))
        return self._ccode

    def convert_to_fortran(self):
        """Returns a list with the fortran source code for the sympy expressions"""
        self._fcode = []
        for iter in self._expr:
            self._fcode.append(fcode(iter))
        return self._fcode
